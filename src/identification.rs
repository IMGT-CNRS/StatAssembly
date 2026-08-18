/*
* IMGT/StatAssembly
Available under EUPL license
Made by: Guilhem Zeitoun and IMGT Team
*/
use bio::io::fasta;
use itertools::Itertools;
use rammap::Aligner;
use reqwest::{StatusCode, header};

use crate::{
    BORNES, LOCUSSEPARATOR, TIMEOUT_IN_MN, getassemblyreader, printpotentialbornes,
    r#struct::{
        Args, Blast, Blastcalc, Blastlevel, Blastmatch, Filecrea, Haplotype, Locus, LocusHaplo,
        LocusInfos, Name, Position, Species, Status, Strand,
    },
    submissions::{
        BORNESLINK, MOTIFLINK, RELEASELINK, REQUESTCLIENT, VQUESTLINK, blastcommand, checkoverlap,
        fileincache, getprogressbar, speciesandorphonfiltering,
    },
};
use std::{
    cmp::{Ordering, max, min},
    collections::{BTreeMap, BTreeSet, HashMap},
    env,
    fs::File,
    io::{
        self, BufReader, Cursor,
        ErrorKind::{self, InvalidFilename},
        Read, Write,
    },
    path::{Path, PathBuf},
    str::FromStr,
    sync::{Arc, mpsc::channel},
    thread::{self, ScopedJoinHandle, sleep},
    time::Duration,
};
#[allow(clippy::type_complexity)]
pub(crate) fn threadlaunch<T>(
    referencepath: Filecrea,
    subject: T,
    args: &Args,
    bornespath: Option<Filecrea>,
) -> io::Result<(Vec<Blast>, Option<Option<Vec<Blastmatch>>>)>
where
    T: AsRef<Path> + Clone + Send + Sync,
{
    thread::scope(|s| {
        let subjectbis = subject.as_ref().to_path_buf();
        let _ = File::open(subject.as_ref())?; //Check subject can be opened
        let blast = s.spawn(move || {
            let b = match blastcommand(referencepath, subjectbis.into(), Blastlevel::Normal) {
                Ok(b) => b,
                Err(e) => {
                    return Err(io::Error::new(ErrorKind::InvalidData, e));
                }
            };
            Ok(b)
        });
        let (sender, receiver) = channel();
        let sender = Arc::new(sender);
        let bornes = bornespath.as_ref().map(|b| {
            s.spawn(move || {
                match receiver.recv_timeout(Duration::from_mins(TIMEOUT_IN_MN)) {
                    Ok("kill") => {
                        eprintln!("aborted by main");
                        return None;
                    }
                    Ok(_) => (),
                    Err(e) => {
                        eprintln!("Error is {e}");
                        return None;
                    }
                };
                let aligner = Aligner::from_fasta(
                    subject.as_ref().display().to_string().as_str(),
                    rammap::Preset::Asm20,
                )
                .map_err(|d| {
                    eprintln!("Unable to set aligner. Error is {d}");
                })
                .ok()?;
                let options = rammap::MapOpts {
                    cs: Some(true),
                    ..Default::default()
                };
                let read = fasta::Reader::from_file(b.getpath()).ok()?;
                let mut aligns: Vec<Blastmatch> = Vec::new();
                for seq in read.records().filter_map(Result::ok) {
                    let name = format!("NCBI|{}", seq.desc().unwrap_or_default());
                    if let Ok("kill") = receiver.try_recv() {
                        eprintln!("aborted by main");
                        return None;
                    };
                    if let a = aligner
                        .map_seq_with(
                            name.as_str(),
                            seq.seq(),
                            options,
                        )
                        .mappings
                        .into_iter()
                        .find(|f| f.is_primary) //Should be only one primary
                        && let Some(first) = a
                        && let (
                            Ok(qseqid),
                            sseqid,
                            sseq,
                            Ok(sstart),
                            Ok(send),
                            complement,
                            status,
                            identity,
                        ) = (
                            Name::from_str(&name),
                            &first.target_name,
                            String::new(),
                            TryInto::<usize>::try_into(first.target_start),
                            TryInto::<usize>::try_into(first.target_end),
                            if first.strand == rammap::Strand::Reverse {
                                Strand::Minus
                            } else {
                                Strand::Plus
                            },
                            Status::New,
                            first.score,
                        )
                    {
                        let matche = Blastmatch::new(
                            qseqid.clone(),
                            sseqid.to_string(),
                            sseq,
                            sstart,
                            send,
                            complement,
                            status,
                            f32::from_str(identity.to_string().as_str()).ok()?,
                        );
                        aligns.push(matche);
                    }
                }
                println!("Bornes finished.");
                if aligns.is_empty() {
                    None
                } else {
                    Some(aligns)
                }
            })
        });
        let threadtimeout = |bornessome: &ScopedJoinHandle<'_, Option<Vec<Blastmatch>>>| {
            thread::scope(|f| {
                let sender = Arc::clone(&sender);
                f.spawn(move || {
                    let mut time = 0;
                    while time < TIMEOUT_IN_MN.saturating_mul(60) {
                        if bornessome.is_finished() {
                            return;
                        }
                        sleep(Duration::new(1, 0));
                        time += 1;
                    }
                    let _ = sender.send("kill");
                });
            });
        };
        println!("Waiting for blast and bornes to finish. It can take several minutes.");
        let blast = if args.lowmemory {
            println!("You have low memory, BLAST would be launched first.");
            let r = blast
                .join()
                .map_err(|_| io::Error::new(ErrorKind::BrokenPipe, "Error with blast thread"))??;
            println!("Starting minimap2 now.");
            let _ = sender.send("start"); //No need to check, it will turn into Err on receiver and kill the process
            if let Some(bornessome) = &bornes {
                threadtimeout(bornessome);
            }
            r
        } else {
            let _ = sender.send("start"); //No need to check, it will turn into Err on receiver and kill the process
            if let Some(bornessome) = &bornes {
                threadtimeout(bornessome);
            }
            blast
                .join()
                .map_err(|_| io::Error::new(ErrorKind::BrokenPipe, "Error with blast thread"))??
        };
        let bornes = bornes.map(|d| {
            d.join()
                .map_err(|_| io::Error::new(ErrorKind::BrokenPipe, "Error with bornes thread"))
                .ok()?
        });
        let _ = sender.send("kill"); // just in case, should be useless
        Ok((blast, bornes))
    })
}
pub(crate) fn locusallposition(
    species: &Species,
    args: &Args,
) -> io::Result<(Vec<LocusInfos>, Vec<Blastmatch>, String)> {
    let assembly = args
        .assembly
        .as_ref()
        .ok_or(io::Error::new(InvalidFilename, "No assembly given"))?;
    let infos = if args.nobornes {
        None
    } else {
        downloadbornes(args)
    };
    let (bornespath, referencepath, releaseversion) = match (
        infos,
        downloadref(true, args.cachebypass, &None).map(|(mut a, b)| {
            (
                speciesandorphonfiltering(&mut a, None, b.clone(), species, true, args.cachebypass),
                b,
            )
        }),
    ) {
        (Some(a), Some((b, c))) => (Some(a), b?, c),
        (None, Some((b, c))) if args.nobornes => (None, b?, c),
        (None, Some((b, c))) if !args.nobornes => {
            eprintln!("Bornes could not downloaded, skipped.");
            (None, b?, c)
            /* return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Bornes from IMGT cannot be downloaded",
            )); */
        }
        _ => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    let info = threadlaunch(referencepath, &assembly, args, bornespath);
    let (mut blast, mut bornes) = match info {
        Ok((a, Some(Some(b)))) => (a, Some(b)),
        Ok((a, None)) => {
            if !args.nobornes {
                eprintln!("No bornes found.");
            }
            (a, None)
        }
        Ok((_, Some(None))) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                "Bornes could not be analyzed".to_string(),
            ));
        }
        Err(a) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!(
                    "Blast and bornes could not be performed. Error are: {:?}",
                    a
                ),
            ));
        }
    };
    if let Some(p) = bornes.as_mut() {
        p.sort_unstable_by(|f, g| match f.getallelename().cmp(g.getallelename()) {
            Ordering::Equal => f.getidentity().total_cmp(&g.getidentity()).reverse(),
            e => e,
        });
    }
    //Filter by locus
    retainbestmatch(&mut blast);
    blast.shrink_to_fit();
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
        } else {
            p.complement = Strand::Plus;
        }
        p.setstatus();
    });
    let mut statusvec: BTreeMap<(String, usize, usize), Blast> = BTreeMap::new();
    for elem in blast.into_iter() {
        if let Some((k, b)) = statusvec.iter().find(|((s, r1, r2), b)| {
            b != &&elem
                && s.as_str() == elem.sseqid.as_str()
                && checkoverlap(&(*r1..=*r2), &(elem.sstart..=elem.send))
        }) {
            let scoreactual = b.length / b.qlen * b.pident.powi(2) as usize;
            let newscore = elem.length / elem.qlen * elem.pident.powi(2) as usize;
            if newscore > scoreactual {
                statusvec.insert(k.clone(), elem); //Erase all value
            }
        } else {
            statusvec.insert((elem.sseqid.clone(), elem.sstart, elem.send), elem);
        }
    }
    let mut data: Vec<Blastmatch> = statusvec.into_values().map(|v| v.into()).collect();
    //Sort by name then starting then ending position
    data.sort_unstable();
    data.dedup();
    data.shrink_to_fit();
    let range = find_global_best_range(&data, &bornes, &args).ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No locus found after BLAST analysis",
    ));
    let infos = getassemblyreader(args)?;
    //Set to maximum length of assembly if the locus is at the end of the sequence
    let range = match range {
        Ok(mut b) => {
            b.iter_mut().for_each(|a| {
                let max = infos
                    .index
                    .sequences()
                    .iter()
                    .find(|c| &c.name == a.getcontig())
                    .map(|b| b.len);
                if let Some(max) = max
                    && let Ok(maxi) = max.try_into()
                    && a.end.getzbasedpos() > maxi
                {
                    a.end = Position::newfromzposition(maxi);
                }
            });
            Ok(b)
        }
        Err(a) => Err(a),
    };
    range.map(|f| (f, data, releaseversion))
}
pub(crate) fn retainbestmatch(blast: &mut Vec<Blast>) {
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
}
pub(crate) fn downloadbornes(args: &Args) -> Option<Filecrea> {
    println!("Downloading bornes sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "bornes.fasta");
    if args.cachebypass || !fileincache(&tempfile) {
        println!("Downloading IMGT bornes, please wait...");
        match sendresultcompressed(&REQUESTCLIENT, BORNESLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    _ => {
                        eprintln!("Cannot write bornes in sequence.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT bornes data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!("IMGT bornes were already downloaded, retrieving...");
    };
    Some(Filecrea::createfrompath(tempfile))
}
pub(crate) fn downloadmotifs(args: &Args) -> Option<Filecrea> {
    println!("Downloading motifs sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "motifs.fasta");
    let file = Filecrea::createfrompath(&tempfile);
    if args.cachebypass || !fileincache(&tempfile) {
        println!("Downloading IMGT motifs, please wait...");
        match sendresultcompressed(&REQUESTCLIENT, MOTIFLINK.as_str()) {
            Ok(e) => {
                match file.setfile().and_then(|mut f| f.write_all(e.as_bytes())) {
                    Ok(_) => (),
                    _ => {
                        eprintln!("Cannot write refseq in sequence.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT motifs data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!("IMGT motifs were already downloaded, retrieving...");
    };
    Some(file)
}
pub(crate) fn find_global_best_range(
    blastcheck: &[Blastmatch],
    bornes: &Option<Vec<Blastmatch>>,
    args: &Args,
) -> Option<Vec<LocusInfos>> {
    let mut groups: HashMap<LocusHaplo, (String, bool, usize, usize, bool)> = HashMap::new();
    let blastcheck = blastcheck.to_vec();
    let mut hash = blastcheck
        .iter()
        .filter_map(|p| p.getlocusname().map(|d| (d, p)))
        .into_group_map_by(|(l, f)| (l.clone(), f.getsubject()));
    let mut bornes = match bornes {
        Some(bornes) => Some(bornes.iter().into_group_map_by(|f| f.getsubject())),
        None => None,
    };
    if let Some(p) = &mut bornes {
        p.values_mut().for_each(|d| d.sort_unstable());
    }
    hash.iter_mut().for_each(|(_, v)| {
        v.sort_unstable_by(|(_, a), (_, b)| match a.getsubject().cmp(b.getsubject()) {
            Ordering::Equal => match a.getpos().0.cmp(&b.getpos().0) {
                Ordering::Equal => b.getpos().1.cmp(&b.getpos().1).reverse(),
                ord => ord,
            },
            ord => ord,
        })
    });
    let mut locusborne: HashMap<
        (Locus, String, u32),
        (
            std::collections::BTreeSet<Blastmatch>,
            (Option<Blastmatch>, Option<Blastmatch>),
        ),
    > = HashMap::new();
    for ((loci, sseqname), data) in hash.iter() {
        let mut index = 0;
        for (_, elem) in data {
            let d = match locusborne.get_mut(&(loci.clone(), sseqname.to_string(), index)) {
                Some(a) => a,
                None => {
                    let mut e = BTreeSet::new();
                    e.insert((*elem).clone());
                    locusborne.insert(
                        (loci.clone(), sseqname.to_string(), index),
                        (e, (None, None)),
                    );
                    continue;
                }
            };
            if let Some(b) = d.0.last()
                && std::cmp::max(b.sstart, b.send).abs_diff(std::cmp::min(elem.send, elem.sstart))
                    <= LOCUSSEPARATOR
            {
                d.0.insert((*elem).clone());
            } else {
                index += 1;
                let mut e = BTreeSet::new();
                e.insert((*elem).clone());
                locusborne.insert(
                    (loci.clone(), sseqname.to_string(), index),
                    (e, (None, None)),
                );
            }
        }
    }
    if let Some(borneblast) = bornes {
        for (key, (val, (bornefirst, bornesecond))) in locusborne.iter_mut() {
            let itera = val.clone();
            let (loc, min, max) = match (key, itera.iter().min(), itera.iter().min()) {
                (c, Some(a), Some(b)) => (c, a, b),
                _ => continue,
            };
            let mut raw = borneblast
                .iter()
                .filter(|(key, _)| **key == loc.1.as_str())
                .flat_map(|(_, hit)| hit.into_iter());
            let findfirstborne = raw.find(|val| {
                val.sstart.abs_diff(min.sstart) <= LOCUSSEPARATOR
                    || val.sstart.abs_diff(max.send) <= LOCUSSEPARATOR
                    || val.send.abs_diff(max.send) <= LOCUSSEPARATOR
                    || val.send.abs_diff(min.sstart) <= LOCUSSEPARATOR
            });
            let findsecondborne = raw.find(|val| {
                findfirstborne.is_none_or(|f| f.getallelename() != val.getallelename()) //Do not take the same borne twice
                    && (val.sstart.abs_diff(min.sstart) <= LOCUSSEPARATOR
                        || val.sstart.abs_diff(max.send) <= LOCUSSEPARATOR
                        || val.send.abs_diff(max.send) <= LOCUSSEPARATOR
                        || val.send.abs_diff(min.sstart) <= LOCUSSEPARATOR)
            });
            let data = |borne: Option<&&Blastmatch>,
                        toput: &mut Option<Blastmatch>,
                        min: &Blastmatch,
                        max: &Blastmatch,
                        val: &mut BTreeSet<Blastmatch>| {
                if let Some(hit) = borne {
                    *toput = Some((*hit).clone());
                    let range = min.sstart..=max.send;
                    if range.contains(&hit.sstart) || range.contains(&hit.send) {
                        if hit.send.saturating_sub(min.sstart) < hit.sstart.saturating_sub(max.send)
                        {
                            //Cutting after the borne
                            //Bornes is at the beggining of the locus
                            val.retain(|p| p.sstart > hit.sstart);
                        } else {
                            //Bornes is at the end of the locus
                            val.retain(|p| p.send < hit.send);
                        }
                    } else if hit.sstart.abs_diff(min.sstart) <= BORNES
                        || hit.sstart.abs_diff(max.send) <= BORNES
                        || hit.send.abs_diff(max.send) <= BORNES
                        || hit.send.abs_diff(min.sstart) <= BORNES
                    {
                        val.insert((*hit).clone());
                    }
                };
            };
            data(findfirstborne, bornefirst, &min, &max, val);
            data(findsecondborne, bornesecond, &min, &max, val);
        }
    }
    if locusborne.is_empty() {
        eprintln!("Cannot find any locus, aborting.");
        return None;
    }
    let real = locusborne.into_iter().sorted_unstable_by(
        |((al, ass, _), (av, _)), ((bl, bs, _), (bv, _))| match al.cmp(bl) {
            Ordering::Equal => match av.iter().count().cmp(&bv.iter().count()) {
                Ordering::Equal => ass.cmp(bs),
                ord2 => ord2.reverse(),
            },
            ord => ord,
        },
    );
    let bornesscript = |bornes: (Option<Blastmatch>, Option<Blastmatch>),
                        locus: &LocusHaplo,
                        bornestoprint: &mut Vec<Blastmatch>| {
        let addlocustoname = |name: &mut Blastmatch| {
            name.qseqid.gene = format!("{}_{}", name.qseqid.gene, locus.to_string());
        };
        match bornes {
            (Some(mut first), Some(mut second)) => {
                addlocustoname(&mut first);
                bornestoprint.push(first);
                addlocustoname(&mut second);
                bornestoprint.push(second);
                true
            }
            (Some(mut first), None) => {
                addlocustoname(&mut first);
                bornestoprint.push(first);
                true
            }
            (None, Some(_)) => unreachable!("First borne cannot be null if second borne"),
            (None, None) => false,
        }
    };
    let mut bornestoprint = Vec::new();
    for ((loci, sseq, _), (blast, bornes)) in real {
        let primaryhap = LocusHaplo::new(loci.clone(), Haplotype::Primary);
        let secondhap = LocusHaplo::new(loci.clone(), Haplotype::Alternate);
        if !groups.contains_key(&primaryhap) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let hasbornes = bornesscript(bornes, &primaryhap, &mut bornestoprint);
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert(primaryhap, (sseq, complement, min, max, hasbornes));
        } else if !groups.contains_key(&LocusHaplo::new(loci.clone(), Haplotype::Alternate)) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let hasbornes = bornesscript(bornes, &secondhap, &mut bornestoprint);
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert(secondhap, (sseq, complement, min, max, hasbornes));
        }
    }
    let mut loc: Vec<LocusInfos> = groups
        .into_iter()
        .map(|(haplo, (sseq, complement, start, end, hasbornes))| {
            LocusInfos::newfromlocushap(
                haplo,
                sseq,
                Position::new(
                    false,
                    start
                        .checked_sub(if hasbornes { 0 } else { BORNES })
                        .unwrap_or(1)
                        .try_into()
                        .unwrap_or(1),
                ),
                Position::new(
                    false,
                    end.checked_add(if hasbornes { 0 } else { BORNES })
                        .unwrap_or(1)
                        .try_into()
                        .unwrap_or(1),
                ),
                if complement {
                    Strand::Minus
                } else {
                    Strand::Plus
                },
            )
        })
        .collect();
    bornestoprint.sort_unstable();
    if let Err(e) = printpotentialbornes(&bornestoprint, args) {
        eprintln!("Cannot print bornes. Error is {e}");
    }
    loc.sort_unstable();
    Some(loc)
}
pub(crate) fn sendresult(request: &reqwest::blocking::Client, url: &str) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                let total_size = match (
                    e.content_length(),
                    e.headers().get(header::CONTENT_LENGTH).map(|a| {
                        a.to_str()
                            .map_err(|_| "Error reading headers".to_string())
                            .map(|a| a.parse::<u64>())
                            .map_err(|_| "Error reading headers".to_string())
                    }),
                ) {
                    (Some(a), _) => a,
                    (_, Some(Ok(Ok(b)))) => b,
                    _ => 0,
                };
                let stream = e.text().map_err(|e| e.to_string())?;
                let mut downloaded = 0;
                // Read the stream in chunks
                let d = Cursor::new(&stream);
                let reader = BufReader::new(d);
                let mut buffer = vec![0; 1024]; // 1KB buffer
                let mut pb = getprogressbar(total_size, reader).map_err(|e| e.to_string())?;
                while let Ok(n) = pb.read(&mut buffer) {
                    if n == 0 {
                        break;
                    } // End of stream
                    let new = min(downloaded + (n as u64), total_size);
                    downloaded = new;
                }
                pb.progress_bar.finish_with_message("Finished");
                Ok(stream)
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
        Err(e) if e.is_timeout() => Err(format!("Error getting URL (timeout): {e}").to_string()),
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn sendresultcompressed(
    request: &reqwest::blocking::Client,
    url: &str,
) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                let mut buf = Vec::new();
                let (contentlength, headers) = (e.content_length(), e.headers().clone());
                let stream = e.bytes().map_err(|e| e.to_string())?;
                let total_size = match (
                    contentlength,
                    headers.get(header::CONTENT_LENGTH).map(|a| {
                        a.to_str()
                            .map_err(|_| "Error reading headers".to_string())
                            .map(|a| a.parse::<u64>())
                            .map_err(|_| "Error reading headers".to_string())
                    }),
                    stream.len(),
                ) {
                    (Some(a), ..) => a,
                    (_, Some(Ok(Ok(b))), _) => b,
                    (.., c) if let Ok(s) = c.try_into() => s,
                    _ => 0,
                };
                let mut downloaded = 0;
                // Read the stream in chunks
                let d = Cursor::new(stream);
                let reader = BufReader::new(d);
                let mut pb = getprogressbar(total_size, reader).map_err(|e| e.to_string())?;
                let mut buffer = vec![0; 2048]; // 2KB buffer
                while let Ok(n) = pb.read(&mut buffer) {
                    if n == 0 {
                        break;
                    } // End of stream
                    buf.extend(&buffer[0..n]);
                    let new = min(downloaded + (n as u64), total_size);
                    downloaded = new;
                }
                // Decode the buffer (assuming it's a gzip stream)
                let mut decoder = flate2::read::GzDecoder::new(buf.as_slice());
                let mut string = String::new();
                if decoder.read_to_string(&mut string).is_err() {
                    return Err("Cannot decrypt data".to_string());
                }

                if string.trim().is_empty() {
                    return Err("Cannot decrypt motifs, empty data.".to_string());
                }
                pb.progress_bar.finish_with_message("Finished");
                string.shrink_to_fit();
                Ok(string)
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
        Err(e) if e.is_timeout() => Err(format!("Error getting URL (timeout): {e}").to_string()),
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn downloadref(
    allowdownload: bool,
    force: bool,
    releaseversion: &Option<String>,
) -> Option<(Filecrea, String)> {
    let releaseversion = if let Some(r) = releaseversion {
        r.to_string()
    } else {
        println!("Checking reference sequence from IMGT/GENE-DB");
        match sendresult(&REQUESTCLIENT, RELEASELINK.as_str()) {
            Ok(e) => e,
            Err(e) => {
                println!("Release fetched failed because: {e}");
                return None;
            }
        }
    };
    if !allowdownload {
        return Some((
            Filecrea::createtemp(None::<PathBuf>, None).ok()?,
            releaseversion,
        ));
    };
    let tempfile = Filecrea::createfrompath(Path::join(
        &env::temp_dir(),
        format!("refseq{}.fasta", releaseversion),
    ));
    if !tempfile.getpath().is_file() || force {
        println!(
            "Downloading IMGT/GENE-DB release {}, please wait...",
            releaseversion
        );
        match sendresult(&REQUESTCLIENT, VQUESTLINK.as_str()) {
            Ok(e) => {
                match tempfile.setfile().map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    Err(e) | Ok(Err(e)) => {
                        eprintln!("Cannot write refseq in sequence. Error is {e}.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT/GENE-DB data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!(
            "Release {} was already downloaded, retrieving...",
            releaseversion
        );
    };
    Some((tempfile, releaseversion))
}

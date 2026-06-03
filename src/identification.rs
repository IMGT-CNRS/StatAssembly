use bio::io::fasta;
use itertools::Itertools;
use minimap2::Aligner;
use reqwest::StatusCode;

use crate::{
    BORNES, printpotentialbornes,
    r#struct::{
        Args, Blast, Blastcalc, Blastlevel, Blastmatch, Haplotype, Locus, LocusInfos, Name,
        Position, Species, Status, Strand,
    },
    submissions::{
        BORNESLINK, LOCUSSEPARATOR, MOTIFLINK, RELEASELINK, REQUESTCLIENT, VQUESTLINK,
        blastcommand, checkoverlap, fileincache, speciesandorphonfiltering,
    },
};
use std::{
    cmp::{Ordering, max, min},
    collections::{BTreeMap, BTreeSet, HashMap},
    env,
    fs::File,
    io::{self, ErrorKind, Read, Write},
    path::{Path, PathBuf},
    str::FromStr,
    thread,
};

pub(crate) fn locusallposition(
    subject: &Path,
    species: &Species,
    args: &Args,
) -> io::Result<(Vec<LocusInfos>, Vec<Blastmatch>)> {
    let infos = if args.nobornes {
        None
    } else {
        downloadbornes()
    };
    let (bornespath, referencepath) = match (
        infos,
        downloadref(true)
            .map(|(a, b)| speciesandorphonfiltering(&a, None, b, species, true, args.cacheerase)),
    ) {
        (Some(a), Some(b)) => (Some(a), b?),
        (None, Some(b)) if args.nobornes => (None, b?),
        (None, ..) if !args.nobornes => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Bornes from IMGT cannot be downloaded",
            ));
        }
        _ => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    /* let file = NamedTempFile::new_in(&env::temp_dir())?;
    let mut merge = File::create(&file)?;
    let read = fs::read_to_string(&reference.0)?;
    let read2 = fs::read_to_string(&reference.1)?;
    merge.write_all(read.as_bytes())?;
    merge.write_all(read2.as_bytes())?;
    println!("Blasting to get position of loci."); */
    let (blast, bornes) = thread::scope(|s| {
        let blast = s.spawn(|| {
            let b = match blastcommand(
                referencepath.clone(),
                subject.to_path_buf(),
                Blastlevel::Normal,
            ) {
                Ok(b) => b,
                Err(e) => {
                    return Err(io::Error::new(ErrorKind::InvalidData, e));
                }
            };
            Ok(b.into_iter().collect())
        });
        let bornes = if let Some(b) = bornespath {
            Some(s.spawn(|| {
                let aligner = Aligner::builder()
                    .asm20()
                    .with_cigar()
                    .with_index_threads(args.threads)
                    .with_index(subject, None)
                    .ok()?;
                let read = fasta::Reader::from_file(b).ok()?;
                let mut aligns: Vec<Blastmatch> = Vec::new();
                for seq in read.records().filter_map(Result::ok) {
                    if let Ok(a) = aligner
                        .map(
                            seq.seq(),
                            false,
                            false,
                            None,
                            None,
                            Some(format!("NCBI|{}", seq.desc().unwrap_or_default()).as_bytes()),
                        )
                        .map(|mut f| {
                            f.retain(|p| p.is_primary);
                            f
                        })
                        && let Some(first) = a.first()
                        && let (
                            Some(Ok(qseqid)),
                            Some(sseqid),
                            sseq,
                            Ok(sstart),
                            Ok(send),
                            complement,
                            status,
                            Some(Some(identity)),
                        ) = (
                            &first
                                .query_name
                                .as_ref()
                                .map(|a| Name::from_str(a.as_str())),
                            &first.target_name,
                            String::new(),
                            TryInto::<usize>::try_into(first.target_start),
                            TryInto::<usize>::try_into(first.target_end),
                            if first.strand == minimap2::Strand::Reverse {
                                Strand::Minus
                            } else {
                                Strand::Plus
                            },
                            Status::New,
                            (&first)
                                .alignment
                                .as_ref()
                                .map(|p| p.alignment_score.map(|d| d as f32)),
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
                            identity,
                        );
                        aligns.push(matche);
                    }
                }
                if aligns.is_empty() {
                    None
                } else {
                    Some(aligns)
                }
            }))
        } else {
            None
        };
        println!("Waiting for blast and bornes to finish.");
        (blast.join(), bornes.map(|d| d.join()))
    });
    let (mut blast, mut bornes) = match (blast, bornes) {
        (Ok(Ok(a)), Some(Ok(Some(b)))) => (a, Some(b)),
        (Ok(Ok(a)), Some(Ok(None))) | (Ok(Ok(a)), None) => {
            if !args.nobornes {
                eprintln!("No bornes found.");
            }
            (a, None)
        }
        (Ok(Err(a)), _) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Blast error. Error is: {}", a),
            ));
        }
        (Err(a), Some(Err(b))) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!(
                    "Blast and bornes could not be performed. Error are: {:?} and {:?}",
                    a, b
                ),
            ));
        }
        (Err(a), _) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Blast could not be performed. Error is: {:?}", a),
            ));
        }
        (_, Some(Err(b))) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Bornes could not be performed. Error is: {:?}", b),
            ));
        }
    };
    bornes.as_mut().map(|p| {
        p.sort_unstable_by(|f, g| match f.getallelename().cmp(&g.getallelename()) {
            Ordering::Equal => f.getidentity().total_cmp(&g.getidentity()).reverse(),
            e => e,
        });
    });
    //Filter by locus
    retainbestmatch(&mut blast);
    blast.shrink_to_fit();
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
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
                statusvec.insert(k.clone(), elem);
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
    let writer = csv::WriterBuilder::new()
        .delimiter(b',')
        .comment(Some(b'#'))
        .has_headers(true)
        .from_path(Path::join(&env::temp_dir(), "test2.txt"));
    if let Ok(mut r) = writer {
        for elem in data.iter() {
            let _ = r.serialize(elem);
        }
    }
    let range = find_global_best_range(&data, &bornes).ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No locus found after BLAST analysis",
    ));
    if let Some(bornes2) = &mut bornes
        && let Ok(d) = &range
    {
        bornes2.retain(|p| d.iter().any(|k| k.contig == p.getsubject()));
        bornes2.dedup_by(|a, b| {
            if a.getallelename().eq_ignore_ascii_case(&b.getallelename()) {
                true
            } else {
                false
            }
        });
        printpotentialbornes(bornes2, &args)?;
    }
    range.map(|f| (f, data))
}
pub(crate) fn retainbestmatch(blast: &mut Vec<Blast>) {
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
}
pub(crate) fn downloadbornes() -> Option<PathBuf> {
    println!("Downloading bornes sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "bornes.fasta");
    if !fileincache(&tempfile) {
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
    Some(tempfile)
}
pub(crate) fn downloadmotifs() -> Option<PathBuf> {
    println!("Downloading motifs sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "motifs.fasta");
    if !fileincache(&tempfile) {
        println!("Downloading IMGT motifs, please wait...");
        match sendresultcompressed(&REQUESTCLIENT, MOTIFLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
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
    Some(tempfile)
}
pub(crate) fn find_global_best_range(
    blastcheck: &[Blastmatch],
    bornes: &Option<Vec<Blastmatch>>,
) -> Option<Vec<LocusInfos>> {
    let mut groups: HashMap<(Locus, Haplotype), (String, bool, usize, usize)> = HashMap::new();
    let blastcheck = blastcheck.to_vec();
    let mut hash = blastcheck
        .iter()
        .filter_map(|p| p.getlocusname().map(|d| (d, p)))
        .into_group_map_by(|(l, f)| (l.clone(), f.getsubject()));
    hash.iter_mut().for_each(|(_, v)| {
        v.sort_unstable_by(|(_, a), (_, b)| match a.getsubject().cmp(&b.getsubject()) {
            Ordering::Equal => match a.getpos().0.cmp(&b.getpos().0) {
                Ordering::Equal => b.getpos().1.cmp(&b.getpos().1).reverse(),
                ord => ord,
            },
            ord => ord,
        })
    });
    let mut locus: HashMap<(Locus, String, u32), std::collections::BTreeSet<Blastmatch>> =
        HashMap::new();
    for ((loci, sseqname), data) in hash.iter() {
        let mut index = 0;
        for (_, elem) in data {
            let d = match locus.get_mut(&(loci.clone(), sseqname.to_string(), index)) {
                Some(a) => a,
                None => {
                    let mut e = BTreeSet::new();
                    e.insert((*elem).clone());
                    locus.insert((loci.clone(), sseqname.to_string(), index), e);
                    continue;
                }
            };
            if let Some(b) = d.last()
                && std::cmp::max(b.sstart, b.send).abs_diff(std::cmp::min(elem.send, elem.sstart))
                    <= LOCUSSEPARATOR
            {
                d.insert((*elem).clone());
            } else {
                index += 1;
                let mut e = BTreeSet::new();
                e.insert((*elem).clone());
                locus.insert((loci.clone(), sseqname.to_string(), index), e);
            }
        }
    }
    if let Some(borneblast) = &bornes {
        'borne: for borneblast in borneblast {
            for (_, values) in locus.iter_mut() {
                for val in values.iter() {
                    let a = borneblast;
                    let b = val;
                    if std::cmp::max(b.sstart, b.send).abs_diff(std::cmp::min(a.send, a.sstart))
                        <= BORNES
                    {
                        values.insert(borneblast.clone().into());
                        break 'borne;
                    }
                }
            }
        }
    }
    let real =
        locus
            .into_iter()
            .sorted_unstable_by(|((al, ass, _), av), ((bl, bs, _), bv)| match al.cmp(bl) {
                Ordering::Equal => match av.iter().count().cmp(&bv.iter().count()) {
                    Ordering::Equal => ass.cmp(bs),
                    ord2 => ord2.reverse(),
                },
                ord => ord,
            });
    for ((loci, sseq, _), blast) in real {
        if !groups.contains_key(&(loci.clone(), Haplotype::Primary)) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert((loci, Haplotype::Primary), (sseq, complement, min, max));
        } else if !groups.contains_key(&(loci.clone(), Haplotype::Alternate)) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert((loci, Haplotype::Alternate), (sseq, complement, min, max));
        }
    }
    Some(
        groups
            .into_iter()
            .map(|((locus, hap), (sseq, complement, start, end))| {
                LocusInfos::new(
                    locus,
                    hap,
                    sseq,
                    Position::new(
                        false,
                        start.saturating_sub(BORNES).try_into().unwrap_or_default(),
                    ),
                    Position::new(
                        false,
                        end.saturating_add(BORNES).try_into().unwrap_or_default(),
                    ),
                    if complement {
                        Strand::Minus
                    } else {
                        Strand::Plus
                    },
                )
            })
            .collect(),
    )
}
pub(crate) fn sendresult(request: &reqwest::blocking::Client, url: &str) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                Ok(e.text().unwrap_or("Error getting data".to_string()))
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
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
                let mut string = String::new();
                if let Ok(bytes) = e.bytes()
                    && let mut decode = flate2::read::GzDecoder::new(&bytes[..])
                    && let Err(e) = decode.read_to_string(&mut string)
                {
                    return Err(format!("Cannot decrypt motifs: {e}"));
                }
                if string.trim().is_empty() {
                    return Err(format!("Cannot decrypt motifs, invalid body."));
                }
                string.shrink_to_fit();
                Ok(string)
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn downloadref(allowdownload: bool) -> Option<(PathBuf, String)> {
    println!("Checking reference sequence from IMGT/GENE-DB");
    let releaseversion = match sendresult(&REQUESTCLIENT, RELEASELINK.as_str()) {
        Ok(e) => e,
        Err(e) => {
            println!("Release fetched failed because: {e}");
            return None;
        }
    };
    if !allowdownload {
        return Some((PathBuf::new(), releaseversion));
    };
    let tempfile = Path::join(&env::temp_dir(), format!("refseq{}.fasta", releaseversion));
    if !tempfile.is_file() {
        println!(
            "Downloading IMGT/GENE-DB release {}, please wait...",
            releaseversion
        );
        match sendresult(&REQUESTCLIENT, VQUESTLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    _ => {
                        eprintln!("Cannot write refseq in sequence.");
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

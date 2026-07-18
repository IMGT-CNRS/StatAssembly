#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
use crate::identification::{downloadmotifs, downloadref, retainbestmatch};
use crate::r#struct::{
    Args, Blast, Blastcalc, Blastlevel, Blastmatch, Filecrea, GeneInfos, Locus, LocusInfos, Name,
    Newfasta, Position, ProgressReader, Seqresult, Species, Strand, safestring,
};
use crate::{EMAIL, NAME, PHYLUMLIMIT, TIMEOUT_IN_MN, VERSION, getassemblyreader, getreaderoffile};
use bio::io::fasta;
use extended_htslib::bam::index::Type;
use extended_htslib::bam::{self, Read};
use flate2::GzBuilder;
use flate2::write::GzDecoder;
use flate2::{Compression, write::GzEncoder};
use indicatif::ProgressBar;
use itertools::Itertools;
use lazy_static::lazy_static;
use quick_xml::events::Event;
use reqwest::blocking::multipart;
use reqwest::{StatusCode, tls};
use std::cmp::{Ordering, max, min};
use std::collections::HashMap;
use std::fmt::Debug;
use std::io::ErrorKind::InvalidData;
use std::io::{BufWriter, IsTerminal, Seek, Write};
use std::path::PathBuf;
use std::str::FromStr;
use std::thread::{self, available_parallelism, sleep};
use std::{
    env::{self},
    error::Error,
    fs::{self, File},
    io::{self, BufRead, BufReader, ErrorKind, Read as _},
    ops::{Not, RangeInclusive},
    path::Path,
    process::{Command, Stdio},
    time::Duration,
};
use strum::IntoEnumIterator;
use tempfile::TempDir;
lazy_static! {
    pub static ref REQUESTCLIENT: reqwest::blocking::Client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::new(15, 0))
        .timeout(Duration::new(600, 0))
        .tcp_user_timeout(Duration::new(400, 0))
        .user_agent(format!("IMGT/StatAssembly version {}", VERSION))
        .referer(false)
        .tls_version_min(tls::Version::TLS_1_2)
        .https_only(true)
        .http3_max_idle_timeout(Duration::new(10,0))
        .http3_prior_knowledge()
        .build()
        .unwrap_or_default();
    pub static ref WEBSERVER: String = obfstr::obfstring!("https://imgt.org");
    pub static ref RELEASELINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstr!("/download/GENE-DB/RELEASE")
    );
    pub static ref LIMITDATE: String = obfstr::obfstring!("2026-08-31T23:59:59Z");
    pub static ref MOTIFLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!("/submissions/newmotif_fusionne.fasta.gz")
    );
    pub static ref INVALIDCOVERAGE: String = obfstr::obfstring!("Coverage is out of the window");
    pub static ref NOTENOUGHMATCHREADS: String = obfstr::obfstring!("Not enough matching reads");
    pub static ref SOFTCLIPTOOMUCH: String = obfstr::obfstring!("Too much softclips");
    pub static ref SUSPICIOUSPOSITIONALERT: String = obfstr::obfstring!("The position");
    /* pub static ref MOTIFLINK: String =
        obfstr::obfstring!("http://localhost:8910/submissions/newmotif_fusionne.fasta.gz");
    pub static ref BORNESLINK: String =
        obfstr::obfstring!("http://localhost:8910/submissions/bornes_mammals.fasta.gz");  */
        pub static ref BORNESLINK: String = format!(
            "{}{}",
            WEBSERVER.as_str(),
            obfstr::obfstring!("/submissions/bornes_mammals.fasta.gz")
        );
    pub static ref VQUESTLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!(
            "/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"
        )
    );
    pub static ref GITHUBVERSION: String = obfstr::obfstring!(
        "https://src.koda.cnrs.fr/api/v4/projects/imgt-igh%2Fstatassembly/repository/tags");
    pub static ref SUBMISSIONLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!("/submissions/")
    );
}

pub(crate) fn readfromterminal(
    yes: &char,
    no: &char,
    nothing: Option<&char>,
    force: bool,
) -> Option<bool> {
    let io = io::stdin();
    if !io.is_terminal() {
        return Some(false);
    }
    let mut data = String::new();
    loop {
        data.clear();
        if io.read_line(&mut data).is_err() {
            return Some(false);
        }
        if data.trim().is_empty() && !force {
            return Some(false);
        }
        if data.to_lowercase().trim().chars().all(|p| &p == yes) {
            return Some(true);
        } else if let Some(c) = nothing
            && data.to_lowercase().trim().chars().all(|p| &p == c)
        {
            return None;
        } else if data.to_lowercase().trim().chars().all(|p| &p == no) || !force {
            return Some(false);
        }
    }
}
#[deprecated]
#[allow(unused)]
pub(crate) fn getnamefromblast<T>(text: T) -> Option<String>
where
    T: AsRef<str>,
{
    text.as_ref()
        .to_ascii_lowercase()
        .split('|')
        .nth(2)
        .map(|f| f.to_string())
}
pub(crate) fn getallelefromblast(text: &Name) -> &str {
    &text.gene
}
pub(crate) fn getchromosomefromblast<T>(text: T) -> Option<String>
where
    T: AsRef<str>,
{
    text.as_ref()
        .to_ascii_uppercase()
        .split('|')
        .nth(3)
        .map(|f| f.to_string())
}
pub(crate) fn getpositionfromblast(text: &str) -> Option<(Position, Position, Strand)> {
    text.to_ascii_uppercase().split("|").nth(4).and_then(|p| {
        let p = p.replace("/", "..");
        let mut split = p.splitn(3, "..");
        match (
            split.next().map(|p| p.parse::<i64>()),
            split.next().map(|p| p.parse::<i64>()),
            split.next().map(|a| a.trim()),
        ) {
            (Some(Ok(a)), Some(Ok(b)), c) => {
                let complement = c.is_some_and(|f| !f.is_empty());
                Some((
                    Position::newfromoposition(a),
                    Position::newfromoposition(b),
                    if complement {
                        Strand::Minus
                    } else {
                        Strand::Plus
                    },
                ))
            }
            _ => None,
        }
    })
}
pub(crate) fn removeallspaces(text: &mut String) {
    *text = text.chars().filter(|p| !p.is_ascii_whitespace()).collect();
}
#[must_use]
pub(crate) fn fastafilter(text: &str, find: &str, present: bool, species: bool) -> String {
    let buffer = BufReader::new(text.as_bytes());
    let parser = fasta::Reader::new(buffer);
    let mut lines: Vec<String> = Vec::with_capacity(text.lines().count());
    for record in parser
        .records()
        .filter_map(Result::ok)
        .filter(|p| p.check().is_ok())
    {
        let val = format!(
            "{}{}",
            record.id(),
            record.desc().map_or(String::new(), |a| format!(" {a}"))
        );
        let keep = match Name::from_str(&val) {
            Ok(mut result) => {
                let mut find = find.to_string();
                removeallspaces(&mut result.species);
                removeallspaces(&mut result.gene);
                removeallspaces(&mut find);
                match (present, species) {
                    (true, true) => result.species.eq_ignore_ascii_case(&find),
                    (false, true) => !result.species.eq_ignore_ascii_case(&find),
                    (false, false) => !result.gene.contains(&find),
                    (true, false) => result.gene.contains(&find),
                }
            }
            Err(_e) => {
                //eprintln!("Error parsing line: {}. Error is {e}", val);
                false
            }
        };
        if keep {
            lines.push(record.to_string());
        }
    }
    lines
        .iter()
        .fold(String::new(), |mut acc, f| {
            acc.push_str(format!("\n{}", f.trim()).as_str());
            acc
        })
        .trim()
        .to_string()
}
pub(crate) fn checkoverlap(a: &RangeInclusive<usize>, b: &RangeInclusive<usize>) -> bool {
    max(a.start(), b.start()) <= min(a.end(), b.end())
}
pub(crate) fn fileincache<T>(tempfile: T) -> bool
where
    T: AsRef<Path>,
{
    if !tempfile.as_ref().is_file() {
        return false;
    }
    tempfile
        .as_ref()
        .metadata()
        .and_then(|p| {
            p.accessed()
                .map(|p| p.elapsed().map(|p| p < Duration::new(3600 * 24 * 7, 0)))
        })
        .is_ok_and(|f| f.unwrap_or_default())
}
pub(crate) fn checkifblastpresent() -> bool {
    println!("Detect if BLAST is operating.");
    let command = Command::new("blastn")
        .arg("-version")
        .stdout(Stdio::null())
        .stdin(Stdio::null())
        .status()
        .is_ok_and(|f| f.success());
    if !command {
        eprintln!("BLAST was not found or has failed. Check if present in PATH.");
        false
    } else {
        println!("BLAST is working. Continuing.");
        true
    }
}
pub(crate) fn positionfiltering<T>(locus: &[LocusInfos], blast: &mut Vec<T>)
where
    T: Blastcalc,
{
    blast.retain(|p| {
        locus.iter().any(|f| {
            let range: RangeInclusive<usize> = f.start.getobasedpos().try_into().unwrap_or_default()
                ..=f.end.getobasedpos().try_into().unwrap_or_default();
            f.contig.as_str() == p.getsubject()
                && checkoverlap(&(p.getpos().0..=p.getpos().1), &range)
        })
    })
}
pub(crate) fn locusfiltering<T>(locus: &Locus, blast: &mut Vec<T>)
where
    T: Blastcalc,
{
    //TRD is inside TRA locus
    if locus == &Locus::TRA {
        blast.retain(|p| {
            p.getlocusname().as_ref() == Some(locus) || p.getallelename().starts_with("TRD")
        });
    } else {
        blast.retain(|p| p.getlocusname().as_ref() == Some(locus));
    }
}
pub(crate) fn checklocus(p: &fasta::Record, locus: &Locus) -> bool {
    match locus {
        Locus::IGH => p.id().contains("IGH") || p.desc().is_some_and(|f| f.contains("IGH")),
        Locus::IGK => p.id().contains("IGK") || p.desc().is_some_and(|f| f.contains("IGK")),
        Locus::IGL => p.id().contains("IGL") || p.desc().is_some_and(|f| f.contains("IGL")),
        Locus::TRA => {
            p.id().contains("TRA")
                || p.desc().is_some_and(|f| f.contains("TRA"))
                || p.id().contains("TRD")
                || p.desc().is_some_and(|f| f.contains("TRD"))
        }
        Locus::TRB => p.id().contains("TRB") || p.desc().is_some_and(|f| f.contains("TRB")),
        Locus::TRG => p.id().contains("TRG") || p.desc().is_some_and(|f| f.contains("TRG")),
    }
}
pub(crate) fn speciesandorphonfiltering(
    tempfile: &mut Filecrea,
    locus: Option<&Locus>,
    releaseversion: String,
    species: &Species,
    orphonfilter: bool,
    force: bool,
) -> io::Result<Filecrea> {
    let outfile = if force {
        Filecrea::createtemp(
            None,
            Some(format!(
                "refseq{}-{}{}.fasta",
                releaseversion.replace(" ", "_"),
                species.getname().replace(" ", "-"),
                locus.map_or("".to_string(), |l| format!("-{}", l))
            )),
        )?
    } else {
        Filecrea::createfrompath(Path::join(
            &env::temp_dir(),
            format!(
                "refseq{}-{}{}.fasta",
                releaseversion.replace(" ", "_"),
                species.getname().replace(" ", "-"),
                locus.map_or("".to_string(), |l| format!("-{}", l))
            ),
        ))
    };
    if !outfile.istemp() && outfile.getpath().try_exists().unwrap_or(false) && !force {
        println!("Filtering already done, retrieving...");
        return Ok(outfile);
    }
    println!(
        "Filtering based on {} {}.",
        species.getrank(),
        species.getname()
    );
    let mut file = String::new();
    tempfile.getfile()?.read_to_string(&mut file)?;
    let info = fastafilter(&file, species.getname(), true, true);
    let (info, file) = if orphonfilter {
        println!("Orphon filtering");
        (
            fastafilter(&info, "/OR", false, false).replace(" ", "_"),
            fastafilter(&file, "/OR", false, false).replace(" ", "_"),
        )
    } else {
        (info, file)
    };
    let finale = {
        //Check if there is match for all loci or the loci
        let tempnew = io::Cursor::new(&info);
        let buf = BufReader::new(tempnew);
        let read = fasta::Reader::from_bufread(buf);
        let allmatch = if let Some(l) = locus {
            read.records()
                .filter_map(Result::ok)
                .any(|p| checklocus(&p, l))
        } else {
            let mut count = Locus::iter().map(|p| (p, 0)).collect_vec();
            read.records().filter_map(Result::ok).for_each(|p| {
                for (locus, count) in count.iter_mut() {
                    if checklocus(&p, locus) {
                        *count += 1;
                        break;
                    }
                }
            });
            count.iter().all(|(_, f)| *f > 0)
        };
        match (info.is_empty(), allmatch) {
            (true, _) => {
                println!(
                    "No match with IMGT/GENE-DB, the {} {} might be new.",
                    species.getrank(),
                    species.getname()
                );
                file
            }
            (false, false) => {
                if let Some(a) = locus {
                    println!(
                        "No match with IMGT/GENE-DB for {} and {} loci.",
                        species.getname(),
                        a
                    );
                } else {
                    println!(
                        "No match with IMGT/GENE-DB with {} and some loci.",
                        species.getname()
                    );
                }
                file
            }
            _ => info,
        }
    };
    let info = if let Some(l) = locus {
        let tempnew = io::Cursor::new(&finale);
        let buf = BufReader::new(tempnew);
        let read = fasta::Reader::from_bufread(buf);
        let records = read
            .records()
            .filter_map(Result::ok)
            .filter(|a| a.check().is_ok())
            .filter(|p| checklocus(p, l));
        let seq = records.fold(String::new(), |mut acc, f| {
            acc.push_str(&format!(
                ">{} {}\n{}\n",
                f.id(),
                f.desc().unwrap_or_default(),
                String::from_utf8_lossy(f.seq())
            ));
            acc
        });
        seq.trim().to_string()
    } else {
        finale
    };
    if info.is_empty() {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "Empty data after filtering",
        ));
    }
    /* let newdata = if newdata.is_empty() {
        println!("New species!!");
        let bufreader = io::Cursor::new(newdata);
        let reader = fasta::Reader::new(bufreader);
        let fastaread = readfastareader(reader);
        return match fastaread {
            Ok(b) => {
                let mut new = Vec::new();
                b.iter().for_each(|f| {
                    if let Ok(v) = Newfasta::new(f) {
                        new.push(v);
                    } else {
                        eprintln!("Line {} is skipped as an error occured", f.name)
                    }
                });
                (Some(new), None)
            }
            Err(e) => {
                eprintln!("{e}");
                (None, None)
            }
        };
        file
    } else {
        newdata
    }; */
    if tempfile.getpath() != outfile.getpath() {
        std::fs::write(outfile.getpath(), &info)?;
    }
    Ok(outfile)
}
pub(crate) fn readfastareader<T>(fasta: fasta::Reader<T>) -> io::Result<Vec<Newfasta>>
where
    T: std::io::BufRead,
{
    let mut seqs = Vec::new();
    for record in fasta.records().filter_map(|p| p.ok()) {
        if let Err(v) = record.check() {
            eprintln!("Seq {} is invalid, skipped. Error is {v}", record.id());
            continue;
        }
        let text = format!(
            "{} {}\n{}",
            record.id(),
            record.desc().unwrap_or_default(),
            String::from_utf8_lossy(record.seq())
        );
        match Newfasta::from_str(&text) {
            Ok(a) => seqs.push(a),
            Err(e) => {
                eprintln!("Seq {} is invalid, skipped. Error is {}", record.id(), e);
                continue;
            }
        };
    }
    if seqs.is_empty() {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "No fasta record valid found",
        ));
    }
    Ok(seqs)
}
#[allow(unused)]
pub(crate) fn readfastafile(seq: &Path) -> io::Result<Vec<Newfasta>> {
    let file = File::open(seq)?;
    let fasta = fasta::Reader::new(file);
    readfastareader(fasta)
}
pub(crate) fn statusblastmotifs(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.sseqid.cmp(&b.sseqid) {
        std::cmp::Ordering::Equal => match a.sstart.cmp(&b.sstart) {
            Ordering::Equal => a.send.cmp(&b.send).reverse(),
            ord => ord,
        },
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && g.getallelename() == f.getallelename()
                    && checkoverlap(&(g.sstart..=g.send), &(f.sstart..=f.send))
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn statusblastvs(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.getsubject().cmp(b.getsubject()) {
        std::cmp::Ordering::Equal => a.getquery().gene.cmp(&b.getquery().gene),
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && g.sseqid == f.sseqid
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
#[allow(unused)]
pub(crate) fn statusblast(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.sstart.cmp(&b.sstart) {
        std::cmp::Ordering::Equal => a.send.cmp(&b.send),
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && checkoverlap(&g.getposrange(), &f.getposrange())
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn matchmotif(
    subject: &Path,
    species: &Species,
    args: &Args,
    locus: Option<&Locus>,
) -> io::Result<Vec<Blastmatch>> {
    let reference = match downloadmotifs(args).map(|mut a| {
        speciesandorphonfiltering(&mut a, locus, "motifs".to_string(), species, false, false)
    }) {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Motifs from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast: Vec<Blast> = match blastcommand(
        reference,
        subject.to_path_buf().into(),
        Blastlevel::default(),
    ) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    if let Some(a) = locus {
        locusfiltering(a, &mut blast);
    }
    blast.iter_mut().for_each(|_p| {
        //Already performed
        //if p.sstart > p.send {
        //    (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
        //}
    });
    statusblastmotifs(&mut blast);
    //let _ = fs::remove_file(name);
    Ok(blast.into_iter().map(|f| f.into()).collect())
}
pub(crate) fn genesblast(
    subject: &[GeneInfos],
    args: &Args,
    releaseversion: &Option<String>,
    species: &Species,
    locus: &Locus,
) -> io::Result<(Vec<Blastmatch>, Option<String>)> {
    let name = Filecrea::createtemp(None, Some("genes_blast.txt"))?;
    let mut fastawriter = fasta::Writer::to_file(name.getpath())?;
    let mut reader = getassemblyreader(args)?;
    let count = subject
        .iter()
        .map(|f| {
            let elem = f
                .extractsequence(&mut reader)
                .map_err(|p| io::Error::new(ErrorKind::InvalidInput, p));
            let bool = elem
                .as_ref()
                .is_ok_and(|p| f.addtosequence(p, &mut fastawriter, species).is_err());

            if bool {
                Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    "Cannot print sequence",
                ))
            } else {
                elem
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    fastawriter.flush()?;
    if count.is_empty() {
        eprintln!("No genes found for this locus, skipped.");
        return Ok((Vec::new(), None));
    }
    let (reference, version) =
        match downloadref(true, args.cacheerase, releaseversion).map(|(mut a, b)| {
            (
                speciesandorphonfiltering(
                    &mut a,
                    Some(locus),
                    b.clone(),
                    species,
                    false,
                    args.cacheerase,
                ),
                b,
            )
        }) {
            Some((a, b)) => (a?, b),
            None => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Reference from IMGT cannot be downloaded",
                ));
            }
        };
    let mut blast = match blastcommand(reference, name, Blastlevel::default()) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    locusfiltering(locus, &mut blast);
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
        }
    });
    statusblastvs(&mut blast);
    Ok((blast.into_iter().map(|f| f.into()).collect(), Some(version)))
}
/*
pub(crate) fn locuspos(
    search: &Locus,
    hap: &Haplotype,
    locus: &[LocusInfos],
    blast: &[Blastmatch],
) -> Option<(LocusInfos, Vec<Blastmatch>)> {
    let opt = locus
        .iter()
        .find(|p| &p.getlocus()== search && &p.haplotype == hap)?;
    let fil = |a: &Blastmatch| {
        let opt = opt.clone();
        a.getlocusname() == Some(opt.locus)
            && a.sseqid == opt.contig
            && (opt.start.getobasedpos()..=opt.end.getobasedpos())
                .contains(&a.sstart.try_into().unwrap_or_default())
    };
    if !blast.iter().any(fil) {
        return None;
    }
    let mut bl = blast.to_vec();
    bl.retain(fil);
    Some((opt.clone(), bl))
}
     */
pub(crate) fn filter_new_alleles<'a, T>(data: &'a [T], motifs: &[T]) -> impl Iterator<Item = &'a T>
where
    T: Blastcalc + Debug,
{
    data.iter().filter(|p| p.onlynewalleles()).filter(move |f| {
        motifs.iter().any(|p| {
            checkoverlap(&p.getposrange(), &f.getposrange())
                && p.getposrange() != f.getposrange()
                && p.getposrange().count() > f.getposrange().count()
        })
    })
}
//if there is space, BLAST would cut the header and parsing would fail, we therefore regenerate a file for reference without this issue
fn checknospaceinheaderandvalidseq(
    reference: &mut Filecrea,
    _subject: &Filecrea,
) -> io::Result<()> {
    let refreader = fasta::Reader::from_file(reference.getpath())
        .map_err(|b| io::Error::new(ErrorKind::NotFound, b.to_string()))?;
    let haspaces = refreader.records().filter_map(Result::ok).any(|f| {
        f.desc().is_some() || f.id().contains(',') || f.desc().is_some_and(|f| f.contains(','))
    }); //One sequence has a space because it has a description
    if haspaces {
        let refnew = Filecrea::createtemp(None::<&Path>, None::<&Path>)?;
        let mut writer = fasta::Writer::to_file(refnew.getpath())?;
        let refreader = fasta::Reader::from_file(reference.getpath())
            .map_err(|b| io::Error::new(ErrorKind::NotFound, b.to_string()))?;
        for seq in refreader.records().filter_map(Result::ok) {
            if let Err(e) = seq.check() {
                eprintln!("Sequence {} has invalid sequence. Error is {}", seq.id(), e);
                continue;
            }
            writer.write(
                &format!(
                    "{}{}",
                    safestring(seq.id()),
                    seq.desc()
                        .map_or(String::new(), |d| format!("_{}", safestring(d)))
                ),
                None,
                seq.seq(),
            )?;
        }
        *reference = refnew;
    }
    Ok(())
}
pub(crate) fn blastcommand<T>(
    reference: T,
    subject: T,
    blastlevel: Blastlevel,
) -> io::Result<Vec<Blast>>
where
    T: Into<Filecrea>,
{
    let (mut reference, subject) = (reference.into(), subject.into());
    checknospaceinheaderandvalidseq(&mut reference, &subject)?;
    let (refpath, subpath) = (reference.getpath(), subject.getpath());
    if !refpath.try_exists().unwrap_or(false) {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Reference file was not found",
        ));
    }
    if !subpath.try_exists().unwrap_or(false) {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Subject file was not found",
        ));
    }
    let reference = &format!("{}", refpath.display());
    let subject = &format!("{}", subpath.display());
    let output = Filecrea::createtemp(None::<&Path>, None::<&Path>)?;
    let output = &format!("{}", output.getpath().display());
    let mut command = Command::new("blastn")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .args([
            "-query",
            reference,
            "-subject",
            subject,
            "-perc_identity",
            "75",
            "-out",
            output,
            "-word_size",
            blastlevel.as_word_size().to_string().as_str(),
            "-max_target_seqs",
            blastlevel.as_max_matches().to_string().as_str(),
            "-max_hsps",
            blastlevel.as_max_matches().to_string().as_str(),
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send qlen length pident gaps sseq",
        ])
        .spawn()?;
    println!("Launching {} against {}", reference, subject);
    println!(
        "BLAST has been launched with id {}. Please wait.",
        command.id()
    );
    let threadoutput = thread::spawn(move || {
        let mut time = 0;
        while time < TIMEOUT_IN_MN.saturating_mul(60) {
            if !command.try_wait().is_ok_and(|f| f.is_none()) {
                return command.wait_with_output().map(Some);
            }
            sleep(Duration::new(1, 0));
            time += 1;
        }
        command.kill().map(|_| None)
    })
    .join();
    let threadoutput = threadoutput.map_err(|_| io::Error::from(ErrorKind::BrokenPipe))?;
    let mut result: Vec<Blast> = Vec::new();
    {
        match threadoutput {
            Err(e) => {
                return Err(io::Error::new(
                    ErrorKind::BrokenPipe,
                    format!("BLAST has failed. Error is {}", e),
                ));
            }
            Ok(None) => {
                return Err(io::Error::new(
                    ErrorKind::ConnectionAborted,
                    "BLAST has timeout. Please retry later.".to_string(),
                ));
            }
            Ok(Some(b)) if !b.status.success() => {
                return Err(io::Error::new(
                    ErrorKind::BrokenPipe,
                    format!(
                        "BLAST has failed. Errorcode is {}. Error is {}",
                        b.status,
                        String::from_utf8_lossy(&b.stderr)
                    ),
                ));
            }
            Ok(Some(_)) => {
                println!("BLAST has been done. Parsing.");
                let file = File::open(output)?;
                let reader = BufReader::new(file);
                let mut csv = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .comment(Some(b'#'))
                    .has_headers(false)
                    .from_reader(reader);
                for record in csv.deserialize::<Blast>() {
                    if let Ok(mut r) = record {
                        if r.sstart > r.send {
                            (r.sstart, r.send, r.complement) = (r.send, r.sstart, Strand::Minus)
                        }
                        result.push(r);
                    } else if let Err(r) = record {
                        eprintln!("Error parsing line. Error is: {r}.");
                        continue;
                        //return Err(io::Error::new(ErrorKind::InvalidData, r.to_string()));
                    }
                }
            }
        };
    };
    //let _ = fs::remove_file(output);
    Ok(result)
}
/// Returns rank, name and id
pub(crate) fn getspeciesfromncbi<T>(
    client: &reqwest::blocking::Client,
    species: &T,
) -> Result<(String, String, usize), Box<dyn Error>>
where
    T: AsRef<str>,
{
    let species = species.as_ref();
    let id = if let Ok(val) = species.parse::<usize>() {
        val
    } else {
        let response = client
            .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?")
            .query(&[
                ("db", "taxonomy"),
                ("term", species),
                ("tool", NAME.as_str()),
                ("email", EMAIL),
                ("format", "xml"),
            ])
            .send()?;
        if response.status().is_server_error() {
            return Err(Box::new(io::Error::new(
                ErrorKind::HostUnreachable,
                "Unsuccessful NCBI response (server error)",
            )));
        } else if !response.status().is_success() {
            return Err(Box::new(io::Error::new(
                ErrorKind::NotFound,
                format!("Unsuccessful NCBI response. Error is {}", response.status()),
            )));
        }
        let responsebody = match response.text() {
            Ok(a) => a,
            Err(_) => {
                return Err(Box::new(io::Error::new(
                    ErrorKind::HostUnreachable,
                    "Badly formatted NCBI output.",
                )));
            }
        };
        let mut xml = quick_xml::Reader::from_str(&responsebody);
        loop {
            match xml.read_event() {
                Ok(Event::Eof) => {
                    return Err(Box::new(io::Error::new(
                        ErrorKind::InvalidInput,
                        "Invalid species from NCBI.",
                    )));
                }
                Ok(Event::Start(c)) if c.name().as_ref().eq_ignore_ascii_case(b"Id") => {
                    match xml
                        .read_text(c.name())
                        .map(|a| String::from_utf8_lossy(&a).parse::<usize>())
                    {
                        Ok(Ok(a)) => break a,
                        _ => continue,
                    }
                }
                _ => continue,
            }
        }
    };
    sleep(Duration::from_millis(500)); //sleep for NCBI
    let response = client
        .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?")
        .query(&[
            ("db", "taxonomy"),
            ("id", &format!("{id}")),
            ("tool", NAME.as_str()),
            ("email", EMAIL),
            ("format", "xml"),
        ])
        .send()?;
    if response.status().is_server_error() {
        return Err(Box::new(io::Error::new(
            ErrorKind::HostUnreachable,
            "Unsuccessful NCBI response (server error)",
        )));
    } else if !response.status().is_success() {
        return Err(Box::new(io::Error::new(
            ErrorKind::NotFound,
            format!(
                "Unsuccessful NCBI response. Code is {}, text is {}",
                response.status(),
                response.text().unwrap_or_default()
            ),
        )));
    }
    let responsebody = match response.text() {
        Ok(a) => a,
        Err(_) => {
            return Err(Box::new(io::Error::new(
                ErrorKind::HostUnreachable,
                "Badly formatted NCBI output.",
            )));
        }
    };
    let mut xml = quick_xml::Reader::from_str(&responsebody);
    let mut lineage = false;
    let mut rank = String::new();
    let mut scientificname = String::new();
    let mut id = 0;
    let goodtaxononomy = loop {
        match xml.read_event() {
            Ok(Event::Eof) => break false,
            Ok(Event::Start(b)) if b.name().as_ref() == b"LineageEx" => {
                lineage = true;
                continue;
            }
            Ok(Event::Start(a)) if a.name().as_ref() == b"TaxId" && lineage => {
                match xml
                    .read_text(a.name())
                    .map(|a| String::from_utf8_lossy(&a).parse::<usize>())
                {
                    Ok(Ok(a)) if a == PHYLUMLIMIT => break true,
                    _ => continue,
                }
            }
            Ok(Event::Start(c))
                if c.name().as_ref().eq_ignore_ascii_case(b"scientificname") && !lineage =>
            {
                match xml.read_text(c.name()) {
                    Ok(a) => scientificname = String::from_utf8_lossy(&a).to_string(),
                    _ => continue,
                }
            }
            Ok(Event::Start(c)) if c.name().as_ref().eq_ignore_ascii_case(b"Rank") && !lineage => {
                match xml.read_text(c.name()) {
                    Ok(a) => rank = String::from_utf8_lossy(&a).to_string(),
                    _ => continue,
                }
            }
            Ok(Event::Start(c)) if c.name().as_ref().eq_ignore_ascii_case(b"TaxId") && !lineage => {
                match xml
                    .read_text(c.name())
                    .map(|a| String::from_utf8_lossy(&a).parse::<usize>())
                {
                    Ok(Ok(a)) => id = a,
                    _ => continue,
                }
            }
            _ => continue,
        }
    };
    match (rank, id, goodtaxononomy, scientificname) {
        (_, 0, ..) => Err(Box::new(io::Error::new(
            ErrorKind::InvalidInput,
            "Species was not found on NCBI.",
        ))),
        (b, .., a) if a.is_empty() || b.is_empty() => Err(Box::new(io::Error::new(
            ErrorKind::InvalidInput,
            "Species was not found on NCBI.",
        ))),
        (rank, _, false, _)
            if rank.eq_ignore_ascii_case("species") || rank.eq_ignore_ascii_case("subspecies") =>
        {
            Err(Box::new(io::Error::new(
                ErrorKind::Unsupported,
                "Species is not a jawed vertebrate.",
            )))
        }
        (rank, id, true, name)
            if !name.is_empty() && rank.eq_ignore_ascii_case("species")
                || rank.eq_ignore_ascii_case("subspecies") =>
        {
            Ok((rank, name, id))
        }
        _ => Err(Box::new(io::Error::new(
            ErrorKind::InvalidInput,
            "The term used is not a species or an subspecies.",
        ))),
    }
}
/*
pub(crate) fn launchblast(
    species: &str,
    locus: &Locus,
    subject: &Path,
    args: &Args,
) -> Result<Vec<Newfasta>, ()> {
    let (realspecies, _) = match getspeciesfromncbi(&REQUESTCLIENT, &species) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    println!("The species is {}.", realspecies.);
    let (path, releaseversion) = match downloadref(true).map(|(a, b)| {
        (
            speciesandorphonfiltering(
                &a,
                Some(locus),
                b.clone(),
                species.as_ref(),
                false,
                args.cacheerase,
            ),
            b,
        )
    }) {
        None => return Err(()),
        Some((Err(a), ..)) => return Err(()),
        Some((Ok(a), b)) => (a, b),
    };
    let filtering = match speciesandorphonfiltering(
        &path,
        Some(locus),
        releaseversion,
        &realspecies,
        false,
        args.cacheerase,
    ) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    let result: Vec<Newfasta> =
        match blastcommand(filtering.as_path(), subject, Blastlevel::default()) {
            Err(e) => {
                eprintln!("{e}");
                return Err(());
            }
            Ok(mut c) => {
                statusblast(&mut c);
                let status: Vec<Newfasta> = c
                    .iter()
                    .map(|c| {
                        let b: &dyn Blastcalc = c;
                        Newfasta::from(b)
                    })
                    .collect();
                let sequence = status.iter().fold(String::new(), |mut acc, f| {
                    let f: &dyn Seqresult = f;
                    acc.push_str(&format!("\n{}", f));
                    acc
                });
                if let Err(e) = fs::write("sequence_finished.fasta", sequence) {
                    eprintln!("An error has occured while priting sequence: {e}.");
                    return Err(());
                }
                status
            }
        };
    if result.is_empty() {
        eprintln!(
            "BLAST result is empty. No new alleles matching IMGT criteria has been identified."
        );
        return Err(());
    }
    Ok(result)
} */
pub(crate) fn submit(
    args: &Args,
    locus: &[crate::LocusInfos],
    c: &[Blastmatch],
    realspecies: &Species,
) -> Result<(), String> {
    //let result: Vec<Newfasta> = c.into_iter().map(Newfasta::newfromblastowner).collect();
    let dir = env::temp_dir();
    /* if dir.is_dir() {
        eprintln!("Archive directory exists, going to be deleted.");
        if let Err(e) = fs::remove_dir_all(&dir) {
            let dir = dir.display();
            return Err(format!("Cannot remove the directory {dir}, error is {e}."));
        }
    } */
    let dirtemp = match TempDir::new_in(dir) {
        Err(e) => return Err(format!("Cannot create archive directory, error is {e}")),
        Ok(p) => p,
    };
    let dir = dirtemp.path();
    let _ = match args.assembly.as_ref().map(|_| getassemblyreader(args)) {
        None => return Err("No assembly provided.".to_string()),
        Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
        Some(Ok(b)) => b,
    };
    /* let locuspos =
        File::create(dir.join("newloc.csv")).map_err(|f| format!("Locus csv, error is {f}"))?;
    locuspos
        .lock()
        .map_err(|f| format!("Error acquiring lock. Error is {f}"))?;
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(locuspos);
    for loci in locus.iter() {
        csv.serialize(loci)
            .map_err(|p| format!("Error serializing locus position: {p}"))?;
    }
    csv.flush()
        .map_err(|p| format!("Error serializing locus position: {p}"))?; */
    let lightbam = dir.join("outlight.bam");
    generatelightbam(args, &lightbam, None, locus)?;
    let sequencefile = generatesequence(args, dir, false, locus)?;
    let sequencefiledecr = decompressseq(&sequencefile).map_err(|d| d.to_string())?;
    let mut motifs = matchmotif(sequencefiledecr.getpath(), realspecies, args, None)
        .map_err(|f| format!("Error matching motifs: {f}").to_string())?;
    motifs.iter_mut().for_each(|p| {
        if let Some(find) = locus
            .iter()
            .find(|k| format!("{}:{}", k.getlocus(), k.contig) == p.sseqid)
            && let Some((newstart, newend, newcomplement)) = find.positioninlocus(
                &Position::newfromoposition(p.sstart.try_into().unwrap_or_default()),
                &Position::newfromoposition(p.send.try_into().unwrap_or_default()),
                &p.complement,
            )
        {
            p.sstart = newstart.getobasedpos().try_into().unwrap_or_default();
            p.send = newend.getobasedpos().try_into().unwrap_or_default();
            p.complement = newcomplement;
        }
    });

    let mut c = c.to_vec();
    c.iter_mut().for_each(|p| {
        if let Some(loc) = p.getlocusname()
            && let Some(find) = locus.iter().find(|fi| {
                p.getchromosomefromsubject().is_some_and(|a| a == fi.contig)
                    && fi.getlocus() == &loc
            })
            && p.sseqid.starts_with("GENE")
            && let Some((start, end, complement)) = p.getpositionfromsubject()
            && let Some((newstart, newend, newcomplement)) = find.locusinposition(
                &start,
                &end,
                &complement,
                &Position::newfromoposition(p.sstart.try_into().unwrap_or_default()),
                &Position::newfromoposition(p.send.try_into().unwrap_or_default()),
                &p.complement,
            )
        {
            p.sstart = newstart.getobasedpos().try_into().unwrap_or_default();
            p.send = newend.getobasedpos().try_into().unwrap_or_default();
            p.qseqid = Name::new(
                Some("GENE".to_string()),
                p.getallelename().to_string(),
                realspecies.to_string(),
                None,
                start,
                end,
                complement,
            );
            p.complement = newcomplement;
        }
    });
    let mut matche = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(dir.join("motifs.txt"))
        .map_err(|e| format!("Error setting motifs match: {e}"))?;
    for matches in motifs.iter() {
        let _ = matche.serialize(matches);
    }
    let _ = matche.flush();
    let mut matche = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(dir.join("motifs2.txt"))
        .map_err(|e| format!("Error setting motifs match: {e}"))?;
    for matches in c.iter() {
        let _ = matche.serialize(matches);
    }
    let _ = matche.flush();
    let file = dir.join("newpotentialalleles.fasta");
    let sequence = filter_new_alleles(&c, &motifs).fold(String::new(), |mut acc, f| {
        let f: &dyn Blastcalc = f;
        let f: &dyn Seqresult = &Newfasta::from(f);
        acc.push_str(&format!("\n{}", f));
        acc
    });
    println!("BLAST results were added.");
    if !sequence.trim().is_empty()
        && let Err(e) = fs::write(file, &sequence)
    {
        eprintln!("An error has occured while priting sequence: {e}.");
        return Ok(());
    }
    if sequence.trim().is_empty() && args.mytoken.is_none() && !asknonewalleles() {
        println!("Exiting");
        return Ok(());
    }
    let link = Path::join(&args.outdir, "submission");
    #[cfg(target_family = "unix")]
    let _ = std::os::unix::fs::symlink(dirtemp.path(), &link);
    #[cfg(target_family = "windows")]
    let _ = std::os::windows::fs::symlink_dir(dirtemp.path(), &link);
    if !args.nosubmit {
        let result = match &args.mytoken {
            Some(_) => Ok(()),
            None => browseropening(),
        };
        if let Err(e) = result {
            let _ = fs::remove_file(link);
            return Err(e.to_string());
        } else {
            preparesubmission(dir, realspecies, args);
        }
    }
    let _ = fs::remove_file(link);
    //let _ = fs::remove_dir_all(dir);
    Ok(())
    //form(&client);
}
pub fn asknonewalleles() -> bool {
    println!(
        "No new alleles has been found. IMGT can still process your data, do you want to continue (Y/n):"
    );
    readfromterminal(&'y', &'n', None, false).is_some_and(|a| a)
}
pub(crate) fn askforsubmission(
    realspecies: &Species,
    locus: &[LocusInfos],
    args: &Args,
    infos: &HashMap<LocusInfos, Vec<Blastmatch>>,
) -> io::Result<()> {
    let quest = REQUESTCLIENT
        .get(SUBMISSIONLINK.as_str())
        .query(&[("ping", "1")]);
    match quest.send() {
        Ok(p) if p.status().is_success() => (),
        Ok(p) if p.status().is_client_error() => {
            /* return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                "IMGT does not currently support automatic submissions.".to_string(),
            )); */

            //return Ok(());
        }
        Ok(p) => {
            return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                format!(
                    "Unable to contact IMGT servers at the moment. Please retry later. Error is {}",
                    p.status().canonical_reason().unwrap_or_default()
                ),
            ));
        }
        Err(e) => {
            return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                format!("Unable to contact IMGT servers. Error is {}", e),
            ));
        }
    };
    if !args.nosubmit {
        if !infos
            .iter()
            .any(|(_, p)| p.iter().any(|p| p.onlynewalleles()))
        {
            if !asknonewalleles() {
                return Ok(());
            }
        } else if args.mytoken.is_none() {
            println!("Do you want to submit your sequences to IMGT (y to yes or n to no)?");
            if !readfromterminal(&'y', &'n', None, false).is_some_and(|a| a) {
                println!("Your sequences won't be submitted.");
                return Ok(());
            }
        }
        let mut blastmatch: Vec<Blastmatch> = Vec::new();
        for data in infos.values() {
            blastmatch.append(&mut data.clone());
        }
        submit(args, locus, &blastmatch, realspecies)
            .map_err(|f| io::Error::new(io::ErrorKind::InvalidInput, f.to_string()))?;
    }
    Ok(())
}
/// Take a compressed GzFile and give a GzFile decompressed
pub(crate) fn decompressseq(infos: &Filecrea) -> io::Result<Filecrea> {
    let mut reader = BufReader::new(infos.getfile()?);
    let new = Filecrea::createtemp(None, Some("sequence.fasta"))?;
    let mut encoder = GzDecoder::new(new.setfile()?);
    let mut count = false;
    loop {
        let mut buf = [0u8; 2048];
        //Read from fasta
        let n = match reader.read(&mut buf)? {
            0 => break,
            n => {
                count = true;
                n
            }
        };
        let newbuf = buf.split_at(n).0;
        encoder.write_all(&newbuf)?;
    }
    if !count {
        return Err(io::Error::new(ErrorKind::InvalidData, "Invalid data"));
    }
    encoder.flush()?;
    let _ = encoder.finish()?;
    Ok(new)
}
pub(crate) fn generatesequenceraw(
    args: &Args,
    locus: &[LocusInfos],
) -> io::Result<io::Cursor<Vec<u8>>> {
    let mut cursor = io::Cursor::new(Vec::new());
    {
        //To free borrow of cursor
        let bufwriter = BufWriter::new(&mut cursor);
        let mut fastawriter = fasta::Writer::from_bufwriter(bufwriter);
        fastawriter.set_linewrap(Some(100));
        let mut assembly = match args.assembly.as_ref().map(|_| getassemblyreader(args)) {
            None => return Err(io::Error::new(InvalidData, "No assembly provided.")),
            Some(Err(e)) => return Err(e),
            Some(Ok(b)) => b,
        };
        for list in locus.iter() {
            let seq = list
                .extractsequence(&mut assembly)
                .unwrap_or("Sequence is unavailable".to_string());
            fastawriter.write(
                &format!("{}:{}", list.getlocus(), list.contig),
                Some(&format!(
                    "{}:{}-{}/{}",
                    list.getlocushaplo(),
                    list.start.getobasedpos(),
                    list.end.getobasedpos(),
                    list.complement
                )),
                seq.as_bytes(),
            )?;
        }
        fastawriter.flush()?;
    }
    cursor.rewind()?;
    Ok(cursor)
}
pub(crate) fn generatesequence(
    args: &Args,
    dir: &Path,
    tmp: bool,
    locus: &[LocusInfos],
) -> Result<Filecrea, String> {
    let sequencefile = if tmp {
        Filecrea::createtemp(None, Some("sequence.fasta.gz")).map_err(|d| d.to_string())?
    } else {
        Filecrea::createfrompath(Path::join(dir, "sequence.fasta.gz"))
    };
    let mut cursor = generatesequenceraw(args, locus).map_err(|d| d.to_string())?;
    let elem = sequencefile
        .setfile()
        .map_err(|d| format!("Cannot create archive: {d}"))?;
    let mut encoder = GzBuilder::default()
        .filename("sequence.fasta")
        .comment("Sequence")
        .write(elem, Compression::fast());
    loop {
        let mut buf = [0u8; 2048];
        let n = match cursor
            .read(&mut buf)
            .map_err(|e| format!("Error making fasta archive: {e}"))?
        {
            0 => break,
            n => n,
        };
        let newbuf = buf.split_at(n).0;
        encoder
            .write_all(&newbuf)
            .map_err(|d| format!("Cannot write sequence. Error is {d}"))?;
    }
    encoder.flush().map_err(|_| format!("Error flushing"))?;
    encoder.finish().map_err(|_| format!("Error gzip"))?;
    Ok(sequencefile)
}
pub(crate) fn getindexforbam(outpath: &Path) -> Option<PathBuf> {
    let mut f = outpath.to_path_buf();
    if f.add_extension("csi") {
        Some(f)
    } else {
        None
    }
}
pub(crate) fn generatelightbam<T>(
    args: &Args,
    light: T,
    lightindex: Option<T>,
    locus: &[LocusInfos],
) -> Result<(), String>
where
    T: AsRef<Path>,
{
    println!("Generating small BAM for submission");
    let bam = if let Ok(r) = getreaderoffile(args) {
        r
    } else {
        return Err("Cannot access BAM file for light bam.".to_string());
    };
    {
        let mut writer = if let Ok(files) = bam::Writer::from_path(
            light.as_ref(),
            &bam::Header::from_template(bam.header()),
            bam::Format::Bam,
        ) {
            files
        } else {
            let file = light.as_ref().display();
            return Err(format!("Cannot create file {file} for light bam."));
        };
        for f in locus.iter() {
            let mut bam = if let Ok(r) = getreaderoffile(args) {
                r
            } else {
                return Err("Cannot access BAM file for light bam.".to_string());
            };
            if bam
                .fetch((
                    f.contig.as_bytes(),
                    f.start.getzbasedpos(),
                    f.end.getzbasedpos().saturating_add(1),
                ))
                .is_err()
            {
                return Err("Cannot read BAM file region for light bam.".to_string());
            }
            for read in bam.rc_records().filter_map(Result::ok) {
                if writer.write(&read).is_err() {
                    return Err("Cannot read BAM file region for light bam.".to_string());
                };
            }
        }
    }
    //Drop writer else there is an issue when making index
    println!("Building index");
    bam::index::build(
        light.as_ref(),
        lightindex.as_ref().map(|a| a.as_ref()),
        Type::Csi(2_u32.pow(14)),
        available_parallelism()
            .map(|a| a.get())
            .unwrap_or(1)
            .try_into()
            .unwrap_or(1),
    )
    .map_err(|e| format!("Cannot build index. Error is {}", e))?;
    Ok(())
}
pub(crate) fn browseropening() -> io::Result<()> {
    let link = SUBMISSIONLINK.as_str();
    println!(
        "Opening web browser to continue submission. Type (Y) to open the web browser, (e) to exit or (n) and go by yourself to {link}."
    );
    match readfromterminal(&'y', &'e', Some(&'n'), false) {
        Some(true) => webbrowser::open(link),
        None => Ok(()),
        Some(false) => Err(io::Error::new(
            ErrorKind::ConnectionAborted,
            "Your sequences won't be submitted",
        )),
    }
}
pub(crate) fn createarchive(args: &Args, dir: &Path) -> io::Result<Filecrea> {
    let temp = Filecrea::createtemp(Some(dir), Some(Path::new("submission.tar.gz")))?;
    let file = File::create(temp.getpath())?;
    let archive = GzEncoder::new(file, Compression::best());
    let mut tar = tar::Builder::new(archive);
    tar.follow_symlinks(false);
    tar.append_dir_all(".", &args.outdir)?;
    tar.finish()?;
    Ok(temp)
}
pub(crate) fn preparesubmission(path: &Path, species: &Species, args: &Args) -> bool {
    let token = if let Some(a) = &args.mytoken {
        a.to_string()
    } else {
        println!(
            "Please give the token provided by the submission. Type (e) to exit. Analysis files and new sequences would be sent to IMGT for submission."
        );
        loop {
            let mut block = io::stdin().lock().take(100);
            let mut string = String::new();
            let val = block.read_line(&mut string);
            match (val, string.trim().to_lowercase()) {
                (Err(_), _) | (Ok(0), _) => {
                    eprintln!("Exiting");
                    return false;
                }
                (Ok(_), v) if v == "e" => {
                    eprintln!("Exiting");
                    return false;
                }
                (Ok(_), v) if v.chars().all(|f| f.is_ascii_hexdigit()) && v.len() == 10 => break v,
                (Ok(_), _) => {
                    eprintln!("Token is 10 hexadecimal characters, please retry.");
                }
            };
        }
    };
    let archive = match createarchive(args, path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error filling archive: {e}");
            return false;
        }
    };
    println!("Submission in progress");
    if let Err(e) = submission(&token, species, archive) {
        eprintln!("An error has occured during submission: {e}. Please retry later.");
        return false;
    }
    println!(
        "Your submission has been made successfully. Thank you for submitting your sequences to IMGT. A confirmation email has been sent."
    );
    true
}
pub(crate) fn getprogressbarspin() -> io::Result<ProgressBar> {
    let pb = ProgressBar::new(0);
    pb.enable_steady_tick(Duration::from_millis(100));
    pb.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos}")
            .map_err(|b| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("issue with progress bar: {b}"),
                )
            })?,
        //.progress_chars("#>-"),
    );
    Ok(pb)
}
pub(crate) fn getprogressbarclassic(total_size: u64) -> io::Result<ProgressBar> {
    let pb = ProgressBar::new(total_size);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) ({eta}) {msg}",
            )
            .map_err(|b| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("issue with progress bar: {b}"),
                )
            })?
        //.progress_chars("#>-"),
    );
    Ok(pb)
}
pub(crate) fn getprogressbar<R>(size: u64, reader: R) -> io::Result<ProgressReader<R>> {
    let pb = ProgressBar::new(size);
    pb.set_style(
    indicatif::ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) {msg}")
        .map_err(|b| io::Error::new(ErrorKind::InvalidData,format!("issue with progress bar: {b}")))?
    //.progress_chars("#>-"),
);

    Ok(ProgressReader {
        reader,
        progress_bar: pb,
        total_bytes: size,
    })
}
pub(crate) fn submission(token: &str, species: &Species, archive: Filecrea) -> io::Result<()> {
    /* let multipart = reqwest::blocking::multipart::Form::new()
    .file("genelist", "submission/genelist.csv")?
    .file("sequences", "submission/sequences.txt")?
    .file("locuspos", "submission/locus.txt")?
    .text("type","submission")
    .file("locus", "submission/locus.bam")?; */
    let file_size = archive.getpath().metadata()?.len();
    let progress_reader = getprogressbar(file_size, archive.getfile()?)?;
    let zip = reqwest::blocking::multipart::Form::new()
        .part("file", multipart::Part::reader(progress_reader))
        .text("version", VERSION)
        .text("type", "submission")
        .text("species", species.to_string())
        .text("validspecies", species.ischecked().to_string());
    match REQUESTCLIENT
        .post(SUBMISSIONLINK.as_str())
        .bearer_auth(token)
        .multipart(zip)
        .send()
    {
        Ok(a) => match a.status() {
            StatusCode::OK | StatusCode::NO_CONTENT => Ok(()),
            StatusCode::UNAUTHORIZED | StatusCode::FORBIDDEN => Err(io::Error::new(
                io::ErrorKind::PermissionDenied,
                format!(
                    "The token is invalid or has expired. Please retry a submission. Response is {}",
                    a.text().unwrap_or_default()
                ),
            )),
            StatusCode::BAD_REQUEST => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "The data sent is invalid, please check submission and retry. Response is {}",
                    a.text().unwrap_or_default()
                ),
            )),
            a if a.is_server_error() => Err(io::Error::new(
                io::ErrorKind::ConnectionAborted,
                "The server is unavailable. Please retry a submission",
            )),
            e => Err(io::Error::new(
                io::ErrorKind::ConnectionReset,
                format!("An unexpected error has occured, please retry later. Error is {e}"),
            )),
        },
        Err(e) if e.is_timeout() => Err(io::Error::new(
            io::ErrorKind::TimedOut,
            "The submission has timed out. Please retry later.",
        )),
        Err(e) => Err(io::Error::new(
            io::ErrorKind::NetworkUnreachable,
            format!("An unexpected error has occured. Please retry later. Error is {e}"),
        )),
    }?;
    //std::fs::remove_file(archive)?;
    Ok(())
}

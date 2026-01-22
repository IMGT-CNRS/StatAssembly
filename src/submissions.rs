#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
use bio::io::fasta;
use flate2::{Compression, write::GzEncoder};
use itertools::Itertools;
use lazy_static::lazy_static;
use reqwest::{StatusCode, tls};
use serde::{Deserialize, Serialize};
use serde_json::{self as json, Value};
use std::cmp::{Ordering, max, min};
use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::str::FromStr;
use std::{
    borrow::Cow,
    env::{self, current_dir, temp_dir},
    error::Error,
    fs::{self, File},
    io::{self, BufRead, BufReader, ErrorKind, Read, Write},
    ops::{Not, Range, RangeInclusive},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{Duration, SystemTime},
};
use strum::IntoEnumIterator;
use tempfile::NamedTempFile;

use crate::r#struct::{
    Args, Blast, Blastcalc, Blastmatch, GeneInfos, Haplotype, Locus, LocusInfos, Newfasta,
    Ourfasta, Position, Seqresult, Status, Strand,
};
use crate::{generatelightbam, getassemblyreader, getreaderoffile};
lazy_static! {
    pub static ref REQUESTCLIENT: reqwest::blocking::Client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::new(10, 0))
        .timeout(Duration::new(300, 0))
        .user_agent(format!("IMGT/StatAssembly version {}", VERSION))
        .referer(false)
        .tls_version_min(tls::Version::TLS_1_2)
        .https_only(true)
        .build()
        .unwrap_or_default();
}
const VERSION: &str = env!("CARGO_PKG_VERSION");
const RELEASELINK: &str = "https://www.imgt.org/download/GENE-DB/RELEASE";
const VQUESTLINK: &str = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP";
const SUBMISSIONLINK: &str = "https://www.imgt.org/submissions/";
pub(crate) const DELIMITERFASTA: char = '/';
const LOCUSSEPARATOR: usize = 1_000_000;

#[must_use]
pub(crate) fn sendresult(request: &reqwest::blocking::Client, url: &str) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                Ok(e.text().unwrap_or("Error getting data".to_string()))
            } else {
                Err(e
                    .error_for_status()
                    .map_or("Error getting URL".to_string(), |f| {
                        format!("{}", f.status())
                    }))
            }
        }
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn getnamefromblast(text: &str) -> Option<String> {
    text.to_ascii_lowercase()
        .split('|')
        .nth(2)
        .map(|f| f.to_string())
}
pub(crate) fn getallelefromblast(text: &str) -> Option<String> {
    text.to_ascii_uppercase()
        .split('|')
        .nth(1)
        .map(|f| f.to_string())
}
#[must_use]
pub(crate) fn fastafilter(text: &str, find: &str, present: bool) -> String {
    let mut lines: Vec<&str> = Vec::with_capacity(text.lines().count());
    let mut keep = true;
    for line in text.lines() {
        let result = line.trim();
        if result.starts_with(">") {
            keep = if present {
                getnamefromblast(result).is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
            } else {
                !getnamefromblast(result).is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
            };
        }
        if keep {
            lines.push(result);
        }
    }
    lines
        .iter()
        .fold(String::new(), |mut acc, f| {
            acc.push_str(format!("\n{}", f).as_str());
            acc
        })
        .trim()
        .to_string()
}
pub(crate) fn checkoverlap(a: &RangeInclusive<usize>, b: &RangeInclusive<usize>) -> bool {
    max(a.start(), b.start()) <= min(a.end(), b.end())
}
pub(crate) fn checkifblastpresent() -> bool {
    println!("Detect if BLAST is operating.");
    let command = Command::new("blastn")
        .arg("-version")
        .stdout(Stdio::piped())
        .stdin(Stdio::piped())
        .status()
        .is_ok_and(|f| f.success());
    if !command {
        eprintln!("BLAST was not found. Check if present in PATH.");
        return false;
    } else {
        println!("BLAST is working. Continuing");
        return true;
    }
}
pub(crate) fn downloadref() -> Option<(PathBuf, String)> {
    println!("Downloading reference sequence from GENE-DB");
    let releaseversion = match sendresult(&REQUESTCLIENT, RELEASELINK) {
        Ok(e) => e,
        Err(e) => {
            println!("Release fetched failed because: {e}");
            return None;
        }
    };
    let tempfile = Path::join(&env::temp_dir(), format!("refseq{}.fasta", releaseversion));
    if !tempfile.is_file() {
        println!(
            "Downloading GENE-DB release {}, please wait...",
            releaseversion
        );
        match sendresult(&REQUESTCLIENT, VQUESTLINK) {
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
                println!("V-QUEST data retrieval failed because: {e}.");
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
pub(crate) fn locusfiltering(locus: &Locus, blast: &mut Vec<Blast>) {
    //TRD is inside TRA locus
    if locus == &Locus::TRA {
        blast.retain(|p| {
            p.qseqid.contains(&format!("{}", locus)) || p.qseqid.contains(&"TRD".to_string())
        });
    } else {
        blast.retain(|p| p.qseqid.contains(&format!("{}", locus)));
    }
}
pub(crate) fn speciesandorphonfiltering(
    tempfile: &Path,
    releaseversion: String,
    species: &str,
    orphonfilter: bool,
) -> io::Result<PathBuf> {
    println!("Filtering based on species {}.", species);
    let file = std::fs::read_to_string(&tempfile)?;
    let info = fastafilter(&file, species, true).replace(" ", "_");
    let info = if orphonfilter {
        fastafilter(&info, "/OR", false)
    } else {
        info
    };
    let info = if info.is_empty() {
        println!("New species!!");
        file
    } else {
        info
    };
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
    let tempfile = Path::join(
        &env::temp_dir(),
        format!(
            "refseq{}-{}.fasta",
            releaseversion.replace(" ", "_"),
            species.replace(" ", "_")
        ),
    );
    if let Err(e) = std::fs::write(&tempfile, &info) {
        return Err(e);
    }
    return Ok(tempfile);
}
pub(crate) fn readfastareader<T>(fasta: fasta::Reader<T>) -> io::Result<Vec<Ourfasta>>
where
    T: std::io::BufRead,
{
    let mut seqs = Vec::new();
    for record in fasta.records().filter_map(|p| p.ok()) {
        if let Err(v) = record.check() {
            eprintln!("Seq {} is invalid, skipped. Error is {v}", record.id());
            continue;
        }
        seqs.push(Ourfasta {
            name: format!("{} {}", record.id(), record.desc().unwrap_or_default()),
            seq: String::from_utf8_lossy(record.seq()).to_string(),
        })
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
pub(crate) fn readfastafile(seq: &Path) -> io::Result<Vec<Ourfasta>> {
    let file = File::open(seq)?;
    let fasta = fasta::Reader::new(file);
    readfastareader(fasta)
}
#[must_use]
#[allow(unused)]
pub(crate) fn selectnewalleles(result: &[Blast]) -> Vec<Ourfasta> {
    let mut fastas = Vec::new();
    for seq in result {
        let name = Ourfasta {
            name: seq.qseqid.clone(),
            seq: seq.sseq.clone(),
        };
        fastas.push(name);
    }
    fastas
}
pub(crate) fn statusblastvs(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.sseqid.cmp(&b.sseqid) {
        std::cmp::Ordering::Equal => a.qseqid.cmp(&b.qseqid),
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
                    && checkoverlap(&g.pos(), &f.pos())
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn find_global_best_range(data: &[Blastmatch]) -> Option<Vec<LocusInfos>> {
    let mut groups: HashMap<(Locus, Haplotype), (String, bool, usize, usize)> = HashMap::new();
    let blastcheck = data.into_iter();
    let blastcheck = blastcheck.sorted_unstable_by(|a, b| match a.sseqid.cmp(&b.sseqid) {
        Ordering::Equal => match a.sstart.cmp(&b.sstart) {
            Ordering::Equal => b.sstart.cmp(&b.sstart),
            ord => ord,
        },
        ord => ord,
    });
    let hash = blastcheck
        .filter_map(|p| {
            match Locus::try_from(
                p.getallelename()
                    .unwrap_or_default()
                    .chars()
                    .take(3)
                    .collect::<String>(),
            ) {
                Ok(d) => Some((d, p)),
                Err(_) => None,
            }
        })
        .into_group_map_by(|(l, f)| (l.clone(), f.sseqid.clone()));
    let mut locus: HashMap<(Locus, String), BinaryHeap<&Blastmatch>> = HashMap::new();
    for ((loci, sseqname), data) in hash.iter() {
        for (_, elem) in data {
            let d = match locus.get_mut(&(loci.clone().clone(), sseqname.clone())) {
                Some(a) => a,
                None => {
                    let mut e = BinaryHeap::new();
                    e.push(elem.clone());
                    locus.insert((loci.clone().clone(), sseqname.clone()), e);
                    continue;
                }
            };
            if let Some(b) = d.peek()
                && max(b.sstart, b.send).abs_diff(min(elem.send, elem.sstart)) <= LOCUSSEPARATOR
            {
                d.push(&elem);
            }
        }
    }
    let real = locus
        .into_iter()
        .sorted_unstable_by(|((al, ass), av), ((bl, bs), bv)| match al.cmp(&bl) {
            Ordering::Equal => match av.len().cmp(&bv.len()) {
                Ordering::Equal => ass.cmp(&bs),
                ord2 => ord2.reverse(),
            },
            ord => ord,
        });
    for ((loci, sseq), blast) in real {
        if !groups.contains_key(&(loci.clone(), Haplotype::Primary)) {
            let (min, max) = match (blast.iter().next(), blast.peek()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| f.sstart > f.send);
            let complement = if split.get(&true).map_or(0, |s| s.len())
                > split.get(&false).map_or(0, |s| s.len())
            {
                true
            } else {
                false
            };
            groups.insert((loci, Haplotype::Primary), (sseq, complement, min, max));
        } else if groups.contains_key(&(loci.clone(), Haplotype::Alternate)) {
            let (min, max) = match (blast.iter().next(), blast.peek()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| f.sstart > f.send);
            let complement = if split.get(&true).map_or(0, |s| s.len())
                > split.get(&false).map_or(0, |s| s.len())
            {
                true
            } else {
                false
            };
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
                    Position::new(false, start.try_into().unwrap_or_default()),
                    Position::new(false, end.try_into().unwrap_or_default()),
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
pub(crate) fn retainbestmatch(blast: &mut Vec<Blast>) {
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
}
pub(crate) fn genesblast<T>(
    subject: &[GeneInfos],
    args: &Args,
    species: T,
    locus: &Locus,
) -> io::Result<Vec<Blastmatch>>
where
    T: AsRef<str>,
{
    let mut reader = getassemblyreader(&args)?;
    let name = temp_dir().join("genes_blast.txt");
    let file = File::create(&name)?;
    let mut fastawriter = fasta::Writer::new(file);
    subject
        .iter()
        .map(|f| {
            let elem = f
                .extractsequence(&mut reader)
                .map_err(|p| io::Error::new(ErrorKind::InvalidInput, p));
            let bool = elem
                .iter()
                .any(|p| f.addtosequence(p, &mut fastawriter).is_err());
            let elem = if bool {
                Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    "Cannot print sequence",
                ))
            } else {
                elem
            };
            elem
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reference = match downloadref()
        .map(|(a, b)| speciesandorphonfiltering(&a, b, species.as_ref(), false))
    {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast: Vec<Blast> = match blastcommand(reference.as_path(), name.as_path()) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    locusfiltering(&locus, &mut blast);
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
    });
    statusblastvs(&mut blast);
    Ok(blast.into_iter().map(|f| f.into()).collect())
}
pub(crate) fn locusposition<T>(
    subject: &Path,
    species: T,
    locus: &Locus,
    full: bool,
) -> io::Result<(Vec<LocusInfos>, Vec<Blastmatch>)>
where
    T: AsRef<str>,
{
    let reference = downloadref();
    let reference =
        match reference.map(|(a, b)| speciesandorphonfiltering(&a, b, species.as_ref(), false)) {
            Some(a) => a?,
            None => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Reference from IMGT cannot be downloaded",
                ));
            }
        };
    let mut blast: Vec<Blast> = match blastcommand(reference.as_path(), subject) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    if !full {
        locusfiltering(&locus, &mut blast);
    }
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
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
    let mut data: Vec<Blastmatch> = statusvec
        .into_iter()
        .map(|(_k, v)| Blastmatch::new(v.qseqid, v.sseqid, v.sseq, v.sstart, v.send, v.status))
        .collect();
    //Sort by name then starting then ending position
    data.sort_unstable();
    data.dedup();
    let range = find_global_best_range(&data).ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No locus found after BLAST analysis",
    ));
    if !full && let Ok(e) = &range {
        data.retain(|p| {
            e.iter().any(|f| {
                &f.locus == locus
                    && f.contig == p.qseqid
                    && (f.start.getobasedpos()..=f.end.getobasedpos())
                        .contains(&p.sstart.try_into().unwrap_or_default())
            })
        })
    }
    range.map(|f| (f, data))
}
#[must_use]
pub(crate) fn filter_new_alleles<T>(data: &[T]) -> impl Iterator<Item = &T>
where
    T: Blastcalc,
{
    data.iter().filter(|p| p.onlynewalleles())
}
pub(crate) fn blastcommand<T>(reference: T, subject: T) -> io::Result<Vec<Blast>>
where
    T: AsRef<Path>,
{
    let (reference, subject) = (reference.as_ref(), subject.as_ref());
    if !reference.exists() {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Reference file was not found",
        ));
    }
    if !subject.exists() {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Subject file was not found",
        ));
    }
    let reference = &format!("{}", reference.display()).replace(" ", "_");
    let subject = &format!("{}", subject.display()).replace(" ", "_");
    let output = Path::join(&env::temp_dir(), "blast.txt");
    let output = &format!("{}", output.display());
    let command = Command::new("blastn")
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
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "5",
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send qlen length pident gaps sseq",
        ])
        .spawn()?;
    println!("Launching {} against {}", reference, subject);
    println!(
        "BLAST has been launched with id {}. Please wait.",
        command.id()
    );
    let outputc = command.wait_with_output()?;
    let mut result: Vec<Blast> = Vec::new();
    {
        if !outputc.status.success() {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!(
                    "BLAST has failed with status {}. Error is {}",
                    outputc.status,
                    &String::from_utf8_lossy(&outputc.stderr)
                ),
            ));
        } else {
            println!("BLAST has been done. Parsing.");
            let file = File::open(output)?;
            let reader = BufReader::new(file);
            let mut csv = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(reader);
            for record in csv.deserialize() {
                if let Ok(r) = record {
                    result.push(r);
                } else if let Err(r) = record {
                    eprintln!("Error in {r}");
                }
            }
        };
    };
    //let _ = fs::remove_file(&output);
    Ok(result)
}
pub(crate) fn getspeciesfromncbi<T>(
    client: &reqwest::blocking::Client,
    species: &T,
) -> Result<String, Box<dyn Error>>
where
    T: AsRef<str>,
{
    let species = species.as_ref();
    let val = if let Ok(val) = species.parse::<usize>() {
        val
    } else {
        let response = client
            .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?")
            .query(&[("db", "taxonomy"), ("term", species), ("format", "json")])
            .send()?;
        let jsone: json::Value = json::from_str(&response.text().unwrap_or(String::new()))?;
        jsone["esearchresult"]["idlist"]
            .as_array()
            .map(|f| f.iter().next())
            .ok_or(io::Error::new(ErrorKind::InvalidInput, "Idlist not found"))?
            .ok_or(io::Error::new(ErrorKind::InvalidInput, "Idlist not found"))?
            .as_str()
            .ok_or(io::Error::new(
                ErrorKind::InvalidInput,
                "No result for the term used",
            ))?
            .parse::<usize>()
            .unwrap_or(0)
    };
    let val = &format!("{}", val);
    let response = client
        .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")
        .query(&[("db", "taxonomy"), ("id", val), ("format", "json")])
        .send()?;
    let jsone: json::Value = json::from_str(&response.text().unwrap_or(String::new()))?;
    let mapping = jsone["result"].as_object().ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No result for the term used",
    ))?;
    if let Some(data) = mapping.get(val) {
        let elem = data.as_object().ok_or(io::Error::new(
            ErrorKind::InvalidInput,
            "No result for the id found",
        ))?;
        let text = String::from("species");
        match (elem.get("rank"), elem.get("scientificname")) {
            (Some(Value::String(text2)), Some(Value::String(name)))
                if text2.as_str() == text.as_str() =>
            {
                Ok(name.to_string())
            }
            (Some(Value::String(rank)), _) => Err(Box::new(io::Error::new(
                ErrorKind::InvalidInput,
                format!("The term used is not a species but a {}", rank),
            ))),
            _ => Err(Box::new(io::Error::new(
                ErrorKind::InvalidInput,
                "Invalid term used",
            ))),
        }
    } else {
        Err(Box::new(io::Error::new(
            ErrorKind::InvalidInput,
            "Invalid mapping",
        )))
    }
}
#[must_use]
pub(crate) fn launchblast(species: &str, subject: &Path) -> Result<Vec<Newfasta>, ()> {
    let realspecies = match getspeciesfromncbi(&REQUESTCLIENT, &species) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    println!("The species is {}.", realspecies);
    let (path, releaseversion) = match downloadref() {
        None => return Err(()),
        Some((a, b)) => (a, b),
    };
    let filtering = match speciesandorphonfiltering(&path, releaseversion, &realspecies, false) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    let result: Vec<Newfasta> = match blastcommand(filtering.as_path(), subject) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(mut c) => {
            statusblast(&mut c);
            let status: Vec<Newfasta> = c.iter().map(|b| Newfasta::newfromblast(b)).collect();
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
}
pub(crate) fn submit(
    args: &Args,
    locus: &[crate::LocusInfos],
    c: &[Blastmatch],
    realspecies: String,
) -> Result<(), String> {
    //let result: Vec<Newfasta> = c.into_iter().map(Newfasta::newfromblastowner).collect();
    let dir = Path::new(&current_dir().unwrap_or(temp_dir())).join("archive");
    if dir.is_dir() {
        eprintln!("Archive directory exists, going to be deleted.");
        if let Err(e) = fs::remove_dir_all(&dir) {
            let dir = dir.display();
            return Err(format!("Cannot remove the directory {dir}, error is {e}."));
        }
    }
    if let Err(e) = fs::create_dir(&dir) {
        return Err(format!("Cannot create archive directory, error is {e}"));
    }
    let _ = match args.assembly.as_ref().map(|_| getassemblyreader(&args)) {
        None => return Err("No assembly provided.".to_string()),
        Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
        Some(Ok(b)) => b,
    };
    let locuspos =
        File::create(dir.join("newloc.csv")).map_err(|f| format!("Locus csv, error is {f}"))?;
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(locuspos);
    for loci in locus.iter() {
        csv.serialize(loci)
            .map_err(|p| format!("Error serializing locus position: {p}"))?;
    }
    let lightbam = dir.join("outlight.bam");
    generatelightbam(args, &lightbam, locus)?;
    let sequencefile = dir.join("sequence.fasta");
    let mut fastawriter = fasta::Writer::to_file(sequencefile).map_err(|f| format!("{f}"))?;
    for list in locus.iter() {
        let mut assembly = match args.assembly.as_ref().map(|_| getassemblyreader(&args)) {
            None => return Err("No assembly provided.".to_string()),
            Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
            Some(Ok(b)) => b,
        };
        let seq = list
            .extractsequence(&mut assembly)
            .map_err(|f| format!("{f}"))?;
        fastawriter
            .write(
                &format!("{}", list.locus),
                Some(&format!(
                    "{}:{}-{}/{}/{}",
                    list.contig,
                    list.start.getobasedpos(),
                    list.end.getobasedpos(),
                    list.complement,
                    list.haplotype
                )),
                seq.as_bytes(),
            )
            .map_err(|f| {
                format!(
                    "Unable to write fasta sequence {}. Error is {f}",
                    list.locus
                )
            })?;
    }
    let file = dir.join("newalleles.fasta");
    let sequence = filter_new_alleles(c).fold(String::new(), |mut acc, f| {
        let f: &dyn Seqresult = f;
        acc.push_str(&format!("\n{}", f));
        acc
    });
    if let Err(e) = fs::write(file, sequence) {
        eprintln!("An error has occured while priting sequence: {e}.");
        return Ok(());
    }
    println!("BLAST results was added.");
    browseropening().map_err(|f| f.to_string())?;
    preparesubmission(&dir, realspecies);
    let _ = fs::remove_dir_all(dir);
    Ok(())
    //form(&client);
}
pub(crate) fn browseropening() -> io::Result<()> {
    println!(
        "Opening web browser to continue submission. Type (Y) to open the web browser, (N) to refuse and go by yourself to {SUBMISSIONLINK} or (e) to exit."
    );
    let mut val = String::new();
    let _ = io::stdin().read_line(&mut val);
    let val = val.trim().to_ascii_lowercase();
    if val == "y" {
        let _ = webbrowser::open(SUBMISSIONLINK);
    } else if val == "e" {
        println!("Exiting");
        return Err(io::Error::from(ErrorKind::InvalidData));
    }
    Ok(())
}
pub(crate) fn createarchive(dir: &Path) -> io::Result<NamedTempFile> {
    let temp = tempfile::NamedTempFile::with_suffix("submission.tar.gz")?;
    let file = File::create(&temp)?;
    let archive = GzEncoder::new(file, Compression::best());
    let mut tar = tar::Builder::new(archive);
    tar.append_dir_all(
        "",
        dir.file_name()
            .map(|f| f.to_str())
            .unwrap_or_default()
            .unwrap_or_default(),
    )?;
    tar.finish()?;
    Ok(temp)
}
pub(crate) fn preparesubmission(path: &Path, species: String) -> bool {
    println!("Please give the token provided by the submission. Type (e) to exit.");
    let token = loop {
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
    };
    println!(
        "Analysis files and new sequences would be sent to IMGT for submission. Type (Y) to validate or something else to exit"
    );
    let mut val = String::new();
    let _ = io::stdin().read_line(&mut val);
    let val = val.trim().to_ascii_lowercase();
    let archive = match createarchive(path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error filling archive: {e}");
            return false;
        }
    };
    if val == "y" {
        if let Err(e) = submission(&token, species, archive) {
            eprintln!("An error has occured during submission: {e}. Please retry later.");
            return false;
        }
        println!(
            "Your submission has been made successfully. Thank you for submitting your sequences to IMGT. A confirmation email has been sent."
        );
    } else {
        println!("Exiting");
        return false;
    }
    true
}
pub(crate) fn submission(token: &str, species: String, archive: NamedTempFile) -> io::Result<()> {
    /* let multipart = reqwest::blocking::multipart::Form::new()
    .file("genelist", "submission/genelist.csv")?
    .file("sequences", "submission/sequences.txt")?
    .file("locuspos", "submission/locus.txt")?
    .text("type","submission")
    .file("locus", "submission/locus.bam")?; */
    let zip = reqwest::blocking::multipart::Form::new()
        .file("archive", archive.path())?
        .text("version", VERSION)
        .text("type", "submission")
        .text("species", species);
    match REQUESTCLIENT
        .post(SUBMISSIONLINK)
        .bearer_auth(token)
        .multipart(zip)
        .send()
    {
        Ok(a) => match a.status() {
            StatusCode::OK | StatusCode::NO_CONTENT => Ok(()),
            StatusCode::UNAUTHORIZED => Err(io::Error::new(
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
                io::ErrorKind::ResourceBusy,
                "The server is unavailable. Please retry a submission.",
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
    std::fs::remove_file(archive)?;
    Ok(())
}

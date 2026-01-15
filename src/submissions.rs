#[deny(clippy::unwrap_used)]
#[deny(clippy::expect_used)]
use bio::io::fasta;
use flate2::{Compression, write::GzEncoder};
use lazy_static::lazy_static;
use reqwest::{StatusCode, tls};
use serde::{Deserialize, Serialize};
use serde_json::{self as json, Value};
use std::cmp::{Ordering, max, min};
use std::collections::{BTreeMap, HashMap};
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
use tempfile::NamedTempFile;

use crate::r#struct::Locus;
lazy_static! {
    pub static ref REQUESTCLIENT: reqwest::blocking::Client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::new(10, 0))
        .timeout(Duration::new(300, 0))
        .user_agent(format!("IMGT/StatAssembly version {}", VERSION))
        .referer(false)
        .tls_version_min(tls::Version::TLS_1_2)
        .https_only(true)
        .build()
        .unwrap();
}
const VERSION: &str = env!("CARGO_PKG_VERSION");
const RELEASELINK: &str = "https://www.imgt.org/download/GENE-DB/RELEASE";
const VQUESTLINK: &str = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP";
const SUBMISSIONLINK: &str = "https://www.imgt.org/submissions/";
const DELIMITERFASTA: char = '/';
#[derive(Clone, Debug, Deserialize)]
pub(crate) struct Blast {
    qseqid: String,
    sseqid: String,
    qstart: usize,
    qend: usize,
    sstart: usize,
    send: usize,
    qlen: usize,
    length: usize,
    pident: f32,
    gaps: usize,
    sseq: String,
    #[serde(skip_deserializing)]
    complement: bool,
    #[serde(skip_deserializing)]
    status: Status,
}
#[derive(Clone, Debug)]
pub(crate) struct Blastmatch {
    qseqid: String,
    sseqid: String,
    sseq: String,
    sstart: usize,
    send: usize,
    status: Status,
}
impl PartialOrd for Blastmatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}
impl Ord for Blastmatch {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.sseqid.cmp(&other.sseqid) {
            std::cmp::Ordering::Equal => match self.sstart.cmp(&other.sstart) {
                std::cmp::Ordering::Equal => self.send.cmp(&other.send),
                ord => ord,
            },
            ord => ord,
        }
    }
}
impl PartialEq for Blastmatch {
    fn eq(&self, other: &Self) -> bool {
        self.qseqid == other.qseqid && self.sstart == other.sstart && self.send == other.send
    }
}
impl Eq for Blastmatch {}
impl Blastmatch {
    fn new(
        qseqid: String,
        sseqid: String,
        sseq: String,
        sstart: usize,
        send: usize,
        status: Status,
    ) -> Self {
        Self {
            qseqid,
            sseqid,
            sseq,
            sstart,
            send,
            status,
        }
    }
}
impl Blast {
    fn getstatus(&self) -> &Status {
        &self.status
    }
    fn setstatus(&mut self) {
        match (self.pident, self.qlen, self.length) {
            (100.0, a, b) if a == b => self.status = Status::Equal,
            (100.0, ..) => self.status = Status::Shorter,
            _ => self.status = Status::New,
        }
    }
}
impl ToString for Blastmatch {
    fn to_string(&self) -> String {
        return format!(
            ">{}{DELIMITERFASTA}{}:{}-{}{DELIMITERFASTA}{}\n{}",
            self.qseqid,
            self.sseqid,
            self.sstart,
            self.send,
            self.status.to_string(),
            self.sseq
        )
        .to_string();
    }
}
impl FromStr for Blastmatch {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if !s.starts_with(">") {
            return Err("No header format".to_string());
        }
        let s = s.trim_start_matches('>');
        let (split1, seq) = match s.split_once('\n') {
            Some((a, b)) => (a, b),
            None => return Err("No status in header format: {s}".to_string()),
        };
        if !seq.chars().all(|p| p.is_ascii_alphabetic()) {
            return Err("Absence of ASCII alphabetic in sequence: {s}".to_string());
        }
        let (qseqid, sseqid, sstart, send, status) = match split1.splitn(3, DELIMITERFASTA) {
            mut a if a.clone().count() == 3 => {
                let (name, infos, status) =
                    (a.next().unwrap(), a.next().unwrap(), a.next().unwrap());
                let status = if let Ok(a) = Status::try_from(status) {
                    a
                } else {
                    return Err("Invalid status in header: {s}".to_string());
                };
                let (sseqid, sstart, send) = if let Some((a, b, c)) =
                    infos.split_once(':').and_then(|p| {
                        if let Some((a, b)) = p.1.split_once('-') {
                            let (start, end) = match (a.parse::<usize>(), b.parse::<usize>()) {
                                (Ok(b), Ok(c)) => (b, c),
                                _ => return None,
                            };
                            Some((p.0, start, end))
                        } else {
                            None
                        }
                    }) {
                    (a, b, c)
                } else {
                    return Err("Invalid header in sequence header: {s}".to_string());
                };
                (name.to_string(), sseqid, sstart, send, status)
            }
            _ => return Err("Lacking info in sequence header: {s}".to_string()),
        };
        Ok(Self {
            qseqid: qseqid.to_string(),
            sseqid: sseqid.to_string(),
            sseq: seq.to_string(),
            sstart,
            send,
            status,
        })
    }
}
impl PartialEq for Blast {
    fn eq(&self, other: &Self) -> bool {
        self.qseqid == other.qseqid && self.sstart == other.sstart && self.send == other.send
    }
}
impl Eq for Blast {}
pub trait Blastcalc {
    fn getseq<'a>(&self) -> Cow<str>;
    fn getstatus(&self) -> &Status;
    /// Only return new alleles
    fn onlynewalleles(&self) -> bool {
        self.getstatus() == &Status::New
    }
}
impl Blastcalc for Blast {
    fn getseq(&self) -> Cow<str> {
        Cow::Borrowed(&self.sseq)
    }
    fn getstatus(&self) -> &Status {
        &self.status
    }
}
impl Blastcalc for Blastmatch {
    fn getseq(&self) -> Cow<str> {
        Cow::Borrowed(&self.sseq)
    }
    fn getstatus(&self) -> &Status {
        &self.status
    }
}
#[derive(Clone, Debug)]
pub(crate) struct Newfasta {
    qseqid: String,
    sseq: String,
    pos: RangeInclusive<usize>,
    status: Status,
}
impl Seqresult for Newfasta {
    fn qseqid(&self) -> &str {
        &self.qseqid
    }

    fn sseq(&self) -> &str {
        &self.sseq
    }

    fn status(&self) -> &Status {
        &self.status
    }
    fn pos(&self) -> RangeInclusive<usize> {
        self.pos.clone()
    }
}
impl Seqresult for Blast {
    fn qseqid(&self) -> &str {
        &self.qseqid
    }

    fn sseq(&self) -> &str {
        &self.sseq
    }

    fn status(&self) -> &Status {
        &self.status
    }
    fn pos(&self) -> RangeInclusive<usize> {
        self.sstart..=self.send
    }
}
impl Newfasta {
    pub(crate) fn new(reader: &Ourfasta) -> io::Result<Self> {
        let split: Vec<&str> = reader.name.splitn(4, DELIMITERFASTA).collect();
        if split.len() != 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "{DELIMITERFASTA} in fasta lacking for {reader.name}",
            ));
        }
        let (a, b, c, d) = (split[0], split[1], split[2], split[3]);
        let b = match Status::try_from(b) {
            Ok(b) => b,
            Err(_) => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Invalid status for {reader.name}",
                ));
            }
        };
        let (c, d) = match (c.parse::<usize>(), d.parse::<usize>()) {
            (Ok(c), Ok(d)) => (c, d),
            _ => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Position are invalid for {reader.name}.",
                ));
            }
        };
        Ok(Self {
            qseqid: a.to_string(),
            sseq: reader.seq.to_string(),
            pos: (c..=d),
            status: b,
        })
    }
    #[allow(unused)]
    pub(crate) fn newfromblast(blast: &Blast) -> Self {
        Self {
            qseqid: blast.qseqid.clone(),
            sseq: blast.sseq.clone(),
            pos: blast.pos().clone(),
            status: blast.status.clone(),
        }
    }
    pub(crate) fn newfromblastowner(blast: Blast) -> Self {
        Self {
            pos: blast.pos().clone(),
            qseqid: blast.qseqid,
            sseq: blast.sseq,
            status: blast.status,
        }
    }
}
pub(crate) trait Seqresult {
    fn qseqid(&self) -> &str;
    fn sseq(&self) -> &str;
    fn status(&self) -> &Status;
    fn pos(&self) -> RangeInclusive<usize>;
}
impl TryFrom<&str> for Status {
    type Error = io::Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_lowercase().as_str() {
            "shorter" | "short" => Ok(Status::Shorter),
            "equal" | "eq" => Ok(Status::Equal),
            "new" => Ok(Self::New),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                "{value} is not valid, expect shorter, new or equal",
            )),
        }
    }
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) enum Status {
    #[default]
    Shorter,
    Equal,
    New,
}
impl Status {
    pub(crate) fn to_serialize(&self) -> Cow<'_, str> {
        Cow::Owned(match self {
            Status::Shorter => "shorter".to_string(),
            Status::Equal => "equal".to_string(),
            Status::New => "new".to_string(),
        })
    }
}
impl ToString for Status {
    fn to_string(&self) -> String {
        return format!("{}", self.to_serialize()).to_string();
    }
}
impl Newfasta {
    pub(crate) fn multipledeserialize(data: String) -> io::Result<Vec<Self>> {
        let mut elements = Vec::new();
        for elem in data.trim().trim_matches('>').split(">") {
            let data: Vec<&str> = elem.split([DELIMITERFASTA]).collect();
            if data.len() == 4
                && let Some((status, seq)) =
                    data.iter().last().map(|p| p.split_once("\n")).flatten()
            {
                let (posstart, posend, status) = match (
                    data[1].parse::<usize>(),
                    data[2].parse::<usize>(),
                    Status::try_from(status),
                ) {
                    (Ok(a), Ok(b), Ok(c)) => (a, b, c),
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Position or status is/are wrong for line {elem}",
                        ));
                    }
                };
                let entry = Newfasta {
                    qseqid: data[0].to_string(),
                    sseq: seq.to_string(),
                    pos: posstart..=posend,
                    status,
                };
                elements.push(entry);
            } else {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Invalid format for line {elem}",
                ));
            }
        }
        Ok(elements)
    }
}
impl std::fmt::Display for &dyn Seqresult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let pos = self.pos();
        write!(
            f,
            ">{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}\n{}",
            self.qseqid(),
            pos.start(),
            pos.end(),
            self.status().to_serialize(),
            self.sseq()
        )
    }
}
#[derive(Clone, Debug)]
pub(crate) struct Ourfasta {
    name: String,
    seq: String,
}
impl Serialize for Ourfasta {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let text = format!(">{}\n{}", self.name, self.seq);
        serializer.serialize_str(&text)
    }
}
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
#[must_use]
pub(crate) fn fastafilter(text: &str, find: &str, present: bool) -> String {
    let mut lines: Vec<&str> = Vec::with_capacity(text.lines().count());
    let mut keep = true;
    for line in text.lines() {
        let result = line.trim();
        if result.starts_with(">") {
            keep = if present {
                result
                    .to_ascii_lowercase()
                    .split('|')
                    .nth(2)
                    .is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
            } else {
                !result
                    .to_ascii_lowercase()
                    .split('|')
                    .nth(2)
                    .is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
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
                let mut file =
                    File::create(&tempfile).expect("Cannot create refseq file in temp dir");
                file.write_all(e.as_bytes())
                    .expect("Cannot create refseq file in temp dir.");
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
    blast.retain(|p| p.qseqid.contains(&format!("{}", locus)));
}
pub(crate) fn speciesandorphonfiltering(
    tempfile: &Path,
    releaseversion: String,
    species: &str,
    orphonfilter: bool,
) -> io::Result<PathBuf> {
    println!("Filtering based on species");
    let file = std::fs::read_to_string(&tempfile).unwrap();
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
        match (
            blastresult.pident,
            blastresult.gaps,
            blastresult.length,
            blastresult.qlen,
        ) {
            (100.0, 0, a, b) if a < b => blastresult.status = Status::Shorter,
            (100.0, 0, ..) => blastresult.status = Status::Equal,
            _ => blastresult.status = Status::New,
        }
    }
}
pub(crate) fn find_global_best_range(data: &[Blastmatch]) -> Option<(String, usize, usize)> {
    let mut groups: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    // Group by string, storing (pos, val) pairs
    for (s, pos, val) in data.iter().map(|f| (f.sseqid.clone(), f.sstart, f.send)) {
        groups.entry(s).or_insert_with(Vec::new).push((pos, val));
    }

    let mut global_best_name = String::new();
    let mut global_best_min = 0;
    let mut global_best_max = 0;
    let mut global_max_count = 0;

    // Process each group
    for (name, mut pairs) in groups {
        // Sort by the second usize (index 1)
        pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let mut max_count = 0;
        let mut best_min = 0;
        let mut best_max = 0;

        for left in 0..pairs.len() {
            let mut count = 0;
            let mut right = left;
            let right = loop {
                right += 1;
                if pairs.get(right).is_none() {
                    break right.saturating_sub(1);
                }
                // Move left to ensure window is valid
                if pairs[right].0.abs_diff(pairs[right.saturating_sub(1)].1) >= 1_000_000 {
                    break right;
                }
                count += 1;
            };
            if count > max_count {
                max_count = count;
                best_min = pairs[left].0;
                best_max = pairs[right].1;
            }
        }

        // Update global best
        if max_count > global_max_count {
            global_max_count = max_count;
            global_best_name = name;
            global_best_min = best_min;
            global_best_max = best_max;
        }
    }

    if global_max_count > 0 {
        Some((
            global_best_name.to_string(),
            global_best_min,
            global_best_max,
        ))
    } else {
        None
    }
}
pub(crate) fn locusposition<T>(
    subject: &Path,
    species: T,
    locus: &Locus,
) -> io::Result<(String, RangeInclusive<usize>)>
where
    T: AsRef<str>,
{
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
    let mut blast: Vec<Blast> = match blastcommand(reference.as_path(), subject) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
    locusfiltering(&locus, &mut blast);
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
        p.setstatus();
    });
    let mut statusvec: BTreeMap<(String, usize, usize), Blast> = BTreeMap::new();
    for elem in blast.into_iter() {
        if let Some((k, b)) = statusvec.clone().iter_mut().find(|((s, r1, r2), b)| {
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
    let range = find_global_best_range(&data)
        .and_then(|(n, pos1, pos2)| Some((n, pos1..=pos2)))
        .ok_or(io::Error::new(
            ErrorKind::InvalidInput,
            "No locus found after BLAST analysis",
        ));
    if let Ok((name, dat)) = &range {
        data.retain(|p| p.sseqid.as_str() == name && dat.contains(&p.sstart))
    }
    let text = data.iter().fold(String::new(), |mut acc, f| {
        acc.push_str(&format!("\n{}", f.to_string()));
        acc
    });
    let text = text.trim();
    //TODO: Separate
    fs::write("/tmp/genes.txt", text);
    range
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
            "-out",
            output,
            "-max_hsps",
            "5",
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send qlen length pident gaps sseq",
        ])
        .spawn()?;
    println!("Launching {} against {}", reference, subject);
    println!("BLAST has been launched with id {}", command.id());
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
pub(crate) fn filternewandsubmit(mut c: Vec<Blast>, realspecies: String) -> () {
    filterblast(&mut c);
    let result: Vec<Newfasta> = c.into_iter().map(Newfasta::newfromblastowner).collect();
    let dir = Path::new(&current_dir().unwrap_or(temp_dir())).join("archive");
    if dir.is_dir() {
        eprintln!("Archive directory exists, remove directory to submit if needed.");
        return;
    }
    if let Err(e) = fs::create_dir(&dir) {
        eprintln!("Cannot create archive directory, error is {e}");
        return;
    }
    if browseropening().is_err() {
        let _ = fs::remove_dir_all(dir);
        return;
    }
    let file = dir.join("sequence.fasta");
    let sequence = result.iter().fold(String::new(), |mut acc, f| {
        let f: &dyn Seqresult = f;
        acc.push_str(&format!("\n{}", f));
        acc
    });
    if let Err(e) = fs::write(file, sequence) {
        eprintln!("An error has occured while priting sequence: {e}.");
        return;
    }
    println!("BLAST results was added.");
    preparesubmission(&dir, realspecies, &REQUESTCLIENT);
    let _ = fs::remove_dir_all(dir);
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
    let archive = GzEncoder::new(file, Compression::default());
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
pub(crate) fn preparesubmission(
    path: &Path,
    species: String,
    release: &reqwest::blocking::Client,
) -> bool {
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
        if let Err(e) = submission(&token, species, archive, release) {
            eprintln!("An error has occured during submission: {e}. Please retry later.");
            return false;
        }
        println!("Your submission has been made. Thank you for submitting your sequences to IMGT.");
    } else {
        println!("Exiting");
        return false;
    }
    true
}
pub(crate) fn submission(
    token: &str,
    species: String,
    archive: NamedTempFile,
    release: &reqwest::blocking::Client,
) -> io::Result<()> {
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
    match release
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
            StatusCode::SERVICE_UNAVAILABLE | StatusCode::INTERNAL_SERVER_ERROR => {
                Err(io::Error::new(
                    io::ErrorKind::ResourceBusy,
                    "The server is unavailable. Please retry a submission.",
                ))
            }
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

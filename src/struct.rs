use bio::io::fasta;
use itertools::Itertools;
use serde::ser::SerializeTupleStruct;
use serde_with::{DefaultOnError, DisplayFromStr, serde_as};
/*
This software allows the analysis of BAM files to identify the confidence on a locus (specifically IG and TR) as well as allele confidence.
It was created and used by IMGT Team (https://www.imgt.org).
Available under EUPL license
Made by: Guilhem Zeitoun
*/
use crate::submissions::{
    DELIMITERFASTA, getallelefromblast, getchromosomefromblast, getpositionfromblast,
};
use crate::{MATCHREADS, SOFTCLIPRATIO, locusisokay};
use clap::{Parser, crate_authors};
use serde::{Deserialize, Serialize, de};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io::ErrorKind;
use std::ops::{Not, RangeInclusive};
use std::str::FromStr;
use std::{borrow::Cow, fmt::Display, fs::File, hash::Hash, io, path::PathBuf};
use strum_macros::EnumIter;
#[derive(Parser, Debug)]
#[clap(
    author = crate_authors!("\n"),
    before_help = "This script analyzes BAM files coming from reads assembled on an assembly to assess IG/TR loci quality.\nYou can also submit new sequences to IMGT.",
    after_help = format!("This code was made by and for IMGT (the international ImMunoGeneTics information system) under {} license.",env!("CARGO_PKG_LICENSE")),
    arg_required_else_help=true,
    display_name="IMGT/StatAssembly",
    help_template = "\
    {name} {version}
    Authors: {author-section}
    {before-help}
    About: {about-with-newline}
    {usage-heading} {usage}

    {all-args}{after-help}
    "
)]
#[command(version, author, about, long_about = None)]
pub(crate) struct Args {
    /// Input file (BAM-indexed file)
    #[arg(short, long)]
    pub(crate) file: PathBuf,
    /// Index file if not default
    #[arg(short, long)]
    pub(crate) index: Option<PathBuf>,
    ///CSV containing locus infos. See example file for blueprint. If unset, will try all loci on both haplotypes
    #[arg(short, long)]
    pub(crate) locuspos: Option<PathBuf>,
    /// Minimal number of reads (included) to declare a break in coverage
    #[arg(short, long, default_value_t = 3)]
    pub(crate) breaks: u32,
    /// Coverage to calculate on CSV
    #[arg(short, long, default_value_t = 10, value_parser=greater_than_0)]
    pub(crate) coverage: u32,
    /// Minimum number of match reads (included) for warning positions
    #[arg(long, default_value_t = 10)]
    pub(crate) minreads: u32,
    /// Percent warning position for mismatch reads (included)
    #[arg(long, default_value_t = 80, value_parser=less_than_100)]
    pub(crate) percentwarning: u8,
    /// Percent alerting position for mismatch reads (included)
    #[arg(long, default_value_t = 60, value_parser=less_than_100)]
    pub(crate) percentalerting: u8,
    /// Force cigar even if no =. Some functionalities would be disabled
    #[arg(long)]
    pub(crate) force: bool,
    /// If the BAM file is truncated, length of the overall extracted sequence (default: 0 meaning full length). This analysis takes between 10 and 15 minutes.
    #[arg(long, conflicts_with = "meancoverage", default_value_t = 0)]
    pub(crate) extractedlength: u64,
    /// The mean coverage is already calculated, else will be calculated at startup (10-15 minutes) if not stored already.
    #[arg(long, value_parser=greater_than_0_64)]
    pub(crate) meancoverage: Option<u64>,
    /// Huge region (more than 10 Mb)
    #[arg(long)]
    pub(crate) hugeregion: bool,
    /// Number of threads to decrypt bgzf files (0 for number of threads up to 12)
    #[arg(long, default_value_t = 0)]
    pub(crate) threads: usize,
    /// Only strand-specific alignments to reference
    #[arg(long)]
    pub(crate) forward: bool,
    /// Assembly fasta file
    #[arg(long, short = 'a')]
    pub(crate) assembly: Option<PathBuf>,
    /// Assembly Index file if not default
    #[arg(long, short = 'j')]
    pub(crate) assemblyindex: Option<PathBuf>,
    /// Query full quality PHRED score (script will be longer to execute)
    #[arg(long)]
    pub(crate) fullquality: bool,
    /// Calculate total reads mismatch
    #[arg(long, short = 't')]
    pub(crate) totalread: bool,
    /// Size of legend axis (default 16)
    #[arg(long, default_value_t = 16, value_parser=greater_than_0)]
    pub(crate) fontlegendsize: u32,
    /// No legend on graphs
    #[arg(long)]
    pub(crate) nolegend: bool,
    /// Get supplementary and secondary alignments on gene graphs as well
    #[arg(long)]
    pub(crate) allreads: bool,
    /// Save as SVG images (create big images)
    #[arg(long)]
    pub(crate) svg: bool,
    ///Species name (must match NCBI taxonomy & for folder creation)
    #[arg(short, long)]
    pub(crate) species: String,
    ///Gene location (csv file). See example file for blueprint.
    #[arg(short, long)]
    pub(crate) geneloc: Option<PathBuf>,
    ///Output directory (created or overwritten)
    #[arg(short, long)]
    pub(crate) outdir: PathBuf,
    ///Output light BAM for submission
    #[arg(short = 'z', long)]
    pub(crate) outlightbam: Option<PathBuf>,
    ///Haploid status (only if locuspos is not set)
    #[arg(long, conflicts_with = "locuspos")]
    pub(crate) haploid: bool,
    /// Do not submit to IMGT
    #[arg(long)]
    pub(crate) nosubmit: bool,
    /// Erase cache files
    #[arg(long)]
    pub(crate) cacheerase: bool,
    /// Automatic token to submit
    #[arg(long, value_parser=checktoken, conflicts_with = "nosubmit")]
    pub(crate) mytoken: Option<String>,
}
fn checktoken(s: &str) -> Result<String, String> {
    let s = s.trim();
    if !s.is_ascii() || s.len() != 24 || !s.chars().any(|p| p.is_ascii_alphanumeric()) {
        Err(String::from(
            "Invalid token, please remove it to process or verify its value.",
        ))
    } else {
        Ok(String::from(s))
    }
}
pub(crate) fn less_than_100(s: &str) -> Result<u8, String> {
    match s.parse::<u8>() {
        Ok(s) if (0..=100).contains(&s) => Ok(s),
        _ => Err(String::from("Bad number, must be greater than 0.")),
    }
}
pub(crate) fn greater_than_0(s: &str) -> Result<u32, String> {
    match s.parse::<u32>() {
        Ok(s) if s != u32::MIN => Ok(s),
        _ => Err(String::from("Bad number, must be greater than 0.")),
    }
}
pub(crate) fn greater_than_0_64(s: &str) -> Result<u64, String> {
    match s.parse::<u64>() {
        Ok(s) if s != u64::MIN => Ok(s),
        _ => Err(String::from("Bad number, must be greater than 0.")),
    }
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize, Hash, EnumIter)]
#[allow(clippy::upper_case_acronyms)]
pub(crate) enum Locus {
    IGH,
    IGK,
    IGL,
    TRA,
    TRB,
    TRG,
}
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Hash)]
pub(crate) enum OkStatus {
    #[default]
    Unknown,
    Accepted,
    Rejected,
}
impl Not for OkStatus {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            OkStatus::Unknown => Self::Unknown,
            OkStatus::Accepted => Self::Rejected,
            OkStatus::Rejected => Self::Accepted,
        }
    }
}
impl Display for OkStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OkStatus::Unknown => write!(f, "Unknown"),
            OkStatus::Accepted => write!(f, "Accepted"),
            OkStatus::Rejected => write!(f, "Rejected"),
        }
    }
}
impl OkStatus {
    pub(crate) fn isvalid(&self) -> bool {
        *self == OkStatus::Accepted
    }
    #[allow(dead_code)]
    pub(crate) fn isinvalid(&self) -> bool {
        *self == OkStatus::Rejected
    }
    #[allow(dead_code)]
    pub(crate) fn isunknown(&self) -> bool {
        *self == OkStatus::Unknown
    }
}
impl TryFrom<String> for Locus {
    type Error = io::Error;
    fn try_from(value: String) -> Result<Self, Self::Error> {
        match value.trim().to_ascii_uppercase().as_str() {
            "IGH" => Ok(Locus::IGH),
            "IGK" => Ok(Locus::IGK),
            "IGL" => Ok(Locus::IGL),
            "TRA" | "TRD" => Ok(Locus::TRA),
            "TRB" => Ok(Locus::TRB),
            "TRG" => Ok(Locus::TRG),
            _ => Err(io::Error::new(
                ErrorKind::InvalidInput,
                format!("{} is not a valid locus.", value),
            )),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(clippy::upper_case_acronyms)]
pub enum Alertpos {
    Valid,
    Warning,
    Suspicious,
}
pub(crate) trait Alerting {
    ///Is a warning position
    fn iswarning(&self) -> bool;
    ///Is a suspicious position
    fn issuspicious(&self) -> bool;
    #[allow(dead_code)]
    /// Is a non-suspicious nor warning position
    fn isvalid(&self) -> bool;
}
impl Alertpos {
    fn new(record: &Posread) -> Self {
        let percent = if record.total > 0 {
            record
                .r#match
                .saturating_mul(100)
                .saturating_div(record.total)
        } else {
            0
        };
        match (
            record.percentalerting,
            record.percentwarning,
            record.minreads,
        ) {
            (d, ..) if percent <= d.into() => Alertpos::Suspicious,
            (_, d, e)
                if record.r#match <= e.try_into().unwrap_or(usize::MAX) || percent <= d.into() =>
            {
                Alertpos::Warning
            }
            _ => Alertpos::Valid,
        }
    }
}
impl Alerting for Alertpos {
    ///Is a warning position
    fn iswarning(&self) -> bool {
        matches!(self, Alertpos::Warning)
    }
    ///Is a suspicious position
    fn issuspicious(&self) -> bool {
        matches!(self, Alertpos::Suspicious)
    }
    /// Is a non-suspicious nor warning position
    #[allow(dead_code)]
    fn isvalid(&self) -> bool {
        matches!(self, Alertpos::Valid)
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]

pub(crate) struct Position {
    zbased: bool,
    position: i64,
}
impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Position {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.getzbasedpos().cmp(&other.getzbasedpos())
    }
}
///Serialize as a 1-based position
impl Serialize for Position {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_i64(self.getobasedpos())
    }
}
impl FromStr for Position {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Position::new(false, s.parse()?))
    }
}
impl<'de> Deserialize<'de> for Position {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        let s: &str = de::Deserialize::deserialize(deserializer)?;

        match s.parse::<i64>() {
            Ok(pos) => Ok(Position::new(false, pos)),
            Err(_) => Err(de::Error::invalid_type(
                de::Unexpected::Str(s),
                &"expected i64",
            )),
        }
    }
}
impl Default for Position {
    fn default() -> Self {
        Self {
            zbased: true,
            position: 0,
        }
    }
}
impl Position {
    pub(crate) fn new(zbased: bool, position: i64) -> Self {
        Position { zbased, position }
    }
    pub(crate) fn length(&self, other: &Self) -> i64 {
        let min = std::cmp::min(self.getzbasedpos(), other.getzbasedpos());
        let max = std::cmp::max(self.getzbasedpos(), other.getzbasedpos());
        max.checked_sub(min)
            .and_then(|f| f.checked_add(1))
            .unwrap_or(0) //Calculate the length
    }
    pub(crate) fn getzbasedpos(&self) -> i64 {
        if self.zbased {
            self.position
        } else {
            self.position.saturating_sub(1)
        }
    }
    pub(crate) fn getobasedpos(&self) -> i64 {
        if self.zbased {
            self.position.saturating_add(1)
        } else {
            self.position
        }
    }
    #[allow(dead_code)]
    pub(crate) fn iszbased(&self) -> bool {
        self.zbased
    }
}
#[derive(Clone, Debug, PartialEq, Copy)]
pub(crate) struct Posread {
    pub(crate) r#match: usize,
    pub(crate) indel: usize,
    pub(crate) total: usize,
    pub(crate) minreads: u32,
    pub(crate) percentwarning: u8,
    pub(crate) percentalerting: u8,
    pub(crate) softclips: f32,
}
impl Alerting for Posread {
    fn iswarning(&self) -> bool {
        self.getstate().iswarning()
    }

    fn issuspicious(&self) -> bool {
        self.getstate().issuspicious()
    }

    fn isvalid(&self) -> bool {
        self.getstate().isvalid()
    }
}
#[derive(Clone, Debug, Deserialize)]
pub(crate) struct Blast {
    pub(crate) qseqid: String,
    pub(crate) sseqid: String,
    #[allow(dead_code)]
    pub(crate) qstart: usize,
    #[allow(dead_code)]
    pub(crate) qend: usize,
    pub(crate) sstart: usize,
    pub(crate) send: usize,
    pub(crate) qlen: usize,
    pub(crate) length: usize,
    pub(crate) pident: f32,
    #[allow(dead_code)]
    pub(crate) gaps: usize,
    pub(crate) sseq: String,
    #[serde(skip_deserializing)]
    pub(crate) complement: bool,
    #[serde(skip_deserializing)]
    pub(crate) status: Status,
}
#[repr(C)]
#[derive(Clone, Debug, Serialize)]
pub(crate) struct Blastmatch {
    pub(crate) qseqid: String,
    pub(crate) sseqid: String,
    pub(crate) sseq: String,
    pub(crate) sstart: usize,
    pub(crate) send: usize,
    pub(crate) complement: Strand,
    pub(crate) status: Status,
    pub(crate) identity: f32,
}
impl PartialOrd for Blastmatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Blastmatch {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.sseqid.cmp(&other.sseqid) {
            std::cmp::Ordering::Equal => match self.sstart.cmp(&other.sstart) {
                std::cmp::Ordering::Equal => self.send.cmp(&other.send).reverse(),
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
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        qseqid: String,
        sseqid: String,
        sseq: String,
        sstart: usize,
        send: usize,
        complement: Strand,
        status: Status,
        identity: f32,
    ) -> Self {
        Self {
            qseqid,
            sseqid,
            sseq,
            sstart,
            send,
            complement,
            status,
            identity,
        }
    }
}
impl Blast {
    #[allow(dead_code)]
    pub(crate) fn getstatus(&self) -> &Status {
        &self.status
    }
    pub(crate) fn setstatus(&mut self) {
        match (self.pident, self.qlen, self.length) {
            (100.0, a, b) if a == b => self.status = Status::Equal,
            (100.0, a, b) if b < a => self.status = Status::Shorter,
            _ => self.status = Status::New,
        }
    }
}
#[allow(clippy::from_over_into)]
//Blastmatch is for matches, should not be converted
impl Into<Blastmatch> for Blast {
    fn into(mut self) -> Blastmatch {
        self.setstatus();
        Blastmatch::new(
            self.qseqid,
            self.sseqid,
            self.sseq,
            self.sstart,
            self.send,
            if self.complement {
                Strand::Minus
            } else {
                Strand::Plus
            },
            self.status,
            self.pident,
        )
    }
}
impl Display for Blastmatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!(
            ">{}{DELIMITERFASTA}{}:{}-{}-{}{DELIMITERFASTA}{}\n{}",
            self.qseqid,
            self.sseqid,
            self.sstart,
            self.send,
            self.complement,
            self.status,
            self.sseq
        ))
    }
}
impl PartialEq for Blast {
    fn eq(&self, other: &Self) -> bool {
        self.qseqid == other.qseqid && self.sstart == other.sstart && self.send == other.send
    }
}
impl Eq for Blast {}
pub trait Blastcalc {
    #[allow(dead_code)]
    fn getseq(&self) -> Cow<'_, str>;
    fn getstatus(&self) -> &Status;
    /// Only return new alleles
    fn onlynewalleles(&self) -> bool {
        self.getstatus() == &Status::New
    }
    /// get subject
    fn getsubject(&self) -> &str;
    fn getpos(&self) -> (usize, usize);
    fn getposrange(&self) -> RangeInclusive<usize> {
        let r = self.getpos();
        r.0..=r.1
    }
    #[allow(dead_code)]
    fn getstrand(&self) -> Strand {
        match self.getpos() {
            (a, b) if a > b => Strand::Minus,
            _ => Strand::Plus,
        }
    }
    fn getquery(&self) -> &str;
    fn getallelename(&self) -> Option<String> {
        getallelefromblast(self.getquery())
    }
    fn getposition(&self) -> Option<(Position, Position, Strand)> {
        getpositionfromblast(self.getsubject())
    }
    fn getchromosomefromblast(&self) -> Option<String> {
        getchromosomefromblast(self.getsubject())
    }
    #[allow(dead_code)]
    fn getlocusname(&self) -> Option<Locus> {
        self.getallelename()
            .and_then(|f| Locus::try_from(f.split_at(3).0.to_string()).ok())
    }
    fn getidentity(&self) -> f32;
}
impl Blastcalc for Blast {
    fn getseq(&self) -> Cow<'_, str> {
        Cow::Borrowed(&self.sseq)
    }
    fn getsubject(&self) -> &str {
        &self.sseqid
    }
    fn getstatus(&self) -> &Status {
        &self.status
    }
    fn getpos(&self) -> (usize, usize) {
        (self.sstart, self.send)
    }
    fn getquery(&self) -> &str {
        &self.qseqid
    }
    fn getidentity(&self) -> f32 {
        self.pident
    }
}
impl Blastcalc for Blastmatch {
    /// Get sequence
    fn getseq(&self) -> Cow<'_, str> {
        Cow::Borrowed(&self.sseq)
    }
    fn getsubject(&self) -> &str {
        &self.sseqid
    }
    fn getstatus(&self) -> &Status {
        &self.status
    }
    fn getpos(&self) -> (usize, usize) {
        (self.sstart, self.send)
    }
    /// Query id
    fn getquery(&self) -> &str {
        &self.qseqid
    }
    /// Get identity
    fn getidentity(&self) -> f32 {
        self.identity
    }
}
#[derive(Clone, Debug)]
pub(crate) struct Newfasta {
    qseqid: String,
    sseqid: String,
    sseq: String,
    pos: RangeInclusive<usize>,
    status: Status,
    identity: f32,
}
impl Blastcalc for Newfasta {
    /// Get sequence
    fn getseq(&self) -> Cow<'_, str> {
        Cow::Borrowed(&self.sseq)
    }
    fn getsubject(&self) -> &str {
        &self.sseqid
    }
    fn getstatus(&self) -> &Status {
        &self.status
    }
    fn getpos(&self) -> (usize, usize) {
        (*self.pos.start(), *self.pos.end())
    }
    /// Query id
    fn getquery(&self) -> &str {
        &self.qseqid
    }
    /// Get identity
    fn getidentity(&self) -> f32 {
        self.identity
    }
}
impl Seqresult for Newfasta {}
impl From<&dyn Blastcalc> for Newfasta {
    fn from(value: &dyn Blastcalc) -> Self {
        Self {
            qseqid: value.getquery().to_string(),
            sseqid: value.getsubject().to_string(),
            sseq: value.getseq().to_string(),
            pos: value.getposrange(),
            status: value.getstatus().clone(),
            identity: value.getidentity(),
        }
    }
}
pub(crate) trait Seqresult: Blastcalc {}
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
impl Serialize for Status {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_serialize())
    }
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
impl Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{}", self.to_serialize()))
    }
}
/* impl Newfasta {
    pub(crate) fn multipledeserialize(data: String) -> io::Result<Vec<Self>> {
        let mut elements = Vec::new();
        for elem in data.trim().trim_matches('>').split(">") {
            let elem =
                Newfasta::from_str(elem).map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
            elements.push(elem);
        }
        Ok(elements)
    }
} */
impl std::fmt::Display for &dyn Seqresult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (start, end) = self.getpos();
        write!(
            f,
            ">{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}{DELIMITERFASTA}{}\n{}",
            self.getquery(),
            self.getsubject(),
            start,
            end,
            self.getstatus().to_serialize(),
            self.getidentity(),
            self.getseq(),
        )
    }
}
impl FromStr for Newfasta {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = if let Some(a) = s.trim().strip_prefix('>') {
            a
        } else {
            return Err("Does not start with fasta format".to_string());
        };
        let (data, seq) = match s.split_once('\n') {
            Some((a, b)) => (a, b),
            None => return Err("No sequence found".to_string()),
        };
        let seq = regex::Regex::new(r"\s+")
            .unwrap_or_else(|_| unreachable!("Regex is ok"))
            .replace_all(seq, "")
            .to_ascii_lowercase();
        if seq
            .chars()
            .any(|p| !p.is_ascii_alphabetic() || !['n', 'c', 'a', 't'].contains(&p))
        {
            return Err("Sequence has invalid nucleotides. Only atcgn are accepted.".to_string());
        }
        let split: Vec<&str> = data.splitn(6, DELIMITERFASTA).collect();
        if split.len() < 6 {
            return Err("Invalid header".to_string());
        }
        let [query, subject, start, end, status, identity] = split[..6] else {
            return Err("Invalid header".to_string());
        };
        let status = match Status::try_from(status) {
            Ok(b) => b,
            Err(_) => return Err("Status is invalid".to_string()),
        };
        let (start, end) = match (start.parse::<usize>(), end.parse::<usize>()) {
            (Ok(a), Ok(b)) if b >= a => (a, b),
            _ => return Err("Position are invalid".to_string()),
        };
        Ok(Newfasta {
            qseqid: query.to_string(),
            sseqid: subject.to_string(),
            sseq: seq,
            pos: (start..=end),
            status,
            identity: identity.parse::<f32>().unwrap_or(0f32),
        })
    }
}
#[derive(Debug)]
pub(crate) struct MyError(String);

impl std::fmt::Display for MyError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl std::error::Error for MyError {}
impl Posread {
    #[allow(dead_code)]
    pub(crate) fn new(
        r#match: usize,
        indel: usize,
        total: usize,
        softclips: f32,
        args: &Args,
    ) -> Result<Self, MyError> {
        if r#match + indel > total {
            return Err(MyError(String::from("Invalid total")));
        }
        Ok(Self {
            r#match,
            indel,
            total,
            minreads: args.minreads,
            percentwarning: args.percentwarning,
            percentalerting: args.percentalerting,
            softclips,
        })
    }
    ///Get the state of the position
    fn getstate(&self) -> Alertpos {
        Alertpos::new(self)
    }
    pub(crate) fn gettotal(&self) -> usize {
        self.total
    }
    pub(crate) fn getmismatchcount(&self) -> usize {
        self.getindel().saturating_sub(self.getmatch())
    }
    pub(crate) fn getindelcount(&self) -> usize {
        self.gettotal().saturating_sub(self.getindel())
    }
    pub(crate) fn addtotal(&mut self, count: usize) {
        self.total += count
    }
    pub(crate) fn getmatch(&self) -> usize {
        self.r#match
    }
    pub(crate) fn addmatch(&mut self, count: usize) {
        self.r#match += count
    }
    pub(crate) fn getindel(&self) -> usize {
        self.indel
    }
    pub(crate) fn addindel(&mut self, count: usize) {
        self.indel += count
    }
}
impl Display for Locus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Locus::IGH => write!(f, "IGH"),
            Locus::IGK => write!(f, "IGK"),
            Locus::IGL => write!(f, "IGL"),
            Locus::TRA => write!(f, "TRA_TRD"),
            Locus::TRB => write!(f, "TRB"),
            Locus::TRG => write!(f, "TRG"),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Default, PartialOrd, Ord, Hash, EnumIter)]
pub(crate) enum Haplotype {
    #[default]
    Primary,
    Alternate,
}
impl FromStr for Haplotype {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "primary" | "pri" | "p" => Ok(Haplotype::Primary),
            "alternate" | "alt" | "a" => Ok(Haplotype::Alternate),
            _ => Err("Invalid haplotype".to_string()),
        }
    }
}
impl<'de> Deserialize<'de> for Haplotype {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        let s: &str = de::Deserialize::deserialize(deserializer)?;
        Haplotype::from_str(s).map_err(|_| {
            de::Error::unknown_variant(s, &["primary or pri or p", "alternate or alt or a"])
        })
    }
}
impl Serialize for Haplotype {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            Haplotype::Primary => serializer.serialize_str("Primary"),
            Haplotype::Alternate => serializer.serialize_str("Alternate"),
        }
    }
}
impl Display for Haplotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Primary => write!(f, "Primary"),
            Self::Alternate => write!(f, "Alternate"),
        }
    }
}
impl Haplotype {
    pub(crate) fn isprimary(&self) -> bool {
        self == &Haplotype::Primary
    }
}
#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct GeneInfos {
    pub(crate) gene: String,
    pub(crate) chromosome: String,
    pub(crate) strand: Strand,
    pub(crate) start: Position,
    pub(crate) end: Position,
    #[serde(skip_deserializing)]
    pub(crate) status: OkStatus,
}
impl GeneInfos {
    pub(crate) fn new(
        gene: String,
        chromosome: String,
        strand: Strand,
        start: Position,
        end: Position,
    ) -> Self {
        GeneInfos {
            gene,
            chromosome,
            strand,
            start,
            end,
            status: OkStatus::default(),
        }
    }
    pub(crate) fn addtosequence<T>(&self, seq: T, fasta: &mut fasta::Writer<File>) -> io::Result<()>
    where
        T: AsRef<[u8]>,
    {
        fasta.write(
            &format!(
                "GENE|{}|{}|{}|{}..{}{}",
                self.getgene(),
                "species",
                self.getchromosome(),
                self.getstart().getobasedpos(),
                self.getend().getobasedpos(),
                if self.getstrand().isrev() { "/rc" } else { "" }
            ),
            None,
            seq.as_ref(),
        )
    }
    pub(crate) fn extractsequence(
        &self,
        fasta: &mut fasta::IndexedReader<File>,
    ) -> io::Result<String> {
        let (startpos, endpos) = match (
            self.getstart().getzbasedpos().try_into(),
            self.getend().getobasedpos().try_into(),
        ) {
            (Ok(a), Ok(b)) => (a, b),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid start and end position",
                ));
            }
        };
        fasta.fetch(self.getchromosome(), startpos, endpos)?;
        let length = endpos.saturating_sub(startpos).saturating_add(1);
        let mut cap = Vec::with_capacity(length.try_into().unwrap_or(0));
        fasta.read(&mut cap)?;
        if self.strand.isrev() {
            cap = bio::alphabets::dna::revcomp(&cap);
        }
        String::from_utf8(cap).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
impl GenesList for GeneInfos {
    fn getgene(&self) -> &str {
        &self.gene
    }

    fn getchromosome(&self) -> &str {
        &self.chromosome
    }

    fn getstrand(&self) -> &Strand {
        &self.strand
    }

    fn getstart(&self) -> &Position {
        &self.start
    }

    fn getend(&self) -> &Position {
        &self.end
    }

    fn getstatus(&self) -> &OkStatus {
        &self.status
    }

    fn alterstatus(&mut self, status: OkStatus) {
        self.status = status;
    }
}
impl LocusInfos {
    pub(crate) fn setstatus(&mut self, mean: u64, pos: &BTreeMap<Position, HashMapinfo>) {
        if locusisokay(mean, &pos.values().collect_vec()) {
            self.status = OkStatus::Accepted
        } else {
            self.status = OkStatus::Rejected
        }
    }
    pub(crate) fn locusinposition(
        &self,
        newstart: &Position,
        newend: &Position,
        complement: &Strand,
        matchstart: &Position,
        matchend: &Position,
        matchcomplement: &Strand,
    ) -> Option<(Position, Position, Strand)> {
        if newend.length(newstart) < matchend.length(matchstart) {
            return None;
        }
        let (start, end) = if self.complement.isrev() {
            let end = newend.getobasedpos() - std::cmp::min(matchstart, matchend).getobasedpos();
            let start = newend.getobasedpos() - std::cmp::max(matchstart, matchend).getobasedpos();
            (start, end)
        } else {
            let start =
                newstart.getobasedpos() + std::cmp::min(matchstart, matchend).getobasedpos();
            let end = newstart.getobasedpos() + std::cmp::max(matchstart, matchend).getobasedpos();
            (start, end)
        };
        let complement = matchcomplement.isrev() ^ complement.isrev();
        Some((
            Position::new(false, start),
            Position::new(false, end),
            if complement {
                Strand::Minus
            } else {
                Strand::Plus
            },
        ))
    }
    pub(crate) fn positioninlocus(
        &self,
        start: &Position,
        end: &Position,
        complement: &Strand,
    ) -> Option<(Position, Position, Strand)> {
        if start > end || end.length(start) > self.getlength() {
            return None;
        }
        if self.complement.isrev() {
            let newend = self
                .end
                .getobasedpos()
                .saturating_sub(std::cmp::min(end.getobasedpos(), start.getobasedpos()));
            let newstart = self
                .end
                .getobasedpos()
                .saturating_sub(std::cmp::max(end.getobasedpos(), start.getobasedpos()));
            let complement = !complement.clone();
            Some((
                Position::new(false, newstart),
                Position::new(false, newend),
                complement,
            ))
        } else {
            let newstart = self
                .start
                .getobasedpos()
                .saturating_add(start.getobasedpos());
            let newend = self.start.getobasedpos().saturating_add(end.getobasedpos());
            let complement = complement.clone();
            Some((
                Position::new(false, newstart),
                Position::new(false, newend),
                complement,
            ))
        }
    }
    pub(crate) fn extractsequence(
        &self,
        fasta: &mut fasta::IndexedReader<File>,
    ) -> io::Result<String> {
        let fake = GeneInfos {
            gene: "FAKE".to_string(),
            chromosome: self.contig.clone(),
            start: self.start,
            end: self.end,
            strand: self.complement.clone(),
            status: OkStatus::default(),
        };
        fake.extractsequence(fasta)
    }
}
impl PartialEq for GeneInfos {
    fn eq(&self, other: &Self) -> bool {
        self.gene == other.gene
    }
}
impl std::hash::Hash for GeneInfos {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.gene.hash(state);
    }
}
impl PartialOrd for GeneInfos {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for GeneInfos {
    //Filter by chromosome, then by start ASC and end DESC
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.chromosome.cmp(&other.chromosome) {
            std::cmp::Ordering::Equal => match self.start.cmp(&other.start) {
                std::cmp::Ordering::Equal => self.end.cmp(&other.end).reverse(),
                ord => ord,
            },
            e => e,
        }
    }
}
impl Eq for GeneInfos {}
#[derive(Clone, Debug, Default, Serialize, PartialEq, Eq, Hash)]
pub(crate) enum Strand {
    #[default]
    Plus,
    Minus,
}
impl Strand {
    pub(crate) fn isrev(&self) -> bool {
        self == &Strand::Minus
    }
    #[allow(unused)]
    pub(crate) fn isfwd(&self) -> bool {
        self == &Strand::Plus
    }
}
impl Not for Strand {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }
}
impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        let s: &str = de::Deserialize::deserialize(deserializer)?;

        match s.to_lowercase().as_str() {
            "1" | "-" | "minus" => Ok(Strand::Minus),
            "0" | "+" | "plus" => Ok(Strand::Plus),
            _ => Err(de::Error::unknown_variant(
                s,
                &["1 or - or minus", "0 or + or plus"],
            )),
        }
    }
}
impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Plus => write!(f, "FWD"),
            Self::Minus => write!(f, "REV"),
        }
    }
}
pub trait GenesList {
    fn getgene(&self) -> &str;
    fn getchromosome(&self) -> &str;
    fn getstrand(&self) -> &Strand;
    fn getstart(&self) -> &Position;
    fn getend(&self) -> &Position;
    #[allow(dead_code)]
    fn getstatus(&self) -> &OkStatus;
    fn alterstatus(&mut self, status: OkStatus);
    fn setstatus(&mut self, reads100m: usize, hash: &BTreeMap<Position, Posread>) {
        if hash
            .iter()
            .all(|(_, f)| f.isvalid() && f.softclips < SOFTCLIPRATIO)
            && reads100m >= MATCHREADS
        {
            self.alterstatus(OkStatus::Accepted);
        } else {
            self.alterstatus(OkStatus::Rejected);
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub(crate) struct GeneInfosFinish {
    pub(crate) gene: String,
    pub(crate) chromosome: String,
    pub(crate) strand: Strand,
    pub(crate) start: Position,
    pub(crate) end: Position,
    pub(crate) length: i64,
    pub(crate) readscoverage: f32,
    pub(crate) reads: usize,
    pub(crate) matchpos: String,
    pub(crate) readsfull: usize,
    pub(crate) reads100: usize,
    pub(crate) reads100m: usize,
    pub(crate) coveragex: usize,
    pub(crate) status: OkStatus,
}
impl GenesList for GeneInfosFinish {
    fn getgene(&self) -> &str {
        &self.gene
    }

    fn getchromosome(&self) -> &str {
        &self.chromosome
    }

    fn getstrand(&self) -> &Strand {
        &self.strand
    }

    fn getstart(&self) -> &Position {
        &self.start
    }

    fn getend(&self) -> &Position {
        &self.end
    }

    fn getstatus(&self) -> &OkStatus {
        &self.status
    }

    fn alterstatus(&mut self, status: OkStatus) {
        self.status = status;
    }
}
impl Ord for dyn GenesList {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.getchromosome().cmp(other.getchromosome()) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        match self.getstart().cmp(other.getstart()) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        match self.getend().cmp(other.getend()) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        std::cmp::Ordering::Equal
    }
}
impl PartialOrd for dyn GenesList {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for dyn GenesList {
    fn eq(&self, other: &Self) -> bool {
        self.getgene() == other.getgene()
    }
}
impl Eq for dyn GenesList {}
impl From<GeneInfosFinish> for GeneInfos {
    fn from(value: GeneInfosFinish) -> Self {
        GeneInfos {
            gene: value.gene,
            chromosome: value.chromosome,
            strand: value.strand,
            start: value.start,
            end: value.end,
            status: value.status,
        }
    }
}
impl GeneInfosFinish {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        gene: GeneInfos,
        reads: usize,
        readsfull: usize,
        matchpos: Option<String>,
        reads100: usize,
        reads100m: usize,
        readscoverage: f32,
        coveragex: usize,
    ) -> Self {
        GeneInfosFinish {
            gene: gene.gene,
            chromosome: gene.chromosome,
            strand: gene.strand,
            length: gene.end.length(&gene.start),
            start: gene.start,
            end: gene.end,
            reads,
            matchpos: matchpos.unwrap_or(String::from("N/A")),
            readsfull,
            reads100,
            reads100m,
            readscoverage,
            coveragex,
            status: gene.status,
        }
    }
    pub(crate) fn make_default(gene: GeneInfos) -> Self {
        Self::new(gene, 0, 0, None, 0, 0, 0.0, 0)
    }
}
impl Ord for GeneInfosFinish {
    fn cmp(&self, other: &Self) -> Ordering {
        let a: &dyn GenesList = self;
        let b: &dyn GenesList = other;
        a.cmp(b)
    }
}
impl PartialOrd for GeneInfosFinish {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for GeneInfosFinish {
    fn eq(&self, other: &Self) -> bool {
        let a: &dyn GenesList = self;
        let b: &dyn GenesList = other;
        a.eq(b)
    }
}
impl Eq for GeneInfosFinish {}
#[serde_as]
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Hash)]
pub(crate) struct FakeLocusinfo {
    pub(crate) locus: Locus,
    #[serde_as(as = "DefaultOnError<Option<DisplayFromStr>>")]
    #[serde(default)]
    pub(crate) haplotype: Option<Haplotype>,
    #[serde_as(as = "DefaultOnError<Option<DisplayFromStr>>")]
    #[serde(default)]
    pub(crate) contig: Option<String>,
    //#[serde_as(as = "DefaultOnError<Option<DisplayFromStr>>")]
    #[serde(default)]
    pub(crate) start: Option<Position>,
    //#[serde_as(as = "DefaultOnError<Option<DisplayFromStr>>")]
    #[serde(default)]
    pub(crate) end: Option<Position>,
    #[serde(skip)]
    pub(crate) complement: Strand,
}
impl From<LocusInfos> for FakeLocusinfo {
    fn from(value: LocusInfos) -> Self {
        FakeLocusinfo {
            locus: value.locus,
            haplotype: Some(value.haplotype),
            contig: Some(value.contig),
            start: Some(value.start),
            end: Some(value.end),
            complement: value.complement,
        }
    }
}
impl FakeLocusinfo {
    pub(crate) fn new(
        locus: Locus,
        haplotype: Option<Haplotype>,
        contig: Option<String>,
        start: Option<Position>,
        end: Option<Position>,
        complement: Option<Strand>,
    ) -> Self {
        FakeLocusinfo {
            locus,
            haplotype,
            contig,
            start,
            end,
            complement: complement.unwrap_or_default(),
        }
    }
    pub(crate) fn intoloc(self) -> io::Result<LocusInfos> {
        match (self.contig, self.start, self.end) {
            (Some(a), Some(b), Some(c)) if a.to_lowercase() != "auto" => Ok(LocusInfos::new(
                self.locus,
                self.haplotype.unwrap_or_default(),
                a,
                b,
                c,
                self.complement,
            )),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No assembly given with auto locus",
            )),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Hash)]
#[repr(C)]
pub(crate) struct LocusInfos {
    pub(crate) locus: Locus,
    pub(crate) haplotype: Haplotype,
    pub(crate) contig: String,
    pub(crate) start: Position,
    pub(crate) end: Position,
    #[serde(skip)]
    pub(crate) complement: Strand,
    #[serde(skip_deserializing)]
    pub(crate) status: OkStatus,
}
impl Serialize for LocusInfos {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut new = self.clone();
        if new.complement.isrev() {
            (new.start, new.end) = (self.end, self.start)
        };
        let mut se = serializer.serialize_tuple_struct("LocusInfos", 6)?;
        se.serialize_field(&new.locus)?;
        se.serialize_field(&new.haplotype)?;
        se.serialize_field(&new.contig)?;
        se.serialize_field(&new.start)?;
        se.serialize_field(&new.end)?;
        se.serialize_field(&new.status)?;
        se.end()
    }
}
impl LocusInfos {
    pub(crate) fn new(
        locus: Locus,
        haplotype: Haplotype,
        contig: String,
        start: Position,
        end: Position,
        complement: Strand,
    ) -> Self {
        Self {
            locus,
            haplotype,
            contig,
            start,
            end,
            complement,
            status: OkStatus::default(),
        }
    }
}
impl Ord for LocusInfos {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.locus.to_string().cmp(&other.locus.to_string()) {
            std::cmp::Ordering::Equal => (),
            ord => return ord,
        };
        self.haplotype.cmp(&other.haplotype)
    }
}
impl PartialOrd for LocusInfos {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl LocusInfos {
    pub(crate) fn getlength(&self) -> i64 {
        self.end.length(&self.start)
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub(crate) struct HashMapinfo {
    pub(crate) locuspos: Position,
    pub(crate) position: Position,
    pub(crate) map60: i64,
    pub(crate) map1: i64,
    pub(crate) map0: i64,
    #[serde(serialize_with = "globalmismatch")]
    pub(crate) globalmismatch: usize,
    pub(crate) overlaps: i64,
    pub(crate) secondary: i64,
    pub(crate) supplementary: i64,
    pub(crate) mismatches: i64,
    pub(crate) misalign: i64,
    #[serde(skip_serializing_if = "iszero")]
    pub(crate) qual: usize,
    #[serde(rename = "percent-softclips")]
    pub(crate) softclips: f32,
}
impl PartialEq for HashMapinfo {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position
    }
}
impl Eq for HashMapinfo {}
impl HashMapinfo {
    pub(crate) fn getmaxvalue(&self) -> i64 {
        let elem = [
            self.map0,
            self.map1,
            self.overlaps,
            self.secondary,
            self.supplementary,
        ];
        elem.into_iter().max().unwrap_or(self.overlaps)
    }
    /// Get the number of primary reads
    pub(crate) fn gettotalmap(&self) -> i64 {
        self.map0 + self.map1 + self.map60
    }
    /// Get the number of primary reads **(except map0)**
    #[allow(dead_code)]
    pub(crate) fn gettotalscore(&self) -> i64 {
        self.map1 + self.map60
    }
}
impl PartialOrd for HashMapinfo {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HashMapinfo {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.getmaxvalue().cmp(&other.getmaxvalue())
    }
}
impl HashMapinfo {
    #[allow(dead_code, clippy::too_many_arguments)]
    pub(crate) fn new(
        locuspos: Position,
        position: Position,
        map60: i64,
        map1: i64,
        map0: i64,
        secondary: i64,
        supplementary: i64,
        globalmismatch: usize,
        overlaps: i64,
        mismatches: i64,
        misalign: i64,
        qual: usize,
        softclips: f32,
    ) -> Self {
        HashMapinfo {
            locuspos,
            position,
            map60,
            map1,
            map0,
            globalmismatch,
            secondary,
            supplementary,
            overlaps,
            mismatches,
            misalign,
            qual,
            softclips,
        }
    }
}
pub(crate) fn iszero(num: &usize) -> bool {
    *num == 0
}
pub(crate) fn globalmismatch<S>(num: &usize, s: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let val = *num as f32 / crate::GLOBALMISMATCHFLOATING as f32;
    s.collect_str(&format!("{val}"))
}

/*
This software allows the analysis of BAM files to identify the confidence on a locus (specifically IG and TR) as well as allele confidence.
It was created and used by IMGT Team (https://www.imgt.org).
Available under EUPL license
*/
use clap::{Parser, crate_authors};
use serde::{Deserialize, Serialize, de};
use std::{fmt::Display, hash::Hash, path::PathBuf};
#[derive(Parser, Debug)]
#[clap(
    author = crate_authors!("\n"),
    before_help = "This script analyzes BAM files coming from reads assembled on an assembly.",
    after_help = "This code was made by and for IMGT (the international ImMunoGeneTics information system).",
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
    ///CSV containing locus infos. See example file for blueprint.
    #[arg(short, long)]
    pub(crate) locuspos: PathBuf,
    /// Minimal number of reads (included) to declare a break in coverage
    #[arg(short, long, default_value_t = 3)]
    pub(crate) breaks: u32,
    /// Coverage to calculate on CSV
    #[arg(short, long, default_value_t = 10)]
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
    /// Huge region
    #[arg(long)]
    pub(crate) hugeregion: bool,
    /// Number of threads to decrypt bgzf files (0 for number of threads up to 12)
    #[arg(long, default_value_t = 0)]
    pub(crate) threads: usize,
    /// Only strand-specific alignments to reference
    #[arg(long)]
    pub(crate) forward: bool,
    /// Query full quality PHRED score (script will be longer to execute)
    #[arg(long)]
    pub(crate) fullquality: bool,
    /// Calculate total reads mismatch
    #[arg(long)]
    pub(crate) totalread: bool,
    /// Size of legend axis (default 16)
    #[arg(long, default_value_t = 16)]
    pub(crate) fontlegendsize: u32,
    /// No legend on graphs
    #[arg(long)]
    pub(crate) nolegend: bool,
    /// Get supplementary and secondary alignments on gene graphs
    #[arg(long)]
    pub(crate) allreads: bool,
    /// Save as SVG images (create big images)
    #[arg(long)]
    pub(crate) svg: bool,
    ///Species name (for folder creation)
    #[arg(short, long)]
    pub(crate) species: String,
    ///Gene location (csv file). See example file for blueprint.
    #[arg(short, long)]
    pub(crate) geneloc: Option<PathBuf>,
    ///Output directory (created or overwritten)
    #[arg(short, long)]
    pub(crate) outdir: PathBuf,
}
pub(crate) fn less_than_100(s: &str) -> Result<u8, String> {
    match s.parse::<u8>() {
        Ok(s) if (0..=100).contains(&s) => Ok(s),
        _ => Err(String::from("Bad number, must be between 0 and 100.")),
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Hash)]
#[allow(clippy::upper_case_acronyms)]
pub(crate) enum Locus {
    IGH,
    IGK,
    IGL,
    TRA,
    TRB,
    TRG,
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
            record.r#match * 100 / record.total
        } else {
            0
        };
        match (
            record.percentalerting,
            record.percentwarning,
            record.minreads,
        ) {
            (d, ..) if percent <= d.into() => Alertpos::Suspicious,
            (_, d, e) if record.r#match <= e.try_into().unwrap() || percent <= d.into() => {
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
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
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
                &"expected i32",
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
        max.checked_sub(min).unwrap().checked_add(1).unwrap() //Calculate the length
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
#[derive(Clone, Debug, Eq, PartialEq, Copy)]
pub(crate) struct Posread {
    pub(crate) r#match: usize,
    pub(crate) indel: usize,
    pub(crate) total: usize,
    pub(crate) minreads: u32,
    pub(crate) percentwarning: u8,
    pub(crate) percentalerting: u8,
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
        self.getindel().checked_sub(self.getmatch()).unwrap()
    }
    pub(crate) fn getindelcount(&self) -> usize {
        self.gettotal().checked_sub(self.getindel()).unwrap()
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
            Locus::TRA => write!(f, "TRA"),
            Locus::TRB => write!(f, "TRB"),
            Locus::TRG => write!(f, "TRG"),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Haplotype {
    Primary,
    Alternate,
}
impl<'de> Deserialize<'de> for Haplotype {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        let s: &str = de::Deserialize::deserialize(deserializer)?;

        match s.to_lowercase().as_str() {
            "primary" | "pri" | "p" => Ok(Haplotype::Primary),
            "alternate" | "alt" | "a" => Ok(Haplotype::Alternate),
            _ => Err(de::Error::unknown_variant(
                s,
                &["primary or pri or p", "alternate or alt or a"],
            )),
        }
    }
}
impl Ord for Haplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (Self::Primary, Self::Primary) | (Self::Alternate, Self::Alternate) => {
                std::cmp::Ordering::Equal
            }
            (Self::Primary, Self::Alternate) => std::cmp::Ordering::Less,
            (Self::Alternate, Self::Primary) => std::cmp::Ordering::Greater,
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
impl PartialOrd for Haplotype {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Haplotype {
    pub(crate) fn isprimary(&self) -> bool {
        self == &Haplotype::Primary
    }
}
#[derive(Clone, Debug, Deserialize)]
pub(crate) struct GeneInfos {
    pub(crate) gene: String,
    pub(crate) chromosome: String,
    pub(crate) strand: Strand,
    pub(crate) start: Position,
    pub(crate) end: Position,
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
#[derive(Clone, Debug, Serialize, PartialEq, Eq, Hash)]
pub(crate) enum Strand {
    Plus,
    Minus,
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
            Self::Minus => write!(f, "RVD"),
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
    length: i64,
    pub(crate) readscoverage: f32,
    pub(crate) reads: usize,
    pub(crate) matchpos: String,
    pub(crate) readsfull: usize,
    pub(crate) reads100: usize,
    pub(crate) reads100m: usize,
    pub(crate) coveragex: usize,
}
impl Ord for GeneInfosFinish {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.chromosome.cmp(&other.chromosome) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        match self.start.cmp(&other.start) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        match self.end.cmp(&other.end) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        std::cmp::Ordering::Equal
    }
}
impl PartialOrd for GeneInfosFinish {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for GeneInfosFinish {
    fn eq(&self, other: &Self) -> bool {
        self.gene == other.gene
    }
}
impl Eq for GeneInfosFinish {}
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
        }
    }
    pub(crate) fn make_default(gene: GeneInfos) -> Self {
        Self::new(gene, 0, 0, None, 0, 0, 0.0, 0)
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Hash)]
pub(crate) struct LocusInfos {
    pub(crate) locus: Locus,
    pub(crate) haplotype: Haplotype,
    pub(crate) contig: String,
    pub(crate) start: Position,
    pub(crate) end: Position,
    #[serde(skip)]
    pub(crate) complement: bool,
}
#[derive(Debug, Clone, Serialize, Default)]
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
        elem.into_iter().max().unwrap()
    }
    pub(crate) fn gettotalmap(&self) -> i64 {
        self.map0 + self.map1 + self.map60
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
    s.collect_str(&format!("{}", val))
}

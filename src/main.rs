use clap::{crate_authors, Parser};
use colors::full_palette::GREY_300;
use plotters::coord::Shift;
use std::collections::BTreeSet;
use std::io::{stderr, stdout};
use std::num::NonZero;
use std::time::Instant;
//use noodles_fasta::{self as fasta, record::Sequence};
use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::prelude::{
    full_palette, AreaSeries, BitMapBackend, DrawingArea, DrawingBackend, IntoDrawingArea,
    IntoSegmentedCoord, PathElement, SVGBackend,
};
use plotters::series::{Histogram, LineSeries};
use plotters::style::*;
use rust_htslib::bam::ext::{BamRecordExtensions, IterAlignedPairs};
use rust_htslib::bam::{self, IndexedReader, Read};
use serde::{de, Deserialize, Serialize};
use std::{
    cmp::{max, min},
    collections::BTreeMap,
    fmt::Display,
    fs::File,
    io::Write,
    path::PathBuf,
};
///Assess quality of an assembly based on reads mapping
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
struct Args {
    /// Input file (SAM or BAM)
    #[arg(short, long)]
    file: PathBuf,
    /// Index file if not default
    #[arg(short, long)]
    index: Option<PathBuf>,
    ///CSV containing locus infos
    #[arg(short, long)]
    locuspos: PathBuf,
    /// Minimal number of reads (included) to declare a break in coverage
    #[arg(short, long, default_value_t = 3)]
    breaks: u32,
    /// Coverage to calculate on CSV
    #[arg(short, long, default_value_t = 10)]
    coverage: u32,
    /// Force cigar even if no =
    #[arg(long)]
    force: bool,
    /// Save as SVG images
    #[arg(long)]
    svg: bool,
    ///Species
    #[arg(short, long)]
    species: String,
    ///Gene location (csv file)
    #[arg(short, long)]
    geneloc: Option<PathBuf>,
    ///Output directory
    #[arg(short, long)]
    outdir: PathBuf,
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
#[allow(clippy::upper_case_acronyms)]
enum Locus {
    IGH,
    IGK,
    IGL,
    TRA,
    TRB,
    TRG,
}
#[derive(Clone, Debug, Eq, PartialEq, Copy, Default)]
struct Posread {
    r#match: usize,
    indel: usize,
    total: usize,
}
impl Posread {
    #[allow(dead_code)]
    fn new(r#match: usize, indel: usize, total: usize) -> Result<Self, &'static str> {
        if r#match + indel > total {
            return Err("Invalid total");
        }
        Ok(Self {
            r#match,
            indel,
            total,
        })
    }
    fn gettotal(&self) -> usize {
        self.total
    }
    fn addtotal(&mut self, count: usize) {
        self.total += count
    }
    fn getmatch(&self) -> usize {
        self.r#match
    }
    fn addmatch(&mut self, count: usize) {
        self.r#match += count
    }
    fn getindel(&self) -> usize {
        self.indel
    }
    fn addindel(&mut self, count: usize) {
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
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
enum Haplotype {
    Primary,
    Alternate,
}
impl Ord for Haplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if (self == &Haplotype::Primary && other == &Haplotype::Primary)
            || (self == &Haplotype::Alternate && other == &Haplotype::Alternate)
        {
            std::cmp::Ordering::Equal
        } else if self == &Haplotype::Primary {
            std::cmp::Ordering::Less
        } else if other == &Haplotype::Primary {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
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
    fn isprimary(&self) -> bool {
        self == &Haplotype::Primary
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
struct GeneInfos {
    gene: String,
    chromosome: String,
    strand: Strand,
    start: i64,
    end: i64,
}
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
enum Strand {
    Plus,
    Minus,
}
impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        let s: &str = de::Deserialize::deserialize(deserializer)?;

        match s {
            "1" => Ok(Strand::Minus),
            "0" => Ok(Strand::Plus),
            _ => Err(de::Error::unknown_variant(s, &["1", "0"])),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize)]
struct GeneInfosFinish {
    gene: String,
    chromosome: String,
    strand: Strand,
    start: i64,
    end: i64,
    length: i64,
    coverageperc: f32,
    reads: usize,
    matchpos: String,
    reads100: usize,
    reads100m: usize,
    coveragex: usize,
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
struct LocusInfos {
    locus: Locus,
    haplotype: Haplotype,
    contig: String,
    start: i64,
    end: i64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct HashMapinfo {
    map60: i64,
    map1: i64,
    map0: i64,
    overlaps: i64,
    secondary: i64,
    supplementary: i64,
    mismatches: i64,
    misalign: i64,
    qual: usize,
}
impl HashMapinfo {
    fn getmaxvalue(&self) -> i64 {
        let elem = [
            self.map0,
            self.map1,
            self.overlaps,
            self.secondary,
            self.supplementary,
        ];
        *elem.iter().max().unwrap()
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
    fn new(
        map60: i64,
        map1: i64,
        map0: i64,
        secondary: i64,
        supplementary: i64,
        overlaps: i64,
        mismatches: i64,
        misalign: i64,
        qual: usize,
    ) -> Self {
        HashMapinfo {
            map60,
            map1,
            map0,
            secondary,
            supplementary,
            overlaps,
            mismatches,
            misalign,
            qual,
        }
    }
    fn default() -> Self {
        HashMapinfo {
            map60: 0,
            map1: 0,
            map0: 0,
            secondary: 0,
            supplementary: 0,
            overlaps: 0,
            mismatches: 0,
            misalign: 0,
            qual: 0,
        }
    }
}
fn getreaderoffile(args: &Args) -> IndexedReader {
    let mut reader = match &args.index {
        Some(d) => bam::IndexedReader::from_path_and_index(&args.file, d),
        None => bam::IndexedReader::from_path(&args.file),
    }
    .unwrap();
    let threads = match std::thread::available_parallelism() {
        Ok(d) => max(d, NonZero::new(12).unwrap()),
        Err(_) => NonZero::new(4).unwrap(),
    };
    reader.set_threads(threads.get()).unwrap();
    reader
}
fn main() {
    let args = Args::parse();
    let mut csv = csv::ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(&args.locuspos)
        .unwrap();
    let mut locus: Vec<LocusInfos> = Vec::new();
    for record in csv.deserialize() {
        let record = record.expect("Invalid CSV format, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend");
        locus.push(record);
    }
    if locus.is_empty() {
        panic!("Invalid CSV format, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend");
    }
    let outputdir = match args.outdir.is_dir() {
        true => &args.outdir,
        false => {
            eprintln!(
                "Folder {} does not exist, attempt to create.",
                args.outdir.display()
            );
            std::fs::create_dir(&args.outdir).unwrap();
            &args.outdir
        }
    };
    locus.sort_unstable_by(|a, b| match a.locus.to_string().cmp(&b.locus.to_string()) {
        std::cmp::Ordering::Equal => a.haplotype.cmp(&b.haplotype),
        o => o,
    });
    let grouped = locus.into_iter().into_group_map_by(|f| f.locus.to_string());
    for (_, locus) in grouped {
        let floci = locus.first().unwrap();
        let haplotype = locus.len();
        if haplotype > 2 {
            panic!("There is more than 2 haplotypes for {}", floci.locus);
        }
        let haplotypebool = haplotype == 1;
        println!(
            "Going for {} locus - {}",
            floci.locus,
            if !haplotypebool { "diploid" } else { "haploid" }
        );
        let mut outputfile1 = PathBuf::new();
        let mut outputfile2 = PathBuf::new();
        let mut outputfile3 = PathBuf::new();
        let mut outputfile4 = PathBuf::new();
        let (
            readgraphtop,
            readgraphbottom,
            mismatchgraphtop,
            mismatchgraphbottom,
            readgraphtop2,
            readgraphbottom2,
            mismatchgraphtop2,
            mismatchgraphbottom2,
        ) = {
            // Second graph with reads
            let fgraph = "readresult";
            let (top, bottom, topb, bottomb) = if args.svg {
                let outputfile = outputdir.join(givename(
                    &args.species,
                    &floci.locus,
                    &floci.contig,
                    haplotypebool,
                    &format!("{}.svg", fgraph),
                    true,
                ));
                outputfile1 = outputfile;
                let root = SVGBackend::new(
                    &outputfile1,
                    (
                        900,
                        500 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
                    ),
                )
                .into_drawing_area();
                if haplotypebool {
                    let (top, bottom) = root.split_vertically((50).percent_height());
                    (Some(top), Some(bottom), None, None)
                    //readgraph(outputfile, locus.first().unwrap(), &pos, &args, root);
                } else {
                    (Some(root), None, None, None)
                }
            } else {
                let outputfile = outputdir.join(givename(
                    &args.species,
                    &floci.locus,
                    &floci.contig,
                    haplotypebool,
                    &format!("{}.png", fgraph),
                    true,
                ));
                outputfile2 = outputfile;
                let root = BitMapBackend::new(
                    &outputfile2,
                    (
                        1800,
                        1000 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
                    ),
                )
                .into_drawing_area();
                if !haplotypebool {
                    let (top, bottom) = root.split_vertically((50).percent_height());
                    (None, None, Some(top), Some(bottom))
                    //readgraph(outputfile, locus.first().unwrap(), &pos, &args, root);
                } else {
                    (None, None, Some(root), None)
                }
            };
            // Second graph with mismatches
            let sgraph = "mismatchresult";
            let (mistop, misbottom, mistopb, misbottomb) = if args.svg {
                let outputfile = outputdir.join(givename(
                    &args.species,
                    &floci.locus,
                    &floci.contig,
                    haplotype == 1,
                    &format!("{}.svg", sgraph),
                    true,
                ));
                outputfile3 = outputfile;
                let root = SVGBackend::new(
                    &outputfile3,
                    (
                        1200,
                        600 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
                    ),
                )
                .into_drawing_area();
                if !haplotypebool {
                    let (top, bottom) = root.split_vertically((50).percent_height());
                    (Some(top), Some(bottom), None, None)
                    //readgraph(outputfile, locus.first().unwrap(), &pos, &args, root);
                } else {
                    (Some(root), None, None, None)
                }
                //mismatchgraph(outputfile, floci, &pos, &args, root);
            } else {
                let outputfile = outputdir.join(givename(
                    &args.species,
                    &floci.locus,
                    &floci.contig,
                    haplotype == 1,
                    &format!("{}.png", sgraph),
                    true,
                ));
                outputfile4 = outputfile;
                let root = BitMapBackend::new(
                    &outputfile4,
                    (
                        1200,
                        600 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
                    ),
                )
                .into_drawing_area();
                if !haplotypebool {
                    let (top, bottom) = root.split_vertically((50).percent_height());
                    (None, None, Some(top), Some(bottom))
                    //readgraph(outputfile, locus.first().unwrap(), &pos, &args, root);
                } else {
                    (None, None, Some(root), None)
                }
                //mismatchgraph(outputfile, floci, &pos, &args, root);
            };
            (
                top, bottom, mistop, misbottom, topb, bottomb, mistopb, misbottomb,
            )
        };
        let mut lock = stdout().lock();
        for loci in locus.iter() {
            let mut reader = getreaderoffile(&args);
            reader
                .fetch((&loci.contig, loci.start, loci.end + 1))
                .unwrap_or_else(|_| {
                    panic!(
                        "The region {}:{}-{} cannot be found, exiting.",
                        loci.contig, loci.start, loci.end
                    )
                });
            let mut nocount = true;
            //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
            //let file = File::create(&filename).unwrap();
            //let mut writer = BufWriter::new(file);
            let mut pos: BTreeMap<i64, HashMapinfo> = BTreeMap::new();
            (loci.start..=loci.end).for_each(|p| {
                pos.insert(p, HashMapinfo::default());
            });
            let mut message = false;
            println!("Region {} fetched, analyzing all reads.", loci.locus);
            let mut count = 0;
            let time = Instant::now();
            //let sep = max((loci.end - loci.start + 1) / 250, 100); //250 points for quality point
            for p in reader.rc_records().filter_map(Result::ok) {
                count += 1;
                if count % 100 == 0 {
                    writeln!(
                        lock,
                        "Process {} reads in {} s",
                        count.to_formatted_string(&Locale::en),
                        Instant::now().saturating_duration_since(time).as_secs_f32()
                    )
                    .unwrap();
                }
                nocount = false;
                let newrange =
                    max(p.reference_start(), loci.start)..min(loci.end + 1, p.reference_end());
                if p.is_secondary() || p.is_supplementary() {
                    newrange.for_each(|i| {
                        let targeting = pos.get_mut(&i).unwrap();
                        if p.is_secondary() {
                            targeting.secondary += 1;
                        } else {
                            targeting.supplementary += 1;
                            /* targeting
                            .supplementary
                            .push(String::from_utf8_lossy(p.qname()).to_string()); */
                        }
                    });
                    continue;
                }
                let overlaprange = max(p.reference_start() + 1, loci.start)
                    ..min(loci.end + 1, p.reference_end() - 1); //Remove 1 because does not overlap on borders
                let matched: BTreeSet<(i64, i64)> = match p.aligned_pairs_match() {
                    Some(a) => a.map(|[a, b]| (a, b)).collect(),
                    None => {
                        let text = "No = CIGAR given";
                        if !args.force {
                            panic!("{}", text);
                        } else if !message {
                            eprintln!("{} but it was forced.", text);
                            message = true
                        }
                        std::iter::empty::<IterAlignedPairs>()
                            .map(|_| (0, 0))
                            .collect()
                    }
                };
                let aligned: BTreeSet<(i64, i64)> =
                    p.aligned_pairs().map(|[a, b]| (a, b)).collect();
                for (i, targeting) in pos.range_mut(newrange) {
                    match p.mapq() {
                        0 => targeting.map0 += 1,
                        1..=59 => targeting.map1 += 1,
                        60 => targeting.map60 += 1,
                        _ => panic!(
                            "MAPQ score is invalid. Got {} for {}",
                            p.mapq(),
                            String::from_utf8_lossy(p.qname())
                        ),
                    };
                    //Get quality of reads depending on genomic position
                    if let Some(d) = p.aligned_pairs().find(|p| p[1] == *i) {
                        let index = d[0] as usize;
                        targeting.qual += *p.qual().get(index).unwrap() as usize;
                    }
                    targeting.misalign += 1;
                }
                //Overlap
                for (_, targeting) in pos.range_mut(overlaprange) {
                    targeting.overlaps += 1;
                }
                //Remove misalign and mismatched
                matched.iter().for_each(|p| {
                    pos.range_mut(p.0..=p.1).for_each(|(_, f)| {
                        f.misalign -= 1;
                    })
                });
                //Aligned reads but without matches ones
                aligned.difference(&matched).for_each(|p| {
                    pos.range_mut(p.0..=p.1).for_each(|(_, f)| {
                        f.misalign -= 1;
                        //Aligns but not correct nt
                        if !args.force {
                            f.mismatches += 1;
                        };
                    })
                });
            }
            if nocount {
                panic!(
                    "The region {}:{}-{} cannot be found, exiting.",
                    loci.contig, loci.start, loci.end
                )
            }
            println!("Making graphs");
            match (loci.haplotype.isprimary(), args.svg) {
                (true, true) => {
                    readgraph(
                        &outputfile1,
                        loci,
                        &pos,
                        &args,
                        readgraphtop.clone().unwrap(),
                    );
                    mismatchgraph(
                        &outputfile3,
                        loci,
                        &pos,
                        &args,
                        mismatchgraphtop.clone().unwrap(),
                    );
                }
                (false, true) => {
                    readgraph(
                        &outputfile1,
                        loci,
                        &pos,
                        &args,
                        readgraphbottom.clone().unwrap(),
                    );
                    mismatchgraph(
                        &outputfile3,
                        loci,
                        &pos,
                        &args,
                        mismatchgraphbottom.clone().unwrap(),
                    );
                }
                (true, false) => {
                    readgraph(
                        &outputfile2,
                        loci,
                        &pos,
                        &args,
                        readgraphtop2.clone().unwrap(),
                    );
                    mismatchgraph(
                        &outputfile4,
                        loci,
                        &pos,
                        &args,
                        mismatchgraphtop2.clone().unwrap(),
                    );
                }
                (false, false) => {
                    readgraph(
                        &outputfile2,
                        loci,
                        &pos,
                        &args,
                        readgraphbottom2.clone().unwrap(),
                    );
                    mismatchgraph(
                        &outputfile4,
                        loci,
                        &pos,
                        &args,
                        mismatchgraphbottom2.clone().unwrap(),
                    );
                }
            }
            //Create CSV from HashMap
            createcsv(outputdir, loci, &pos, &args);
        }
        //Create gene CSV
        if args.geneloc.is_some() {
            println!("Gene list starting!");
            genelist(outputdir, floci, &args);
            println!("Gene list finished");
        }
        println!("Locus {} is done!", &floci.locus);
    }
}
fn genelist(outputdir: &std::path::Path, floci: &LocusInfos, args: &Args) {
    let outputfile = outputdir.join(givename(
        &args.species,
        &floci.locus,
        &floci.contig,
        floci.haplotype.isprimary(),
        "geneanalysis.csv",
        false,
    ));
    let mut lock: std::io::StderrLock<'_> = stderr().lock();
    let mut csv = csv::ReaderBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b',')
        .from_path(args.geneloc.as_ref().unwrap())
        .unwrap();
    let mut genes: Vec<GeneInfos> = Vec::new();
    for record in csv.deserialize() {
        let record = record
            .expect("Invalid CSV format, waiting gene,chromosome,strand,start,end case sensitive");
        genes.push(record);
    }
    if genes.is_empty() {
        panic!("Invalid CSV format, waiting gene,chromosome,strand,start,end case sensitive");
    }
    let mut finale: Vec<GeneInfosFinish> = Vec::with_capacity(genes.len());
    for mut gene in genes {
        /* if gene.chromosome != loci.contig
                || !(loci.start..=loci.end).contains(&gene.start)
                || !(loci.start..=loci.end).contains(&gene.end)
        {
            continue; //gene not on this loci
        } */
        let mut reader = getreaderoffile(args);
        let (mut reads, mut reads100, mut reads100m) = (0, 0, 0);
        if gene.start > gene.end {
            (gene.end, gene.start) = (gene.start, gene.end) //Swap position
        }
        let range = ranges::Ranges::from(vec![gene.start..=gene.end]);
        reader
            .fetch((&gene.chromosome, gene.start + 1, gene.end + 1))
            .unwrap();
        let records = reader.records();
        let mut hash: BTreeMap<i64, Posread> = BTreeMap::new(); //Match and full match and total
        range.clone().into_iter().for_each(|p| {
            hash.insert(p, Posread::default());
        });
        let mut coverageperc = 0;
        let mut empty = true;
        for record in records
            .filter_map(Result::ok)
            .filter(|p| !p.is_secondary() && !p.is_supplementary())
        {
            empty = false;
            reads += 1;
            coverageperc +=
                ranges::Ranges::from(record.reference_start() + 1..record.reference_end() + 1)
                    .into_iter()
                    .count();
            for p in record.reference_start() + 1..record.reference_end() + 1 {
                match hash.get_mut(&p) {
                    Some(d) => d.addtotal(1),
                    None => {
                        if record.reference_start() > gene.end {
                            break;
                        }
                    } //Outside coverage of gene
                }
            }
            'outer: for [start, end] in record.aligned_blocks() {
                for p in start..end {
                    match hash.get_mut(&p) {
                        Some(d) => d.addindel(1),
                        None => {
                            if start > gene.end {
                                break 'outer;
                            }
                        } //Outside coverage of gene
                    }
                }
            }
            if !args.force {
                'outer: for [start, end] in record.aligned_blocks_match().unwrap() {
                    for p in start..end {
                        match hash.get_mut(&p) {
                            Some(d) => d.addmatch(1),
                            None => {
                                if start > gene.end {
                                    break 'outer;
                                }
                            } //Outside coverage of gene
                        }
                    }
                }
            }
            if record
                .aligned_blocks()
                .any(|p| p[0] < gene.start && p[1] + 1 > gene.end)
            {
                reads100 += 1;
            }
            if !args.force
                && record
                    .aligned_blocks_match()
                    .unwrap()
                    .any(|p| p[0] < gene.start && p[1] + 1 > gene.end)
            {
                reads100m += 1;
            }
        }
        if empty {
            writeln!(lock, "Empty records for gene {}", gene.gene).unwrap();
            continue;
        }
        let coverage = hash
            .clone()
            .into_values()
            .filter(|p| p.gettotal() >= args.coverage.try_into().unwrap())
            .count();
        let text = hash.iter().fold(String::new(), |mut acc, (_, f)| {
            acc.push_str(&format!(
                "{}/{}({})-",
                f.getindel(),
                f.gettotal(),
                f.getmatch()
            ));
            acc
        });
        let text = String::from(text.trim_end_matches('-'));
        let genename = gene.gene.clone();
        let plots = outputdir
            .join(format!("gene_{}", args.species))
            .join(floci.haplotype.to_string().as_str().to_lowercase());
        if !std::fs::exists(&plots).unwrap() {
            println!("Creating the folder {}", plots.display());
            std::fs::create_dir_all(&plots).unwrap();
        };
        let mut output = plots.join(&genename);
        output.set_extension("png");
        let root = BitMapBackend::new(&output, (700, 400)).into_drawing_area();
        //Gene graph
        genegraph(&hash, &gene, floci, root);
        let elem = GeneInfosFinish {
            gene: gene.gene,
            chromosome: gene.chromosome,
            strand: gene.strand,
            start: gene.start,
            end: gene.end,
            length: gene
                .end
                .checked_sub(gene.start)
                .unwrap()
                .checked_add(1)
                .unwrap(),
            reads,
            matchpos: text,
            reads100,
            reads100m,
            coverageperc: ((coverageperc / reads / range.into_iter().count()) as f32 * 1000.0)
                .round()
                / 1000.0,
            coveragex: coverage,
        };
        finale.push(elem);
    }
    let mut csv = csv::WriterBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(&outputfile)
        .unwrap();
    for gene in finale {
        csv.serialize(gene).unwrap();
    }
    csv.flush().unwrap();
    println!("Gene analysis has been saved to {}", outputfile.display());
}
fn genegraph<T>(
    hash: &BTreeMap<i64, Posread>,
    gene: &GeneInfos,
    loci: &LocusInfos,
    root: DrawingArea<T, Shift>,
) where
    T: DrawingBackend,
{
    let genename = gene.gene.to_string();
    let text_style = ("sans-serif", 14, &BLACK).into_text_style(&root);
    let _ = root.fill(&plotters::prelude::WHITE);
    let max = hash.values().map(|p| p.gettotal()).max().unwrap() + 5;
    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .right_y_label_area_size(40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            format!("Reads alignment for {} ({})", genename, loci.haplotype),
            ("sans-serif", 22),
        )
        .build_cartesian_2d(0..hash.len(), 0..max)
        .unwrap();
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos, val.getmatch())),
            full_palette::RED_200.mix(0.8),
        ))
        .unwrap()
        .label("Equal")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::RED_200));
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos, val.getindel())),
            full_palette::GREEN_600.mix(0.8),
        ))
        .unwrap()
        .label("Match")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::GREEN_600));
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos, val.gettotal())),
            full_palette::LIGHTBLUE_300.mix(0.8),
        ))
        .unwrap()
        .label("Total")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::LIGHTBLUE_300));
    chart
        .configure_mesh()
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        .x_desc("Position in sequence (bp)")
        .y_desc("Reads count")
        .x_label_style(text_style.clone())
        .y_label_style(text_style)
        .light_line_style(GREY_300)
        .x_max_light_lines(5)
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .unwrap();
    chart
        .configure_series_labels()
        .position(plotters::chart::SeriesLabelPosition::LowerRight)
        .background_style(WHITE.mix(0.6))
        .border_style(BLACK.mix(0.8))
        .draw()
        .unwrap();
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
}
fn givename(
    species: &str,
    locus: &Locus,
    contig: &str,
    haplo: bool,
    suffix: &str,
    image: bool,
) -> String {
    format!(
        "{}_{}_{}{}_{}",
        species,
        locus,
        if haplo {
            format!("{}-", contig)
        } else {
            String::new()
        },
        if haplo {
            "primary"
        } else if image {
            "full"
        } else {
            "alternate"
        },
        suffix
    )
}
fn createcsv(
    outputdir: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
) {
    let outputfile = outputdir.join(givename(
        &args.species,
        &loci.locus,
        &loci.contig,
        loci.haplotype.isprimary(),
        "positionresult.csv",
        false,
    ));
    let outputfile = outputfile.as_path();
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .flexible(false)
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(true)
        .from_path(outputfile)
        .unwrap();
    csv.write_record([
        "Position",
        "map60",
        "map1",
        "map0",
        "overlaps",
        "secondary",
        "mismatches",
        "misalign",
        "qual",
    ])
    .unwrap();
    for (pos, record) in pos.iter() {
        csv.write_field(format!("{}", pos)).unwrap();
        csv.serialize(record).unwrap();
    }
    csv.flush().unwrap();
    println!("CSV analysis has been saved to {}", outputfile.display());
}
fn mismatchgraph<T>(
    _outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    _args: &Args,
    root: DrawingArea<T, Shift>,
) where
    T: DrawingBackend,
{
    let _ = root.fill(&plotters::prelude::WHITE);
    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Mismatches rate and quality over the locus {} ({}-{})",
                loci.locus, loci.contig, loci.haplotype
            ),
            ("sans-serif", 28),
        )
        .build_cartesian_2d(loci.start..loci.end, 0..100)
        .unwrap();
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .style(full_palette::BROWN.mix(0.8).filled())
                .margin(0)
                .data(pos.iter().map(|p| {
                    (
                        *p.0,
                        (p.1.mismatches as f64 / (p.1.map0 + p.1.map1 + p.1.map60) as f64 * 100.0)
                            .round() as i32,
                    )
                })),
        )
        .unwrap()
        .label("Mismatches (%)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BROWN));
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .style(full_palette::RED_400.mix(0.8).filled())
                .margin(0)
                .data(pos.iter().map(|p| {
                    (
                        *p.0,
                        (p.1.misalign as f64 / (p.1.map0 + p.1.map1 + p.1.map60) as f64 * 100.0)
                            .round() as i32,
                    )
                })),
        )
        .unwrap()
        .label("Misalign (%)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::RED_400));
    let mut secondary = chart.set_secondary_coord(loci.start..loci.end, 0..100usize);
    secondary
        .configure_mesh()
        .y_label_formatter(&|f| format!("{}%", f))
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        .x_desc("Genomic position (bp)")
        .y_desc("Mismatch rate (%)")
        .disable_x_mesh()
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .unwrap();
    //let mut second = chart.set_secondary_coord(loci.start..loci.end, 0..max);
    secondary
        .draw_secondary_series(LineSeries::new(
            pos.iter().filter_map(|p| {
                if p.1.qual > 0 {
                    Some((
                        *p.0,
                        p.1.qual
                            / std::convert::TryInto::<usize>::try_into(
                                p.1.map0 + p.1.map1 + p.1.map60,
                            )
                            .unwrap(),
                    ))
                } else {
                    None
                }
            }),
            full_palette::BLACK.mix(0.4),
        ))
        .unwrap()
        .label("Quality")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLACK));
    secondary
        .configure_secondary_axes()
        .y_desc("Quality (PHRED score)")
        //.disable_y_mesh()
        .draw()
        .unwrap();
    secondary
        .configure_series_labels()
        .position(plotters::chart::SeriesLabelPosition::UpperRight)
        .background_style(WHITE)
        .border_style(BLACK.mix(0.8))
        .draw()
        .unwrap();
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
}
fn readgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
    root: DrawingArea<T, Shift>,
) where
    T: DrawingBackend,
{
    let max = pos.values().map(|max| max.getmaxvalue()).max().unwrap() + 5;
    let _ = root.fill(&plotters::prelude::WHITE);
    let (top, bottom) = root.split_vertically((80).percent_height());
    let mut chart = ChartBuilder::on(&top)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Reads mapping quality over the locus {} ({}-{})",
                loci.locus, loci.contig, loci.haplotype
            ),
            ("sans-serif", 28),
        )
        .build_cartesian_2d(loci.start..loci.end, 0..max)
        .unwrap();
    let _ = chart
        .configure_mesh()
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        //.x_desc("Genomic position (bp)")
        .y_desc("Coverage")
        .disable_x_mesh()
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw();
    chart
        .draw_series(AreaSeries::new(
            pos.iter().map(|p| (*p.0, p.1.map0)),
            0,
            full_palette::RED_300.mix(0.6),
        ))
        .unwrap()
        .label("MAPQ: 0")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::RED_300.filled(),
            )
        });
    chart
        .draw_series(AreaSeries::new(
            pos.iter().map(|p| (*p.0, p.1.map1)),
            0,
            full_palette::YELLOW_900.mix(0.5),
        ))
        .unwrap()
        .label("MAPQ: 1-59")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::YELLOW_900.filled(),
            )
        });
    chart
        .draw_series(AreaSeries::new(
            pos.iter().map(|p| (*p.0, p.1.map60)),
            0,
            full_palette::GREEN_400.mix(0.4),
        ))
        .unwrap()
        .label("MAPQ: 60")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::GREEN_400.filled(),
            )
        });
    chart
        .draw_series(LineSeries::new(
            pos.iter().map(|p| (*p.0, p.1.secondary)),
            full_palette::BLACK,
        ))
        .unwrap()
        .label("Secondary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLACK));
    chart
        .draw_series(LineSeries::new(
            pos.iter().map(|p| (*p.0, p.1.supplementary)),
            full_palette::BLUE_700,
        ))
        .unwrap()
        .label("Supplementary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLUE_700));
    chart
        .draw_series(LineSeries::new(
            pos.iter().map(|p| (*p.0, p.1.overlaps)),
            full_palette::ORANGE_300,
        ))
        .unwrap()
        .label("Overlapping reads")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::ORANGE_300));
    chart
        .configure_series_labels()
        .position(plotters::chart::SeriesLabelPosition::UpperRight)
        .background_style(WHITE)
        .border_style(BLACK.mix(0.8))
        .draw()
        .unwrap();
    //Bottom graph
    let mut chart = ChartBuilder::on(&bottom)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        /*.caption(
            format!(
                "Break in coverage {} ({})",
                loci.locus, loci.contig
            ),
            ("sans-serif", 40),
        )  */
        .build_cartesian_2d((loci.start..loci.end).into_segmented(), 0..1i64)
        .unwrap();
    let text_style = ("sans-serif", 14, &BLACK).into_text_style(&root);
    let _ = chart
        .configure_mesh()
        .x_desc("Genomic position (bp)")
        .x_label_style(text_style)
        .disable_x_axis()
        .draw();
    let breaks = pos.iter().filter_map(|(pos, elem)| {
        if elem.overlaps <= args.breaks.into() {
            Some((*pos, 1))
        } else {
            None
        }
    });
    let mut breakfile = File::create(outputfile.parent().unwrap().join(givename(
        &args.species,
        &loci.locus,
        &loci.contig,
        loci.haplotype.isprimary(),
        "break.txt",
        false,
    )))
    .unwrap();
    let mut result = Vec::new();
    let mut prev = None;
    breaks.clone().fold((), |_, (num, _)| {
        if let Some(prev_num) = prev {
            if num - prev_num != 1 {
                result.push(format!("{}:{}..{}", loci.contig, prev_num, num));
            }
        } else {
            result.push(format!("{}:{}", loci.contig, num));
        }
        prev = Some(num);
    });
    let breakcode = result.into_iter().fold(String::new(), |mut acc, f| {
        acc.push_str(&format!("Break: {}\n", f));
        acc
    });
    breakfile.write_all(breakcode.trim().as_bytes()).unwrap();
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .style(full_palette::RED.filled())
                .data(breaks)
                .margin(0),
        )
        .unwrap()
        .label("Coverage break")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
    chart
        .configure_series_labels()
        .position(plotters::chart::SeriesLabelPosition::UpperRight)
        .background_style(WHITE)
        .border_style(BLACK.mix(0.8))
        .draw()
        .unwrap();
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
}

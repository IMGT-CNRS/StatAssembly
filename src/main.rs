use bio_types::genome::AbstractInterval;
use clap::{crate_authors, Parser};
use csv;
use noodles_core::Region;
use noodles_fasta::{self as fasta, record::Sequence};
use plotters::{self, coord};
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::prelude::{full_palette, AreaSeries, BitMapBackend, IntoDrawingArea, PathElement};
use plotters::series::LineSeries;
use plotters::style::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read};
use serde::Deserialize;
use num_format::{ToFormattedString,Locale};
use core::num;
use std::{
    cmp::{max, min},
    collections::BTreeMap,
    fmt::Display,
    fs::File,
    io::{BufReader, BufWriter, Write},
    ops::{Range, RangeInclusive},
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
    ///Species
    #[arg(short, long)]
    species: String,
    ///Assembly file (FASTA)
    #[arg(short, long)]
    assembly: PathBuf,
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
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
struct LocusInfos {
    locus: Locus,
    haplotype: Haplotype,
    contig: String,
    start: i64,
    end: i64,
}
#[derive(Debug, Clone,PartialEq, Eq)]
struct HashMapinfo {
    map60: i64,
    map1: i64,
    map0: i64,
    secondary: i64
}
impl PartialOrd for HashMapinfo {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HashMapinfo {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let a = std::cmp::max(std::cmp::max(std::cmp::max(self.map0,self.map1),self.map60),self.secondary);
        let b = std::cmp::max(std::cmp::max(std::cmp::max(other.map0,other.map1),other.map60),other.secondary);
        a.cmp(&b)
    }
}
impl HashMapinfo {
    fn new(map60: i64, map1: i64, map0: i64, secondary: i64) -> Self {
        HashMapinfo { map60, map1, map0, secondary }
    }
    fn default() -> Self {
        HashMapinfo { map60: 0, map1: 0, map0: 0, secondary: 0 }
    }
}
fn main() {
    let args = Args::parse();
    let path = &args.file;
    let mut csv = csv::ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(args.locuspos)
        .unwrap();
    let mut locus: Vec<LocusInfos> = Vec::new();
    for record in csv.deserialize() {
        let record = record.expect("Invalid CSV format, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend");
        locus.push(record);
    }
    if locus.is_empty() {
        panic!("Invalid CSV format, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend");
    }
    let outputdir = &match args.outdir.is_dir() {
        true => args.outdir,
        false => {
            eprintln!(
                "Folder {} does not exist, attempt to create.",
                args.outdir.display()
            );
            std::fs::create_dir(&args.outdir).unwrap();
            args.outdir
        }
    };
    for loci in locus {
        let mut reader = match &args.index {
            Some(d) => bam::IndexedReader::from_path_and_index(path, d),
            None => bam::IndexedReader::from_path(path),
        }
        .unwrap();
        reader
            .fetch((&loci.contig, loci.start, loci.end + 1))
            .unwrap();
        //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
        //let file = File::create(&filename).unwrap();
        //let mut writer = BufWriter::new(file);
        let mut pos: BTreeMap<i64, HashMapinfo> = BTreeMap::new();
        (loci.start..=loci.end).for_each(|p| {
            pos.insert(p, HashMapinfo::default());
        });
        for p in reader.rc_records() {
            let p = p.unwrap();
            let newrange =
                max(p.reference_start(), loci.start)..min(loci.end + 1, p.reference_end());
            if p.is_secondary() {
                newrange.for_each(|i| {
                    pos.get_mut(&i).unwrap().secondary += 1;
                });
                continue;
            }
            newrange.for_each(|i| match p.mapq() {
                0 => pos.get_mut(&i).unwrap().map0 += 1,
                1..=59 => pos.get_mut(&i).unwrap().map1 += 1,
                60 => pos.get_mut(&i).unwrap().map60 += 1,
                _ => panic!(
                    "MAPQ score is invalid. Got {} for {}",
                    p.mapq(),
                    String::from_utf8_lossy(p.qname())
                ),
            });
        }
        let max = pos.values().map(|max| {
            std::cmp::max(std::cmp::max(std::cmp::max(max.map0,max.map1),max.map60),max.secondary)
        }).max().unwrap()+5;
        //writer.flush().unwrap();
        let outputfile = "result.png";
        let root = BitMapBackend::new(outputfile, (1200, 600)).into_drawing_area();
        let _ = root.fill(&plotters::prelude::WHITE);
        let mut chart = ChartBuilder::on(&root)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            .caption(format!("Reads mapping quality over the locus {}",loci.locus), ("sans-serif", 40))
            .build_cartesian_2d(loci.start..loci.end, 0..max).unwrap();
        let _ = chart
            .configure_mesh()
            .x_label_formatter(&|f| { f.to_formatted_string(&Locale::en).to_string() })
            .disable_x_mesh()
            .y_max_light_lines(3)
            //.disable_y_mesh()
            .draw();
        chart.draw_series(
            AreaSeries::new(
                pos.iter().map(|p| (*p.0,p.1.map0)),0,full_palette::RED_300.mix(0.6)
            ),

        ).unwrap().label("MAPQ: 0").legend(|(x,y)| plotters::element::Rectangle::new([(x,y),(x+10,y+5)],full_palette::RED_300));
        chart.draw_series(
            AreaSeries::new(
                pos.iter().map(|p| (*p.0,p.1.map1)),0,full_palette::YELLOW_900.mix(0.5)
            )
        ).unwrap().label("MAPQ: 1-59").legend(|(x,y)| plotters::element::Rectangle::new([(x,y),(x+10,y+5)],full_palette::YELLOW_900));
        chart.draw_series(
            AreaSeries::new(
                pos.iter().map(|p| (*p.0,p.1.map60)),0,full_palette::GREEN_300.mix(0.4)
            )
        ).unwrap().label("MAPQ: 60").legend(|(x,y)| plotters::element::Rectangle::new([(x,y),(x+10,y+5)],full_palette::GREEN_300));
        chart.draw_series(
            LineSeries::new(
                pos.iter().map(|p| (*p.0,p.1.secondary)),full_palette::BLACK
            )
        ).unwrap().label("Secondary alignments").legend(|(x,y)| PathElement::new(vec![(x,y), (x+20,y)],full_palette::BLACK));
        chart.configure_series_labels().background_style(WHITE).border_style(BLACK.mix(0.8)).draw().unwrap();
        // To avoid the IO failure being ignored silently, we manually call the present function
        root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
        println!("Result has been saved to {}", outputfile);
    }
}

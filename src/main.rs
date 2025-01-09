use bio_types::genome::AbstractInterval;
use clap::{crate_authors, Parser};
use core::num;
use csv;
use noodles_core::Region;
use noodles_fasta::{self as fasta, record::Sequence};
use num_format::{Locale, ToFormattedString};
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::prelude::{full_palette, AreaSeries, BitMapBackend, IntoDrawingArea, PathElement};
use plotters::series::{Histogram, LineSeries};
use plotters::style::*;
use plotters::{self, coord};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read};
use serde::Deserialize;
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
#[derive(Debug, Clone, PartialEq, Eq)]
struct HashMapinfo {
    map60: i64,
    map1: i64,
    map0: i64,
    overlaps: i64,
    secondary: i64,
    mismatches: i64,
    misalign: i64,
    qual: usize,
}
impl HashMapinfo {
    fn getmaxvalue(&self) -> i64 {
        assert!(std::mem::size_of::<HashMapinfo>() <= isize::MAX as _);
        let elem = unsafe { std::slice::from_raw_parts(self as *const Self as *const i64, 7) };
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
    fn new(
        map60: i64,
        map1: i64,
        map0: i64,
        secondary: i64,
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
            overlaps: 0,
            mismatches: 0,
            misalign: 0,
            qual: 0,
        }
    }
}
fn main() {
    let args = Args::parse();
    let path = &args.file;
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
    for loci in locus {
        let mut reader = match &args.index {
            Some(d) => bam::IndexedReader::from_path_and_index(path, d),
            None => bam::IndexedReader::from_path(path),
        }
        .unwrap();
        reader.set_threads(4).unwrap();
        reader
            .fetch((&loci.contig, loci.start, loci.end + 1))
            .unwrap_or_else(|_| {
                panic!(
                    "The region {}:{}-{} cannot be found, exiting.",
                    loci.contig, loci.start, loci.end
                )
            });
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
            let overlaprange =
                max(p.reference_start(), loci.start + 1)..min(loci.end, p.reference_end()); //Remove 1 because does not overlap on borders
            let mut matched = p.aligned_pairs_match().unwrap().map(|[a, b]| a..=b);
            let mut aligned = p.aligned_pairs().map(|[a, b]| a..=b);
            if p.is_secondary() {
                newrange.for_each(|i| {
                    pos.get_mut(&i).unwrap().secondary += 1;
                });
                continue;
            }
            newrange.for_each(|i| {
                match p.mapq() {
                    0 => pos.get_mut(&i).unwrap().map0 += 1,
                    1..=59 => pos.get_mut(&i).unwrap().map1 += 1,
                    60 => pos.get_mut(&i).unwrap().map60 += 1,
                    _ => panic!(
                        "MAPQ score is invalid. Got {} for {}",
                        p.mapq(),
                        String::from_utf8_lossy(p.qname())
                    ),
                };
                if i % 455 == 0 { //Check every 455 nt
                    if let Some(d) = p.aligned_pairs().find(|p| p[1]==i) {
                        let index = d[0] as usize;
                        pos.get_mut(&i).unwrap().qual += *p.qual().get(index).unwrap() as usize;
                    }
                }
                if overlaprange.contains(&i) {
                    pos.get_mut(&i).unwrap().overlaps += 1;
                }
                if matched.any(|f| f.contains(&i)) {
                    //Match skipped
                } else if aligned.any(|f| f.contains(&i)) {
                    //Aligns but not correct nt
                    pos.get_mut(&i).unwrap().mismatches += 1;
                } else {
                    //No alignment (deletion in read probably)
                    pos.get_mut(&i).unwrap().misalign += 1;
                }
            });
        }
        //writer.flush().unwrap();
        readgraph(outputdir, &loci, &pos, &args);
        // Second graph with mismatches
        mismatchgraph(outputdir, &loci, &pos, &args);
    }
}
fn mismatchgraph(
    outputdir: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
) {
    let outputfile = outputdir.join(format!(
        "{}_{}_mismatchresult.png",
        &args.species, &loci.locus
    ));
    let outputfile = outputfile.as_path();
    let root = BitMapBackend::new(outputfile, (1200, 600)).into_drawing_area();
    let _ = root.fill(&plotters::prelude::WHITE);
    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Mismatches rate and quality over the locus {} ({})",
                loci.locus, loci.contig
            ),
            ("sans-serif", 40),
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
            pos.iter().filter_map(|p| if p.1.qual > 0 { Some((*p.0, p.1.qual / std::convert::TryInto::<usize>::try_into(p.1.map0 + p.1.map1 + p.1.map60).unwrap())) } else { None }),
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
    println!("Result has been saved to {}", outputfile.display());
}
fn readgraph(
    outputdir: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
) {
    let max = pos.values().map(|max| max.getmaxvalue()).max().unwrap() + 5;
    let outputfile = outputdir.join(format!("{}_{}_readresult.png", &args.species, &loci.locus));
    let outputfile = outputfile.as_path();
    let root = BitMapBackend::new(outputfile, (1200, 600)).into_drawing_area();
    let _ = root.fill(&plotters::prelude::WHITE);
    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Reads mapping quality over the locus {} ({})",
                loci.locus, loci.contig
            ),
            ("sans-serif", 40),
        )
        .build_cartesian_2d(loci.start..loci.end, 0..max)
        .unwrap();
    let _ = chart
        .configure_mesh()
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        .x_desc("Genomic position (bp)")
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
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", outputfile.display());
}

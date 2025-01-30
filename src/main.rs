use bio_types::genome::AbstractInterval;
use clap::{crate_authors, Parser};
use plotters::coord::Shift;
use serde::ser::SerializeMap;
use core::num;
use csv;
use noodles_core::Region;
//use noodles_fasta::{self as fasta, record::Sequence};
use num_format::{Locale, ToFormattedString};
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::prelude::{full_palette, AreaSeries, BitMapBackend, DrawingArea, DrawingBackend, IntoDrawingArea, IntoSegmentedCoord, PathElement, SVGBackend};
use plotters::series::{Histogram, LineSeries};
use plotters::style::*;
use plotters::{self, coord};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, IndexedReader, Read};
use serde::{de, Deserialize, Serialize};
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
    /// Minimal number of reads (included) to declare a break in coverage
    #[arg(short,long, default_value_t=3)]
    breaks: u32,
    /// Coverage to calculate on CSV
    #[arg(short,long, default_value_t=10)]
    coverage: u32,
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
    Minus
}
impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de> {
            let s: &str = de::Deserialize::deserialize(deserializer)?;
            
                match s {
                    "1" => Ok(Strand::Minus),
                    "0" => Ok(Strand::Plus),
                    _ => Err(de::Error::unknown_variant(s, &["1","0"])),
                }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct GeneInfosFinish {
    gene: String,
    chromosome: String,
    strand: Strand,
    start: i64,
    end: i64,
    length: i64,
    reads: usize,
    matchpos: String,
    reads100: usize,
    reads100m: usize,
    coverage10x: usize
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
        let mut nocount=true;
        //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
        //let file = File::create(&filename).unwrap();
        //let mut writer = BufWriter::new(file);
        let mut pos: BTreeMap<i64, HashMapinfo> = BTreeMap::new();
        (loci.start..=loci.end).for_each(|p| {
            pos.insert(p, HashMapinfo::default());
        });
        for p in reader.rc_records() {
            nocount=false;
            let p = p.unwrap();
            let newrange =
                max(p.reference_start(), loci.start)..min(loci.end + 1, p.reference_end());
            if p.is_secondary() {
                newrange.for_each(|i| {
                    let targeting = pos.get_mut(&i).unwrap();
                    targeting.secondary += 1;
                });
                continue;
            }
            let overlaprange =
                max(p.reference_start()+1, loci.start)..min(loci.end + 1, p.reference_end() - 1); //Remove 1 because does not overlap on borders
            let mut matched = p.aligned_pairs_match().expect("No CIGAR = given.").map(|[a, b]| a..=b);
            let mut aligned = p.aligned_pairs().map(|[a, b]| a..=b);
            newrange.for_each(|i| {
                let targeting = pos.get_mut(&i).unwrap();
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
                if i % 455 == 0 {
                    //Check every 455 nt
                    if let Some(d) = p.aligned_pairs().find(|p| p[1] == i) {
                        let index = d[0] as usize;
                        targeting.qual += *p.qual().get(index).unwrap() as usize;
                    }
                }
                if overlaprange.contains(&i) {
                    targeting.overlaps += 1;
                }
                if matched.any(|f| f.contains(&i)) {
                    //Match skipped
                } else if aligned.any(|f| f.contains(&i)) {
                    //Aligns but not correct nt
                    targeting.mismatches += 1;
                } else {
                    //No alignment (deletion in read probably)
                    targeting.misalign += 1;
                }
            });
        }
        if nocount {
            panic!(
                "The region {}:{}-{} cannot be found, exiting.",
                loci.contig, loci.start, loci.end
            )
        }
        // Second graph with reads
        let fgraph = "readresult";
        if args.svg {
            let outputfile = outputdir.join(givename(
                &args.species, &loci.locus,&loci.contig,&format!("{}.svg",fgraph)
            ));
            let outputfile = outputfile.as_path();
            let root = SVGBackend::new(outputfile, (900, 500)).into_drawing_area();
            readgraph(outputfile, &loci, &pos, &args, root);
        } else {
            let outputfile = outputdir.join(givename(
                &args.species, &loci.locus,&loci.contig,&format!("{}.png",fgraph)
            ));
            let outputfile = outputfile.as_path();
            let root = BitMapBackend::new(outputfile, (1800, 1000)).into_drawing_area();
            readgraph(outputfile, &loci, &pos, &args, root);
        }
        // Second graph with mismatches
        let sgraph = "mismatchresult";
        if args.svg {
            let outputfile = outputdir.join(givename(
                &args.species, &loci.locus,&loci.contig,&format!("{}.svg",sgraph)
            ));
            let outputfile = outputfile.as_path();
            let root = SVGBackend::new(outputfile, (1200, 600)).into_drawing_area();
            mismatchgraph(outputfile, &loci, &pos, &args,root);
        } else {
            let outputfile = outputdir.join(givename(
                &args.species, &loci.locus,&loci.contig,&format!("{}.png",sgraph)
            ));
            let outputfile = outputfile.as_path();
            let root = BitMapBackend::new(outputfile, (1200, 600)).into_drawing_area();
            mismatchgraph(outputfile, &loci, &pos, &args, root);
        }
        //Create CSV from HashMap
        createcsv(outputdir, &loci, &pos, &args);
        //Create gene CSV
        if args.geneloc.is_some() {
            genelist(outputdir,&loci,reader.rc_records(),&args);
        }
    }
}
fn genelist(outputdir: &std::path::Path,
    loci: &LocusInfos,
    records: bam::RcRecords<'_, IndexedReader>,
    args: &Args) {
        let records: Vec<Result<std::rc::Rc<rust_htslib::bam::Record>, rust_htslib::errors::Error>> = records.collect();
        if records.is_empty() {
            panic!("Error with records1");
        }
        let records: Vec<std::rc::Rc<rust_htslib::bam::Record>> = records.into_iter().filter_map(Result::ok).collect();
        if records.is_empty() {
            panic!("Error with records");
        }
        let outputfile = outputdir.join(givename(
            &args.species, &loci.locus,&loci.contig,"geneanalysis.csv"
        ));
        let mut csv = csv::ReaderBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b',')
        .from_path(&args.geneloc.as_ref().unwrap())
        .unwrap();
    let mut genes: Vec<GeneInfos> = Vec::new();
    for record in csv.deserialize() {
        let record = record.expect("Invalid CSV format, waiting gene,chromosome,strand,start,end");
        genes.push(record);
    }
    if genes.is_empty() {
        panic!("Invalid CSV format, waiting gene,chromosome,strand,start,end");
    }
    let mut finale: Vec<GeneInfosFinish> = Vec::with_capacity(genes.len());
    for mut gene in genes {
        let (mut reads, mut reads100,mut reads100m) = (0,0,0);
        if gene.start > gene.end {
            (gene.end,gene.start) = (gene.start,gene.end) //Swap position
        }
        let mut hash: BTreeMap<i64,(usize,usize)> = BTreeMap::new(); //Match and full match
        let range = ranges::Ranges::from(vec![gene.start..=gene.end]);
        range.clone().into_iter().for_each(|p| { hash.insert(p, (0,0)); });
        let records = records.iter().filter(|p| { 
            let firstrange = ranges::Ranges::from(vec![p.reference_start()+1..p.reference_end()+1]);
            !firstrange.intersect(range.clone()).is_empty()
        });
        if records.clone().count() == 0 {
            panic!("Empty records");
        }
        for record in records {
            reads += 1;
            'outer: for [start,end] in record.aligned_blocks() {
                for p in start..end {
                    match hash.get_mut(&p) {
                        Some((d,_)) => { *d +=1 },
                        None => { if start > gene.end {
                            break 'outer;
                        } } //Outside coverage of gene
                    }
                }
            }
            'outer: for [start,end] in record.aligned_blocks_match().unwrap() {
                for p in start..end {
                    match hash.get_mut(&p) {
                        Some((d,_)) => { *d +=1 },
                        None => { if start > gene.end {
                            break 'outer;
                        } } //Outside coverage of gene
                    }
                }
            }
            let range= gene.start..=gene.end;
            if record.aligned_blocks().any(|p| p[0] <= *range.start() && p[1] >= *range.end()) {
                reads100 += 1;
            }
            if record.aligned_blocks_match().unwrap().any(|p| p[0] <= *range.start() && p[1] >= *range.end()) {
                reads100m += 1;
            }
        }
        let coverage = hash.clone().into_values().filter(|p| p.0 >= args.coverage.try_into().unwrap()).count();
        let text = hash.into_values().fold(String::new(), |mut acc, f| {
                acc.push_str(&format!("{}({})-",f.0,f.1));
                acc
        });
        let text = String::from(text.trim_end_matches('-'));
        let elem = GeneInfosFinish {
            gene: gene.gene,
            chromosome: gene.chromosome,
            strand: gene.strand,
            start: gene.start,
            end: gene.end,
            length: gene.end.checked_sub(gene.start).unwrap().checked_add(1).unwrap(),
            reads,
            matchpos: text,
            reads100,
            reads100m,
            coverage10x: coverage
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
fn givename(species: &str, locus: &Locus, contig: &str, suffix: &str) -> String {
    format!("{}_{}_{}_{}",species,locus,contig,suffix)
}
fn createcsv(outputdir: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
) {
    let outputfile = outputdir.join(givename(
        &args.species, &loci.locus,&loci.contig,"positionresult.csv"
    ));
    let outputfile = outputfile.as_path();
    let mut csv = csv::WriterBuilder::new().comment(Some(b'#')).flexible(false).has_headers(false).delimiter(b'\t').from_path(outputfile).unwrap();
    csv.write_record(["Position","map60","map1","map0","overlaps","secondary","mismatches","misalign","qual"]).unwrap();
    for (pos,record) in pos.iter() {
        csv.write_field(format!("{}",pos)).unwrap();
        csv.serialize(record).unwrap();
    }
    csv.flush().unwrap();
    println!("CSV analysis has been saved to {}", outputfile.display());
}
fn mismatchgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
    root: DrawingArea<T, Shift>
) where T: DrawingBackend {
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
    println!("Result has been saved to {}", outputfile.display());
}
fn readgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<i64, HashMapinfo>,
    args: &Args,
    root: DrawingArea<T, Shift>
) where T: DrawingBackend {
    let max = pos.values().map(|max| max.getmaxvalue()).max().unwrap() + 5;
    let _ = root.fill(&plotters::prelude::WHITE);
    let (top, bottom) = root.split_vertically((80).percent_height()); 
    let mut chart = ChartBuilder::on(&top)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Reads mapping quality over the locus {} ({})",
                loci.locus, loci.contig
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
    let breaks = pos.iter().filter_map(|(pos,elem)| {
        if elem.overlaps <= args.breaks.into() {
            Some((*pos,1))
        } else {
            None
        }
    });
    chart
        .draw_series(Histogram::vertical(&chart)
        .style(full_palette::RED.filled())
        .data(breaks)
        .margin(0)).unwrap()
        .label("Coverage break")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
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

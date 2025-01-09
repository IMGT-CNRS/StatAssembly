use bio_types::genome::AbstractInterval;
use clap::{crate_authors, Parser};
use csv;
use rust_htslib::bam::ext::BamRecordExtensions;
use noodles_fasta::{self as fasta, record::Sequence};
use rust_htslib::bam::{self, Read};
use serde::Deserialize;
use std::{
    fmt::Display, fs::File, io::{BufReader, BufWriter, Write}, ops::{Range, RangeInclusive}, path::PathBuf
};
use noodles_core::Region;
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
        /* if reader.records().count() == 0 {
            eprintln!("Locus {} on {}:{}-{} has no entry, skipped.",loci.locus,loci.contig,loci.start,loci.end);
            continue;
        } */
        let filename = outputdir.join(format!("{}.pileup", &loci.locus));
        let file = File::create(&filename).unwrap();
        let mut writer = BufWriter::new(file);
        /* if reader.pileup().count() == 0 {
            eprintln!("Locus {} on {}:{}-{} has no pileups, skipped.",loci.locus,loci.contig,loci.start,loci.end);
            continue;
        } */
        /* let pileups = reader.pileup().filter_map(Result::ok).filter(|f| f.alignments().any(|f| {
            f.record().contig() == loci.contig && f.record().reference_start() == loci.start && f.record().reference_end() == loci.end
        })); */
        let mut fasta = fasta::io::indexed_reader::Builder::default().build_from_path(&args.assembly).unwrap();
        for p in reader.pileup() {
            let p = p.unwrap();
            let elem = noodles_core::Position::new(p.pos().try_into().unwrap()).unwrap();
            let region = Region::new::<std::string::String,RangeInclusive<noodles_core::Position>>(loci.contig.clone(),elem..=elem);
            let fastachar = fasta.query(&region).unwrap();
            let seq = String::from_utf8_lossy(fastachar.sequence().as_ref());
            write!(writer, "{}\t{}\t{}\t{}\t", p.tid(), p.pos(), p.depth(),seq).unwrap();
            let (mut score, total) = (0,p.alignments().len());
            p.alignments().for_each(|f| {
                if !f.is_del() && !f.is_refskip() {
                    let byte = [f.record().seq()[f.qpos().unwrap()]];
                    let char = String::from_utf8_lossy(&byte);
                    if seq.to_ascii_lowercase() == char.to_ascii_lowercase() {
                        score += 1;
                    }
                    write!(
                        writer,
                        "{}",
                        char
                    )
                    .unwrap();
                } else {
                    write!(
                        writer,
                        "*"
                    ).unwrap();
                };
                match f.indel() {
                    bam::pileup::Indel::Ins(len) => write!(writer, "+{}", len).unwrap(),
                    bam::pileup::Indel::Del(len) => write!(writer, "-{}", len).unwrap(),
                    bam::pileup::Indel::None => (),
                }
            });
            write!(writer, "\t").unwrap();
            let score = score / total * 100;
            write!(writer, "{}\t",score).unwrap();
            p.alignments().for_each(|f| {
                if !f.is_del() {
                    write!(writer, "{}-", f.record().qual()[f.qpos().unwrap()]).unwrap();
                } else {
                    write!(writer, "-").unwrap();
                }
            });
            write!(writer, "\t").unwrap();
            p.alignments().for_each(|f| {
                write!(writer, "{}-", f.record().mapq()).unwrap();
            });
            writeln!(writer).unwrap();
        }
        writer.flush().unwrap();
        println!("Finished, file created {}", &filename.display());
    }
}

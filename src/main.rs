/*
This software allows the analysis of BAM files to identify the confidence on a locus (specifically IG and TR) as well as allele confidence.
It was created and used by IMGT Team (https://www.imgt.org).
Available under EUPL license
*/
///Assess quality of an assembly based on reads mapping
use clap::Parser;
use colors::full_palette::GREY_400;
use itertools::Itertools;
use plotters::coord::Shift;
use std::io::{stderr, stdout};
use std::num::NonZero;
use std::ops::RangeInclusive;
use std::process::ExitCode;
use std::time::Instant;
//use noodles_fasta::{self as fasta, record::Sequence};
use crate::r#struct::*;
use extended_htslib::bam::ext::{BamRecordExtensions, CsValue, IterAlignedPairs};
use extended_htslib::bam::{self, IndexedReader, Read};
use lazy_static::lazy_static;
use num_format::{Locale, ToFormattedString};
use plotters::chart::{ChartBuilder, LabelAreaPosition};
use plotters::prelude::{
    AreaSeries, BitMapBackend, DrawingArea, DrawingBackend, IntoDrawingArea, IntoSegmentedCoord,
    PathElement, SVGBackend, full_palette,
};
use plotters::series::{Histogram, LineSeries};
use plotters::style::*;
use std::{
    cmp::{max, min},
    collections::BTreeMap,
    fs::File,
    io::Write,
    path::PathBuf,
};
mod r#struct;
const VERSION: &str = env!("CARGO_PKG_VERSION");
const AUTHOR: &str = "IMGT";
const GLOBALMISMATCHFLOATING: usize = 10_000;
const ALERTLOCUSSIZE: i64 = 10_000_000;
lazy_static! {
    static ref fontstyle: (&'static str, u32, &'static RGBColor) = {
        let args = Args::parse();
        ("sans-serif", args.fontlegendsize, &BLACK)
    };
    static ref NAME: String = env!("CARGO_PKG_NAME").replacen('_', "/", 1);
    static ref regexpword: regex::Regex = regex::Regex::new(r"[^-\w()]").unwrap();
}
//Return block of positions thanks to CS/MD tag or CIGAR = (preferred if existing)
fn iterblock(record: &bam::Record) -> Option<Vec<[i64; 2]>> {
    match (record.getcsaligned(), record.aligned_blocks_match()) {
        //There is a CIGAR =
        (_, Some(d)) => Some(d.collect()),
        //There is a MD/CS tag
        (Some(d), None) => Some(
            d.into_iter()
                .filter_map(|p| {
                    if let CsValue::Same(d) = p.state {
                        let pos = p.getgenomepos().unwrap();
                        Some([
                            pos,
                            pos.checked_add(d.try_into().unwrap())
                                .unwrap()
                                .checked_sub(1)
                                .unwrap(),
                        ])
                    } else {
                        None
                    }
                })
                .collect(),
        ),
        //We have nothing
        (None, None) => None,
    }
}
#[allow(clippy::type_complexity)]
fn iteralert(
    args: &Args,
    mut message: bool,
    record: &bam::Record,
) -> (
    bool,
    Option<Vec<RangeInclusive<i64>>>,
    Vec<RangeInclusive<i64>>,
) {
    let aligned: Vec<RangeInclusive<i64>> = record.aligned_blocks().map(|[a, b]| a..=b).collect();
    match iterblock(record) {
        Some(a) => {
            if args.force && !message {
                eprintln!(
                    "Force used but = CIGAR or MD/CS given. Remove force to have full results. Ctrl+C to quit or wait to continue."
                );
                message = true;
                std::thread::sleep(std::time::Duration::new(5, 0));
            }
            (
                message,
                Some(a.into_iter().map(|[a, b]| a..=b).collect()),
                aligned,
            )
        }
        //Check if forced or not, if yes, force software
        None => {
            let text = "No = CIGAR given";
            if !args.force {
                eprintln!(
                    "{}. Add --force to force even without = or MD/CS tag (some results won't be available).",
                    text
                );
                return (false, None, aligned);
            } else if !message {
                eprintln!("{} but it was forced... Continuing...", text);
                message = true;
            }
            (
                message,
                Some(
                    std::iter::empty::<IterAlignedPairs>()
                        .map(|_| 0..=0)
                        .collect(),
                ),
                aligned,
            )
        }
    }
}
//Filter reads thanks to args provided, remove reverse or supplementary/secondary alignments
fn filterread(args: &Args, record: &bam::Record) -> bool {
    if args.forward && record.is_reverse() {
        return false;
    }
    if !args.allreads && (record.is_supplementary() || record.is_secondary()) {
        return false;
    }
    true
}
//Check we can read BAM file and return the reader with desired threads
fn getreaderoffile(args: &Args) -> Result<IndexedReader, extended_htslib::errors::Error> {
    let mut reader = match &args.index {
        Some(d) => bam::IndexedReader::from_path_and_index(&args.file, d),
        None => bam::IndexedReader::from_path(&args.file),
    }?;
    //Set threads if given or thanks to parralelism. If nothing available, set to 4 by default.
    let threads = match (
        NonZero::new(args.threads),
        std::thread::available_parallelism(),
    ) {
        (Some(d), _) => d,
        (_, Ok(d)) => min(d, NonZero::new(12).unwrap()),
        _ => NonZero::new(4).unwrap(),
    };
    reader.set_threads(threads.get()).unwrap();
    Ok(reader)
}
//Check there is one alternate for one primary.
fn mergelocus(locus: Vec<LocusInfos>) -> Option<Vec<Vec<LocusInfos>>> {
    let mut elem: Vec<Vec<LocusInfos>> = Vec::with_capacity(locus.len());
    let mut alternate = false;
    let mut actual: Vec<LocusInfos> = Vec::new();
    for loci in locus {
        match elem
            .iter()
            .find(|p| p.iter().any(|f| f.locus == loci.locus))
        {
            Some(d) if d != elem.last().unwrap() => {
                eprintln!("Locus {} is splited! Aborted.", loci.locus);
                return None;
            }
            _ => (),
        };
        if !loci.haplotype.isprimary() {
            if actual.is_empty() {
                eprintln!("Empty without a corresponding primary!");
                return None;
            } else {
                alternate = true;
            }
        }
        match actual.first() {
            Some(e) if e.locus == loci.locus && alternate && actual.len() >= 2 => {
                eprintln!("Only one alternate is allowed!");
                return None;
            }
            Some(e) if e.locus == loci.locus && alternate => actual.push(loci),
            Some(e) if e.locus != loci.locus && !loci.haplotype.isprimary() => {
                eprintln!("Alternate without a corresponding primary!");
                return None;
            }
            _ => {
                if !actual.is_empty() {
                    elem.push(actual.clone())
                };
                actual.clear();
                actual.push(loci);
                alternate = false;
            }
        }
    }
    if !actual.is_empty() {
        elem.push(actual)
    };
    Some(elem)
}
//Get number of mismatches for the record (x 10_000 to get as an integer)
fn getglobalmismatch(args: &Args, record: &bam::Record) -> usize {
    let length = if record.seq_len() != 0 {
        record.seq_len()
    } else {
        1
    };
    match (args.totalread, record.aux(b"NM")) {
        (true, Ok(extended_htslib::bam::record::Aux::U8(d))) => {
            if d == 0 {
                0
            } else {
                (d as usize * GLOBALMISMATCHFLOATING) / length
            }
        }
        _ => 0,
    }
}
//Parse the location csv with locus infos
fn locusposparser(args: &Args) -> std::io::Result<Vec<LocusInfos>> {
    let mut csv = match csv::ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(&args.locuspos)
    {
        Ok(c) => c,
        Err(e) => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CSV file of location position cannot be found. Error is {}",
                    e
                ),
            ));
        }
    };
    let mut locus: Vec<LocusInfos> = Vec::new();
    for record in csv.deserialize() {
        let record = match record {
            Ok(r) => r,
            Err(e) => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Invalid CSV format (tabular fields) required for file {}, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend\n{}",
                        args.locuspos.display(),
                        e
                    ),
                ));
            }
        };
        locus.push(record);
    }
    if locus.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Empty CSV format for file {}, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend. Nothing found.",
                args.locuspos.display()
            ),
        ));
    }
    //At least one duplicate line
    if let Some(d) = locus.iter().duplicates().next() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "The locus {} ({}) on region {}:{}-{} appears more than once. Please provide a unique gene name.",
                d.locus,
                d.haplotype,
                d.contig,
                d.start.getobasedpos(),
                d.end.getobasedpos()
            ),
        ));
    }
    //make complement if locus is complement
    locus.iter_mut().for_each(|r| {
        if r.start >= r.end {
            (r.end, r.start) = (r.start.clone(), r.end.clone());
            r.complement = true;
        }
    });
    if !args.hugeregion {
        if let Some(big) = locus
            .iter()
            .find(|p| p.end.length(&p.start) >= ALERTLOCUSSIZE)
        {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "The region {}-{} is more than {}. Check your input file ({}) is correct.n\nIf wanted, please add --hugeregion parameters. Be careful software might use a lot of memory",
                    big.start.getobasedpos(),
                    big.end.getobasedpos(),
                    ALERTLOCUSSIZE,
                    args.locuspos.display()
                ),
            ));
        }
    }
    Ok(locus)
}
//Check BAM file exists and outputdir is created and return it
fn checklocusandoutput(args: &Args) -> std::io::Result<&PathBuf> {
    //Check bam file exists
    if let Err(e) = getreaderoffile(args) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!(
                "Cannot read bam file ({}). Error is: {}. Exiting.",
                args.file.display(),
                e
            ),
        ));
    }
    let outputdir = match args.outdir.is_dir() {
        true => {
            eprintln!(
                "Output folder exists: {}. Will be overwritten.",
                &args.outdir.display()
            );
            &args.outdir
        }
        false => {
            eprintln!(
                "Folder {} does not exist, attempt to create.",
                args.outdir.display()
            );
            std::fs::create_dir(&args.outdir)?;
            &args.outdir
        }
    };
    Ok(outputdir)
}
//Process countings for each read
fn processcounting(
    args: &Args,
    pos: &mut BTreeMap<Position, HashMapinfo>,
    newrange: std::ops::Range<Position>,
    record: &bam::Record,
    sep: i64,
    matched: &[RangeInclusive<i64>],
    aligned: &[RangeInclusive<i64>],
) {
    for (i, targeting) in pos.range_mut(newrange) {
        let i = &i.getzbasedpos();
        let _time = Instant::now();
        targeting.globalmismatch += getglobalmismatch(args, record);
        match record.mapq() {
            0 => targeting.map0 += 1,
            1..=59 => targeting.map1 += 1,
            60 => targeting.map60 += 1,
            _ => {
                eprintln!(
                    "MAPQ score is invalid. Got {} for {}",
                    record.mapq(),
                    String::from_utf8_lossy(record.qname())
                );
                return;
            }
        };
        if i % sep == 0 {
            //Get quality of reads depending on genomic position
            if let Some(d) = record.aligned_pairs().find(|p| p[1] == *i) {
                let index = d[0] as usize;
                targeting.qual += *record.qual().get(index).unwrap() as usize;
            }
        }
        //If not a match, add to mismatches or misalign (if none is found)
        if !matched
            .iter()
            .skip_while(|p| p.end() <= i)
            .take_while(|p| p.start() <= i)
            .any(|s| s.contains(i))
        {
            if aligned
                .iter()
                .skip_while(|p| p.end() <= i)
                .take_while(|p| p.start() <= i)
                .any(|s| s.contains(i))
            {
                //Aligns but not correct nt
                if !args.force {
                    targeting.mismatches += 1;
                };
            } else {
                //No alignment (deletion in read probably)
                targeting.misalign += 1;
            }
        }
        //Overlap reads
        if record.reference_start() != *i && record.reference_end() != *i {
            targeting.overlaps += 1;
        }
    }
}
fn main() -> ExitCode {
    let firstinstant = Instant::now();
    let args = Args::parse();
    if args.percentalerting >= args.percentwarning {
        eprintln!("Percent warning must be greater than percent alerting.");
        return ExitCode::FAILURE;
    }
    //Get locus, geneloc and outputdir, print errors if we have
    let (outputdir, locus) = match (checklocusandoutput(&args), locusposparser(&args)) {
        (Err(e), _) => {
            eprintln!("{}", e);
            return ExitCode::FAILURE;
        }
        (_, Err(f)) => {
            eprintln!("{}", f);
            return ExitCode::FAILURE;
        }
        (Ok(a), Ok(b)) => (a, b),
    };
    //Group between primary and alternate
    let grouped = match mergelocus(locus) {
        Some(g) => g,
        None => {
            eprintln!(
                "Check order of loci in the file {}.",
                args.locuspos.display()
            );
            return ExitCode::FAILURE;
        }
    };
    for locus in grouped {
        let floci = locus.first().unwrap();
        let haplotype = locus.len();
        if haplotype > 2 {
            eprintln!("There is more than 2 haplotypes for {}", floci.locus);
            return ExitCode::FAILURE;
        }
        let haplotypebool = haplotype == 1; //IS there one or two haplotypes?
        println!(
            "Going for {} locus - {}",
            floci.locus,
            if !haplotypebool { "diploid" } else { "haploid" }
        );
        let readresultsize = (
            1800,
            1000 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
        );
        let mismatchsize = (
            1200,
            600 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap(),
        );
        //Get infos for graph
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
                let root = SVGBackend::new(&outputfile1, readresultsize).into_drawing_area();
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
                let root = BitMapBackend::new(&outputfile2, readresultsize).into_drawing_area();
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
                    haplotypebool,
                    &format!("{}.svg", sgraph),
                    true,
                ));
                outputfile3 = outputfile;
                let root = SVGBackend::new(&outputfile3, mismatchsize).into_drawing_area();
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
                    haplotypebool,
                    &format!("{}.png", sgraph),
                    true,
                ));
                outputfile4 = outputfile;
                let root = BitMapBackend::new(&outputfile4, mismatchsize).into_drawing_area();
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
        //For each individual haplotype inside locus
        for loci in locus.iter() {
            let mut reader = match getreaderoffile(&args) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("Cannot read bam file. Error is {}. Exiting", e);
                    return ExitCode::FAILURE;
                }
            };
            //0-based except end because end is exclusive
            if reader
                .fetch((
                    &loci.contig,
                    loci.start.getzbasedpos(),
                    loci.end.getobasedpos(),
                ))
                .is_err()
            {
                eprintln!(
                    "The region {}:{}-{} cannot be found, exiting.",
                    loci.contig,
                    loci.start.getobasedpos(),
                    loci.end.getobasedpos()
                );
                return ExitCode::FAILURE;
            };
            let mut nocount = true;
            //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
            //let file = File::create(&filename).unwrap();
            //let mut writer = BufWriter::new(file);
            let locusrange = loci.start.getzbasedpos()..=loci.end.getzbasedpos();
            let mut pos: BTreeMap<Position, HashMapinfo> = BTreeMap::new();
            //Populate all B-Tree position, 0-based, invert if locus complement
            let iterator: Vec<(usize, i64)> = if loci.complement {
                locusrange.rev().enumerate().collect()
            } else {
                locusrange.enumerate().collect()
            };
            iterator.into_iter().for_each(|(i, p)| {
                let locuspos: i64 = i.try_into().unwrap();
                let default = HashMapinfo {
                    locuspos: Position::new(true, locuspos),
                    position: Position::new(true, p),
                    ..Default::default()
                };
                pos.insert(Position::new(true, p), default);
            });
            let mut message = false;
            println!("Region {} fetched, analyzing all reads.", loci.locus);
            let mut count = 0;
            let time = Instant::now();
            let sep = if args.fullquality {
                1
            } else {
                max(
                    (loci.end.getobasedpos() - loci.start.getobasedpos() + 1) / 250,
                    100,
                ) //250 points for quality point or 100nt break
            };
            for p in reader
                .rc_records()
                .filter_map(Result::ok)
                .filter(|p| !(args.forward && p.is_reverse()))
            {
                count += 1;
                //Print every 100 reads done
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
                //Get range to put the reads inclusive pos
                let newrange =
                    Position::new(true, max(p.reference_start(), loci.start.getzbasedpos()))
                        ..Position::new(true, min(loci.end.getobasedpos(), p.reference_end()));
                if p.is_secondary() || p.is_supplementary() {
                    for (_, targeting) in pos.range_mut(newrange) {
                        if p.is_secondary() {
                            targeting.secondary += 1;
                        } else {
                            targeting.supplementary += 1;
                            /* targeting
                            .supplementary
                            .push(String::from_utf8_lossy(p.qname()).to_string()); */
                        }
                    }
                    continue;
                }
                let (matched, aligned) = match iteralert(&args, message, &p) {
                    (_, None, _) => {
                        return ExitCode::FAILURE;
                    } //Kill software, errors sent by iteralert
                    (newmessage, Some(p), aligned) => {
                        message = newmessage;
                        (p, aligned)
                    }
                };
                processcounting(&args, &mut pos, newrange, &p, sep, &matched, &aligned);
            }
            if nocount {
                eprintln!(
                    "The region {}:{}-{} cannot be found, exiting.",
                    loci.contig,
                    loci.start.getobasedpos(),
                    loci.end.getobasedpos()
                );
                return ExitCode::FAILURE;
            }
            //Quality and mismatch is the sum of reads so dividing to get real results
            pos.iter_mut().for_each(|(_, p)| {
                if p.qual > 0 {
                    p.qual /= max(
                        std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap(),
                        1,
                    )
                }
                p.globalmismatch /= max(
                    std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap(),
                    1,
                );
            });
            println!("Making graphs");
            match (loci.haplotype.isprimary(), args.svg) {
                (true, true) => {
                    let f = readgraph(
                        &outputfile1,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        readgraphtop.clone().unwrap(),
                    );
                    if let Err(e) = f {
                        eprintln!("Cannot create read graph. Error is {}", e);
                        return ExitCode::FAILURE;
                    }
                    mismatchgraph(
                        &outputfile3,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        mismatchgraphtop.clone().unwrap(),
                    );
                }
                (false, true) => {
                    let f = readgraph(
                        &outputfile1,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        readgraphbottom.clone().unwrap(),
                    );
                    if let Err(e) = f {
                        eprintln!("Cannot create read graph. Error is {}", e);
                        return ExitCode::FAILURE;
                    }
                    mismatchgraph(
                        &outputfile3,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        mismatchgraphbottom.clone().unwrap(),
                    );
                }
                (true, false) => {
                    let f = readgraph(
                        &outputfile2,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        readgraphtop2.clone().unwrap(),
                    );
                    if let Err(e) = f {
                        eprintln!("Cannot create read graph. Error is {}", e);
                        return ExitCode::FAILURE;
                    }
                    mismatchgraph(
                        &outputfile4,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        mismatchgraphtop2.clone().unwrap(),
                    );
                }
                (false, false) => {
                    let f = readgraph(
                        &outputfile2,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        readgraphbottom2.clone().unwrap(),
                    );
                    if let Err(e) = f {
                        eprintln!("Cannot create read graph. Error is {}", e);
                        return ExitCode::FAILURE;
                    }
                    mismatchgraph(
                        &outputfile4,
                        loci,
                        pos.values().collect_vec().as_slice(),
                        &args,
                        mismatchgraphbottom2.clone().unwrap(),
                    );
                }
            }
            println!("Graphs finished.");
            //Create CSV from HashMap
            if let Err(e) = createcsv(
                outputdir,
                loci,
                pos.values().collect_vec().as_slice(),
                &args,
            ) {
                eprintln!("Cannot create csv file. Error is {}", e);
                return ExitCode::FAILURE;
            }
            //Create gene CSV
            if args.geneloc.is_some() {
                println!("Gene list starting!");
                if let Err(e) = genelist(outputdir, loci, &args) {
                    eprintln!("Cannot create gene list. Error is {}", e);
                    return ExitCode::FAILURE;
                }
                println!("Gene list finished.");
            } else {
                eprintln!(
                    "You have not provided a gene list, skipped. Provide one to get more datas."
                );
            }
        }
        println!("Locus {} is done!", &floci.locus);
    }
    println!(
        "{} done in {} s",
        NAME.as_str(),
        firstinstant.elapsed().as_secs_f32()
    );
    ExitCode::SUCCESS
}
fn extractgenelist(
    args: &Args,
    loci: &LocusInfos,
) -> Result<Vec<GeneInfos>, Box<dyn std::error::Error>> {
    let geneloc = match &args.geneloc {
        Some(l) => l,
        None => {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "No localization given",
            )));
        }
    };
    let mut genes: Vec<GeneInfos> = Vec::new();
    let mut csv = csv::ReaderBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b',')
        .from_path(geneloc)?;
    for record in csv.deserialize() {
        let record = match record {
            Ok(r) => r,
            Err(e) => {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Invalid CSV format (comma separated) for file {}, waiting gene,chromosome,strand,start,end case sensitive. Have you kept the header?\n{}",
                        geneloc.display(),
                        e
                    ),
                )));
            }
        };
        genes.push(record);
    }
    if genes.is_empty() {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Empty CSV format for file {}, waiting gene,chromosome,strand,start,end case sensitive. Nothing found.",
                geneloc.display()
            ),
        )));
    }
    //Invert if start >= end because strand is given and won't work
    genes.iter_mut().filter(|p| p.start > p.end).for_each(|p| {
        (p.start, p.end) = (p.end.clone(), p.start.clone());
    });
    //Retain genes inside the correct loci
    genes.retain(|gene| {
        gene.chromosome == loci.contig
            && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                .contains(&gene.start.getobasedpos())
            && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                .contains(&gene.end.getobasedpos())
    });
    if genes.is_empty() {
        println!("No gene identified for locus {}, skipped.", loci.locus);
        return Ok(Vec::new());
    }
    //At least one duplicate line
    if let Some(d) = genes
        .iter()
        .duplicates_by(|g| (g.gene.as_str(), loci.haplotype.isprimary()))
        .next()
    {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "The gene {} appears more than once (same gene and haplotype). Please provide a unique gene name or consider splitting gene list.",
                d.gene
            ),
        )));
    }
    Ok(genes)
}
fn genelist(
    outputdir: &std::path::Path,
    loci: &LocusInfos,
    args: &Args,
) -> Result<(), Box<dyn std::error::Error>> {
    let genes = extractgenelist(args, loci)?;
    if genes.is_empty() {
        return Ok(());
    }
    let outputfile = outputdir.join(givename(
        &args.species,
        &loci.locus,
        &loci.contig,
        loci.haplotype.isprimary(),
        "geneanalysis.csv",
        false,
    ));
    let mut lock: std::io::StderrLock<'_> = stderr().lock();
    let mut finale: Vec<GeneInfosFinish> = Vec::with_capacity(genes.len());
    //For each gene, list of alerting positions, bbool said suspicious or warning position
    let mut alertingpositions: BTreeMap<GeneInfos, Vec<(bool, usize)>> = BTreeMap::new();
    for gene in genes {
        let mut reader = getreaderoffile(args)?;
        let (mut reads, mut readsfull, mut reads100, mut reads100m) = (0, 0, 0, 0);
        let (mut hash, records) = {
            //O position is exclusive
            let genegenericrange = gene.start.getzbasedpos()..gene.end.getobasedpos();
            //As gene start is 1-ranged, put it as 0-range with -1. End is exclusive so -1/+1 = 0
            reader.fetch((
                &gene.chromosome,
                genegenericrange.start,
                genegenericrange.end,
            ))?;
            let mut hash: BTreeMap<i64, Posread> = BTreeMap::new(); //Match and full match and total
            //Hash contains 1-based positions
            genegenericrange.for_each(|p| {
                hash.insert(p, Posread::new(0, 0, 0, args).unwrap());
            });
            (hash, reader.records())
        };
        let mut coverageperc = 0;
        let mut empty = true;
        for record in records
            .filter_map(Result::ok)
            .filter(|p| filterread(args, p))
        {
            empty = false;
            reads += 1;
            let range = record.reference_start()..record.reference_end();
            coverageperc += ranges::Ranges::from(range.clone()).into_iter().count();
            'outer: for [start, end] in record.aligned_blocks() {
                for p in start..end {
                    match hash.get_mut(&p) {
                        Some(d) => d.addindel(1),
                        None => {
                            if start > gene.end.getzbasedpos() {
                                break 'outer;
                            }
                        } //Outside coverage of gene
                    }
                }
            }
            if !args.force {
                'outer: for [start, end] in iterblock(&record).unwrap() {
                    for p in start..end {
                        match hash.get_mut(&p) {
                            Some(d) => d.addmatch(1),
                            None => {
                                if start > gene.end.getzbasedpos() {
                                    break 'outer;
                                }
                            } //Outside coverage of gene
                        }
                    }
                }
            }
            //Check 0-based gene is inside the inclusive range
            let validrange = |p: [i64; 2], genestart: &Position, geneend: &Position| {
                let range = p[0]..=p[1];
                range.contains(&genestart.getzbasedpos()) && range.contains(&geneend.getzbasedpos())
            };
            if record
                .aligned_blocks()
                .any(|p| validrange(p, &gene.start, &gene.end))
            {
                reads100 += 1;
            }
            if range.contains(&gene.start.getzbasedpos())
                && range.contains(&gene.end.getzbasedpos())
            {
                readsfull += 1;
            }
            if !args.force
                && iterblock(&record)
                    .unwrap()
                    .into_iter()
                    .any(|p| validrange(p, &gene.start, &gene.end))
            {
                reads100m += 1;
            }
            for p in range {
                match hash.get_mut(&p) {
                    Some(d) => d.addtotal(1),
                    None => {
                        if record.reference_start() > gene.end.getzbasedpos() {
                            break;
                        }
                    } //Outside coverage of gene
                }
            }
        }
        if empty {
            writeln!(lock, "Empty records for gene {}", gene.gene).unwrap();
            //PUT 0 value on the CSV
            let elem = GeneInfosFinish::make_default(gene);
            finale.push(elem);
            continue;
        }
        //Coverage calculus
        let coverage = hash
            .iter()
            .filter(|(_, p)| p.gettotal() >= args.coverage.try_into().unwrap())
            .count();
        //Reverse if complement
        let iterator: Vec<(&i64, &Posread)> = match gene.strand {
            Strand::Plus => hash.iter().collect(),
            Strand::Minus => hash.iter().rev().collect(),
        };
        //Merging data
        let text = iterator.into_iter().fold(String::new(), |mut acc, (_, f)| {
            acc.push_str(&format!(
                "{}({}={}X{}ID)-",
                f.gettotal(),
                f.getmatch(),
                f.getmismatchcount(),
                f.getindelcount()
            ));
            acc
        });
        let text = String::from(text.trim_end_matches('-'));
        let plots = outputdir
            .join(format!("gene_{}", args.species))
            .join(loci.haplotype.to_string().as_str().to_lowercase());
        if !std::fs::exists(&plots)? {
            println!("Creating the folder {}", plots.display());
            std::fs::create_dir_all(&plots)?;
        };
        let mut output = plots.join(regexpword.replace_all(&gene.gene, "_").to_uppercase());
        if !args.svg {
            output.set_extension("png");
            let root = BitMapBackend::new(&output, (700, 400)).into_drawing_area();
            //Gene graph
            genegraph(
                args,
                &hash,
                &gene,
                loci,
                root,
                &mut alertingpositions,
                reads100m,
            );
        } else {
            output.set_extension("svg");
            let root = SVGBackend::new(&output, (700, 400)).into_drawing_area();
            //Gene graph
            genegraph(
                args,
                &hash,
                &gene,
                loci,
                root,
                &mut alertingpositions,
                reads100m,
            );
        }
        let coverageperc = ((coverageperc * 1_000
            / reads
            / usize::try_from(gene.end.length(&gene.start)).unwrap())
            as f32)
            .round()
            / 1_000.0;
        let elem = GeneInfosFinish::new(
            gene,
            reads,
            readsfull,
            Some(text),
            reads100,
            reads100m,
            coverageperc,
            coverage,
        );
        finale.push(elem);
    }
    finale.sort_unstable(); //Sort the table
    let mut csv = csv::WriterBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(&outputfile)?;
    for gene in finale {
        csv.serialize(gene)?;
    }
    csv.flush()?;
    //Sorting the positions
    alertingpositions.iter_mut().for_each(|(_, vec)| {
        vec.sort_unstable_by(|a, b| match a.1.cmp(&b.1) {
            std::cmp::Ordering::Equal => a.0.cmp(&b.0),
            ord => ord,
        });
    });
    //Print position if nothing was forced
    if !args.force {
        printpossus(args, loci, outputdir, &alertingpositions)?;
    }
    println!("Gene analysis has been saved to {}", outputfile.display());
    Ok(())
}
fn printbreaks<T>(
    args: &Args,
    finalpos: i64,
    breaks: T,
    loci: &LocusInfos,
    outputfile: &std::path::Path,
) -> std::io::Result<()>
where
    T: Iterator<Item = (i64, i64)> + Clone,
{
    let mut breakfile = File::create(outputfile.parent().unwrap().join(givename(
        &args.species,
        &loci.locus,
        &loci.contig,
        loci.haplotype.isprimary(),
        "break.txt",
        false,
    )))?;
    let mut first = None;
    let mut prev = None;
    //Might be none if no breaks
    let finalbreak = match breaks.clone().map(|p| p.0).max() {
        Some(d) => d,
        None => {
            breakfile.write_all("No breaks.".as_bytes())?;
            return Ok(());
        }
    };
    let mut acc = breaks.fold(String::new(), |mut acc, (num, _)| {
        if first.is_none() {
            first = Some(num);
            prev = Some(num);
        } else if let (Some(mut prev_num), Some(f)) = (prev, first) {
            if num - prev_num != 1 || num == finalpos || num == finalbreak {
                if num == finalpos {
                    prev_num = finalpos;
                } else if num == finalbreak {
                    prev_num = finalbreak;
                }
                if f == prev_num {
                    acc.push_str(&format!("{}:{}\n", loci.contig, f));
                } else {
                    acc.push_str(&format!("{}:{}..{}\n", loci.contig, f, prev_num));
                }
                if num != finalpos && num != finalbreak {
                    first = Some(num);
                    prev = Some(num);
                } else {
                    first = None;
                }
            } else {
                prev = Some(num);
            }
        }
        acc
    });
    if let Some(d) = first {
        acc.push_str(&format!("{}:{}\n", loci.contig, d + 1));
    }
    breakfile.write_all(acc.trim().as_bytes())?;
    breakfile.flush()?;
    Ok(())
}
fn printpossus(
    args: &Args,
    loci: &LocusInfos,
    outputdir: &std::path::Path,
    data: &BTreeMap<GeneInfos, Vec<(bool, usize)>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let outputfile = outputdir.join(givename(
        &args.species,
        &loci.locus,
        &loci.contig,
        loci.haplotype.isprimary(),
        "allele_confidence.csv",
        false,
    ));
    let mut csv = csv::WriterBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(&outputfile)?;
    /* if data.is_empty() {
        csv.write_record(["No warning or alerting positions found."])?;
    } else { */
    csv.write_record(["Gene", "Positions (! for alerting, ~ for warning)"])?;
    for (gene, vec) in data {
        csv.write_field(&gene.gene)?;
        let infos = vec.iter().fold(String::new(), |mut acc, f| {
            acc.push_str(&format!("-{}({})", f.1, if f.0 { "!" } else { "~" }));
            acc
        });
        let infos = infos.trim_matches('-');
        csv.write_field(infos)?;
        csv.write_record(None::<&[u8]>)?;
    }
    //}
    csv.flush()?;
    println!(
        "Gene suspicious position has been saved to {}",
        outputfile.display()
    );
    Ok(())
}
fn genegraph<T>(
    args: &Args,
    hash: &BTreeMap<i64, Posread>,
    gene: &GeneInfos,
    loci: &LocusInfos,
    root: DrawingArea<T, Shift>,
    alerting: &mut BTreeMap<GeneInfos, Vec<(bool, usize)>>,
    reads100m: usize,
) where
    T: DrawingBackend,
{
    let genename = gene.gene.to_string();
    let text_style = fontstyle.into_text_style(&root);
    let _ = root.fill(&plotters::prelude::WHITE);
    let max = hash.values().map(|p| p.gettotal()).max().unwrap() + 5;
    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .right_y_label_area_size(40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            format!(
                "Reads coverage for {} ({}-{})",
                genename, loci.haplotype, gene.strand
            ),
            ("sans-serif", 22),
        )
        .build_cartesian_2d(1..hash.len(), 0..max)
        .unwrap();
    //Reverse complement genes
    let hash: Vec<(&i64, &Posread)> = match gene.strand {
        Strand::Plus => hash.iter().collect(),
        Strand::Minus => hash.iter().rev().collect(),
    };
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos + 1, val.gettotal())),
            full_palette::BLUE_300.mix(0.8).stroke_width(4),
        ))
        .unwrap()
        .label("Total reads")
        .legend(|(x, y)| {
            PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLUE_300.mix(0.8))
        });
    //Enumerate to get position relative to the gene (+1 because 0-related)
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos + 1, val.getindel())),
            full_palette::ORANGE_300.mix(0.7).stroke_width(3),
        ))
        .unwrap()
        .label("Sequence match")
        .legend(|(x, y)| {
            PathElement::new(vec![(x, y), (x + 15, y)], full_palette::ORANGE_300.mix(0.8))
        });
    //Three levels if not forced
    if !args.force {
        chart
            .draw_series(LineSeries::new(
                hash.iter()
                    .enumerate()
                    .map(|(pos, (_, val))| (pos + 1, val.getmatch())),
                full_palette::GREEN_400.mix(0.5).stroke_width(2),
            ))
            .unwrap()
            .label("Sequence equal")
            .legend(|(x, y)| {
                PathElement::new(vec![(x, y), (x + 15, y)], full_palette::GREEN_400.mix(0.5))
            });
        //Line of all full reads
        chart
            .draw_series(LineSeries::new(
                hash.iter()
                    .enumerate()
                    .map(|(pos, ..)| (pos + 1, reads100m)),
                BLACK.mix(0.7).stroke_width(3),
            ))
            .unwrap()
            .label("Reads 100% match")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLACK.mix(0.7)));
        chart
            .draw_series(
                Histogram::vertical(&chart)
                    .style(full_palette::ORANGE_300.mix(0.3).filled())
                    .data(hash.iter().enumerate().filter_map(|(pos, (_, val))| {
                        let pos1 = pos + 1;
                        if val.iswarning() {
                            match alerting.get_mut(gene) {
                                Some(d) => d.push((false, pos1)),
                                None => {
                                    alerting.insert(gene.clone(), vec![(false, pos1)]);
                                }
                            };
                            Some((pos1, max))
                        } else {
                            None
                        }
                    }))
                    .margin(0),
            )
            .unwrap()
            .label("Warning positions")
            .legend(|(x, y)| {
                plotters::element::Rectangle::new(
                    [(x, y), (x + 15, y + 5)],
                    full_palette::ORANGE_300.filled(),
                )
            });
        chart
            .draw_series(
                Histogram::vertical(&chart)
                    .style(full_palette::RED_400.mix(0.3).filled())
                    .data(hash.iter().enumerate().filter_map(|(pos, (_, val))| {
                        let pos1 = pos + 1;
                        if val.issuspicious() {
                            match alerting.get_mut(gene) {
                                Some(d) => d.push((true, pos1)),
                                None => {
                                    alerting.insert(gene.clone(), vec![(true, pos1)]);
                                }
                            };
                            Some((pos1, max))
                        } else {
                            None
                        }
                    }))
                    .margin(0),
            )
            .unwrap()
            .label("Suspicious positions")
            .legend(|(x, y)| {
                plotters::element::Rectangle::new(
                    [(x, y), (x + 15, y + 5)],
                    full_palette::RED_400.filled(),
                )
            });
    }
    //Continue graph
    chart
        .configure_mesh()
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        .x_desc("Position in sequence (bp)")
        .y_desc("Reads count")
        .label_style(text_style.clone())
        .light_line_style(GREY_400.mix(0.6))
        .x_max_light_lines(5)
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .unwrap();
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::LowerRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style)
            .border_style(BLACK.mix(0.8))
            .draw()
            .unwrap();
    }
    // To avoid the IO failure being ignored silently, we manually call the present function
    drawnoticetext(&root);
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
        if image && !haplo {
            String::new()
        } else {
            format!("{}-", contig)
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
    pos: &[&HashMapinfo],
    args: &Args,
) -> Result<(), Box<dyn std::error::Error>> {
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
        .has_headers(true)
        .delimiter(b'\t')
        .flexible(true)
        .from_path(outputfile)?;
    for record in pos {
        csv.serialize(record)?;
    }
    csv.flush()?;
    println!("CSV analysis has been saved to {}", outputfile.display());
    Ok(())
}
fn drawnoticetext<T>(root: &DrawingArea<T, Shift>)
where
    T: DrawingBackend,
{
    let text = format!(
        "Graph made by {} version {} ({})",
        NAME.as_str(),
        VERSION,
        AUTHOR
    );
    let text_style = ("Georgia", 11, FontStyle::Oblique, &BLACK).into_text_style(root);
    let size = root.estimate_text_size(&text, &text_style).unwrap();
    let size = (
        (root.dim_in_pixel().0 - size.0 - 10) as i32,
        (root.dim_in_pixel().1 - size.1 - 10) as i32,
    );
    root.draw_text(&text, &text_style, size).unwrap();
}
fn mismatchgraph<T>(
    _outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    args: &Args,
    root: DrawingArea<T, Shift>,
) where
    T: DrawingBackend,
{
    let _ = root.fill(&plotters::prelude::WHITE);
    let text_style = fontstyle.into_text_style(&root);
    let (top, bottom) = match args.totalread {
        true => {
            let (top, bottom) = root.split_vertically((50).percent_height());
            (top, Some(bottom))
        }
        false => (root.clone(), None),
    };
    let mut chart = ChartBuilder::on(&top)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Mismatches rate and quality on the locus {} ({}-{})",
                loci.locus, loci.contig, loci.haplotype
            ),
            ("sans-serif", 28),
        )
        .build_cartesian_2d(loci.start.getobasedpos()..loci.end.getobasedpos(), 0..100)
        .unwrap();
    if !args.force {
        chart
            .draw_series(LineSeries::new(
                pos.iter().filter_map(|p| {
                    let val: i32 = (p.mismatches * 100 / max(1, p.gettotalmap())) as i32;
                    if val > 0 && p.gettotalmap() > 0 {
                        Some((p.position.getobasedpos(), val))
                    } else {
                        None
                    }
                }),
                full_palette::DEEPPURPLE_400.mix(0.8).filled(),
            ))
            .unwrap()
            .label("Mismatches (%)")
            .legend(|(x, y)| {
                PathElement::new(vec![(x, y), (x + 15, y)], full_palette::DEEPPURPLE_400)
            });
    }
    chart
        .draw_series(LineSeries::new(
            pos.iter().filter_map(|p| {
                let val = (p.misalign * 100 / max(1, p.gettotalmap())) as i32;
                if val > 0 && p.gettotalmap() > 0 {
                    Some((p.position.getobasedpos(), val))
                } else {
                    None
                }
            }),
            full_palette::RED_400.mix(0.8).filled(),
        ))
        .unwrap()
        .label("Misalign (%)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::RED_400));
    let mut secondary = chart.set_secondary_coord(
        loci.start.getobasedpos()..loci.end.getobasedpos(),
        0..100usize,
    );
    secondary
        .configure_mesh()
        .y_label_formatter(&|f| format!("{}%", f))
        .x_label_formatter(&|f| {
            format!(
                "{} ({})",
                f.to_formatted_string(&Locale::en),
                pos.iter()
                    .find(|p| { p.position.getobasedpos() == *f })
                    .unwrap()
                    .locuspos
                    .getobasedpos()
                    .to_formatted_string(&Locale::en)
            )
        })
        .x_desc("Genomic position (bp)")
        .y_desc("Mismatch rate (%)")
        .disable_x_mesh()
        .label_style(text_style.clone())
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .unwrap();
    //let mut second = chart.set_secondary_coord(loci.start..loci.end, 0..max);
    secondary
        .draw_secondary_series(LineSeries::new(
            pos.iter().filter_map(|p| {
                if p.qual > 0 {
                    Some((p.position.getobasedpos(), p.qual))
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
    if !args.nolegend {
        secondary
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .unwrap();
    }
    //Bottom graph
    if let Some(bottom) = bottom {
        let max = pos.iter().map(|f| f.globalmismatch).max().unwrap();
        let mut chart = ChartBuilder::on(&bottom)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .right_y_label_area_size(60)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            /*.caption(
                format!(
                    "Break in coverage {} ({})",
                    loci.locus, loci.contig
                ),
                ("sans-serif", 40),
            )  */
            .build_cartesian_2d(
                loci.start.getobasedpos()..loci.end.getobasedpos(),
                0.0..max as f64 * 1.1 / 10_000.0,
            )
            .unwrap();
        let _ = chart
            .configure_mesh()
            .y_label_formatter(&|f| format!("{}%", (f * 10_000.0).round() / 100.0))
            .x_label_formatter(&|f| {
                format!(
                    "{} ({})",
                    f.to_formatted_string(&Locale::en),
                    pos.iter()
                        .find(|p| { p.position.getobasedpos() == *f })
                        .unwrap()
                        .locuspos
                        .getobasedpos()
                        .to_formatted_string(&Locale::en)
                )
            })
            .x_desc("Genomic position (bp)")
            .y_desc("Mismatch full rate (%)")
            .x_label_style(text_style.clone())
            .disable_x_mesh()
            .y_max_light_lines(2)
            .draw();
        chart
            .draw_series(AreaSeries::new(
                pos.iter().filter_map(|p| {
                    let score = p.globalmismatch as f64 / GLOBALMISMATCHFLOATING as f64;
                    if score.is_finite() && score != 0.0 && p.gettotalmap() > 0 {
                        Some((p.position.getobasedpos(), score))
                    } else {
                        None
                    }
                }),
                0.0,
                full_palette::DEEPORANGE_200.mix(0.8).filled(),
            ))
            .unwrap()
            .label("Mismatch full rate (%)")
            .legend(|(x, y)| {
                PathElement::new(vec![(x, y), (x + 15, y)], full_palette::DEEPORANGE_200)
            });
        if !args.nolegend {
            chart
                .configure_series_labels()
                .position(plotters::chart::SeriesLabelPosition::UpperRight)
                .label_font(text_style)
                .background_style(WHITE.mix(0.6))
                .border_style(BLACK.mix(0.8))
                .draw()
                .unwrap();
        }
    }
    // To avoid the IO failure being ignored silently, we manually call the present function
    drawnoticetext(&root);
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
}
fn readgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    args: &Args,
    root: DrawingArea<T, Shift>,
) -> Result<(), Box<dyn std::error::Error>>
where
    T: DrawingBackend,
{
    let text_style = fontstyle.into_text_style(&root);
    let max = pos.iter().map(|max| max.getmaxvalue()).max().unwrap() + 5;
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
        .build_cartesian_2d(loci.start.getobasedpos()..loci.end.getobasedpos(), 0..max)
        .unwrap();
    let _ = chart
        .configure_mesh()
        .x_label_formatter(&|f| {
            format!(
                "{} ({})",
                f.to_formatted_string(&Locale::en),
                pos.iter()
                    .find(|p| { p.position.getobasedpos() == *f })
                    .unwrap()
                    .locuspos
                    .getobasedpos()
                    .to_formatted_string(&Locale::en)
            )
        })
        //.x_desc("Genomic position (bp)")
        .y_desc("Coverage")
        .label_style(text_style.clone())
        .disable_x_mesh()
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw();
    chart
        .draw_series(AreaSeries::new(
            pos.iter().map(|p| (p.position.getobasedpos(), p.map0)),
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
            pos.iter().map(|p| (p.position.getobasedpos(), p.map1)),
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
            pos.iter().map(|p| (p.position.getobasedpos(), p.map60)),
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
            pos.iter().map(|p| (p.position.getobasedpos(), p.secondary)),
            full_palette::BLACK,
        ))
        .unwrap()
        .label("Secondary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLACK));
    chart
        .draw_series(LineSeries::new(
            pos.iter()
                .map(|p| (p.position.getobasedpos(), p.supplementary)),
            full_palette::BLUE_700,
        ))
        .unwrap()
        .label("Supplementary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLUE_700));
    chart
        .draw_series(LineSeries::new(
            pos.iter().map(|p| (p.position.getobasedpos(), p.overlaps)),
            full_palette::ORANGE_300,
        ))
        .unwrap()
        .label("Overlapping reads")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::ORANGE_300));
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .border_style(BLACK.mix(0.8))
            .label_font(text_style)
            .draw()
            .unwrap();
    }
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
        .build_cartesian_2d(
            (loci.start.getobasedpos()..loci.end.getobasedpos()).into_segmented(),
            0..1i64,
        )
        .unwrap();
    let text_style = fontstyle.into_text_style(&root);
    let _ = chart
        .configure_mesh()
        .x_desc("Genomic position (bp)")
        .x_label_style(text_style.clone())
        .disable_x_axis()
        .draw();
    let breaks = pos.iter().filter_map(|elem| {
        if elem.overlaps <= args.breaks.into() {
            Some((elem.position.getobasedpos(), 1))
        } else {
            None
        }
    });
    let finalpos = pos.iter().last().unwrap().position.getobasedpos();
    printbreaks(args, finalpos, breaks.clone(), loci, outputfile)?;
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
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .border_style(BLACK.mix(0.8))
            .label_font(text_style)
            .draw()
            .unwrap();
    }
    drawnoticetext(&root);
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    Ok(())
}

#![warn(clippy::unwrap_used)]
#![warn(clippy::expect_used)]
use bio::io::fasta;
use bio_types::sequence::SequenceRead;
/*
This software allows the analysis of BAM files to identify the confidence on a locus (specifically IG and TR) as well as allele confidence.
It was created and used by IMGT Team (https://www.imgt.org).
Available under EUPL license
Made by: Guilhem Zeitoun
*/
//TODO: Soft clips dans nouveau tableau et vérifier les valeurs.
use crate::identification::{downloadref, locusallposition};
///Assess quality of an assembly based on reads mapping, pourquoi la fin c'est 9 overlaps?, dû au samtools view
use clap::Parser;
use itertools::Itertools;
use plotters::coord::Shift;
use std::collections::HashMap;
use std::io::ErrorKind::InvalidInput;
use std::io::{stderr, stdout};
use std::num::NonZero;
use std::ops::RangeInclusive;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;
use std::{fs, io};
use strum::IntoEnumIterator;
//use noodles_fasta::{self as fasta, record::Sequence};
use crate::r#struct::*;
use crate::submissions::{
    REQUESTCLIENT, askforsubmission, checkifblastpresent, generatelightbam, genesblast,
    getspeciesfromncbi, positionfiltering,
};
use extended_htslib::bam::ext::{BamRecordExtensions, CsValue, IterAlignedPairs};
use extended_htslib::bam::{self, FetchDefinition, IndexedReader, Read};
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
mod identification;
mod r#struct;
mod submissions;
const VERSION: &str = env!("CARGO_PKG_VERSION");
const AUTHOR: &str = "IMGT®";
const GLOBALMISMATCHFLOATING: usize = 10_000;
const ALERTLOCUSSIZE: i64 = 10_000_000;
const MIN_READLENGTH: u64 = 1_000;
const READGAPMESSAGE: u64 = 200;
const MINIMUMCOVERAGE: usize = 10;
const MAXCOVERAGERATIO: f32 = 2.0;
const MATCHREADS: usize = 10;
const BORNES: usize = 10_000;
const SOFTCLIPRATIO: f32 = 0.4;
lazy_static! {
    static ref fontstyle: (&'static str, u32, &'static RGBColor) = {
        let args = Args::parse();
        ("sans-serif", args.fontlegendsize, &BLACK)
    };
    static ref NAME: String = env!("CARGO_PKG_NAME").replacen('_', "/", 1);
    #[allow(clippy::unwrap_used)]
    static ref regexpword: regex::Regex = regex::Regex::new(r"[^-\w()]").unwrap_or_else(|_| unreachable!("Regex issue"));
}
//Return block of positions thanks to CS/MD tag or CIGAR = (preferred if existing)
fn iterblock(record: &bam::Record) -> Option<Vec<[i64; 2]>> {
    match (record.getcsaligned(), record.aligned_blocks_match()) {
        //There is a CIGAR =
        (_, Some(d)) => Some(d.collect()),
        //There is a MD/CS tag
        (Some(d), None) => Some(
            d.into_iter()
                .filter_map(|p| match (&p.state, p.getgenomepos()) {
                    (CsValue::Same(d), Some(pos)) => Some([
                        pos,
                        pos.checked_add((*d).try_into().unwrap_or(0))
                            .and_then(|f| f.checked_sub(1))
                            .unwrap_or(0),
                    ]),
                    _ => None,
                })
                .collect(),
        ),
        //We have nothing
        (None, None) => None,
    }
}
/* fn blasttogenelist(list: &[Blastmatch], new: bool) -> Vec<GeneInfos> {
    let mut vec = Vec::new();
    for elem in list.iter().filter(|p| !new || p.onlynewalleles()) {
        let (posa, posb) = elem.getpos();
        let (posa, posb) = match (posa.try_into(), posb.try_into()) {
            (Ok(a), Ok(b)) => (a, b),
            _ => (0, 0),
        };
        let info = GeneInfos {
            gene: getnamefromblast(elem.qseqid()).unwrap_or("unknown".to_string()),
            chromosome: elem.sseqid.clone(),
            strand: elem.getstrand(),
            start: Position::new(false, min(posa, posb)),
            end: Position::new(false, max(posa, posb)),
            status: OkStatus::default(),
        };
        vec.push(info);
    }
    vec
} */
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
                    "{text}. Add --force to force even without = or MD/CS tag (some results won't be available)."
                );
                return (false, None, aligned);
            } else if !message {
                eprintln!("{text} but it was forced... Continuing...");
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
fn getassemblyreader(args: &Args) -> io::Result<fasta::IndexedReader<File>> {
    let (assembly, index) = match (args.assembly.as_ref(), args.assemblyindex.as_ref()) {
        (Some(a), Some(b)) => (Some(File::open(a)?), Some(File::open(b)?)),
        (Some(a), None) => (Some(File::open(a)?), None),
        (None, _) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No assembly provided",
            ));
        }
    };
    let elem = match (args.assembly.as_ref(), assembly, index) {
        (_, Some(a), Some(b)) => fasta::IndexedReader::new(a, b)
            .map_err(|p| io::Error::new(io::ErrorKind::InvalidInput, p.to_string())),
        (Some(p), ..) => fasta::IndexedReader::from_file(&p)
            .map_err(|p| io::Error::new(io::ErrorKind::InvalidInput, p.to_string())),
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No assembly provided",
            ));
        }
    };
    let elem = match elem {
        Ok(d) => Ok(d),
        Err(e) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Assembly error, maybe index is missing (create it with samtools): {e}"),
        )),
    }?;
    Ok(elem)
}
//Check we can read BAM file and return the reader with desired threads
fn getreaderoffile(args: &Args) -> Result<IndexedReader, extended_htslib::errors::Error> {
    let mut reader = match (&args.file, &args.index) {
        (Some(file), Some(d)) => bam::IndexedReader::from_path_and_index(file, d),
        (Some(file), _) => bam::IndexedReader::from_path(file),
        _ => {
            return Err(extended_htslib::errors::Error::FileSeek);
        }
    }?;
    //Set threads if given or thanks to parralelism. If nothing available, set to 4 by default.
    let threads = match (
        NonZero::new(args.threads),
        std::thread::available_parallelism(),
    ) {
        (Some(d), _) => d,
        #[allow(clippy::unwrap_used)]
        (_, Ok(d)) => min(d, NonZero::new(12).unwrap()),
        #[allow(clippy::unwrap_used)]
        _ => NonZero::new(4).unwrap(),
    };
    let _ = reader.set_threads(threads.get());
    Ok(reader)
}
//Check there is one alternate for one primary.
fn mergelocus(mut locus: Vec<LocusInfos>) -> Option<Vec<Vec<LocusInfos>>> {
    let mut elem: Vec<Vec<LocusInfos>> = Vec::with_capacity(locus.len());
    let mut actual: Vec<LocusInfos> = Vec::new();
    locus.sort_unstable();
    for loci in locus {
        match elem
            .iter()
            .find(|p| p.iter().any(|f| f.locus == loci.locus))
        {
            Some(d) if Some(d) != elem.last() => {
                eprintln!("Locus {} is splited! Aborted.", loci.locus);
                return None;
            }
            _ => (),
        };
        let alternate = if !loci.haplotype.isprimary() {
            if actual.is_empty() {
                eprintln!("Empty without a corresponding primary!");
                return None;
            } else {
                true
            }
        } else {
            false
        };
        match actual.first() {
            Some(e) if e.locus == loci.locus && alternate && actual.len() >= 2 => {
                eprintln!("Only one alternate is allowed!");
                return None;
            }
            Some(e) if e.locus == loci.locus && alternate => actual.push(loci.clone()),
            Some(e) if e.locus != loci.locus && alternate => {
                eprintln!("Alternate without a corresponding primary!");
                return None;
            }
            _ => {
                if !actual.is_empty() {
                    elem.push(actual.clone())
                };
                actual.clear();
                actual.push(loci.clone());
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
fn locusposparser(
    args: &Args,
    realspecies: &Species,
    blastpresent: bool,
) -> std::io::Result<(Vec<LocusInfos>, Option<Vec<Blastmatch>>)> {
    let realspecies = realspecies.to_string();
    let mut records = Vec::new();
    if let Some(arg) = &args.locuspos {
        let detectheaders = csv::ReaderBuilder::new()
            .flexible(true)
            .from_path(arg)
            .is_ok_and(|mut f| {
                f.headers()
                    .is_ok_and(|p| p.as_slice().to_ascii_lowercase().contains("locus"))
            });
        let mut csv = match csv::ReaderBuilder::new()
            .has_headers(detectheaders)
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .flexible(true)
            .from_path(arg)
        {
            Ok(c) => c,
            Err(e) => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("CSV file of location position cannot be found. Error is {e}"),
                ));
            }
        };
        for record in csv.deserialize::<FakeLocusinfo>() {
            let recorda = match record {
                Ok(r) => r,
                Err(e) => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Invalid CSV auto format. Cannot compute: {e}"),
                    ));
                }
            };
            records.push(recorda);
        }
    } else {
        let elem = itertools::iproduct!(Haplotype::iter(), Locus::iter());
        for (haplo, locu) in elem {
            records.push(FakeLocusinfo::new(
                locu,
                Some(haplo),
                None,
                None,
                None,
                None,
            ));
        }
        if args.haploid {
            records.retain(|p| p.haplotype.as_ref().map_or(false, |f| f.isprimary()));
        }
    }
    let mut locusrecord: (Vec<LocusInfos>, Option<Vec<Blastmatch>>) = (Vec::new(), None);
    match (
        args.assembly.as_ref(),
        records
            .iter_mut()
            .filter(|p| p.contig.is_none())
            .collect::<Vec<&mut FakeLocusinfo>>(),
    ) {
        (Some(path), b) if !b.is_empty() && blastpresent => {
            let (locus, blast) = locusallposition(path, realspecies, args)?;
            b.into_iter().for_each(|p| {
                let find = locus.iter().find(|r| {
                    r.haplotype == p.haplotype.clone().unwrap_or_default() && r.locus == p.locus
                });
                if let Some(a) = find {
                    *p = a.clone().into();
                }
            });
            locusrecord.1 = Some(blast);
        }
        (Some(_), b) if !b.is_empty() && !blastpresent => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Cannot perform automatic detection without BLAST. Please provide position."
                    .to_string(),
            ));
        }
        (None, b) if !b.is_empty() => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Automatic detection cannot be performed without an assembly. Provide position or provide assembly.".to_string(),
            ));
        }
        _ => (),
    }
    for elem in records.iter().filter(|p| p.contig.is_none()) {
        eprintln!(
            "Locus {} for {} could not be identified and would be skipped",
            &elem.locus,
            elem.haplotype.as_ref().unwrap_or(&Haplotype::default())
        );
    }
    records.retain(|p| p.contig.is_some());
    for record in records {
        let elem = record.intoloc()?;
        locusrecord.0.push(elem);
    }
    if locusrecord.0.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Empty CSV format for file {}, waiting locus\thaplotype (Primary or Alternate)\tcontig\tstart\tend. Nothing found.",
                args.locuspos
                    .as_ref()
                    .map(|f| format!("{}", f.display()))
                    .unwrap_or("no loc given".to_string())
            ),
        ));
    }
    //At least one duplicate line
    if let Some(d) = locusrecord.0.iter().duplicates().next() {
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
    locusrecord.0.iter_mut().for_each(|r| {
        if r.start >= r.end {
            (r.end, r.start) = (r.start, r.end);
            r.complement = Strand::Minus;
        }
    });
    if !args.hugeregion
        && let Some(big) = locusrecord
            .0
            .iter()
            .find(|p| p.end.length(&p.start) >= ALERTLOCUSSIZE)
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "The region {}-{} ({}) is more than {} bp and might be incorrect. Check your input file ({}) is correct.\nIf wanted, please add --hugeregion parameters. Be careful software might use a lot of memory.",
                big.start.getobasedpos(),
                big.end.getobasedpos(),
                big.locus,
                ALERTLOCUSSIZE.to_formatted_string(&Locale::en),
                &args
                    .locuspos
                    .as_ref()
                    .map(|f| format!("{}", f.display()))
                    .unwrap_or("no loc given".to_string())
            ),
        ));
    }
    Ok(locusrecord)
}
fn generategenelist<T>(
    locushashresult: Option<&HashMap<Locus, Vec<Blastmatch>>>,
    speciesblast: &Species,
    locus: &[LocusInfos],
    assembly: T,
    args: &mut Args,
) -> io::Result<Option<Vec<Blastmatch>>>
where
    T: AsRef<Path>,
{
    let speciesblast = speciesblast.to_string();
    let func = {
        eprintln!("You have not provided a gene list, BLASTING to get one.");
        locusallposition(assembly.as_ref(), &speciesblast, &args).map(|(_, b)| b)?
    };
    let mut data = match locushashresult {
        None => func,
        Some(a) if a.is_empty() => func,
        Some(b) => {
            eprintln!("Generating a gene list.");
            b.clone().into_values().flatten().collect_vec()
        }
    };
    //locusfiltering(&loci.locus, &mut data);
    positionfiltering(locus, &mut data);
    if data.is_empty() {
        eprintln!("No data after filtering gene list. Skipped.");
    } else {
        let mut finish: Vec<GeneInfos> = data
            .iter()
            .map(|p| {
                let strand = if p.send < p.sstart {
                    Strand::Minus
                } else {
                    Strand::Plus
                };
                GeneInfos::new(
                    p.getquery().gene.clone(),
                    p.getsubject().to_string(),
                    strand,
                    Position::new(false, p.sstart.try_into().unwrap_or_default()),
                    Position::new(false, p.send.try_into().unwrap_or_default()),
                )
            })
            .collect();
        checkandcorrectgenelistduplicate(&mut finish);
        let genenamefile = args.outdir.join("genelist_new.csv");
        let genefile = File::create(&genenamefile)?;
        let mut csv = csv::WriterBuilder::new()
            .delimiter(b',')
            .has_headers(true)
            .quote_style(csv::QuoteStyle::NonNumeric)
            .comment(Some(b'#'))
            .from_writer(genefile);
        for gene in finish.iter() {
            if let Err(e) = csv.serialize(gene) {
                return Err(io::Error::new(
                    InvalidInput,
                    format!("Cannot print gene list. Error is {e}"),
                ));
            }
        }
        csv.flush()?;
        args.geneloc = Some(genenamefile);
        return Ok(Some(data));
    }
    Ok(None)
}
fn printnewloc(args: &Args, locus: &[LocusInfos]) -> io::Result<()> {
    let file = File::create(args.outdir.join("newloc.csv"))?;
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(file);
    csv.write_record(["locus", "haplotype", "contig", "start", "end", "status"])?;
    for loci in locus.iter() {
        csv.serialize(loci)?;
    }
    csv.flush()?;
    Ok(())
}
//Check BAM file exists and outputdir is created and return it
fn checklocusandoutput(args: &Args) -> std::io::Result<&PathBuf> {
    //Check bam file exists
    if let Some(a) = &args.file
        && let Err(e) = getreaderoffile(args)
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!(
                "Cannot read bam file ({}). Error is: {}. Exiting.",
                a.display(),
                e
            ),
        ));
    }
    //Check assembly is okay if existing
    if args.assembly.is_some() {
        getassemblyreader(args)?;
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
    mut message: bool,
    locus: &LocusInfos,
    record: &bam::Record,
    sep: i64,
) -> io::Result<bool> {
    let start = &Position::new(true, record.reference_start());
    let end = &Position::new(true, record.reference_end());
    //Get range to put the reads inclusive pos
    let newrange = Position::new(
        true,
        max(record.reference_start(), locus.start.getzbasedpos()),
    )
        ..Position::new(true, min(locus.end.getobasedpos(), record.reference_end()));
    if pos.contains_key(start)
        && let Some(d) = pos.get_mut(start)
        && record.cigar().leading_softclips() > 0
    {
        d.softclips += 1.0;
    } else if pos.contains_key(end)
        && let Some(d) = pos.get_mut(end)
        && record.cigar().trailing_softclips() > 0
    {
        d.softclips += 1.0;
    }
    let (message, matched, aligned) = match iteralert(args, message, record) {
        (_, None, _) => {
            return Ok(false);
        } //Kill software, errors sent by iteralert
        (newmessage, Some(p), aligned) => {
            message = newmessage;
            (message, p, aligned)
        }
    };
    for (i, targeting) in pos.range_mut(newrange) {
        if record.is_secondary() || record.is_supplementary() {
            if record.is_secondary() {
                targeting.secondary += 1;
            } else {
                targeting.supplementary += 1;
                /* targeting
                .supplementary
                .push(String::from_utf8_lossy(p.qname()).to_string()); */
            }
            return Ok(false);
        }
        let i = &i.getzbasedpos();
        let _time = Instant::now();
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
                return Ok(false);
            }
        };
        targeting.globalmismatch += getglobalmismatch(args, record);
        /* match record.mapq() {
            0 => (),
            _ => targeting.globalmismatch += getglobalmismatch(args, record),
        }; */
        if i % sep == 0 {
            //Get quality of reads depending on genomic position
            if let Some(d) = record.aligned_pairs().find(|p| p[1] == *i) {
                let index = d[0] as usize;
                targeting.qual += record.qual().get(index).map_or(0, |f| (*f).into());
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
    Ok(message)
}
fn getmeancoverageandlength(args: &Args) -> std::io::Result<(u64, u64)> {
    let time = Instant::now();
    let mut reader: IndexedReader = getreaderoffile(args)
        .map_err(|f| std::io::Error::new(std::io::ErrorKind::InvalidInput, f))?;
    reader
        .fetch(FetchDefinition::All)
        .map_err(|f| std::io::Error::new(std::io::ErrorKind::InvalidInput, f))?;
    let mut readlength = 0;
    let mut count = 0;
    let length = reader
        .rc_records()
        .filter_map(Result::ok)
        .filter_map(|f| {
            //Remove reads with 0 as MAPQ
            if f.mapq() != 0 && !f.is_unmapped() && !f.is_secondary() && !f.is_supplementary() {
                readlength += f.len();
                count += 1;
                Some(usize::try_from(f.reference_end() - f.reference_start() + 1).unwrap_or(0))
            } else {
                None
            }
        })
        .sum::<usize>();
    let total_mapped = match args.extractedlength {
        0 => {
            println!("Going for stats");
            let stats = reader
                .index_stats()
                .map_err(|f| std::io::Error::new(std::io::ErrorKind::InvalidInput, f))?;
            stats.iter().map(|f| f.1).sum::<u64>()
        }
        e => {
            println!("Forced length used");
            e
        }
    };
    let values = u64::try_from(length)
        .map_err(|f| std::io::Error::new(std::io::ErrorKind::InvalidInput, f))?;
    let readvalues = u64::try_from(readlength)
        .map_err(|f| std::io::Error::new(std::io::ErrorKind::InvalidInput, f))?;
    let mean = values / total_mapped;
    let readavg = readvalues / count;
    println!("Calculated in {} s", time.elapsed().as_secs());
    Ok((mean, readavg))
}
fn calculatemean(meanpath: &Path, args: &Args) -> io::Result<u64> {
    match (
        args.meancoverage,
        meanpath.try_exists(),
        fs::read_to_string(meanpath).map(|p| p.parse::<u64>()),
    ) {
        //Mean is given
        (Some(mean), ..) => {
            println!(
                "Mean was provided, getting given value and setting the value for further usage."
            );
            let _ = fs::write(meanpath, mean.to_string().as_bytes());
            Ok(mean)
        }
        //The mean was calculated with a file
        (None, Ok(true), Ok(Ok(mean))) => {
            println!("Mean was already calculated, retrieving...");
            Ok(mean)
        }
        _ => {
            println!("Getting mean coverage from calculations, it might take some minutes.");
            let (mean, average) = match getmeancoverageandlength(args) {
                Ok((0, _)) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Probably truncated BAM, you should give extracted length with the argument --extractedlength.",
                    ));
                }
                Ok((_, b)) if b < MIN_READLENGTH => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "The average read length does not exceed {MIN_READLENGTH} bp. For accurate results, please provide long-reads to the software.",
                    ));
                }
                Ok((a, b)) => (a, b),
                Err(e) => {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, format! {"{e}"}));
                }
            };
            println!(
                "Mean coverage is {} and average length of reads is {}",
                mean, average
            );
            let _ = fs::write(meanpath, mean.to_string().as_bytes());
            Ok(mean)
        }
    }
}
#[allow(clippy::complexity)]
fn paintgraph<T>(
    outputfileread: &std::path::Path,
    outputfilemismatch: &std::path::Path,
    loci: &LocusInfos,
    pos: &BTreeMap<Position, HashMapinfo>,
    args: &Args,
    readgraphelem: DrawingArea<T, Shift>,
    mismatchgraphelem: DrawingArea<T, Shift>,
    mean: u64,
) -> io::Result<()>
where
    T: DrawingBackend,
{
    if let Err(e) = readgraph(
        outputfileread,
        loci,
        pos.values().collect_vec().as_slice(),
        args,
        readgraphelem,
        mean,
    )
    .and_then(|_| {
        mismatchgraph(
            outputfilemismatch,
            loci,
            pos.values().collect_vec().as_slice(),
            args,
            mismatchgraphelem,
        )
    }) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{}", e),
        ));
    };
    Ok(())
}
fn main() -> ExitCode {
    let mut args = Args::parse();
    let firstinstant = Instant::now();
    let speciesblast = match Species::new(&args.species) {
        Ok(b) => b,
        Err(_) => return ExitCode::FAILURE,
    };
    println!(
        "Species is {} (taxon: {}).",
        speciesblast.to_string(),
        speciesblast
            .getid()
            .map_or("Unknown".to_string(), |f| f.to_string())
    );
    if args.percentalerting >= args.percentwarning {
        eprintln!("Percent warning must be greater or equal than percent alerting.");
        return ExitCode::FAILURE;
    }
    let blastpresent = checkifblastpresent();
    if !blastpresent {
        eprintln!("BLAST is not found and won't be used, some analysis won't be performed");
        args.nosubmit = true;
    }
    //Get locus, geneloc and outputdir, print errors if we have
    let outputdir = match checklocusandoutput(&args) {
        Ok(a) => a.clone(),
        Err(e) => {
            eprintln!("Error with locus file or output: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (locus, blastcheck) = match locusposparser(&args, &speciesblast, blastpresent) {
        Err(f) => {
            eprintln!("Error locus parser: {f}");
            return ExitCode::FAILURE;
        }
        Ok(b) => b,
    };
    if args.geneloc.is_some()
        && let Err(e) = checkgenelistformat(&args)
    {
        eprintln!("Error gene list: {e}");
        return ExitCode::FAILURE;
    }
    if args.command == Command::Find {
        if !blastpresent {
            eprintln!("To find the locus, BLAST must be present");
            return ExitCode::FAILURE;
        }
        if let Err(e) = printnewloc(&args, &locus) {
            eprintln!("Error setting new locus result: {e}");
            return ExitCode::FAILURE;
        };
        let r = args.assembly.clone();
        let infos: HashMap<Locus, Vec<Blastmatch>> = if let Some(a) = blastcheck {
            a.into_iter()
                .filter_map(|f| {
                    if let Some(a) = locus.iter().find(|loci| {
                        f.sseqid == loci.contig
                            && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                                .contains(&f.sstart.try_into().unwrap_or_default())
                    }) {
                        Some((a.clone(), f))
                    } else {
                        None
                    }
                })
                .into_group_map_by(|c| c.0)
                .values_mut()
                .map(|p| p.into_iter().map(|a| a.1).collect_vec())
        } else {
            HashMap::new()
        };
        if infos.is_empty() {
            eprintln!("Error setting gene list result. Empty list.");
            return ExitCode::FAILURE;
        }
        if let Some(b) = r
            && let Err(e) = generategenelist(Some(infos), &speciesblast, &locus, b, &mut args)
        {
            eprintln!("Error setting gene list result: {e}");
            return ExitCode::FAILURE;
        }
        endmessage(firstinstant, &args);
        return ExitCode::SUCCESS;
    }
    let initiallocus = &locus;
    //Group between primary and alternate
    let mut grouped = match mergelocus(locus.clone()) {
        Some(g) => g,
        None => {
            eprintln!(
                "Check order of loci in the file {}.",
                args.locuspos
                    .map(|f| format!("{}", f.display()))
                    .unwrap_or("no loc given".to_string())
            );
            return ExitCode::FAILURE;
        }
    };
    let meanpath = outputdir.join(".mean");
    let mean = match calculatemean(&meanpath, &args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::FAILURE;
        }
    };
    let mut locushashresult: HashMap<Locus, Vec<Blastmatch>> = HashMap::new();
    for locus in grouped.iter_mut() {
        let haplotype = locus.len();
        let nlocus = locus.clone();
        let floci = match nlocus.first() {
            Some(f) => f,
            None => {
                eprintln!("There is no locus after grouping found.");
                return ExitCode::FAILURE;
            }
        };
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
            1000 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap_or(1),
        );
        let mismatchsize = (
            1200,
            600 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap_or(1),
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
                    &format!("{fgraph}.svg"),
                    true,
                ));
                outputfile1 = outputfile;
                let root = SVGBackend::new(&outputfile1, readresultsize).into_drawing_area();
                if haplotypebool {
                    let (top, bottom) = root.split_vertically((50).percent_height());
                    (Some(top), Some(bottom), None, None)
                    //readgraph(outpufile, locus.first().unwrap(), &pos, &args, root);
                } else {
                    (Some(root), None, None, None)
                }
            } else {
                let outputfile = outputdir.join(givename(
                    &args.species,
                    &floci.locus,
                    &floci.contig,
                    haplotypebool,
                    &format!("{fgraph}.png"),
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
                    &format!("{sgraph}.svg"),
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
                    &format!("{sgraph}.png"),
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
        for loci in locus.iter_mut() {
            let mut reader = match getreaderoffile(&args) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("Cannot read bam file. Error is {e}. Exiting");
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
                    "The region {}:{}-{} cannot be found, skipped.",
                    loci.contig,
                    loci.start.getobasedpos(),
                    loci.end.getobasedpos()
                );
                continue;
                //return ExitCode::FAILURE;
            };
            let mut nocount = true;
            //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
            //let file = File::create(&filename).unwrap();
            //let mut writer = BufWriter::new(file);
            let locusrange = loci.start.getzbasedpos()..=loci.end.getzbasedpos();
            let mut pos: BTreeMap<Position, HashMapinfo> = BTreeMap::new();
            //Populate all B-Tree position, 0-based, invert if locus complement
            let iterator: Vec<(usize, i64)> = if loci.complement.isrev() {
                locusrange.rev().enumerate().collect()
            } else {
                locusrange.enumerate().collect()
            };
            let iterator: Vec<(i64, i64)> = iterator
                .into_iter()
                .filter_map(|(a, b)| match i64::try_from(a) {
                    Ok(d) => Some((d, b)),
                    Err(_) => {
                        eprintln!("Position {a} is too big and is not currently supported");
                        None
                    }
                })
                .collect();
            iterator.into_iter().for_each(|(i, p)| {
                let default = HashMapinfo {
                    locuspos: Position::new(true, i),
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
                //Print every x reads done
                if count % READGAPMESSAGE == 0 {
                    let _ = writeln!(
                        lock,
                        "Processed {} reads in {:.3} s",
                        count.to_formatted_string(&Locale::en),
                        Instant::now().saturating_duration_since(time).as_secs_f32()
                    );
                }
                nocount = false;
                match processcounting(&args, &mut pos, message, loci, &p, sep) {
                    Err(e) => {
                        eprintln!("Error is {e}");
                        return ExitCode::FAILURE;
                    }
                    Ok(b) => {
                        if !message && b {
                            message = b;
                        }
                    }
                }
            }
            loci.setstatus(mean, &pos);
            if nocount {
                eprintln!(
                    "The region {}:{}-{} has no data, skipped.",
                    loci.contig,
                    loci.start.getobasedpos(),
                    loci.end.getobasedpos()
                );
                continue;
                //return ExitCode::FAILURE;
            }
            //Quality and mismatch is the sum of reads so dividing to get real results
            pos.iter_mut().for_each(|(_, p)| {
                if p.qual > 0 {
                    p.qual /= max(
                        std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap_or(1),
                        1,
                    )
                }
                if p.softclips.is_normal() {
                    p.softclips =
                        (p.softclips * 100f32 / max(p.gettotalmap(), 1) as f32).round() / 100f32
                }
                p.globalmismatch /= max(
                    std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap_or(1),
                    1,
                );
            });
            println!("Making graphs");
            match (
                loci.haplotype.isprimary(),
                args.svg,
                readgraphtop.clone(),
                mismatchgraphtop.clone(),
                readgraphbottom.clone(),
                mismatchgraphbottom.clone(),
                readgraphtop2.clone(),
                mismatchgraphtop2.clone(),
                readgraphbottom2.clone(),
                mismatchgraphbottom2.clone(),
            ) {
                (true, true, Some(readgraphtop), Some(mismatchgraphtop), ..)
                    if !outputfile1.file_name().is_none_or(|p| p.is_empty())
                        && !outputfile3.file_name().is_none_or(|p| p.is_empty()) =>
                {
                    if let Err(e) = paintgraph(
                        &outputfile1,
                        &outputfile3,
                        loci,
                        &pos,
                        &args,
                        readgraphtop,
                        mismatchgraphtop,
                        mean,
                    ) {
                        eprintln!("Cannot create read graphs. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                }
                (false, true, _, _, Some(readgraphbottom), Some(mismatchgraphbottom), ..)
                    if !outputfile1.file_name().is_none_or(|p| p.is_empty())
                        && !outputfile3.file_name().is_none_or(|p| p.is_empty()) =>
                {
                    if let Err(e) = paintgraph(
                        &outputfile1,
                        &outputfile3,
                        loci,
                        &pos,
                        &args,
                        readgraphbottom,
                        mismatchgraphbottom,
                        mean,
                    ) {
                        eprintln!("Cannot create read graphs. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                }
                (true, false, _, _, _, _, Some(readgraphtop2), Some(mismatchgraphtop2), ..)
                    if !outputfile2.file_name().is_none_or(|p| p.is_empty())
                        && !outputfile4.file_name().is_none_or(|p| p.is_empty()) =>
                {
                    if let Err(e) = paintgraph(
                        &outputfile2,
                        &outputfile4,
                        loci,
                        &pos,
                        &args,
                        readgraphtop2,
                        mismatchgraphtop2,
                        mean,
                    ) {
                        eprintln!("Cannot create read graphs. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                }
                (false, false, .., Some(readgraphbottom2), Some(mismatchgraphbottom2))
                    if !outputfile2.file_name().is_none_or(|p| p.is_empty())
                        && !outputfile4.file_name().is_none_or(|p| p.is_empty()) =>
                {
                    if let Err(e) = paintgraph(
                        &outputfile2,
                        &outputfile4,
                        loci,
                        &pos,
                        &args,
                        readgraphbottom2,
                        mismatchgraphbottom2,
                        mean,
                    ) {
                        eprintln!("Cannot create read graphs. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                }
                _ => {
                    eprintln!(
                        "Graphs for {} ({}) could not be created, skipped.",
                        loci.locus, loci.haplotype
                    );
                    continue;
                }
            }
            println!("Graphs finished.");
            //Create CSV from HashMap
            if let Err(e) = createcsv(
                outputdir.as_path(),
                loci,
                pos.values().collect_vec().as_slice(),
                &args,
            ) {
                eprintln!("Cannot create csv file. Error is {e}");
                return ExitCode::FAILURE;
            }
            println!(
                "Locus {} ({}) is {}.",
                loci.locus, loci.haplotype, loci.status
            );
            //Create gene CSV
            if let Some(blast) = &blastcheck {
                let element: Vec<Blastmatch> = blast
                    .iter()
                    .filter_map(|f| {
                        if f.sseqid == loci.contig
                            && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                                .contains(&f.sstart.try_into().unwrap_or_default())
                        {
                            Some(f.clone())
                        } else {
                            None
                        }
                    })
                    .collect();
                locushashresult.insert(loci.locus.clone(), element);
            }
            if args.geneloc.is_some() {
                println!("Gene list starting!");
                match genelist(loci, &args, false) {
                    Err(e) => {
                        eprintln!("Cannot create gene list. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                    Ok(b) if b.is_empty() => {
                        eprintln!("No gene list for locus {}. Skipped.", loci.locus);
                    }
                    Ok(b) => {
                        if args.assembly.as_ref().is_some() && !args.nosubmit && blastpresent {
                            println!("Blasting gene list");
                            let geneinfo: Vec<GeneInfos> =
                                b.iter().map(|f| f.clone().into()).collect();
                            match genesblast(&geneinfo, &args, &speciesblast, &loci.locus) {
                                Ok(blast) => {
                                    if loci.status.isvalid() && b.iter().any(|f| f.status.isvalid())
                                    {
                                        locushashresult.insert(loci.locus.clone(), blast);
                                    }
                                }
                                Err(e) => {
                                    eprintln!("Cannot blast gene list. Error is {e}");
                                    continue;
                                }
                            };
                        }
                    }
                }
                println!("Gene list finished.");
            } else {
                let r = args.assembly.clone();
                if !args.nosubmit
                    && let Some(b) = r
                {
                    let mut data = match generategenelist(
                        Some(&locushashresult),
                        &speciesblast,
                        &initiallocus,
                        b,
                        &mut args,
                    ) {
                        Ok(Some(a)) => a,
                        Ok(None) => continue,
                        Err(e) => {
                            eprintln!(
                                "Gene list generation for locus {} has an error: {}.",
                                loci.locus, e
                            );
                            continue;
                        }
                    };
                    let locivec = vec![loci.clone()];
                    positionfiltering(&locivec, &mut data);
                    let result = match genelist(loci, &args, false) {
                        Err(e) => {
                            eprintln!(
                                "Gene list generation for locus {} has an error: {}.",
                                loci.locus, e
                            );
                            continue;
                        }
                        Ok(b) => b,
                    };
                    if loci.status.isvalid() && result.iter().any(|f| f.status.isvalid()) {
                        locushashresult.insert(loci.locus.clone(), data);
                    }
                } else if !args.nosubmit {
                    eprintln!("No assembly to check gene list.");
                }
            }
        }
        println!("Locus {} is done!", &floci.locus);
    }
    let mergedloci: Vec<LocusInfos> = grouped.iter().flatten().cloned().collect();
    if let Some(light) = &args.outlightbam
        && let Err(e) = generatelightbam(&args, light, &mergedloci)
    {
        eprintln!("{e}");
        return ExitCode::FAILURE;
    }
    if let Err(e) = printnewloc(&args, &mergedloci) {
        eprintln!("Error setting new locus result: {e}");
    }
    if let Err(e) = printpotentialalleles(
        &args,
        downloadref(false).map(|(_, release)| release),
        &locushashresult,
    ) {
        eprintln!("Error setting alleles result: {e}");
    }
    if !args.nosubmit
    /* && locushashresult
    .iter()
    .any(|(_, f)| f.iter().any(|g| g.onlynewalleles())) */
    {
        match askforsubmission(&speciesblast, &mergedloci, &args, &locushashresult) {
            Ok(_) => (),
            Err(e) => eprintln!("Error submitting sequences: {e}"),
        }
    }
    endmessage(firstinstant, &args);
    ExitCode::SUCCESS
}
fn endmessage(firstinstant: Instant, args: &Args) {
    println!(
        "{} (version {}) done sucessfully in {:.3} seconds. Output files are located in {}.",
        NAME.as_str(),
        VERSION,
        firstinstant.elapsed().as_secs_f32(),
        args.outdir.display()
    );
}
fn printpotentialbornes<T>(bornes: &[T], args: &Args) -> io::Result<()>
where
    T: Blastcalc,
{
    let path = args.outdir.join("potentialbornes.csv");
    let writer = File::create(path)?;
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b',')
        .comment(Some(b'#'))
        .flexible(true)
        .quote_style(csv::QuoteStyle::Never)
        .from_writer(writer);
    for borne in bornes {
        csv.write_record(&[
            borne.getallelename(),
            borne.getsubject(),
            &format!("{}", borne.getpos().0),
            &format!("{}", borne.getpos().1),
            borne.getstrand().to_string().as_str(),
        ])?;
    }
    csv.flush()?;
    Ok(())
}
fn printpotentialalleles<T>(
    args: &Args,
    release: Option<T>,
    locushash: &HashMap<Locus, Vec<Blastmatch>>,
) -> io::Result<()>
where
    T: AsRef<str>,
{
    let path = args.outdir.join("potentialalleles.fasta");
    let writer = File::create(path)?;
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b',')
        .comment(Some(b'#'))
        .flexible(true)
        .quote_style(csv::QuoteStyle::Never)
        .from_writer(writer);
    csv.write_record([
        "name",
        "position",
        "locus",
        "strand",
        "bestmatch",
        "identity",
    ])?;
    csv.write_record([b"#Name of the gene given or detected, Locus, Strand, Best match based on IMGT-GENE-DB/Identity (in %)"])?;
    csv.write_record(&[format!(
        "#IMGT/GENE-DB release {}",
        release.map_or("Unknown".to_string(), |p| String::from(p.as_ref()))
    )])?;
    for (locus, elem) in locushash {
        for matches in elem {
            let name = matches.getallelename();
            let (start, end) = matches.getpos();
            csv.write_record(&[
                name.to_string(),
                format!("{}-{}", start, end),
                locus.to_string(),
                matches.complement.to_string(),
                name.to_string(),
                matches.getidentity().to_string(),
            ])?;
        }
    }
    Ok(())
}
fn checkgenelistformat(args: &Args) -> Result<Vec<GeneInfos>, Box<dyn std::error::Error>> {
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
        .from_path(geneloc)
        .map_err(|_| {
            Box::new(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("The file given {} is not found", geneloc.display()),
            ))
        })?;
    for record in csv.deserialize() {
        let record = match record {
            Ok(r) => r,
            Err(e) => {
                println!("Line {e} is passed.");
                continue;
                /* return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Invalid CSV format (comma separated) for file {}, waiting gene,chromosome,strand,start,end case sensitive. Have you kept the header?\n{}",
                        geneloc.display(),
                        e
                    ),
                ))); */
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
    checkandcorrectgenelistduplicate(&mut genes);
    Ok(genes)
}
/// Add underscore in case gene list has duplicates
fn checkandcorrectgenelistduplicate(genes: &mut [GeneInfos]) {
    //Invert if start >= end because strand is given and won't work
    genes.iter_mut().filter(|p| p.start > p.end).for_each(|p| {
        (p.start, p.end) = (p.end, p.start);
    });
    let geneclone = genes.to_vec();
    let finish: Vec<&GeneInfos> = geneclone
        .iter()
        .duplicates_by(|g| (g.gene.as_str(), g.chromosome.as_str()))
        .collect();
    if !finish.is_empty() {
        println!(
            "Some genes has the same name on the same chromosome. Underscore (_) would be added"
        );
        let mut count = 0;
        for name in genes.iter_mut() {
            name.gene = if finish.iter().any(|g| g.gene.as_str() == name.gene.as_str()) {
                count += 1;
                format!("{}_{}", name.gene, count)
                    .replace(",", "_")
                    .trim()
                    .to_string()
            } else {
                name.gene.to_string().replace(",", "_").trim().to_string()
            };
            name.chromosome = name.chromosome.replace(",", "_").trim().to_string()
        }
    }
    genes.sort_unstable();
}
fn extractgenelist(
    args: &Args,
    loci: &LocusInfos,
    full: bool,
) -> Result<Vec<GeneInfos>, Box<dyn std::error::Error>> {
    let mut genes = checkgenelistformat(args)?;
    //Retain genes inside the correct loci
    if !full {
        genes.retain(|gene| {
            gene.chromosome == loci.contig
                && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                    .contains(&gene.start.getobasedpos())
                && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                    .contains(&gene.end.getobasedpos())
        });
    }
    if genes.is_empty() {
        println!("No gene identified for locus {}, skipped.", loci.locus);
        return Ok(Vec::new());
    }
    //At least one duplicate line
    checkandcorrectgenelistduplicate(&mut genes);
    Ok(genes)
}
fn genelist(
    loci: &LocusInfos,
    args: &Args,
    full: bool,
) -> Result<Vec<GeneInfosFinish>, Box<dyn std::error::Error>> {
    let genes = extractgenelist(args, loci, full)?;
    if genes.is_empty() {
        return Ok(Vec::new());
    }
    let outputdir = &args.outdir;
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
    for mut gene in genes {
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
            let mut hash: BTreeMap<Position, Posread> = BTreeMap::new(); //Match and full match and total
            //Hash contains 1-based positions
            genegenericrange.for_each(|p| {
                hash.insert(
                    Position::new(true, p),
                    //Default should not trigger as no error possible
                    Posread::new(0, 0, 0, 0f32, args)
                        .unwrap_or_else(|_| unreachable!("Error on Posread")),
                );
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
            if let Some(d) = hash.get_mut(&Position::new(true, record.reference_start()))
                && record.cigar().leading_softclips() > 0
            {
                d.softclips += 1f32
            } else if let Some(d) = hash.get_mut(&Position::new(true, record.reference_end()))
                && record.cigar().trailing_softclips() > 0
            {
                d.softclips += 1f32
            }
            let range = record.reference_start()..record.reference_end();
            coverageperc += ranges::Ranges::from(range.clone()).into_iter().count();
            'outer: for [start, end] in record.aligned_blocks() {
                for p in start..end {
                    match hash.get_mut(&Position::new(true, p)) {
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
                'outer: for [start, end] in iterblock(&record).unwrap_or_default() {
                    for p in start..end {
                        match hash.get_mut(&Position::new(true, p)) {
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
                    .is_some_and(|f| f.into_iter().any(|p| validrange(p, &gene.start, &gene.end)))
            {
                reads100m += 1;
            }
            for p in range {
                match hash.get_mut(&Position::new(true, p)) {
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
            let _ = writeln!(lock, "Empty records for gene {}", gene.gene);
            //PUT 0 value on the CSV
            let elem = GeneInfosFinish::make_default(gene);
            finale.push(elem);
            continue;
        }
        hash.iter_mut().for_each(|(_, p)| {
            if p.softclips.is_normal() {
                p.softclips = (p.softclips * 100f32 / max(p.gettotal(), 1) as f32).round() / 100f32
            }
        });
        //Coverage calculus
        let coverage = hash
            .iter()
            .filter(|(_, p)| p.gettotal() >= args.coverage.try_into().unwrap_or(usize::MAX))
            .count();
        //Reverse if complement
        let text = {
            let iterator: Vec<(&Position, &Posread)> = match gene.strand {
                Strand::Plus => hash.iter().collect(),
                Strand::Minus => hash.iter().rev().collect(),
            };
            //Merging data
            iterator.into_iter().fold(String::new(), |mut acc, (_, f)| {
                acc.push_str(&format!(
                    "{}({}={}X{}ID)-",
                    f.gettotal(),
                    f.getmatch(),
                    f.getmismatchcount(),
                    f.getindelcount()
                ));
                acc
            })
        };
        let text = String::from(text.trim_end_matches('-'));
        let plots = outputdir.join(format!(
            "gene_{}",
            loci.haplotype.to_string().as_str().to_lowercase()
        ));
        if !std::fs::exists(&plots)? {
            println!("Creating the folder {}", plots.display());
            std::fs::create_dir_all(&plots)?;
        };
        let mut output = plots.join(regexpword.replace_all(&gene.gene, "_").to_uppercase());
        gene.setstatus(reads100m, &hash);
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
            )?;
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
            )?;
        }
        let coverageperc = ((coverageperc * 1_000
            / reads
            / usize::try_from(gene.end.length(&gene.start)).unwrap_or(usize::MAX))
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
    for gene in finale.iter() {
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
    Ok(finale)
}
fn printbreaks(
    args: &Args,
    finalpos: i64,
    breaks: &[(i64, i64)],
    loci: &LocusInfos,
    outputfile: &std::path::Path,
) -> std::io::Result<()> {
    let outputdir = match outputfile.parent() {
        Some(d) => d,
        None => {
            return Err(io::Error::new(
                io::ErrorKind::NotADirectory,
                "Invalid root for directory, cannot create break file.",
            ));
        }
    };
    let mut breakfile = File::create(outputdir.join(givename(
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
    let finalbreak = match breaks.iter().map(|p| p.0).max() {
        Some(d) => d,
        None => {
            breakfile.write_all("No breaks.".as_bytes())?;
            return Ok(());
        }
    };
    let mut acc = breaks.iter().fold(String::new(), |mut acc, (num, _)| {
        if first.is_none() {
            first = Some(num);
            prev = Some(num);
        } else if let (Some(mut prev_num), Some(f)) = (prev, first) {
            if num - prev_num != 1 || *num == finalpos || *num == finalbreak {
                if *num == finalpos {
                    prev_num = &finalpos;
                } else if *num == finalbreak {
                    prev_num = &finalbreak;
                }
                if f == prev_num {
                    acc.push_str(&format!("{}:{}\n", loci.contig, f));
                } else {
                    acc.push_str(&format!("{}:{}..{}\n", loci.contig, f, prev_num));
                }
                if *num != finalpos && *num != finalbreak {
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
    hash: &BTreeMap<Position, Posread>,
    gene: &GeneInfos,
    loci: &LocusInfos,
    root: DrawingArea<T, Shift>,
    alerting: &mut BTreeMap<GeneInfos, Vec<(bool, usize)>>,
    reads100m: usize,
) -> Result<(), Box<dyn std::error::Error>>
where
    T: DrawingBackend,
{
    let genename = gene.gene.to_string();
    let text_style = fontstyle.into_text_style(&root);
    let _ = root.fill(&plotters::prelude::WHITE);
    let max = hash.values().map(|p| p.gettotal()).max().unwrap_or(0) + 5;
    let colorgene = if gene.status.isvalid() {
        full_palette::GREEN
    } else {
        full_palette::RED
    };

    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            format!(
                "Reads coverage for {} ({}-{})",
                genename, loci.haplotype, gene.strand
            ),
            ("sans-serif", 22, &colorgene),
        )
        .build_cartesian_2d(1..hash.len(), 0..max)
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    //Reverse complement genes
    let hash: Vec<(&Position, &Posread)> = match gene.strand {
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
            .label("Suspicious positions")
            .legend(|(x, y)| {
                plotters::element::Rectangle::new(
                    [(x, y), (x + 15, y + 5)],
                    full_palette::RED_400.filled(),
                )
            });
    }
    //Continue graph
    /* chart
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
    } */
    //Secondary
    /* let softclipmax = ((softclipmax as f32 /max as f32 * 10.0).ceil() / 10.0 * max as f32) as i64;
    //let softclipmax = core::cmp::max(calc,max);
     */
    let mut secondary = chart.set_secondary_coord(
        usize::try_from(loci.start.getobasedpos()).unwrap_or(0)
            ..usize::try_from(loci.end.getobasedpos()).unwrap_or(0),
        0..100usize,
    );
    secondary
        .configure_mesh()
        .x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        .x_desc("Genomic position (bp)")
        .y_desc("Reads count")
        .disable_x_mesh()
        .label_style(text_style.clone())
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    //let mut second = chart.set_secondary_coord(loci.start..loci.end, 0..max);
    secondary
        .draw_secondary_series(
            Histogram::vertical(&secondary)
                .baseline(0)
                .margin(3)
                .data(hash.iter().filter_map(|(pos, p)| {
                    if p.softclips.is_normal() {
                        Some((
                            usize::try_from(pos.getobasedpos()).unwrap_or(0),
                            (p.softclips * 100f32) as usize,
                        ))
                    } else {
                        None
                    }
                }))
                .style(full_palette::BLACK.mix(0.4).filled()),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Average softclips")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::BLACK.filled(),
            )
        });
    secondary
        .configure_secondary_axes()
        .y_desc("Average softclips")
        .y_label_formatter(&|f| format!("{f}%"))
        .label_style(text_style.clone())
        //.disable_y_mesh()
        .draw()
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    if !args.nolegend {
        secondary
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    }
    // To avoid the IO failure being ignored silently, we manually call the present function
    drawnoticetext(&root)?;
    root.present().map_err(|_e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir")))?;
    Ok(())
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
        locus.safestring(),
        if image && !haplo {
            String::new()
        } else {
            format!("{contig}-")
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
fn drawnoticetext<T>(root: &DrawingArea<T, Shift>) -> Result<(), Box<dyn std::error::Error>>
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
    let size = root
        .estimate_text_size(&text, &text_style)
        .unwrap_or_default();
    let size = (
        (root.dim_in_pixel().0 - size.0 - 10) as i32,
        (root.dim_in_pixel().1 - size.1 - 10) as i32,
    );
    Ok(root
        .draw_text(&text, &text_style, size)
        .map_err(|p| Box::new(io::Error::other(p.to_string())))?)
}
fn mismatchgraph<T>(
    _outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    args: &Args,
    root: DrawingArea<T, Shift>,
) -> Result<(), Box<dyn std::error::Error>>
where
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Misalign (%)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::RED_400));
    let mut secondary = chart.set_secondary_coord(
        loci.start.getobasedpos()..loci.end.getobasedpos(),
        0..100usize,
    );
    secondary
        .configure_mesh()
        .y_label_formatter(&|f| format!("{f}%"))
        .x_label_formatter(&|f| {
            format!(
                "{} ({})",
                f.to_formatted_string(&Locale::en),
                pos.iter()
                    .find(|p| { p.position.getobasedpos() == *f })
                    .map(|p| p.locuspos.getobasedpos().to_formatted_string(&Locale::en))
                    .unwrap_or("Unknown".to_string())
            )
        })
        .x_desc("Genomic position (bp)")
        .y_desc("Mismatch rate (%)")
        .disable_x_mesh()
        .label_style(text_style.clone())
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Quality")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLACK));
    secondary
        .configure_secondary_axes()
        .y_desc("Quality (PHRED score)")
        //.disable_y_mesh()
        .label_style(text_style.clone())
        .draw()
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    if !args.nolegend {
        secondary
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    }
    //Bottom graph
    if let Some(bottom) = bottom {
        let max = pos
            .iter()
            .map(|f| f.globalmismatch)
            .max()
            .unwrap_or_default();
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
        let _ = chart
            .configure_mesh()
            .y_label_formatter(&|f| format!("{}%", (f * 10_000.0).round() / 100.0))
            .x_label_formatter(&|f| {
                format!(
                    "{} ({})",
                    f.to_formatted_string(&Locale::en),
                    pos.iter()
                        .find(|p| { p.position.getobasedpos() == *f })
                        .map(|p| p.locuspos.getobasedpos().to_formatted_string(&Locale::en))
                        .unwrap_or("Unknown".to_string())
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
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
                .map_err(|e| {
                    Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))
                })?;
        }
    }
    // To avoid the IO failure being ignored silently, we manually call the present function
    drawnoticetext(&root)?;
    root.present().map_err(|_| Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir")))?;
    Ok(())
}
pub(crate) fn locusisokay(mean: u64, graph: &[&HashMapinfo]) -> bool {
    //Between a minimum and a maximum number of reads
    graph.iter().all(|f| {
        f.overlaps >= MINIMUMCOVERAGE.try_into().unwrap_or(i64::MAX)
            && f.overlaps < ((mean as f32 * MAXCOVERAGERATIO).round() as i64)
    })
}
fn readgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    args: &Args,
    root: DrawingArea<T, Shift>,
    mean: u64,
) -> Result<(), Box<dyn std::error::Error>>
where
    T: DrawingBackend,
{
    let text_style = fontstyle.into_text_style(&root);
    let max = max(
        i64::try_from(mean).unwrap_or(i64::MAX),
        pos.iter().map(|max| max.getmaxvalue()).max().unwrap_or(0),
    ) + 5;
    let _ = root.fill(&plotters::prelude::WHITE);
    let (top, bottom) = root.split_vertically((80).percent_height());
    let tenlines = loci.getlength() / 9;
    let colorgene = if loci.status.isvalid() {
        full_palette::GREEN
    } else {
        full_palette::RED
    };
    let mut chart = ChartBuilder::on(&top)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(
            format!(
                "Reads mapping quality over the locus {} ({}-{})",
                loci.locus, loci.contig, loci.haplotype
            ),
            ("sans-serif", 28, &colorgene),
        )
        .build_cartesian_2d(loci.start.getobasedpos()..loci.end.getobasedpos(), 0..max)
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    let _ = chart
        .configure_mesh()
        .x_label_formatter(&|f| {
            format!(
                "{} ({})",
                f.to_formatted_string(&Locale::en),
                pos.iter()
                    .find(|p| { p.position.getobasedpos() == *f })
                    .map(|f| f.locuspos.getobasedpos().to_formatted_string(&Locale::en))
                    .unwrap_or_default()
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
        .draw_series(
            LineSeries::new(
                pos.iter().filter_map(|p| {
                    if (p.position.getobasedpos() - loci.start.getobasedpos())
                        .rem_euclid(std::cmp::max(tenlines, 1))
                        == 0
                    {
                        Some((
                            p.position.getobasedpos(),
                            mean.try_into().unwrap_or(i64::MAX),
                        ))
                    } else {
                        None
                    }
                }),
                full_palette::BROWN.mix(0.7).filled(),
            )
            .point_size(3),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Mean coverage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BROWN));
    chart
        .draw_series(
            AreaSeries::new(
                pos.iter().filter_map(|p| {
                    if (p.position.getobasedpos() - loci.start.getobasedpos())
                        .rem_euclid(std::cmp::max(tenlines, 1))
                        == 0
                    {
                        Some((
                            p.position.getobasedpos(),
                            (mean as f32 * MAXCOVERAGERATIO).round() as i64,
                        ))
                    } else {
                        None
                    }
                }),
                MINIMUMCOVERAGE.try_into().unwrap_or(i64::MAX),
                full_palette::GREY_A700.mix(0.4).filled(),
            )
            .border_style(full_palette::BLACK.mix(0.8)),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Coverage boundaries")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::GREY_A700.filled(),
            )
        });
    chart
        .draw_series(AreaSeries::new(
            pos.iter().map(|p| (p.position.getobasedpos(), p.map0)),
            0,
            full_palette::RED_300.mix(0.6),
        ))
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
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
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Secondary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLACK));
    chart
        .draw_series(LineSeries::new(
            pos.iter()
                .map(|p| (p.position.getobasedpos(), p.supplementary)),
            full_palette::BLUE_700,
        ))
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Supplementary alignments")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::BLUE_700));
    chart
        .draw_series(LineSeries::new(
            pos.iter().map(|p| (p.position.getobasedpos(), p.overlaps)),
            full_palette::ORANGE_300,
        ))
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Overlapping reads")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], full_palette::ORANGE_300));
    //Secondary
    /* let softclipmax = ((softclipmax as f32 /max as f32 * 10.0).ceil() / 10.0 * max as f32) as i64;
    //let softclipmax = core::cmp::max(calc,max);
     */
    /* let mut secondary = chart.set_secondary_coord(
        loci.start.getobasedpos()..loci.end.getobasedpos(),
        0..softclipmax,
    );
    secondary
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
        .x_desc("Genomic position (bp)")
        .y_desc("Average softclip")
        .disable_x_mesh()
        .label_style(text_style.clone())
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .unwrap();
    //let mut second = chart.set_secondary_coord(loci.start..loci.end, 0..max);
    secondary
        .draw_secondary_series(
            Histogram::vertical(&secondary)
                .baseline(0)
                .margin(3)
                .data(pos.iter().filter_map(|p| {
                    if p.softclips.is_normal() {
                        Some((
                            usize::try_from(p.position.getobasedpos()).unwrap(),
                            (p.softclips * 100f32) as i64,
                        ))
                    } else {
                        None
                    }
                }))
                .style(full_palette::BLACK.mix(0.4).filled()),
        )
        .unwrap()
        .label("Average softclip")
        .legend(|(x, y)| {
            plotters::element::Rectangle::new(
                [(x, y), (x + 15, y + 5)],
                full_palette::BLACK.filled(),
            )
        });
    secondary
        .configure_secondary_axes()
        .y_desc("Soft clips")
        .label_style(text_style.clone())
        //.disable_y_mesh()
        .draw()
        .unwrap();     */
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
    }
    //Bottom graph
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
            (loci.start.getobasedpos()..loci.end.getobasedpos()).into_segmented(),
            0..100i64,
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    let text_style = fontstyle.into_text_style(&root);
    let _ = chart
        .configure_mesh()
        .x_desc("Genomic position (bp)")
        .y_label_formatter(&|f| format!("{f}%"))
        .x_label_style(text_style.clone())
        .disable_x_axis()
        .disable_y_mesh()
        .draw();
    let breaks: Vec<(i64, i64)> = pos
        .iter()
        .filter_map(|elem| {
            if elem.overlaps <= args.breaks.into() {
                Some((elem.position.getobasedpos(), 100))
            } else {
                None
            }
        })
        .collect();
    let finalpos = pos
        .iter()
        .last()
        .unwrap_or_else(|| unreachable!("Invalid pos variable"))
        .position
        .getobasedpos();
    printbreaks(args, finalpos, &breaks, loci, outputfile)?;
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .style(full_palette::RED.filled())
                .data(breaks)
                .baseline(0)
                .margin(0),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Coverage break")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .baseline(0)
                .margin(3)
                .data(pos.iter().filter_map(|p| {
                    if p.softclips.is_normal() {
                        Some((p.position.getobasedpos(), (p.softclips * 100f32) as i64))
                    } else {
                        None
                    }
                }))
                .style(full_palette::BLACK.mix(0.8).filled()),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Softclips percent")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.6))
            .border_style(BLACK.mix(0.8))
            .label_font(text_style)
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    }
    drawnoticetext(&root)?;
    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().map_err(|_| Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir")))?;
    Ok(())
}

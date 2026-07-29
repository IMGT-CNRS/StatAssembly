/*
This software allows the analysis of BAM files to identify the confidence on a locus (specifically IG and TR) as well as allele confidence.
It was created and used by IMGT Team (https://www.imgt.org).
Available under EUPL license
Made by: Guilhem Zeitoun
*/
#![warn(clippy::unwrap_used)]
#![warn(clippy::expect_used)]
use crate::identification::{downloadref, locusallposition};
#[cfg(feature = "pdf")]
use crate::pdf::generatepdf;
#[cfg(feature = "bam")]
use crate::pilebam::pileup;
use crate::r#struct::AcceptedStatus::Rejected;
use bio::io::fasta;
use bio_types::sequence::SequenceRead;
use clap::Parser;
use itertools::Itertools;
use plotters::coord::Shift;
use plotters::coord::ranged1d::SegmentValue;
use plotters::style::full_palette::{BROWN_500, GREY_400};
use regex::Regex;
use std::collections::HashMap;
use std::io::ErrorKind::{InvalidData, InvalidInput};
use std::num::NonZero;
use std::ops::{Add, Div, Mul, RangeInclusive};
use std::path::Path;
use std::process::ExitCode;
use std::str::FromStr;
use std::time::Instant;
use std::{env, fs, io};
use strum::IntoEnumIterator;
//use noodles_fasta::{self as fasta, record::Sequence};
use crate::r#struct::*;
use crate::submissions::{
    GITHUBVERSION, INVALIDCOVERAGE, LIMITDATE, NOTENOUGHMATCHREADS, REQUESTCLIENT, SOFTCLIPTOOMUCH,
    SUSPICIOUSPOSITIONALERT, askforsubmission, checkifblastpresent, generatelightbam,
    generatesequence, genesblast, getindexforbam, getprogressbarclassic, getprogressbarspin,
    positionfiltering,
};
use extended_htslib::bam::ext::{BamRecordExtensions, IterAlignedPairs};
use extended_htslib::bam::record::{Cigar, CsValue};
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
#[cfg(feature = "pdf")]
mod pdf;
#[cfg(feature = "bam")]
mod pilebam;
mod r#struct;
mod submissions;
#[cfg(test)]
mod test;
const VERSION: &str = env!("CARGO_PKG_VERSION");
const AUTHOR: &str = "IMGT®";
const EMAIL: &str = "guilhem.zeitoun@cnrs.fr";
const GLOBALMISMATCHFLOATING: usize = 10_000;
const ALERTLOCUSSIZE: i64 = 10_000_000;
const MIN_READLENGTH: u64 = 1_000;
const MIN_PHREDSCORE: u64 = 30;
const TIMEOUT_IN_MN: u64 = 60;
const READGAPMESSAGE: u64 = 200;
const TELOMERESEP: u64 = 1_000; //Length when coverage is not taken into account.
const MINIMUMCOVERAGE: usize = 10;
const MAXCOVERAGERATIO: f32 = 2.0;
const PHYLUMLIMIT: usize = 7776; //obfstr::obf!("Gnathostomata");
const MATCHREADS: u32 = 10;
const BORNES: usize = 10_000;
const WARNINGPERC: u8 = 80;
const ALERTPERC: u8 = 60;
const SOFTCLIPRATIO: f32 = 0.4;
pub(crate) const DELIMITERFASTA: char = '/';
pub(crate) const LOCUSSEPARATOR: usize = 1_000_000;
lazy_static! {
    static ref fontstyle: (&'static str, u32, &'static RGBColor) = {
        let args = Args::parse();
        ("sans-serif", args.fontlegendsize.get(), &BLACK)
    };
    static ref NAME: String = env!("CARGO_PKG_NAME").replacen('_', "/", 1);
    #[allow(clippy::unwrap_used)]
    static ref regexpword: regex::Regex = regex::Regex::new(r"[^-\w()]").unwrap_or_else(|_| unreachable!("Regex issue"));
}
//Return block of positions thanks to CS/MD tag or CIGAR = (preferred if existing)
fn iterblock(record: &bam::Record) -> Option<Vec<[i64; 2]>> {
    match (
        record.aligned_blocks_match(),
        record.getorgeneratecsblock(false),
    ) {
        //There is a CIGAR =
        (Some(d), _) => Some(d.collect()),
        //There is a MD/CS tag
        (None, Ok(d)) => Some(
            d.into_iter()
                .filter_map(|p| match (&p.state, p.getgenomepos()) {
                    (CsValue::Same(d), Some(pos)) => Some([
                        pos,
                        pos.checked_add(d.getsize().checked_sub(1)?.try_into().ok()?)
                            .unwrap_or(0),
                    ]),
                    _ => None,
                })
                .collect(),
        ),
        //We have nothing
        (None, Err(e)) => {
            eprintln!("No cigar X/= and MD/CS tag cannot be parsed. Error is {e}");
            None
        }
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
            start: Position::newfromoposition( min(posa, posb)),
            end: Position::newfromoposition( max(posa, posb)),
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
    if !args.allreads && !record.is_primary() {
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
            format!("Assembly error, maybe index is missing (create it with samtools faidx): {e}"),
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
            .find(|p| p.iter().any(|f| f.getlocus() == loci.getlocus()))
        {
            Some(d) if Some(d) != elem.last() => {
                eprintln!("Locus {} is splited! Aborted.", loci.getlocus());
                return None;
            }
            _ => (),
        };
        let alternate = if !loci.gethaplotype().isprimary() {
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
            Some(e) if e.getlocus() == loci.getlocus() && alternate && actual.len() >= 2 => {
                eprintln!("Only one alternate is allowed!");
                return None;
            }
            Some(e) if e.getlocus() == loci.getlocus() && alternate => actual.push(loci.clone()),
            Some(e) if e.getlocus() != loci.getlocus() && alternate => {
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
) -> std::io::Result<(Vec<LocusInfos>, Option<Vec<Blastmatch>>, String)> {
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
                Ok(mut r) => {
                    r.checkandresetifauto();
                    r
                }
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
            records.retain(|p| p.haplotype.as_ref().is_some_and(|f| f.isprimary()));
        }
    }
    let mut locusrecord: (Vec<LocusInfos>, Option<Vec<Blastmatch>>) = (Vec::new(), None);
    let mut release = String::new();
    match (
        args.assembly.as_ref(),
        records
            .iter_mut()
            .filter(|p| p.contig.is_none())
            .collect::<Vec<&mut FakeLocusinfo>>(),
    ) {
        (Some(path), b) if !b.is_empty() && blastpresent => {
            let (locus, blast, release4) = locusallposition(path, realspecies, args)?;
            release = release4;
            b.into_iter().for_each(|p| {
                let find = locus.iter().find(|r| {
                    r.gethaplotype() == &p.haplotype.clone().unwrap_or_default()
                        && r.getlocus() == &p.locus
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
            elem.locus,
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
                d.getlocus(),
                d.gethaplotype(),
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
                "The region {}-{} ({}) from {} is more than {} bp (actually: {}) and might be incorrect. Check your input file ({}) is correct.\nIf wanted, please add --hugeregion parameters. Be careful software might use a lot of memory.",
                big.start.getobasedpos(),
                big.end.getobasedpos(),
                big.contig,
                big.getlocus(),
                ALERTLOCUSSIZE.to_formatted_string(&Locale::en),
                big.getlength().to_formatted_string(&Locale::en),
                args.locuspos
                    .as_ref()
                    .map(|f| format!("{}", f.display()))
                    .unwrap_or("no loc given".to_string())
            ),
        ));
    }
    Ok((locusrecord.0, locusrecord.1, release))
}
fn generategenelist<T>(
    locushashresult: &Option<Vec<Blastmatch>>,
    speciesblast: &Species,
    locus: &[LocusInfos],
    assembly: T,
    args: &Args,
) -> io::Result<(Option<Vec<Blastmatch>>, Option<String>)>
where
    T: AsRef<Path>,
{
    let func = || {
        eprintln!("You have not provided a gene list, BLASTING to get one.");
        locusallposition(assembly.as_ref(), speciesblast, args).map(|(_, b, c)| (b, c))
    };
    let (mut data, release) = match locushashresult {
        None => {
            let a = func()?;
            (a.0, Some(a.1))
        }
        Some(a) if a.is_empty() => {
            let a = func()?;
            (a.0, Some(a.1))
        }
        Some(b) => {
            eprintln!("Generating a gene list.");
            (b.to_vec(), None)
        }
    };
    //locusfiltering(&loci.getlocus(), &mut data);
    positionfiltering(locus, &mut data);
    if data.is_empty() {
        eprintln!("No data after filtering gene list. Skipped.");
    } else {
        return Ok((Some(data), release));
    }
    Ok((None, None))
}
fn getgeneinfos(data: &[Blastmatch]) -> Vec<GeneInfos> {
    data.iter()
        .filter_map(|p| {
            let strand = if p.send < p.sstart {
                Strand::Minus
            } else {
                Strand::Plus
            };
            let gene = match Genename::new(&p.getquery().gene, p.getquery().label.clone()) {
                Ok(a) => a,
                Err(e) => {
                    eprintln!("Error with gene {}: {e}", p.getquery().gene);
                    return None;
                }
            };
            let subject = match Name::from_str(p.getsubject()) {
                Ok(a) => a.numacc.unwrap_or_default(),
                Err(_) => p.getsubject().to_string(),
            };
            Some(GeneInfos::new(
                gene,
                subject,
                strand,
                Position::newfromoposition(p.sstart.try_into().unwrap_or_default()),
                Position::newfromoposition(p.send.try_into().unwrap_or_default()),
            ))
        })
        .collect()
}
fn printgenelist(data: &[Blastmatch], args: &mut Args, tmp: bool) -> io::Result<Option<Filecrea>> {
    let mut finish: Vec<GeneInfos> = getgeneinfos(data);
    checkandcorrectgenelistduplicate(&mut finish);
    let genenamefile = if tmp {
        Filecrea::createtemp(None, Some("genelist_new.csv"))?
    } else {
        Filecrea::createfrompath(args.outdir.join("genelist_new.csv"))
    };
    let genefile = File::create(genenamefile.getpath())?;
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
    args.geneloc = Some(genenamefile.getpath().to_path_buf());
    Ok(Some(genenamefile))
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
        getassemblyreader(args)
            .map_err(|p| io::Error::new(p.kind(), format!("Cannot open assembly. Error is {p}")))?;
    }
    let outputdir = match args.outdir.is_dir() {
        true => {
            eprintln!(
                "Output folder exists: {}. Will be overwritten.",
                args.outdir.display()
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
    let start = &Position::newfromzposition(record.reference_start());
    let end = &Position::newfromzposition(record.reference_end());
    //Get range to put the reads inclusive pos
    let newrange = Position::new(
        true,
        max(record.reference_start(), locus.start.getzbasedpos()),
    )
        ..Position::newfromzposition(min(locus.end.getobasedpos(), record.reference_end()));
    let hitter = if let Some(d) = pos.get_mut(start) {
        Some(d)
    } else if let Some(d) = pos.get_mut(end) {
        Some(d)
    } else {
        None
    };
    match ( //Leading and trailing are only present at the start or end so no need to check that.
        hitter,
        record.cigar().leading_softclips() > 0
            || record.cigar().trailing_softclips() > 0
            || record.cigar().leading_hardclips() > 0
            || record.cigar().trailing_hardclips() > 0,
        record.is_primary(),
    ) {
        (_, false, _) | (None, ..) => (),
        (Some(d), true, true) => d.psoftclips += 1.0,
        (Some(d), true, false) => d.osoftclips += 1.0,
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
        targeting.total += 1;
        if !record.is_primary() {
            if record.is_secondary() {
                targeting.secondary += 1;
            } else if record.is_supplementary() {
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
                targeting.qual = Some(
                    targeting.qual.unwrap_or_default()
                        + record.qual().get(index).map_or(0, |f| (*f).into()),
                );
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
fn getmeancoveragelengthandphred(args: &Args) -> io::Result<Params> {
    let time = Instant::now();
    let mut reader: IndexedReader =
        getreaderoffile(args).map_err(|f| io::Error::new(io::ErrorKind::InvalidInput, f))?;
    reader
        .fetch(FetchDefinition::All)
        .map_err(|f| io::Error::new(io::ErrorKind::InvalidInput, f))?;
    let mut readlength = 0;
    let mut count = 0;
    let mut phred = 0;
    let length = reader
        .rc_records()
        .filter_map(Result::ok)
        .filter_map(|f| {
            //Remove reads with 0 as MAPQ
            if f.mapq() != 0 && !f.is_unmapped() && f.is_primary() {
                readlength += f.len();
                count += 1;
                phred += f
                    .qual()
                    .iter()
                    .map(|a| <u8 as Into<u64>>::into(*a))
                    .sum::<u64>()
                    / f.len().try_into().unwrap_or(1);
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
    let phred = phred / count;
    println!("Calculated in {:.1} s", time.elapsed().as_secs());
    Ok(Params::new(readavg, mean, phred))
}
fn getorsetparams(meanpath: &Path, args: &Args) -> io::Result<Params> {
    let params = match (
        args.paramsfile
            .as_ref()
            .map(|b| fs::read_to_string(b).map(|p| Params::from_str(&p))),
        args.extractedlength,
        fs::read_to_string(meanpath).map(|p| Params::from_str(&p)),
    ) {
        /* //Mean is given
        (Some(mean), ..) => {
            println!(
                "Mean was provided, getting given value and setting the value for further usage."
            );
            let _ = fs::write(meanpath, mean.to_string().as_bytes());
            Ok(mean)
        } */
        //The mean was calculated with a file
        (Some(Ok(Ok(params))), b, _) | (.., b, Ok(Ok(params))) if b == 0 => {
            println!("Params were already calculated, retrieving...");
            params
        }
        _ => {
            println!("Getting mean coverage from calculations, it might take some minutes.");
            let params = getmeancoveragelengthandphred(args)?;
            println!(
                "Mean coverage is {} and average length of reads is {}. PHRED score is {}",
                params.getmean(),
                params.getavg(),
                params.getphred()
            );
            let _ = fs::write(meanpath, params.to_string().as_bytes());
            params
        }
    };
    params.goodparams().map(|_| params)
}
#[allow(clippy::complexity)]
fn paintgraph<T>(
    loci: &LocusInfos,
    pos: &BTreeMap<Position, HashMapinfo>,
    args: &Args,
    species: &Species,
    readgraphelem: DrawingArea<T, Shift>,
    mismatchgraphelem: DrawingArea<T, Shift>,
    mean: u64,
) -> io::Result<()>
where
    T: DrawingBackend,
{
    if let Err(e) = mismatchgraph(
        &args.outdir,
        loci,
        pos.values().collect_vec().as_slice(),
        args,
        mismatchgraphelem,
    )
    .and_then(|_| {
        readgraph(
            &args.outdir,
            loci,
            &pos.iter().map(|a| a.1).collect_vec(),
            args,
            species,
            readgraphelem,
            mean,
        )
    }) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{}", e),
        ));
    };
    Ok(())
}
fn available_memory() -> Option<usize> {
    let contents = std::fs::read_to_string("/proc/meminfo").ok()?;
    let mem_info = contents
        .lines()
        .find(|line| line.starts_with("MemAvailable"))?;
    let size = mem_info.split(" ").nth(4)?;
    let available_mem: usize = size.parse().ok()?;
    Some(available_mem)
}
/*
fn setgraphbitmap<'a>(
    haplotype: usize,
    haplotypebool: bool,
    args: &'a Args,
    outputdir: &Path,
    speciesblast: &Species,
    floci: &LocusInfos,
) -> Result<HashMap<Haplotype, Vec<(String, ImageType)>>, String> {
    let readresultsize = (
        1800,
        1000 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap_or(1),
    );
    let mismatchsize = (
        1200,
        600 * std::convert::TryInto::<u32>::try_into(haplotype).unwrap_or(1),
    );
    let elem = if haplotypebool {
        vec![Haplotype::Primary]
    } else {
        Haplotype::iter().collect()
    };
    let mut list: HashMap<Haplotype, Vec<(String, ImageType)>> = HashMap::new();
    //Need to peuplate the HashMap before else Image has bad lifetime
    for haplo in elem.iter() {
        for fgraph in ["readresult", "mismatchgraph"] {
            let ext = if args.svg { "svg" } else { "png" };
            let outputfile = outputdir.join(givename(
                &speciesblast,
                &floci.getlocus(),
                &floci.contig,
                haplotypebool,
                &format!("{fgraph}.{ext}"),
                true,
            ));
            let size = if fgraph == "readresult" {
                readresultsize
            } else {
                mismatchsize
            };
            let root = if args.svg {
                ImageType::Svg((outputfile, size))
            } else {
                ImageType::Png((outputfile, size))
            };
            if let Some(b) = list.get_mut(&haplo) {
                b.push((fgraph.to_string(), root));
            } else {
                list.insert(haplo.clone(), vec![(fgraph.to_string(), root)]);
            }
        }
    }
    Ok(list)
}
*/
fn posread(args: &Args, loci: &LocusInfos) -> io::Result<BTreeMap<Position, HashMapinfo>> {
    let mut pos: BTreeMap<Position, HashMapinfo> = BTreeMap::new();
    let mut reader = match getreaderoffile(&args) {
        Ok(r) => r,
        Err(e) => {
            return Err(io::Error::new(
                InvalidData,
                format!("Cannot read bam file. Error is {e}. Exiting"),
            ));
        }
    };
    //0-based except end because end is exclusive
    if loci.fetch(&mut reader).is_err() {
        eprintln!(
            "The region {}:{}-{} cannot be found, skipped.",
            loci.contig,
            loci.start.getobasedpos(),
            loci.end.getobasedpos()
        );
        return Ok(pos);
        //return ExitCode::FAILURE;
    }
    let mut nocount = true;
    //let filename = outputdir.join(format!("{}.pileup", &loci.locus));
    //let file = File::create(&filename).unwrap();
    //let mut writer = BufWriter::new(file);
    let locusrange = loci.start.getzbasedpos()..=loci.end.getzbasedpos();
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
        let default = HashMapinfo::smalldefault(Position::newfromzposition(i),Position::newfromzposition(p));
        pos.insert(Position::newfromzposition(p), default);
    });
    let mut message = false;
    println!("Region {} fetched, analyzing all reads.", loci.getlocus());
    let mut count = 0;
    //let time = Instant::now();
    let sep = if args.fullquality {
        1
    } else {
        max(
            (loci.end.getobasedpos() - loci.start.getobasedpos() + 1) / 250,
            100,
        ) //250 points for quality point or 100nt break
    };
    let progress = getprogressbarspin();
    for p in reader
        .rc_records()
        .filter_map(Result::ok)
        .filter(|p| !(args.forward && p.is_reverse()))
    {
        count += 1;
        //Print every x reads done
        if count % READGAPMESSAGE == 0
            && let Ok(a) = &progress
        {
            a.set_position(count);
            /* let _ = writeln!(
                lock,
                "Processed {} reads in {:.3} s",
                count.to_formatted_string(&Locale::en),
                Instant::now().saturating_duration_since(time).as_secs_f32()
            ); */
        }
        nocount = false;
        match processcounting(&args, &mut pos, message, loci, &p, sep) {
            Err(e) => {
                return Err(io::Error::new(io::ErrorKind::InvalidData, e));
            }
            Ok(b) => {
                if !message && b {
                    message = b;
                }
            }
        }
    }
    if let Ok(a) = progress {
        a.finish_and_clear();
    }
    if nocount {
        eprintln!(
            "The region {}:{}-{} has no data, skipped.",
            loci.contig,
            loci.start.getobasedpos(),
            loci.end.getobasedpos()
        );
        return Ok(pos);
        //return ExitCode::FAILURE;
    }
    //Quality and mismatch is the sum of reads so dividing to get real results
    pos.iter_mut().for_each(|(_, p)| {
        if let Some(q) = p.qual
            && q > 0
        {
            p.qual = Some(
                q / max(
                    std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap_or(1),
                    1,
                ),
            )
        }
        if p.psoftclips.is_normal() {
            p.psoftclips = (p.psoftclips * 100f32 / max(p.gettotalmap(), 1) as f32).round() / 100f32
        }
        if p.osoftclips.is_normal() {
            p.osoftclips = (p.psoftclips * 100f32 / max(p.getfullsecondary().unwrap_or(1), 1) as f32).round() / 100f32
        }
        p.globalmismatch /= max(
            std::convert::TryInto::<usize>::try_into(p.gettotalmap()).unwrap_or(1),
            1,
        );
    });
    Ok(pos)
}
fn main() -> ExitCode {
    /*
    let mainpa = env::current_dir().unwrap();
    let testo = tempfile::TempDir::new().unwrap();
    let path = &Path::new(&mainpa).join("full_cs.bam").display().to_string();
    let path2 = Path::new(&mainpa)
        .join("assembly.fasta")
        .display()
        .to_string();
    let path3 = Path::new(&mainpa).join("locus.txt").display().to_string();
    let fake_args = vec![
        "IMGT_StatAssembly",
        "-f",
        path,
        "-a",
        &path2,
        "-l",
        &path3,
        "-s",
        "human",
        "-o",
        testo.path().to_str().unwrap(),
        "analyze",
    ];
    let human = Species::new("Homo sapiens").unwrap();
    let mut args4 = Args::try_parse_from(fake_args)
        .map_err(|f| f.to_string())
        .unwrap();
    let (locus, ..) = locusposparser(&args4, &human, false).unwrap();
    let mut append = std::fs::OpenOptions::new()
        .append(true)
        .open("/home/gzeitoun/historic.txt")
        .unwrap();
    let string = locus.iter().fold(String::new(), |mut acc, loci| {
        acc.push_str(&format!(
            "{} - Locus {} is invalid at {}\n",
            mainpa.display().to_string(),
            loci.getlocushaplo(),
            loci.contig
        ));
        acc
    });
    //let s = string.trim();
    let _ = append.lock(); //Would fail on Windows if append so reject the error
    append.write_all(string.as_bytes()).unwrap();
    let _ = append.unlock();
    let (hits, _) = generategenelist(
        &None,
        &human,
        &locus,
        args4.assembly.as_ref().unwrap(),
        &args4,
    )
    .unwrap();
    let listtokkep = printgenelist(&hits.unwrap(), &mut args4, true).unwrap(); //So the file is not dropped.
    for loci in locus {
        let mut genes = extractgenelist(&args4, &loci, true).unwrap();
        let invalidgenes = genes
            .iter_mut()
            .map(|mut a| generategeneinfos(&args4, &mut a).unwrap())
            .filter(|p| !p.0.getstatus().getstatus().isvalid());
        let string = invalidgenes.fold(String::new(), |mut acc, (gene, _)| {
            acc.push_str(&format!(
                "{} - {} is invalid at {}:{}-{}\n",
                mainpa.display().to_string(),
                gene.getgene(),
                gene.getchromosome(),
                gene.getstart().getobasedpos(),
                gene.getend().getobasedpos()
            ));
            acc
        });
        //let s = string.trim();
        let _ = append.lock(); //Would fail on Windows if append so reject the error
        append.write_all(string.as_bytes()).unwrap();
        let _ = append.unlock();
    }
    return ExitCode::SUCCESS;
     */
    let mut args = Args::parse();
    if checknewversion() {
        return ExitCode::FAILURE;
    }
    if !args.lowmemory
        && let Some(b) = available_memory()
        && b < 12_usize.saturating_mul(10_usize.saturating_pow(6))
    {
        println!("You have low memory available.");
        args.lowmemory = true; //Force low memory if less than 12 Go of free memory (in kb normally in /proc/meminfo).
    }
    let firstinstant = Instant::now();
    let speciesblast = match Species::new(&args.species) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("Species {} has failed: {e}", args.species);
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{} is {} (taxon: {}).",
        speciesblast
            .getrank()
            .chars()
            .take(1)
            .map(|p| p.to_ascii_uppercase())
            .take_while(|_| true)
            .collect::<String>(),
        speciesblast.getname(),
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
            eprintln!("Error with locus file, assembly or output: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (locus, blastcheck, releaseversion) =
        match locusposparser(&args, &speciesblast, blastpresent) {
            Err(f) => {
                eprintln!("Error locus parser: {f}");
                return ExitCode::FAILURE;
            }
            Ok(b) => b,
        };
    let mut releaseversion = if releaseversion.trim().is_empty() {
        None
    } else {
        Some(releaseversion.trim().to_string())
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
        /* todo!("To delete");
        let infos: HashMap<Locus, Vec<Blastmatch>> = if let Some(blastmatches) = blastcheck {
            blastmatches
                .into_iter()
                .filter_map(|f| {
                    if let Some(a) = locus.iter().find(|loci| {
                        f.sseqid == loci.contig
                            && (loci.start.getobasedpos()..=loci.end.getobasedpos())
                                .contains(&f.sstart.try_into().unwrap_or_default())
                    }) {
                        Some((a.locus.clone(), f))
                    } else {
                        None
                    }
                })
                .into_group_map_by(|(locus, _)| locus.clone())
                .into_iter()
                .map(|(locus, matches)| {
                    (
                        locus,
                        matches
                            .into_iter()
                            .map(|(_, blast_match)| blast_match)
                            .collect(),
                    )
                })
                .collect()
        } else {
            HashMap::new()
        }; */
        if blastcheck.as_ref().is_none_or(|a| a.is_empty()) {
            eprintln!("Error setting gene list result. Empty list.");
            return ExitCode::FAILURE;
        }
        if let Some(b) = r {
            let v = match generategenelist(&blastcheck, &speciesblast, &locus, b, &args) {
                Ok((Some(a), b)) => {
                    #[allow(unused_assignments)]
                    if let Some(release) = b {
                        releaseversion = Some(release);
                    }
                    printgenelist(&a, &mut args, false).map(|_| ())
                }
                Ok((None, _)) => Ok(()),
                Err(a) => Err(a),
            };
            if let Err(e) = v {
                eprintln!("Error setting gene list result: {e}");
                return ExitCode::FAILURE;
            }
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
    let mean = match getorsetparams(&meanpath, &args) {
        Ok(a) => a.getmean(),
        Err(e) => {
            eprintln!("Error making parameters: {e}.");
            return ExitCode::FAILURE;
        }
    };
    let mut locushashresult: HashMap<LocusInfos, Vec<Blastmatch>> = HashMap::new();
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
            eprintln!("There is more than 2 haplotypes for {}", floci.getlocus());
            return ExitCode::FAILURE;
        }
        let haplotypebool = haplotype == 1; //IS there one or two haplotypes?
        println!(
            "Going for {} locus - {}",
            floci.getlocus(),
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
        let loc = if let Some(a) = locus.first() {
            a
        } else {
            eprintln!("Error with some locus, continuing");
            continue;
        };
        let readname = Path::join(
            &args.outdir,
            givename(
                &speciesblast,
                loc.getlocus(),
                &loc.contig,
                haplotypebool,
                &format!("readresult.{}", if args.svg { "svg" } else { "png" }),
                true,
            ),
        );
        let mismatchname = Path::join(
            &args.outdir,
            givename(
                &speciesblast,
                loc.getlocus(),
                &loc.contig,
                haplotypebool,
                &format!("mismatchresult.{}", if args.svg { "svg" } else { "png" }),
                true,
            ),
        );
        let mut drawingmaps = {
            if args.svg {
                let (top, bottom) = if !haplotypebool {
                    let r = SVGBackend::new(readname.as_path(), readresultsize)
                        .into_drawing_area()
                        .split_vertically((50).percent_height());
                    (r.0, Some(r.1))
                } else {
                    let r = SVGBackend::new(readname.as_path(), readresultsize).into_drawing_area();
                    (r, None)
                };
                let (mtop, mbottom) = if !haplotypebool {
                    let r = SVGBackend::new(mismatchname.as_path(), mismatchsize)
                        .into_drawing_area()
                        .split_vertically((50).percent_height());
                    (r.0, Some(r.1))
                } else {
                    let r =
                        SVGBackend::new(mismatchname.as_path(), mismatchsize).into_drawing_area();
                    (r, None)
                };
                DrawingImage {
                    readresulttop: Some(Backend::Svg(top)),
                    readresultbottom: bottom.map(Backend::Svg),
                    mismatchtop: Some(Backend::Svg(mtop)),
                    mismatchbottom: mbottom.map(Backend::Svg),
                }
            } else {
                let (top, bottom) = if !haplotypebool {
                    let r = BitMapBackend::new(readname.as_path(), readresultsize)
                        .into_drawing_area()
                        .split_vertically((50).percent_height());
                    (r.0, Some(r.1))
                } else {
                    let r =
                        BitMapBackend::new(readname.as_path(), readresultsize).into_drawing_area();
                    (r, None)
                };
                let (mtop, mbottom) = if !haplotypebool {
                    let r = BitMapBackend::new(mismatchname.as_path(), mismatchsize)
                        .into_drawing_area()
                        .split_vertically((50).percent_height());
                    (r.0, Some(r.1))
                } else {
                    let r = BitMapBackend::new(mismatchname.as_path(), mismatchsize)
                        .into_drawing_area();
                    (r, None)
                };
                DrawingImage {
                    readresulttop: Some(Backend::Png(top)),
                    readresultbottom: bottom.map(Backend::Png),
                    mismatchtop: Some(Backend::Png(mtop)),
                    mismatchbottom: mbottom.map(Backend::Png),
                }
            }
        };
        //let mut lock = stdout().lock();
        //For each individual haplotype inside locus
        for loci in locus.iter_mut() {
            let pos = match posread(&args, loci) {
                Ok(b) if b.is_empty() => continue,
                Ok(a) => a,
                Err(e) => {
                    eprintln!("Error setting counts: {e}");
                    return ExitCode::FAILURE;
                }
            };
            loci.setstatus(mean, &args, &pos);
            //let _ = progress.map(|d| d.finish());
            println!("Making graphs");
            {
                /*
                let imagepng = |path: &PathBuf, size: (u32, u32)| -> BitMapBackend<'a, _> {
                    if !haplotypebool {
                        let (a, b) = BitMapBackend::new(path.as_path(), size)
                            .into_drawing_area()
                            .split_vertically(50);
                        if loci.gethaplotype().isprimary() { a } else { b }
                    } else {
                        let i = BitMapBackend::new(path.as_path(), size).into_drawing_area();
                        i
                    }
                };
                let imagesvg = |path: &PathBuf, size: (u32, u32)| {
                    if !haplotypebool {
                        let (a, b) = SVGBackend::new(path.as_path(), size)
                            .into_drawing_area()
                            .split_vertically(50);
                        if loci.gethaplotype().isprimary() { a } else { b }
                    } else {
                        let i = SVGBackend::new(path.as_path(), size).into_drawing_area();
                        i
                    }
                };
                */
                let (readmap, mismatchmap, new) = if loci.gethaplotype().isprimary() {
                    drawingmaps.gettop()
                } else {
                    drawingmaps.getbottom()
                };
                drawingmaps = new;
                match (readmap, mismatchmap) {
                    (Some(a), Some(b)) => {
                        if a.issvg() {
                            if let (Some(a), Some(b)) = (a.getsvg(), b.getsvg()) {
                                if let Err(e) =
                                    paintgraph(loci, &pos, &args, &speciesblast, a, b, mean)
                                {
                                    eprintln!("Error drawing graphs. Error is {e}");
                                    continue;
                                }
                            } else {
                                eprintln!("Error drawing graphs.");
                                continue;
                            }
                        } else {
                            if let (Some(a), Some(b)) = (a.getpng(), b.getpng()) {
                                if let Err(e) =
                                    paintgraph(loci, &pos, &args, &speciesblast, a, b, mean)
                                {
                                    eprintln!("Error drawing graphs. Error is {e}");
                                    continue;
                                }
                            } else {
                                eprintln!("Error drawing graphs.");
                                continue;
                            }
                        }
                    }
                    _ => {
                        eprintln!("Error making graphs.");
                        continue;
                    }
                }
            }
            println!("Graphs finished.");
            //Create CSV from HashMap
            if let Err(e) = createcsv(
                outputdir.as_path(),
                &speciesblast,
                loci,
                pos.values().collect_vec().as_slice(),
                &args,
            ) {
                eprintln!("Cannot create csv file. Error is {e}");
                return ExitCode::FAILURE;
            }
            println!(
                "Locus {} ({}) is {}.",
                loci.getlocus(),
                loci.gethaplotype(),
                loci.status
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
                locushashresult.insert(loci.clone(), element);
            }
            if args.geneloc.is_some() {
                println!("Gene list starting!");
                let _filepath = if let Some(a) = blastcheck.as_ref()
                    && !args
                        .geneloc
                        .as_ref()
                        .is_some_and(|f| f.try_exists().unwrap_or(false))
                {
                    match printgenelist(a, &mut args, true) {
                        Ok(a) => a,
                        Err(e) => {
                            eprintln!("Cannot create gene list. Error is {e}");
                            return ExitCode::FAILURE;
                        }
                    }
                } else {
                    None
                };
                match genelist(loci, &speciesblast, &args, false) {
                    Err(e) => {
                        eprintln!("Cannot create gene list. Error is {e}");
                        return ExitCode::FAILURE;
                    }
                    Ok(b) if b.is_empty() => {
                        eprintln!("No gene list for locus {}. Skipped.", loci.getlocus());
                    }
                    Ok(b) => {
                        if args.assembly.as_ref().is_some() && !args.nosubmit && blastpresent {
                            println!("Blasting gene list");
                            let geneinfo: Vec<GeneInfos> =
                                b.iter().map(|f| f.clone().into()).collect();
                            match genesblast(
                                &geneinfo,
                                &args,
                                &releaseversion,
                                &speciesblast,
                                loci.getlocus(),
                            ) {
                                Ok((blast, r)) => {
                                    if releaseversion.is_none() {
                                        releaseversion = r;
                                    }
                                    if loci.status.getstatus().isvalid()
                                        && b.iter().any(|f| f.status.getstatus().isvalid())
                                    {
                                        locushashresult.insert(loci.clone(), blast);
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
                        &blastcheck,
                        &speciesblast,
                        initiallocus,
                        b,
                        &args,
                    ) {
                        Ok((Some(a), b)) => {
                            if let Some(release) = b {
                                releaseversion = Some(release);
                            }
                            a
                        }
                        Ok((None, _)) => continue,
                        Err(e) => {
                            eprintln!(
                                "Gene list generation for locus {} has an error: {}.",
                                loci.getlocus(),
                                e
                            );
                            continue;
                        }
                    };
                    let locivec = vec![loci.clone()];
                    positionfiltering(&locivec, &mut data);
                    let result = match (
                        printgenelist(&data, &mut args, true),
                        genelist(loci, &speciesblast, &args, false),
                    ) {
                        (Err(e), _) => {
                            eprintln!(
                                "Gene list generation for locus {} has an error: {}.",
                                loci.getlocus(),
                                e
                            );
                            continue;
                        }
                        (Ok(_a), Err(e)) => {
                            eprintln!(
                                "Gene list generation for locus {} has an error: {}.",
                                loci.getlocus(),
                                e
                            );
                            continue;
                        }
                        (Ok(_), Ok(b)) => b,
                    };
                    data.retain(|a| {
                        result
                            .iter()
                            .find(|c| {
                                Genename::new(a.qseqid.gene.clone(), a.qseqid.label.clone())
                                    .is_ok_and(|a| a == c.gene)
                            })
                            .is_some_and(|b| b.status.getstatus().isvalid())
                    });
                    if loci.status.getstatus().isvalid()
                        && result.iter().any(|f| f.status.getstatus().isvalid())
                    {
                        locushashresult.insert(loci.clone(), data);
                    }
                } else if !args.nosubmit {
                    eprintln!("No assembly to check gene list.");
                }
            }
            #[cfg(feature = "bam")]
            if let Err(e) = pileup(loci, &args) {
                eprintln!("Cannot make pileup. Error is {e}");
            }
        }
        println!("Locus {} is done!", floci.getlocus());
    }
    let mergedloci: Vec<LocusInfos> = grouped.iter().flatten().cloned().collect();
    if let Some(light) = &args.outlightbam
        && let Err(e) = generatelightbam(&args, light, getindexforbam(&light).as_ref(), &mergedloci)
    {
        eprintln!("{e}");
        return ExitCode::FAILURE;
    }
    if let Err(e) = printnewloc(&args, &mergedloci) {
        eprintln!("Error setting new locus result: {e}");
    }
    if let Some(a) = blastcheck
        && let Err(e) = printgenelist(&a, &mut args, false)
    {
        eprintln!("Error setting new gene list: {e}");
    }
    if !locushashresult.is_empty()
        && let Err(e) = printvalidatedalleles(
            &args,
            downloadref(false, args.cacheerase, &releaseversion).map(|(_, release)| release),
            &locushashresult,
        )
    {
        eprintln!("Error setting alleles result: {e}");
    }
    if let Err(e) = generatesequence(&args, &args.outdir, false, &mergedloci) {
        eprintln!("{e}");
    }
    #[cfg(feature = "pdf")]
    if let Err(e) = generatepdf(&mergedloci, &locushashresult, &args.outdir) {
        eprintln!("{e}");
    }
    if !args.nosubmit
    /* && locushashresult
    .iter()
    .any(|(_, f)| f.iter().any(|g| g.onlynewalleles())) */
    {
        locushashresult
            .values_mut()
            .for_each(|a| a.retain(|p| p.onlynewalleles()));
        locushashresult.retain(|_, v| !v.is_empty());
        if !locushashresult.is_empty() {
            match askforsubmission(&speciesblast, &mergedloci, &args, &locushashresult) {
                Ok(_) => (),
                Err(e) => eprintln!("Error submitting sequences: {e}"),
            }
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
        csv.write_record([
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
fn printvalidatedalleles<T>(
    args: &Args,
    release: Option<T>,
    locushash: &HashMap<LocusInfos, Vec<Blastmatch>>,
) -> io::Result<()>
where
    T: AsRef<str>,
{
    let path = args.outdir.join("validatedalleles.fasta");
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
        "subject",
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
            let (subject, start, end, strand) = if let Ok(r) = Name::from_str(&matches.sseqid) //Check if it is a full composite or not
                && let Some(a) = r.numacc && let Some((start,end,strand)) = locus.fullposition(matches)
            {
                (a, start, end, strand)
            } else {
                (
                    matches.sseqid.clone(),
                    Position::newfromoposition(matches.sstart.try_into().unwrap_or_default()),
                    Position::newfromoposition(matches.send.try_into().unwrap_or_default()),
                    matches.complement.clone(),
                )
            };
            csv.write_record(&[
                name.to_string(),
                format!(
                    "{}-{}{}",
                    start.getobasedpos(),
                    end.getobasedpos(),
                    if strand.isrev() { "/rc" } else { "" }
                ),
                subject,
                LocusHaplo::from(locus.clone()).to_string(),
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
        .duplicates_by(|g| (&g.gene, &g.chromosome))
        .collect();
    if !finish.is_empty() {
        println!(
            "Some genes has the same name on the same chromosome. Underscore (_) would be added"
        );
        let mut count = 0;
        for name in genes.iter_mut() {
            name.gene.name = if finish.iter().any(|g| g.gene.eq(&name.gene)) {
                count += 1;
                format!("{}_{}", name.gene.name, count)
                    .replace(",", "_")
                    .trim()
                    .to_string()
            } else {
                name.gene.name.replace(",", "_").trim().to_string()
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
        println!("No gene identified for locus {}, skipped.", loci.getlocus());
        return Ok(Vec::new());
    }
    //At least one duplicate line
    checkandcorrectgenelistduplicate(&mut genes);
    Ok(genes)
}
fn generategeneinfos(
    args: &Args,
    gene: &mut GeneInfos,
) -> io::Result<(GeneInfosFinish, BTreeMap<Position, Posread>)> {
    let mut lock = io::stdout().lock();
    let mut reader =
        getreaderoffile(args).map_err(|e| io::Error::new(io::ErrorKind::BrokenPipe, e))?;
    let (mut reads, mut readsfull, mut reads100, mut reads100m, mut realreads100m, mut phredscore) =
        (0, 0, 0, 0, 0f32, Vec::new());
    let (mut hash, records) = {
        //O position is exclusive
        let genegenericrange = gene.start.getzbasedpos()..gene.end.getobasedpos();
        //As gene start is 1-ranged, put it as 0-range with -1. End is exclusive so -1/+1 = 0
        reader
            .fetch((
                &gene.chromosome,
                genegenericrange.start,
                genegenericrange.end,
            ))
            .map_err(|e| io::Error::new(io::ErrorKind::BrokenPipe, e))?;
        let mut hash: BTreeMap<Position, Posread> = BTreeMap::new(); //Match and full match and total
        //Hash contains 1-based positions
        genegenericrange.for_each(|p| {
            hash.insert(
                Position::newfromzposition(p),
                //Default should not trigger as no error possible
                Posread::new(0, 0, 0, 0, 0f32, args)
                    .unwrap_or_else(|_| unreachable!("Error on Posread")),
            );
        });
        (hash, reader.records())
    };
    let mut coverageperc = 0;
    let mut empty = true;
    for mut record in records
        .filter_map(Result::ok)
        .filter(|p| filterread(args, p))
    {
        empty = false;
        reads += 1;
        record.cache_cigar();
        if let Some(d) = hash.get_mut(&Position::newfromzposition(record.reference_start()))
            && (record.cigar_cached().unwrap().leading_softclips() > 0 || record.cigar_cached().unwrap().leading_softclips() > 0)
        {
            d.addsoftclip(1);
        } else if let Some(d) = hash.get_mut(&Position::newfromzposition(record.reference_end()))
            && (record.cigar_cached().unwrap().leading_softclips() > 0 || record.cigar_cached().unwrap().leading_hardclips() > 0)
        {
            d.addsoftclip(1);
        }
        let range = record.reference_start()..record.reference_end();
        coverageperc += ranges::Ranges::from(range.clone()).into_iter().count();
        let mut oldcigar = None;
        'outer: for (_, refpos, cigar) in record.aligned_pairs_full() {
            match (cigar, oldcigar) {
                (Cigar::Ins(_), Some(old)) => {
                    match hash.get_mut(&Position::newfromzposition(old)) {
                        Some(d) => d.addinsertion(1),
                        None => {
                            if old > gene.end.getzbasedpos() {
                                break 'outer;
                            }
                        } //Outside coverage of gene
                    }
                }
                _ => {
                    oldcigar = refpos;
                }
            }
        }
        'outer: for [start, end] in record.aligned_blocks() {
            for p in start..end {
                match hash.get_mut(&Position::newfromzposition(p)) {
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
                    match hash.get_mut(&Position::newfromzposition(p)) {
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
            let range = p[0]..p[1];
            range.contains(&genestart.getzbasedpos()) && range.contains(&geneend.getzbasedpos())
        };
        if record
            .aligned_blocks()
            .any(|p| validrange(p, &gene.start, &gene.end))
        {
            reads100 += 1;
        }
        if validrange([range.start, range.end], &gene.start, &gene.end) {
            readsfull += 1;
        }
        if !args.force
            && iterblock(&record)
                .is_some_and(|f| f.into_iter().any(|p| validrange(p, &gene.start, &gene.end)))
        {
            let mut aligned = record.aligned_pairs();
            let realstart = aligned.find(|p| p.last() == Some(&gene.start.getobasedpos()));
            let realend = aligned.find(|p| p.last() == Some(&gene.end.getobasedpos()));
            if let (Some(Some(start)), Some(Some(end))) = (
                realstart.map(|p| p.first().copied()),
                realend.map(|p| p.first().copied()),
            ) {
                let newphredscore = record
                    .qual()
                    .iter()
                    .skip(start.try_into().unwrap_or_default())
                    .take(
                        end.saturating_sub(start)
                            .saturating_sub(1)
                            .try_into()
                            .unwrap_or_default(),
                    )
                    .copied()
                    .collect::<Vec<u8>>();
                let minscore = newphredscore
                    .iter()
                    .filter(|a| **a > 0u8)
                    .nth(newphredscore.len() / 2);
                //Ponderation based on lowest quality of a read at a gene position
                let fscore = match minscore {
                    Some(0..=10) | None => 0,
                    Some(11..=20) => 1,
                    Some(21..=30) => 3,
                    Some(31..=40) => 7,
                    Some(41..=50) => 9,
                    Some(51..) => 10,
                } as f32;
                phredscore.push(newphredscore);
                realreads100m = realreads100m.mul(10f32).add(fscore).div(10f32);
            }
            reads100m += 1;
        }
        for p in range {
            match hash.get_mut(&Position::newfromzposition(p)) {
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
        let elem = GeneInfosFinish::make_default(gene.clone());
        return Ok((elem, hash));
    }
    hash.iter_mut().for_each(|(_, p)| {
        if p.softclips.is_normal() {
            p.softclips = (p.softclips * 100f32 / max(p.gettotal(), 1) as f32).round() / 100f32
        }
    });
    //Coverage calculus
    let coverage = hash
        .iter()
        .filter(|(_, p)| p.gettotal() >= args.coverage.get().try_into().unwrap_or(usize::MAX))
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
                "{}({}={}X{}D{}I{}S)-",
                f.gettotal(),
                f.getmatch(),
                f.getmismatchcount(),
                f.getindelcount(),
                f.getinsertion(),
                f.getsoftclip()
            ));
            acc
        })
    };
    let text = String::from(text.trim_end_matches('-'));
    gene.setstatus(reads100m, &hash);
    let coverageperc = ((coverageperc * 1_000
        / reads
        / usize::try_from(gene.end.length(&gene.start)).unwrap_or(usize::MAX))
        as f32)
        .round()
        / 1_000.0;
    let elem = GeneInfosFinish::new(
        gene.clone(),
        reads,
        readsfull,
        Some(text),
        reads100,
        reads100m,
        realreads100m,
        coverageperc,
        phredscore,
        coverage,
    );
    Ok((elem, hash))
}
fn genelist(
    loci: &LocusInfos,
    species: &Species,
    args: &Args,
    full: bool,
) -> Result<Vec<GeneInfosFinish>, Box<dyn std::error::Error>> {
    let genes = extractgenelist(args, loci, full)?;
    if genes.is_empty() {
        return Ok(Vec::new());
    }
    let outputdir = &args.outdir;
    let outputfile = outputdir.join(givename(
        species,
        loci.getlocus(),
        &loci.contig,
        loci.gethaplotype().isprimary(),
        "geneanalysis.csv",
        false,
    ));
    let mut finale: Vec<GeneInfosFinish> = Vec::with_capacity(genes.len());
    //For each gene, list of alerting positions, bbool said suspicious or warning position
    let mut alertingpositions: BTreeMap<GeneInfos, Vec<(bool, usize)>> = BTreeMap::new();
    let progressbar = getprogressbarclassic(genes.len().try_into().unwrap_or_default())?;
    for (pos, mut gene) in genes.into_iter().enumerate() {
        let (elem, hash) = generategeneinfos(args, &mut gene)?;
        let plots = outputdir.join(format!(
            "gene_{}",
            loci.gethaplotype().to_string().as_str().to_lowercase()
        ));
        if !std::fs::exists(&plots)? {
            println!("Creating the folder {}", plots.display());
            std::fs::create_dir_all(&plots)?;
        };
        let mut output = plots.join(
            regexpword
                .replace_all(&gene.gene.to_string(), "_")
                .to_uppercase(),
        );
        let size = (700, 700);
        if !args.svg {
            output.set_extension("png");
            let root = BitMapBackend::new(&output, size).into_drawing_area();
            //Gene graph
            genegraph(args, &hash, loci, root, &mut alertingpositions, &elem)?;
        } else {
            output.set_extension("svg");
            let root = SVGBackend::new(&output, size).into_drawing_area();
            //Gene graph
            genegraph(args, &hash, loci, root, &mut alertingpositions, &elem)?;
        }
        finale.push(elem);
        progressbar.set_position(pos.saturating_add(1).try_into().unwrap_or_default());
    }
    progressbar.finish();
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
        printpossus(args, species, loci, outputdir, &alertingpositions)?;
    }
    println!("Gene analysis has been saved to {}", outputfile.display());
    Ok(finale)
}
fn printbreaks(
    species: &Species,
    finalpos: i64,
    breaks: &[(i64, i64)],
    loci: &LocusInfos,
    outputdir: &std::path::Path,
) -> std::io::Result<()> {
    let mut breakfile = File::create(outputdir.join(givename(
        species,
        loci.getlocus(),
        &loci.contig,
        loci.gethaplotype().isprimary(),
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
    _args: &Args,
    species: &Species,
    loci: &LocusInfos,
    outputdir: &std::path::Path,
    data: &BTreeMap<GeneInfos, Vec<(bool, usize)>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let outputfile = outputdir.join(givename(
        species,
        loci.getlocus(),
        &loci.contig,
        loci.gethaplotype().isprimary(),
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
        csv.write_field(gene.gene.to_string())?;
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
    loci: &LocusInfos,
    root: DrawingArea<T, Shift>,
    alerting: &mut BTreeMap<GeneInfos, Vec<(bool, usize)>>,
    gene: &GeneInfosFinish,
) -> Result<(), Box<dyn std::error::Error>>
where
    T: DrawingBackend,
{
    let genename = gene.gene.to_string();
    let text_style = fontstyle.into_text_style(&root);
    let _ = root.fill(&plotters::prelude::WHITE);
    let max = hash.values().map(|p| p.gettotal()).max().unwrap_or(0) + 5;
    let colorgene = if gene.status.getstatus().isvalid() {
        full_palette::GREEN
    } else {
        full_palette::RED
    };
    let (top, bottom) = root.split_vertically((70).percent_height());
    let mut chart = ChartBuilder::on(&top)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            format!(
                "Reads coverage for {} ({}-{})",
                genename,
                loci.gethaplotype(),
                gene.strand
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
    chart
        .draw_series(LineSeries::new(
            hash.iter()
                .enumerate()
                .map(|(pos, (_, val))| (pos + 1, val.getinsertion())),
            full_palette::BLUEGREY_600.mix(0.7).stroke_width(3),
        ))
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Insertion")
        .legend(|(x, y)| {
            PathElement::new(
                vec![(x, y), (x + 15, y)],
                full_palette::BLUEGREY_600.mix(0.8),
            )
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
                    .map(|(pos, ..)| (pos + 1, gene.reads100m)),
                BLACK.mix(0.7).stroke_width(3),
            ))
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
            .label("Reads 100% match")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLACK.mix(0.7)));
        //Line of all full reads
        chart
            .draw_series(LineSeries::new(
                hash.iter()
                    .enumerate()
                    .map(|(pos, _)| (pos + 1, gene.realreads100m.round() as usize)),
                BROWN_500.mix(0.7).stroke_width(3),
            ))
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
            .label("Reads 100% match (recalculated)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BROWN_500.mix(0.7)));
        chart
            .draw_series(
                Histogram::vertical(&chart)
                    .style(full_palette::ORANGE_300.mix(0.3).filled())
                    .data(hash.iter().enumerate().filter_map(|(pos, (_, val))| {
                        let pos1 = pos + 1;
                        if val.iswarning() {
                            match alerting.get_mut(&GeneInfos::from(gene.clone())) {
                                Some(d) => d.push((false, pos1)),
                                None => {
                                    alerting
                                        .insert(GeneInfos::from(gene.clone()), vec![(false, pos1)]);
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
                            match alerting.get_mut(&GeneInfos::from(gene.clone())) {
                                Some(d) => d.push((true, pos1)),
                                None => {
                                    alerting
                                        .insert(GeneInfos::from(gene.clone()), vec![(true, pos1)]);
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
        0..max,
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
                .data(hash.iter().enumerate().filter_map(|(index, (_, p))| {
                    if p.softclips.is_normal() {
                        Some((index + 1, (p.softclips * 100f32) as usize))
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
            .position(plotters::chart::SeriesLabelPosition::LowerLeft)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    }
    let max = 100;
    let mut chart = ChartBuilder::on(&bottom)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .right_y_label_area_size(60)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Average PHRED score for matching reads", ("sans-serif", 18))
        .build_cartesian_2d(1..hash.len(), 0..max)
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    chart
        .draw_series(LineSeries::new(
            hash.iter().enumerate().filter_map(|(index, _)| {
                let iter = gene.phredscore.getscore();
                let infos: Vec<&u8> = iter.iter().filter_map(|b| b.get(index)).collect();
                infos
                    .get(infos.len() / 2)
                    .map(|b| (index + 1, min(100, (**b).into())))
            }),
            full_palette::BLUEGREY.mix(0.4),
        ))
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("PHRED score")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLACK.mix(0.7)));
    chart
        .draw_series(
            Histogram::vertical(&secondary)
                .style(full_palette::ORANGE_300.mix(0.3).filled())
                .data(hash.iter().enumerate().filter_map(|(pos, (_, val))| {
                    let pos1 = pos + 1;
                    if val.iswarning() {
                        Some((pos1, max))
                    } else {
                        None
                    }
                }))
                .margin(0),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    chart
        .draw_series(
            Histogram::vertical(&secondary)
                .style(full_palette::ORANGE_300.mix(0.3).filled())
                .data(hash.iter().enumerate().filter_map(|(pos, (_, val))| {
                    let pos1 = pos + 1;
                    if val.issuspicious() {
                        Some((pos1, max))
                    } else {
                        None
                    }
                }))
                .margin(0),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    chart
        .configure_mesh()
        //.x_label_formatter(&|f| f.to_formatted_string(&Locale::en).to_string())
        //.x_desc("Position in sequence (bp)")
        .y_desc("PHRED score")
        .disable_x_mesh()
        .label_style(text_style.clone())
        .light_line_style(GREY_400.mix(0.6))
        .x_max_light_lines(5)
        .y_max_light_lines(2)
        //.disable_y_mesh()
        .draw()
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    if !args.nolegend {
        chart
            .configure_series_labels()
            .position(plotters::chart::SeriesLabelPosition::LowerRight)
            .background_style(WHITE.mix(0.6))
            .label_font(text_style.clone())
            .border_style(BLACK.mix(0.8))
            .draw()
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?;
    }
    // To avoid the IO failure being ignored silently, we manually call the present function
    drawnoticetext(&bottom)?;
    bottom.present().map_err(|_e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir")))?;
    top.present().map_err(|_e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, "Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir")))?;
    Ok(())
}
fn givename(
    species: &Species,
    locus: &Locus,
    contig: &str,
    haplo: bool,
    suffix: &str,
    image: bool,
) -> String {
    format!(
        "{}_{}_{}{}_{}",
        species.safestring(),
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
    species: &Species,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    _args: &Args,
) -> Result<(), Box<dyn std::error::Error>> {
    let outputfile = outputdir.join(givename(
        species,
        loci.getlocus(),
        &loci.contig,
        loci.gethaplotype().isprimary(),
        "positionresult.csv",
        false,
    ));
    let outputfile = outputfile.as_path();
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .has_headers(true)
        .delimiter(b'\t')
        .flexible(false)
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
                loci.getlocus(),
                loci.contig,
                loci.gethaplotype()
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
                if let Some(a) = p.qual
                    && a > 0
                {
                    Some((p.position.getobasedpos(), a))
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
                    loci.getlocus(), loci.contig
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
pub(crate) fn checknewversion() -> bool {
    let r = match REQUESTCLIENT
        .get(GITHUBVERSION.to_string())
        .send()
        .map(|a| a.json::<Vec<GitLabTag>>())
    {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error checking new version: {e}");
            return false;
        }
    };
    let r = if let Ok(r) = r {
        r
    } else {
        return false;
    };
    if let Some(a) = r.first() {
        let version = &a.name;
        let regexp = if let Ok(a) = Regex::new(r"([0-9]+)\.([0-9]+)\.([0-9]+)") {
            a
        } else {
            return false;
        };
        let cap = regexp.captures(version);
        let cap2 = regexp.captures(VERSION);
        let (major, minor, patch, major2, minor2, patch2) = match (cap, cap2) {
            (Some(b), Some(c))
                if let Some(maj) = b.get(1)
                    && let Some(min) = b.get(2)
                    && let Some(patch) = b.get(3)
                    && let Some(maj2) = c.get(1)
                    && let Some(min2) = c.get(2)
                    && let Some(patch2) = c.get(3) =>
            {
                (
                    maj.as_str(),
                    min.as_str(),
                    patch.as_str(),
                    maj2.as_str(),
                    min2.as_str(),
                    patch2.as_str(),
                )
            }
            _ => {
                eprintln!("Invalid version");
                return false;
            }
        };
        let (major, minor, patch, major2, minor2, patch2) = match (
            major.parse::<usize>(),
            minor.parse::<usize>(),
            patch.parse::<usize>(),
            major2.parse::<usize>(),
            minor2.parse::<usize>(),
            patch2.parse::<usize>(),
        ) {
            (Ok(a), Ok(b), Ok(c), Ok(d), Ok(e), Ok(f)) => (a, b, c, d, e, f),
            _ => {
                eprintln!("Invalid version: {}", version);
                return false;
            }
        };
        let name = env!("CARGO_CRATE_NAME").replace("_", "/");
        if major > major2 || (major == major2 && minor > minor2) {
            eprintln!("Another version of {name} is available. Please update.");
            true
        } else if major == major2 && minor == minor2 && patch > patch2 {
            eprintln!("A patch version of {name} is available. Please update.");
            false
        } else if major2 > major || (major == major2 && minor2 > minor) {
            let now = chrono::Local::now();
            let limit = chrono::DateTime::parse_from_str(LIMITDATE.as_str(), "%+");
            if let Ok(date) = limit
                && date.timestamp() < now.timestamp()
            {
                eprintln!(
                    "You have a test version of {name}. Please use the public version or contact {AUTHOR} to reactive the test version."
                );
                true
            } else {
                false
            }
        } else {
            false
        }
    } else {
        false
    }
}
pub(crate) fn geneisokay(reads100m: usize, hash: &BTreeMap<Position, Posread>) -> OkStatus {
    if let Some(a) = hash
        .iter()
        .find(|(_, f)| !f.isvalid() || f.softclips >= SOFTCLIPRATIO)
    {
        match a {
            (f, a) if !a.isvalid() => OkStatus::new(
                AcceptedStatus::Rejected,
                Some(format!(
                    "{} {} is a {}",
                    *SUSPICIOUSPOSITIONALERT,
                    f.getobasedpos(),
                    if a.iswarning() {
                        "warning"
                    } else {
                        "suspicious"
                    }
                )),
            ),
            (f, _) => OkStatus::new(
                AcceptedStatus::Rejected,
                Some(format!(
                    "{} at position {}",
                    *SOFTCLIPTOOMUCH,
                    f.getobasedpos()
                )),
            ),
        }
    } else if reads100m < MATCHREADS.try_into().unwrap_or_default() {
        OkStatus::new(
            AcceptedStatus::Rejected,
            Some(NOTENOUGHMATCHREADS.to_string()),
        )
    } else {
        OkStatus::new(AcceptedStatus::Accepted, None)
    }
}
pub(crate) fn telomereposition(p: &Position, lengthloci: &Position) -> bool {
    lengthloci.getobasedpos().saturating_sub(p.getobasedpos())
        < TELOMERESEP.try_into().unwrap_or(i64::MAX)
        || p.getobasedpos() < TELOMERESEP.try_into().unwrap_or(i64::MAX)
}
pub(crate) fn locusisokay(mean: u64, lengthloci: &Position, graph: &[&HashMapinfo]) -> OkStatus {
    //Between a minimum and a maximum number of reads
    if let Some(a) = graph
        .iter()
        .filter(|d| !telomereposition(&d.position, lengthloci))
        .find(|f| !coveragewindow(mean).contains(&f.overlaps.try_into().unwrap_or_default()))
    {
        match a {
            a if a.overlaps < usize::try_from(*coveragewindow(mean).start()).unwrap_or_default() => OkStatus::new(
                Rejected,
                Some(format!(
                    "{} at position {} ({})",
                    *INVALIDCOVERAGE,
                    a.position.getobasedpos(),
                    "too low",
                )),
            ),
            _ => OkStatus::new(
                Rejected,
                Some(format!(
                    "{} at position {} ({})",
                    *INVALIDCOVERAGE,
                    a.position.getobasedpos(),
                    "too high",
                )),
            ),
        }
    } else if let Some(a) = graph.iter().find(|f| f.psoftclips > SOFTCLIPRATIO) {
        OkStatus::new(
            Rejected,
            Some(format!(
                "{} at position {}",
                *SOFTCLIPTOOMUCH,
                a.position.getobasedpos()
            )),
        )
    } else {
        OkStatus::new(AcceptedStatus::Accepted, None)
    }
}
pub(crate) fn coveragewindow(mean: u64) -> RangeInclusive<u64> {
    max(
        MINIMUMCOVERAGE.try_into().unwrap_or_default(),
        (mean as f32 / MAXCOVERAGERATIO).round() as u64,
    )..=((mean as f32 * MAXCOVERAGERATIO).round() as u64)
}
fn readgraph<T>(
    outputfile: &std::path::Path,
    loci: &LocusInfos,
    pos: &[&HashMapinfo],
    args: &Args,
    species: &Species,
    root: DrawingArea<T, Shift>,
    mean: u64,
) -> Result<(), Box<dyn std::error::Error>>
where
    T: DrawingBackend,
{
    let text_style = fontstyle.into_text_style(&root);
    let max = max(
        usize::try_from(mean).unwrap_or(usize::MAX),
        pos.iter().map(|max| max.getmaxvalue()).max().unwrap_or(0).try_into().unwrap_or_default(),
    ) + 5;
    let _ = root.fill(&plotters::prelude::WHITE);
    let (top, bottom) = root.split_vertically((80).percent_height());
    let (start, end, tenlines) = if let Some(full) = loci.getfulllength(args)
        && telomereposition(&loci.start, &full)
    {
        (
            None,
            Some(full),
            loci.getlength().saturating_sub(full.length(&loci.end)) / 9,
        )
    } else if let Some(full) = loci.getfulllength(args)
        && (telomereposition(&loci.end, &full))
    {
        (
            Some(
                TELOMERESEP
                    .saturating_sub(loci.start.getobasedpos().try_into().unwrap_or_default()),
            ),
            Some(full),
            loci.getlength().saturating_sub(full.length(&loci.end)) / 9,
        )
    } else {
        println!("Length for {} and {} is {:?}",loci.getlocushaplo(),loci.contig,loci.getfulllength(args));
        (None, None, loci.getlength() / 9)
    };
    let colorgene = if loci.status.getstatus().isvalid() {
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
                loci.getlocus(),
                loci.contig,
                loci.gethaplotype()
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
                        || p.position == loci.end
                    {
                        Some((
                            p.position.getobasedpos(),
                            mean.try_into().unwrap_or(usize::MAX),
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
    let coveragewindows = coveragewindow(mean);
    chart
        .draw_series(
            AreaSeries::new(
                pos.iter().filter_map(|p| {
                    if (p.position.getobasedpos() - loci.start.getobasedpos())
                        .rem_euclid(std::cmp::max(tenlines, 1))
                        == 0
                        || p.position == loci.end
                    {
                        Some((p.position.getobasedpos(), usize::try_from(*coveragewindows.end()).unwrap_or_default()))
                    } else {
                        None
                    }
                }),
                usize::try_from(*coveragewindows.start()).unwrap_or_default(),
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
                loci.getlocus(), loci.contig
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
            if elem.overlaps <= usize::try_from(args.breaks).unwrap_or_default() {
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
    printbreaks(species, finalpos, &breaks, loci, outputfile)?;
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
    if start.is_some() || end.is_some() {
        chart
            .draw_series(AreaSeries::new(
                pos.iter().filter_map(|elem| {
                    if telomereposition(
                        &elem.position,
                        &end.unwrap_or(Position::new(true, i64::MAX)),
                    ) {
                        Some((SegmentValue::CenterOf(elem.position.getobasedpos()), 100))
                    } else {
                        None
                    }
                }),
                0,
                full_palette::BLUEGREY_A200.mix(0.8).filled(),
            ))
            .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
            .label("Telomere region")
            .legend(|(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], full_palette::BLUEGREY_A200)
            });
    }
    chart
        .draw_series(
            Histogram::vertical(&chart)
                .baseline(0)
                .margin(3)
                .data(pos.iter().filter_map(|p| {
                    if p.psoftclips.is_normal() {
                        Some((p.position.getobasedpos(), (p.psoftclips * 100f32) as i64))
                    } else {
                        None
                    }
                }))
                .style(full_palette::BLACK.mix(0.7).filled()),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Not primary softclips percent")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));
        chart
        .draw_series(
            Histogram::vertical(&chart)
                .baseline(0)
                .margin(3)
                .data(pos.iter().filter_map(|p| {
                    if p.psoftclips.is_normal() {
                        Some((p.position.getobasedpos(), (p.psoftclips * 100f32) as i64))
                    } else {
                        None
                    }
                }))
                .style(full_palette::DEEPPURPLE_600.mix(0.8).filled()),
        )
        .map_err(|e| Box::new(io::Error::new(io::ErrorKind::InvalidInput, e.to_string())))?
        .label("Primary softclips percent")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], full_palette::DEEPPURPLE_600));
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

#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
use crate::identification::{downloadmotifs, downloadref, retainbestmatch};
use crate::r#struct::{
    Args, Blast, Blastcalc, Blastlevel, Blastmatch, GeneInfos, Haplotype, Locus, LocusInfos, Name,
    Newfasta, Position, Seqresult, Species, Strand,
};
use crate::{getassemblyreader, getreaderoffile};
use bio::io::fasta;
use extended_htslib::bam::{self, Read};
use flate2::{Compression, write::GzEncoder};
use itertools::Itertools;
use lazy_static::lazy_static;
use reqwest::{StatusCode, tls};
use serde_json::{self as json, Value};
use std::cmp::{Ordering, max, min};
use std::collections::HashMap;
use std::fmt::Debug;
use std::io::IsTerminal;
use std::str::FromStr;
use std::{
    env::{self, temp_dir},
    error::Error,
    fs::{self, File},
    io::{self, BufRead, BufReader, ErrorKind, Read as _},
    ops::{Not, RangeInclusive},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::Duration,
};
use strum::IntoEnumIterator;
use tempfile::{NamedTempFile, TempDir};
lazy_static! {
    pub static ref REQUESTCLIENT: reqwest::blocking::Client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::new(15, 0))
        .timeout(Duration::new(300, 0))
        .tcp_user_timeout(Duration::new(400, 0))
        .user_agent(format!("IMGT/StatAssembly version {}", VERSION))
        .referer(false)
        .tls_version_min(tls::Version::TLS_1_2)
        .https_only(true)
        .http3_max_idle_timeout(Duration::new(10,0))
        .http3_prior_knowledge()
        .build()
        .unwrap_or_default();
    pub static ref WEBSERVER: String = obfstr::obfstring!("https://imgt.org");
    pub static ref RELEASELINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstr!("/download/GENE-DB/RELEASE")
    );
    pub static ref MOTIFLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!("/submissions/newmotif_fusionne.fasta.gz")
    );
    /* pub static ref MOTIFLINK: String =
        obfstr::obfstring!("http://localhost:8910/submissions/newmotif_fusionne.fasta.gz");
    pub static ref BORNESLINK: String =
        obfstr::obfstring!("http://localhost:8910/submissions/bornes_mammals.fasta.gz");  */
        pub static ref BORNESLINK: String = format!(
            "{}{}",
            WEBSERVER.as_str(),
            obfstr::obfstring!("/submissions/bornes_mammals.fasta.gz")
        );
    pub static ref VQUESTLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!(
            "/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"
        )
    );
    pub static ref SUBMISSIONLINK: String = format!(
        "{}{}",
        WEBSERVER.as_str(),
        obfstr::obfstring!("/submissions/")
    );
}
pub(crate) const VERSION: &str = env!("CARGO_PKG_VERSION");
pub(crate) const DELIMITERFASTA: char = '/';
pub(crate) const LOCUSSEPARATOR: usize = 1_000_000;

pub(crate) fn readfromterminal(yes: &char, no: &char, force: bool) -> bool {
    let io = io::stdin();
    if !io.is_terminal() {
        return false;
    }
    let mut data = String::new();
    loop {
        data.clear();
        if io.read_line(&mut data).is_err() {
            return false;
        }
        if data.trim().is_empty() && !force {
            return false;
        }
        if data.to_lowercase().trim().chars().all(|p| &p == yes) {
            return true;
        } else if data.to_lowercase().trim().chars().all(|p| &p == no) || !force {
            return false;
        }
    }
}
#[deprecated]
pub(crate) fn getnamefromblast<T>(text: T) -> Option<String>
where
    T: AsRef<str>,
{
    text.as_ref()
        .to_ascii_lowercase()
        .split('|')
        .nth(2)
        .map(|f| f.to_string())
}
pub(crate) fn getallelefromblast(text: &Name) -> &str {
    &text.gene
}
pub(crate) fn getchromosomefromblast<T>(text: T) -> Option<String>
where
    T: AsRef<str>,
{
    text.as_ref()
        .to_ascii_uppercase()
        .split('|')
        .nth(3)
        .map(|f| f.to_string())
}
pub(crate) fn getpositionfromblast(text: &str) -> Option<(Position, Position, Strand)> {
    text.to_ascii_uppercase().split("|").nth(4).and_then(|p| {
        let p = p.replace("/", "..");
        let mut split = p.splitn(3, "..");
        match (
            split.next().map(|p| p.parse::<i64>()),
            split.next().map(|p| p.parse::<i64>()),
            split.next().map(|a| a.trim()),
        ) {
            (Some(Ok(a)), Some(Ok(b)), c) => {
                let complement = c.is_some_and(|f| !f.is_empty());
                Some((
                    Position::new(false, a),
                    Position::new(false, b),
                    if complement {
                        Strand::Minus
                    } else {
                        Strand::Plus
                    },
                ))
            }
            _ => None,
        }
    })
}
#[must_use]
pub(crate) fn fastafilter(text: &str, find: &str, present: bool, species: bool) -> String {
    let buffer = BufReader::new(text.as_bytes());
    let parser = fasta::Reader::new(buffer);
    let mut lines: Vec<String> = Vec::with_capacity(text.lines().count());
    for record in parser
        .records()
        .filter_map(Result::ok)
        .filter(|p| p.check().is_ok())
    {
        let val = format!(
            "{}{}",
            record.id(),
            record.desc().map_or(String::new(), |a| format!(" {a}"))
        );
        let keep = match Name::from_str(&val) {
            Ok(result) => match (present, species) {
                (true, true) => result.species.eq_ignore_ascii_case(find),
                (false, true) => !result.species.eq_ignore_ascii_case(find),
                (false, false) => !result.gene.eq_ignore_ascii_case(find),
                (true, false) => result.gene.eq_ignore_ascii_case(find),
            },
            Err(_e) => {
                //eprintln!("Error parsing line: {}. Error is {e}", val);
                false
            }
        };
        if keep {
            lines.push(record.to_string());
        }
    }
    lines
        .iter()
        .fold(String::new(), |mut acc, f| {
            acc.push_str(format!("\n{}", f.trim()).as_str());
            acc
        })
        .trim()
        .to_string()
}
pub(crate) fn checkoverlap(a: &RangeInclusive<usize>, b: &RangeInclusive<usize>) -> bool {
    max(a.start(), b.start()) <= min(a.end(), b.end())
}
pub(crate) fn fileincache<T>(tempfile: T) -> bool
where
    T: AsRef<Path>,
{
    if !tempfile.as_ref().is_file() {
        return false;
    }
    tempfile
        .as_ref()
        .metadata()
        .and_then(|p| {
            p.accessed()
                .map(|p| p.elapsed().map(|p| p < Duration::new(3600 * 24 * 7, 0)))
        })
        .is_ok_and(|f| f.unwrap_or_default())
}
pub(crate) fn checkifblastpresent() -> bool {
    println!("Detect if BLAST is operating.");
    let command = Command::new("blastn")
        .arg("-version")
        .stdout(Stdio::null())
        .stdin(Stdio::null())
        .status()
        .is_ok_and(|f| f.success());
    if !command {
        eprintln!("BLAST was not found or has failed. Check if present in PATH.");
        false
    } else {
        println!("BLAST is working. Continuing.");
        true
    }
}
pub(crate) fn positionfiltering<T>(locus: &[LocusInfos], blast: &mut Vec<T>)
where
    T: Blastcalc,
{
    blast.retain(|p| {
        locus.iter().any(|f| {
            let range: RangeInclusive<usize> = f.start.getobasedpos().try_into().unwrap_or_default()
                ..=f.end.getobasedpos().try_into().unwrap_or_default();
            f.contig.as_str() == p.getsubject()
                && checkoverlap(&(p.getpos().0..=p.getpos().1), &range)
        })
    })
}
pub(crate) fn locusfiltering<T>(locus: &Locus, blast: &mut Vec<T>)
where
    T: Blastcalc,
{
    //TRD is inside TRA locus
    if locus == &Locus::TRA {
        blast.retain(|p| {
            p.getlocusname().as_ref() == Some(locus) || p.getallelename().starts_with("TRD")
        });
    } else {
        blast.retain(|p| p.getlocusname().as_ref() == Some(locus));
    }
}
pub(crate) fn checklocus(p: &fasta::Record, locus: &Locus) -> bool {
    match locus {
        Locus::IGH => p.id().contains("IGH") || p.desc().is_some_and(|f| f.contains("IGH")),
        Locus::IGK => p.id().contains("IGK") || p.desc().is_some_and(|f| f.contains("IGK")),
        Locus::IGL => p.id().contains("IGL") || p.desc().is_some_and(|f| f.contains("IGL")),
        Locus::TRA => {
            p.id().contains("TRA")
                || p.desc().is_some_and(|f| f.contains("TRA"))
                || p.id().contains("TRD")
                || p.desc().is_some_and(|f| f.contains("TRD"))
        }
        Locus::TRB => p.id().contains("TRB") || p.desc().is_some_and(|f| f.contains("TRB")),
        Locus::TRG => p.id().contains("TRG") || p.desc().is_some_and(|f| f.contains("TRG")),
    }
}
pub(crate) fn speciesandorphonfiltering(
    tempfile: &Path,
    locus: Option<&Locus>,
    releaseversion: String,
    species: &str,
    orphonfilter: bool,
    force: bool,
) -> io::Result<PathBuf> {
    let outfile = if force {
        NamedTempFile::with_prefix_in(
            format!(
                "refseq{}-{}{}.fasta",
                releaseversion.replace(" ", "_"),
                species.replace(" ", "-"),
                locus.map_or("".to_string(), |l| format!("-{}", l))
            ),
            env::temp_dir(),
        )?
        .into_temp_path()
        .to_path_buf()
    } else {
        Path::join(
            &env::temp_dir(),
            format!(
                "refseq{}-{}{}.fasta",
                releaseversion.replace(" ", "_"),
                species.replace(" ", "-"),
                locus.map_or("".to_string(), |l| format!("-{}", l))
            ),
        )
    };
    if outfile.try_exists().unwrap_or(false) && !force {
        println!("Filtering already done, retrieving...");
        return Ok(outfile);
    }
    println!("Filtering based on species {}.", species);
    let file = std::fs::read_to_string(tempfile)?;
    let info = fastafilter(&file, species, true, true).replace(" ", "_");
    let info = if orphonfilter {
        println!("Orphon filtering");
        fastafilter(&info, "/OR", false, false)
    } else {
        info
    };
    let tempnew = NamedTempFile::new_in(temp_dir())?;
    fs::write(&tempnew, &info)?;
    let read = fasta::Reader::from_file(&tempnew)
        .map_err(|f| io::Error::new(ErrorKind::InvalidInput, f))?;
    let allmatch = if let Some(l) = locus {
        read.records()
            .filter_map(Result::ok)
            .any(|p| checklocus(&p, l))
    } else {
        let mut count = Locus::iter().map(|p| (p, 0)).collect_vec();
        read.records().filter_map(Result::ok).for_each(|p| {
            for (locus, count) in count.iter_mut() {
                if checklocus(&p, locus) {
                    *count += 1;
                    break;
                }
            }
        });
        count.iter().all(|(_, f)| *f > 0)
    };
    let finale = match (info.is_empty(), allmatch) {
        (true, _) => {
            println!(
                "No match with IMGT/GENE-DB, the species {} might be a new species.",
                species
            );
            file
        }
        (false, false) => {
            if let Some(a) = locus {
                println!(
                    "No match with IMGT/GENE-DB with {} and {} loci.",
                    species, a
                );
            } else {
                println!("No match with IMGT/GENE-DB with {} and some loci.", species);
            }
            file
        }
        _ => info,
    };
    let info = if let Some(l) = locus {
        let read = fasta::Reader::from_file(&tempnew)
            .map_err(|f| io::Error::new(ErrorKind::InvalidInput, f))?;
        let records = read
            .records()
            .filter_map(Result::ok)
            .filter(|p| checklocus(p, l));
        let seq = records.fold(String::new(), |mut acc, f| {
            acc.push_str(&format!(
                ">{} {}\n{}\n",
                f.id(),
                f.desc().unwrap_or_default(),
                String::from_utf8_lossy(f.seq())
            ));
            acc
        });
        seq.trim().to_string()
    } else {
        finale
    };
    if info.is_empty() {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "Empty data after filtering",
        ));
    }
    /* let newdata = if newdata.is_empty() {
        println!("New species!!");
        let bufreader = io::Cursor::new(newdata);
        let reader = fasta::Reader::new(bufreader);
        let fastaread = readfastareader(reader);
        return match fastaread {
            Ok(b) => {
                let mut new = Vec::new();
                b.iter().for_each(|f| {
                    if let Ok(v) = Newfasta::new(f) {
                        new.push(v);
                    } else {
                        eprintln!("Line {} is skipped as an error occured", f.name)
                    }
                });
                (Some(new), None)
            }
            Err(e) => {
                eprintln!("{e}");
                (None, None)
            }
        };
        file
    } else {
        newdata
    }; */
    if tempfile != outfile {
        std::fs::write(&outfile, &info)?;
    }
    Ok(outfile)
}
pub(crate) fn readfastareader<T>(fasta: fasta::Reader<T>) -> io::Result<Vec<Newfasta>>
where
    T: std::io::BufRead,
{
    let mut seqs = Vec::new();
    for record in fasta.records().filter_map(|p| p.ok()) {
        if let Err(v) = record.check() {
            eprintln!("Seq {} is invalid, skipped. Error is {v}", record.id());
            continue;
        }
        let text = format!(
            "{} {}\n{}",
            record.id(),
            record.desc().unwrap_or_default(),
            String::from_utf8_lossy(record.seq()).to_string()
        );
        match Newfasta::from_str(&text) {
            Ok(a) => seqs.push(a),
            Err(e) => {
                eprintln!("Seq {} is invalid, skipped. Error is {}", record.id(), e);
                continue;
            }
        };
    }
    if seqs.is_empty() {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "No fasta record valid found",
        ));
    }
    Ok(seqs)
}
#[allow(unused)]
pub(crate) fn readfastafile(seq: &Path) -> io::Result<Vec<Newfasta>> {
    let file = File::open(seq)?;
    let fasta = fasta::Reader::new(file);
    readfastareader(fasta)
}
pub(crate) fn statusblastmotifs(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.sseqid.cmp(&b.sseqid) {
        std::cmp::Ordering::Equal => match a.sstart.cmp(&b.sstart) {
            Ordering::Equal => a.send.cmp(&b.send).reverse(),
            ord => ord,
        },
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && g.getallelename() == f.getallelename()
                    && checkoverlap(&(g.sstart..=g.send), &(f.sstart..=f.send))
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn statusblastvs(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.getsubject().cmp(&b.getsubject()) {
        std::cmp::Ordering::Equal => a.getquery().gene.cmp(&b.getquery().gene),
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && g.sseqid == f.sseqid
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn statusblast(data: &mut Vec<Blast>) {
    data.sort_unstable_by(|a, b| match a.sstart.cmp(&b.sstart) {
        std::cmp::Ordering::Equal => a.send.cmp(&b.send),
        ord => ord,
    });
    let mut other = Vec::new();
    data.clone_into(&mut other);
    data.retain(|f| {
        other
            .iter()
            .any(|g| {
                g != f
                    && checkoverlap(&g.getposrange(), &f.getposrange())
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        blastresult.setstatus();
    }
}
pub(crate) fn matchmotif<T>(
    subject: &Path,
    species: T,
    locus: Option<&Locus>,
) -> io::Result<Vec<Blastmatch>>
where
    T: AsRef<str>,
{
    let reference = match downloadmotifs().map(|a| {
        speciesandorphonfiltering(
            &a,
            locus,
            "motifs".to_string(),
            species.as_ref(),
            false,
            false,
        )
    }) {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Motifs from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast: Vec<Blast> =
        match blastcommand(reference.as_path(), subject, Blastlevel::default()) {
            Ok(b) => b,
            Err(e) => {
                return Err(io::Error::new(ErrorKind::InvalidData, e));
            }
        }
        .into_iter()
        .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    if let Some(a) = locus {
        locusfiltering(a, &mut blast);
    }
    blast.iter_mut().for_each(|p| {
        //Already performed
        //if p.sstart > p.send {
        //    (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
        //}
    });
    statusblastmotifs(&mut blast);
    //let _ = fs::remove_file(name);
    Ok(blast.into_iter().map(|f| f.into()).collect())
}
pub(crate) fn genesblast(
    subject: &[GeneInfos],
    args: &Args,
    species: &Species,
    locus: &Locus,
) -> io::Result<Vec<Blastmatch>> {
    let mut reader = getassemblyreader(args)?;
    let name = NamedTempFile::with_suffix("genes_blast.txt")?;
    let file = File::create(&name)?;
    let mut fastawriter = fasta::Writer::new(file);
    subject
        .iter()
        .map(|f| {
            let elem = f
                .extractsequence(&mut reader)
                .map_err(|p| io::Error::new(ErrorKind::InvalidInput, p));
            let bool = elem
                .as_ref()
                .map_or(false, |p| f.addtosequence(p, &mut fastawriter).is_err());

            if bool {
                Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    "Cannot print sequence",
                ))
            } else {
                elem
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reference = match downloadref(true).map(|(a, b)| {
        speciesandorphonfiltering(
            &a,
            Some(locus),
            b,
            &species.to_string(),
            false,
            args.cacheerase,
        )
    }) {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast = match blastcommand(
        reference.as_path(),
        &name.into_temp_path(),
        Blastlevel::default(),
    ) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    retainbestmatch(&mut blast);
    locusfiltering(locus, &mut blast);
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, Strand::Minus);
        }
    });
    statusblastvs(&mut blast);
    Ok(blast.into_iter().map(|f| f.into()).collect())
}
pub(crate) fn locuspos(
    search: &Locus,
    hap: &Haplotype,
    locus: &[LocusInfos],
    blast: &[Blastmatch],
) -> Option<(LocusInfos, Vec<Blastmatch>)> {
    let opt = locus
        .iter()
        .find(|p| &p.locus == search && &p.haplotype == hap)?;
    let fil = |a: &Blastmatch| {
        let opt = opt.clone();
        a.getlocusname() == Some(opt.locus)
            && a.sseqid == opt.contig
            && (opt.start.getobasedpos()..=opt.end.getobasedpos())
                .contains(&a.sstart.try_into().unwrap_or_default())
    };
    if !blast.iter().any(fil) {
        return None;
    }
    let mut bl = blast.to_vec();
    bl.retain(fil);
    Some((opt.clone(), bl))
}
pub(crate) fn filter_new_alleles<'a, T>(data: &'a [T], motifs: &[T]) -> impl Iterator<Item = &'a T>
where
    T: Blastcalc + Debug,
{
    data.iter().filter(|p| p.onlynewalleles()).filter(move |f| {
        motifs.iter().any(|p| {
            checkoverlap(&p.getposrange(), &f.getposrange())
                && p.getposrange() != f.getposrange()
                && p.getposrange().count() > f.getposrange().count()
        })
    })
}
pub(crate) fn blastcommand<T>(
    reference: T,
    subject: T,
    blastlevel: Blastlevel,
) -> io::Result<Vec<Blast>>
where
    T: AsRef<Path>,
{
    let (reference, subject) = (reference.as_ref(), subject.as_ref());
    if !reference.try_exists().unwrap_or(false) {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Reference file was not found",
        ));
    }
    if !subject.try_exists().unwrap_or(false) {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Subject file was not found",
        ));
    }
    let reference = &format!("{}", reference.display()).replace(" ", "_");
    let subject = &format!("{}", subject.display()).replace(" ", "_");
    let output = NamedTempFile::new()?;
    let output = &format!("{}", output.path().display());
    let command = Command::new("blastn")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .args([
            "-query",
            reference,
            "-subject",
            subject,
            "-perc_identity",
            "75",
            "-out",
            output,
            "-word_size",
            blastlevel.into_word_size().to_string().as_str(),
            "-max_target_seqs",
            blastlevel.into_max_matches().to_string().as_str(),
            "-max_hsps",
            blastlevel.into_max_matches().to_string().as_str(),
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send qlen length pident gaps sseq",
        ])
        .spawn()?;
    println!("Launching {} against {}", reference, subject);
    println!(
        "BLAST has been launched with id {}. Please wait.",
        command.id()
    );
    let outputc = command.wait_with_output()?;
    let mut result: Vec<Blast> = Vec::new();
    {
        if !outputc.status.success() {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!(
                    "BLAST has failed with status {}. Error is {}",
                    outputc.status,
                    &String::from_utf8_lossy(&outputc.stderr)
                ),
            ));
        } else {
            println!("BLAST has been done. Parsing.");
            let file = File::open(output)?;
            let reader = BufReader::new(file);
            let mut csv = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(reader);
            for record in csv.deserialize::<Blast>() {
                if let Ok(mut r) = record {
                    if r.sstart > r.send {
                        (r.sstart, r.send, r.complement) = (r.send, r.sstart, Strand::Minus)
                    }
                    result.push(r);
                } else if let Err(r) = record {
                    eprintln!("Error in {r}");
                    return Err(io::Error::new(ErrorKind::InvalidData, "error"));
                }
            }
        };
    };
    //let _ = fs::remove_file(output);
    Ok(result)
}
pub(crate) fn getspeciesfromncbi<T>(
    client: &reqwest::blocking::Client,
    species: &T,
) -> Result<(String, usize), Box<dyn Error>>
where
    T: AsRef<str>,
{
    let species = species.as_ref();
    let id = if let Ok(val) = species.parse::<usize>() {
        val
    } else {
        let response = client
            .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?")
            .query(&[("db", "taxonomy"), ("term", species), ("format", "json")])
            .send()?;
        if response.status().is_server_error() {
            return Err(Box::new(io::Error::new(
                ErrorKind::HostUnreachable,
                "Unsuccessful NCBI response",
            )));
        } else if !response.status().is_success() {
            return Err(Box::new(io::Error::new(
                ErrorKind::NotFound,
                "Unsuccessful NCBI response",
            )));
        }
        let jsone: json::Value = json::from_str(&response.text().unwrap_or(String::new()))?;
        jsone["esearchresult"]["idlist"]
            .as_array()
            .map(|f| f.iter().next())
            .ok_or(io::Error::new(
                ErrorKind::InvalidInput,
                "Species not found on NCBI Taxonomy.",
            ))?
            .ok_or(io::Error::new(
                ErrorKind::InvalidInput,
                "Species not found on NCBI Taxonomy.",
            ))?
            .as_str()
            .ok_or(io::Error::new(
                ErrorKind::InvalidInput,
                "No result for the term used",
            ))?
            .parse::<usize>()
            .unwrap_or(0)
    };
    let val = &format!("{}", id);
    let response = client
        .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")
        .query(&[("db", "taxonomy"), ("id", val), ("format", "json")])
        .send()?;
    let jsone: json::Value = json::from_str(&response.text().unwrap_or(String::new()))?;
    let mapping = jsone["result"].as_object().ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No result for the term used",
    ))?;
    if let Some(data) = mapping.get(val) {
        println!("Data is {:?}", data);
        let elem = data.as_object().ok_or(io::Error::new(
            ErrorKind::InvalidInput,
            "No result for the id found",
        ))?;
        let text = String::from("species");
        let subspecies = String::from("subspecies");
        match (elem.get("rank"), elem.get("scientificname")) {
            (Some(Value::String(text2)), Some(Value::String(name)))
                if text2.as_str() == text.as_str() || text2.as_str() == subspecies.as_str() =>
            {
                Ok((name.to_string(), id))
            }
            (Some(Value::String(rank)), _) => Err(Box::new(io::Error::new(
                ErrorKind::InvalidInput,
                format!("The term used is not a (sub)species but a {}", rank),
            ))),
            _ => Err(Box::new(io::Error::new(
                ErrorKind::InvalidInput,
                "Invalid term used",
            ))),
        }
    } else {
        Err(Box::new(io::Error::new(
            ErrorKind::InvalidInput,
            "Invalid mapping",
        )))
    }
}
/*
pub(crate) fn launchblast(
    species: &str,
    locus: &Locus,
    subject: &Path,
    args: &Args,
) -> Result<Vec<Newfasta>, ()> {
    let (realspecies, _) = match getspeciesfromncbi(&REQUESTCLIENT, &species) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    println!("The species is {}.", realspecies.);
    let (path, releaseversion) = match downloadref(true).map(|(a, b)| {
        (
            speciesandorphonfiltering(
                &a,
                Some(locus),
                b.clone(),
                species.as_ref(),
                false,
                args.cacheerase,
            ),
            b,
        )
    }) {
        None => return Err(()),
        Some((Err(a), ..)) => return Err(()),
        Some((Ok(a), b)) => (a, b),
    };
    let filtering = match speciesandorphonfiltering(
        &path,
        Some(locus),
        releaseversion,
        &realspecies,
        false,
        args.cacheerase,
    ) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    let result: Vec<Newfasta> =
        match blastcommand(filtering.as_path(), subject, Blastlevel::default()) {
            Err(e) => {
                eprintln!("{e}");
                return Err(());
            }
            Ok(mut c) => {
                statusblast(&mut c);
                let status: Vec<Newfasta> = c
                    .iter()
                    .map(|c| {
                        let b: &dyn Blastcalc = c;
                        Newfasta::from(b)
                    })
                    .collect();
                let sequence = status.iter().fold(String::new(), |mut acc, f| {
                    let f: &dyn Seqresult = f;
                    acc.push_str(&format!("\n{}", f));
                    acc
                });
                if let Err(e) = fs::write("sequence_finished.fasta", sequence) {
                    eprintln!("An error has occured while priting sequence: {e}.");
                    return Err(());
                }
                status
            }
        };
    if result.is_empty() {
        eprintln!(
            "BLAST result is empty. No new alleles matching IMGT criteria has been identified."
        );
        return Err(());
    }
    Ok(result)
} */
pub(crate) fn submit(
    args: &Args,
    locus: &[crate::LocusInfos],
    c: &[Blastmatch],
    realspecies: &Species,
) -> Result<(), String> {
    //let result: Vec<Newfasta> = c.into_iter().map(Newfasta::newfromblastowner).collect();
    let dir = Path::new(&args.outdir);
    /* if dir.is_dir() {
        eprintln!("Archive directory exists, going to be deleted.");
        if let Err(e) = fs::remove_dir_all(&dir) {
            let dir = dir.display();
            return Err(format!("Cannot remove the directory {dir}, error is {e}."));
        }
    } */
    let dirtemp = match TempDir::new_in(dir) {
        Err(e) => return Err(format!("Cannot create archive directory, error is {e}")),
        Ok(p) => p,
    };
    let dir = dirtemp.path();
    let _ = match args.assembly.as_ref().map(|_| getassemblyreader(args)) {
        None => return Err("No assembly provided.".to_string()),
        Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
        Some(Ok(b)) => b,
    };
    /* let locuspos =
        File::create(dir.join("newloc.csv")).map_err(|f| format!("Locus csv, error is {f}"))?;
    locuspos
        .lock()
        .map_err(|f| format!("Error acquiring lock. Error is {f}"))?;
    let mut csv = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(locuspos);
    for loci in locus.iter() {
        csv.serialize(loci)
            .map_err(|p| format!("Error serializing locus position: {p}"))?;
    }
    csv.flush()
        .map_err(|p| format!("Error serializing locus position: {p}"))?; */
    let lightbam = dir.join("outlight.bam");
    generatelightbam(args, &lightbam, locus)?;
    let sequencefile = dir.join("sequence.fasta");
    let mut fastawriter = fasta::Writer::to_file(&sequencefile).map_err(|f| format!("{f}"))?;
    for list in locus.iter() {
        let mut assembly = match args.assembly.as_ref().map(|_| getassemblyreader(args)) {
            None => return Err("No assembly provided.".to_string()),
            Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
            Some(Ok(b)) => b,
        };
        let seq = list
            .extractsequence(&mut assembly)
            .unwrap_or("Sequence is unavailable".to_string());
        fastawriter
            .write(
                &format!("{}:{}", list.locus, list.contig),
                Some(&format!(
                    "{}:{}-{}/{}/{}",
                    list.locus,
                    list.start.getobasedpos(),
                    list.end.getobasedpos(),
                    list.complement,
                    list.haplotype
                )),
                seq.as_bytes(),
            )
            .map_err(|f| {
                format!(
                    "Unable to write fasta sequence {}. Error is {f}",
                    list.locus
                )
            })?;
    }
    let mut motifs = matchmotif(&sequencefile, &realspecies.to_string(), None)
        .map_err(|f| format!("Error matching motifs: {f}").to_string())?;
    motifs.iter_mut().for_each(|p| {
        if let Some(find) = locus
            .iter()
            .find(|k| format!("{}:{}", k.locus, k.contig) == p.sseqid)
            && let Some((newstart, newend, newcomplement)) = find.positioninlocus(
                &Position::new(false, p.sstart.try_into().unwrap_or_default()),
                &Position::new(false, p.send.try_into().unwrap_or_default()),
                &p.complement,
            )
        {
            p.sstart = newstart.getobasedpos().try_into().unwrap_or_default();
            p.send = newend.getobasedpos().try_into().unwrap_or_default();
            p.complement = newcomplement;
        }
    });

    let mut c = c.to_vec();
    c.iter_mut().for_each(|p| {
        if let Some(loc) = p.getlocusname()
            && let Some(find) = locus.iter().find(|fi| {
                p.getchromosomefromsubject().is_some_and(|a| a == fi.contig) && fi.locus == loc
            })
            && p.sseqid.starts_with("GENE")
            && let Some((start, end, complement)) = p.getpositionfromsubject()
            && let Some((newstart, newend, newcomplement)) = find.locusinposition(
                &start,
                &end,
                &complement,
                &Position::new(false, p.sstart.try_into().unwrap_or_default()),
                &Position::new(false, p.send.try_into().unwrap_or_default()),
                &p.complement,
            )
        {
            p.sstart = newstart.getobasedpos().try_into().unwrap_or_default();
            p.send = newend.getobasedpos().try_into().unwrap_or_default();
            p.qseqid = Name::new(
                Some("GENE".to_string()),
                p.getallelename().to_string(),
                realspecies.to_string(),
                None,
                start,
                end,
                complement,
            );
            p.complement = newcomplement;
        }
    });
    let mut matche = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(dir.join("motifs.txt"))
        .map_err(|e| format!("Error setting motifs match: {e}"))?;
    for matches in motifs.iter() {
        let _ = matche.serialize(matches);
    }
    let _ = matche.flush();
    let mut matche = csv::WriterBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_path(dir.join("motifs2.txt"))
        .map_err(|e| format!("Error setting motifs match: {e}"))?;
    for matches in c.iter() {
        let _ = matche.serialize(matches);
    }
    let _ = matche.flush();
    let file = dir.join("newpotentialalleles.fasta");
    let sequence = filter_new_alleles(&c, &motifs).fold(String::new(), |mut acc, f| {
        let f: &dyn Blastcalc = f;
        let f: &dyn Seqresult = &Newfasta::from(f);
        acc.push_str(&format!("\n{}", f));
        acc
    });
    println!("BLAST results were added.");
    if let Err(e) = fs::write(file, &sequence) {
        eprintln!("An error has occured while priting sequence: {e}.");
        return Ok(());
    }
    if !args.nosubmit {
        if sequence.trim().is_empty() && args.mytoken.is_none() && !asknonewalleles() {
            println!("Exiting");
            return Ok(());
        }
        if args.mytoken.is_none() {
            browseropening().map_err(|f| f.to_string())?;
        }
        preparesubmission(dir, realspecies, args);
    }
    //let _ = fs::remove_dir_all(dir);
    Ok(())
    //form(&client);
}
pub fn asknonewalleles() -> bool {
    println!(
        "No new alleles has been found. IMGT can still process your data, do you want to continue (Y/n):"
    );
    readfromterminal(&'y', &'n', false)
}
pub(crate) fn askforsubmission(
    realspecies: &Species,
    locus: &[LocusInfos],
    args: &Args,
    infos: &HashMap<Locus, Vec<Blastmatch>>,
) -> io::Result<()> {
    let quest = REQUESTCLIENT
        .get(SUBMISSIONLINK.as_str())
        .query(&[("ping", "1")]);
    match quest.send() {
        Ok(p) if p.status().is_success() => (),
        Ok(p) if p.status().is_client_error() => {
            /* return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                "IMGT does not currently support automatic submissions.".to_string(),
            )); */

            //return Ok(());
        }
        Ok(p) => {
            return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                format!(
                    "Unable to contact IMGT servers at the moment. Please retry later. Error is {}",
                    p.status().canonical_reason().unwrap_or_default()
                ),
            ));
        }
        Err(e) => {
            return Err(io::Error::new(
                io::ErrorKind::NetworkUnreachable,
                format!("Unable to contact IMGT servers. Error is {}", e),
            ));
        }
    };
    if !args.nosubmit {
        if !infos
            .iter()
            .any(|(_, p)| p.iter().any(|p| p.onlynewalleles()))
        {
            if !asknonewalleles() {
                return Ok(());
            }
        } else if args.mytoken.is_none() {
            println!("Do you want to submit your sequences to IMGT (y to yes or n to no)?");
            if !readfromterminal(&'y', &'n', false) {
                println!("Your sequences won't be submitted.");
                return Ok(());
            }
        }
        let mut blastmatch: Vec<Blastmatch> = Vec::new();
        for (_, data) in infos.iter() {
            blastmatch.append(&mut data.clone());
        }
        submit(args, locus, &blastmatch, realspecies)
            .map_err(|f| io::Error::new(io::ErrorKind::InvalidInput, f.to_string()))?;
    }
    Ok(())
}
pub(crate) fn generatelightbam(
    args: &Args,
    light: &Path,
    locus: &[LocusInfos],
) -> Result<(), String> {
    println!("Generating small BAM for submission");
    let bam = if let Ok(r) = getreaderoffile(args) {
        r
    } else {
        return Err("Cannot access BAM file for light bam.".to_string());
    };
    let mut writer = if let Ok(files) = bam::Writer::from_path(
        light,
        &bam::Header::from_template(bam.header()),
        bam::Format::Bam,
    ) {
        files
    } else {
        let file = light.display();
        return Err(format!("Cannot create file {file} for light bam."));
    };
    for f in locus.iter() {
        let mut bam = if let Ok(r) = getreaderoffile(args) {
            r
        } else {
            return Err("Cannot access BAM file for light bam.".to_string());
        };
        if bam
            .fetch((
                f.contig.as_bytes(),
                f.start.getzbasedpos(),
                f.end.getzbasedpos().saturating_add(1),
            ))
            .is_err()
        {
            return Err("Cannot read BAM file region for light bam.".to_string());
        }
        for read in bam.rc_records().filter_map(Result::ok) {
            if writer.write(&read).is_err() {
                return Err("Cannot read BAM file region for light bam.".to_string());
            };
        }
    }
    Ok(())
}
pub(crate) fn browseropening() -> io::Result<()> {
    let link = SUBMISSIONLINK.as_str();
    println!(
        "Opening web browser to continue submission. Type (Y) to open the web browser, (e) to exit or (n) and go by yourself to {link}."
    );
    if readfromterminal(&'y', &'e', false) {
        let _ = webbrowser::open(link);
    } else {
        return Err(io::Error::new(
            ErrorKind::ConnectionAborted,
            "Your sequences won't be submitted",
        ));
    }
    Ok(())
}
pub(crate) fn createarchive(args: &Args, dir: &Path) -> io::Result<NamedTempFile> {
    let temp = tempfile::NamedTempFile::with_suffix_in("submission.tar.gz", dir)?;
    let file = File::create(&temp)?;
    let archive = GzEncoder::new(file, Compression::best());
    let mut tar = tar::Builder::new(archive);
    tar.append_dir_all(".", &args.outdir)?;
    tar.finish()?;
    Ok(temp)
}
pub(crate) fn preparesubmission(path: &Path, species: &Species, args: &Args) -> bool {
    let token = if let Some(a) = &args.mytoken {
        a.to_string()
    } else {
        println!(
            "Please give the token provided by the submission. Type (e) to exit. Analysis files and new sequences would be sent to IMGT for submission."
        );
        loop {
            let mut block = io::stdin().lock().take(100);
            let mut string = String::new();
            let val = block.read_line(&mut string);
            match (val, string.trim().to_lowercase()) {
                (Err(_), _) | (Ok(0), _) => {
                    eprintln!("Exiting");
                    return false;
                }
                (Ok(_), v) if v == "e" => {
                    eprintln!("Exiting");
                    return false;
                }
                (Ok(_), v) if v.chars().all(|f| f.is_ascii_hexdigit()) && v.len() == 10 => break v,
                (Ok(_), _) => {
                    eprintln!("Token is 10 hexadecimal characters, please retry.");
                }
            };
        }
    };
    let archive = match createarchive(args, path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error filling archive: {e}");
            return false;
        }
    };
    if let Err(e) = submission(&token, species, archive) {
        eprintln!("An error has occured during submission: {e}. Please retry later.");
        return false;
    }
    println!(
        "Your submission has been made successfully. Thank you for submitting your sequences to IMGT. A confirmation email has been sent."
    );
    true
}
pub(crate) fn submission(token: &str, species: &Species, archive: NamedTempFile) -> io::Result<()> {
    /* let multipart = reqwest::blocking::multipart::Form::new()
    .file("genelist", "submission/genelist.csv")?
    .file("sequences", "submission/sequences.txt")?
    .file("locuspos", "submission/locus.txt")?
    .text("type","submission")
    .file("locus", "submission/locus.bam")?; */
    let zip = reqwest::blocking::multipart::Form::new()
        .file("archive", archive.path())?
        .text("version", VERSION)
        .text("type", "submission")
        .text("species", species.to_string())
        .text("validspecies", species.ischecked().to_string());
    match REQUESTCLIENT
        .post(SUBMISSIONLINK.as_str())
        .bearer_auth(token)
        .multipart(zip)
        .send()
    {
        Ok(a) => match a.status() {
            StatusCode::OK | StatusCode::NO_CONTENT => Ok(()),
            StatusCode::UNAUTHORIZED | StatusCode::FORBIDDEN => Err(io::Error::new(
                io::ErrorKind::PermissionDenied,
                format!(
                    "The token is invalid or has expired. Please retry a submission. Response is {}",
                    a.text().unwrap_or_default()
                ),
            )),
            StatusCode::BAD_REQUEST => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "The data sent is invalid, please check submission and retry. Response is {}",
                    a.text().unwrap_or_default()
                ),
            )),
            a if a.is_server_error() => Err(io::Error::new(
                io::ErrorKind::ConnectionAborted,
                "The server is unavailable. Please retry a submission",
            )),
            e => Err(io::Error::new(
                io::ErrorKind::ConnectionReset,
                format!("An unexpected error has occured, please retry later. Error is {e}"),
            )),
        },
        Err(e) if e.is_timeout() => Err(io::Error::new(
            io::ErrorKind::TimedOut,
            "The submission has timed out. Please retry later.",
        )),
        Err(e) => Err(io::Error::new(
            io::ErrorKind::NetworkUnreachable,
            format!("An unexpected error has occured. Please retry later. Error is {e}"),
        )),
    }?;
    std::fs::remove_file(archive)?;
    Ok(())
}

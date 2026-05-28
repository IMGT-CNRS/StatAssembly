#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
use crate::r#struct::{
    Args, Blast, Blastcalc, Blastlevel, Blastmatch, FakeLocusinfo, GeneInfos, Haplotype, Locus,
    LocusInfos, Newfasta, Position, Seqresult, Status, Strand,
};
use crate::{BORNES, getassemblyreader, getreaderoffile, printpotentialbornes};
use bio::io::fasta;
use extended_htslib::bam::{self, Read};
use flate2::{Compression, write::GzEncoder};
use itertools::{Either, Itertools};
use lazy_static::lazy_static;
use minimap2::{Aligner, Mapping};
use reqwest::{StatusCode, tls};
use serde_json::{self as json, Value};
use std::cmp::{Ordering, max, min};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Debug;
use std::io::IsTerminal;
use std::str::FromStr;
use std::thread::{self};
use std::{
    env::{self, temp_dir},
    error::Error,
    fs::{self, File},
    io::{self, BufRead, BufReader, ErrorKind, Read as _, Write},
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
        obfstr::obfstring!("/statassembly/downloadmotif.txt")
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
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub(crate) const DELIMITERFASTA: char = '/';
pub(crate) const LOCUSSEPARATOR: usize = 1_000_000;

pub(crate) fn sendresult(request: &reqwest::blocking::Client, url: &str) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                Ok(e.text().unwrap_or("Error getting data".to_string()))
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn sendresultcompressed(
    request: &reqwest::blocking::Client,
    url: &str,
) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                let mut string = String::new();
                if let Ok(bytes) = e.bytes()
                    && let mut decode = flate2::read::GzDecoder::new(&bytes[..])
                    && let Err(e) = decode.read_to_string(&mut string)
                {
                    return Err(format!("Cannot decrypt motifs: {e}"));
                }
                if string.trim().is_empty() {
                    return Err(format!("Cannot decrypt motifs, invalid body."));
                }
                string.shrink_to_fit();
                Ok(string)
            } else {
                Err(format!("Error getting URL. Code is {}", e.status()))
            }
        }
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
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
pub(crate) fn getallelefromblast<T>(text: T) -> Option<String>
where
    T: AsRef<str>,
{
    text.as_ref()
        .to_ascii_uppercase()
        .split('|')
        .nth(1)
        .map(|f| f.to_string())
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
    let mut lines: Vec<&str> = Vec::with_capacity(text.lines().count());
    let mut keep = true;
    for line in text.lines() {
        let result = line.trim();
        if result.starts_with(">") {
            keep = match (present, species) {
                (true, true) => {
                    getnamefromblast(result).is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
                }
                (false, true) => !getnamefromblast(result)
                    .is_some_and(|f| f.to_lowercase().contains(&find.to_ascii_lowercase())),
                (false, false) => !getallelefromblast(result)
                    .is_some_and(|f| f.to_lowercase().contains(&find.to_ascii_lowercase())),
                (true, false) => getallelefromblast(result)
                    .is_some_and(|f| f.to_lowercase().contains(&find.to_ascii_lowercase())),
            };
        }
        if keep {
            lines.push(result);
        }
    }
    lines
        .iter()
        .fold(String::new(), |mut acc, f| {
            acc.push_str(format!("\n{}", f).as_str());
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
pub(crate) fn downloadbornes() -> Option<PathBuf> {
    println!("Downloading bornes sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "bornes.fasta");
    if !fileincache(&tempfile) {
        println!("Downloading IMGT bornes, please wait...");
        match sendresultcompressed(&REQUESTCLIENT, BORNESLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    _ => {
                        eprintln!("Cannot write bornes in sequence.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT bornes data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!("IMGT bornes were already downloaded, retrieving...");
    };
    Some(tempfile)
}
pub(crate) fn downloadmotifs() -> Option<PathBuf> {
    println!("Downloading motifs sequence from IMGT");
    let tempfile = Path::join(&env::temp_dir(), "motifs.fasta");
    if !fileincache(&tempfile) {
        println!("Downloading IMGT motifs, please wait...");
        match sendresultcompressed(&REQUESTCLIENT, MOTIFLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    _ => {
                        eprintln!("Cannot write refseq in sequence.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT motifs data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!("IMGT motifs were already downloaded, retrieving...");
    };
    Some(tempfile)
}
pub(crate) fn downloadref(allowdownload: bool) -> Option<(PathBuf, String)> {
    println!("Checking reference sequence from IMGT/GENE-DB");
    let releaseversion = match sendresult(&REQUESTCLIENT, RELEASELINK.as_str()) {
        Ok(e) => e,
        Err(e) => {
            println!("Release fetched failed because: {e}");
            return None;
        }
    };
    if !allowdownload {
        return Some((PathBuf::new(), releaseversion));
    };
    let tempfile = Path::join(&env::temp_dir(), format!("refseq{}.fasta", releaseversion));
    if !tempfile.is_file() {
        println!(
            "Downloading IMGT/GENE-DB release {}, please wait...",
            releaseversion
        );
        match sendresult(&REQUESTCLIENT, VQUESTLINK.as_str()) {
            Ok(e) => {
                match File::create(&tempfile).map(|mut f| f.write_all(e.as_bytes())) {
                    Ok(Ok(_)) => (),
                    _ => {
                        eprintln!("Cannot write refseq in sequence.");
                        return None;
                    }
                }
                println!("Success.")
            }
            Err(e) => {
                println!("IMGT/GENE-DB data retrieval failed because: {e}.");
                return None;
            }
        }
    } else {
        println!(
            "Release {} was already downloaded, retrieving...",
            releaseversion
        );
    };
    Some((tempfile, releaseversion))
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
            p.getquery().contains(&format!("{}", locus))
                || p.getquery().contains(&"TRD".to_string())
        });
    } else {
        blast.retain(|p| p.getquery().contains(&format!("{}", locus)));
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
    if outfile.exists() && !force {
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
    data.sort_unstable_by(|a, b| match a.sseqid.cmp(&b.sseqid) {
        std::cmp::Ordering::Equal => a.qseqid.cmp(&b.qseqid),
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
pub(crate) fn find_global_best_range(
    blastcheck: &[Blastmatch],
    bornes: &Option<Vec<Blastmatch>>,
) -> Option<Vec<LocusInfos>> {
    let mut groups: HashMap<(Locus, Haplotype), (String, bool, usize, usize)> = HashMap::new();
    let blastcheck = blastcheck.to_vec();
    let mut hash = blastcheck
        .iter()
        .filter_map(|p| p.getlocusname().map(|d| (d, p)))
        .into_group_map_by(|(l, f)| (l.clone(), f.getsubject()));
    hash.iter_mut().for_each(|(_, v)| {
        v.sort_unstable_by(|(_, a), (_, b)| match a.getsubject().cmp(&b.getsubject()) {
            Ordering::Equal => match a.getpos().0.cmp(&b.getpos().0) {
                Ordering::Equal => b.getpos().1.cmp(&b.getpos().1).reverse(),
                ord => ord,
            },
            ord => ord,
        })
    });
    let mut locus: HashMap<(Locus, String, u32), std::collections::BTreeSet<Blastmatch>> =
        HashMap::new();
    for ((loci, sseqname), data) in hash.iter() {
        let mut index = 0;
        for (_, elem) in data {
            let d = match locus.get_mut(&(loci.clone(), sseqname.to_string(), index)) {
                Some(a) => a,
                None => {
                    let mut e = BTreeSet::new();
                    e.insert((*elem).clone());
                    locus.insert((loci.clone(), sseqname.to_string(), index), e);
                    continue;
                }
            };
            if let Some(b) = d.last()
                && max(b.sstart, b.send).abs_diff(min(elem.send, elem.sstart)) <= LOCUSSEPARATOR
            {
                d.insert((*elem).clone());
            } else {
                index += 1;
                let mut e = BTreeSet::new();
                e.insert((*elem).clone());
                locus.insert((loci.clone(), sseqname.to_string(), index), e);
            }
        }
    }
    if let Some(borneblast) = &bornes {
        'borne: for borneblast in borneblast {
            for (_, values) in locus.iter_mut() {
                for val in values.iter() {
                    let a = borneblast;
                    let b = val;
                    if max(b.sstart, b.send).abs_diff(min(a.send, a.sstart)) <= BORNES {
                        values.insert(borneblast.clone().into());
                        break 'borne;
                    }
                }
            }
        }
    }
    let real =
        locus
            .into_iter()
            .sorted_unstable_by(|((al, ass, _), av), ((bl, bs, _), bv)| match al.cmp(bl) {
                Ordering::Equal => match av.iter().count().cmp(&bv.iter().count()) {
                    Ordering::Equal => ass.cmp(bs),
                    ord2 => ord2.reverse(),
                },
                ord => ord,
            });
    for ((loci, sseq, _), blast) in real {
        if !groups.contains_key(&(loci.clone(), Haplotype::Primary)) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert((loci, Haplotype::Primary), (sseq, complement, min, max));
        } else if !groups.contains_key(&(loci.clone(), Haplotype::Alternate)) {
            let (min, max) = match (blast.first(), blast.last()) {
                (Some(a), Some(b)) => (min(a.sstart, a.send), max(b.sstart, b.send)),
                _ => unreachable!("vec is not empty"),
            };
            let split = blast.iter().into_group_map_by(|f| &f.complement);
            let complement = split.get(&Strand::Minus).map_or(0, |s| s.len())
                > split.get(&Strand::Plus).map_or(0, |s| s.len());
            groups.insert((loci, Haplotype::Alternate), (sseq, complement, min, max));
        }
    }
    Some(
        groups
            .into_iter()
            .map(|((locus, hap), (sseq, complement, start, end))| {
                LocusInfos::new(
                    locus,
                    hap,
                    sseq,
                    Position::new(
                        false,
                        start.saturating_sub(BORNES).try_into().unwrap_or_default(),
                    ),
                    Position::new(
                        false,
                        end.saturating_add(BORNES).try_into().unwrap_or_default(),
                    ),
                    if complement {
                        Strand::Minus
                    } else {
                        Strand::Plus
                    },
                )
            })
            .collect(),
    )
}
pub(crate) fn retainbestmatch(blast: &mut Vec<Blast>) {
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
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
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
    });
    statusblastmotifs(&mut blast);
    //let _ = fs::remove_file(name);
    Ok(blast.into_iter().map(|f| f.into()).collect())
}
pub(crate) fn genesblast<T>(
    subject: &[GeneInfos],
    args: &Args,
    species: T,
    locus: &Locus,
) -> io::Result<Vec<Blastmatch>>
where
    T: AsRef<str>,
{
    let mut reader = getassemblyreader(args)?;
    let name = temp_dir().join("genes_blast.txt");
    let file = File::create(&name)?;
    let mut fastawriter = fasta::Writer::new(file);
    subject
        .iter()
        .map(|f| {
            let elem = f
                .extractsequence(&mut reader)
                .map_err(|p| io::Error::new(ErrorKind::InvalidInput, p));
            let bool = elem
                .iter()
                .any(|p| f.addtosequence(p, &mut fastawriter).is_err());

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
        speciesandorphonfiltering(&a, Some(locus), b, species.as_ref(), false, args.cacheerase)
    }) {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast: Vec<Blast> =
        match blastcommand(reference.as_path(), name.as_path(), Blastlevel::default()) {
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
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
    });
    statusblastvs(&mut blast);
    let _ = fs::remove_file(name);
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
pub(crate) fn locusallposition<T>(
    subject: &Path,
    species: T,
    args: &Args,
) -> io::Result<(Vec<LocusInfos>, Vec<Blastmatch>)>
where
    T: AsRef<str>,
{
    let (bornespath, referencepath) = match (
        downloadbornes(),
        downloadref(true).map(|(a, b)| {
            speciesandorphonfiltering(&a, None, b, species.as_ref(), true, args.cacheerase)
        }),
    ) {
        (Some(a), Some(b)) => (a, b?),
        (None, ..) => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Bornes from IMGT cannot be downloaded",
            ));
        }
        (.., None) => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    println!("Merging bornes and IMGT/GENE-DB data.");
    /* let file = NamedTempFile::new_in(&env::temp_dir())?;
    let mut merge = File::create(&file)?;
    let read = fs::read_to_string(&reference.0)?;
    let read2 = fs::read_to_string(&reference.1)?;
    merge.write_all(read.as_bytes())?;
    merge.write_all(read2.as_bytes())?;
    println!("Blasting to get position of loci."); */
    let (blast, bornes) = thread::scope(|s| {
        let blast = s.spawn(|| {
            let b = match blastcommand(
                referencepath.clone(),
                subject.to_path_buf(),
                Blastlevel::Normal,
            ) {
                Ok(b) => b,
                Err(e) => {
                    return Err(io::Error::new(ErrorKind::InvalidData, e));
                }
            };
            Ok(b.into_iter().collect())
        });
        let bornes = if args.nobornes {
            None
        } else {
            Some(s.spawn(|| {
                let aligner = Aligner::builder()
                    .asm20()
                    .with_cigar()
                    .with_index_threads(args.threads)
                    .with_index(subject, None)
                    .ok()?;
                let read = fasta::Reader::from_file(bornespath.clone()).ok()?;
                let mut aligns: Vec<Blastmatch> = Vec::new();
                for seq in read.records().filter_map(Result::ok) {
                    if let Ok(a) = aligner
                        .map(
                            seq.seq(),
                            false,
                            false,
                            None,
                            None,
                            Some(format!("NCBI|{}", seq.desc().unwrap_or_default()).as_bytes()),
                        )
                        .map(|mut f| {
                            f.retain(|p| p.is_primary);
                            f
                        })
                        && let Some(first) = a.first()
                        && let (
                            Some(qseqid),
                            Some(sseqid),
                            sseq,
                            Ok(sstart),
                            Ok(send),
                            complement,
                            status,
                            Some(Some(identity)),
                        ) = (
                            &first.query_name,
                            &first.target_name,
                            String::new(),
                            TryInto::<usize>::try_into(first.target_start),
                            TryInto::<usize>::try_into(first.target_end),
                            if first.strand == minimap2::Strand::Reverse {
                                Strand::Minus
                            } else {
                                Strand::Plus
                            },
                            Status::New,
                            (&first)
                                .alignment
                                .as_ref()
                                .map(|p| p.alignment_score.map(|d| d as f32)),
                        )
                    {
                        let matche = Blastmatch::new(
                            qseqid.to_string(),
                            sseqid.to_string(),
                            sseq,
                            sstart,
                            send,
                            complement,
                            status,
                            identity,
                        );
                        aligns.push(matche);
                    }
                }
                if aligns.is_empty() {
                    None
                } else {
                    Some(aligns)
                }
            }))
        };
        println!("Waiting for blast and bornes to finish.");
        (blast.join(), bornes.map(|d| d.join()))
    });
    let (mut blast, mut bornes) = match (blast, bornes) {
        (Ok(Ok(a)), Some(Ok(Some(b)))) => (a, Some(b)),
        (Ok(Ok(a)), Some(Ok(None))) | (Ok(Ok(a)), None) => {
            eprintln!("No bornes found.");
            (a, None)
        }
        (Ok(Err(a)), _) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Blast error. Error is: {}", a),
            ));
        }
        (Err(a), Some(Err(b))) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!(
                    "Blast and bornes could not be performed. Error are: {:?} and {:?}",
                    a, b
                ),
            ));
        }
        (Err(a), _) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Blast could not be performed. Error is: {:?}", a),
            ));
        }
        (_, Some(Err(b))) => {
            return Err(io::Error::new(
                ErrorKind::BrokenPipe,
                format!("Bornes could not be performed. Error is: {:?}", b),
            ));
        }
    };
    bornes.as_mut().map(|p| {
        p.sort_unstable_by(|f, g| match f.getallelename().cmp(&g.getallelename()) {
            Ordering::Equal => f.getidentity().total_cmp(&g.getidentity()).reverse(),
            e => e,
        });
        p.dedup_by(|f, g| match (f.getallelename(), g.getallelename()) {
            (Some(a), Some(b)) if a.eq_ignore_ascii_case(&b) => true,
            _ => false,
        });
        p.sort_unstable();
        p.shrink_to_fit();
    });
    if let Some(bornes2) = &bornes {
        printpotentialbornes(bornes2, &args)?;
    }
    //Filter by locus
    retainbestmatch(&mut blast);
    blast.shrink_to_fit();
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
        p.setstatus();
    });
    let mut statusvec: BTreeMap<(String, usize, usize), Blast> = BTreeMap::new();
    for elem in blast.into_iter() {
        if let Some((k, b)) = statusvec.iter().find(|((s, r1, r2), b)| {
            b != &&elem
                && s.as_str() == elem.sseqid.as_str()
                && checkoverlap(&(*r1..=*r2), &(elem.sstart..=elem.send))
        }) {
            let scoreactual = b.length / b.qlen * b.pident.powi(2) as usize;
            let newscore = elem.length / elem.qlen * elem.pident.powi(2) as usize;
            if newscore > scoreactual {
                statusvec.insert(k.clone(), elem);
            }
        } else {
            statusvec.insert((elem.sseqid.clone(), elem.sstart, elem.send), elem);
        }
    }
    let mut data: Vec<Blastmatch> = statusvec.into_values().map(|v| v.into()).collect();
    //Sort by name then starting then ending position
    data.sort_unstable();
    data.dedup();
    let writer = csv::WriterBuilder::new()
        .delimiter(b',')
        .comment(Some(b'#'))
        .has_headers(true)
        .from_path(Path::join(&temp_dir(), "test2.txt"));
    if let Ok(mut r) = writer {
        for elem in data.iter() {
            let _ = r.serialize(elem);
        }
    }
    let range = find_global_best_range(&data, &bornes).ok_or(io::Error::new(
        ErrorKind::InvalidInput,
        "No locus found after BLAST analysis",
    ));
    range.map(|f| (f, data))
}
pub(crate) fn filter_new_alleles<'a, T>(data: &'a [T], motifs: &[T]) -> impl Iterator<Item = &'a T>
where
    T: Blastcalc + Debug,
{
    data.iter().filter(|p| p.onlynewalleles()).filter(move |f| {
        motifs.iter().any(|p| {
            checkoverlap(&p.getposrange(), &f.getposrange()) && p.getposrange() != f.getposrange()
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
    if !reference.exists() {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "Reference file was not found",
        ));
    }
    if !subject.exists() {
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
            for record in csv.deserialize() {
                if let Ok(r) = record {
                    result.push(r);
                } else if let Err(r) = record {
                    eprintln!("Error in {r}");
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
) -> Result<String, Box<dyn Error>>
where
    T: AsRef<str>,
{
    let species = species.as_ref();
    let val = if let Ok(val) = species.parse::<usize>() {
        val
    } else {
        let response = client
            .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?")
            .query(&[("db", "taxonomy"), ("term", species), ("format", "json")])
            .send()?;
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
    let val = &format!("{}", val);
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
                Ok(name.to_string())
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
pub(crate) fn launchblast(
    species: &str,
    locus: &Locus,
    subject: &Path,
    args: &Args,
) -> Result<Vec<Newfasta>, ()> {
    let realspecies = match getspeciesfromncbi(&REQUESTCLIENT, &species) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    println!("The species is {}.", realspecies);
    let (path, releaseversion) = match downloadref(true) {
        None => return Err(()),
        Some((a, b)) => (a, b),
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
}
pub(crate) fn submit(
    args: &Args,
    locus: &[crate::LocusInfos],
    c: &[Blastmatch],
    realspecies: String,
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
    let mut motifs = matchmotif(&sequencefile, &realspecies, None)
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
                p.getchromosomefromblast().is_some_and(|a| a == fi.contig) && fi.locus == loc
            })
            && p.sseqid.starts_with("GENE")
            && let Some((start, end, complement)) = p.getposition()
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
            p.qseqid = p.getchromosomefromblast().unwrap_or_default();
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
    realspecies: &str,
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
        submit(args, locus, &blastmatch, realspecies.to_string())
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
pub(crate) fn preparesubmission(path: &Path, species: String, args: &Args) -> bool {
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
pub(crate) fn submission(token: &str, species: String, archive: NamedTempFile) -> io::Result<()> {
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
        .text("species", species);
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

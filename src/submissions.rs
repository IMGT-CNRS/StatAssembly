#[deny(clippy::unwrap_used)]
#[deny(clippy::expect_used)]
use bio::io::fasta;
use flate2::{Compression, write::GzEncoder};
use lazy_static::lazy_static;
use reqwest::{StatusCode, tls};
use serde::{Deserialize, Serialize};
use serde_json::{self as json, Value};
use std::cmp::{Ordering, max, min};
use std::collections::{BTreeMap, HashMap};
use std::str::FromStr;
use std::{
    borrow::Cow,
    env::{self, current_dir, temp_dir},
    error::Error,
    fs::{self, File},
    io::{self, BufRead, BufReader, ErrorKind, Read, Write},
    ops::{Not, Range, RangeInclusive},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{Duration, SystemTime},
};
use tempfile::NamedTempFile;

use crate::generatelightbam;
use crate::r#struct::{
    Args, Blast, Blastcalc, Blastmatch, Locus, Newfasta, Ourfasta, Seqresult, Status,
};
lazy_static! {
    pub static ref REQUESTCLIENT: reqwest::blocking::Client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::new(10, 0))
        .timeout(Duration::new(300, 0))
        .user_agent(format!("IMGT/StatAssembly version {}", VERSION))
        .referer(false)
        .tls_version_min(tls::Version::TLS_1_2)
        .https_only(true)
        .build()
        .unwrap();
}
const VERSION: &str = env!("CARGO_PKG_VERSION");
const RELEASELINK: &str = "https://www.imgt.org/download/GENE-DB/RELEASE";
const VQUESTLINK: &str = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP";
const SUBMISSIONLINK: &str = "https://www.imgt.org/submissions/";
pub(crate) const DELIMITERFASTA: char = '/';
const LOCUSSEPARATOR: usize = 1_000_000;

#[must_use]
pub(crate) fn sendresult(request: &reqwest::blocking::Client, url: &str) -> Result<String, String> {
    match request.get(url).send() {
        Ok(e) => {
            if e.status() == StatusCode::OK {
                Ok(e.text().unwrap_or("Error getting data".to_string()))
            } else {
                Err(e
                    .error_for_status()
                    .map_or("Error getting URL".to_string(), |f| {
                        format!("{}", f.status())
                    }))
            }
        }
        Err(e) => Err(format!("Error getting URL: {e}").to_string()),
    }
}
pub(crate) fn getnamefromblast(text: &str) -> Option<String> {
    text.to_ascii_lowercase()
        .split('|')
        .nth(2)
        .map(|f| f.to_string())
}
#[must_use]
pub(crate) fn fastafilter(text: &str, find: &str, present: bool) -> String {
    let mut lines: Vec<&str> = Vec::with_capacity(text.lines().count());
    let mut keep = true;
    for line in text.lines() {
        let result = line.trim();
        if result.starts_with(">") {
            keep = if present {
                getnamefromblast(result).is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
            } else {
                !getnamefromblast(result).is_some_and(|f| f.contains(&find.to_ascii_lowercase()))
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
pub(crate) fn checkifblastpresent() -> bool {
    println!("Detect if BLAST is operating.");
    let command = Command::new("blastn")
        .arg("-version")
        .stdout(Stdio::piped())
        .stdin(Stdio::piped())
        .status()
        .is_ok_and(|f| f.success());
    if !command {
        eprintln!("BLAST was not found. Check if present in PATH.");
        return false;
    } else {
        println!("BLAST is working. Continuing");
        return true;
    }
}
pub(crate) fn downloadref() -> Option<(PathBuf, String)> {
    println!("Downloading reference sequence from GENE-DB");
    let releaseversion = match sendresult(&REQUESTCLIENT, RELEASELINK) {
        Ok(e) => e,
        Err(e) => {
            println!("Release fetched failed because: {e}");
            return None;
        }
    };
    let tempfile = Path::join(&env::temp_dir(), format!("refseq{}.fasta", releaseversion));
    if !tempfile.is_file() {
        println!(
            "Downloading GENE-DB release {}, please wait...",
            releaseversion
        );
        match sendresult(&REQUESTCLIENT, VQUESTLINK) {
            Ok(e) => {
                let mut file =
                    File::create(&tempfile).expect("Cannot create refseq file in temp dir");
                file.write_all(e.as_bytes())
                    .expect("Cannot create refseq file in temp dir.");
                println!("Success.")
            }
            Err(e) => {
                println!("V-QUEST data retrieval failed because: {e}.");
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
pub(crate) fn locusfiltering(locus: &Locus, blast: &mut Vec<Blast>) {
    //TRD is inside TRA locus
    if locus == &Locus::TRA {
        blast.retain(|p| {
            p.qseqid.contains(&format!("{}", locus)) || p.qseqid.contains(&format!("TRD"))
        });
    } else {
        blast.retain(|p| p.qseqid.contains(&format!("{}", locus)));
    }
}
pub(crate) fn speciesandorphonfiltering(
    tempfile: &Path,
    releaseversion: String,
    species: &str,
    orphonfilter: bool,
) -> io::Result<PathBuf> {
    println!("Filtering based on species");
    let file = std::fs::read_to_string(&tempfile).unwrap();
    let info = fastafilter(&file, species, true).replace(" ", "_");
    let info = if orphonfilter {
        fastafilter(&info, "/OR", false)
    } else {
        info
    };
    let info = if info.is_empty() {
        println!("New species!!");
        file
    } else {
        info
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
    let tempfile = Path::join(
        &env::temp_dir(),
        format!(
            "refseq{}-{}.fasta",
            releaseversion.replace(" ", "_"),
            species.replace(" ", "_")
        ),
    );
    if let Err(e) = std::fs::write(&tempfile, &info) {
        return Err(e);
    }
    return Ok(tempfile);
}
pub(crate) fn readfastareader<T>(fasta: fasta::Reader<T>) -> io::Result<Vec<Ourfasta>>
where
    T: std::io::BufRead,
{
    let mut seqs = Vec::new();
    for record in fasta.records().filter_map(|p| p.ok()) {
        if let Err(v) = record.check() {
            eprintln!("Seq {} is invalid, skipped. Error is {v}", record.id());
            continue;
        }
        seqs.push(Ourfasta {
            name: format!("{} {}", record.id(), record.desc().unwrap_or_default()),
            seq: String::from_utf8_lossy(record.seq()).to_string(),
        })
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
pub(crate) fn readfastafile(seq: &Path) -> io::Result<Vec<Ourfasta>> {
    let file = File::open(seq)?;
    let fasta = fasta::Reader::new(file);
    readfastareader(fasta)
}
#[must_use]
#[allow(unused)]
pub(crate) fn selectnewalleles(result: &[Blast]) -> Vec<Ourfasta> {
    let mut fastas = Vec::new();
    for seq in result {
        let name = Ourfasta {
            name: seq.qseqid.clone(),
            seq: seq.sseq.clone(),
        };
        fastas.push(name);
    }
    fastas
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
                    && checkoverlap(&g.pos(), &f.pos())
                    && (g.length as f32 / g.qlen as f32 * g.pident.powf(1.1))
                        >= (f.length as f32 / f.qlen as f32 * f.pident.powf(1.1))
            })
            .not()
    });
    for blastresult in data {
        match (
            blastresult.pident,
            blastresult.gaps,
            blastresult.length,
            blastresult.qlen,
        ) {
            (100.0, 0, a, b) if a < b => blastresult.status = Status::Shorter,
            (100.0, 0, ..) => blastresult.status = Status::Equal,
            _ => blastresult.status = Status::New,
        }
    }
}
pub(crate) fn find_global_best_range(data: &[Blastmatch]) -> Option<(String, usize, usize)> {
    let mut groups: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    // Group by string, storing (pos, val) pairs
    for (s, pos, val) in data.iter().map(|f| (f.sseqid.clone(), f.sstart, f.send)) {
        groups.entry(s).or_insert_with(Vec::new).push((pos, val));
    }

    let mut global_best_name = String::new();
    let mut global_best_min = 0;
    let mut global_best_max = 0;
    let mut global_max_count = 0;

    // Process each group
    for (name, mut pairs) in groups {
        // Sort by the second usize (index 1)
        pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let mut max_count = 0;
        let mut best_min = 0;
        let mut best_max = 0;

        for left in 0..pairs.len() {
            let mut count = 0;
            let mut right = left;
            let right = loop {
                right += 1;
                if pairs.get(right).is_none() {
                    break right.saturating_sub(1);
                }
                // Move left to ensure window is valid
                if pairs[right].0.abs_diff(pairs[right.saturating_sub(1)].1) >= LOCUSSEPARATOR {
                    break right;
                }
                count += 1;
            };
            if count > max_count {
                max_count = count;
                best_min = pairs[left].0;
                best_max = pairs[right].1;
            }
        }

        // Update global best
        if max_count > global_max_count {
            global_max_count = max_count;
            global_best_name = name;
            global_best_min = best_min;
            global_best_max = best_max;
        }
    }

    if global_max_count > 0 {
        Some((
            global_best_name.to_string(),
            global_best_min,
            global_best_max,
        ))
    } else {
        None
    }
}
pub(crate) fn locusposition<T>(
    subject: &Path,
    species: T,
    locus: &Locus,
) -> io::Result<(String, RangeInclusive<usize>, Vec<Blastmatch>)>
where
    T: AsRef<str>,
{
    let reference = match downloadref()
        .map(|(a, b)| speciesandorphonfiltering(&a, b, species.as_ref(), false))
    {
        Some(a) => a?,
        None => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Reference from IMGT cannot be downloaded",
            ));
        }
    };
    let mut blast: Vec<Blast> = match blastcommand(reference.as_path(), subject) {
        Ok(b) => b,
        Err(e) => {
            return Err(io::Error::new(ErrorKind::InvalidData, e));
        }
    }
    .into_iter()
    .collect();
    //Filter by locus
    blast.retain(|f| f.length * 100 / f.qlen > 80 && f.pident >= 75.0);
    locusfiltering(&locus, &mut blast);
    blast.iter_mut().for_each(|p| {
        if p.sstart > p.send {
            (p.sstart, p.send, p.complement) = (p.send, p.sstart, true);
        }
        p.setstatus();
    });
    let mut statusvec: BTreeMap<(String, usize, usize), Blast> = BTreeMap::new();
    for elem in blast.into_iter() {
        if let Some((k, b)) = statusvec.clone().iter_mut().find(|((s, r1, r2), b)| {
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
    let mut data: Vec<Blastmatch> = statusvec
        .into_iter()
        .map(|(_k, v)| Blastmatch::new(v.qseqid, v.sseqid, v.sseq, v.sstart, v.send, v.status))
        .collect();
    //Sort by name then starting then ending position
    data.sort_unstable();
    data.dedup();
    let range = find_global_best_range(&data)
        .and_then(|(n, pos1, pos2)| Some((n, pos1..=pos2)))
        .ok_or(io::Error::new(
            ErrorKind::InvalidInput,
            "No locus found after BLAST analysis",
        ));
    if let Ok((name, dat)) = &range {
        data.retain(|p| p.sseqid.as_str() == name && dat.contains(&p.sstart))
    }
    range.map(|f| (f.0, f.1, data))
}
#[must_use]
pub(crate) fn filter_new_alleles<T>(data: &[T]) -> impl Iterator<Item = &T>
where
    T: Blastcalc,
{
    data.iter().filter(|p| p.onlynewalleles())
}
pub(crate) fn blastcommand<T>(reference: T, subject: T) -> io::Result<Vec<Blast>>
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
    let output = Path::join(&env::temp_dir(), "blast.txt");
    let output = &format!("{}", output.display());
    let command = Command::new("blastn")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .args([
            "-query",
            reference,
            "-subject",
            subject,
            "-out",
            output,
            "-max_hsps",
            "5",
            "-outfmt",
            "6 qseqid sseqid qstart qend sstart send qlen length pident gaps sseq",
        ])
        .spawn()?;
    println!("Launching {} against {}", reference, subject);
    println!("BLAST has been launched with id {}", command.id());
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
    let _ = fs::remove_file(&output);
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
            .ok_or(io::Error::new(ErrorKind::InvalidInput, "Idlist not found"))?
            .ok_or(io::Error::new(ErrorKind::InvalidInput, "Idlist not found"))?
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
        match (elem.get("rank"), elem.get("scientificname")) {
            (Some(Value::String(text2)), Some(Value::String(name)))
                if text2.as_str() == text.as_str() =>
            {
                Ok(name.to_string())
            }
            (Some(Value::String(rank)), _) => Err(Box::new(io::Error::new(
                ErrorKind::InvalidInput,
                format!("The term used is not a species but a {}", rank),
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
#[must_use]
pub(crate) fn launchblast(species: &str, subject: &Path) -> Result<Vec<Newfasta>, ()> {
    let realspecies = match getspeciesfromncbi(&REQUESTCLIENT, &species) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    println!("The species is {}.", realspecies);
    let (path, releaseversion) = match downloadref() {
        None => return Err(()),
        Some((a, b)) => (a, b),
    };
    let filtering = match speciesandorphonfiltering(&path, releaseversion, &realspecies, false) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(b) => b,
    };
    let result: Vec<Newfasta> = match blastcommand(filtering.as_path(), subject) {
        Err(e) => {
            eprintln!("{e}");
            return Err(());
        }
        Ok(mut c) => {
            statusblast(&mut c);
            let status: Vec<Newfasta> = c.iter().map(|b| Newfasta::newfromblast(b)).collect();
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
    let dir = Path::new(&current_dir().unwrap_or(temp_dir())).join("archive");
    if dir.is_dir() {
        return Err("Archive directory exists, remove directory to submit if needed.".to_string());
    }
    if let Err(e) = fs::create_dir(&dir) {
        return Err(format!("Cannot create archive directory, error is {e}"));
    }
    let _ = match args
        .assembly
        .as_ref()
        .map(|p| fasta::IndexedReader::from_file(&p))
    {
        None => return Err("No assembly provided.".to_string()),
        Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
        Some(Ok(b)) => b,
    };
    let lightbam = dir.join("outlight.bam");
    generatelightbam(args, &lightbam, locus)?;
    let sequencefile = dir.join("sequence.fasta");
    let mut fastawriter = fasta::Writer::to_file(sequencefile).map_err(|f| format!("{f}"))?;
    for list in locus {
        let assembly = match args
            .assembly
            .as_ref()
            .map(|p| fasta::IndexedReader::from_file(&p))
        {
            None => return Err("No assembly provided.".to_string()),
            Some(Err(e)) => return Err(format!("Error with assembly: {e}")),
            Some(Ok(b)) => b,
        };
        let seq = list.extractsequence(assembly).map_err(|f| format!("{f}"))?;
        fastawriter
            .write(
                &format!("{}", list.locus),
                Some(&format!(
                    "{}:{}-{}/{}/{}",
                    list.contig,
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
    let file = dir.join("newalleles.fasta");
    let sequence = result.fold(String::new(), |mut acc, f| {
        let f: &dyn Seqresult = f;
        acc.push_str(&format!("\n{}", f));
        acc
    });
    if let Err(e) = fs::write(file, sequence) {
        eprintln!("An error has occured while priting sequence: {e}.");
        return Ok(());
    }
    println!("BLAST results was added.");
    preparesubmission(&dir, realspecies);
    let _ = fs::remove_dir_all(dir);
    return Ok(());
    //form(&client);
}
pub(crate) fn browseropening() -> io::Result<()> {
    println!(
        "Opening web browser to continue submission. Type (Y) to open the web browser, (N) to refuse and go by yourself to {SUBMISSIONLINK} or (e) to exit."
    );
    let mut val = String::new();
    let _ = io::stdin().read_line(&mut val);
    let val = val.trim().to_ascii_lowercase();
    if val == "y" {
        let _ = webbrowser::open(SUBMISSIONLINK);
    } else if val == "e" {
        println!("Exiting");
        return Err(io::Error::from(ErrorKind::InvalidData));
    }
    Ok(())
}
pub(crate) fn createarchive(dir: &Path) -> io::Result<NamedTempFile> {
    let temp = tempfile::NamedTempFile::with_suffix("submission.tar.gz")?;
    let file = File::create(&temp)?;
    let archive = GzEncoder::new(file, Compression::default());
    let mut tar = tar::Builder::new(archive);
    tar.append_dir_all(
        "",
        dir.file_name()
            .map(|f| f.to_str())
            .unwrap_or_default()
            .unwrap_or_default(),
    )?;
    tar.finish()?;
    Ok(temp)
}
pub(crate) fn preparesubmission(path: &Path, species: String) -> bool {
    println!("Please give the token provided by the submission. Type (e) to exit.");
    let token = loop {
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
    };
    println!(
        "Analysis files and new sequences would be sent to IMGT for submission. Type (Y) to validate or something else to exit"
    );
    let mut val = String::new();
    let _ = io::stdin().read_line(&mut val);
    let val = val.trim().to_ascii_lowercase();
    let archive = match createarchive(path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error filling archive: {e}");
            return false;
        }
    };
    if val == "y" {
        if let Err(e) = submission(&token, species, archive) {
            eprintln!("An error has occured during submission: {e}. Please retry later.");
            return false;
        }
        println!("Your submission has been made. Thank you for submitting your sequences to IMGT.");
    } else {
        println!("Exiting");
        return false;
    }
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
        .post(SUBMISSIONLINK)
        .bearer_auth(token)
        .multipart(zip)
        .send()
    {
        Ok(a) => match a.status() {
            StatusCode::OK | StatusCode::NO_CONTENT => Ok(()),
            StatusCode::UNAUTHORIZED => Err(io::Error::new(
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
            StatusCode::SERVICE_UNAVAILABLE | StatusCode::INTERNAL_SERVER_ERROR => {
                Err(io::Error::new(
                    io::ErrorKind::ResourceBusy,
                    "The server is unavailable. Please retry a submission.",
                ))
            }
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

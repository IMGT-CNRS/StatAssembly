#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use std::{
        fs::File,
        io::{BufReader, ErrorKind::UnexpectedEof, Read as read2},
        path::PathBuf,
        process::Command,
        thread::sleep,
        time::Duration,
    };

    use crate::{
        extractgenelist, generategeneinfos, getorsetparams,
        identification::{sendresult, sendresultcompressed},
        locusposparser, posread,
        r#struct::{Args, Filecrea, HashMapinfo, Locus, Species, SpeciesError},
        submissions::{
            BORNESLINK, RELEASELINK, REQUESTCLIENT, checkifblastpresent, decompressseq,
            generatelightbam, generatesequence, generatesequenceraw, getprogressbarclassic,
        },
    };
    use clap::Parser;
    use extended_htslib::bam::{self, Read};
    use itertools::Itertools;
    use tempfile::{NamedTempFile, TempDir};
    #[test]
    fn checkblast() {
        assert!(checkifblastpresent())
    }
    #[test]
    #[allow(non_snake_case)]
    fn checkIMGTservers() {
        assert!(
            sendresult(&REQUESTCLIENT, &RELEASELINK).is_ok(),
            "Error with release link"
        );
        assert!(
            sendresultcompressed(&REQUESTCLIENT, &BORNESLINK).is_ok(),
            "Error with bornes link"
        );
    }
    #[test]
    fn generatesequencetest() {
        let testo = TempDir::new().unwrap();
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "-l",
            "example_files/CHM13v2.0loc.csv",
            "-s",
            "human",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "analyze",
        ];
        let human = Species::new("Homo sapiens").unwrap();
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let (mut locus, ..) = locusposparser(&args4, &human, false).unwrap();
        locus.sort_unstable(); //Because sorted in mergedlocus function
        let seq = generatesequence(&args4, testo.path(), true, &locus).unwrap();
        let mut bufone =
            BufReader::new(std::fs::File::open("example_files/results/sequence.fasta.gz").unwrap());
        let mut buftwo = BufReader::new(seq.getfile().unwrap());
        loop {
            let mut buf1 = [0; 2048];
            let mut buf2 = [0; 2048];
            match (bufone.read_exact(&mut buf1), buftwo.read_exact(&mut buf2)) {
                (Err(a), Err(b)) if b.kind() == UnexpectedEof && a.kind() == UnexpectedEof => break,
                (Ok(_), Ok(_)) if buf1.eq(&buf2) => (),
                (Ok(_), Ok(_)) => {
                    panic!("Buffer does not have same value.");
                }
                _ => panic!("Fails with buffer"),
            }
        }
    }
    #[test]
    fn checkspecies() {
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let human = Species::new("Homo sapiens").unwrap();
        assert!(human.ischecked(), "Species could not be validated");
        assert_eq!(human.getid(), Some(9606), "Species taxon is invalid");
        assert!(
            human.getrank().eq_ignore_ascii_case("species"),
            "Species has not a valid rank"
        );
        assert!(
            human.getname().eq_ignore_ascii_case("Homo sapiens"),
            "Homo sapiens is not valid"
        );
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let dog = Species::new("dog").unwrap();
        assert!(dog.ischecked(), "Species could not be validated");
        assert_eq!(dog.getid(), Some(9615), "Species taxon is invalid");
        assert!(
            dog.getrank().eq_ignore_ascii_case("subspecies"),
            "Subspecies has not a valid rank"
        );
        assert!(
            dog.getname().eq_ignore_ascii_case("Canis lupus familiaris"),
            "Dog is not valid"
        );
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let lamprey = Species::new("least brook lamprey");
        assert!(
            matches!(lamprey, Err(SpeciesError::Blocked)),
            "Lamprey should be blocked"
        );
    }
    #[test]
    fn testlocusparsing() {
        let mut val = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(
                "example_files/results/Homo_sapiens_IGH_NC_060938.1-primary_positionresult.csv",
            )
            .unwrap();
        let mut data: Vec<HashMapinfo> = Vec::with_capacity(5000);
        let mut count = 0;
        for line in val.deserialize() {
            count += 1;
            data.push(line.unwrap());
            if count >= 5000 {
                break;
            }
        }
        assert_eq!(data.iter().nth(3069).unwrap().qual, Some(80), "Issue");
        assert!(data.iter().nth(1022).unwrap().qual.is_none(), "Issue");
    }
    #[test]
    fn testlocuspositiononlyandbam() {
        let testo = TempDir::new().unwrap();
        let a = Filecrea::createtemp(None, Some("test.bam")).unwrap();
        let b = Filecrea::createtemp(None, Some("test.bam.csi")).unwrap();
        // Fake args as a Vec<&str>
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "-l",
            "example_files/CHM13v2.0loc.csv",
            "-s",
            "human",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "-z",
            a.getpath().to_str().unwrap(),
            "analyze",
        ];
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::new(&args4.species).unwrap();
        let (mut result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        let meanpath = PathBuf::from("example_files/results").join(".mean");
        assert!(meanpath.exists(), "No .mean in example_files");
        let mean = getorsetparams(&meanpath, &args4).unwrap();
        result.iter_mut().for_each(|p| {
            let pos = posread(&args4, &p).unwrap();
            p.setstatus(mean.mean, &args4, &pos);
        });
        assert!(
            result.iter().all(|f| f.status.getstatus().isvalid()),
            "One locus is not valid"
        );
        generatelightbam(
            &args4,
            args4.outlightbam.as_ref().unwrap(),
            Some(&b.getpath().to_path_buf()),
            &result,
        )
        .unwrap();
        let mut br = bam::Reader::from_path(args4.file.unwrap()).unwrap();
        let mut cr = bam::Reader::from_path(args4.outlightbam.unwrap().as_path()).unwrap();
        let mut activ = false;
        for (record, record2) in br
            .records()
            .filter_map(Result::ok)
            .zip_eq(cr.records().filter_map(Result::ok))
        {
            activ = true;
            assert_eq!(
                record,
                record2,
                "Record4 {} is different",
                String::from_utf8_lossy(record.qname())
            );
        }
        if !activ {
            panic!("No records found");
        }
    }
    #[test]
    fn compressanddecompress() {
        let testo = TempDir::new().unwrap();
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "-l",
            "example_files/CHM13v2.0loc.csv",
            "-s",
            "human",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "analyze",
        ];
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::new(&args4.species).unwrap();
        let (locus, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        // Generate sequence changes the reader so we need to compress and decompress.
        let compress = generatesequenceraw(&args4, &locus).unwrap();
        let compressreal = generatesequence(&args4, testo.path(), true, &locus).unwrap();
        let decompress = decompressseq(&compressreal).unwrap();
        let mut bufone = BufReader::new(compress);
        let mut buftwo = BufReader::new(decompress.getfile().unwrap());
        loop {
            let mut buf1 = [0; 4096];
            let mut buf2 = [0; 4096];
            match (bufone.read_exact(&mut buf1), buftwo.read_exact(&mut buf2)) {
                (Err(a), Err(b)) if b.kind() == UnexpectedEof && a.kind() == UnexpectedEof => break,
                (Ok(_), Ok(_)) if buf1.eq(&buf2) => (),
                (Ok(_), Ok(_)) => {
                    panic!("Buffer does not have same value.");
                }
                _ => panic!("Fails with buffer"),
            }
        }
    }
    #[ignore]
    #[test]
    fn testlocuspositionandgenes() {
        let testo = TempDir::new().unwrap();
        let a = NamedTempFile::with_suffix("test4.bam").unwrap();
        // Fake args as a Vec<&str>
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "-l",
            "example_files/CHM13v2.0loc.csv",
            "-s",
            "human",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "-z",
            a.path().to_str().unwrap(),
            "analyze",
        ];
        // Parse using `try_parse_from`
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::new(&args4.species).unwrap();
        let (result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        assert!(
            result.iter().any(|f| !f.status.status.isvalid()),
            "One locus is not valid"
        );
        let progress = getprogressbarclassic(result.len().try_into().unwrap_or_default()).unwrap();
        let mut i = 0;
        for loci in result {
            let genelist = extractgenelist(&args4, &loci, false).unwrap();
            for mut gene in genelist {
                let v = generategeneinfos(&args4, &mut gene).unwrap().0;
                assert!(v.status.status.isvalid(), "Gene {} is not valid", gene.gene);
            }
            i += 1;
            progress.set_position(i);
        }
        progress.with_finish(indicatif::ProgressFinish::AndClear);
    }
}

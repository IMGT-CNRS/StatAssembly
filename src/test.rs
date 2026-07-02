#[cfg(test)]
mod tests {
    use std::{
        io::{BufReader, ErrorKind::UnexpectedEof, Read as read2, Write},
        thread::sleep,
        time::Duration,
    };

    use crate::{
        extractgenelist, generategeneinfos,
        identification::{sendresult, sendresultcompressed},
        locusposparser,
        r#struct::{Args, Filecrea, GenesList, Species, SpeciesError},
        submissions::{
            BORNESLINK, RELEASELINK, REQUESTCLIENT, checkifblastpresent, generatelightbam,
            generatesequence, getprogressbarclassic,
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
    #[ignore]
    fn validfolder() {
        let testo = TempDir::new().unwrap();
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "full_cs.bam",
            "-a",
            "assembly.fasta",
            "-l",
            "locus.csv",
            "-s",
            "human",
            "-g",
            "geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "analyze",
        ];
        let human = Species::new("Homo sapiens").unwrap();
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let (locus, ..) = locusposparser(&args4, &human, false).unwrap();
        let mut append = std::fs::OpenOptions::new()
            .append(true)
            .open("~guilhem/historic.txt")
            .unwrap();
        let string = locus.iter().fold(String::new(), |mut acc, loci| {
            acc.push_str(&format!(
                "Locus {} is invalid at {}\n",
                loci.gethaplotype(),
                loci.contig
            ));
            acc
        });
        let s = string.trim();
        let _ = append.lock(); //Would fail on Windows if append so reject the error
        append.write_all(s.as_bytes()).unwrap();
        let _ = append.unlock();
        for loci in locus {
            let mut genes = extractgenelist(&args4, &loci, false).unwrap();
            let invalidgenes = genes
                .iter_mut()
                .map(|mut a| generategeneinfos(&args4, &mut a).unwrap())
                .filter(|p| !p.0.getstatus().getstatus().isvalid());
            let string = invalidgenes.fold(String::new(), |mut acc, (gene, _)| {
                acc.push_str(&format!(
                    "{} is invalid at {}:{}-{}\n",
                    gene.getgene(),
                    gene.getchromosome(),
                    gene.getstart().getobasedpos(),
                    gene.getend().getobasedpos()
                ));
                acc
            });
            let s = string.trim();
            let _ = append.lock(); //Would fail on Windows if append so reject the error
            append.write_all(s.as_bytes()).unwrap();
            let _ = append.unlock();
        }
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
        let seq = generatesequence(&args4, testo.path(), &locus).unwrap();
        let mut bufone =
            BufReader::new(std::fs::File::open("example_files/results/sequence.fasta.gz").unwrap());
        let mut buftwo = BufReader::new(std::fs::File::open(seq.getpath()).unwrap());
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
        assert_eq!(human.ischecked(), true, "Species could not be validated");
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
        assert_eq!(dog.ischecked(), true, "Species could not be validated");
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
        let (result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        assert!(
            result.iter().all(|f| f.status.status.isvalid()),
            "One locus is not valid"
        );
        generatelightbam(
            &args4,
            args4.outlightbam.as_ref().unwrap(),
            Some(b.getpath()),
            &result,
        )
        .unwrap();
        let mut br = bam::Reader::from_path(args4.file.unwrap()).unwrap();
        let mut cr = bam::Reader::from_path(&args4.outlightbam.unwrap().as_path()).unwrap();
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

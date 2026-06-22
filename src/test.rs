#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{
        extractgenelist, genelist, generategeneinfos, locusposparser,
        r#struct::{Args, Command, Species},
        submissions::{checkifblastpresent, generatelightbam},
    };
    use clap::Parser;
    use extended_htslib::bam::{self, Read};
    use indicatif::ProgressBar;
    use itertools::Itertools;
    use tempfile::{NamedTempFile, TempDir, env};
    #[test]
    fn checkblast() {
        assert!(checkifblastpresent())
    }
    #[test]
    fn testlocuspositiononlyandbam() {
        let testo = TempDir::new().unwrap();
        let a = NamedTempFile::with_suffix("test.bam").unwrap();
        let b = NamedTempFile::with_suffix("test.bam.csi").unwrap();
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
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::new(&args4.species).unwrap();
        let (result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        assert!(
            result.iter().all(|f| !f.status.status.isvalid()),
            "One locus is not valid"
        );
        generatelightbam(
            &args4,
            args4.outlightbam.as_ref().unwrap(),
            Some(b.path()),
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
        println!("Species is {}", spec);
        let (result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        assert!(
            result.iter().any(|f| !f.status.status.isvalid()),
            "One locus is not valid"
        );
        let progress = ProgressBar::new(result.len().try_into().unwrap_or_default());
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

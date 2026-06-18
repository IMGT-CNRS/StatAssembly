#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{
        genelist, locusposparser,
        r#struct::{Args, Command, Species},
        submissions::checkifblastpresent,
    };
    use clap::Parser;
    use tempfile::env;
    #[test]
    fn checkblast() {
        assert!(checkifblastpresent())
    }
    #[test]
    fn testlocuspositionandgenes() {
        let testo = Path::join(&env::temp_dir(), "test");
        let a = Path::join(&env::temp_dir(), "test.bam");
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
            testo.to_str().unwrap(),
            "-z",
            a.to_str().unwrap(),
            "analyze",
        ];
        let test4 = fake_args.iter().fold(String::new(), |mut acc, f| {
            acc.push_str(&format!(" {}", f));
            acc
        });
        let test4 = test4.trim();
        // Parse using `try_parse_from`
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        println!("Args are {:?}", args4);
        let spec = Species::new(&args4.species).unwrap();
        println!("Species is {}", spec);
        let (result, _blast) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        assert!(
            result.iter().any(|f| !f.status.status.isvalid()),
            "One locus is not valid"
        );
        for loci in result {
            let v = genelist(&loci, &spec, &args4, false).unwrap();
            assert!(
                v.iter().all(|p| p.status.status.isvalid()),
                "One gene is not valid"
            );
        }
    }
}

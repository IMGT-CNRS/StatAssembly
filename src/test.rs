/*
* IMGT/StatAssembly
Available under EUPL license
Made by: Guilhem Zeitoun and IMGT Team
*/
#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use std::{
        collections::BTreeMap,
        fs,
        io::{BufReader, Cursor, ErrorKind::UnexpectedEof, Read as read2, Seek},
        path::{Path, PathBuf},
        thread::sleep,
        time::Duration,
    };

    use crate::{
        extractgenelist, filtervalidatedallelesrecreate, generategeneinfos, getorsetparams,
        getreaderoffile,
        identification::{sendresult, sendresultcompressed},
        locusposparser, posread, printvalidatedallelesraw,
        r#struct::{
            Args, Filecrea, Genehit, Genename, GenesList, Haplotype::Primary, HashMapinfo,
            LocusHaplo, LocusInfos, Positionstrand, Species, SpeciesError,
        },
        submissions::{
            BORNESLINK, RELEASELINK, REQUESTCLIENT, checkifblastpresent, decompressseq,
            generatelightbam, generatesequence, generatesequenceraw, getprogressbarclassic,
        },
    };
    use clap::Parser;
    use extended_htslib::bam::{
        self, FetchDefinition, Read,
        ext::BamRecordExtensions,
        pileup::{RustPileupConfig, RustPileups},
    };
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
    fn checkparams() {
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
            "--extractedlength",
            "4772468",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            "analyze",
        ];
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let params = getorsetparams(&testo.path().join("test.txt"), &args4).unwrap();
        assert_eq!(params.getavg(), 14009, "Average read not the same");
        assert_eq!(params.getmean(), 56, "Mean coverage is not the same");
        assert_eq!(params.getphred(), 78, "Phred score is not the same");
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
        let human = Species::newunchecked("Homo sapiens");
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
            "Species should be valid rank"
        );
        assert!(
            human.getname().eq_ignore_ascii_case("Homo sapiens"),
            "Homo sapiens should be valid"
        );
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let dog = Species::new("dog").unwrap();
        assert!(dog.ischecked(), "Species could not be validated");
        assert_eq!(dog.getid(), Some(9615), "Species taxon is invalid");
        assert!(
            dog.getrank().eq_ignore_ascii_case("subspecies"),
            "Subspecies should be a valid rank"
        );
        assert!(
            dog.getname().eq_ignore_ascii_case("Canis lupus familiaris"),
            "Dog should be valid"
        );
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let lamprey = Species::new("least brook lamprey");
        assert!(
            matches!(lamprey, Err(SpeciesError::Blocked)),
            "Lamprey should be blocked"
        );
        sleep(Duration::new(3, 0)); //Sleep for NCBI
        let pongo = Species::new("pongo");
        let expected = "The term used is not a species or an subspecies.".to_string();
        if let Err(e) = &pongo
            && e == &SpeciesError::Invalid(expected)
        {
            ()
        } else {
            panic!("Pongo should be blocked. Got {:?}", pongo.as_ref());
        }
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
        for line in val.deserialize().take(5000) {
            data.push(line.unwrap());
        }
        assert_eq!(data.iter().nth(3069).unwrap().qual, Some(80), "Issue");
        assert!(data.iter().nth(1022).unwrap().qual.is_none(), "Issue");
    }
    #[test]
    fn testbam() {
        let testo = TempDir::new().unwrap();
        let a = Filecrea::createtemp(None, Some("test.bam")).unwrap();
        let b = Filecrea::createtemp(None, Some("test.bam.csi")).unwrap();
        // Fake args as a Vec<&str>
        let output = format!("-z={}", a.getpath().to_str().unwrap());
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "--totalread",
            "-l",
            "example_files/CHM13v2.0loc.csv",
            "-s",
            "human",
            "-g",
            "example_files/CHM13v2.0geneloc.csv",
            "-o",
            testo.path().to_str().unwrap(),
            &output,
            "analyze",
        ];
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::newunchecked(&args4.species);
        let (result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
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
            .rc_records()
            .filter_map(Result::ok)
            .zip_eq(cr.rc_records().filter_map(Result::ok))
        {
            activ = true;
            assert_eq!(
                record.reference_start(),
                record2.reference_start(),
                "Initial record {} is different because {}",
                String::from_utf8_lossy(record.qname()),
                String::from_utf8_lossy(record2.qname())
            );
        }
        if !activ {
            panic!("No records found");
        }
    }
    #[test]
    fn testlocusposition() {
        let testo = TempDir::new().unwrap();
        // Fake args as a Vec<&str>
        let fake_args = vec![
            "IMGT_StatAssembly",
            "-f",
            "example_files/CHM13v2.0.bam",
            "-a",
            "example_files/assembly.fasta",
            "--totalread",
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
        let spec = Species::newunchecked(&args4.species);
        let (mut result, _blast, _release) = match locusposparser(&args4, &spec, true) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        let meanpath = PathBuf::from("example_files/results").join(".mean");
        assert!(meanpath.exists(), "No .mean in example_files");
        let mean = getorsetparams(&meanpath, &args4).unwrap();
        let mut done = false;
        result.iter_mut().for_each(|p| {
            let pos = posread(&args4, &p).unwrap();
            if p.getlocushaplo() == &LocusHaplo::new(crate::r#struct::Locus::IGK, Primary) {
                let mut elem = pos.iter().skip(1675043).take(11).map(|(_, b)| b);
                let first = elem.next().unwrap();
                let file = Filecrea::createtemp(None::<PathBuf>, None::<PathBuf>).unwrap();
                let mut csv = csv::WriterBuilder::new()
                    .comment(Some(b'#'))
                    .has_headers(false)
                    .delimiter(b'\t')
                    .flexible(false)
                    .from_path(file.getpath())
                    .unwrap();
                for record in elem {
                    csv.serialize(record).unwrap();
                }
                csv.flush().unwrap();
                let text = fs::read_to_string(file.getpath()).unwrap();
                assert_eq!(
                    text.to_string().trim(),
                    "269534	90531414	49	0	0	0.0013	49	34	0	83	0	0	NC	0.0	0.0
269533	90531415	49	0	0	0.0013	49	34	0	83	0	0	NC	0.0	0.0
269532	90531416	49	0	0	0.0013	49	34	0	83	0	0	NC	0.0	0.0
269531	90531417	49	0	0	0.0013	49	34	0	83	0	0	NC	0.0	0.0
269530	90531418	50	0	0	0.0013	49	34	0	84	0	0	NC	0.02	0.0
269529	90531419	50	0	0	0.0013	50	34	0	84	0	0	NC	0.0	0.0
269528	90531420	50	0	0	0.0013	50	34	0	84	0	0	NC	0.0	0.0
269527	90531421	50	0	0	0.0013	50	34	0	84	0	0	NC	0.0	0.0
269526	90531422	50	0	0	0.0013	50	34	0	84	0	0	NC	0.0	0.0
269525	90531423	50	0	0	0.0013	50	34	0	84	0	0	NC	0.0	0.0",
                    "Invalid posresult"
                );
                let mut reader = getreaderoffile(&args4).unwrap();
                reader
                    .fetch(FetchDefinition::RegionString(
                        p.getcontig().to_string().as_bytes(),
                        first.position.getzbasedpos(),
                        first.position.getobasedpos(),
                    ))
                    .unwrap();
                let (mut map60, mut map1, mut map0, mut total, mut sup, mut secondary) =
                    (0, 0, 0, 0, 0, 0);
                for record in reader.records().filter_map(Result::ok) {
                    match (
                        record.is_primary(),
                        record.is_secondary(),
                        record.is_supplementary(),
                        record.mapq(),
                    ) {
                        (true, .., 60) => map60 += 1,
                        (true, .., 1..=59) => map1 += 1,
                        (true, .., 0) => map0 += 1,
                        (false, true, ..) => secondary += 1,
                        (false, false, true, _) => sup += 1,
                        _ => (),
                    }
                    total += 1;
                }
                let pos = pos
                    .iter()
                    .find(|(_, value)| value.position == first.position)
                    .unwrap()
                    .1;
                assert_eq!(map60, pos.map60);
                assert_eq!(map1, pos.map1);
                assert_eq!(map0, pos.map0);
                assert_eq!(total, pos.total);
                assert_eq!(secondary, pos.secondary);
                assert_eq!(sup, pos.supplementary);
                done = true;
            }
            p.setstatus(mean.mean, &args4, &pos);
        });
        assert!(done, "Haplotype IGK not found");
        assert!(
            result.iter().all(|f| f.status.isvalidated()),
            "One locus is not valid"
        );
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
        let spec = Species::newunchecked(&args4.species);
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
        let output = format!("-z={}", a.path().to_str().unwrap());
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
            &output,
            "analyze",
        ];
        // Parse using `try_parse_from`
        let args4 = Args::try_parse_from(fake_args)
            .map_err(|f| f.to_string())
            .unwrap();
        let spec = Species::newunchecked(&args4.species);
        let (mut result, _blast, _release) = match locusposparser(&args4, &spec, false) {
            Err(e) => panic!("Error is {e}"),
            Ok(a) => a,
        };
        let progress = getprogressbarclassic(result.len().try_into().unwrap_or_default()).unwrap();
        let mut i = 0;
        let mut done = false;
        let mean = getorsetparams(Path::new("example_files/results/.mean"), &args4).unwrap();
        let mut locushash: BTreeMap<LocusInfos, Vec<Genehit>> = BTreeMap::new();
        for loci in result.iter_mut() {
            let info = posread(&args4, &loci).unwrap();
            loci.setstatus(mean.getmean(), &args4, &info);
            let genelist = extractgenelist(&args4, &loci, false).unwrap();
            locushash.insert(loci.clone(), Vec::new());
            for mut gene in genelist {
                let (v, hashmap) = generategeneinfos(&args4, &mut gene).unwrap();
                let entry = locushash.get_mut(&loci).unwrap();
                entry.push(Genehit::new(v.clone(), None));
                if gene.gene
                    == Genename::new("IGHV(III)-82*01", Some("V-REGION".to_string())).unwrap()
                {
                    let info = "IGHV(III)-82*01_V-REGION	NC_060938.1	Minus	101154845	101155136	292	47.242	35	34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(32=0X2D0I0S)-34(33=0X1D2I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=1X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D1I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=1X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(32=0X2D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(34=0X0D0I0S)-34(32=0X2D0I0S)-34(33=0X1D0I0S)-34(34=0X0D0I0S)-35(35=0X0D0I0S)-35(34=0X1D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(34=0X1D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(34=0X1D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(34=0X1D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(34=0X1D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)-35(35=0X0D0I0S)	34	24	24	24.0	292	Accepted";
                    let file = Filecrea::createtemp(None::<PathBuf>, None::<PathBuf>).unwrap();
                    let mut csv = csv::WriterBuilder::new()
                        .comment(Some(b'#'))
                        .has_headers(false)
                        .delimiter(b'\t')
                        .flexible(false)
                        .from_path(file.getpath())
                        .unwrap();
                    csv.serialize(&v).unwrap();
                    csv.flush().unwrap();
                    assert_eq!(fs::read_to_string(file.getpath()).unwrap().trim(), info);
                    let pileup = RustPileups::new(
                        "example_files/CHM13v2.0.bam",
                        None,
                        Some("example_files/assembly.fasta"),
                        None,
                        FetchDefinition::RegionString(
                            gene.chromosome.as_bytes(),
                            gene.getpositionstart().getobasedpos(),
                            gene.getpositionend().getzbasedpos(),
                        ),
                        RustPileupConfig::default(),
                    )
                    .unwrap();
                    for hit in pileup {
                        let line = hashmap
                            .get(&crate::Position::newfromoposition(
                                hit.getpos().try_into().unwrap(),
                            ))
                            .unwrap();
                        if hit.rbases.get(0).unwrap().contains("*")
                            || hit.rbases.get(0).unwrap().contains("#")
                        {
                            assert!(line.getindel() < line.gettotal());
                        } else {
                            assert!(
                                line.getindel() == line.gettotal(),
                                "Should equal {} but is {}",
                                line.gettotal(),
                                line.getindel()
                            );
                        }
                        if hit.rbases.get(0).unwrap().contains("+") {
                            assert!(line.getinsertion() > 0);
                        } else {
                            assert!(line.getinsertion() == 0);
                        }
                    }
                    //Check validatedalleles output
                    {
                        let entry = locushash
                            .iter()
                            .flat_map(|(l, p)| {
                                if let Some(b) =
                                    p.iter().find(|e| gene.getgene() == e.geneinfo.getgene())
                                {
                                    Some((l, b))
                                } else {
                                    None
                                }
                            })
                            .next()
                            .unwrap();
                        let mut fakehash = BTreeMap::new();
                        fakehash.insert(entry.0.clone(), vec![entry.1.clone()]);
                        let mut cursor = Cursor::new(Vec::new());
                        {
                            printvalidatedallelesraw(
                                &mut cursor,
                                None::<String>,
                                &args4,
                                &fakehash,
                            )
                            .unwrap();
                        }
                        cursor.rewind().unwrap();
                        let mut val = Vec::new();
                        cursor.read_to_end(&mut val).unwrap();
                        assert_eq!(
                            String::from_utf8_lossy(&val).trim(),
                            "name,position,subject,locus,strand,bestmatch,identity
\"#IMGT/GENE-DB release Unknown\"
IGHV(III)-82*01_V-REGION,101154845-101155136,NC_060938.1,IGH (Primary),REV,No hits,Unknown"
                        );
                    }
                    done = true;
                }
                assert!(v.status.isvalidated(), "Gene {} is not valid", gene.gene);
            }
            i += 1;
            progress.set_position(i);
        }
        assert!(
            result.iter().all(|f| f.status.isvalidated()),
            "One locus is not valid"
        );
        let new = filtervalidatedallelesrecreate(&locushash);
        assert_eq!(
            new.values().map(|p| p.len()).sum::<usize>(),
            locushash.values().map(|p| p.len()).sum::<usize>(),
            "Validated alleles issue"
        );
        progress.with_finish(indicatif::ProgressFinish::AndClear);
        assert!(done, "IGHV(III)-82 not analyzed");
        testo.close().unwrap()
    }
}

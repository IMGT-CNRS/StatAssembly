use extended_htslib::bam::{
    self,
    pileup::{RustPileupConfig, RustPileups},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    io::{self, ErrorKind::InvalidInput, Write},
    path::{Path, PathBuf},
};

use crate::{
    givename,
    r#struct::{Args, Filecrea, LocusInfos, Species},
    submissions::getprogressbarclassic,
};
pub(crate) fn printpileup(species: &Species, args: &Args, locus: &[LocusInfos]) -> io::Result<()> {
    println!("Generating pileups...");
    let progressbar = getprogressbarclassic(locus.len().try_into().unwrap_or_default());
    let path = args
        .file
        .as_ref()
        .ok_or_else(|| io::Error::new(InvalidInput, "No BAM given"))?;
    locus
        .par_iter()
        .map(|loc| {
            let mut file = Filecrea::createfrompath(Path::join(
                &args.outdir,
                PathBuf::from(
                    givename(
                        &species,
                        loc.getlocus(),
                        loc.getcontig(),
                        loc.gethaplotype().isprimary(),
                        "pileup.txt",
                        false,
                    )
                    .to_string(),
                ),
            ))
            .setfile()
            .unwrap();
            let d = RustPileups::new(
                path,
                args.index.as_ref(),
                args.assembly.as_ref(),
                args.assemblyindex.as_ref(),
                bam::FetchDefinition::RegionString(
                    loc.contig.as_bytes(),
                    loc.start.getzbasedpos(),
                    loc.end.getzbasedpos(),
                ),
                RustPileupConfig::default(),
            )
            .map_err(|d| io::Error::new(InvalidInput, d))?;
            file.write_all(d.to_string().as_bytes())?;
            file.flush()?;
            if let Ok(p) = &progressbar {
                p.set_length(p.length().unwrap_or_default().saturating_add(1));
            }
            Ok(())
        })
        .collect::<io::Result<()>>()?;
    if let Ok(p) = progressbar {
        p.with_finish(indicatif::ProgressFinish::AndLeave);
    }
    println!("Pileups finished");
    Ok(())
}

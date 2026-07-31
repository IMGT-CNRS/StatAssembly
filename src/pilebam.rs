use std::{
    io::{self, ErrorKind::InvalidInput},
    num::NonZero,
    thread::available_parallelism,
};

use crate::{
    r#struct::{Args, Filecrea, LocusInfos},
    submissions::{generatelightbam, generatesequence},
};
use bam as pileup;
pub(crate) fn pileup(loci: &LocusInfos, args: &Args) -> io::Result<()> {
    todo!("Change structure");
    let filename = Filecrea::createfrompath(
        args.outdir
            .join(format!("{}_pileup.txt", loci.getlocushaplo())),
    );
    let mut writer = filename.setfile()?;
    let light = Filecrea::createtemp(None, Some("bam.bam"))?;
    let lightindex = Filecrea::createtemp(None, Some("bam.bam.csi"))?;
    let _ = generatelightbam(
        args,
        light.getpath(),
        Some(lightindex.getpath()),
        &[loci.clone()],
    )
    .map_err(|d| io::Error::new(InvalidInput, d))?;
    let mut smallsequence = generatesequence(args, &std::env::temp_dir(), true, &[loci.clone()])
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    if let Some(num) = NonZero::new(1) {
        let mut opts = pileup::MpileupOpts::default();
        opts.min_baseq = 0;
        opts.no_baq = true;
        opts.output_all = 2;
        let _ = pileup::mpileup(
            light.getpath(),
            Some(smallsequence.getpath()),
            &mut writer,
            &opts,
            available_parallelism().unwrap_or(num),
        )
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    }
    smallsequence.flagtodelete();
    Ok(())
}

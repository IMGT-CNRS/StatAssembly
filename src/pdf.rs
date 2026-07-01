use std::io;
// main.rs
use std::io::Write;
use std::path::Path;
use std::{collections::HashMap, fs::File};

use itertools::Itertools;
use pdfrs::security::PdfPermissions;
use pdfrs::{elements, pdf_generator};
use serde::Serialize;

use crate::r#struct::{Blastmatch, Filecrea, GeneInfos, LocusInfos};
use crate::{NAME, VERSION, getgeneinfos};

/// Convertit un Vec<Vec<String>> en tableau Markdown
fn vec_to_markdown_table(data: &[Vec<String>]) -> String {
    if data.is_empty() {
        return String::new();
    }

    let mut table = String::new();

    // Ligne d'en-tête
    table.push_str("| ");
    table.push_str(&data[0].join(" | "));
    table.push_str(" |\n");

    // Séparateur (|---|---|...)
    table.push_str("|");
    for _ in 0..data[0].len() {
        table.push_str("---|");
    }
    table.push_str("\n");

    // Lignes de données
    for row in &data[1..] {
        table.push_str("| ");
        table.push_str(&row.join(" | "));
        table.push_str(" |\n");
    }

    table
}

pub(crate) fn generatepdf(
    infos: &[LocusInfos],
    hash: &HashMap<LocusInfos, Vec<Blastmatch>>,
    outdir: &Path,
) -> io::Result<Filecrea> {
    let locusdata: Vec<Vec<String>> = infos
        .into_iter()
        .map(|a| {
            let mut b = Vec::new();
            b.push(a.locusinfo.to_string());
            b.push(a.contig.to_string());
            b.push(a.start.getobasedpos().to_string());
            b.push(a.end.getobasedpos().to_string());
            b.push(a.complement.to_string());
            b.push(a.status.to_string());
            b
        })
        .collect();
    let hashdata: Vec<(&LocusInfos, Vec<GeneInfos>)> = hash
        .into_iter()
        .map(|(a, f)| (a, getgeneinfos(&f)))
        .collect();
    let hashdata: Vec<GeneInfos> = hashdata.into_iter().map(|(a, b)| b).flatten().collect();
    let elem: Vec<Vec<String>> = hashdata
        .into_iter()
        .map(|a| {
            let mut b = Vec::new();
            b.push(a.gene.to_string());
            b.push(a.chromosome.to_string());
            b.push(a.start.getobasedpos().to_string());
            b.push(a.end.getobasedpos().to_string());
            b.push(a.strand.to_string());
            b.push(a.status.to_string());
            b
        })
        .collect();
    let layout = format!(
        "# {} version {}\n![Logo software](images/logo_imgt.png)\n## List of locus {} ## List of validated alleles {}",
        NAME.as_str(),
        VERSION,
        vec_to_markdown_table(&locusdata),
        vec_to_markdown_table(&elem)
    );
    let elements = elements::parse_markdown(&layout);
    let layout = pdf_generator::PageLayout::portrait();
    let file = Filecrea::createfrompath(Path::join(&outdir, "summary.pdf"));
    pdf_generator::create_pdf_from_elements_with_layout(
        &file.getpath().display().to_string().as_str(),
        &elements,
        "Helvetica",
        12.0,
        layout,
    )
    .map_err(|d| io::Error::new(io::ErrorKind::InvalidInput, d))?;
    let permissions = PdfPermissions::read_only();
    Ok(file)
}

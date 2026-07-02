use std::io;
use std::io::ErrorKind::InvalidInput;
// main.rs
use std::io::Write;
use std::path::Path;
use std::{collections::HashMap, fs::File};

use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use oxidize_pdf::document::DocumentEncryption;
use oxidize_pdf::encryption::{OwnerPassword, UserPassword};
use oxidize_pdf::{Color, Document, Font, Image, Page, PageTables, TableOptions, encryption};
use serde::Serialize;

use crate::r#struct::{Blastmatch, Filecrea, GeneInfos, LocusInfos};
use crate::{AUTHOR, NAME, VERSION, getgeneinfos};

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
            b.push(a.start.getobasedpos().to_formatted_string(&Locale::en));
            b.push(a.end.getobasedpos().to_formatted_string(&Locale::en));
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
            b.push(a.start.getobasedpos().to_formatted_string(&Locale::en));
            b.push(a.end.getobasedpos().to_formatted_string(&Locale::en));
            b.push(a.strand.to_string());
            b.push(a.status.to_string());
            b
        })
        .collect();
    // Create a new document
    let mut doc = Document::new();
    doc.set_title("IMGT/StatAssembly summary");

    // Create a page
    let mut page = Page::a4();
    let black = Color::black();
    // Add text
    page.text()
        .set_font(Font::Helvetica, 24.0)
        .set_fill_color(black)
        .at(20.0, 700.0)
        .write(&format!("{} version {}", NAME.as_str(), VERSION))
        .map_err(|d| io::Error::new(InvalidInput, d))?;

    // Add graphics
    let softwareimage = Image::from_file("images/logo_software.png")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    page.add_image("Software logo", softwareimage);
    let image =
        Image::from_file("images/logo_imgt.png").map_err(|d| io::Error::new(InvalidInput, d))?;
    page.add_image("IMGT logo", image);
    page.text()
        .set_font(Font::Helvetica, 20.0)
        .set_fill_color(black)
        .at(20.0, 600.0)
        .write("Locus found")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    let tableoptions = TableOptions::default();
    page.add_quick_table(
        locusdata,
        5.0,
        400.0,
        page.content_width() - 10.0,
        Some(tableoptions),
    )
    .map_err(|d| io::Error::new(InvalidInput, d))?;
    doc.add_page(page);
    // Create a page
    let mut page = Page::a4();
    page.text()
        .set_font(Font::Helvetica, 20.0)
        .set_fill_color(black)
        .at(2.0, 700.0)
        .write("Validated alleles")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    let tableoptions = TableOptions::default();
    page.add_quick_table(
        elem,
        5.0,
        620.0,
        page.content_width() - 10.0,
        Some(tableoptions),
    )
    .map_err(|d| io::Error::new(InvalidInput, d))?;
    // Save the document
    doc.add_page(page);
    doc.set_producer(AUTHOR);
    doc.set_author(NAME.as_str());
    /* let mut permissions = encryption::Permissions::new();
    permissions.set_copy(true).set_print_high_quality(true);
    doc.set_encryption(DocumentEncryption {
        user_password: UserPassword(String::new()),
        owner_password: OwnerPassword(String::new()),
        permissions: permissions,
        strength: oxidize_pdf::document::EncryptionStrength::Aes256,
    }); */

    let file = Filecrea::createfrompath(Path::join(&outdir, "summary.pdf"));
    doc.save(file.getpath())
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    Ok(file)
}

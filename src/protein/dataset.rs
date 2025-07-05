use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use flate2::read::GzDecoder;
use reqwest::blocking::Client;

#[derive(Debug, Clone)]
pub struct SmallProtein {
    pub species: String,
    pub id: String,
    pub rna_seq: String,
    pub aa_seq: String,
    pub length: usize,
    pub chromosome: String,
    pub start: usize,
    pub stop: usize,
    pub strand: String,
    pub blocks: String,
    pub start_codon: String,
    pub phylo_csf_mean: f64,
}

pub fn download_and_parse_small_protein_dataset() -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    let url = "http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt2_LiteratureMining.txt.gz";
    let temp_file = "small_protein_dataset.txt.gz";
    let extracted_file = "small_protein_dataset.txt";

    if !Path::new(extracted_file).exists() {
        if !Path::new(temp_file).exists() {
            println!("Downloading small protein dataset...");
            let client = Client::new();
            let mut response = client.get(url).send()?;
            let mut file = File::create(temp_file)?;
            io::copy(&mut response, &mut file)?;
        }

        println!("Extracting small protein dataset...");
        let compressed_file = File::open(temp_file)?;
        let decoder = GzDecoder::new(compressed_file);
        let mut reader = BufReader::new(decoder);
        let mut extracted_content = String::new();
        reader.read_to_string(&mut extracted_content)?;

        std::fs::write(extracted_file, extracted_content)?;
    }

    let file = File::open(extracted_file)?;
    let reader = BufReader::new(file);
    let mut proteins = Vec::new();

    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            continue;
        }

        let protein = SmallProtein {
            species: fields[0].to_string(),
            id: fields[1].to_string(),
            rna_seq: fields[2].to_string(),
            aa_seq: fields[3].to_string(),
            length: fields[4].parse().unwrap_or(0),
            chromosome: fields[5].to_string(),
            start: fields[6].parse().unwrap_or(0),
            stop: fields[7].parse().unwrap_or(0),
            strand: fields[8].to_string(),
            blocks: fields[9].to_string(),
            start_codon: fields[10].to_string(),
            phylo_csf_mean: fields[11].parse().unwrap_or(0.0),
        };

        proteins.push(protein);
    }

    Ok(proteins)
}
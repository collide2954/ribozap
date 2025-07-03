use std::io::{self, BufRead, BufReader, Read};
use std::error::Error;
use std::path::Path;
use std::fs::File;
use reqwest::blocking::Client;
use flate2::read::GzDecoder;

#[derive(Debug, Clone)]
struct SmallProtein {
    species: String,
    id: String,
    rna_seq: String,
    aa_seq: String,
    length: usize,
    chromosome: String,
    start: usize,
    stop: usize,
    strand: String,
    blocks: String,
    start_codon: String,
    phylo_csf_mean: f64,
}

fn download_and_parse_small_protein_dataset() -> Result<Vec<SmallProtein>, Box<dyn Error>> {
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

    println!("Parsing small protein dataset...");
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

    println!("Loaded {} small proteins", proteins.len());
    Ok(proteins)
}

fn calculate_dna_similarity(seq1: &str, seq2: &str) -> f64 {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let min_len = seq1.len().min(seq2.len());
    if min_len == 0 {
        return 0.0;
    }

    let mut matches = 0;
    for (c1, c2) in seq1.chars().zip(seq2.chars()) {
        if c1 == c2 {
            matches += 1;
        }
    }

    (matches as f64 / min_len as f64) * 100.0
}

fn main() -> Result<(), Box<dyn Error>> {
    let proteins = download_and_parse_small_protein_dataset()?;

    println!("First few proteins:");
    for (i, protein) in proteins.iter().take(5).enumerate() {
        println!("{}. {} ({}): RNA length = {}, AA length = {}", 
            i + 1, 
            protein.id, 
            protein.species, 
            protein.rna_seq.len(), 
            protein.aa_seq.len()
        );
    }

    if !proteins.is_empty() {
        let test_seq = "ATGAAAAACCCCAGTTGGATTAGAAAGAACTGGCTTCTTGTGGCTGGGGTGACTTTCATAGGCGTCCATCTTGGAACATACTTTATACAGAGAGTTGCAAAAGAGTCTGTGAGGTCTGAGGCCAGAGGCAGACAAAAGAATATTGAAGAATGA";

        let mut best_match = None;
        let mut best_similarity = 0.0;

        for protein in &proteins {
            let similarity = calculate_dna_similarity(test_seq, &protein.rna_seq);
            if similarity > best_similarity {
                best_similarity = similarity;
                best_match = Some(protein);
            }
        }

        if let Some(protein) = best_match {
            println!("\nClosest protein to test sequence:");
            println!("ID: {}", protein.id);
            println!("Species: {}", protein.species);
            println!("Similarity: {:.2}%", best_similarity);
            println!("RNA Seq: {}", protein.rna_seq);
            println!("AA Seq: {}", protein.aa_seq);
        }
    }

    Ok(())
}

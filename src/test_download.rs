//! Test utility for downloading and parsing protein dataset

use std::error::Error;
use ribozap::protein::{download_and_parse_small_protein_dataset, calculate_dna_similarity};

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

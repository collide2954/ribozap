//! Protein analysis and matching module

pub mod dataset;
pub mod matching;
pub mod molecular_weights;

// Re-export commonly used types and functions
pub use dataset::*;
pub use matching::*;
pub use molecular_weights::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_protein_dataset_and_similarity() -> Result<(), Box<dyn std::error::Error>> {
        // Test downloading and parsing protein dataset
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

    #[test]
    fn test_dna_similarity_basic() {
        let seq1 = "ATCG";
        let seq2 = "ATCG";
        assert_eq!(calculate_dna_similarity(seq1, seq2), 100.0);

        let seq1 = "ATCG";
        let seq2 = "TTCG";
        assert_eq!(calculate_dna_similarity(seq1, seq2), 75.0);

        let seq1 = "ATCG";
        let seq2 = "TTTT";
        assert_eq!(calculate_dna_similarity(seq1, seq2), 25.0);

        let seq1 = "";
        let seq2 = "ATCG";
        assert_eq!(calculate_dna_similarity(seq1, seq2), 0.0);
    }
}
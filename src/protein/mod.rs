pub mod dataset;
pub mod matching;
pub mod molecular_weights;

pub use dataset::*;
pub use dataset::DatasetProgress;
pub use matching::*;
pub use molecular_weights::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_protein_dataset_and_similarity() -> Result<(), Box<dyn std::error::Error>> {
        let proteins = download_and_parse_small_protein_dataset()?;

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
            if let Some(_protein) = best_match {
            }
        }
        Ok(())
    }
}
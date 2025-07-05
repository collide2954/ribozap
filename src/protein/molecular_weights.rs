//! Molecular weights of amino acids
//!
//! This module contains the exact molecular weights of all 20 standard amino acids
//! in Daltons (Da). These weights are for the amino acid residues in a protein chain,
//! which means they exclude the weight of water that is removed during peptide bond formation.

use std::collections::HashMap;

/// Get the molecular weight of an amino acid by its single-letter code
pub fn get_amino_acid_molecular_weight(amino_acid: char) -> f64 {
    match amino_acid {
        'A' => 71.04,   // Alanine
        'R' => 156.10,  // Arginine
        'N' => 114.04,  // Asparagine
        'D' => 115.03,  // Aspartic acid
        'C' => 103.01,  // Cysteine
        'E' => 129.04,  // Glutamic acid
        'Q' => 128.06,  // Glutamine
        'G' => 57.02,   // Glycine
        'H' => 137.06,  // Histidine
        'I' => 113.08,  // Isoleucine
        'L' => 113.08,  // Leucine
        'K' => 128.09,  // Lysine
        'M' => 131.04,  // Methionine
        'F' => 147.07,  // Phenylalanine
        'P' => 97.05,   // Proline
        'S' => 87.03,   // Serine
        'T' => 101.05,  // Threonine
        'W' => 186.08,  // Tryptophan
        'Y' => 163.06,  // Tyrosine
        'V' => 99.07,   // Valine
        '*' => 0.0,     // Stop codon
        _ => 0.0,       // Unknown amino acid
    }
}

/// Get a HashMap of all amino acid molecular weights
pub fn get_all_molecular_weights() -> HashMap<char, f64> {
    let mut weights = HashMap::new();

    weights.insert('A', 71.04);   // Alanine
    weights.insert('R', 156.10);  // Arginine
    weights.insert('N', 114.04);  // Asparagine
    weights.insert('D', 115.03);  // Aspartic acid
    weights.insert('C', 103.01);  // Cysteine
    weights.insert('E', 129.04);  // Glutamic acid
    weights.insert('Q', 128.06);  // Glutamine
    weights.insert('G', 57.02);   // Glycine
    weights.insert('H', 137.06);  // Histidine
    weights.insert('I', 113.08);  // Isoleucine
    weights.insert('L', 113.08);  // Leucine
    weights.insert('K', 128.09);  // Lysine
    weights.insert('M', 131.04);  // Methionine
    weights.insert('F', 147.07);  // Phenylalanine
    weights.insert('P', 97.05);   // Proline
    weights.insert('S', 87.03);   // Serine
    weights.insert('T', 101.05);  // Threonine
    weights.insert('W', 186.08);  // Tryptophan
    weights.insert('Y', 163.06);  // Tyrosine
    weights.insert('V', 99.07);   // Valine

    weights
}

/// Calculate the exact molecular weight of a protein sequence
pub fn calculate_protein_molecular_weight(amino_acid_sequence: &str) -> f64 {
    let mut total_weight = 0.0;

    for amino_acid in amino_acid_sequence.chars() {
        total_weight += get_amino_acid_molecular_weight(amino_acid);
    }

    // Add the weight of water (18.015 Da) to account for the N-terminus and C-terminus
    // The protein has one additional H2O compared to the sum of residue weights
    total_weight + 18.015
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_individual_amino_acid_weights() {
        assert_eq!(get_amino_acid_molecular_weight('A'), 71.04);
        assert_eq!(get_amino_acid_molecular_weight('G'), 57.02);
        assert_eq!(get_amino_acid_molecular_weight('W'), 186.08);
        assert_eq!(get_amino_acid_molecular_weight('*'), 0.0);
        assert_eq!(get_amino_acid_molecular_weight('X'), 0.0);
    }

    #[test]
    fn test_protein_molecular_weight() {
        // Test with a simple dipeptide: Ala-Gly
        // Should be 71.04 + 57.02 + 18.015 = 146.075
        let weight = calculate_protein_molecular_weight("AG");
        assert!((weight - 146.075).abs() < 0.001);
    }

    #[test]
    fn test_empty_sequence() {
        assert_eq!(calculate_protein_molecular_weight(""), 18.015);
    }
}

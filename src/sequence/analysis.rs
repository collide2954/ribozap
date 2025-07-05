//! DNA/RNA sequence analysis functions

use crate::sequence::codon::codon_to_single_letter_amino_acid;
use crate::sequence::conversion::dna_sequence_to_mrna;
use crate::protein::molecular_weights::get_amino_acid_molecular_weight;

/// Calculate GC content as a percentage
pub fn calculate_gc_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let gc_count = dna.chars()
        .filter(|&c| c == 'G' || c == 'g' || c == 'C' || c == 'c')
        .count();

    (gc_count as f64 / dna.len() as f64) * 100.0
}

/// Calculate AT content as a percentage
pub fn calculate_at_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let at_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'T' || c == 't')
        .count();

    (at_count as f64 / dna.len() as f64) * 100.0
}

/// Calculate purine content (A, G) as a percentage
pub fn calculate_purine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let purine_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'G' || c == 'g')
        .count();

    (purine_count as f64 / dna.len() as f64) * 100.0
}

/// Calculate pyrimidine content (C, T) as a percentage
pub fn calculate_pyrimidine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let pyrimidine_count = dna.chars()
        .filter(|&c| c == 'C' || c == 'c' || c == 'T' || c == 't')
        .count();

    (pyrimidine_count as f64 / dna.len() as f64) * 100.0
}

/// Calculate amino acid length from DNA sequence
pub fn calculate_amino_acid_length(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }
    dna.len() / 3
}

/// Calculate precise molecular weight based on actual amino acid composition
pub fn estimate_molecular_weight(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    // Convert DNA to mRNA
    let mrna = dna_sequence_to_mrna(dna);
    let mut total_weight = 0.0;

    // Process each codon
    for i in (0..mrna.len()).step_by(3) {
        if i + 2 < mrna.len() {
            let codon = &mrna[i..i+3];
            let amino_acid = codon_to_single_letter_amino_acid(codon);

            // Skip stop codons in the weight calculation
            if amino_acid != '*' {
                total_weight += get_amino_acid_molecular_weight(amino_acid);
            }
        }
    }

    // Add the weight of water for the N-terminus and C-terminus
    total_weight + 18.015
}

/// Calculate hydrophobicity index
pub fn calculate_hydrophobicity_index(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    let mut hydrophobic_count = 0;
    let mut total_amino_acids = 0;

    for i in (0..dna.len()).step_by(3) {
        if i + 2 < dna.len() {
            let codon = &dna[i..i+3];
            let amino_acid = super::codon::codon_to_amino_acid(codon);
            
            total_amino_acids += 1;

            // Hydrophobic amino acids
            if matches!(amino_acid, "Phe" | "Leu" | "Ile" | "Met" | "Val" | "Ala" | "Trp") {
                hydrophobic_count += 1;
            }
        }
    }

    if total_amino_acids > 0 {
        (hydrophobic_count as f64 / total_amino_acids as f64) * 100.0
    } else {
        0.0
    }
}

/// Count charged residues (positive and negative)
pub fn count_charged_residues(dna: &str) -> (usize, usize) {
    if dna.len() < 3 {
        return (0, 0);
    }

    let mut positive_count = 0;
    let mut negative_count = 0;

    for i in (0..dna.len()).step_by(3) {
        if i + 2 < dna.len() {
            let codon = &dna[i..i+3];
            let amino_acid = super::codon::codon_to_amino_acid(codon);
            
            match amino_acid {
                "Lys" | "Arg" => positive_count += 1,
                "Asp" | "Glu" => negative_count += 1,
                _ => {}
            }
        }
    }

    (positive_count, negative_count)
}

/// Count open reading frames (ORFs)
pub fn count_orfs(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }

    let mut orf_count = 0;
    let mut in_orf = false;

    for i in (0..dna.len()).step_by(3) {
        if i + 2 < dna.len() {
            let codon = &dna[i..i+3].to_uppercase();
            
            if codon == "ATG" && !in_orf {
                in_orf = true;
            } else if (codon == "TAA" || codon == "TAG" || codon == "TGA") && in_orf {
                orf_count += 1;
                in_orf = false;
            }
        }
    }

    orf_count
}
use bio_seq::prelude::*;
use bio_seq::translation::{TranslationTable, STANDARD};
use crate::sequence::conversion::dna_sequence_to_mrna;
use crate::protein::molecular_weights::get_amino_acid_molecular_weight;

pub fn calculate_gc_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let gc_count = dna.chars()
        .filter(|&c| c == 'G' || c == 'g' || c == 'C' || c == 'c')
        .count();

    (gc_count as f64 / dna.len() as f64) * 100.0
}

pub fn calculate_at_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let at_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'T' || c == 't')
        .count();

    (at_count as f64 / dna.len() as f64) * 100.0
}

pub fn calculate_purine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let purine_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'G' || c == 'g')
        .count();

    (purine_count as f64 / dna.len() as f64) * 100.0
}

pub fn calculate_pyrimidine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let pyrimidine_count = dna.chars()
        .filter(|&c| c == 'C' || c == 'c' || c == 'T' || c == 't')
        .count();

    (pyrimidine_count as f64 / dna.len() as f64) * 100.0
}

pub fn calculate_amino_acid_length(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }
    dna.len() / 3
}

pub fn estimate_molecular_weight(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    if let Ok(translation) = crate::sequence::translation::translate_dna_to_amino(dna) {
        let mut total_weight = 0.0;

        for amino_char in translation.chars() {
            if amino_char != '*' {
                total_weight += get_amino_acid_molecular_weight(amino_char);
            }
        }

        total_weight + 18.015
    } else {
        let mrna = dna_sequence_to_mrna(dna);
        let mut total_weight = 0.0;

        for i in (0..mrna.len()).step_by(3) {
            if i + 2 < mrna.len() {
                let codon = &mrna[i..i+3];
                if let Ok(codon_seq) = codon.parse::<Seq<Dna>>() {
                    if codon_seq.len() == 3 {
                        let amino = STANDARD.to_amino(&codon_seq);
                        let amino_str = amino.to_string();
                        if let Some(amino_char) = amino_str.chars().next() {
                            if amino_char != '*' {
                                total_weight += get_amino_acid_molecular_weight(amino_char);
                            }
                        }
                    }
                }
            }
        }

        total_weight + 18.015
    }
}

pub fn calculate_hydrophobicity_index(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    if let Ok(translation) = crate::sequence::translation::translate_dna_to_amino(dna) {
        let mut hydrophobic_count = 0;
        let total_amino_acids = translation.len();

        for amino_char in translation.chars() {
            if matches!(amino_char, 'F' | 'L' | 'I' | 'M' | 'V' | 'A' | 'W') {
                hydrophobic_count += 1;
            }
        }

        if total_amino_acids > 0 {
            (hydrophobic_count as f64 / total_amino_acids as f64) * 100.0
        } else {
            0.0
        }
    } else {
        0.0
    }
}

pub fn count_charged_residues(dna: &str) -> (usize, usize) {
    if dna.len() < 3 {
        return (0, 0);
    }

    if let Ok(translation) = crate::sequence::translation::translate_dna_to_amino(dna) {
        let mut positive_count = 0;
        let mut negative_count = 0;

        for amino_char in translation.chars() {
            match amino_char {
                'K' | 'R' => positive_count += 1,
                'D' | 'E' => negative_count += 1,
                _ => {}
            }
        }

        (positive_count, negative_count)
    } else {
        (0, 0)
    }
}

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
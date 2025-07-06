use bio_seq::prelude::*;
use bio_seq::translation::{TranslationTable, STANDARD};

pub fn translate_dna_to_amino(dna: &str) -> Result<String, String> {
    if dna.len() % 3 != 0 {
        return Err("DNA sequence length must be divisible by 3".to_string());
    }

    if dna.chars().all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'T' | 'G' | 'C')) {
        if let Ok(seq) = dna.parse::<Seq<Dna>>() {
            let mut amino_acids = Vec::new();

            for codon_chunk in seq.chunks(3) {
                if codon_chunk.len() == 3 {
                    let amino = STANDARD.to_amino(codon_chunk);
                    amino_acids.push(amino.to_string());
                }
            }

            return Ok(amino_acids.join(""));
        }
    }

    Err("Invalid DNA sequence".to_string())
}

pub fn translate_all_reading_frames(dna: &str) -> Result<Vec<String>, String> {
    let mut translations = Vec::new();

    for offset in 0..3 {
        if offset < dna.len() {
            let frame_dna = &dna[offset..];
            if let Ok(translation) = translate_dna_to_amino(frame_dna) {
                translations.push(translation);
            }
        }
    }

    let revcomp = crate::sequence::conversion::get_reverse_complement(dna);
    for offset in 0..3 {
        if offset < revcomp.len() {
            let frame_dna = &revcomp[offset..];
            if let Ok(translation) = translate_dna_to_amino(frame_dna) {
                translations.push(translation);
            }
        }
    }

    Ok(translations)
}

pub fn find_longest_orf(dna: &str) -> Result<(String, usize, usize), String> {
    let mut longest_orf = String::new();
    let mut longest_start = 0;
    let mut longest_end = 0;

    let frames = find_reading_frames_simple(dna);

    for frame_dna in &frames {
        let mut current_orf = String::new();
        let mut in_orf = false;
        let mut orf_start = 0;

        for i in (0..frame_dna.len()).step_by(3) {
            if i + 2 < frame_dna.len() {
                let codon = &frame_dna[i..i+3].to_uppercase();

                if codon == "ATG" && !in_orf {
                    in_orf = true;
                    orf_start = i;
                    current_orf = "M".to_string();
                } else if matches!(codon.as_str(), "TAA" | "TAG" | "TGA") && in_orf {
                    if current_orf.len() > longest_orf.len() {
                        longest_orf = current_orf.clone();
                        longest_start = orf_start;
                        longest_end = i + 3;
                    }
                    current_orf.clear();
                    in_orf = false;
                } else if in_orf {
                    if let Ok(codon_seq) = codon.parse::<Seq<Dna>>() {
                        if codon_seq.len() == 3 {
                            let amino = STANDARD.to_amino(&codon_seq);
                            let amino_str = amino.to_string();
                            if amino_str != "*" {
                                current_orf.push_str(&amino_str);
                            }
                        }
                    }
                }
            }
        }
    }

    Ok((longest_orf, longest_start, longest_end))
}

pub fn calculate_codon_usage(dna: &str) -> Result<std::collections::HashMap<String, usize>, String> {
    let mut codon_counts = std::collections::HashMap::new();

    for i in (0..dna.len()).step_by(3) {
        if i + 2 < dna.len() {
            let codon = dna[i..i+3].to_uppercase();
            *codon_counts.entry(codon).or_insert(0) += 1;
        }
    }

    Ok(codon_counts)
}

fn find_reading_frames_simple(dna: &str) -> Vec<String> {
    let mut frames = Vec::new();

    for offset in 0..3 {
        if offset < dna.len() {
            frames.push(dna[offset..].to_string());
        }
    }

    let revcomp = crate::sequence::conversion::get_reverse_complement(dna);
    for offset in 0..3 {
        if offset < revcomp.len() {
            frames.push(revcomp[offset..].to_string());
        }
    }

    frames
}

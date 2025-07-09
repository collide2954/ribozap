use bio_seq::prelude::*;
use bio_seq::translation::{TranslationTable, STANDARD};

pub fn dna_codon_to_amino_acid(codon: &str) -> String {
    if let Ok(codon_seq) = codon.parse::<Seq<Dna>>() {
        if codon_seq.len() == 3 {
            STANDARD.to_amino(&codon_seq).to_string()
        } else {
            "?".to_string()
        }
    } else {
        "?".to_string()
    }
}

pub fn count_total_codons(dna: &str) -> usize {
    dna.len() / 3
}

pub fn count_complete_incomplete_codons(dna: &str) -> (usize, usize) {
    let complete = dna.len() / 3;
    let incomplete = if dna.len() % 3 > 0 { 1 } else { 0 };
    (complete, incomplete)
}

pub fn count_start_codons(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }

    let mut count = 0;
    for i in 0..=(dna.len() - 3) {
        let codon = &dna[i..i+3].to_uppercase();
        if codon == "ATG" || codon == "AUG" {
            count += 1;
        }
    }
    count
}

pub fn count_stop_codons(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }

    let mut count = 0;
    for i in 0..=(dna.len() - 3) {
        let codon = &dna[i..i+3].to_uppercase();
        if codon == "TAA" || codon == "TAG" || codon == "TGA" {
            count += 1;
        }
    }
    count
}

pub fn find_reading_frames(dna: &str) -> Vec<String> {
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

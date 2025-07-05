//! Codon-related functions and analysis

/// Convert a codon to its corresponding amino acid
pub fn codon_to_amino_acid(codon: &str) -> &'static str {
    match codon {
        "UUU" | "UUC" => "Phe",
        "UUA" | "UUG" | "CUU" | "CUC" | "CUA" | "CUG" => "Leu",
        "AUU" | "AUC" | "AUA" => "Ile",
        "AUG" => "Met",
        "GUU" | "GUC" | "GUA" | "GUG" => "Val",
        "UCU" | "UCC" | "UCA" | "UCG" | "AGU" | "AGC" => "Ser",
        "CCU" | "CCC" | "CCA" | "CCG" => "Pro",
        "ACU" | "ACC" | "ACA" | "ACG" => "Thr",
        "GCU" | "GCC" | "GCA" | "GCG" => "Ala",
        "UAU" | "UAC" => "Tyr",
        "CAU" | "CAC" => "His",
        "CAA" | "CAG" => "Gln",
        "AAU" | "AAC" => "Asn",
        "AAA" | "AAG" => "Lys",
        "GAU" | "GAC" => "Asp",
        "GAA" | "GAG" => "Glu",
        "UGU" | "UGC" => "Cys",
        "UGG" => "Trp",
        "CGU" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => "Arg",
        "GGU" | "GGC" | "GGA" | "GGG" => "Gly",
        "UAA" | "UAG" | "UGA" => "Stop",
        _ => "???",
    }
}

/// Count total possible codons in a DNA sequence
pub fn count_total_codons(dna: &str) -> usize {
    dna.len() / 3
}

/// Count complete and incomplete codons
pub fn count_complete_incomplete_codons(dna: &str) -> (usize, usize) {
    let complete = dna.len() / 3;
    let incomplete = if dna.len() % 3 > 0 { 1 } else { 0 };
    (complete, incomplete)
}

/// Count start codons (ATG/AUG) in the sequence
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

/// Count stop codons in the sequence
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

/// Convert three-letter amino acid code to single-letter code
pub fn three_letter_to_single_letter(three_letter: &str) -> char {
    match three_letter {
        "Ala" => 'A',
        "Arg" => 'R',
        "Asn" => 'N',
        "Asp" => 'D',
        "Cys" => 'C',
        "Glu" => 'E',
        "Gln" => 'Q',
        "Gly" => 'G',
        "His" => 'H',
        "Ile" => 'I',
        "Leu" => 'L',
        "Lys" => 'K',
        "Met" => 'M',
        "Phe" => 'F',
        "Pro" => 'P',
        "Ser" => 'S',
        "Thr" => 'T',
        "Trp" => 'W',
        "Tyr" => 'Y',
        "Val" => 'V',
        "Stop" => '*',
        _ => '?',
    }
}

/// Convert a codon to its corresponding single-letter amino acid code
pub fn codon_to_single_letter_amino_acid(codon: &str) -> char {
    let three_letter = codon_to_amino_acid(codon);
    three_letter_to_single_letter(three_letter)
}

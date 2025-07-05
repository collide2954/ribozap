//! Base conversion functions for DNA/RNA sequences

/// Convert a DNA base to its complementary base
pub fn get_complementary_base(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

/// Convert a DNA base to its corresponding mRNA base
pub fn dna_to_mrna(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'U',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

/// Convert a complete DNA sequence to mRNA sequence
pub fn dna_sequence_to_mrna(dna: &str) -> String {
    dna.chars()
        .map(|base| dna_to_mrna(base))
        .collect()
}

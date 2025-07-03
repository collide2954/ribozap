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

/// Convert an mRNA base to its corresponding tRNA base
pub fn mrna_to_trna(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'U',
        'U' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}
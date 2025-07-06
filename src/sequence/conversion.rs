pub fn get_complementary_base(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

pub fn dna_to_mrna(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'U',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

pub fn dna_sequence_to_mrna(dna: &str) -> String {
    dna.chars()
        .map(dna_to_mrna)
        .collect()
}

pub fn get_reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(get_complementary_base)
        .collect()
}

pub fn get_complement(dna: &str) -> String {
    dna.chars()
        .map(get_complementary_base)
        .collect()
}

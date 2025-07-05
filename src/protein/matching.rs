pub fn calculate_dna_similarity(seq1: &str, seq2: &str) -> f64 {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let min_len = seq1.len().min(seq2.len());
    if min_len == 0 {
        return 0.0;
    }

    let mut matches = 0;
    for (c1, c2) in seq1.chars().zip(seq2.chars()) {
        if c1 == c2 {
            matches += 1;
        }
    }

    (matches as f64 / min_len as f64) * 100.0
}

pub fn identify_matching_positions(seq1: &str, seq2: &str) -> Vec<bool> {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let mut matches = Vec::new();
    for (c1, c2) in seq1.chars().zip(seq2.chars()) {
        matches.push(c1 == c2);
    }

    matches
}
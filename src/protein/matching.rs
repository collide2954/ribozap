pub fn calculate_dna_similarity(seq1: &str, seq2: &str) -> f64 {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let min_len = seq1.len().min(seq2.len());
    if min_len == 0 {
        return 0.0;
    }

    let matches = seq1.chars()
        .zip(seq2.chars())
        .take(min_len)
        .filter(|(a, b)| a == b)
        .count();

    (matches as f64 / min_len as f64) * 100.0
}

pub fn identify_matching_positions(seq1: &str, seq2: &str) -> Vec<bool> {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    seq1.chars()
        .zip(seq2.chars())
        .map(|(a, b)| a == b)
        .collect()
}

pub fn calculate_kmer_similarity<const K: usize>(seq1: &str, seq2: &str) -> f64 {
    if seq1.len() < K || seq2.len() < K {
        return 0.0;
    }

    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let kmers1: std::collections::HashSet<_> = (0..=seq1.len() - K)
        .map(|i| &seq1[i..i + K])
        .collect();

    let kmers2: std::collections::HashSet<_> = (0..=seq2.len() - K)
        .map(|i| &seq2[i..i + K])
        .collect();

    let intersection = kmers1.intersection(&kmers2).count();
    let union = kmers1.union(&kmers2).count();

    if union == 0 {
        0.0
    } else {
        (intersection as f64 / union as f64) * 100.0
    }
}

pub fn find_longest_common_subsequence(seq1: &str, seq2: &str) -> String {
    let seq1: Vec<char> = seq1.to_uppercase().chars().collect();
    let seq2: Vec<char> = seq2.to_uppercase().chars().collect();

    let len1 = seq1.len();
    let len2 = seq2.len();
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];

    for i in 1..=len1 {
        for j in 1..=len2 {
            if seq1[i-1] == seq2[j-1] {
                dp[i][j] = dp[i-1][j-1] + 1;
            } else {
                dp[i][j] = dp[i-1][j].max(dp[i][j-1]);
            }
        }
    }

    let mut lcs = Vec::new();
    let mut i = len1;
    let mut j = len2;

    while i > 0 && j > 0 {
        if seq1[i-1] == seq2[j-1] {
            lcs.push(seq1[i-1]);
            i -= 1;
            j -= 1;
        } else if dp[i-1][j] > dp[i][j-1] {
            i -= 1;
        } else {
            j -= 1;
        }
    }

    lcs.reverse();
    lcs.iter().collect()
}

pub fn calculate_amino_acid_similarity(seq1: &str, seq2: &str) -> f64 {
    let min_len = seq1.len().min(seq2.len());
    if min_len == 0 {
        return 0.0;
    }

    let matches = seq1.chars()
        .zip(seq2.chars())
        .take(min_len)
        .filter(|(a, b)| a == b)
        .count();

    (matches as f64 / min_len as f64) * 100.0
}
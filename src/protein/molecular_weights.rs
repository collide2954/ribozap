use std::collections::HashMap;

pub fn get_amino_acid_molecular_weight(amino_acid: char) -> f64 {
    match amino_acid {
        'A' => 71.04,
        'R' => 156.10,
        'N' => 114.04,
        'D' => 115.03,
        'C' => 103.01,
        'E' => 129.04,
        'Q' => 128.06,
        'G' => 57.02,
        'H' => 137.06,
        'I' => 113.08,
        'L' => 113.08,
        'K' => 128.09,
        'M' => 131.04,
        'F' => 147.07,
        'P' => 97.05,
        'S' => 87.03,
        'T' => 101.05,
        'W' => 186.08,
        'Y' => 163.06,
        'V' => 99.07,
        '*' => 0.0,
        _ => 0.0,
    }
}

pub fn get_all_molecular_weights() -> HashMap<char, f64> {
    let mut weights = HashMap::new();
    weights.insert('A', 71.04);
    weights.insert('R', 156.10);
    weights.insert('N', 114.04);
    weights.insert('D', 115.03);
    weights.insert('C', 103.01);
    weights.insert('E', 129.04);
    weights.insert('Q', 128.06);
    weights.insert('G', 57.02);
    weights.insert('H', 137.06);
    weights.insert('I', 113.08);
    weights.insert('L', 113.08);
    weights.insert('K', 128.09);
    weights.insert('M', 131.04);
    weights.insert('F', 147.07);
    weights.insert('P', 97.05);
    weights.insert('S', 87.03);
    weights.insert('T', 101.05);
    weights.insert('W', 186.08);
    weights.insert('Y', 163.06);
    weights.insert('V', 99.07);
    weights.insert('*', 0.0);
    weights
}

pub fn calculate_protein_molecular_weight(amino_acid_sequence: &str) -> f64 {
    let mut total_weight = 0.0;

    for amino_acid in amino_acid_sequence.chars() {
        total_weight += get_amino_acid_molecular_weight(amino_acid);
    }

    total_weight + 18.015
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_individual_amino_acid_weights() {
        assert_eq!(get_amino_acid_molecular_weight('A'), 71.04);
        assert_eq!(get_amino_acid_molecular_weight('G'), 57.02);
        assert_eq!(get_amino_acid_molecular_weight('W'), 186.08);
        assert_eq!(get_amino_acid_molecular_weight('*'), 0.0);
        assert_eq!(get_amino_acid_molecular_weight('X'), 0.0);
    }

    #[test]
    fn test_protein_molecular_weight() {
        let weight = calculate_protein_molecular_weight("AG");
        assert!((weight - 146.075).abs() < 0.001);
    }

    #[test]
    fn test_empty_sequence() {
        assert_eq!(calculate_protein_molecular_weight(""), 18.015);
    }
}

//! Main application state and logic

use ratatui::style::Color;
use crate::protein::{SmallProtein, download_and_parse_small_protein_dataset, calculate_dna_similarity, identify_matching_positions};
use crate::sequence::{get_complementary_base, dna_to_mrna, mrna_to_trna, codon_to_amino_acid};
use crate::ui::colors::get_amino_acid_color;

/// Main application state
pub struct App {
    pub input: String,
    pub complementary: String,
    pub mrna: String,
    pub trna: String,
    pub amino_acids: String,
    pub amino_acids_colored: Vec<(String, Color)>,
    pub current_codon_position: usize,
    pub small_proteins: Vec<SmallProtein>,
    pub closest_protein: Option<SmallProtein>,
    pub is_loading_proteins: bool,
    pub loading_error: Option<String>,
    pub loaded_proteins_count: usize,
    pub is_positive_strand: bool,
    pub matching_positions: Vec<bool>,
    pub current_strand_confidence: f64,
    pub opposite_strand_confidence: f64,
    pub last_input_length: usize,
    pub protein_match_needed: bool,
}

impl App {
    /// Create a new App instance
    pub fn new() -> App {
        let mut app = App {
            input: String::new(),
            complementary: String::new(),
            mrna: String::new(),
            trna: String::new(),
            amino_acids: String::new(),
            amino_acids_colored: Vec::new(),
            current_codon_position: 0,
            small_proteins: Vec::new(),
            closest_protein: None,
            is_loading_proteins: false,
            loading_error: None,
            loaded_proteins_count: 0,
            is_positive_strand: true,
            matching_positions: Vec::new(),
            current_strand_confidence: 0.0,
            opposite_strand_confidence: 0.0,
            last_input_length: 0,
            protein_match_needed: false,
        };

        match download_and_parse_small_protein_dataset() {
            Ok(proteins) => {
                app.loaded_proteins_count = proteins.len();
                app.small_proteins = proteins;
            },
            Err(e) => {
                app.loading_error = Some(format!("Error loading proteins: {}", e));
            }
        }

        app
    }

    /// Find the closest protein match
    pub fn find_closest_protein(&mut self) {
        if self.input.is_empty() || self.small_proteins.is_empty() {
            self.closest_protein = None;
            self.matching_positions.clear();
            self.current_strand_confidence = 0.0;
            self.opposite_strand_confidence = 0.0;
            return;
        }

        let mut best_match = None;
        let mut best_similarity = 0.0;
        let mut best_matching_positions = Vec::new();

        let mut positive_strand_similarities = Vec::new();
        let mut negative_strand_similarities = Vec::new();

        for protein in &self.small_proteins {
            let positive_similarity = calculate_dna_similarity(&self.input, &protein.rna_seq);
            let negative_similarity = calculate_dna_similarity(&self.complementary, &protein.rna_seq);

            // Collect ALL similarities for both strands regardless of protein strand annotation
            positive_strand_similarities.push(positive_similarity);
            negative_strand_similarities.push(negative_similarity);

            let (compare_seq, protein_seq, similarity) = if protein.strand == "-" {
                (&self.complementary, &protein.rna_seq, negative_similarity)
            } else {
                (&self.input, &protein.rna_seq, positive_similarity)
            };

            if similarity > best_similarity {
                best_similarity = similarity;
                best_match = Some(protein.clone());
                best_matching_positions = identify_matching_positions(compare_seq, protein_seq);
            }
        }

        // Always calculate confidence for current vs opposite strand
        self.current_strand_confidence = self.calculate_strand_confidence(&positive_strand_similarities);
        self.opposite_strand_confidence = self.calculate_strand_confidence(&negative_strand_similarities);

        self.closest_protein = best_match;
        self.matching_positions = best_matching_positions;
    }

    /// Calculate strand confidence based on similarities
    pub fn calculate_strand_confidence(&self, similarities: &[f64]) -> f64 {
        if similarities.is_empty() {
            return 0.0;
        }

        let mut sorted_similarities = similarities.to_vec();
        sorted_similarities.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

        let top_n = sorted_similarities.len().min(5);
        if top_n == 0 {
            return 0.0;
        }

        let sum: f64 = sorted_similarities.iter().take(top_n).sum();
        sum / top_n as f64
    }

    /// Toggle between positive and negative strand modes
    pub fn toggle_strand_mode(&mut self) {
        std::mem::swap(&mut self.input, &mut self.complementary);
        self.is_positive_strand = !self.is_positive_strand;
        self.update_sequences();
        self.find_closest_protein();
        self.protein_match_needed = false;
    }

    /// Update all derived sequences based on current input
    pub fn update_sequences(&mut self) {
        if self.is_positive_strand {
            self.complementary = self.input
                .chars()
                .map(get_complementary_base)
                .collect();

            self.mrna = self.input
                .chars()
                .map(dna_to_mrna)
                .collect();
        } else {
            let positive_strand: String = self.complementary
                .chars()
                .map(get_complementary_base)
                .collect();

            self.input = positive_strand;

            self.mrna = self.input
                .chars()
                .map(dna_to_mrna)
                .collect();
        }

        self.trna = self.mrna
            .chars()
            .map(mrna_to_trna)
            .collect();

        self.current_codon_position = self.mrna.len() % 3;

        self.update_amino_acids();

        let current_length = self.input.len();
        if current_length < 10 || 
           self.last_input_length == 0 || 
           current_length.abs_diff(self.last_input_length) >= 3 {
            self.protein_match_needed = true;
        }
        self.last_input_length = current_length;
    }

    /// Perform protein matching if needed
    pub fn perform_protein_matching_if_needed(&mut self) {
        if self.protein_match_needed {
            self.find_closest_protein();
            self.protein_match_needed = false;
        }
    }

    /// Get the current partial codon for completion display
    pub fn get_current_partial_codon(&self) -> String {
        if self.mrna.is_empty() {
            return String::new();
        }

        let mrna_str = self.mrna.to_uppercase();
        let codon_start = (self.mrna.len() / 3) * 3;

        if codon_start >= mrna_str.len() {
            return String::new();
        }

        mrna_str[codon_start..].to_string()
    }

    /// Update amino acid sequence with colors
    fn update_amino_acids(&mut self) {
        self.amino_acids = String::new();
        self.amino_acids_colored.clear();

        let mrna_str = self.mrna.to_uppercase();
        let mut i = 0;

        while i + 2 < mrna_str.len() {
            let codon = &mrna_str[i..i+3];
            let amino = codon_to_amino_acid(codon);
            let color = get_amino_acid_color(amino);

            if !self.amino_acids.is_empty() {
                self.amino_acids.push(' ');
            }
            self.amino_acids.push_str(amino);

            self.amino_acids_colored.push((amino.to_string(), color));

            i += 3;
        }

        if i < mrna_str.len() {
            if !self.amino_acids.is_empty() {
                self.amino_acids.push(' ');
            }
            self.amino_acids.push_str("???");

            self.amino_acids_colored.push(("???".to_string(), Color::White));
        }
    }

    /// Handle key input
    pub fn on_key(&mut self, c: char) {
        if self.is_positive_strand {
            self.input.push(c);
        } else {
            self.complementary.push(c);
        }
        self.update_sequences();
    }

    /// Handle backspace
    pub fn on_backspace(&mut self) {
        if self.is_positive_strand {
            self.input.pop();
        } else {
            self.complementary.pop();
        }
        self.update_sequences();
    }
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}
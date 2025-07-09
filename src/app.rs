use ratatui::style::Color;
use bio_seq::prelude::*;
use bio_seq::translation::{TranslationTable, STANDARD};
use crate::protein::{SmallProtein, calculate_dna_similarity, identify_matching_positions, DatasetProgress};
use crate::sequence::{get_complementary_base, dna_to_mrna};
use std::collections::HashMap;
use std::sync::mpsc::{self, Receiver};
use std::thread;
use log::{info, warn, error, debug, trace};



pub struct App {
    pub input: String,
    pub complementary: String,
    pub mrna: String,
    pub amino_acids: String,
    pub amino_acids_colored: Vec<(String, Color)>,
    pub current_codon_position: usize,
    pub small_proteins: Vec<SmallProtein>,
    pub closest_protein: Option<SmallProtein>,
    pub is_loading_proteins: bool,
    pub loading_error: Option<String>,
    pub loaded_proteins_count: usize,
    pub dataset_progress: Option<DatasetProgress>,
    pub is_positive_strand: bool,
    pub matching_positions: Vec<bool>,
    pub current_strand_confidence: f64,
    pub opposite_strand_confidence: f64,
    pub last_input_length: usize,
    pub protein_match_needed: bool,
    pub show_protein_searcher: bool,
    pub searcher_input: String,
    pub searcher_field: SearchField,
    pub filtered_proteins: Vec<SmallProtein>,
    pub selected_protein_index: usize,
    pub selected_search_field: usize,
    pub search_filters: HashMap<SearchField, String>,
    pub multi_search_mode: bool,
    pub show_protein_detail: bool,
    pub detailed_protein: Option<SmallProtein>,
    pub progress_receiver: Option<Receiver<DatasetProgress>>,
    pub protein_receiver: Option<Receiver<Result<Vec<SmallProtein>, String>>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SearchField {
    Species,
    Id,
    Chromosome,
    Strand,
    StartCodon,
    MinLength,
    MaxLength,
    MinPhyloCSF,
    MaxPhyloCSF,
}

impl App {
    pub fn new() -> App {
        debug!("Creating new App instance");
        App {
            input: String::new(),
            complementary: String::new(),
            mrna: String::new(),
            amino_acids: String::new(),
            amino_acids_colored: Vec::new(),
            current_codon_position: 0,
            small_proteins: Vec::new(),
            closest_protein: None,
            is_loading_proteins: true,
            loading_error: None,
            loaded_proteins_count: 0,
            dataset_progress: Some(DatasetProgress::CheckingCache),
            is_positive_strand: true,
            matching_positions: Vec::new(),
            current_strand_confidence: 0.0,
            opposite_strand_confidence: 0.0,
            last_input_length: 0,
            protein_match_needed: false,
            show_protein_searcher: false,
            searcher_input: String::new(),
            searcher_field: SearchField::Species,
            filtered_proteins: Vec::new(),
            selected_protein_index: 0,
            selected_search_field: 0,
            search_filters: HashMap::new(),
            multi_search_mode: false,
            show_protein_detail: false,
            detailed_protein: None,
            progress_receiver: None,
            protein_receiver: None,
        }
    }

    pub fn load_datasets(&mut self) {
        use crate::protein::download_and_parse_small_protein_dataset_with_progress;
        
        info!("Starting dataset loading");
        self.is_loading_proteins = true;
        self.loading_error = None;
        self.dataset_progress = Some(DatasetProgress::CheckingCache);
        
        match download_and_parse_small_protein_dataset_with_progress(None) {
            Ok(proteins) => {
                self.loaded_proteins_count = proteins.len();
                self.small_proteins = proteins;
                self.is_loading_proteins = false;
                self.dataset_progress = Some(DatasetProgress::Complete);
                info!("Successfully loaded {} proteins", self.loaded_proteins_count);
            },
            Err(e) => {
                error!("Failed to load protein datasets: {e}");
                self.loading_error = Some(format!("Error loading proteins: {e}"));
                self.is_loading_proteins = false;
                self.dataset_progress = Some(DatasetProgress::Error(e.to_string()));
            }
        }
    }

    pub fn update_progress(&mut self, progress: DatasetProgress) {
        self.dataset_progress = Some(progress);
    }

    fn check_and_load_dataset_file(&mut self) -> bool {
        use crate::protein::dataset::get_data_dir;
        if let Ok(data_dir) = get_data_dir() {
            let extracted_file = data_dir.join("small_protein_dataset.txt");
            if extracted_file.exists() {
                self.dataset_progress = Some(DatasetProgress::Parsing { lines_parsed: 0 });
                self.load_datasets();
                return true;
            }
        }
        self.dataset_progress = Some(DatasetProgress::CheckingCache);
        false
    }

    pub fn check_and_load_datasets(&mut self) -> bool {
        self.check_and_load_dataset_file()
    }

    pub fn start_loading_if_needed(&mut self) -> bool {
        if self.is_loading_proteins && self.small_proteins.is_empty() && self.loading_error.is_none() {
            return self.check_and_load_dataset_file();
        }
        false
    }

    pub fn find_closest_protein(&mut self) {
        if self.input.is_empty() || self.small_proteins.is_empty() {
            self.closest_protein = None;
            self.matching_positions.clear();
            self.current_strand_confidence = 0.0;
            self.opposite_strand_confidence = 0.0;
            debug!("No input or proteins available for matching");
            return;
        }

        trace!("Finding closest protein match for input length: {}", self.input.len());
        let mut best_match = None;
        let mut best_similarity = 0.0;
        let mut best_matching_positions = Vec::new();

        let mut positive_strand_similarities = Vec::new();
        let mut negative_strand_similarities = Vec::new();

        for protein in &self.small_proteins {
            let positive_similarity = calculate_dna_similarity(&self.input, &protein.rna_seq);
            let negative_similarity = calculate_dna_similarity(&self.complementary, &protein.rna_seq);

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

        self.current_strand_confidence = self.calculate_strand_confidence(&positive_strand_similarities);
        self.opposite_strand_confidence = self.calculate_strand_confidence(&negative_strand_similarities);

        if let Some(ref protein) = best_match {
            debug!("Found best protein match: {} (similarity: {:.4})", protein.id, best_similarity);
        }

        self.closest_protein = best_match;
        self.matching_positions = best_matching_positions;
    }

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

    pub fn toggle_strand_mode(&mut self) {
        info!("Toggling strand mode from {} to {}",
              if self.is_positive_strand { "positive" } else { "negative" },
              if self.is_positive_strand { "negative" } else { "positive" });

        std::mem::swap(&mut self.input, &mut self.complementary);
        self.is_positive_strand = !self.is_positive_strand;
        self.update_sequences();
        self.find_closest_protein();
        self.protein_match_needed = false;
    }

    pub fn update_sequences(&mut self) {
        trace!("Updating sequences for {} strand",
               if self.is_positive_strand { "positive" } else { "negative" });

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

        self.current_codon_position = self.mrna.len() % 3;

        self.update_amino_acids();

        let current_length = self.input.len();
        if current_length < 10 || 
           self.last_input_length == 0 ||
           current_length.abs_diff(self.last_input_length) >= 3 {
            self.protein_match_needed = true;
            debug!("Protein matching needed due to sequence length change: {} -> {}",
                   self.last_input_length, current_length);
        }
        self.last_input_length = current_length;
    }

    pub fn perform_protein_matching_if_needed(&mut self) {
        if self.protein_match_needed {
            self.find_closest_protein();
            self.protein_match_needed = false;
        }
    }

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

    fn update_amino_acids(&mut self) {
        self.amino_acids = String::new();
        self.amino_acids_colored.clear();

        let mrna_str = self.mrna.to_uppercase();
        let mut i = 0;

        while i + 2 < mrna_str.len() {
            let codon = &mrna_str[i..i+3];
            let dna_codon = codon.replace('U', "T");

            let amino = if let Ok(codon_seq) = dna_codon.parse::<Seq<Dna>>() {
                if codon_seq.len() == 3 {
                    STANDARD.to_amino(&codon_seq).to_string()
                } else {
                    "?".to_string()
                }
            } else {
                "?".to_string()
            };
            let color = Color::White;

            if !self.amino_acids.is_empty() {
                self.amino_acids.push(' ');
            }
            self.amino_acids.push_str(&amino);

            self.amino_acids_colored.push((amino, color));

            i += 3;
        }

        if i < mrna_str.len() {
            if !self.amino_acids.is_empty() {
                self.amino_acids.push(' ');
            }
            self.amino_acids.push('_');

            self.amino_acids_colored.push(("_".to_string(), Color::White));
        }
    }

    pub fn on_key(&mut self, c: char) {
        trace!("Adding nucleotide '{}' to {} strand", c,
               if self.is_positive_strand { "positive" } else { "negative" });

        if self.is_positive_strand {
            self.input.push(c);
        } else {
            self.complementary.push(c);
        }
        self.update_sequences();
    }

    pub fn on_backspace(&mut self) {
        trace!("Removing nucleotide from {} strand",
               if self.is_positive_strand { "positive" } else { "negative" });

        if self.is_positive_strand {
            self.input.pop();
        } else {
            self.complementary.pop();
        }
        self.update_sequences();
    }

    pub fn toggle_protein_searcher(&mut self) {
        self.show_protein_searcher = !self.show_protein_searcher;
        debug!("Protein searcher toggled: {}", self.show_protein_searcher);

        if self.show_protein_searcher {
            self.searcher_input.clear();
            self.selected_protein_index = 0;
            self.selected_search_field = 0;
            self.filter_proteins();
            info!("Protein searcher opened with {} total proteins", self.small_proteins.len());
        }
    }

    pub fn searcher_on_key(&mut self, c: char) {
        if self.show_protein_searcher {
            self.searcher_input.push(c);
            self.filter_proteins();
        }
    }

    pub fn searcher_on_backspace(&mut self) {
        if self.show_protein_searcher {
            self.searcher_input.pop();
            self.filter_proteins();
        }
    }

    pub fn searcher_next_field(&mut self) {
        if self.show_protein_searcher {
            self.selected_search_field = (self.selected_search_field + 1) % 9;
            self.update_search_field();
            self.searcher_input.clear();
            self.filter_proteins();
        }
    }

    pub fn searcher_prev_field(&mut self) {
        if self.show_protein_searcher {
            self.selected_search_field = if self.selected_search_field == 0 { 8 } else { self.selected_search_field - 1 };
            self.update_search_field();
            self.searcher_input.clear();
            self.filter_proteins();
        }
    }

    pub fn searcher_next_protein(&mut self) {
        if self.show_protein_searcher && !self.filtered_proteins.is_empty() {
            self.selected_protein_index = (self.selected_protein_index + 1) % self.filtered_proteins.len();
        }
    }

    pub fn searcher_prev_protein(&mut self) {
        if self.show_protein_searcher && !self.filtered_proteins.is_empty() {
            self.selected_protein_index = if self.selected_protein_index == 0 {
                self.filtered_proteins.len() - 1
            } else {
                self.selected_protein_index - 1
            };
        }
    }

    pub fn select_current_protein(&mut self) {
        if self.show_protein_searcher
            && !self.filtered_proteins.is_empty()
            && self.selected_protein_index < self.filtered_proteins.len() {
            self.detailed_protein = Some(self.filtered_proteins[self.selected_protein_index].clone());
            self.show_protein_detail = true;
        }
    }

    pub fn return_to_search(&mut self) {
        self.show_protein_detail = false;
        self.detailed_protein = None;
    }

    pub fn select_detailed_protein(&mut self) {
        if let Some(protein) = &self.detailed_protein {
            self.closest_protein = Some(protein.clone());
            self.show_protein_searcher = false;
            self.show_protein_detail = false;
            self.detailed_protein = None;
        }
    }

    fn update_search_field(&mut self) {
        self.searcher_field = match self.selected_search_field {
            0 => SearchField::Species,
            1 => SearchField::Id,
            2 => SearchField::Chromosome,
            3 => SearchField::Strand,
            4 => SearchField::StartCodon,
            5 => SearchField::MinLength,
            6 => SearchField::MaxLength,
            7 => SearchField::MinPhyloCSF,
            8 => SearchField::MaxPhyloCSF,
            _ => SearchField::Species,
        };
    }

    pub fn toggle_multi_search_mode(&mut self) {
        self.multi_search_mode = !self.multi_search_mode;
        if self.multi_search_mode {
            if !self.searcher_input.is_empty() {
                self.search_filters.insert(self.searcher_field, self.searcher_input.clone());
            }
            self.searcher_input.clear();
        } else {
            self.search_filters.clear();
        }
        self.filter_proteins();
    }

    pub fn add_current_filter(&mut self) {
        if !self.searcher_input.is_empty() {
            self.search_filters.insert(self.searcher_field, self.searcher_input.clone());
            self.searcher_input.clear();
            self.filter_proteins();
        }
    }

    pub fn clear_current_filter(&mut self) {
        self.search_filters.remove(&self.searcher_field);
        self.filter_proteins();
    }

    pub fn clear_all_filters(&mut self) {
        self.search_filters.clear();
        self.searcher_input.clear();
        self.filter_proteins();
    }

    pub fn get_active_filters(&self) -> Vec<(SearchField, String)> {
        self.search_filters.iter()
            .map(|(field, value)| (*field, value.clone()))
            .collect()
    }

    fn filter_proteins(&mut self) {
        let initial_count = self.small_proteins.len();

        if self.multi_search_mode {
            self.filtered_proteins = self.small_proteins.iter()
                .filter(|protein| {
                    for (field, value) in &self.search_filters {
                        if !self.matches_field_criteria(protein, *field, value) {
                            return false;
                        }
                    }
                    if !self.searcher_input.is_empty()
                        && !self.matches_field_criteria(protein, self.searcher_field, &self.searcher_input) {
                        return false;
                    }
                    true
                })
                .cloned()
                .collect();
        } else if self.searcher_input.is_empty() {
            self.filtered_proteins = self.small_proteins.clone();
        } else {
            self.filtered_proteins = self.small_proteins.iter()
                .filter(|protein| self.matches_field_criteria(protein, self.searcher_field, &self.searcher_input))
                .cloned()
                .collect();
        }

        if self.selected_protein_index >= self.filtered_proteins.len() {
            self.selected_protein_index = 0;
        }

        debug!("Protein search filtered from {} to {} results", initial_count, self.filtered_proteins.len());
    }

    fn matches_field_criteria(&self, protein: &SmallProtein, field: SearchField, value: &str) -> bool {
        let search_term = value.to_lowercase();

        match field {
            SearchField::Species => protein.species.to_lowercase().contains(&search_term),
            SearchField::Id => protein.id.to_lowercase().contains(&search_term),
            SearchField::Chromosome => protein.chromosome.to_lowercase().contains(&search_term),
            SearchField::Strand => protein.strand.to_lowercase().contains(&search_term),
            SearchField::StartCodon => protein.start_codon.to_lowercase().contains(&search_term),
            SearchField::MinLength => {
                if let Ok(min_length) = value.parse::<usize>() {
                    protein.length >= min_length
                } else {
                    true
                }
            },
            SearchField::MaxLength => {
                if let Ok(max_length) = value.parse::<usize>() {
                    protein.length <= max_length
                } else {
                    true
                }
            },
            SearchField::MinPhyloCSF => {
                if let Ok(min_phylo) = value.parse::<f64>() {
                    protein.phylo_csf_mean >= min_phylo
                } else {
                    true
                }
            },
            SearchField::MaxPhyloCSF => {
                if let Ok(max_phylo) = value.parse::<f64>() {
                    protein.phylo_csf_mean <= max_phylo
                } else {
                    true
                }
            },
        }
    }

    pub fn get_search_field_name(&self) -> &'static str {
        match self.searcher_field {
            SearchField::Species => "Species",
            SearchField::Id => "ID",
            SearchField::Chromosome => "Chromosome",
            SearchField::Strand => "Strand",
            SearchField::StartCodon => "Start Codon",
            SearchField::MinLength => "Min Length",
            SearchField::MaxLength => "Max Length",
            SearchField::MinPhyloCSF => "Min PhyloCSF",
            SearchField::MaxPhyloCSF => "Max PhyloCSF",
        }
    }

    pub fn simulate_loading_step(&mut self) -> bool {
        // Check if file already exists first
        use crate::protein::dataset::get_data_dir;
        if let Ok(data_dir) = get_data_dir() {
            let extracted_file = data_dir.join("small_protein_dataset.txt");
            if extracted_file.exists() {
                // File exists, load directly without showing progress
                self.load_datasets();
                return true;
            }
        }

        match &self.dataset_progress {
            Some(DatasetProgress::CheckingCache) => {
                // Move to downloading state
                self.dataset_progress = Some(DatasetProgress::Downloading { 
                    bytes_downloaded: 0, 
                    total_bytes: Some(15_000_000) // Simulate ~15MB file
                });
                false // Continue loading
            },
            Some(DatasetProgress::Downloading { bytes_downloaded, total_bytes }) => {
                let new_downloaded = bytes_downloaded + 2_000_000; // Simulate 2MB chunks
                if let Some(total) = total_bytes {
                    if new_downloaded >= *total {
                        // Download complete, move to extracting (restart progress)
                        self.dataset_progress = Some(DatasetProgress::Extracting);
                    } else {
                        // Continue downloading
                        self.dataset_progress = Some(DatasetProgress::Downloading {
                            bytes_downloaded: new_downloaded,
                            total_bytes: *total_bytes,
                        });
                    }
                }
                false // Continue loading
            },
            Some(DatasetProgress::Extracting) => {
                // Extraction complete, load the actual data
                self.load_datasets();
                true // Loading complete
            },
            _ => true // Already complete or error
        }
    }

    pub fn start_threaded_loading(&mut self) {
        if self.progress_receiver.is_some() || self.protein_receiver.is_some() {
            warn!("Attempted to start loading while already in progress");
            return;
        }

        info!("Starting threaded protein dataset loading");

        // Create channels for progress and result communication
        let (progress_tx, progress_rx) = mpsc::channel();
        let (result_tx, result_rx) = mpsc::channel();

        self.progress_receiver = Some(progress_rx);
        self.protein_receiver = Some(result_rx);

        // Spawn background thread for dataset loading
        thread::spawn(move || {
            use crate::protein::download_and_parse_small_protein_dataset_with_progress;

            debug!("Background thread started for dataset loading");

            // Create progress callback that sends updates through channel
            let progress_callback = Box::new(move |progress: DatasetProgress| {
                let _ = progress_tx.send(progress);
            });

            // Load dataset with progress callback
            let result = download_and_parse_small_protein_dataset_with_progress(Some(progress_callback));

            // Send final result
            let final_result = match result {
                Ok(proteins) => {
                    info!("Background loading completed successfully with {} proteins", proteins.len());
                    Ok(proteins)
                },
                Err(e) => {
                    error!("Background loading failed: {e}");
                    Err(e.to_string())
                }
            };

            let _ = result_tx.send(final_result);
            debug!("Background thread completed");
        });
    }

    pub fn check_loading_progress(&mut self) {
        // Check for progress updates
        if let Some(ref progress_rx) = self.progress_receiver {
            while let Ok(progress) = progress_rx.try_recv() {
                trace!("Loading progress update: {progress:?}");
                self.dataset_progress = Some(progress);
            }
        }

        // Check for final result
        if let Some(ref result_rx) = self.protein_receiver {
            if let Ok(result) = result_rx.try_recv() {
                match result {
                    Ok(proteins) => {
                        self.loaded_proteins_count = proteins.len();
                        self.small_proteins = proteins;
                        self.is_loading_proteins = false;
                        self.dataset_progress = Some(DatasetProgress::Complete);
                        info!("Loading completed successfully. {} proteins loaded", self.loaded_proteins_count);
                    },
                    Err(e) => {
                        error!("Loading failed with error: {e}");
                        self.loading_error = Some(format!("Error loading proteins: {e}"));
                        self.is_loading_proteins = false;
                        self.dataset_progress = Some(DatasetProgress::Error(e));
                    }
                }
                
                // Clear receivers as loading is complete
                self.progress_receiver = None;
                self.protein_receiver = None;
            }
        }
    }
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}
use std::io::{self, BufRead, BufReader, Read};
use std::error::Error;
use std::path::Path;
use std::fs::File;
use reqwest::blocking::Client;
use flate2::read::GzDecoder;
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout},
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph, BarChart},
    Terminal,
};

fn get_complementary_base(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

fn dna_to_mrna(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'U',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

fn mrna_to_trna(base: char) -> char {
    match base.to_uppercase().next().unwrap_or(' ') {
        'A' => 'U',
        'U' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => '?',
    }
}

fn codon_to_amino_acid(codon: &str) -> &'static str {
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

fn format_triplets(sequence: &str) -> String {
    let mut result = String::new();
    let mut count = 0;

    for c in sequence.chars() {
        result.push(c);
        count += 1;
        if count % 3 == 0 && count < sequence.len() {
            result.push(' ');
        }
    }

    result
}

fn create_codon_completion_display(partial_codon: &str) -> Vec<Line<'static>> {
    let mut lines = Vec::new();

    let mut current_codon_text = vec![
        Span::raw("Current codon: "),
    ];

    if partial_codon.is_empty() {
        current_codon_text.push(Span::styled("None", Style::default().fg(Color::DarkGray)));
    } else {
        current_codon_text.push(Span::styled(partial_codon.to_string(), Style::default().fg(Color::Green)));

        for _ in 0..(3 - partial_codon.len()) {
            current_codon_text.push(Span::styled("_", Style::default().fg(Color::DarkGray)));
        }
    }

    lines.push(Line::from(current_codon_text));

    lines.push(Line::from(vec![Span::raw("")]));

    let nucleotides = ['U', 'C', 'A', 'G'];

    match partial_codon.len() {
        0 => {
            lines.push(Line::from(vec![
                Span::styled("Start a new codon with any base:", Style::default().fg(Color::White)),
            ]));

            let mut first_options = Vec::new();
            for &base in &nucleotides {
                first_options.push(Span::styled(format!("{base} "), Style::default().fg(Color::Cyan)));
            }
            lines.push(Line::from(first_options));
        },
        1 => {
            let first_base = partial_codon.chars().next().unwrap();
            lines.push(Line::from(vec![
                Span::styled(format!("With first base {first_base}, add second base:"), Style::default().fg(Color::White)),
            ]));

            for &second_base in &nucleotides {
                let mut row = Vec::new();
                row.push(Span::styled(format!("{first_base}{second_base}_ → "), Style::default().fg(Color::Cyan)));

                let mut possible_aminos = Vec::new();
                for &third_base in &nucleotides {
                    let codon = format!("{first_base}{second_base}{third_base}");
                    let amino = codon_to_amino_acid(&codon);
                    possible_aminos.push(amino);
                }

                possible_aminos.sort();
                possible_aminos.dedup();
                let mut colored_amino_list = Vec::new();
                for (i, amino) in possible_aminos.iter().enumerate() {
                    if i > 0 {
                        colored_amino_list.push(Span::raw("/"));
                    }
                    colored_amino_list.push(Span::styled(*amino, Style::default().fg(get_amino_acid_color(amino))));
                }
                row.extend(colored_amino_list);

                lines.push(Line::from(row));
            }
        },
        2 => {
            let first_base = partial_codon.chars().next().unwrap();
            let second_base = partial_codon.chars().nth(1).unwrap();

            lines.push(Line::from(vec![
                Span::styled(format!("With bases {first_base}{second_base}, complete codon with:"), Style::default().fg(Color::White)),
            ]));

            for &third_base in &nucleotides {
                let codon = format!("{first_base}{second_base}{third_base}");
                let amino = codon_to_amino_acid(&codon);

                let color = get_amino_acid_color(amino);

                lines.push(Line::from(vec![
                    Span::styled(format!("{first_base}{second_base}{third_base} → "), Style::default().fg(Color::Cyan)),
                    Span::styled(amino, Style::default().fg(color)),
                ]));
            }
        },
        _ => {
            lines.push(Line::from(vec![
                Span::styled("Ready for next codon", Style::default().fg(Color::White)),
            ]));
        }
    }

    lines
}

fn calculate_gc_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let gc_count = dna.chars()
        .filter(|&c| c == 'G' || c == 'g' || c == 'C' || c == 'c')
        .count();

    (gc_count as f64 / dna.len() as f64) * 100.0
}

fn calculate_at_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let at_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'T' || c == 't')
        .count();

    (at_count as f64 / dna.len() as f64) * 100.0
}

fn calculate_purine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let purine_count = dna.chars()
        .filter(|&c| c == 'A' || c == 'a' || c == 'G' || c == 'g')
        .count();

    (purine_count as f64 / dna.len() as f64) * 100.0
}

fn calculate_pyrimidine_content(dna: &str) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }

    let pyrimidine_count = dna.chars()
        .filter(|&c| c == 'C' || c == 'c' || c == 'T' || c == 't')
        .count();

    (pyrimidine_count as f64 / dna.len() as f64) * 100.0
}

fn count_total_codons(dna: &str) -> usize {
    dna.len() / 3
}

fn count_complete_incomplete_codons(dna: &str) -> (usize, usize) {
    let complete = dna.len() / 3;
    let incomplete = if dna.len() % 3 > 0 { 1 } else { 0 };
    (complete, incomplete)
}

fn count_start_codons(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }

    let mut count = 0;
    for i in 0..=(dna.len() - 3) {
        let codon = &dna[i..i+3].to_uppercase();
        if codon == "ATG" {
            count += 1;
        }
    }
    count
}

fn count_stop_codons(dna: &str) -> usize {
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

fn calculate_amino_acid_length(dna: &str) -> usize {
    if dna.len() < 3 {
        return 0;
    }

    dna.len() / 3
}

fn estimate_molecular_weight(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    let amino_acid_count = dna.len() / 3;
    amino_acid_count as f64 * 110.0
}

fn calculate_hydrophobicity_index(dna: &str) -> f64 {
    if dna.len() < 3 {
        return 0.0;
    }

    let mut hydrophobic_count = 0;
    let mut total_amino_acids = 0;

    for i in (0..dna.len()).step_by(3) {
        if i + 3 <= dna.len() {
            let codon = &dna[i..i+3].to_uppercase();
            let amino = codon_to_amino_acid(codon);

            if amino == "Ala" || amino == "Ile" || amino == "Leu" || 
               amino == "Met" || amino == "Phe" || amino == "Val" || 
               amino == "Pro" || amino == "Gly" {
                hydrophobic_count += 1;
            }

            total_amino_acids += 1;
        }
    }

    if total_amino_acids == 0 {
        return 0.0;
    }

    (hydrophobic_count as f64 / total_amino_acids as f64) * 100.0
}

fn count_charged_residues(dna: &str) -> (usize, usize) {
    if dna.len() < 3 {
        return (0, 0);
    }

    let mut positive_count = 0;
    let mut negative_count = 0;

    for i in (0..dna.len()).step_by(3) {
        if i + 3 <= dna.len() {
            let codon = &dna[i..i+3].to_uppercase();
            let amino = codon_to_amino_acid(codon);

            if amino == "Arg" || amino == "Lys" || amino == "His" {
                positive_count += 1;
            }

            if amino == "Asp" || amino == "Glu" {
                negative_count += 1;
            }
        }
    }

    (positive_count, negative_count)
}

fn count_orfs(dna: &str) -> usize {
    if dna.len() < 6 {
        return 0;
    }

    let mut orf_count = 0;
    let mut in_orf = false;

    for i in (0..dna.len()-2).step_by(3) {
        let codon = &dna[i..i+3].to_uppercase();

        if codon == "ATG" && !in_orf {
            in_orf = true;
        } else if (codon == "TAA" || codon == "TAG" || codon == "TGA") && in_orf {
            orf_count += 1;
            in_orf = false;
        }
    }

    orf_count
}

fn get_amino_acid_color(amino: &str) -> Color {
    match amino {
        "Phe" => Color::Red,
        "Leu" => Color::Green,
        "Ile" => Color::Yellow,
        "Met" => Color::Blue,
        "Val" => Color::Magenta,
        "Ser" => Color::Cyan,
        "Pro" => Color::Gray,
        "Thr" => Color::DarkGray,
        "Ala" => Color::LightRed,
        "Tyr" => Color::LightGreen,
        "His" => Color::LightYellow,
        "Gln" => Color::LightBlue,
        "Asn" => Color::LightMagenta,
        "Lys" => Color::LightCyan,
        "Asp" => Color::White,
        "Glu" => Color::Red,
        "Cys" => Color::Green,
        "Trp" => Color::Yellow,
        "Arg" => Color::Blue,
        "Gly" => Color::Magenta,
        "Stop" => Color::Red,
        _ => Color::White,
    }
}

#[derive(Debug, Clone)]
struct SmallProtein {
    species: String,
    id: String,
    rna_seq: String,
    aa_seq: String,
    length: usize,
    chromosome: String,
    start: usize,
    stop: usize,
    strand: String,
    blocks: String,
    start_codon: String,
    phylo_csf_mean: f64,
}

fn download_and_parse_small_protein_dataset() -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    let url = "http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt2_LiteratureMining.txt.gz";
    let temp_file = "small_protein_dataset.txt.gz";
    let extracted_file = "small_protein_dataset.txt";

    if !Path::new(extracted_file).exists() {
        if !Path::new(temp_file).exists() {
            println!("Downloading small protein dataset...");
            let client = Client::new();
            let mut response = client.get(url).send()?;
            let mut file = File::create(temp_file)?;
            io::copy(&mut response, &mut file)?;
        }

        println!("Extracting small protein dataset...");
        let compressed_file = File::open(temp_file)?;
        let decoder = GzDecoder::new(compressed_file);
        let mut reader = BufReader::new(decoder);
        let mut extracted_content = String::new();
        reader.read_to_string(&mut extracted_content)?;

        std::fs::write(extracted_file, extracted_content)?;
    }

    let file = File::open(extracted_file)?;
    let reader = BufReader::new(file);
    let mut proteins = Vec::new();

    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            continue;
        }

        let protein = SmallProtein {
            species: fields[0].to_string(),
            id: fields[1].to_string(),
            rna_seq: fields[2].to_string(),
            aa_seq: fields[3].to_string(),
            length: fields[4].parse().unwrap_or(0),
            chromosome: fields[5].to_string(),
            start: fields[6].parse().unwrap_or(0),
            stop: fields[7].parse().unwrap_or(0),
            strand: fields[8].to_string(),
            blocks: fields[9].to_string(),
            start_codon: fields[10].to_string(),
            phylo_csf_mean: fields[11].parse().unwrap_or(0.0),
        };

        proteins.push(protein);
    }

    Ok(proteins)
}

fn calculate_dna_similarity(seq1: &str, seq2: &str) -> f64 {
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

fn identify_matching_positions(seq1: &str, seq2: &str) -> Vec<bool> {
    let seq1 = seq1.to_uppercase();
    let seq2 = seq2.to_uppercase();

    let mut matches = Vec::new();
    for (c1, c2) in seq1.chars().zip(seq2.chars()) {
        matches.push(c1 == c2);
    }

    matches
}

struct App {
    input: String,
    complementary: String,
    mrna: String,
    trna: String,
    amino_acids: String,
    amino_acids_colored: Vec<(String, Color)>,
    current_codon_position: usize,
    small_proteins: Vec<SmallProtein>,
    closest_protein: Option<SmallProtein>,
    is_loading_proteins: bool,
    loading_error: Option<String>,
    loaded_proteins_count: usize,
    is_positive_strand: bool,
    matching_positions: Vec<bool>,
    current_strand_confidence: f64,
    opposite_strand_confidence: f64,
    histogram_data: Vec<(SmallProtein, f64)>,
    max_histogram_entries: usize,
    last_input_length: usize,
    protein_match_needed: bool,
}

impl App {
    fn new() -> App {
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
            histogram_data: Vec::new(),
            max_histogram_entries: 10,
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

    fn find_closest_protein(&mut self) {
        if self.input.is_empty() || self.small_proteins.is_empty() {
            self.closest_protein = None;
            self.matching_positions.clear();
            self.current_strand_confidence = 0.0;
            self.opposite_strand_confidence = 0.0;
            self.histogram_data.clear();
            return;
        }

        let mut best_match = None;
        let mut best_similarity = 0.0;
        let mut best_matching_positions = Vec::new();

        let mut positive_strand_similarities = Vec::new();
        let mut negative_strand_similarities = Vec::new();

        self.histogram_data.clear();

        let mut all_similarities = Vec::new();

        for protein in &self.small_proteins {
            let positive_similarity = calculate_dna_similarity(&self.input, &protein.rna_seq);
            let negative_similarity = calculate_dna_similarity(&self.complementary, &protein.rna_seq);

            if protein.strand == "+" {
                positive_strand_similarities.push(positive_similarity);
            } else if protein.strand == "-" {
                negative_strand_similarities.push(negative_similarity);
            }
            let (compare_seq, protein_seq, similarity) = if protein.strand == "-" {
                (&self.complementary, &protein.rna_seq, negative_similarity)
            } else {
                (&self.input, &protein.rna_seq, positive_similarity)
            };

            all_similarities.push((protein.clone(), similarity));

            if similarity > best_similarity {
                best_similarity = similarity;
                best_match = Some(protein.clone());
                best_matching_positions = identify_matching_positions(compare_seq, protein_seq);
            }
        }

        all_similarities.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        self.histogram_data = all_similarities.into_iter()
            .take(self.max_histogram_entries)
            .collect();

        if self.is_positive_strand {
            self.current_strand_confidence = self.calculate_strand_confidence(&positive_strand_similarities);
            self.opposite_strand_confidence = self.calculate_strand_confidence(&negative_strand_similarities);
        } else {
            self.current_strand_confidence = self.calculate_strand_confidence(&negative_strand_similarities);
            self.opposite_strand_confidence = self.calculate_strand_confidence(&positive_strand_similarities);
        }

        self.closest_protein = best_match;
        self.matching_positions = best_matching_positions;
    }

    fn calculate_strand_confidence(&self, similarities: &[f64]) -> f64 {
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

    fn toggle_strand_mode(&mut self) {
        std::mem::swap(&mut self.input, &mut self.complementary);

        self.is_positive_strand = !self.is_positive_strand;

        self.update_sequences();

        self.find_closest_protein();
        self.protein_match_needed = false;
    }

    fn update_sequences(&mut self) {
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

    fn perform_protein_matching_if_needed(&mut self) {
        if self.protein_match_needed {
            self.find_closest_protein();
            self.protein_match_needed = false;
        }
    }

    fn get_current_partial_codon(&self) -> String {
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

    fn on_key(&mut self, c: char) {
        if self.is_positive_strand {
            self.input.push(c);
        } else {
            self.complementary.push(c);
        }
        self.update_sequences();
    }

    fn on_backspace(&mut self) {
        if self.is_positive_strand {
            self.input.pop();
        } else {
            self.complementary.pop();
        }
        self.update_sequences();
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let mut app = App::new();
    loop {
        app.perform_protein_matching_if_needed();

        terminal.draw(|f| {
            let main_horizontal_split = Layout::default()
                .direction(Direction::Horizontal)
                .margin(2)
                .constraints(
                    [
                        Constraint::Percentage(70),
                        Constraint::Percentage(30),
                    ]
                    .as_ref(),
                )
                .split(f.size());

            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints(
                    [
                        Constraint::Length(3),
                        Constraint::Length(5),
                        Constraint::Length(5),
                        Constraint::Length(5),
                        Constraint::Length(25),
                        Constraint::Min(3),
                    ]
                    .as_ref(),
                )
                .split(main_horizontal_split[0]);

            let title = "DNA to Protein Sequence Converter";

            let mut spans = Vec::new();

            spans.push(Span::styled(title, Style::default().fg(Color::Cyan)));

            spans.push(Span::raw("   "));

            spans.push(Span::styled(
                format!("Loaded {} small proteins", app.loaded_proteins_count),
                Style::default().fg(Color::Green)
            ));

            spans.push(Span::raw("   "));
            let strand_mode_spans = if app.is_positive_strand {
                vec![
                    Span::raw("Strand: "),
                    Span::styled("[+] Positive", Style::default().fg(Color::Green)),
                    Span::raw(" / "),
                    Span::styled("[-] Negative", Style::default().fg(Color::DarkGray)),
                ]
            } else {
                vec![
                    Span::raw("Strand: "),
                    Span::styled("[+] Positive", Style::default().fg(Color::DarkGray)),
                    Span::raw(" / "),
                    Span::styled("[-] Negative", Style::default().fg(Color::Yellow)),
                ]
            };

            spans.extend(strand_mode_spans);

            let title_text = Line::from(spans);

            let title = Paragraph::new(vec![title_text])
                .block(Block::default().borders(Borders::ALL));
            f.render_widget(title, chunks[0]);

            let formatted_input = format_triplets(&app.input);
            let input_text = Line::from(vec![
                Span::raw("Positive Strand: "),
                Span::styled(&formatted_input, Style::default().fg(Color::Green)),
            ]);
            let input = Paragraph::new(vec![input_text])
                .block(Block::default().borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(input, chunks[1]);

            let formatted_complementary = format_triplets(&app.complementary);
            let complementary_text = Line::from(vec![
                Span::raw("Negative Strand: "),
                Span::styled(&formatted_complementary, Style::default().fg(Color::Yellow)),
            ]);
            let complementary = Paragraph::new(vec![complementary_text])
                .block(Block::default().borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(complementary, chunks[2]);
            let formatted_mrna = format_triplets(&app.mrna);
            let mrna_text = Line::from(vec![
                Span::raw("mRNA:           "),
                Span::styled(&formatted_mrna, Style::default().fg(Color::Magenta)),
            ]);
            let mrna = Paragraph::new(vec![mrna_text])
                .block(Block::default().borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(mrna, chunks[3]);

            let amino_chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(30), Constraint::Percentage(30), Constraint::Percentage(40)].as_ref())
                .split(chunks[4]);
            let mut amino_spans = vec![Span::raw("Amino Acids: ")];

            for (i, (amino, color)) in app.amino_acids_colored.iter().enumerate() {
                if i > 0 {
                    amino_spans.push(Span::raw(" "));
                }

                amino_spans.push(Span::styled(amino, Style::default().fg(*color)));
            }

            let amino_text = Line::from(amino_spans);
            let amino = Paragraph::new(vec![amino_text])
                .block(Block::default().title("Amino Acid Sequence").borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(amino, amino_chunks[0]);

            let partial_codon = app.get_current_partial_codon();
            let codon_completion = create_codon_completion_display(&partial_codon);
            let codon_completion_widget = Paragraph::new(codon_completion)
                .block(Block::default().title("Codon Completion Guide").borders(Borders::ALL));
            f.render_widget(codon_completion_widget, amino_chunks[1]);
            let closest_protein_text = if app.is_loading_proteins {
                vec![Line::from(vec![
                    Span::styled("Loading protein database...", Style::default().fg(Color::Yellow)),
                ])]
            } else if let Some(error) = &app.loading_error {
                vec![Line::from(vec![
                    Span::styled(error, Style::default().fg(Color::Red)),
                ])]
            } else if let Some(protein) = &app.closest_protein {
                let rna_seq_clone = protein.rna_seq.clone();

                let mut rna_seq_spans = Vec::new();

                let mut triplet_count = 0;
                for (i, c) in rna_seq_clone.chars().enumerate() {
                    let is_match = i < app.matching_positions.len() && app.matching_positions[i];

                    let style = if is_match {
                        Style::default().fg(Color::Green).bg(Color::DarkGray)
                    } else {
                        Style::default().fg(Color::Cyan)
                    };

                    rna_seq_spans.push(Span::styled(c.to_string(), style));

                    triplet_count += 1;
                    if triplet_count == 3 && i < rna_seq_clone.len() - 1 {
                        rna_seq_spans.push(Span::raw(" "));
                        triplet_count = 0;
                    }
                }

                vec![
                    Line::from(vec![
                        Span::raw("Species: "),
                        Span::styled(&protein.species, Style::default().fg(Color::Green)),
                    ]),
                    Line::from(vec![
                        Span::raw("ID: "),
                        Span::styled(&protein.id, Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from(vec![
                        Span::raw("Length: "),
                        Span::styled(protein.length.to_string(), Style::default().fg(Color::Blue)),
                    ]),
                    Line::from(vec![
                        Span::raw("Chromosome: "),
                        Span::styled(&protein.chromosome, Style::default().fg(Color::Cyan)),
                    ]),
                    Line::from(vec![
                        Span::raw("Start: "),
                        Span::styled(protein.start.to_string(), Style::default().fg(Color::Green)),
                    ]),
                    Line::from(vec![
                        Span::raw("Stop: "),
                        Span::styled(protein.stop.to_string(), Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from(vec![
                        Span::raw("Strand: "),
                        Span::styled(&protein.strand, Style::default().fg(Color::Blue)),
                    ]),
                    Line::from(vec![
                        Span::raw("Blocks: "),
                        Span::styled(&protein.blocks, Style::default().fg(Color::Cyan)),
                    ]),
                    Line::from(vec![
                        Span::raw("Start Codon: "),
                        Span::styled(&protein.start_codon, Style::default().fg(Color::Green)),
                    ]),
                    Line::from(vec![
                        Span::raw("PhyloCSF Mean: "),
                        Span::styled(protein.phylo_csf_mean.to_string(), Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from({
                        let mut spans = vec![Span::raw("RNA Seq: ")];
                        spans.extend(rna_seq_spans);
                        spans
                    }),
                    Line::from(vec![
                        Span::raw("AA Seq: "),
                        Span::styled(&protein.aa_seq, Style::default().fg(Color::Magenta)),
                    ]),
                ]
            } else {
                vec![Line::from(vec![
                    Span::styled("No matching protein found", Style::default().fg(Color::DarkGray)),
                ])]
            };

            let title_spans = vec![Span::raw("Closest Small Protein Match")];
            let closest_protein_widget = Paragraph::new(closest_protein_text)
                .block(Block::default().title(Line::from(title_spans)).borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(closest_protein_widget, amino_chunks[2]);

            {
                let match_colors = [
                    Color::Red,
                    Color::Green,
                    Color::Yellow,
                    Color::Blue,
                    Color::Magenta,
                    Color::Cyan,
                    Color::LightRed,
                    Color::LightGreen,
                    Color::LightYellow,
                    Color::LightBlue,
                ];

                let top_matches_title = Line::from(vec![
                    Span::styled("General Statistics", Style::default().fg(Color::Cyan))
                ]);

                let top_matches_block = Block::default()
                    .title(top_matches_title)
                    .borders(Borders::ALL);

                f.render_widget(top_matches_block.clone(), main_horizontal_split[1]);

                let inner_area = top_matches_block.inner(main_horizontal_split[1]);

                let content_layout = Layout::default()
                    .direction(Direction::Vertical)
                    .constraints(
                        [
                            Constraint::Length(12),
                            Constraint::Min(1),
                        ]
                        .as_ref(),
                    )
                    .split(inner_area);

                let stats_layout = Layout::default()
                    .direction(Direction::Vertical)
                    .constraints(
                        [
                            Constraint::Length(3),
                            Constraint::Length(4),
                            Constraint::Length(5),
                        ]
                        .as_ref(),
                    )
                    .split(content_layout[0]);

                let sequence_length = app.input.len();
                let gc_content = calculate_gc_content(&app.input);
                let at_content = calculate_at_content(&app.input);
                let purine_content = calculate_purine_content(&app.input);
                let pyrimidine_content = calculate_pyrimidine_content(&app.input);

                let total_codons = count_total_codons(&app.input);
                let (complete_codons, incomplete_codons) = count_complete_incomplete_codons(&app.input);
                let start_codon_count = count_start_codons(&app.input);
                let stop_codon_count = count_stop_codons(&app.input);

                let aa_length = calculate_amino_acid_length(&app.input);
                let mol_weight = estimate_molecular_weight(&app.input);
                let hydrophobicity = calculate_hydrophobicity_index(&app.input);
                let (positive_residues, negative_residues) = count_charged_residues(&app.input);
                let orf_count = count_orfs(&app.input);

                let seq_comp_lines = vec![
                    Line::from(vec![
                        Span::styled("Sequence: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{} bp", sequence_length), Style::default().fg(Color::Green)),
                    ]),
                    Line::from(vec![
                        Span::styled("AT: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1}%", at_content), Style::default().fg(Color::Red)),
                        Span::raw("  "),
                        Span::styled("GC: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1}%", gc_content), Style::default().fg(Color::Blue)),
                        Span::raw("  "),
                        Span::styled("Pur: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1}%", purine_content), Style::default().fg(Color::Yellow)),
                        Span::raw("  "),
                        Span::styled("Pyr: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1}%", pyrimidine_content), Style::default().fg(Color::Magenta)),
                    ]),
                ];

                let codon_lines = vec![
                    Line::from(vec![
                        Span::styled("Codons: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{} total", total_codons), Style::default().fg(Color::Green)),
                        Span::raw("  "),
                        Span::styled(format!("{} complete", complete_codons), Style::default().fg(Color::Blue)),
                        Span::raw("  "),
                        Span::styled(format!("{} incomplete", incomplete_codons), Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from(vec![
                        Span::styled("Start (ATG): ", Style::default().fg(Color::White)),
                        Span::styled(format!("{}", start_codon_count), Style::default().fg(Color::Green)),
                        Span::raw("  "),
                        Span::styled("Stop: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{}", stop_codon_count), Style::default().fg(Color::Red)),
                    ]),
                ];

                let protein_lines = vec![
                    Line::from(vec![
                        Span::styled("Amino acids: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{}", aa_length), Style::default().fg(Color::Green)),
                        Span::raw("  "),
                        Span::styled("MW: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1} kDa", mol_weight / 1000.0), Style::default().fg(Color::Blue)),
                    ]),
                    Line::from(vec![
                        Span::styled("Hydrophobicity: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{:.1}%", hydrophobicity), Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from(vec![
                        Span::styled("Charged: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{} (+)", positive_residues), Style::default().fg(Color::Blue)),
                        Span::raw("  "),
                        Span::styled(format!("{} (-)", negative_residues), Style::default().fg(Color::Red)),
                        Span::raw("  "),
                        Span::styled("ORFs: ", Style::default().fg(Color::White)),
                        Span::styled(format!("{}", orf_count), Style::default().fg(Color::Green)),
                    ]),
                ];

                let seq_comp_widget = Paragraph::new(seq_comp_lines)
                    .block(Block::default().title("Sequence Composition").borders(Borders::ALL));

                let codon_widget = Paragraph::new(codon_lines)
                    .block(Block::default().title("Codon Analysis").borders(Borders::ALL));

                let protein_widget = Paragraph::new(protein_lines)
                    .block(Block::default().title("Protein Analysis").borders(Borders::ALL));

                f.render_widget(seq_comp_widget, stats_layout[0]);
                f.render_widget(codon_widget, stats_layout[1]);
                f.render_widget(protein_widget, stats_layout[2]);

                let mut match_lines = Vec::new();

                if !app.histogram_data.is_empty() {
                    let mut display_entries = app.max_histogram_entries;
                    let mut has_low_threshold = false;

                    for (_, similarity) in &app.histogram_data {
                        if *similarity < 0.5 {
                            has_low_threshold = true;
                            break;
                        }
                    }

                    if has_low_threshold {
                        display_entries = display_entries.min(app.histogram_data.len());
                        display_entries = display_entries / 2;
                        display_entries = display_entries.max(1);
                    }

                    for (i, (protein, similarity)) in app.histogram_data.iter().enumerate().take(display_entries) {
                        let percentage = *similarity;

                        let color_index = i % match_colors.len();
                        let color = match_colors[color_index];

                        let match_line = Line::from(vec![
                            Span::styled(
                                format!("{:>2}. ", i + 1),
                                Style::default().fg(Color::White)
                            ),
                            Span::styled(
                                format!("{:<15}", protein.id),
                                Style::default().fg(color)
                            ),
                            Span::raw(" - "),
                            Span::styled(
                                format!("{:.1}%", percentage),
                                Style::default().fg(if percentage >= 50.0 { Color::Green } else { Color::Yellow })
                            ),
                            Span::raw(" match"),
                        ]);

                        match_lines.push(match_line);
                    }
                }

                if match_lines.is_empty() {
                    match_lines.push(Line::from(vec![
                        Span::styled(
                            "No significant matches found",
                            Style::default().fg(Color::DarkGray)
                        ),
                    ]));
                }

                let matches_widget = Paragraph::new(match_lines)
                    .wrap(ratatui::widgets::Wrap { trim: true });

                f.render_widget(matches_widget, content_layout[1]);
            }

            let current_confidence_color = if app.current_strand_confidence > app.opposite_strand_confidence {
                Color::Green
            } else if app.current_strand_confidence < app.opposite_strand_confidence {
                Color::Red
            } else {
                Color::Yellow
            };

            let opposite_confidence_color = if app.opposite_strand_confidence > app.current_strand_confidence {
                Color::Green
            } else if app.opposite_strand_confidence < app.current_strand_confidence {
                Color::Red
            } else {
                Color::Yellow
            };

            let mut confidence_spans = Vec::new();

            if app.closest_protein.is_some() {
                confidence_spans.push(Span::raw("Confidence: "));

                confidence_spans.push(Span::styled(
                    format!("Current {:.1}%", app.current_strand_confidence),
                    Style::default().fg(current_confidence_color)
                ));

                confidence_spans.push(Span::raw(" / "));

                confidence_spans.push(Span::styled(
                    format!("Opposite {:.1}%", app.opposite_strand_confidence),
                    Style::default().fg(opposite_confidence_color)
                ));

                if app.current_strand_confidence > app.opposite_strand_confidence {
                    confidence_spans.push(Span::styled(
                        " (Current strand more likely)",
                        Style::default().fg(Color::Green)
                    ));
                } else if app.current_strand_confidence < app.opposite_strand_confidence {
                    confidence_spans.push(Span::styled(
                        " (Opposite strand more likely)",
                        Style::default().fg(Color::Red)
                    ));
                }
            } else {
                confidence_spans.push(Span::styled(
                    "No protein match found yet",
                    Style::default().fg(Color::DarkGray)
                ));
            }

            let confidence_text = Line::from(confidence_spans);
            let confidence_widget = Paragraph::new(vec![confidence_text])
                .block(Block::default().title("Confidence").borders(Borders::ALL))
                .alignment(ratatui::layout::Alignment::Center);
            f.render_widget(confidence_widget, chunks[5]);
        })?;

        if let Event::Key(key) = event::read()? {
            match key.code {
                KeyCode::Char('q') => break,
                KeyCode::Char(' ') => {},
                KeyCode::Char('s') => app.toggle_strand_mode(),
                KeyCode::Char(c) => app.on_key(c),
                KeyCode::Backspace => app.on_backspace(),
                KeyCode::Esc => break,
                _ => {}
            }
        }
    }

    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    Ok(())
}

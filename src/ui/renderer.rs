use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph, Gauge},
    Frame,
};

use crate::{
    App,
    protein::DatasetProgress,
    sequence::*,
    ui::{format_triplets, create_codon_completion_display},
};

// Helper functions to eliminate code duplication
fn create_conditional_style(condition: bool, true_color: Color, false_color: Color) -> Style {
    if condition {
        Style::default().fg(true_color)
    } else {
        Style::default().fg(false_color)
    }
}

fn create_selection_style(is_selected: bool) -> Style {
    if is_selected {
        Style::default().fg(Color::Black).bg(Color::Yellow)
    } else {
        Style::default().fg(Color::White)
    }
}

fn create_match_style(is_match: bool) -> Style {
    if is_match {
        Style::default().fg(Color::Green).bg(Color::DarkGray)
    } else {
        Style::default().fg(Color::Cyan)
    }
}

fn create_labeled_span(label: &str, value: String, color: Color) -> Vec<Span<'_>> {
    vec![
        Span::raw(label.to_string()),
        Span::styled(value, Style::default().fg(color)),
    ]
}

fn create_strand_mode_spans(is_positive_strand: bool) -> Vec<Span<'static>> {
    vec![
        Span::raw("Strand: "),
        Span::styled("[+] Positive", create_conditional_style(is_positive_strand, Color::Green, Color::DarkGray)),
        Span::raw(" / "),
        Span::styled("[-] Negative", create_conditional_style(!is_positive_strand, Color::Yellow, Color::DarkGray)),
    ]
}

fn create_help_widget(help_lines: Vec<Line>) -> Paragraph {
    Paragraph::new(help_lines)
        .block(Block::default()
            .title("Help")
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::Cyan)))
}

pub fn render_ui(f: &mut Frame, app: &App) {
    // Show loading screen if datasets are being loaded
    if app.is_loading_proteins {
        render_loading_screen(f, app);
        return;
    }

    let main_horizontal_split = Layout::default()
        .direction(Direction::Horizontal)
        .margin(2)
        .constraints([
            Constraint::Percentage(70),
            Constraint::Percentage(30),
        ])
        .split(f.area());

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(5),
            Constraint::Length(5),
            Constraint::Length(5),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(main_horizontal_split[0]);

    render_title(f, app, chunks[0]);
    render_sequence_strands(f, app, &chunks[1..4]);
    render_amino_acid_section(f, app, chunks[4]);
    render_right_panel(f, app, main_horizontal_split[1]);
    render_status_bar(f, app, chunks[5]);

    if app.show_protein_searcher {
        if app.show_protein_detail {
            render_protein_detail(f, app);
        } else {
            render_protein_searcher(f, app);
        }
    }
}

fn render_title(f: &mut Frame, app: &App, area: Rect) {
    let mut spans = vec![
        Span::styled("Ribozap", Style::default().fg(Color::Cyan)),
        Span::raw("   "),
        Span::styled(
            format!("Loaded {} small proteins", app.loaded_proteins_count),
            Style::default().fg(Color::Green)
        ),
        Span::raw("   "),
    ];

    spans.extend(create_strand_mode_spans(app.is_positive_strand));

    let title_widget = Paragraph::new(vec![Line::from(spans)])
        .block(Block::default().borders(Borders::ALL));
    f.render_widget(title_widget, area);
}

fn render_sequence_strands(f: &mut Frame, app: &App, areas: &[Rect]) {
    let formatted_input = format_triplets(&app.input);
    let input_text = Line::from(vec![
        Span::raw("Positive Strand: "),
        Span::styled(&formatted_input, Style::default().fg(Color::Green)),
    ]);
    let input_widget = Paragraph::new(vec![input_text])
        .block(Block::default().borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(input_widget, areas[0]);

    let formatted_complementary = format_triplets(&app.complementary);
    let complementary_text = Line::from(vec![
        Span::raw("Negative Strand: "),
        Span::styled(&formatted_complementary, Style::default().fg(Color::Yellow)),
    ]);
    let complementary_widget = Paragraph::new(vec![complementary_text])
        .block(Block::default().borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(complementary_widget, areas[1]);

    let formatted_mrna = format_triplets(&app.mrna);
    let mrna_text = Line::from(vec![
        Span::raw("mRNA:           "),
        Span::styled(&formatted_mrna, Style::default().fg(Color::Magenta)),
    ]);
    let mrna_widget = Paragraph::new(vec![mrna_text])
        .block(Block::default().borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(mrna_widget, areas[2]);
}

fn render_amino_acid_section(f: &mut Frame, app: &App, area: Rect) {
    let amino_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(30),
            Constraint::Percentage(30),
            Constraint::Percentage(40),
        ])
        .split(area);

    render_amino_acid_sequence(f, app, amino_chunks[0]);
    render_codon_completion(f, app, amino_chunks[1]);
    render_protein_match(f, app, amino_chunks[2]);
}

fn render_amino_acid_sequence(f: &mut Frame, app: &App, area: Rect) {
    let mut amino_spans = vec![Span::raw("Amino Acids: ")];

    for (amino, color) in app.amino_acids_colored.iter() {
        amino_spans.push(Span::styled(amino, Style::default().fg(*color)));
    }

    let amino_widget = Paragraph::new(vec![Line::from(amino_spans)])
        .block(Block::default().title("Amino Acid Sequence").borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(amino_widget, area);
}

fn render_codon_completion(f: &mut Frame, app: &App, area: Rect) {
    let partial_codon = app.get_current_partial_codon();
    let codon_completion = create_codon_completion_display(&partial_codon);
    let codon_completion_widget = Paragraph::new(codon_completion)
        .block(Block::default().title("Codon Completion Guide").borders(Borders::ALL));
    f.render_widget(codon_completion_widget, area);
}

fn render_protein_match(f: &mut Frame, app: &App, area: Rect) {
    let protein_text = if app.is_loading_proteins {
        vec![Line::from(vec![
            Span::styled("Loading protein database...", Style::default().fg(Color::Yellow)),
        ])]
    } else if let Some(error) = &app.loading_error {
        vec![Line::from(vec![
            Span::styled(error, Style::default().fg(Color::Red)),
        ])]
    } else if let Some(protein) = &app.closest_protein {
        build_protein_info_lines(protein, &app.matching_positions)
    } else {
        vec![Line::from(vec![
            Span::styled("No matching protein found", Style::default().fg(Color::DarkGray)),
        ])]
    };

    let protein_widget = Paragraph::new(protein_text)
        .block(Block::default().title("Closest Small Protein Match").borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(protein_widget, area);
}

fn build_protein_info_lines(protein: &crate::SmallProtein, matching_positions: &[bool]) -> Vec<Line<'static>> {
    let mut rna_seq_spans = Vec::new();
    let mut triplet_count = 0;

    for (i, c) in protein.rna_seq.chars().enumerate() {
        let is_match = i < matching_positions.len() && matching_positions[i];
        let style = create_match_style(is_match);

        rna_seq_spans.push(Span::styled(c.to_string(), style));

        triplet_count += 1;
        if triplet_count == 3 && i < protein.rna_seq.len() - 1 {
            rna_seq_spans.push(Span::raw(" "));
            triplet_count = 0;
        }
    }

    vec![
        Line::from(create_labeled_span("Species: ", protein.species.clone(), Color::Green)),
        Line::from(create_labeled_span("ID: ", protein.id.clone(), Color::Yellow)),
        Line::from(create_labeled_span("Length: ", protein.length.to_string(), Color::Blue)),
        Line::from(create_labeled_span("Chromosome: ", protein.chromosome.clone(), Color::Cyan)),
        Line::from(create_labeled_span("Start: ", protein.start.to_string(), Color::Green)),
        Line::from(create_labeled_span("Stop: ", protein.stop.to_string(), Color::Yellow)),
        Line::from(create_labeled_span("Strand: ", protein.strand.clone(), Color::Blue)),
        Line::from(create_labeled_span("Blocks: ", protein.blocks.clone(), Color::Cyan)),
        Line::from(create_labeled_span("Start Codon: ", protein.start_codon.clone(), Color::Green)),
        Line::from(create_labeled_span("PhyloCSF Mean: ", protein.phylo_csf_mean.to_string(), Color::Yellow)),
        Line::from({
            let mut spans = vec![Span::raw("RNA Seq: ")];
            spans.extend(rna_seq_spans);
            spans
        }),
        Line::from(create_labeled_span("AA Seq: ", protein.aa_seq.clone(), Color::Magenta)),
    ]
}

fn render_right_panel(f: &mut Frame, app: &App, area: Rect) {
    let right_panel_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(20),
            Constraint::Min(15),
        ])
        .split(area);

    render_sequence_analysis(f, app, right_panel_chunks[0]);
    render_protein_analysis(f, app, right_panel_chunks[1]);
}

fn render_sequence_analysis(f: &mut Frame, app: &App, area: Rect) {
    let confidence_diff = app.current_strand_confidence - app.opposite_strand_confidence;
    let likely_strand = if confidence_diff > 5.0 {
        ("Current strand more likely", Color::Green)
    } else if confidence_diff < -5.0 {
        ("Opposite strand more likely", Color::Yellow)
    } else if confidence_diff.abs() < 0.1 && app.current_strand_confidence == 0.0 {
        ("No confidence data", Color::DarkGray)
    } else {
        ("Similar confidence", Color::White)
    };

    let composition_lines = vec![
        Line::from(vec![
            Span::styled("Composition Analysis", Style::default().fg(Color::Cyan)),
        ]),
        Line::from(vec![Span::raw("")]),
        Line::from(vec![
            Span::raw("GC Content: "),
            Span::styled(format!("{:.1}%", calculate_gc_content(&app.input)), Style::default().fg(Color::Green)),
        ]),
        Line::from(vec![
            Span::raw("AT Content: "),
            Span::styled(format!("{:.1}%", calculate_at_content(&app.input)), Style::default().fg(Color::Yellow)),
        ]),
        Line::from(vec![
            Span::raw("Purine Content: "),
            Span::styled(format!("{:.1}%", calculate_purine_content(&app.input)), Style::default().fg(Color::Blue)),
        ]),
        Line::from(vec![
            Span::raw("Pyrimidine Content: "),
            Span::styled(format!("{:.1}%", calculate_pyrimidine_content(&app.input)), Style::default().fg(Color::Magenta)),
        ]),
        Line::from(vec![Span::raw("")]),
        Line::from(vec![
            Span::raw("Total Codons: "),
            Span::styled(count_total_codons(&app.input).to_string(), Style::default().fg(Color::White)),
        ]),
        Line::from(vec![
            Span::raw("Complete/Incomplete: "),
            Span::styled(format!("{}/{}", count_complete_incomplete_codons(&app.input).0, count_complete_incomplete_codons(&app.input).1), Style::default().fg(Color::White)),
        ]),
        Line::from(vec![
            Span::raw("Start Codons: "),
            Span::styled(count_start_codons(&app.input).to_string(), Style::default().fg(Color::Green)),
        ]),
        Line::from(vec![
            Span::raw("Stop Codons: "),
            Span::styled(count_stop_codons(&app.input).to_string(), Style::default().fg(Color::Red)),
        ]),
        Line::from(vec![
            Span::raw("ORFs: "),
            Span::styled(count_orfs(&app.input).to_string(), Style::default().fg(Color::Cyan)),
        ]),
        Line::from(vec![Span::raw("")]),
        Line::from(vec![
            Span::styled("Strand Confidence:", Style::default().fg(Color::Cyan)),
        ]),
        Line::from(vec![
            Span::raw("Current: "),
            Span::styled(format!("{:.1}%", app.current_strand_confidence), Style::default().fg(Color::Green)),
            Span::raw(" Opposite: "),
            Span::styled(format!("{:.1}%", app.opposite_strand_confidence), Style::default().fg(Color::Yellow)),
        ]),
        Line::from(vec![
            Span::raw("Assessment: "),
            Span::styled(likely_strand.0, Style::default().fg(likely_strand.1)),
        ]),
    ];

    let composition_widget = Paragraph::new(composition_lines)
        .block(Block::default().title("Sequence Analysis").borders(Borders::ALL));
    f.render_widget(composition_widget, area);
}

fn render_protein_analysis(f: &mut Frame, app: &App, area: Rect) {
    let (positive_charges, negative_charges) = count_charged_residues(&app.input);
    let protein_lines = vec![
        Line::from(vec![
            Span::styled("Protein Properties", Style::default().fg(Color::Cyan)),
        ]),
        Line::from(vec![Span::raw("")]),
        Line::from(vec![
            Span::raw("Amino Acid Length: "),
            Span::styled(calculate_amino_acid_length(&app.input).to_string(), Style::default().fg(Color::White)),
        ]),
        Line::from(vec![
            Span::raw("Est. Molecular Weight: "),
            Span::styled(format!("{:.1} Da", estimate_molecular_weight(&app.input)), Style::default().fg(Color::Yellow)),
        ]),
        Line::from(vec![
            Span::raw("Hydrophobicity Index: "),
            Span::styled(format!("{:.1}%", calculate_hydrophobicity_index(&app.input)), Style::default().fg(Color::Blue)),
        ]),
        Line::from(vec![
            Span::raw("Positive Charges: "),
            Span::styled(positive_charges.to_string(), Style::default().fg(Color::Green)),
        ]),
        Line::from(vec![
            Span::raw("Negative Charges: "),
            Span::styled(negative_charges.to_string(), Style::default().fg(Color::Red)),
        ]),
        Line::from(vec![
            Span::raw("Net Charge: "),
            Span::styled((positive_charges - negative_charges).to_string(), Style::default().fg(Color::White)),
        ]),
    ];

    let protein_widget = Paragraph::new(protein_lines)
        .block(Block::default().title("Protein Analysis").borders(Borders::ALL));
    f.render_widget(protein_widget, area);
}

fn render_status_bar(f: &mut Frame, app: &App, area: Rect) {
    let status_text = if app.input.is_empty() {
        "Enter DNA sequence (A, T, G, C). Press 'q' to quit, 's' to toggle strand mode, 'p' for protein searcher."
    } else {
        "Continue typing or press 'q' to quit, 's' to toggle strand mode, 'p' for protein searcher."
    };

    let status_widget = Paragraph::new(vec![Line::from(vec![
        Span::styled(status_text, Style::default().fg(Color::White)),
    ])])
    .block(Block::default().title("Status").borders(Borders::ALL));
    f.render_widget(status_widget, area);
}

fn render_loading_screen(f: &mut Frame, app: &App) {
    let area = f.area();

    // Center the loading dialog
    let loading_area = Rect::new(
        area.width / 4,
        area.height / 3,
        area.width / 2,
        area.height / 3,
    );

    // Clear the background
    f.render_widget(ratatui::widgets::Clear, loading_area);

    // Main loading container
    let loading_block = Block::default()
        .title("Loading Dataset")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan));
    f.render_widget(loading_block, loading_area);

    let inner_area = Rect::new(
        loading_area.x + 1,
        loading_area.y + 1,
        loading_area.width - 2,
        loading_area.height - 2,
    );

    let loading_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Min(3),
        ])
        .split(inner_area);

    // Status text
    let (status_text, progress_ratio) = match &app.dataset_progress {
        Some(DatasetProgress::CheckingCache) => ("Checking local cache...".to_string(), 0.1),
        Some(DatasetProgress::Downloading { bytes_downloaded, total_bytes }) => {
            if let Some(total) = total_bytes {
                let ratio = (*bytes_downloaded as f64) / (*total as f64);
                (format!("Downloading... {:.1} MB / {:.1} MB",
                    *bytes_downloaded as f64 / 1_048_576.0,
                    *total as f64 / 1_048_576.0), ratio * 0.7 + 0.1) // 10% to 80%
            } else {
                (format!("Downloading... {:.1} MB",
                    *bytes_downloaded as f64 / 1_048_576.0), 0.4)
            }
        },
        Some(DatasetProgress::Extracting) => ("Extracting compressed file...".to_string(), 1.0), // 100% for extracting
        Some(DatasetProgress::Parsing { .. }) => ("Loading complete!".to_string(), 1.0), // Treat parsing as complete
        Some(DatasetProgress::Complete) => ("Loading complete!".to_string(), 1.0),
        Some(DatasetProgress::Error(err)) => (format!("Error: {err}"), 0.0),
        None => ("Initializing...".to_string(), 0.0),
    };

    let status_widget = Paragraph::new(vec![Line::from(vec![
        Span::styled(status_text, Style::default().fg(Color::White)),
    ])])
    .block(Block::default().borders(Borders::ALL).title("Status"));
    f.render_widget(status_widget, loading_chunks[0]);

    // Progress bar using ratatui's built-in Gauge
    let progress_percentage = (progress_ratio * 100.0) as u16;
    let gauge = Gauge::default()
        .block(Block::default().borders(Borders::ALL).title("Progress"))
        .gauge_style(Style::default().fg(Color::Green))
        .percent(progress_percentage)
        .label(format!("{progress_percentage}%"));
    f.render_widget(gauge, loading_chunks[1]);

    // Data location info
    let data_dir = crate::protein::dataset::get_data_dir().unwrap_or_else(|_| std::path::PathBuf::from("Unable to determine data directory"));
    let location_text = format!("Data stored in: {}", data_dir.display());
    let location_widget = Paragraph::new(vec![Line::from(vec![
        Span::styled(location_text, Style::default().fg(Color::DarkGray)),
    ])])
    .block(Block::default().borders(Borders::ALL).title("Storage Location"))
    .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(location_widget, loading_chunks[2]);
}

fn render_protein_searcher(f: &mut Frame, app: &App) {
    let area = f.area();
    let popup_area = Rect::new(
        area.width / 6,
        area.height / 6,
        area.width * 2 / 3,
        area.height * 2 / 3,
    );

    f.render_widget(
        ratatui::widgets::Clear,
        popup_area,
    );

    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::White))
            .title("Protein Searcher"),
        popup_area,
    );

    let inner_area = Rect::new(
        popup_area.x + 1,
        popup_area.y + 1,
        popup_area.width - 2,
        popup_area.height - 2,
    );

    let searcher_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Length(4),
            Constraint::Min(5),
            Constraint::Length(8),
            Constraint::Length(3),
        ])
        .split(inner_area);

    let field_name = app.get_search_field_name();
    let mode_text = if app.multi_search_mode {
        format!("Multi-Search Mode | Current: {field_name}")
    } else {
        format!("Single Search | Field: {field_name}")
    };

    let field_selector = Paragraph::new(vec![Line::from(vec![
        Span::styled(mode_text, create_conditional_style(app.multi_search_mode, Color::Green, Color::Yellow)),
        Span::raw(" (Tab/Shift+Tab to change, Ctrl+T to toggle mode)"),
    ])])
    .block(Block::default()
        .title("Search Mode")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan)));
    f.render_widget(field_selector, searcher_chunks[0]);

    let search_input = Paragraph::new(vec![Line::from(vec![
        Span::styled(&app.searcher_input, Style::default().fg(Color::White)),
        Span::styled("█", Style::default().fg(Color::Yellow)),
    ])])
    .block(Block::default()
        .title("Search Query")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan)));
    f.render_widget(search_input, searcher_chunks[1]);

    let active_filters = app.get_active_filters();
    let filter_lines = if active_filters.is_empty() {
        vec![Line::from(vec![
            Span::styled("No active filters", Style::default().fg(Color::DarkGray)),
        ])]
    } else {
        active_filters.iter().map(|(field, value)| {
            let field_name = match field {
                crate::SearchField::Species => "Species",
                crate::SearchField::Id => "ID",
                crate::SearchField::Chromosome => "Chromosome",
                crate::SearchField::Strand => "Strand",
                crate::SearchField::StartCodon => "Start Codon",
                crate::SearchField::MinLength => "Min Length",
                crate::SearchField::MaxLength => "Max Length",
                crate::SearchField::MinPhyloCSF => "Min PhyloCSF",
                crate::SearchField::MaxPhyloCSF => "Max PhyloCSF",
            };
            Line::from(vec![
                Span::styled(field_name, Style::default().fg(Color::Yellow)),
                Span::raw(": "),
                Span::styled(value, Style::default().fg(Color::White)),
            ])
        }).collect()
    };

    let filters_widget = Paragraph::new(filter_lines)
        .block(Block::default()
            .title("Active Filters")
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::Cyan)))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(filters_widget, searcher_chunks[2]);

    let results_lines: Vec<Line> = app.filtered_proteins.iter().enumerate().map(|(i, protein)| {
        let is_selected = i == app.selected_protein_index;
        Line::from(vec![
            Span::styled(
                format!("{}: {} ({})", protein.id, protein.species, protein.length),
                create_selection_style(is_selected),
            ),
        ])
    }).collect();

    let results_widget = Paragraph::new(results_lines)
        .block(Block::default()
            .title(format!("Results ({}/{})", app.filtered_proteins.len(), app.small_proteins.len()))
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::Cyan)))
        .wrap(ratatui::widgets::Wrap { trim: true });
    f.render_widget(results_widget, searcher_chunks[3]);

    if !app.filtered_proteins.is_empty() && app.selected_protein_index < app.filtered_proteins.len() {
        let selected_protein = &app.filtered_proteins[app.selected_protein_index];
        let details_lines = vec![
            Line::from(vec![
                Span::raw("ID: "),
                Span::styled(&selected_protein.id, Style::default().fg(Color::Yellow)),
            ]),
            Line::from(vec![
                Span::raw("Species: "),
                Span::styled(&selected_protein.species, Style::default().fg(Color::Green)),
            ]),
            Line::from(vec![
                Span::raw("Chr: "),
                Span::styled(&selected_protein.chromosome, Style::default().fg(Color::Cyan)),
                Span::raw(" Start: "),
                Span::styled(selected_protein.start.to_string(), Style::default().fg(Color::Blue)),
                Span::raw(" Stop: "),
                Span::styled(selected_protein.stop.to_string(), Style::default().fg(Color::Blue)),
            ]),
            Line::from(vec![
                Span::raw("Strand: "),
                Span::styled(&selected_protein.strand, Style::default().fg(Color::Magenta)),
                Span::raw(" Length: "),
                Span::styled(selected_protein.length.to_string(), Style::default().fg(Color::White)),
                Span::raw(" PhyloCSF: "),
                Span::styled(format!("{:.2}", selected_protein.phylo_csf_mean), Style::default().fg(Color::Yellow)),
            ]),
            Line::from(vec![
                Span::raw("Start Codon: "),
                Span::styled(&selected_protein.start_codon, Style::default().fg(Color::Green)),
            ]),
            Line::from(vec![
                Span::raw("AA Seq: "),
                Span::styled(
                    if selected_protein.aa_seq.len() > 50 {
                        format!("{}...", &selected_protein.aa_seq[..50])
                    } else {
                        selected_protein.aa_seq.clone()
                    },
                    Style::default().fg(Color::Magenta),
                ),
            ]),
        ];

        let details_widget = Paragraph::new(details_lines)
            .block(Block::default()
                .title("Selected Protein Details")
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::Cyan)))
            .wrap(ratatui::widgets::Wrap { trim: true });
        f.render_widget(details_widget, searcher_chunks[4]);
    } else {
        let no_selection = Paragraph::new(vec![Line::from(vec![
            Span::styled("No protein selected", Style::default().fg(Color::DarkGray)),
        ])])
        .block(Block::default()
            .title("Selected Protein Details")
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::Cyan)));
        f.render_widget(no_selection, searcher_chunks[4]);
    }

    let help_lines = if app.multi_search_mode {
        vec![
            Line::from(vec![
                Span::styled("Ctrl+T: Toggle Mode | Ctrl+A: Add Filter | Ctrl+C: Clear Current | Ctrl+X: Clear All", Style::default().fg(Color::Green)),
            ]),
            Line::from(vec![
                Span::styled("↑/↓: Navigate | Enter: Select | Tab: Change field | Esc: Close", Style::default().fg(Color::White)),
            ]),
        ]
    } else {
        vec![
            Line::from(vec![
                Span::styled("Ctrl+T: Multi-Search Mode | ↑/↓: Navigate | Enter: Select | Tab: Change field | Esc: Close", Style::default().fg(Color::White)),
            ]),
        ]
    };

    let help_widget = create_help_widget(help_lines);
    f.render_widget(help_widget, searcher_chunks[5]);
}

fn render_protein_detail(f: &mut Frame, app: &App) {
    let area = f.area();
    let popup_area = Rect::new(
        area.width / 6,
        area.height / 6,
        area.width * 2 / 3,
        area.height * 2 / 3,
    );

    f.render_widget(
        ratatui::widgets::Clear,
        popup_area,
    );

    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::White))
            .title("Protein Detail"),
        popup_area,
    );

    let inner_area = Rect::new(
        popup_area.x + 1,
        popup_area.y + 1,
        popup_area.width - 2,
        popup_area.height - 2,
    );

    let detail_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(1),
            Constraint::Length(2),
        ])
        .split(inner_area);

    let header_text = if let Some(protein) = &app.detailed_protein {
        format!("{} - {}", protein.id, protein.species)
    } else {
        "Protein Detail".to_string()
    };

    let header = Paragraph::new(vec![Line::from(vec![
        Span::styled(&header_text, Style::default().fg(Color::Cyan)),
    ])])
    .block(Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan)));
    f.render_widget(header, detail_chunks[0]);

    if let Some(protein) = &app.detailed_protein {
        let sequence_lines = vec![
            Line::from(vec![
                Span::raw("ID: "),
                Span::styled(&protein.id, Style::default().fg(Color::Yellow)),
            ]),
            Line::from(vec![
                Span::raw("Species: "),
                Span::styled(&protein.species, Style::default().fg(Color::Green)),
            ]),
            Line::from(vec![
                Span::raw("Chromosome: "),
                Span::styled(&protein.chromosome, Style::default().fg(Color::Cyan)),
            ]),
            Line::from(vec![
                Span::raw("Strand: "),
                Span::styled(&protein.strand, Style::default().fg(Color::Magenta)),
            ]),
            Line::from(vec![
                Span::raw("Start: "),
                Span::styled(protein.start.to_string(), Style::default().fg(Color::Green)),
                Span::raw(" Stop: "),
                Span::styled(protein.stop.to_string(), Style::default().fg(Color::Red)),
            ]),
            Line::from(vec![
                Span::raw("Length: "),
                Span::styled(protein.length.to_string(), Style::default().fg(Color::Blue)),
            ]),
            Line::from(vec![
                Span::raw("Blocks: "),
                Span::styled(protein.blocks.clone(), Style::default().fg(Color::Cyan)),
            ]),
            Line::from(vec![
                Span::raw("Start Codon: "),
                Span::styled(protein.start_codon.clone(), Style::default().fg(Color::Green)),
            ]),
            Line::from(vec![
                Span::raw("PhyloCSF Mean: "),
                Span::styled(protein.phylo_csf_mean.to_string(), Style::default().fg(Color::Yellow)),
            ]),
            Line::from(vec![
                Span::raw("RNA Seq: "),
                Span::styled(protein.rna_seq.clone(), Style::default().fg(Color::White)),
            ]),
            Line::from(vec![
                Span::raw("AA Seq: "),
                Span::styled(protein.aa_seq.clone(), Style::default().fg(Color::Magenta)),
            ]),
        ];

        let sequence_widget = Paragraph::new(sequence_lines)
            .block(Block::default().title("Sequence Details").borders(Borders::ALL))
            .wrap(ratatui::widgets::Wrap { trim: true });
        f.render_widget(sequence_widget, detail_chunks[1]);
    } else {
        let no_detail = Paragraph::new(vec![Line::from(vec![
            Span::styled("No protein detail available", Style::default().fg(Color::DarkGray)),
        ])])
        .block(Block::default().title("Sequence Details").borders(Borders::ALL));
        f.render_widget(no_detail, detail_chunks[1]);
    }

    let help_lines = vec![
        Line::from(vec![
            Span::styled("Enter: Select & Close | Esc: Return to Search | ↑/↓: Scroll", Style::default().fg(Color::White)),
        ]),
    ];

    let help_widget = create_help_widget(help_lines);
    f.render_widget(help_widget, detail_chunks[2]);
}

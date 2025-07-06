use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph},
    Frame,
};

use crate::{
    App,
    sequence::{
        calculate_gc_content, calculate_at_content, calculate_purine_content,
        calculate_pyrimidine_content, count_total_codons, count_complete_incomplete_codons,
        count_start_codons, count_stop_codons, calculate_amino_acid_length,
        estimate_molecular_weight, calculate_hydrophobicity_index,
        count_charged_residues, count_orfs,
    },
    ui::{format_triplets, create_codon_completion_display},
};

pub fn render_ui(f: &mut Frame, app: &App) {
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
        let style = if is_match {
            Style::default().fg(Color::Green).bg(Color::DarkGray)
        } else {
            Style::default().fg(Color::Cyan)
        };

        rna_seq_spans.push(Span::styled(c.to_string(), style));

        triplet_count += 1;
        if triplet_count == 3 && i < protein.rna_seq.len() - 1 {
            rna_seq_spans.push(Span::raw(" "));
            triplet_count = 0;
        }
    }

    vec![
        Line::from(vec![
            Span::raw("Species: "),
            Span::styled(protein.species.clone(), Style::default().fg(Color::Green)),
        ]),
        Line::from(vec![
            Span::raw("ID: "),
            Span::styled(protein.id.clone(), Style::default().fg(Color::Yellow)),
        ]),
        Line::from(vec![
            Span::raw("Length: "),
            Span::styled(protein.length.to_string(), Style::default().fg(Color::Blue)),
        ]),
        Line::from(vec![
            Span::raw("Chromosome: "),
            Span::styled(protein.chromosome.clone(), Style::default().fg(Color::Cyan)),
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
            Span::styled(protein.strand.clone(), Style::default().fg(Color::Blue)),
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
        Line::from({
            let mut spans = vec![Span::raw("RNA Seq: ")];
            spans.extend(rna_seq_spans);
            spans
        }),
        Line::from(vec![
            Span::raw("AA Seq: "),
            Span::styled(protein.aa_seq.clone(), Style::default().fg(Color::Magenta)),
        ]),
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
            Span::styled((positive_charges as i32 - negative_charges as i32).to_string(), Style::default().fg(Color::White)),
        ]),
    ];

    let protein_widget = Paragraph::new(protein_lines)
        .block(Block::default().title("Protein Analysis").borders(Borders::ALL));
    f.render_widget(protein_widget, area);
}

fn render_status_bar(f: &mut Frame, app: &App, area: Rect) {
    let status_text = if app.input.is_empty() {
        "Enter DNA sequence (A, T, G, C). Press 'q' to quit, 's' to toggle strand mode."
    } else {
        "Continue typing or press 'q' to quit, 's' to toggle strand mode."
    };

    let status_widget = Paragraph::new(vec![Line::from(vec![
        Span::styled(status_text, Style::default().fg(Color::White)),
    ])])
    .block(Block::default().title("Status").borders(Borders::ALL));
    f.render_widget(status_widget, area);
}

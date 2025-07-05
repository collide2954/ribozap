//! Main entry point for the Ribozap application

use std::error::Error;
use std::io;
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
    widgets::{Block, Borders, Paragraph},
    Terminal,
};

use ribozap::{
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
                        Constraint::Length(3),   // Title - fixed
                        Constraint::Length(5),   // Positive strand - fixed
                        Constraint::Length(5),   // Negative strand - fixed
                        Constraint::Length(5),   // mRNA - fixed
                        Constraint::Min(12),     // Amino acid sections - use remaining space
                        Constraint::Length(3),   // Status - absolutely fixed at 3 lines
                    ]
                    .as_ref(),
                )
                .split(main_horizontal_split[0]);

            let title = "Ribozap";

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

            // Remove the old standalone strand analysis box - it's no longer needed

            // Right panel with detailed analysis (now only 2 sections)
            let right_panel_chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(20),
                    Constraint::Min(15),
                ].as_ref())
                .split(main_horizontal_split[1]);

            // Sequence composition analysis with strand confidence
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
            f.render_widget(composition_widget, right_panel_chunks[0]);

            // Protein properties
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
            f.render_widget(protein_widget, right_panel_chunks[1]);

            // Status bar at bottom - use the full chunk area directly
            let status_text = if app.input.is_empty() {
                "Enter DNA sequence (A, T, G, C). Press 'q' to quit, 's' to toggle strand mode."
            } else {
                "Continue typing or press 'q' to quit, 's' to toggle strand mode."
            };

            let status = Paragraph::new(vec![Line::from(vec![
                Span::styled(status_text, Style::default().fg(Color::White)),
            ])])
            .block(Block::default().title("Status").borders(Borders::ALL));
            f.render_widget(status, chunks[5]);
        })?;

        if let Event::Key(key) = event::read()? {
            match key.code {
                KeyCode::Char('q') => break,
                KeyCode::Char('s') => {
                    app.toggle_strand_mode();
                }
                KeyCode::Char(c) if c.is_ascii_alphabetic() => {
                    let upper_c = c.to_uppercase().next().unwrap();
                    if upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C' {
                        app.on_key(upper_c);
                    }
                }
                KeyCode::Backspace => {
                    app.on_backspace();
                }
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

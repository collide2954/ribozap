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
    layout::{Alignment, Constraint, Direction, Layout},
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

            let title = "Ribozap";
            let mut spans = vec![
                Span::styled(title, Style::default().fg(Color::Yellow)),
                Span::raw(" - "),
                Span::styled(
                    format!("Proteins Loaded: {}", app.loaded_proteins_count),
                    Style::default().fg(Color::Green)
                ),
            ];

            if app.is_positive_strand {
                spans.extend([
                    Span::raw(" - "),
                    Span::styled("Mode: ", Style::default().fg(Color::White)),
                    Span::styled("[+] Positive", Style::default().fg(Color::Green)),
                    Span::raw(" / "),
                    Span::styled("[-] Negative", Style::default().fg(Color::Yellow)),
                ]);
            } else {
                spans.extend([
                    Span::raw(" - "),
                    Span::styled("Mode: ", Style::default().fg(Color::White)),
                    Span::styled("[+] Positive", Style::default().fg(Color::Yellow)),
                    Span::raw(" / "),
                    Span::styled("[-] Negative", Style::default().fg(Color::Green)),
                ]);
            }

            let strand_mode_spans = if app.is_positive_strand {
                [
                    Span::raw(" - "),
                    Span::styled("[+] Positive", Style::default().fg(Color::Green)),
                    Span::raw(" / "),
                    Span::styled("[-] Negative", Style::default().fg(Color::Yellow)),
                ]
            } else {
                [
                    Span::raw(" - "),
                    Span::styled("[+] Positive", Style::default().fg(Color::Yellow)),
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
                Span::styled(&formatted_complementary, Style::default().fg(Color::Red)),
            ]);
            let complementary = Paragraph::new(vec![complementary_text])
                .block(Block::default().borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(complementary, chunks[2]);

            let formatted_mrna = format_triplets(&app.mrna);
            let mrna_text = Line::from(vec![
                Span::raw("mRNA: "),
                Span::styled(&formatted_mrna, Style::default().fg(Color::Blue)),
            ]);
            let mrna = Paragraph::new(vec![mrna_text])
                .block(Block::default().borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(mrna, chunks[3]);

            let amino_chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints(
                    [
                        Constraint::Length(6),
                        Constraint::Length(10),
                        Constraint::Min(1),
                    ]
                    .as_ref(),
                )
                .split(chunks[4]);

            let amino_acid_spans = app.amino_acids_colored.iter()
                .enumerate()
                .flat_map(|(i, (amino, color))| {
                    let mut spans = Vec::new();
                    if i > 0 {
                        spans.push(Span::raw(" "));
                    }
                    spans.push(Span::styled(amino.clone(), Style::default().fg(*color)));
                    spans
                })
                .collect::<Vec<Span>>();

            let amino_text = Line::from(amino_acid_spans);
            let amino_acids = Paragraph::new(vec![amino_text])
                .block(Block::default().title("Amino Acids").borders(Borders::ALL))
                .wrap(ratatui::widgets::Wrap { trim: true });
            f.render_widget(amino_acids, amino_chunks[0]);

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
                    if i < app.matching_positions.len() && app.matching_positions[i] {
                        rna_seq_spans.push(Span::styled(c.to_string(), Style::default().fg(Color::Green)));
                    } else {
                        rna_seq_spans.push(Span::styled(c.to_string(), Style::default().fg(Color::Red)));
                    }

                    triplet_count += 1;
                    if triplet_count % 3 == 0 && i < rna_seq_clone.len() - 1 {
                        rna_seq_spans.push(Span::raw(" "));
                    }
                }

                vec![
                    Line::from(vec![
                        Span::styled("Closest Match: ", Style::default().fg(Color::Cyan)),
                        Span::styled(&protein.id, Style::default().fg(Color::Yellow)),
                    ]),
                    Line::from(vec![
                        Span::styled("Species: ", Style::default().fg(Color::Cyan)),
                        Span::styled(&protein.species, Style::default().fg(Color::White)),
                    ]),
                    Line::from(vec![
                        Span::styled("RNA Sequence: ", Style::default().fg(Color::Cyan)),
                    ]),
                    Line::from(rna_seq_spans),
                ]
            } else {
                vec![Line::from(vec![
                    Span::styled("No protein matches found", Style::default().fg(Color::Gray)),
                ])]
            };

            let title_spans = vec![
                Span::styled("Closest Protein Match", Style::default().fg(Color::Cyan)),
            ];

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
                } else {
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
                .alignment(Alignment::Center);
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
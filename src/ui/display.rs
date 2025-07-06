use ratatui::{
    style::{Color, Style},
    text::{Line, Span},
};
use bio_seq::prelude::*;
use bio_seq::translation::{TranslationTable, STANDARD};
use crate::ui::colors::get_amino_acid_color;

pub fn format_triplets(sequence: &str) -> String {
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

pub fn create_codon_completion_display(partial_codon: &str) -> Vec<Line<'static>> {
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

    let nucleotides = ['T', 'C', 'A', 'G'];
    let dna_partial_codon = partial_codon.replace('U', "T");

    match dna_partial_codon.len() {
        0 => {
            lines.push(Line::from(vec![
                Span::styled("Start a new codon with any base:", Style::default().fg(Color::White)),
            ]));

            let mut first_options = Vec::new();
            for &base in &nucleotides {
                let display_base = if base == 'T' { 'U' } else { base };
                first_options.push(Span::styled(format!("{display_base} "), Style::default().fg(Color::Cyan)));
            }
            lines.push(Line::from(first_options));
        },
        1 => {
            let first_base = dna_partial_codon.chars().next().unwrap();
            let display_first = if first_base == 'T' { 'U' } else { first_base };
            lines.push(Line::from(vec![
                Span::styled(format!("With first base {display_first}, add second base:"), Style::default().fg(Color::White)),
            ]));

            for &second_base in &nucleotides {
                let mut row = Vec::new();
                let display_second = if second_base == 'T' { 'U' } else { second_base };
                row.push(Span::styled(format!("{display_first}{display_second}_ → "), Style::default().fg(Color::Cyan)));

                let mut possible_aminos = Vec::new();
                for &third_base in &nucleotides {
                    let codon = format!("{first_base}{second_base}{third_base}");
                    if let Ok(codon_seq) = codon.parse::<Seq<Dna>>() {
                        if codon_seq.len() == 3 {
                            let amino = STANDARD.to_amino(&codon_seq);
                            possible_aminos.push(amino.to_string());
                        }
                    }
                }

                possible_aminos.sort();
                possible_aminos.dedup();
                let mut colored_amino_list = Vec::new();
                for (i, amino) in possible_aminos.iter().enumerate() {
                    if i > 0 {
                        colored_amino_list.push(Span::raw("/"));
                    }
                    colored_amino_list.push(Span::styled(amino.clone(), Style::default().fg(get_amino_acid_color(amino))));
                }
                row.extend(colored_amino_list);

                lines.push(Line::from(row));
            }
        },
        2 => {
            let first_base = dna_partial_codon.chars().next().unwrap();
            let second_base = dna_partial_codon.chars().nth(1).unwrap();
            let display_first = if first_base == 'T' { 'U' } else { first_base };
            let display_second = if second_base == 'T' { 'U' } else { second_base };

            lines.push(Line::from(vec![
                Span::styled(format!("With bases {display_first}{display_second}, complete codon with:"), Style::default().fg(Color::White)),
            ]));

            for &third_base in &nucleotides {
                let codon = format!("{first_base}{second_base}{third_base}");
                let display_third = if third_base == 'T' { 'U' } else { third_base };

                let amino = if let Ok(codon_seq) = codon.parse::<Seq<Dna>>() {
                    if codon_seq.len() == 3 {
                        STANDARD.to_amino(&codon_seq).to_string()
                    } else {
                        "?".to_string()
                    }
                } else {
                    "?".to_string()
                };
                let color = get_amino_acid_color(&amino);

                lines.push(Line::from(vec![
                    Span::styled(format!("{display_first}{display_second}{display_third} → "), Style::default().fg(Color::Cyan)),
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
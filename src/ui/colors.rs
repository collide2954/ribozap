//! Color definitions for amino acids and UI elements

use ratatui::style::Color;

/// Get the display color for an amino acid
pub fn get_amino_acid_color(amino: &str) -> Color {
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
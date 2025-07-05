//! Ribozap - DNA/RNA sequence analysis tool
//! 
//! This library provides functionality for analyzing DNA and RNA sequences,
//! including conversion between different nucleic acid forms, codon analysis,
//! and protein matching.

pub mod app;
pub mod protein;
pub mod sequence;
pub mod ui;

// Re-export main types for convenience
pub use app::App;
pub use protein::SmallProtein;
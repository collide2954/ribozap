use std::error::Error;
use std::io;
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    Terminal,
};

use ribozap::{App, ui::render_ui};

fn main() -> Result<(), Box<dyn Error>> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let mut app = App::new();

    // Start threaded loading immediately
    app.start_threaded_loading();

    loop {
        // Check for loading progress updates from background thread
        if app.is_loading_proteins {
            app.check_loading_progress();
        }

        if !app.is_loading_proteins {
            app.perform_protein_matching_if_needed();
        }

        terminal.draw(|f| render_ui(f, &app))?;

        // Use a timeout for event reading to allow progress updates
        if crossterm::event::poll(std::time::Duration::from_millis(100))? {
            if let Event::Key(key) = event::read()? {
                match key.code {
                    KeyCode::Char('q') => break,
                    KeyCode::Char('r') if app.loading_error.is_some() => {
                        // Retry loading on error
                        app.start_threaded_loading();
                    },
                    _ if app.is_loading_proteins => {
                        // Don't process other keys while loading
                        continue;
                    },
                    _ if app.show_protein_searcher => {
                        match key.code {
                            KeyCode::Char('t') if key.modifiers.contains(crossterm::event::KeyModifiers::CONTROL) => {
                                app.toggle_multi_search_mode();
                            },
                            KeyCode::Char('a') if key.modifiers.contains(crossterm::event::KeyModifiers::CONTROL) => {
                                app.add_current_filter();
                            },
                            KeyCode::Char('c') if key.modifiers.contains(crossterm::event::KeyModifiers::CONTROL) => {
                                app.clear_current_filter();
                            },
                            KeyCode::Char('x') if key.modifiers.contains(crossterm::event::KeyModifiers::CONTROL) => {
                                app.clear_all_filters();
                            },
                            KeyCode::Char(c) if c.is_ascii_alphanumeric() || c == ' ' || c == '.' || c == '-' => {
                                app.searcher_on_key(c);
                            },
                            KeyCode::Backspace => {
                                app.searcher_on_backspace();
                            },
                            KeyCode::Tab => {
                                app.searcher_next_field();
                            },
                            KeyCode::BackTab => {
                                app.searcher_prev_field();
                            },
                            KeyCode::Down => {
                                app.searcher_next_protein();
                            },
                            KeyCode::Up => {
                                app.searcher_prev_protein();
                            },
                            KeyCode::Enter => {
                                if app.show_protein_detail {
                                    app.select_detailed_protein();
                                } else {
                                    app.select_current_protein();
                                }
                            },
                            KeyCode::Esc => {
                                if app.show_protein_detail {
                                    app.return_to_search();
                                } else {
                                    app.show_protein_searcher = false;
                                }
                            },
                            _ => {}
                        }
                    },
                    KeyCode::Char('p') => app.toggle_protein_searcher(),
                    KeyCode::Char('s') => {
                        app.toggle_strand_mode();
                    },
                    KeyCode::Char(c) if c.is_ascii_alphabetic() => {
                        let upper_c = c.to_uppercase().next().unwrap();
                        if matches!(upper_c, 'A' | 'T' | 'G' | 'C') {
                            app.on_key(upper_c);
                        }
                    }
                    KeyCode::Backspace => {
                        app.on_backspace();
                    },
                    _ => {}
                }
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

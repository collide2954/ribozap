use std::error::Error;
use std::io;
use std::panic;
use std::path::PathBuf;
use log::{info, warn, error, debug};
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    Terminal,
};

use ribozap::{App, ui::render_ui, logging};

fn setup_logging() -> Result<PathBuf, Box<dyn Error>> {
    // Set log level from environment or default
    logging::set_log_level();

    // Initialize comprehensive logging
    let log_file = logging::init_logging()?;

    // Log system information
    logging::log_system_info();

    Ok(log_file)
}

fn setup_panic_handler() {
    // Set up human-panic for user-friendly panic messages
    human_panic::setup_panic!();

    // Set up custom panic hook for logging
    let default_panic = panic::take_hook();
    panic::set_hook(Box::new(move |panic_info| {
        // Log the panic using our logging system
        let payload = panic_info.payload();
        let panic_msg = if let Some(s) = payload.downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = payload.downcast_ref::<String>() {
            s.clone()
        } else {
            "Unknown panic occurred".to_string()
        };

        let location = if let Some(location) = panic_info.location() {
            format!("{}:{}:{}", location.file(), location.line(), location.column())
        } else {
            "unknown location".to_string()
        };

        // Log critical error with context
        logging::log_critical_error(
            &format!("PANIC at {location}: {panic_msg}"),
            Some("PANIC_HANDLER")
        );

        // Attempt to restore terminal state before panicking
        let _ = disable_raw_mode();
        let _ = execute!(
            io::stdout(),
            LeaveAlternateScreen,
            DisableMouseCapture
        );

        warn!("Terminal state restored after panic");

        // Call the default panic handler
        default_panic(panic_info);
    }));

    info!("Panic handler configured successfully");
}

fn cleanup_terminal() -> Result<(), Box<dyn Error>> {
    debug!("Cleaning up terminal state");
    disable_raw_mode()?;
    execute!(
        io::stdout(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    info!("Terminal cleanup completed successfully");
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize comprehensive logging first
    let log_file = match setup_logging() {
        Ok(file) => {
            info!("Logging initialized successfully");
            Some(file)
        },
        Err(e) => {
            eprintln!("Failed to initialize logging: {e}");
            logging::log_critical_error(&format!("Logging initialization failed: {e}"), Some("STARTUP"));
            None
        }
    };

    // Set up panic handler
    setup_panic_handler();

    info!("RiboZap application starting");
    if let Some(ref file) = log_file {
        info!("Logs are being written to: {file:?}");
    }
    debug!("Setting up terminal");

    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    info!("Terminal initialized successfully");

    let mut app = App::new();
    debug!("App instance created");

    // Start threaded loading immediately
    info!("Starting protein data loading");
    app.start_threaded_loading();

    // Main application loop
    info!("Entering main application loop");
    let loop_result = run_main_loop(&mut terminal, &mut app);

    // Cleanup
    match cleanup_terminal() {
        Ok(_) => {
            info!("Application terminated gracefully");
            logging::log_shutdown();
        },
        Err(e) => {
            error!("Error during terminal cleanup: {e}");
            logging::log_critical_error(&format!("Terminal cleanup failed: {e}"), Some("CLEANUP"));
        }
    }

    // Handle any errors from the main loop
    if let Err(ref e) = loop_result {
        error!("Main loop error: {e}");
        logging::log_critical_error(&format!("Main loop failed: {e}"), Some("MAIN_LOOP"));
    }

    loop_result
}

fn run_main_loop(terminal: &mut Terminal<CrosstermBackend<io::Stdout>>, app: &mut App) -> Result<(), Box<dyn Error>> {
    loop {
        // Check for loading progress updates from background thread
        if app.is_loading_proteins {
            app.check_loading_progress();
        }

        if !app.is_loading_proteins {
            app.perform_protein_matching_if_needed();
        }

        terminal.draw(|f| render_ui(f, app))?;

        // Use a timeout for event reading to allow progress updates
        if event::poll(std::time::Duration::from_millis(100))? {
            if let Event::Key(key) = event::read()? {
                debug!("Key event received: {key:?}");
                match key.code {
                    KeyCode::Char('q') => {
                        info!("Quit command received");
                        break;
                    },
                    KeyCode::Char('r') if app.loading_error.is_some() => {
                        // Retry loading on error
                        warn!("Retrying protein data loading after error");
                        app.start_threaded_loading();
                    },
                    _ if app.is_loading_proteins => {
                        // Don't process other keys while loading
                        continue;
                    },
                    _ if app.show_protein_searcher => {
                        handle_protein_searcher_keys(&key, app)?;
                    },
                    KeyCode::Char('p') => {
                        debug!("Toggling protein searcher");
                        app.toggle_protein_searcher();
                    },
                    KeyCode::Char('s') => {
                        debug!("Toggling strand mode");
                        app.toggle_strand_mode();
                    },
                    KeyCode::Char(c) if c.is_ascii_alphabetic() => {
                        let upper_c = c.to_uppercase().next().unwrap();
                        if matches!(upper_c, 'A' | 'T' | 'G' | 'C') {
                            debug!("Processing nucleotide input: {upper_c}");
                            app.on_key(upper_c);
                        }
                    }
                    KeyCode::Backspace => {
                        debug!("Processing backspace");
                        app.on_backspace();
                    },
                    _ => {
                        debug!("Unhandled key event: {key:?}");
                    }
                }
            }
        }
    }

    Ok(())
}

fn handle_protein_searcher_keys(key: &event::KeyEvent, app: &mut App) -> Result<(), Box<dyn Error>> {
    match key.code {
        KeyCode::Char('t') if key.modifiers.contains(event::KeyModifiers::CONTROL) => {
            debug!("Toggling multi-search mode");
            app.toggle_multi_search_mode();
        },
        KeyCode::Char('a') if key.modifiers.contains(event::KeyModifiers::CONTROL) => {
            debug!("Adding current filter");
            app.add_current_filter();
        },
        KeyCode::Char('c') if key.modifiers.contains(event::KeyModifiers::CONTROL) => {
            debug!("Clearing current filter");
            app.clear_current_filter();
        },
        KeyCode::Char('x') if key.modifiers.contains(event::KeyModifiers::CONTROL) => {
            debug!("Clearing all filters");
            app.clear_all_filters();
        },
        KeyCode::Char(c) if c.is_ascii_alphanumeric() || c == ' ' || c == '.' || c == '-' => {
            app.searcher_on_key(c);
        },
        KeyCode::Backspace => {
            app.searcher_on_backspace();
        },
        KeyCode::Tab => {
            debug!("Moving to next search field");
            app.searcher_next_field();
        },
        KeyCode::BackTab => {
            debug!("Moving to previous search field");
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
                debug!("Selecting detailed protein");
                app.select_detailed_protein();
            } else {
                debug!("Selecting current protein");
                app.select_current_protein();
            }
        },
        KeyCode::Esc => {
            if app.show_protein_detail {
                debug!("Returning to search from detail view");
                app.return_to_search();
            } else {
                debug!("Closing protein searcher");
                app.show_protein_searcher = false;
            }
        },
        _ => {}
    }
    Ok(())
}

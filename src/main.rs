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
    
    loop {
        app.perform_protein_matching_if_needed();
        terminal.draw(|f| render_ui(f, &app))?;
        if let Event::Key(key) = event::read()? {
            match key.code {
                KeyCode::Char('q') => break,
                KeyCode::Char('s') => app.toggle_strand_mode(),
                KeyCode::Char(c) if c.is_ascii_alphabetic() => {
                    let upper_c = c.to_uppercase().next().unwrap();
                    if matches!(upper_c, 'A' | 'T' | 'G' | 'C') {
                        app.on_key(upper_c);
                    }
                }
                KeyCode::Backspace => app.on_backspace(),
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

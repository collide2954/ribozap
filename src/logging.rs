use std::error::Error;
use std::path::PathBuf;
use chrono::Utc;
use log::{info, error, debug};
use env_logger::{Builder, Target};

/// Initialize logging with comprehensive configuration
pub fn init_logging() -> Result<PathBuf, Box<dyn Error>> {
    // Create logs directory if it doesn't exist
    let log_dir = dirs::home_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(".ribozap")
        .join("logs");

    std::fs::create_dir_all(&log_dir)?;

    // Create timestamped log file
    let log_file = log_dir.join(format!("ribozap_{}.log", Utc::now().format("%Y%m%d_%H%M%S")));

    // Set up environment logger with custom format
    Builder::from_default_env()
        .target(Target::Pipe(Box::new(std::fs::File::create(&log_file)?)))
        .format(|buf, record| {
            use std::io::Write;
            writeln!(buf,
                "{} [{}] [{}:{}] [{}] {}",
                Utc::now().format("%Y-%m-%d %H:%M:%S%.3f UTC"),
                record.level(),
                record.module_path().unwrap_or("unknown"),
                record.line().unwrap_or(0),
                std::thread::current().name().unwrap_or("main"),
                record.args()
            )
        })
        .init();

    info!("Logging system initialized");
    info!("Log file: {log_file:?}");
    debug!("Log directory: {log_dir:?}");

    Ok(log_file)
}

/// Set logging level based on environment variable or default
pub fn set_log_level() {
    let level = std::env::var("RIBOZAP_LOG_LEVEL")
        .unwrap_or_else(|_| "info".to_string())
        .to_lowercase();

    let env_filter = match level.as_str() {
        "trace" => "trace",
        "debug" => "debug",
        "info" => "info",
        "warn" => "warn",
        "error" => "error",
        _ => {
            eprintln!("Invalid log level '{level}', defaulting to 'info'");
            "info"
        }
    };

    std::env::set_var("RUST_LOG", format!("ribozap={env_filter}"));
}

/// Log system information at startup
pub fn log_system_info() {
    info!("=== RiboZap Application Starting ===");
    info!("Version: {}", env!("CARGO_PKG_VERSION"));
    info!("Build target: {}", std::env::consts::ARCH);
    info!("Operating system: {}", std::env::consts::OS);

    if let Ok(hostname) = std::env::var("HOSTNAME") {
        info!("Hostname: {hostname}");
    }

    if let Ok(user) = std::env::var("USER") {
        info!("User: {user}");
    }

    info!("Current working directory: {:?}", std::env::current_dir().unwrap_or_default());
    info!("Available CPU cores: {}", num_cpus::get());

    debug!("Environment variables:");
    for (key, value) in std::env::vars() {
        if key.starts_with("RIBOZAP_") || key == "RUST_LOG" {
            debug!("  {key}: {value}");
        }
    }
}

/// Log application shutdown
pub fn log_shutdown() {
    info!("=== RiboZap Application Shutting Down ===");
    info!("Application terminated at {}", Utc::now().format("%Y-%m-%d %H:%M:%S UTC"));
}

/// Create a panic-safe logging function for critical errors
pub fn log_critical_error(error: &str, context: Option<&str>) {
    // Try to log using the standard logger first
    if let Some(ctx) = context {
        error!("CRITICAL ERROR [{ctx}]: {error}");
    } else {
        error!("CRITICAL ERROR: {error}");
    }

    // Also write to stderr as a fallback
    if let Some(ctx) = context {
        eprintln!("[{}] CRITICAL ERROR [{}]: {}",
                 Utc::now().format("%Y-%m-%d %H:%M:%S UTC"), ctx, error);
    } else {
        eprintln!("[{}] CRITICAL ERROR: {}",
                 Utc::now().format("%Y-%m-%d %H:%M:%S UTC"), error);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_logging_initialization() {
        let temp_dir = tempdir().unwrap();
        std::env::set_var("HOME", temp_dir.path());

        let result = init_logging();
        assert!(result.is_ok());

        let log_file = result.unwrap();
        assert!(log_file.exists());
    }
}

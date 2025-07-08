use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use reqwest::blocking::Client;
use log::{info, warn, error, debug, trace};

#[derive(Debug, Clone)]
pub struct SmallProtein {
    pub species: String,
    pub id: String,
    pub rna_seq: String,
    pub aa_seq: String,
    pub length: usize,
    pub chromosome: String,
    pub start: usize,
    pub stop: usize,
    pub strand: String,
    pub blocks: String,
    pub start_codon: String,
    pub phylo_csf_mean: f64,
}

#[derive(Debug, Clone)]
pub enum DatasetProgress {
    CheckingCache,
    Downloading { bytes_downloaded: u64, total_bytes: Option<u64> },
    Extracting,
    Parsing { lines_parsed: usize },
    Complete,
    Error(String),
}

pub fn get_data_dir() -> Result<PathBuf, Box<dyn Error>> {
    let data_dir = dirs::data_dir()
        .ok_or("Could not determine data directory")?
        .join("ribozap");

    debug!("Data directory path: {data_dir:?}");
    fs::create_dir_all(&data_dir)?;
    trace!("Data directory created/verified: {data_dir:?}");

    Ok(data_dir)
}

pub fn download_and_parse_small_protein_dataset() -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    info!("Starting protein dataset download and parsing (without progress callback)");
    download_and_parse_small_protein_dataset_with_progress(None)
}

pub fn download_and_parse_small_protein_dataset_with_progress(
    progress_callback: Option<Box<dyn Fn(DatasetProgress)>>
) -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    let url = "http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt2_LiteratureMining.txt.gz";

    info!("Starting protein dataset download and parsing with progress tracking");
    debug!("Dataset URL: {url}");

    let data_dir = get_data_dir()?;
    let temp_file = data_dir.join("small_protein_dataset.txt.gz");
    let extracted_file = data_dir.join("small_protein_dataset.txt");

    debug!("Temp file path: {temp_file:?}");
    debug!("Extracted file path: {extracted_file:?}");

    if let Some(ref callback) = progress_callback {
        callback(DatasetProgress::CheckingCache);
    }

    if !extracted_file.exists() {
        info!("Extracted file does not exist, checking for compressed file");
        
        if !temp_file.exists() {
            info!("Compressed file does not exist, starting download");
            
            if let Some(ref callback) = progress_callback {
                callback(DatasetProgress::Downloading { bytes_downloaded: 0, total_bytes: None });
            }

            let client = Client::new();
            debug!("HTTP client created, sending request to: {url}");
            
            let mut response = client.get(url).send()
                .map_err(|e| {
                    error!("Failed to send HTTP request: {e}");
                    e
                })?;

            let total_size = response.content_length();
            debug!("Response received, content length: {total_size:?}");
            
            let mut downloaded = 0u64;

            let mut file = File::create(&temp_file)
                .map_err(|e| {
                    error!("Failed to create temp file {temp_file:?}: {e}");
                    e
                })?;
            
            info!("Starting file download to {temp_file:?}");
            let mut buffer = [0; 8192];

            loop {
                let bytes_read = response.read(&mut buffer)
                    .map_err(|e| {
                        error!("Error reading from HTTP response: {e}");
                        e
                    })?;
                
                if bytes_read == 0 {
                    break;
                }

                file.write_all(&buffer[..bytes_read])
                    .map_err(|e| {
                        error!("Error writing to temp file: {e}");
                        e
                    })?;
                
                downloaded += bytes_read as u64;

                if downloaded % (1024 * 1024) == 0 { // Log every MB
                    trace!("Downloaded {downloaded} bytes");
                }

                if let Some(ref callback) = progress_callback {
                    callback(DatasetProgress::Downloading {
                        bytes_downloaded: downloaded,
                        total_bytes: total_size,
                    });
                }
            }
            
            info!("Download completed successfully. Total bytes: {downloaded}");
        } else {
            info!("Compressed file already exists, skipping download");
        }

        info!("Starting file extraction");
        if let Some(ref callback) = progress_callback {
            callback(DatasetProgress::Extracting);
        }

        let compressed_file = File::open(&temp_file)
            .map_err(|e| {
                error!("Failed to open compressed file {temp_file:?}: {e}");
                e
            })?;
        
        let decoder = GzDecoder::new(compressed_file);
        let mut reader = BufReader::new(decoder);
        let mut extracted_content = String::new();
        
        reader.read_to_string(&mut extracted_content)
            .map_err(|e| {
                error!("Failed to decompress file: {e}");
                e
            })?;

        debug!("Decompressed content size: {} bytes", extracted_content.len());

        std::fs::write(&extracted_file, extracted_content)
            .map_err(|e| {
                error!("Failed to write extracted file {extracted_file:?}: {e}");
                e
            })?;
        
        info!("File extraction completed successfully");
    } else {
        info!("Extracted file already exists, proceeding to parsing");
    }

    info!("Starting protein data parsing");
    let file = File::open(&extracted_file)
        .map_err(|e| {
            error!("Failed to open extracted file {extracted_file:?}: {e}");
            e
        })?;
    
    let reader = BufReader::new(file);
    let mut proteins = Vec::new();
    let mut lines_parsed = 0;
    let mut errors_encountered = 0;

    for (line_num, line) in reader.lines().enumerate().skip(1) {
        let line = line
            .map_err(|e| {
                error!("Error reading line {}: {}", line_num + 1, e);
                e
            })?;
        
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            warn!("Line {} has insufficient fields ({}), skipping", line_num + 1, fields.len());
            errors_encountered += 1;
            continue;
        }

        let protein = SmallProtein {
            species: fields[0].to_string(),
            id: fields[1].to_string(),
            rna_seq: fields[2].to_string(),
            aa_seq: fields[3].to_string(),
            length: parse_usize_field(fields[4], line_num + 1, "length", &mut errors_encountered),
            chromosome: fields[5].to_string(),
            start: parse_usize_field(fields[6], line_num + 1, "start", &mut errors_encountered),
            stop: parse_usize_field(fields[7], line_num + 1, "stop", &mut errors_encountered),
            strand: fields[8].to_string(),
            blocks: fields[9].to_string(),
            start_codon: fields[10].to_string(),
            phylo_csf_mean: parse_float_field(fields[11], line_num + 1, "phylo_csf_mean", &mut errors_encountered),
        };

        proteins.push(protein);
        lines_parsed += 1;

        if lines_parsed % 1000 == 0 {
            trace!("Parsed {lines_parsed} lines");
            if let Some(ref callback) = progress_callback {
                callback(DatasetProgress::Parsing { lines_parsed });
            }
        }
    }

    if errors_encountered > 0 {
        warn!("Parsing completed with {errors_encountered} errors encountered");
    }

    info!("Protein data parsing completed successfully. {} proteins loaded", proteins.len());
    
    if let Some(ref callback) = progress_callback {
        callback(DatasetProgress::Complete);
    }

    Ok(proteins)
}

fn parse_float_field(field: &str, line_num: usize, field_name: &str, errors_encountered: &mut usize) -> f64 {
    let trimmed = field.trim();

    // Handle common invalid float values
    match trimmed {
        "" | "NA" | "NULL" | "null" | "N/A" | "n/a" | "-" | "." | "NaN" | "nan" => {
            debug!("Line {}: {} is '{}', using default value 0.0", line_num, field_name, trimmed);
            0.0
        },
        _ => {
            trimmed.parse().unwrap_or_else(|e| {
                // Only warn for truly unexpected parsing errors, not common placeholders
                if !trimmed.is_empty() && !matches!(trimmed, "NA" | "NULL" | "null" | "N/A" | "n/a" | "-" | "." | "NaN" | "nan") {
                    warn!("Line {}: Failed to parse {} '{}': {}", line_num, field_name, trimmed, e);
                    *errors_encountered += 1;
                } else {
                    debug!("Line {}: {} is '{}', using default value 0.0", line_num, field_name, trimmed);
                }
                0.0
            })
        }
    }
}

fn parse_usize_field(field: &str, line_num: usize, field_name: &str, errors_encountered: &mut usize) -> usize {
    let trimmed = field.trim();

    // Handle common invalid usize values
    match trimmed {
        "" | "NA" | "NULL" | "null" | "N/A" | "n/a" | "-" | "." => {
            debug!("Line {}: {} is '{}', using default value 0", line_num, field_name, trimmed);
            0
        },
        _ => {
            trimmed.parse().unwrap_or_else(|e| {
                // Only warn for truly unexpected parsing errors, not common placeholders
                if !trimmed.is_empty() && !matches!(trimmed, "NA" | "NULL" | "null" | "N/A" | "n/a" | "-" | "." ) {
                    warn!("Line {}: Failed to parse {} '{}': {}", line_num, field_name, trimmed, e);
                    *errors_encountered += 1;
                } else {
                    debug!("Line {}: {} is '{}', using default value 0", line_num, field_name, trimmed);
                }
                0
            })
        }
    }
}

use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use reqwest::blocking::Client;

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

    fs::create_dir_all(&data_dir)?;

    Ok(data_dir)
}

pub fn download_and_parse_small_protein_dataset() -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    download_and_parse_small_protein_dataset_with_progress(None)
}

pub fn download_and_parse_small_protein_dataset_with_progress(
    progress_callback: Option<Box<dyn Fn(DatasetProgress)>>
) -> Result<Vec<SmallProtein>, Box<dyn Error>> {
    let url = "http://bigdata.ibp.ac.cn/SmProt/datadownload/SmProt2_LiteratureMining.txt.gz";

    let data_dir = get_data_dir()?;
    let temp_file = data_dir.join("small_protein_dataset.txt.gz");
    let extracted_file = data_dir.join("small_protein_dataset.txt");

    if let Some(ref callback) = progress_callback {
        callback(DatasetProgress::CheckingCache);
    }

    if !extracted_file.exists() {
        if !temp_file.exists() {
            if let Some(ref callback) = progress_callback {
                callback(DatasetProgress::Downloading { bytes_downloaded: 0, total_bytes: None });
            }

            let client = Client::new();
            let mut response = client.get(url).send()?;

            let total_size = response.content_length();
            let mut downloaded = 0u64;

            let mut file = File::create(&temp_file)?;
            let mut buffer = [0; 8192];

            loop {
                let bytes_read = response.read(&mut buffer)?;
                if bytes_read == 0 {
                    break;
                }

                file.write_all(&buffer[..bytes_read])?;
                downloaded += bytes_read as u64;

                if let Some(ref callback) = progress_callback {
                    callback(DatasetProgress::Downloading {
                        bytes_downloaded: downloaded,
                        total_bytes: total_size,
                    });
                }
            }
        }

        if let Some(ref callback) = progress_callback {
            callback(DatasetProgress::Extracting);
        }

        let compressed_file = File::open(&temp_file)?;
        let decoder = GzDecoder::new(compressed_file);
        let mut reader = BufReader::new(decoder);
        let mut extracted_content = String::new();
        reader.read_to_string(&mut extracted_content)?;

        std::fs::write(&extracted_file, extracted_content)?;
    }

    let file = File::open(&extracted_file)?;
    let reader = BufReader::new(file);
    let mut proteins = Vec::new();
    let mut lines_parsed = 0;

    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            continue;
        }

        let protein = SmallProtein {
            species: fields[0].to_string(),
            id: fields[1].to_string(),
            rna_seq: fields[2].to_string(),
            aa_seq: fields[3].to_string(),
            length: fields[4].parse().unwrap_or(0),
            chromosome: fields[5].to_string(),
            start: fields[6].parse().unwrap_or(0),
            stop: fields[7].parse().unwrap_or(0),
            strand: fields[8].to_string(),
            blocks: fields[9].to_string(),
            start_codon: fields[10].to_string(),
            phylo_csf_mean: fields[11].parse().unwrap_or(0.0),
        };

        proteins.push(protein);
        lines_parsed += 1;

        if lines_parsed % 1000 == 0 {
            if let Some(ref callback) = progress_callback {
                callback(DatasetProgress::Parsing { lines_parsed });
            }
        }
    }

    if let Some(ref callback) = progress_callback {
        callback(DatasetProgress::Complete);
    }

    Ok(proteins)
}
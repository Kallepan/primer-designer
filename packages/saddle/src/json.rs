use serde::{Deserialize, Serialize};
// JSON object structs
#[derive(Deserialize, Serialize, Clone)]
pub struct Primer {
    pub id: String,
    pub sequence: String,
    pub tm: f64,
    pub gc_percent: f64,
    pub hairpin_th: f64,
    pub length: usize,
    pub badness: f64,
}
#[derive(Deserialize, Serialize, Clone)]
pub struct Amplicon {
    pub amplicon_name: String,
    pub forward_primers: Vec<Primer>,
    pub reverse_primers: Vec<Primer>
}
#[derive(Deserialize, Serialize)]
pub struct Region {
    pub amplicons: Vec<Amplicon>,
    pub region_name: String,
}
#[derive(Deserialize, Serialize)]
pub struct Pool {
    pub regions: Vec<Region>,
    pub pool_id: String,
}

// SADDLE structs
#[derive(Clone, Serialize)]
pub struct PrimerPair {
    pub region_name: String,
    pub amplicon_name: String,
    pub forward_primer: Primer,
    pub reverse_primer: Primer,
}
#[derive(Clone, Serialize)]
pub struct Set {
    pub primer_pairs: Vec<PrimerPair>,
    pub pool_id: String,
    pub loss: f64,
}

use std::error::Error;

pub fn load_json_from_file(input_file_path: &str) -> Result<Pool, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file_path)?;
    let pools: Pool = serde_json::from_str(&file_content)?;

    Ok(pools)
}

pub fn write_set_to_file(file_path: &str, set: Set) -> Result<(), Box<dyn Error>>{
    let json = serde_json::to_string_pretty(&set).unwrap();
    std::fs::write(file_path, json)?;

    Ok(())
}
pub fn write_losses_to_file(file_path: &str, losses: Vec<f64>) -> Result<(), Box<dyn Error>>{
    let json = serde_json::to_string_pretty(&losses).unwrap();
    std::fs::write(file_path, json)?;

    Ok(())
}
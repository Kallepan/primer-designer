use serde::{Deserialize, Serialize};
// JSON object structs
#[derive(Deserialize, Serialize, Clone)]
pub struct Primer {
    pub primer_sequence: String,
    pub tm: f32,
    pub gc_percent: f32,
    pub hairpin_th: f32,
    pub primer_length: usize,
    pub badness: f32,
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
    pub pool_id: usize,
}

// SADDLE structs
#[derive(Clone, Serialize)]
pub struct AmpliconPrimerPair {
    pub region_name: String,
    pub amplicon_name: String,
    pub forward_primer: Primer,
    pub reverse_primer: Primer,
}
#[derive(Clone, Serialize)]
pub struct Set {
    pub amplicon_primer_pairs: Vec<AmpliconPrimerPair>,
    pub loss: f64,
}
//
use std::error::Error;

pub fn load_json_from_file(input_file_path: &str) -> Result<Vec<Pool>, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file_path)?;
    let pools: Vec<Pool> = serde_json::from_str(&file_content)?;

    Ok(pools)
}

pub fn write_sets_to_file(file_path: &str, set: Vec<Set>) -> Result<(), Box<dyn Error>>{
    let json = serde_json::to_string_pretty(&set).unwrap();
    std::fs::write(file_path, json)?;

    Ok(())
}
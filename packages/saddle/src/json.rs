use serde::{Deserialize, Serialize};
use std::error::Error;

#[derive(Deserialize, Serialize)]
pub struct Primer {
    pub sequence: String,
    pub tm: f32,
    pub gc_percent: f32,
    pub hairpin_th: f32,
    pub length: usize,
}

#[derive(Deserialize, Serialize)]
pub struct Amplicon {
    pub name: String,
    pub forward_primers: Vec<Primer>,
    pub reverse_primers: Vec<Primer>,
}

#[derive(Deserialize, Serialize)]
pub struct Region {
    pub amplicons: Vec<Amplicon>,
    pub name: String, 
}

#[derive(Deserialize, Serialize)]
pub struct Regions {
    pub sequence_id: String,
    pub regions: Vec<Region>,
}

pub fn load_json_from_file(input_file: &str) -> Result<Regions, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file)?;
    let proto_primers: Regions = serde_json::from_str(&file_content)?;

    Ok(proto_primers)
}
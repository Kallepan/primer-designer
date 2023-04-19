use serde::{Deserialize, Serialize};
use serde_json::Result;
use std;

#[derive(Deserialize, Serialize)]
pub struct Primer {
    sequence: String,
    id: String,
    tm: f32,
    gc_percent: f32,
    hairpin_th: f32,
}

#[derive(Deserialize, Serialize)]
pub struct Amplicon {
    id: String,
    forward_primers: Vec<Primer>,
    reverse_primers: Vec<Primer>,
}

#[derive(Deserialize, Serialize)]
pub struct Json {
    amplicons: Vec<Amplicon>,
}

pub fn load_json(input_file: &str) -> Result<Json>{
    let file_content = std::fs::read_to_string(&input_file).expect("Unable to read file");

    let proto_primers: Json = serde_json::from_str(&file_content).expect("Unable to parse JSON");

    Ok(proto_primers)
}
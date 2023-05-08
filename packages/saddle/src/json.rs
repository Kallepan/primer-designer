use crate::types::{Pool, Set};

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
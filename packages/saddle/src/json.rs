use crate::types::{Pool};

use std::error::Error;

pub fn load_json_from_file(input_file_path: &str) -> Result<Vec<Pool>, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file_path)?;
    let pools: Vec<Pool> = serde_json::from_str(&file_content)?;

    Ok(pools)
}
use std::error::Error;
use std::ops::{Bound, RangeBounds};
use crate::json::{Pool, Data, Loss};

pub struct SubsequenceInfo {
    pub seq: String,
    pub start_index: usize,
    pub end_index: usize,
}

pub trait RoundUtils {
    fn round_dp(&self, decimals: usize) -> f64;
}

impl RoundUtils for f64 {
    fn round_dp(&self, decimals: usize) -> f64 {
        let multiplier = 10.0_f64.powi(decimals as i32);
        (self * multiplier).round() / multiplier
    }
}

// Primer (e.g.: string) utility functions
pub trait PrimerUtils {
    /*
        Returns a substring of the string, between the first argument (inclusive) and the second argument (exclusive).
    */
    fn substring(&self, start: usize, len: usize) -> &str;
    fn slice(&self, range: impl RangeBounds<usize>) -> &str;
    fn get_all_substrings_between(&self, min_size: usize, max_size: usize) -> Vec<SubsequenceInfo>;
    fn num_gc(&self) -> usize;
}

impl PrimerUtils for str {
    fn num_gc(&self) -> usize {
        let mut num_gc = 0;
        for nucleotide in self.chars() {
            match nucleotide {
                'G' => num_gc += 1,
                'C' => num_gc += 1,
                _ => continue,
            }
        }
        num_gc
    }

    fn get_all_substrings_between(&self, min_size: usize, max_size: usize) -> Vec<SubsequenceInfo> {
        let mut substrings = Vec::new();
        
        let mut first_pointer = 0;
        let end = self.len();
        let mut second_pointer = first_pointer + min_size;
    
        loop {
            if first_pointer > end { break; }
            if second_pointer > end { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }
            if second_pointer - first_pointer > max_size { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }
            
            let subsequence_info = SubsequenceInfo {
                seq: self.slice(first_pointer..second_pointer).to_string(),
                start_index: first_pointer,
                end_index: second_pointer,
            };
            
            substrings.push(subsequence_info);
            second_pointer += 1;
        }
    
        substrings
    }

    fn substring(&self, start: usize, len: usize) -> &str {
        let mut char_pos = 0;
        let mut byte_start = 0;
        let mut it = self.chars();
        loop {
            if char_pos == start { break; }
            if let Some(c) = it.next() {
                char_pos += 1;
                byte_start += c.len_utf8();
            } else { break; }
        }
        char_pos = 0;
        let mut byte_end = byte_start;
        loop {
            if char_pos == len { break; }
            if let Some(c) = it.next() {
                char_pos += 1;
                byte_end += c.len_utf8();
            } else { break; }
        }
        &self[byte_start..byte_end]
    }

    fn slice(&self, range: impl RangeBounds<usize>) -> &str {
        let start = match range.start_bound() {
            Bound::Included(bound) | Bound::Excluded(bound) => *bound,
            Bound::Unbounded => 0,
        };
        let len = match range.end_bound() {
            Bound::Included(bound) => *bound + 1,
            Bound::Excluded(bound) => *bound,
            Bound::Unbounded => self.len(),
        } - start;
        self.substring(start, len)
    }
}

// JSON utility functions
pub fn load_json_from_file(input_file_path: &str) -> Result<Pool, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file_path)?;
    let pools: Pool = serde_json::from_str(&file_content)?;

    Ok(pools)
}

pub fn write_data_set_to_file(file_path: &str, data_set: &Data) -> Result<(), Box<dyn Error>>{
    let json = serde_json::to_string_pretty(data_set).unwrap();
    std::fs::write(file_path, json)?;

    Ok(())
}
pub fn write_loss_set_to_file(file_path: &str, loss_set: &Loss) -> Result<(), Box<dyn Error>>{
    let json = serde_json::to_string_pretty(loss_set).unwrap();
    std::fs::write(file_path, json)?;

    Ok(())
}
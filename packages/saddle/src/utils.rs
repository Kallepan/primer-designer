use std::ops::{Bound, RangeBounds};

pub struct SubsequenceInfo {
    pub seq: String,
    pub start_index: usize,
    pub end_index: usize,
}

pub trait PrimerUtils {
    /*
        Returns a substring of the string, between the first argument (inclusive) and the second argument (exclusive).
    */
    fn substring(&self, start: usize, len: usize) -> &str;
    fn slice(&self, range: impl RangeBounds<usize>) -> &str;
    fn get_all_substrings_between(&self, min_size: usize, max_size: usize) -> Vec<SubsequenceInfo>;
}

impl PrimerUtils for str {
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
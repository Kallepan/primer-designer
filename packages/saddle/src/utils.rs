use std::ops::{Bound, RangeBounds};

pub fn combinations(start: usize, end: usize, min_size: usize, max_size: usize) -> Vec<[usize; 2]> {
    /*
        Calculates all possible combinations of all numbers between start and end 
        with max_size being the maximum size of the combination and min_size being the minimum size of the combination.
        EG:
        start = 0
        end = 3
        min_size = 1
        max_size = 2
        -> [[0, 1], [0, 2], [1, 2], [1, 3], [2, 3]]

        Usage: combinations(0, len(str)+1, 4, 8)
    */

    let mut combinations = Vec::new();
    let mut first_pointer = start;
    let mut second_pointer = first_pointer + min_size;

    loop {
        if first_pointer > end { break; }
        if second_pointer > end { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }
        if second_pointer - first_pointer > max_size { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }

        combinations.push([first_pointer, second_pointer]);
        second_pointer += 1;
    }

    combinations
}

pub trait StringUtils {
    /*
        Returns a substring of the string, between the first argument (inclusive) and the second argument (exclusive).
    */
    fn substring(&self, start: usize, len: usize) -> &str;
    fn slice(&self, range: impl RangeBounds<usize>) -> &str;
    fn get_all_substrings_between(&self, min_size: usize, max_size: usize) -> Vec<String>;
}

impl StringUtils for str {
    fn get_all_substrings_between(&self, min_size: usize, max_size: usize) -> Vec<String> {
        let mut substrings = Vec::new();
        
        let mut first_pointer = 0;
        let end = self.len();
        let mut second_pointer = first_pointer + min_size;
    
        loop {
            if first_pointer > end { break; }
            if second_pointer > end { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }
            if second_pointer - first_pointer > max_size { first_pointer += 1; second_pointer = first_pointer+min_size; continue; }
    
            substrings.push(self.slice(first_pointer..second_pointer).to_string());
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
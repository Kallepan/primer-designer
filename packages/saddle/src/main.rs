use std::sync::{ Mutex };
use once_cell::sync::Lazy;

mod saddle;
mod json;
mod errors;
mod utils;
mod types;
mod parser;

static DEBUG: Lazy<Mutex<bool>> = Lazy::new(|| Mutex::new(false));

pub fn is_debug() -> bool {
    *DEBUG.lock().unwrap()
}

fn main() {
    let args = parser::get_args();
    *DEBUG.lock().unwrap() = args.debug;

    // Global Debug variable

    let max_iterations = args.max_iterations.unwrap_or(100);
    let min_subsequence_size = args.min_subsequence_size.unwrap_or(4);
    let max_subsequence_size = args.max_subsequence_size.unwrap_or(8);
    let loss_acceptance_probability = args.loss_acceptance_probability.unwrap_or(0.1);

    saddle::run(
        &args.input_file, 
        &args.output_file,
        max_iterations,
        min_subsequence_size,
        max_subsequence_size,
        loss_acceptance_probability)
}

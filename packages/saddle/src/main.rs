mod saddle;
mod json;
mod errors;
mod utils;
mod parser;
mod logger;

use log;
use std::time;

fn init(debug: bool) {
    let _ = match log::set_logger(&logger::CONSOLE_LOGGER) {
        Ok(_) => {},
        Err(e) => {
            println!("Failed to set logger: {}", e);
            std::process::exit(1);
        }
    };
    log::set_max_level(if debug { log::LevelFilter::Debug } else { log::LevelFilter::Info });
}

fn main() {
    // Setup Parser
    let args = parser::get_args();

    init(args.debug);

    let min_subsequence_size = args.min_subsequence_size.unwrap_or(4);
    let max_subsequence_size = args.max_subsequence_size.unwrap_or(8);
    let max_iterations = args.max_iterations.unwrap_or(40_000);
    let amplicons_weight = args.amplicons_weight.unwrap_or(5.0);
    let primers_weight = args.primers_weight.unwrap_or(10.0);
    let num_primers_to_replace = args.num_primers_to_replace.unwrap_or(1);

    let start = time::Instant::now();
    saddle::run(
        &args.input_file,
        &args.output_file_set,
        &args.output_file_loss,
        min_subsequence_size,
        max_subsequence_size,
        num_primers_to_replace,
        max_iterations,
        amplicons_weight,
        primers_weight,
    );
    println!("Script used {:?}", start.elapsed()); 
}

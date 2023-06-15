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

    let start = time::Instant::now();
    saddle::run(
        &args.input_file, 
        &args.output_folder,
        min_subsequence_size,
        max_subsequence_size,
    );
    println!("Script used {:?}", start.elapsed()); 
}
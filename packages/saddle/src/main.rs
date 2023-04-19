mod saddle;
mod json;
mod errors;
use std::env;
use std::process;

fn get_args() -> Vec<String> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("Usage: saddle <file>");
        process::exit(1);
    }

    return args
}


fn main() {
    let args = get_args();

    let input_file = &args[1];

    println!("Using File: {}", &input_file);

    saddle::run(&input_file)
}

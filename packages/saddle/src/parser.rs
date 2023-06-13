use clap::Parser;

#[derive(Parser)]
#[clap(
    name = "SADDLE",
    version = "1.0.0",
    author = "Kandeepan K.",
    about = "A tool for Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE)."
)]
pub struct Opts {
    // Input file
    #[clap(short, long)]
    pub input_file: String,

    // Output folder
    #[clap(short, long)]
    pub output_folder: String,

    // Debug
    #[clap(short, long)]
    pub debug: bool,

    // SADDLE PARAMS
    #[clap(long)]
    pub max_iterations: Option<usize>,
    #[clap(long)]
    pub min_subsequence_size: Option<usize>,
    #[clap(long)]
    pub max_subsequence_size: Option<usize>,
}

pub fn get_args() -> Opts {
    Opts::parse()
}
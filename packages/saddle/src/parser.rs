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

    // Output files
    #[clap(long)]
    pub output_file_set: String,
    #[clap(long)]
    pub output_file_loss: String,
    
    // Debug
    #[clap(short, long)]
    pub debug: bool,

    // SADDLE PARAMS
    #[clap(long)]
    pub min_subsequence_size: Option<usize>,
    #[clap(long)]
    pub max_subsequence_size: Option<usize>,
    #[clap(long)]
    pub max_iterations: Option<usize>,
    #[clap(long)]
    pub num_primers_to_replace: Option<usize>,
    #[clap(long)]
    pub amplicons_weight: Option<f64>,
    #[clap(long)]
    pub primers_weight: Option<f64>,
}

pub fn get_args() -> Opts {
    Opts::parse()
}
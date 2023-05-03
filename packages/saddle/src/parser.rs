use clap::Parser;

#[derive(Parser)]
#[clap(
    name = "SADDLE",
    version = "1.0.0",
    author = "Kandeepan K.",
    about = "A tool for finding saddle points in time series data."
)]
pub struct Opts {
    /// Input file
    #[clap(short, long)]
    pub input_file: String,

    /// Output file
    #[clap(short, long)]
    pub output_file: String,

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
    #[clap(long)]
    pub simulated_annealing_stop_generation: Option<usize>,

}

pub fn get_args() -> Opts {
    Opts::parse()
}
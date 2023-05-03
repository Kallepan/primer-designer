use clap::Parser;

#[derive(Parser)]
#[clap(
    name = "saddle",
    version = "0.1.0",
    author = "J. R. Ullmann ",
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
    pub loss_acceptance_probability: Option<f64>,

}

pub fn get_args() -> Opts {
    Opts::parse()
}
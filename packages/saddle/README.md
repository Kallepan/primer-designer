# SADDLE (Simulated Annealing Design using Dimer Likelihood Estimation)

The Rust Primer Design Package, utilizing the SADDLE (Simulated Annealing Design using Dimer Likelihood Estimation) algorithm, is a high-performance tool for generating optimized primer sets from a sequence pool. Written in Rust to leverage its performance benefits, the package takes a JSON file with specified formatting and applies SADDLE's intelligent search and optimization techniques. It returns a JSON file containing the most optimal primer sequences, ensuring high amplification efficiency, specificity, and target coverage. Researchers and molecular biologists can rely on this package for efficient and accurate primer design in various applications like PCR and sequencing.

## Compilation

- The necessary configurations are present in the `Cargo.toml` file.

```bash
cargo build --release
```

## Usage

```bash
./saddle \
    --input-file <input_file> \
    --output-folder <output_folder> \
    # Optional arguments, defaults are 4, 8 and 40_000 respectively
    --min-subsequence-size <min_subsequence_size> \ 
    --max-subsequence-size <max_subsequence_size> \
    --optimal-iterations <optimal_iterations>
```

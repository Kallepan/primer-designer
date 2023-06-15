use crate::{json};
use crate::json::{Set, Pool, AmpliconPrimerPair, Region};
use crate::utils::PrimerUtils;

use rand::random;
use std::collections::HashMap;
use std::path::Path;

use std::time::Instant; // TODO for benchmarking

fn calculate_loss(hash_map: &HashMap<String, f64>, amplicon_primer_pairs: &Vec<AmpliconPrimerPair>, subsequence_min_size: usize, subsequence_max_size: usize) -> f64 {
    /*
        - For every sequence in the hash table
        - find the reverse complement in all primers with the distance to 3' End
        - calculate the loss for every subsequence SUM(1 / (distance_of_revcomp + 1)) * 2^length_of_subsequence * 2^num_gc_of_subsequence * Hashvalue
            - if none are found -> loss = 0
    */
    fn calculate_distance(length: usize, subsequence_start_index: usize, subsequence_length: usize) -> f64 {
        1.0 / (length - (subsequence_start_index + subsequence_length) + 1) as f64
    }

    let mut loss = 0.0;

    let time = Instant::now();

    for primer_pair in amplicon_primer_pairs {
        let left_primer = &primer_pair.forward_primer.sequence;
        let right_primer = &primer_pair.reverse_primer.sequence;

        // Add self badness from each primer
        loss += &primer_pair.forward_primer.badness;
        loss += &primer_pair.reverse_primer.badness;

        for subsequence_info in left_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let rev_comp = calculate_reverse_complement(&subsequence_info.seq);
            let hash_value = match hash_map.get(&rev_comp) {
                None => continue,
                Some(value) => value,
            };
            
            loss += hash_value *
            calculate_distance(left_primer.len(), subsequence_info.start_index, subsequence_info.seq.len()) *
                2.0_f64.powi(subsequence_info.seq.len() as i32) * 2.0_f64.powi(subsequence_info.seq.num_gc() as i32);
        }

        for subsequence_info in right_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let rev_comp = calculate_reverse_complement(&subsequence_info.seq);
            let hash_value = match hash_map.get(&rev_comp) {
                None => continue,
                Some(value) => value,
            };
            
            loss += hash_value *
            calculate_distance(right_primer.len(), subsequence_info.start_index, subsequence_info.seq.len()) *
                2.0_f64.powi(subsequence_info.seq.len() as i32) * 2.0_f64.powi(subsequence_info.seq.num_gc() as i32);
       }
    }
    // TODO Debug
    println!("Loss: {}", loss);
    println!("Time for one loss: {:?}", time.elapsed());

    loss
}

fn replace_primer_in_set(primer_pool: &Pool, primer_set: &mut Vec<AmpliconPrimerPair>, hash_map: &mut HashMap<String, f64>, subsequence_min_size: usize, subsequence_max_size: usize) {
    /* 
        Replace one primer pair from the given set, by another pair from the list of proto-primers and update the hash map.
    */

    fn calculate_distance_score(primer_len: usize, end_index: usize) -> f64 {
        let distance = primer_len - end_index;
        let distance_score = 1.0 / (distance + 1) as f64;
        distance_score
    }

    fn remove_primer_from_hash_map(hash_map: &mut HashMap<String, f64>, sequence: &String, length: usize, subsequence_min_size: usize, subsequence_max_size: usize) {
        for subsequence_info in sequence.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(length, subsequence_info.end_index);
            
            let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
            *entry -= distance_score;
        }
    }
    fn add_primer_to_hash_map(hash_map: &mut HashMap<String, f64>, sequence: &String, length: usize, subsequence_min_size: usize, subsequence_max_size: usize) {
        for subsequence_info in sequence.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(length, subsequence_info.end_index);
            
            let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
            *entry += distance_score;
        }
    }

    let random_primer_pair_index = random::<usize>() % primer_set.len();
    let forward_primer = &primer_set[random_primer_pair_index].forward_primer;
    let reverse_primer = &primer_set[random_primer_pair_index].reverse_primer;
    let region_name = &primer_set[random_primer_pair_index].region_name;
    let amplicon_name = &primer_set[random_primer_pair_index].amplicon_name;
    remove_primer_from_hash_map(hash_map, &forward_primer.sequence, forward_primer.length, subsequence_min_size, subsequence_max_size);
    remove_primer_from_hash_map(hash_map, &reverse_primer.sequence, reverse_primer.length, subsequence_min_size, subsequence_max_size);

    // get random new primer from proto-primers
    let potential_primers = primer_pool.regions.iter().find(|region| region.region_name == *region_name).unwrap().amplicons.iter().find(|amplicon| amplicon.amplicon_name == *amplicon_name).unwrap();
    let random_int_forward = random::<usize>() % potential_primers.forward_primers.len();
    let random_int_reverse = random::<usize>() % potential_primers.reverse_primers.len();

    primer_set[random_primer_pair_index].forward_primer = potential_primers.forward_primers[random_int_forward].clone();
    primer_set[random_primer_pair_index].reverse_primer = potential_primers.reverse_primers[random_int_reverse].clone();

    // Populate HashMap with new values
    add_primer_to_hash_map(hash_map, &primer_set[random_primer_pair_index].forward_primer.sequence, primer_set[random_primer_pair_index].forward_primer.length, subsequence_min_size, subsequence_max_size);
    add_primer_to_hash_map(hash_map, &primer_set[random_primer_pair_index].reverse_primer.sequence, primer_set[random_primer_pair_index].reverse_primer.length, subsequence_min_size, subsequence_max_size);
}

fn pick_random_primer_set(pool: &Vec<Region>) -> Vec<AmpliconPrimerPair> {
    /*
    Pick a Set of primers from the given proto primers
    Here: each amplicon should always have a forward and reverse primer.
    Therefore, the simplest solution to keep track of it is to 
    generate a list of dictionaries, where each dictionary contains
    the left and right primer for each amplicon along with its name.
    */
        
    let mut picked_primer_pairs = Vec::new();
    
    for region in pool {
        for amplicon in &region.amplicons {
            let random_forward_selector = random::<usize>() % amplicon.forward_primers.len();
            let random_reverse_selector = random::<usize>() % amplicon.reverse_primers.len();

            let primer_set_entry =  AmpliconPrimerPair {
                region_name: region.region_name.clone(),
                amplicon_name: amplicon.amplicon_name.clone(),
                forward_primer: amplicon.forward_primers[random_forward_selector].clone(),
                reverse_primer: amplicon.reverse_primers[random_reverse_selector].clone(),
            };
            
            picked_primer_pairs.push(primer_set_entry);
        }
    }
    
   picked_primer_pairs
}
    
fn initialize_hash_map(hash_map: &mut HashMap<String, f64>, primer_set: &Vec<AmpliconPrimerPair>, subsequence_min_size: usize, subsequence_max_size: usize) {
    /*
    Populate the hash table with all distance values for every subsequence in every in the primer set.

    For 5' -> 3' sequences, the distance is calculated from the end of the sequence. 
    Since all left and right primer are given 5' -> 3', the distance score for the primer can be calculated as is.
    */

    fn calculate_distance_score(length: usize, end_index: usize) -> f64 {
        let distance = length - end_index;
        let distance_score = 1.0 / (distance + 1) as f64;
        distance_score
    }

    fn add_primer_to_hash_map(hash_map: &mut HashMap<String, f64>, sequence: &String, length: usize, subsequence_min_size: usize, subsequence_max_size: usize) {
        for subsequence_info in sequence.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(length, subsequence_info.end_index);
            
            let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
            *entry += distance_score;
        }
    }

    for primer_pair in primer_set.iter() {
        let forward_primer = &primer_pair.forward_primer;
        let reverse_primer = &primer_pair.reverse_primer;

        add_primer_to_hash_map(hash_map, &forward_primer.sequence, forward_primer.length, subsequence_min_size, subsequence_max_size);
        add_primer_to_hash_map(hash_map, &reverse_primer.sequence, reverse_primer.length, subsequence_min_size, subsequence_max_size);
    }
}

fn calculate_reverse_complement(primer: &str) -> String {
    let mut reverse_complement = String::new();
    
    for nucleotide in primer.chars().rev() {
        match nucleotide {
            'A' => reverse_complement.push('T'),
            'T' => reverse_complement.push('A'),
            'G' => reverse_complement.push('C'),
            'C' => reverse_complement.push('G'),
            'N' => reverse_complement.push('N'),
            _ => panic!("Invalid Nucleotide") // It should not panic, since we do not expect any other nucleotide
        }
    }
    reverse_complement
}

pub fn run(
    input_file_path: &String,
    output_folder_path: &String,
    max_iterations: usize, 
    subsequence_min_size: usize,
    subsequence_max_size: usize)
{
    
    let pool = match json::load_json_from_file(&input_file_path) {
        Ok(data) => data,
        Err(e) => panic!("Failed to parse JSON file. Error: {}", e)
    };
    let mut final_sets = Vec::new();
    let output_folder = Path::new(output_folder_path);

    let mut hash_map: HashMap<String, f64> = HashMap::new();
    let mut past_sets = Vec::new();

    // Simulated Annealing Parameters
    let number_of_primers_in_pool = pool.regions.iter().fold(0, |acc, region| acc + region.amplicons.iter().fold(0, |amp_acc: usize, amplicon| amp_acc + amplicon.forward_primers.len() + amplicon.reverse_primers.len()));
    let sa_temp_initial = (1000 + 10*number_of_primers_in_pool) as f64;
    let numsteps = 10.0 + number_of_primers_in_pool as f64/10.0;

    println!("Primers in Pool: {}", number_of_primers_in_pool);
    println!("Initial Temperature: {}", sa_temp_initial);
    println!("Number of Steps: {}", numsteps);

    // TODO Multi Threading
    let mut iteration = 1;
    let mut sa_temp = sa_temp_initial;

    let amplicon_primer_pairs = pick_random_primer_set(&pool.regions);
    initialize_hash_map(&mut hash_map, &amplicon_primer_pairs, subsequence_min_size, subsequence_max_size);

    let loss = calculate_loss(&hash_map, &amplicon_primer_pairs, subsequence_min_size, subsequence_max_size);

    let mut current_set = Set {
        amplicon_primer_pairs,
        loss,
    };
    while iteration <= max_iterations {
        println!("Pool: {}, Iteration: {}", pool.pool_id, iteration);
        /*
            Generate a temp set and calculate the loss by recalculating the hash map. Store the old hash map in case the temp set is not accepted.
            If the temp set is accepted, store the temp set and continue.
            If the temp set is not accepted, restore the old hash map, store the temp set, and continue.
        */
        let mut temp_set = current_set.clone();
        let old_hash_map = hash_map.clone();

        replace_primer_in_set(&pool, &mut temp_set.amplicon_primer_pairs, &mut hash_map, subsequence_min_size, subsequence_max_size);
        temp_set.loss = calculate_loss(&hash_map, &temp_set.amplicon_primer_pairs, subsequence_min_size, subsequence_max_size);
        
        // Simulated Annealing
        let accept;
        if temp_set.loss <= current_set.loss {
            // Better set -> accept
            accept = true;
        } else {
            // Worse set -> Calculate acceptance
            let acceptance_prob = ((current_set.loss - temp_set.loss) / sa_temp as f64 ).exp();
            let random_number: f64 = random::<f64>();
            println!("Acceptance Probability: {}, Random Number: {}", acceptance_prob, random_number);
            if random_number < acceptance_prob {
                accept = true;
            } else {
                accept = false;
            }
        }
        
        // Implementation of (non)-acceptment of current set
        if accept {
            past_sets.push(current_set);
            current_set = temp_set;
        } else {
            past_sets.push(temp_set);
            hash_map = old_hash_map;
        }

        // Update SA Temp and iteration index
        sa_temp -= f64::max(0.0, sa_temp - (sa_temp_initial / numsteps));
        iteration += 1;
    }
    
    final_sets.push(current_set);

    // Write set history to file
    let file_name = output_folder.join(format!("pool_{}_history.json", pool.pool_id));
    match json::write_set_to_file(file_name.to_str().unwrap(), past_sets) {
        Ok(_) => println!("Successfully wrote set history to file {}_history.json", pool.pool_id),
        Err(e) => println!("Failed to write to set history. Error: {}", e)
    }
    

    let file_name = output_folder.join("final_set.json");
    match json::write_set_to_file(file_name.to_str().unwrap(), final_sets) {
        Ok(_) => println!("Successfully wrote output to file {}.", file_name.to_str().unwrap()),
        Err(e) => println!("Failed to write to file. Error: {}", e)
    }
}
use crate::json;
use crate::utils::PrimerUtils;
use rand::random;
use std::collections::HashMap;

pub struct AmpliconPrimer {
    pub amplicon_name: String,
    pub left_primer: String,
    pub right_primer: String,
}

fn calculate_loss(has_table: HashMap<String, f64>, primer_set: Vec<AmpliconPrimer>) {
    /*
        Calculate the loss for the given primer set.
        - Iterate over every primer again.
        - For every subsequence in the primer, calculate the reverse complement.
        - Look up the reverse complement in the hash table.
        - Calculate the distance score for the subsequence.
        - Multiply distance with the 2^(number of GC)*2^(len).
        - Add the value to the loss.
        - Return the loss.
    */
    
    let mut loss = 0.0;

    for primer_pair in primer_set {
        let left_primer = primer_pair.left_primer;
        let right_primer = primer_pair.right_primer;
    }    
}

fn pick_random_primer_set(proto_primers: &json::Regions) -> Vec<AmpliconPrimer> {
    /*
    Pick a Set of primers from the given proto primers
    Here: each amplicon should always have a left and right primer.
    Therefore, the simplest solution to keep track of it is to 
    generate a list of dictionaries, where each dictionary contains
    the left and right primer for each amplicon along with its name.
    Like so:
    [
        {
            "amplicon_name": "amplicon_1",
            "left_primer": "primer_1",
            "right_primer": "primer_2"
        },
        {
            "amplicon_name": "amplicon_2",
            "left_primer": "primer_3",
            "right_primer": "primer_4"
        }
        ]
        */
        
        let mut primer_set = Vec::new();
        for region in &proto_primers.regions {
            for amplicon in &region.amplicons {
                let random_int_forward = random::<usize>() % amplicon.forward_primers.len();
                let random_int_reverse = random::<usize>() % amplicon.reverse_primers.len();
                
                let primer_set_entry =  AmpliconPrimer {
                    amplicon_name: amplicon.name.clone(),
                    left_primer: amplicon.forward_primers[random_int_forward].sequence.clone(),
                    right_primer: amplicon.reverse_primers[random_int_reverse].sequence.clone(),
                };
                
                primer_set.push(primer_set_entry);
            }            
        }
        
        println!("Primer Set: {:?}", primer_set[0].amplicon_name);
        println!("Primer Set: {:?}", primer_set[0].left_primer);
        println!("Primer Set: {:?}", primer_set[0].right_primer);
        
        primer_set
    }
    
fn calculate_hash_table(hash_table: &mut HashMap<String, f64>, primer_set: Vec<AmpliconPrimer>, subsequence_min_size: usize, subsequence_max_size: usize) {
    /*
    Populate the hash table with all distance values for every subsequence in every in the primer set.

    For 5' -> 3' sequences, the distance is calculated from the end of the sequence. 
    Since all left and right primer are given 5' -> 3', the distance score for the primer can be calculated as is.
    */

    fn calculate_distance_score(primer_len: usize, end_index: usize) -> f64 {
        let distance = primer_len - end_index;
        let distance_score = 1.0 / (distance + 1) as f64;

        distance_score
    }

    for primer_pair in primer_set {
        let left_primer = primer_pair.left_primer;
        let right_primer = primer_pair.right_primer;
        
        for subsequence_info in left_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(left_primer.len(), subsequence_info.end_index);

            let reverse_complement_sequence = calculate_reverse_complement(&subsequence_info.seq);
            
            let entry = hash_table.entry(reverse_complement_sequence).or_insert(0.0);
            *entry += distance_score;
        }
        
        for subsequence_info in right_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(right_primer.len(), subsequence_info.end_index);

            let reverse_complement_sequence = calculate_reverse_complement(&subsequence_info.seq);

            let entry = hash_table.entry(reverse_complement_sequence).or_insert(0.0);
            *entry += distance_score;
        }
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
            _ => panic!("Invalid Nucleotide")
        }
    }
    reverse_complement
}

pub fn run(input_file_path: &str) {
    let mut hash_table: HashMap<String, f64> = HashMap::new();
    
    let data = match json::load_json_from_file(&input_file_path) {
        Ok(data) => data,
        Err(e) => panic!("Failed to parse JSON file. Error: {}", e)
    };
    
    let mut current_primer_set = pick_random_primer_set(&data);
    
    calculate_hash_table(&mut hash_table, current_primer_set, 4, 8); // TODO: Hard coded value change to modifiable value
    // calculate_loss(hash_table, current_primer_set);
}
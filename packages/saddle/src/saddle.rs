use crate::json;
use crate::json::Regions;
use crate::utils::PrimerUtils;
use rand::random;
use std::collections::HashMap;

pub struct AmpliconPrimer {
    pub region: String,
    pub amplicon_name: String,
    pub left_primer: String,
    pub right_primer: String,
}

fn calculate_loss(hash_map: &HashMap<String, f64>, primer_set: &Vec<AmpliconPrimer>) {
    /*
        - For every in the hash table
        - find the reverse complement in all primers with the distance to 3' End
        - calculate the loss for every subsequence SUM(1 / (distance_of_revcomp + 1)) * 2^length_of_subsequence * 2^numGC_of_subsequence * Hashvalue
            - if none are found -> loss = 0
    */

    fn calc_distance(primer_length: usize, found_index: usize, sequence_length: usize) -> f64 {
        1.0 / (primer_length - (found_index + sequence_length) + 1) as f64
    }
    
    let mut loss = 0.0;

    for (sequence, hash_value) in hash_map {
        let mut inverse_distance_sum = 0.0;
        let reverse_complement = calculate_reverse_complement(&sequence);
        for primer_pair in primer_set {
            let left_primer = &primer_pair.left_primer;
            let right_primer = &primer_pair.right_primer;

            match left_primer.find(&reverse_complement) {
                None => {},
                Some(found_index) => inverse_distance_sum += calc_distance(left_primer.len(), found_index, sequence.len())
            };
            match right_primer.find(&reverse_complement) {
                None => {},
                Some(found_index) => inverse_distance_sum += calc_distance(right_primer.len(), found_index, sequence.len())
             };
        }
        let badness = inverse_distance_sum*hash_value*2.0_f64.powi(sequence.len() as i32)*2.0_f64.powi(sequence.numGC() as i32);
        loss += badness;
    }
    println!("Loss: {}", loss);
}

fn replace_primer(proto_primers: Regions, primer_set: &mut Vec<AmpliconPrimer>, hash_map: &mut HashMap<String, f64>, subsequence_min_size: usize, subsequence_max_size: usize) {
    /* 
        Replace one primer pair from the given set, by another pair from the list of proto-primers and update the hash map.
    */
    fn calculate_distance_score(primer_len: usize, end_index: usize) -> f64 {
        let distance = primer_len - end_index;
        let distance_score = 1.0 / (distance + 1) as f64;

        distance_score
    }

    let random_primer_pair_index = random::<usize>() % primer_set.len();
    let left_primer = &primer_set[random_primer_pair_index].left_primer;
    let right_primer = &primer_set[random_primer_pair_index].right_primer;
    let region_name = &primer_set[random_primer_pair_index].region;
    let amplicon_name = &primer_set[random_primer_pair_index].amplicon_name;
    for subsequence_info in left_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
        let distance_score = calculate_distance_score(left_primer.len(), subsequence_info.end_index);
        
        let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
        *entry -= distance_score;
    }

    for subsequence_info in right_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
        let distance_score = calculate_distance_score(right_primer.len(), subsequence_info.end_index);
        
        let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
        *entry -= distance_score;
    }

    // get random new primer from proto-primers
    let potential_primers = proto_primers.regions.iter().find(|region| region.region_name == *region_name).unwrap().amplicons.iter().find(|amplicon| amplicon.amplicon_name == *amplicon_name).unwrap();
    let random_int_forward = random::<usize>() % potential_primers.forward_primers.len();
    let random_int_reverse = random::<usize>() % potential_primers.reverse_primers.len();

    primer_set[random_primer_pair_index].left_primer = potential_primers.forward_primers[random_int_forward].primer_sequence.clone();
    primer_set[random_primer_pair_index].right_primer = potential_primers.reverse_primers[random_int_reverse].primer_sequence.clone();

    let left_primer = &primer_set[random_primer_pair_index].left_primer;
    let right_primer = &primer_set[random_primer_pair_index].right_primer;

    // Populate HashMap with new values
    for subsequence_info in left_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
        let distance_score = calculate_distance_score(left_primer.len(), subsequence_info.end_index);
        
        let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
        *entry += distance_score;
    }

    for subsequence_info in right_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
        let distance_score = calculate_distance_score(right_primer.len(), subsequence_info.end_index);
        
        let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
        *entry += distance_score;
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
            "region": "region_1",
            "amplicon_name": "amplicon_1",
            "left_primer": "primer_1",
            "right_primer": "primer_2"
        },
        {
            "region": "region_2",
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
                    region: region.region_name.clone(),
                    amplicon_name: amplicon.amplicon_name.clone(),
                    left_primer: amplicon.forward_primers[random_int_forward].primer_sequence.clone(),
                    right_primer: amplicon.reverse_primers[random_int_reverse].primer_sequence.clone(),
                };
                
                primer_set.push(primer_set_entry);
            }
        }
        primer_set
    }
    
fn calculate_hash_map(hash_map: &mut HashMap<String, f64>, primer_set: &Vec<AmpliconPrimer>, subsequence_min_size: usize, subsequence_max_size: usize) {
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
        let left_primer = &primer_pair.left_primer;
        let right_primer = &primer_pair.right_primer;
        
        for subsequence_info in left_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(left_primer.len(), subsequence_info.end_index);
            
            let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
            *entry += distance_score;
        }
        
        for subsequence_info in right_primer.get_all_substrings_between(subsequence_min_size, subsequence_max_size) {
            let distance_score = calculate_distance_score(right_primer.len(), subsequence_info.end_index);

            let entry = hash_map.entry(subsequence_info.seq).or_insert(0.0);
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
    let mut hash_map: HashMap<String, f64> = HashMap::new();
    hash_map.insert("GTCGTT".to_string(), 1.58);
    let proto_primers = match json::load_json_from_file(&input_file_path) {
        Ok(data) => data,
        Err(e) => panic!("Failed to parse JSON file. Error: {}", e)
    };
    
    // let mut current_primer_set = pick_random_primer_set(&data);
    
    let mut test_set = Vec::<AmpliconPrimer>::new();

    test_set.push(
        AmpliconPrimer { 
            region: "test".to_string(),
            amplicon_name: "p1".to_string(),
            left_primer: "GATGCTCGATGCATCGGTCGTT".to_string(), 
            right_primer: "TGAATTCTCTACTCTGCAACGACCC".to_string() 
        }
    );
    test_set.push(
        AmpliconPrimer {
            region: "test".to_string(),
            amplicon_name: "p2".to_string(),
            left_primer: "CCGTTTACGACTGTTGTCGTTTT".to_string(),
            right_primer: "CGATCGCGACGTGCCGAACGACA".to_string(),
        }
    );
    test_set.push(
        AmpliconPrimer {
            region: "test".to_string(),
            amplicon_name: "p3".to_string(),
            left_primer: "ACTGAGTACGTAAGTCGTTAGC".to_string(), 
            right_primer: "CGACGCGATCATCCGATAACGAC".to_string() 
        }
    );

    let mut current_set = pick_random_primer_set(&proto_primers);
    calculate_hash_map(&mut hash_map, &test_set, 4, 8); // TODO: Hard coded min and max subsequence value change to modifiable value
    
    calculate_loss(&hash_map, &test_set);
}
use crate::json;
use rand::random;
use std::collections::HashMap;

pub struct AmpliconPrimer {
    pub amplicon_name: String,
    pub left_primer: String,
    pub right_primer: String,
}

fn pick_random_primer_set(proto_primers: &json::ProtoPrimers) -> Vec<AmpliconPrimer> {
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

    for amplicon in &proto_primers.amplicons {
        let random_int_forward = random::<usize>() % amplicon.forward_primers.len();
        let random_int_reverse = random::<usize>() % amplicon.reverse_primers.len();

        let primer_set_entry =  AmpliconPrimer {
            amplicon_name: amplicon.id.clone(),
            left_primer: amplicon.forward_primers[random_int_forward].sequence.clone(),
            right_primer: amplicon.reverse_primers[random_int_reverse].sequence.clone(),
        };

        primer_set.push(primer_set_entry);
    }

    println!("Primer Set: {:?}", primer_set[0].amplicon_name);
    println!("Primer Set: {:?}", primer_set[0].left_primer);
    println!("Primer Set: {:?}", primer_set[0].right_primer);

    primer_set
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
    let mut rev_comp_hash_table: HashMap<String, f64> = HashMap::new();

    let data = match json::load_json_from_file(&input_file_path) {
        Ok(data) => data,
        Err(e) => panic!("Failed to parse JSON file. Error: {}", e)
    };

    let mut current_primer_set = pick_random_primer_set(&data);

}
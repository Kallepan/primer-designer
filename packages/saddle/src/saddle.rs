use crate::json;
use rand::random;

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

fn calculate_badness_for_primerset(primer_set: &Vec<AmpliconPrimer>) -> f32 {
    /*
        Calculate the badness of a primer set.
        The badness is a float value, where a lower value means a better primer set.
        The badness is calculated by summing up the badness of each primer paired against each other in the set.
        The badness of a given primer pair is calculated by the formula:

    */

    let mut badness = 0.0;


    badness
}
pub fn run(input_file_path: &str) {
    let data = match json::load_json_from_file(&input_file_path) {
        Ok(data) => data,
        Err(e) => panic!("Failed to parse JSON file. Error: {}", e)
    };

    let mut primer_set = pick_random_primer_set(&data);

}
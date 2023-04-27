use crate::types::{Regions, Pool};

use std::error::Error;

pub fn load_json_from_file(input_file: &str) -> Result<Regions, Box<dyn Error>> {
    let file_content = std::fs::read_to_string(&input_file)?;
    let proto_primers: Regions = serde_json::from_str(&file_content)?;

    Ok(proto_primers)
}

pub fn extract_pools_from_json(regions_result: &Regions) -> Result<Vec<Pool>, Box<dyn Error>> {
    let mut first_pool = Pool{
        amplicon_primers: Vec::new(),
        pool_id: 1,
    };
    let mut second_pool = Pool{
        amplicon_primers: Vec::new(),
        pool_id: 2,
    };

    for region in &regions_result.regions {
        for amplicon in &region.amplicons {
            let mut forward_primers: Vec<String> = Vec::new();
            let mut reverse_primers: Vec<String> = Vec::new();

            for primer in &amplicon.forward_primers {
                forward_primers.push(primer.primer_sequence.clone());
            }

            for primer in &amplicon.reverse_primers {
                reverse_primers.push(primer.primer_sequence.clone());
            }

            let amplicon_primer = crate::types::AmpliconPrimers {
                amplicon_name: amplicon.amplicon_name.clone(),
                region_name: region.region_name.clone(),
                forward_primers: forward_primers,
                reverse_primers: reverse_primers,
            };

            if amplicon.pool == 1 {
                first_pool.amplicon_primers.push(amplicon_primer);
            } else if amplicon.pool == 2 {
                second_pool.amplicon_primers.push(amplicon_primer);
            }
        }
    }

    Ok(vec![first_pool, second_pool])
}
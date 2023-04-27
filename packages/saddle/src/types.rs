use serde::{Deserialize, Serialize};

// json object structs
#[derive(Deserialize, Serialize)]
pub struct Primer {
    pub primer_sequence: String,
    pub tm: f32,
    pub gc_percent: f32,
    pub hairpin_th: f32,
    pub primer_length: usize,
    pub amplicon_length: usize,
    pub amplicon_sequence: String,
}

#[derive(Deserialize, Serialize)]
pub struct Amplicon {
    pub amplicon_name: String,
    pub forward_primers: Vec<Primer>,
    pub reverse_primers: Vec<Primer>,
    pub pool: usize
}

#[derive(Deserialize, Serialize)]
pub struct Region {
    pub amplicons: Vec<Amplicon>,
    pub region_name: String, 
}

#[derive(Deserialize, Serialize)]
pub struct Regions {
    pub sequence_id: String,
    pub regions: Vec<Region>,
}


// Private Storate Struct
pub struct AmpliconPrimers {
    pub amplicon_name: String,
    pub region_name: String,
    pub forward_primers: Vec<String>,
    pub reverse_primers: Vec<String>,
}

pub struct Pool {
    pub amplicon_primers: Vec<AmpliconPrimers>,
    pub pool_id: usize,
}

// SADDLE Set structs
pub struct AmpliconPrimerPair {
    pub region_name: String,
    pub amplicon_name: String,
    pub forward_primer: String,
    pub reverse_primer: String,
}

pub struct Set {
    pub amplicon_primer_pairs: Vec<AmpliconPrimerPair>,
    pub loss: f32,
}
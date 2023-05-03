use serde::{Deserialize, Serialize};
// json object structs
#[derive(Deserialize, Serialize, Clone)]
pub struct Primer {
    pub primer_sequence: String,
    pub tm: f32,
    pub gc_percent: f32,
    pub hairpin_th: f32,
    pub primer_length: usize,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct Amplicon {
    pub amplicon_name: String,
    pub forward_primers: Vec<Primer>,
    pub reverse_primers: Vec<Primer>
}

#[derive(Deserialize, Serialize)]
pub struct Region {
    pub amplicons: Vec<Amplicon>,
    pub region_name: String, 
}

#[derive(Deserialize, Serialize)]
pub struct Pool {
    pub regions: Vec<Region>,
    pub pool_id: usize,
}

// SADDLE structs
#[derive(Clone)]
pub struct AmpliconPrimerPair {
    pub region_name: String,
    pub amplicon_name: String,
    pub forward_primer: Primer,
    pub reverse_primer: Primer,
}
#[derive(Clone)]
pub struct Set {
    pub amplicon_primer_pairs: Vec<AmpliconPrimerPair>,
    pub loss: f64,
}
// Contains all structs for JSON objects

use serde::{Deserialize, Serialize};
#[derive(Deserialize, Serialize, Clone)]
pub struct Primer {
    pub id: String,
    pub sequence: String,
    pub tm: f64,
    pub gc_percent: f64,
    pub hairpin_th: f64,
    pub length: usize,
    pub position: isize,
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
    pub pool_id: String,
}

// SADDLE structs
#[derive(Clone, Serialize)]
pub struct PrimerPair {
    pub region_name: String,
    pub amplicon_name: String,
    pub forward_primer: Primer,
    pub reverse_primer: Primer,
}
#[derive(Clone, Serialize)]
pub struct Data {
    pub primer_pairs: Vec<PrimerPair>,
    pub pool_id: String,
    pub loss: f64,
}

#[derive(Clone, Serialize)]
pub struct Loss {
    pub losses: Vec<f64>,
    pub pool_id: String,
}
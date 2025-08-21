//! This is the data structure required for the parameters yaml file

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct FilterParams {
    pub parameters: Option<HashMap<String, f32>>,
    pub ref_parameters: Option<HashMap<String, f32>>,
    pub snp_parameters: Option<HashMap<String, f32>>,
    pub indel_parameters: Option<HashMap<String, f32>>,
    pub minor_allele_params: Option<MinorAlleleParams>,
    pub fix_gt: Option<bool>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct MinorAlleleParams {
    pub threshold: i32,
    pub strand_bias: Option<f32>,
    pub min_frs: Option<f32>,
}

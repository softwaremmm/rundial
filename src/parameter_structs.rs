//! This is the data structure required for the parameters yaml file

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
enum FilterValue {
    Int(u32),
    Float(f32),
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
struct FilterParamsUnit {
    pub min_depth: Option<i32>,
    pub min_hq_depth: Option<i32>,
    pub min_qual: Option<f32>,
    pub min_strand_bias: Option<f32>,
    pub strand_mismatch_threshold: Option<i32>,
    pub min_mq: Option<i32>,
    pub min_frs: Option<f32>,
    pub het_threshold: Option<f32>,
    pub min_vbd: Option<f32>,
    pub min_idv: Option<i32>,
    pub min_imf: Option<f32>,
    pub invalid_indel: Option<bool>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct FilterParams {
    pub parameters: Option<HashMap<String, f32>>,
    pub ref_parameters: Option<HashMap<String, f32>>,
    pub snp_parameters: Option<HashMap<String, f32>>,
    pub indel_parameters: Option<HashMap<String, f32>>,
    pub fix_gt: Option<bool>,
}

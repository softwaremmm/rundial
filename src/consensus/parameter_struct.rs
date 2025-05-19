//! This is the data structure required for the parameters yaml file

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum HetOption {
    Mask, // Will set bases to Z
    Ref,  // will use ref if possible, else highest depth
    Alt,  // will use alt if possible, else highest depth
    Best, // will use highest depth allele, use first allele if tie
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ConsensusParams {
    pub skip_indels: bool, // Will skip any vcf row with indels
    pub mask: Option<String>,
    pub mask_missing_sites: bool, // If true set sites to N if not in vcf
    pub het_snp_option: HetOption,
    pub het_indel_option: HetOption,
    pub use_filters: bool, // If true will mask sites which have filters
    pub filter_ignore_list: Option<Vec<String>>, // If set will allow these filters
    pub het_pc_threshold: Option<f32>,
    pub minor_pop_threshold: Option<i32>, // min depth of non-GT allele to be considered a minor population
    pub support_minor_pop_threshold: Option<i32>,
    pub main_caller: Option<String>, // If set will use this Caller on records from main vcf
    pub support_caller: Option<String>, // If set will use this Caller on records from support vcf
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GenomeCreationReport {
    #[serde(rename = "Sequencing Quality")]
    pub sequencing_quality: SequencingQuality,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SequencingQuality {
    #[serde(rename = "Reference genome length")]
    pub genome_length: i32,
    #[serde(rename = "Null calls")]
    pub null_calls: i32,
    #[serde(rename = "Mixed calls")]
    pub mixed_calls: i32,
    #[serde(rename = "Fixed coverage")]
    pub fixed_coverage: f32,
    #[serde(rename = "Null Genotype calls")]
    pub null_genotype_calls: i32,
    #[serde(rename = "Filtered calls")]
    pub filtered_calls: i32,
    #[serde(rename = "Masked calls")]
    pub masked_calls: i32,
    #[serde(rename = "Deleted calls")]
    pub deleted_calls: i32,
}

//! This is the data structure required for the parameters yaml file

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct ConsensusParams {
    pub skip_indels: bool,        // Will skip any vcf row with indels
    pub mask_missing_sites: bool, // If true set sites to N if not in vcf
    pub use_filters: bool,        // If true will mask sites which have filters
    pub filter_ignore_list: Option<Vec<String>>, // If set will allow these filters

    pub minor_pop_thresholds: Option<MinorPopParams>, // requirements for non-GT allele to be considered a minor population
    pub remove_sub_minor_pops: bool, // If true will remove minor alleles which are below thresholds. Applies to main and support vcf
    pub call_snps_in_support: bool,  // If false will filter these rows
    pub remove_minor_pops_in_support: bool,

    pub main_caller: Option<String>, // If set will use this Caller on records from main vcf
    pub support_caller: Option<String>, // If set will use this Caller on records from support vcf
}

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct MinorPopParams {
    pub threshold: i32,
    pub strand_bias: Option<f32>,
    pub min_frs: Option<f32>,
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
    #[serde(rename = "Deleted calls")]
    pub deleted_calls: i32,
    #[serde(rename = "Fixed coverage")]
    pub fixed_coverage: f32,
    #[serde(rename = "Mixed calls")]
    pub mixed_calls: i32,
    #[serde(rename = "Mixed snps")]
    pub mixed_snps: i32,
    #[serde(rename = "Mixed snp clusters")]
    pub mixed_snp_clusters: i32,
}

use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;

use crate::vcf::variant_record::Genotype;
use crate::vcf::vcf_header::{FilterHeader, HeaderLine};
use crate::vcf::{RecordValue, VCFHeader, VCFReader, VCFWriter, VariantRecord};

use phf::phf_map;

use crate::parameter_structs::FilterParams;

const MIN_DP: &str = "MIN_DP";
const MIN_HQ_DP: &str = "MIN_HQ_DP";
const MIN_QUAL: &str = "MIN_QUAL";
const INVALID_INDEL: &str = "INVALID_INDEL";
const STRAND_BIAS: &str = "STRAND_BIAS";
const STRAND_MISMATCH: &str = "STRAND_MISMATCH";
const MIN_FRS: &str = "MIN_FRS";
const MIN_MQ: &str = "MIN_MQ";
const MIN_VDB: &str = "MIN_VDB";
const MIN_IDV: &str = "MIN_IDV";
const MIN_IMF: &str = "MIN_IMF";

static DESCRIPTIONS: phf::Map<&'static str, &'static str> = phf_map! {
    "MIN_DP" => "Basic read depth is less than ?",
    "MIN_HQ_DP" => "High quality read depth is less than ?",
    "MIN_QUAL" => "Quality is less than ?",
    "INVALID_INDEL" => "Indel record with no alternate alleles",
    "STRAND_BIAS" => "Strand bias. One strand is more than ? times the other",
    "STRAND_MISMATCH" => "Strand mismatch. Top alleles (within ? of max depth) on forward and reverse strands are different",
    "MIN_FRS" => "Fraction of reads supporting the main allele is less than ?",
    "MIN_MQ" => "Minimum mapping quality. MQ is less than ?",
    "MIN_VDB" => "Variant distance bias is less than ?",
    "MIN_IDV" => "Minimum number of raw reads supporting an indel is less than ?",
    "MIN_IMF" => "Maximum fraction of raw reads supporting an indel is less than ?",
};

static RELEVANT_INFO_TAGS: phf::Map<&'static str, &'static str> = phf_map! {
    "MIN_MQ" => "MQ",
    "MIN_VDB" => "VDB",
    "MIN_IDV" => "IDV",
    "MIN_IMF" => "IMF",
};

type FilterFunction = Box<dyn Fn(&VariantRecord) -> Option<String>>;

/// A struct that filters a VariantRecord based on provided parameters
struct Filterer {
    filters: Vec<FilterFunction>,
}

impl Filterer {
    /// Create a new Filterer from a set of parameters
    /// 
    /// The params hashmap should contain the flags as the key and the thresholds as the value
    pub fn create_from_params(params: &Option<HashMap<String, f32>>) -> Self {
        let mut filterer = Self {
            filters: Vec::new(),
        };
        let Some(params) = params else {
            return filterer;
        };

        for (flag, threshold) in params.clone().into_iter() {
            match flag.as_str() {
                MIN_DP => filterer
                    .filters
                    .push(Box::new(move |record| is_low_depth(record, threshold))),
                MIN_HQ_DP => filterer
                    .filters
                    .push(Box::new(move |record| is_low_hq_depth(record, threshold))),
                MIN_QUAL => filterer
                    .filters
                    .push(Box::new(move |record| is_low_qual(record, threshold))),
                INVALID_INDEL => filterer
                    .filters
                    .push(Box::new(move |record| is_invalid_indel(record, threshold))),
                STRAND_BIAS => filterer
                    .filters
                    .push(Box::new(move |record| is_strand_bias(record, threshold))),
                STRAND_MISMATCH => filterer
                    .filters
                    .push(Box::new(move |record| is_strand_mismatch(record, threshold))),
                MIN_FRS => filterer
                    .filters
                    .push(Box::new(move |record| is_low_support(record, threshold))),
                MIN_MQ | MIN_VDB | MIN_IDV | MIN_IMF => filterer
                    .filters
                    .push(Box::new(move |record| {
                        is_low_tag(record, threshold, RELEVANT_INFO_TAGS[&flag], &flag)
                    })),
                _ => {println!("Unknown flag: {}", flag);},
            }
        }

        filterer
    }

    pub fn filter(&self, record: &VariantRecord) -> HashSet<String> {
        self.filters
            .iter()
            .filter_map(|f| f(record))
            .collect::<HashSet<String>>()
    }
}

// Filter functions

/// Check if a record has a low value for a specific tag in INFO
fn is_low_tag(record: &VariantRecord, threshold: f32, tag: &str, flag: &str) -> Option<String> {
    if let Some(value) = record.info.get(tag) {
        match value {
            RecordValue::Float(val) => {
                if *val < threshold {
                    return Some(flag.to_string());
                }
            }
            RecordValue::Integer(val) => {
                if *val < threshold as i32 {
                    return Some(flag.to_string());
                }
            }
            _ => {}
        }
    }
    return None;
}

/// Check if a record has a low quality score.
fn is_low_qual(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let Some(quality_score) = record.qual {
        if quality_score < threshold {
            return Some(MIN_QUAL.to_string());
        }
    }
    return None;
}

/// Check if a record has a low depth.
fn is_low_depth(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let Some(depth) = record.depth() {
        if *depth < threshold as i32 {
            return Some(MIN_DP.to_string());
        }
    }
    return None;
}

/// Check if a record has a low high quality depth, based on AD or, ADF and ADR.
fn is_low_hq_depth(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let Some(depths) = record.allele_depths() {
        if depths.iter().sum::<i32>() < threshold as i32 {
            return Some(MIN_HQ_DP.to_string());
        }
    }
    return None;
}

/// Check if a record has low fraction of support for the main allele.
fn is_low_support(record: &VariantRecord, threshold: f32) -> Option<String> {
    let main_allele = record.main_allele() as usize;
    if let Some(depths) = record.allele_depths() {
        let total_depth = depths.iter().sum::<i32>();
        if total_depth == 0 {
            return None;
        }
        if (depths[main_allele] as f32 / total_depth as f32) < threshold {
            return Some(MIN_FRS.to_string());
        }
    }
    return None;
}


/// Check if a record claims to be an indel but has no alternate alleles.
fn is_invalid_indel(record: &VariantRecord, _threshold: f32) -> Option<String> {
    let has_indel_label =
        record.info.contains_key("INDEL") || record.ref_bases.chars().count() != 1;
    if has_indel_label && record.alt.is_empty() {
        return Some(INVALID_INDEL.to_string());
    }
    return None;
}

/// Check if a record has many more reads on one strand than the other.
fn is_strand_bias(record: &VariantRecord, threshold: f32) -> Option<String> {
    let main_allele = record.main_allele();
    if main_allele == -1 {
        return None;
    }
    let main_allele = main_allele as usize;
    if let Some((forward, reverse)) = record.strand_depths() {
        let forward_depth = max(1, forward[main_allele]) as f32;
        let reverse_depth = max(1, reverse[main_allele]) as f32;

        if forward_depth / reverse_depth > threshold || reverse_depth / forward_depth > threshold {
            return Some(STRAND_BIAS.to_string());
        }
    }
    return None;
}

/// Check if a record has different top alleles on forward and reverse strands.
/// 
/// Top alleles are defined as the alleles with depth within `threshold` of the max depth
/// on each strand.
fn is_strand_mismatch(record: &VariantRecord, threshold: f32) -> Option<String> {
    // This checks if the two strands support different alleles
    if let Some((forward, reverse)) = record.strand_depths() {
        // unwrap should be fine as ref should always be present
        let max_forward: &i32 = forward.iter().max().unwrap();
        let max_reverse: &i32 = reverse.iter().max().unwrap();

        let top_forward_indexes: HashSet<usize> = forward
            .iter()
            .enumerate()
            .filter(|(_, &x)| x >= max_forward - threshold as i32)
            .map(|(i, _)| i)
            .collect();
        let top_reverse_indexes: HashSet<usize> = reverse
            .iter()
            .enumerate()
            .filter(|(_, &x)| x >= max_reverse - threshold as i32)
            .map(|(i, _)| i)
            .collect();

        // check if there is any overlap between the two lists
        if top_forward_indexes
            .intersection(&top_reverse_indexes)
            .count()
            == 0
        {
            return Some(STRAND_MISMATCH.to_string());
        }
    }
    return None;
}

// End of filter functions

/// Set the genotype to be the allele with the highest depth (x/x).
fn set_gt_to_highest_depth(record: &mut VariantRecord) {
    if record.alt.is_empty() {
        return;
    }

    if let Some(depths) = record.allele_depths() {
        let max_depth = depths.iter().max().unwrap();
        let max_index = depths.iter().position(|x| x == max_depth).unwrap();

        record.set_genotype(Genotype {
            allele1: max_index as i32,
            allele2: max_index as i32,
        });
    }
}

fn add_filter_to_header(header: &mut VCFHeader, key: &str, desc: &str) {
    let filter = FilterHeader {
        id: key.to_string(),
        desc: desc.to_string(),
    };
    header.add_header_line(HeaderLine::Filter(filter));
}

/// Add all filters to the header.
/// 
/// Filters which are in multiple parameter sets will have their thresholds listed in the description.
fn add_filters_to_header(header: &mut VCFHeader, params: &FilterParams) {
    let all_params: Vec<HashMap<String, f32>> = vec![
        params.parameters.clone(),
        params.ref_parameters.clone(),
        params.snp_parameters.clone(),
        params.indel_parameters.clone(),
    ]
    .into_iter()
    .flatten()
    .collect();

    let all_keys = all_params
        .iter()
        .flat_map(|x| x.keys())
        .collect::<HashSet<&String>>();

    for key in all_keys {
        let thresholds: String = all_params
            .iter()
            .filter(|x| x.contains_key(key))
            .map(|x| x[key].to_string())
            .collect::<Vec<String>>()
            .join(", ");

        let desc: String = DESCRIPTIONS
            .get(key)
            .map(|&desc| desc.to_string())
            .unwrap_or_else(|| format!("{} - thresholds: ?", key))
            .replace('?', &thresholds);

        add_filter_to_header(header, key, &desc);
    }
}

/// Filter a VCF file based on a set of parameters.
pub fn filter_vcf(
    in_vcf: PathBuf,
    out_vcf: PathBuf,
    params: PathBuf,
    overwrite: bool,
    verbose: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let params: FilterParams = serde_yaml::from_reader(File::open(params)?)?;

    let vcf_reader = VCFReader::new(BufReader::new(File::open(in_vcf)?))?;
    let header = vcf_reader.header();

    let mut new_header = header.clone();
    add_filters_to_header(&mut new_header, &params);
    let mut vcf_writer = VCFWriter::new(BufWriter::new(File::create(out_vcf)?), &new_header)?;

    let std_filterer = Filterer::create_from_params(&params.parameters);
    let ref_filterer = Filterer::create_from_params(&params.ref_parameters);
    let snp_filterer = Filterer::create_from_params(&params.snp_parameters);
    let indel_filterer = Filterer::create_from_params(&params.indel_parameters);
    let fix_gt: bool = params.fix_gt.unwrap_or(false);
    
    if verbose {
        println!("Filtering VCF file");
    }
    let mut filter_counter: HashMap<String, i32> = HashMap::new();
    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 {
            println!("Processed {} records", count);
        }
        let mut record = record?;
        if fix_gt {
            set_gt_to_highest_depth(&mut record);
        }

        let mut new_filters = std_filterer.filter(&record);

        if record.is_indel() {
            new_filters.extend(indel_filterer.filter(&record));
        } else if record.is_snp() && ! record.genotype().map(|x| x.is_hom_ref()).unwrap_or(true) {
            new_filters.extend(snp_filterer.filter(&record));
        } else {
            new_filters.extend(ref_filterer.filter(&record));
        }

        if new_filters.contains(INVALID_INDEL) {
            continue;
        }

        let mut new_filters: Vec<String> = new_filters.into_iter().collect();
        new_filters.sort();

        if verbose {
            for filter in new_filters.clone() {
                let count = filter_counter.entry(filter).or_insert(0);
                *count += 1;
            }
        }

        if overwrite {
            record.filter = new_filters;
        } else {
            record.filter.extend(new_filters);
        }
        vcf_writer.write_record(&record)?;
    }

    if verbose {
        println!("Filter counts:");
        for (filter, count) in filter_counter {
            println!("{}: {}", filter, count);
        }
    }

    Ok(())
}

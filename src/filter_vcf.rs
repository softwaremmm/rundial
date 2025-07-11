use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::PathBuf;

use crate::vcf::vcf_header::{FilterHeader, HeaderLine};
use crate::vcf::{Genotype, RecordValue, VCFHeader, VCFReader, VCFWriter, VariantRecord};

use phf::phf_map;

pub mod parameter_struct;
pub use parameter_struct::FilterParams;

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
const MIN_AF: &str = "MIN_AF";

static DESCRIPTIONS: phf::Map<&'static str, &'static str> = phf_map! {
    "MIN_DP" => "Basic read depth is less than ?",
    "MIN_HQ_DP" => "High quality read depth is less than ?",
    "MIN_QUAL" => "Quality is less than ?",
    "INVALID_INDEL" => "Indel record with no alternate alleles",
    "STRAND_BIAS" => "Strand bias. One strand is more than ? times the other",
    "STRAND_MISMATCH" => "Strand mismatch. Top alleles (defined as those within ? of max depth) on forward and reverse strands are different",
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

        // Always filter out invalid indels
        // They are more of a quirk in BCFTools
        filterer
            .filters
            .push(Box::new(move |record| is_invalid_indel(record, 0.0)));

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
                STRAND_BIAS => filterer
                    .filters
                    .push(Box::new(move |record| is_strand_bias(record, threshold))),
                STRAND_MISMATCH => filterer.filters.push(Box::new(move |record| {
                    is_strand_mismatch(record, threshold)
                })),
                MIN_FRS => filterer
                    .filters
                    .push(Box::new(move |record| is_low_support(record, threshold))),
                MIN_AF => filterer.filters.push(Box::new(move |record| {
                    is_low_allele_frequency(record, threshold)
                })),
                MIN_MQ | MIN_VDB | MIN_IDV | MIN_IMF => {
                    filterer.filters.push(Box::new(move |record| {
                        is_low_tag(record, threshold, RELEVANT_INFO_TAGS[&flag], &flag)
                    }))
                }
                _ => {
                    println!("Unknown flag: {}", flag);
                }
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
                if (*val as f32) < threshold {
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
        if (*depth as f32) < threshold {
            return Some(MIN_DP.to_string());
        }
    }
    return None;
}

/// Check if a record has a low high quality depth, based on AD or, ADF and ADR.
fn is_low_hq_depth(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let Some(depths) = record.allele_depths() {
        if (depths.iter().sum::<i32>() as f32) < threshold {
            return Some(MIN_HQ_DP.to_string());
        }
    }
    return None;
}

/// Check if a record has low fraction of support for the main allele.
fn is_low_support(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let (Some(main_allele), Some(depths)) = (record.main_allele(), record.allele_depths()) {
        let main_allele = main_allele as usize;
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

/// Check if a record has low fraction of support for the main allele compared to overall depth.
/// Used in clair3.
fn is_low_allele_frequency(record: &VariantRecord, threshold: f32) -> Option<String> {
    if let (Some(main_allele), Some(depths), &Some(total_depth)) =
        (record.main_allele(), record.allele_depths(), record.depth())
    {
        let main_allele = main_allele as usize;
        if total_depth == 0 {
            return None;
        }
        if (depths[main_allele] as f32 / total_depth as f32) < threshold {
            // This is a special case for clair3 so still returns MIN_FRS
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
    if let Some(main_allele) = record.main_allele() {
        let main_allele = main_allele as usize;
        if let Some((forward, reverse)) = record.strand_depths() {
            let forward_depth = max(1, forward[main_allele]) as f32;
            let reverse_depth = max(1, reverse[main_allele]) as f32;

            if forward_depth / reverse_depth >= threshold
                || reverse_depth / forward_depth >= threshold
            {
                return Some(STRAND_BIAS.to_string());
            }
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
    let threshold = threshold as i32;
    if let Some((forward, reverse)) = record.strand_depths() {
        // unwrap should be fine as ref should always be present
        let max_forward: &i32 = forward.iter().max().unwrap();
        let max_reverse: &i32 = reverse.iter().max().unwrap();

        let top_forward_indexes: HashSet<usize> = forward
            .iter()
            .enumerate()
            .filter(|(_, &x)| x >= max_forward - threshold)
            .map(|(i, _)| i)
            .collect();
        let top_reverse_indexes: HashSet<usize> = reverse
            .iter()
            .enumerate()
            .filter(|(_, &x)| x >= max_reverse - threshold)
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
        record.set_genotype(Genotype {
            allele1: 0,
            allele2: 0,
        });
        return;
    }

    if let Some(depths) = record.allele_depths() {
        let max_depth = depths.iter().max().unwrap();
        let max_index = depths.iter().rposition(|x| x == max_depth).unwrap();

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

    let mut all_keys = all_params
        .iter()
        .flat_map(|x| x.keys())
        .collect::<Vec<&String>>();
    all_keys.sort();

    for key in all_keys {
        if key == MIN_AF || key == MIN_FRS {
            continue;
        }
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

    // MIN_AL is a special kind of MIN_FRS so will get just one header
    let mut min_frs_descs = Vec::new();
    let min_frs_thresholds = all_params
        .iter()
        .filter(|x| x.contains_key(MIN_FRS))
        .map(|x| x[MIN_FRS].to_string())
        .collect::<Vec<String>>();
    if !min_frs_thresholds.is_empty() {
        min_frs_descs.push(format!(
            "{} (compared to other alleles)",
            min_frs_thresholds.join(", ")
        ));
    }

    let min_af_thresholds = all_params
        .iter()
        .filter(|x| x.contains_key(MIN_AF))
        .map(|x| x[MIN_AF].to_string())
        .collect::<Vec<String>>();
    if !min_af_thresholds.is_empty() {
        min_frs_descs.push(format!(
            "{} (compared to other overall depth)",
            min_af_thresholds.join(", ")
        ));
    }

    if !min_frs_descs.is_empty() {
        let min_frs_desc = DESCRIPTIONS
            .get(MIN_FRS)
            .map(|&desc| desc.to_string())
            .unwrap_or_else(|| format!("{} - thresholds: ?", MIN_FRS))
            .replace('?', &min_frs_descs.join(", "));
        add_filter_to_header(header, MIN_FRS, &min_frs_desc);
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
    let params: FilterParams = serde_yaml::from_reader(
        File::open(params).map_err(|e| format!("Failed to read params file. Error: {}", e))?,
    )?;

    let vcf_reader = VCFReader::from_path(in_vcf)?;

    let header = vcf_reader.header();

    let mut new_header = header.clone();
    add_filters_to_header(&mut new_header, &params);
    let mut vcf_writer = VCFWriter::to_path(out_vcf, new_header)?;

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
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {} records", count);
        }
        let mut record = record?;
        if fix_gt {
            set_gt_to_highest_depth(&mut record);
        }

        let mut new_filters = std_filterer.filter(&record);

        if record.is_indel() {
            new_filters.extend(indel_filterer.filter(&record));
        } else if record.is_snp() && !record.genotype().map(|x| x.is_hom_ref()).unwrap_or(true) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::variant_record::tests::{clair3_header, standard_header};

    #[test]
    fn test_is_low_tag() {
        let std_header = standard_header();
        let record_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7".to_string();
        let record: VariantRecord =
            VariantRecord::from_string(&std_header, &record_string).unwrap();

        assert!(is_low_tag(&record, 30.0, "DP", "FLAG") == Some("FLAG".to_string()));
        assert!(is_low_tag(&record, 28.1, "DP", "FLAG") == Some("FLAG".to_string()));
        assert!(is_low_tag(&record, 28.0, "DP", "FLAG").is_none());

        assert!(is_low_tag(&record, 28.0, "MISSING", "FLAG").is_none());
    }

    #[test]
    fn test_is_low_qual() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7",
        ).unwrap();

        assert!(is_low_qual(&record, 250.0) == Some(MIN_QUAL.to_string()));
        assert!(is_low_qual(&record, 244.0).is_none());

        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t.\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7",
        ).unwrap();
        assert!(is_low_qual(&record, 244.0).is_none());
    }

    #[test]
    fn test_is_low_depth() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7",
        ).unwrap();

        assert!(is_low_depth(&record, 30.0) == Some(MIN_DP.to_string()));
        assert!(is_low_depth(&record, 28.1) == Some(MIN_DP.to_string()));
        assert!(is_low_depth(&record, 28.0).is_none());

        let record: VariantRecord = VariantRecord::from_string(
            &clair3_header(),
            "ref\t1\tid\tT\t.\t244.589\tF1;F2\t.\tGT:DP:AD\t0/0:10:5",
        )
        .unwrap();
        assert!(is_low_depth(&record, 15.0) == Some(MIN_DP.to_string()));
        assert!(is_low_depth(&record, 7.0).is_none());
    }

    #[test]
    fn test_is_low_hq_depth() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7",
        ).unwrap();

        assert!(is_low_hq_depth(&record, 19.0) == Some(MIN_HQ_DP.to_string()));
        assert!(is_low_hq_depth(&record, 18.5) == Some(MIN_HQ_DP.to_string()));
        assert!(is_low_hq_depth(&record, 18.0).is_none());
    }

    #[test]
    fn test_is_low_support() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,5,0",
        ).unwrap();

        assert!(is_low_support(&record, 0.51) == Some(MIN_FRS.to_string()));
        assert!(is_low_support(&record, 0.50).is_none());
    }

    #[test]
    fn test_is_invalid_indel() {
        let std_header = standard_header();

        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\t.\t244.589\tF1;F2\tINDEL;DP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,5",
        ).unwrap();
        assert!(is_invalid_indel(&record, 0.0) == Some(INVALID_INDEL.to_string()));

        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tTT\t.\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,5",
        ).unwrap();
        assert!(is_invalid_indel(&record, 0.0) == Some(INVALID_INDEL.to_string()));

        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\t.\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,5",
        ).unwrap();
        assert!(is_invalid_indel(&record, 0.0).is_none());

        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tC,CG\t244.589\tF1;F2\tINDEL;DP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,5",
        ).unwrap();
        assert!(is_invalid_indel(&record, 0.0).is_none());
    }

    #[test]
    fn test_is_strand_bias() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=10,8,1;ADR=1,15,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,6,7",
        ).unwrap();

        assert!(is_strand_bias(&record, 10.0) == Some(STRAND_BIAS.to_string()));
        assert!(is_strand_bias(&record, 10.1).is_none());

        // Each strand is considered to have at least 1 read
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=10,8,1;ADR=0,15,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/0:5,6,7",
        ).unwrap();

        assert!(is_strand_bias(&record, 10.0) == Some(STRAND_BIAS.to_string()));
        assert!(is_strand_bias(&record, 10.1).is_none());
    }

    #[test]
    fn test_is_strand_mismatch() {
        let std_header = standard_header();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=10,8,1;ADR=1,15,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7",
        ).unwrap();

        assert!(is_strand_mismatch(&record, 1.0) == Some(STRAND_MISMATCH.to_string()));
        assert!(is_strand_mismatch(&record, 2.0).is_none());
    }

    #[test]
    fn test_set_gt_to_highest_depth() {
        let std_header = standard_header();

        let mut record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tMQ=53.0\tGT:AD\t0/1:5,6,7",
        )
        .unwrap();
        set_gt_to_highest_depth(&mut record);
        assert_eq!(
            record.genotype().unwrap(),
            &Genotype {
                allele1: 2,
                allele2: 2
            }
        );

        // With a draw will use the latest allele
        let mut record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tMQ=53.0\tGT:AD\t0/1:7,6,7",
        )
        .unwrap();
        set_gt_to_highest_depth(&mut record);
        assert_eq!(
            record.genotype().unwrap(),
            &Genotype {
                allele1: 2,
                allele2: 2
            }
        );

        // Will replace missing GT with 0/0
        let mut record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\t.\t244.589\tF1;F2\tMQ=53.0\tGT:AD\t.:0",
        )
        .unwrap();
        set_gt_to_highest_depth(&mut record);
        assert_eq!(
            record.genotype().unwrap(),
            &Genotype {
                allele1: 0,
                allele2: 0
            }
        );
        let mut record: VariantRecord = VariantRecord::from_string(
            &std_header,
            "ref\t1\tid\tT\tC\t244.589\tF1;F2\tMQ=53.0\tGT:AD\t.:0,0",
        )
        .unwrap();
        set_gt_to_highest_depth(&mut record);
        assert_eq!(
            record.genotype().unwrap(),
            &Genotype {
                allele1: 1,
                allele2: 1
            }
        );
    }

    #[test]
    fn test_add_filter_to_header() {
        let mut header = VCFHeader::new();
        add_filter_to_header(&mut header, "FLAG", "Description of FLAG");
        assert_eq!(header.filters().len(), 1);
        assert!(header.filters().contains_key("FLAG"));
    }

    #[test]
    fn test_add_filters_to_header() {
        let mut header = VCFHeader::new();
        let params = FilterParams {
            parameters: Some(HashMap::from([
                (MIN_DP.to_string(), 10.0),
                (MIN_QUAL.to_string(), 20.0),
            ])),
            ref_parameters: Some(HashMap::from([
                (MIN_DP.to_string(), 15.0),
                (MIN_QUAL.to_string(), 25.0),
            ])),
            snp_parameters: None,
            indel_parameters: None,
            fix_gt: None,
        };
        add_filters_to_header(&mut header, &params);
        assert_eq!(header.filters().len(), 2);
        assert!(header.filters().contains_key(MIN_DP));
        assert!(header.filters().contains_key(MIN_QUAL));
        println!("{:?}", header.filters());
        assert!(header
            .filters()
            .get(MIN_DP)
            .unwrap()
            .desc
            .contains("10, 15"));
        assert!(header
            .filters()
            .get(MIN_QUAL)
            .unwrap()
            .desc
            .contains("20, 25"));
    }
}

use core::panic;
use indexmap::IndexMap;
use ordered_float::OrderedFloat;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
};

use crate::vcf::vcf_header::{HeaderLine, HeaderNumber, HeaderType};
use crate::vcf::{Genotype, RecordValue, VCFHeader, VCFReader, VCFWriter, VariantRecord};

pub mod parameter_struct;
pub use parameter_struct::{ConsensusParams, GenomeCreationReport, SequencingQuality};
pub mod bed;
pub use bed::Bed;
mod classifier;
use classifier::{Change, Classification, Classifier, repeat_char};
mod contig_set;
use contig_set::ContigSet;
mod fasta_tools;
use fasta_tools::{read_fasta, save_fasta};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
type HashMapSet<T> = HashMap<String, HashSet<T>>;

const NULL: char = 'N';
const DELETED: char = '-';

const OVERLAP_FILTER: &str = "OVERLAP_BETTER_VARIANT";
const OVERLAP_FILTER_DESC: &str = "This variant overlaps with a higher priority variant.";
const SNP_IN_SUPPORT_VCF: &str = "SNP_IN_SUPPORT_VCF";
const SNP_IN_SUPPORT_VCF_DESC: &str =
    "This variant is a SNP in the support VCF, which should only be used for information of nulls.";
const CALLER: &str = "CALLER";
const CALLER_DESC: &str = "The variant caller that made the call.";
const MIXED: &str = "MIXED";
const MIXED_DESC: &str = "This variant is counted as a mixed call.";
const HET_GT: &str = "HET_GT";
const HET_GT_DESC: &str = "Het genotypes are treated like MIN_FRS.";
const REMOVED_MINORS: &str = "REMOVED_MINOR_ALLELES";
const REMOVED_MINORS_DESC: &str = "Minor alleles which were removed due to being below thresholds.";

/// Apply variant to consensus sequence, overwriting any existing changes
fn apply_variant_simple(
    chrom: &str,
    classification: &Classification,
    chrom_seq: &mut HashMap<String, Vec<char>>,
) {
    let pos = classification.pos;

    // select the correct chromosome sequence
    let chrom_seq: &mut Vec<char> = chrom_seq
        .get_mut(chrom)
        .expect("Chrom not found in consensus");

    match classification.change {
        Change::Null => {
            assert!(classification.ref_bases.len() == classification.new_bases.len());
            for (i, base) in classification.new_bases.chars().enumerate() {
                chrom_seq[pos + i] = base;
            }
        }
        Change::Ref => {
            for (i, base) in classification.ref_bases.chars().enumerate() {
                chrom_seq[pos + i] = base;
            }
        }
        Change::Snp => {
            assert!(classification.ref_bases.len() == 1 && classification.new_bases.len() == 1);
            chrom_seq[pos] = classification
                .new_bases
                .chars()
                .next()
                .expect("Snp new bases empty");
        }
        Change::Del | Change::ComplexDel | Change::Mnp => {
            // We assume that ref and alt have already been simplified
            assert!(classification.ref_bases.len() >= classification.new_bases.len());
            let replacement_str: String = classification.new_bases.clone()
                + &repeat_char(
                    DELETED,
                    classification.ref_bases.len() - classification.new_bases.len(),
                );
            for (i, base) in replacement_str.chars().enumerate() {
                chrom_seq[pos + i] = base;
            }
        }
        Change::Ins | Change::ComplexIns => {
            // We assume that ref and alt have already been simplified
            let new_bases = classification.new_bases.chars();
            let ref_length = classification.ref_bases.len();
            chrom_seq.splice(pos..pos + ref_length, new_bases);
        }
    }
}

// Score variants for marking overlaps. Overlaps are caused by indels
fn score_variant(
    record: &VariantRecord,
    classification: &Classification,
) -> (i32, OrderedFloat<f32>, i32) {
    let type_score = match (
        classification.is_filtered,
        record.genotype(),
        classification.has_indel_alleles,
    ) {
        (_, None, _) => 0,
        (_, Some(gt), _) if gt.is_null() => 0,
        (false, Some(gt), true) if !gt.is_hom_ref() => 10, // clear indel
        (true, Some(gt), true) if !gt.is_hom_ref() => 9,   // clear indel filtered
        (false, _, false) => 8,                            // snp/ref call
        (false, _, true) => 7,                             // minor indel
        (true, _, true) => 6,                              // minor indel filtered
        (true, _, false) => 5,                             // snp/ref call filtered
    };

    let qual_score = match record.qual {
        Some(f) => OrderedFloat(f),
        _ => OrderedFloat(0.0),
    };
    let depth_score = match record.depth() {
        Some(d) => *d,
        _ => 0,
    };
    return (type_score, qual_score, depth_score);
}

fn overlaps(this: &Classification, other: &Classification) -> bool {
    fn range(c: &Classification) -> (f32, f32) {
        // Position changes on float scale where x.0 is the middle of a base x
        // so snp at base 2 would be 2.0-2.0
        // insertion just after base 2 would be 2.5-2.5
        // deletion at base 2-3 would be 2.0-3.0
        return match c.change {
            Change::Null
            | Change::Ref
            | Change::Snp
            | Change::Mnp
            | Change::Del
            | Change::ComplexDel => (c.pos as f32, (c.pos + c.ref_bases.len()) as f32 - 1.0),
            Change::Ins | Change::ComplexIns => {
                (c.pos as f32 - 0.5, (c.pos + c.ref_bases.len()) as f32 - 0.5)
            }
        };
    }

    let (this_start, this_end) = range(this);
    let (other_start, other_end) = range(other);

    // normally the changes just need to not be touching,
    // but for ins and dels we need more of a buffer
    // as we don't want a del to occur just before or after a insertion.
    match (&this.change, &other.change) {
        (Change::Ins | Change::ComplexIns, Change::Del | Change::ComplexDel)
        | (Change::Del | Change::ComplexDel, Change::Ins | Change::ComplexIns) => {
            if this_end < other_start - 0.6 || this_start > other_end + 0.6 {
                return false;
            }
            return true;
        }
        _ => {
            if this_end < other_start || this_start > other_end {
                return false;
            }
            return true;
        }
    }
}

/// Adds "OverlapsBetterVariant" to Filter column to variants overlapping a higher priority variant
/// Priority:
/// 1. Clear indels (gt!=0/0)
/// 2. snps/refs calls
/// 3. Filtered indels
/// 4. Filtered snps/refs calls
fn mark_overlaps(records: &mut [(VariantRecord, Classification)], classifier: &Classifier) {
    records.sort_by_key(|(r, c)| (r.chrom.clone(), c.pos));

    let mut current_variants: Vec<(&mut VariantRecord, &mut Classification)> = Vec::new();

    for (record, classification) in records.iter_mut() {
        // Discard variants which are too far back to possibly overlap with current variant
        current_variants.retain(|(other_r, _)| {
            other_r.chrom == record.chrom
                && other_r.pos_idx() + other_r.ref_bases.len() >= record.pos_idx() - 1
        });

        // Check for clashes and filter accordingly
        let score = score_variant(record, classification);
        for (other_r, other_c) in current_variants.iter_mut() {
            if !overlaps(classification, other_c) {
                continue;
            }

            let other_score = score_variant(other_r, other_c);
            // arbitrarily we favour the earlier variant in a tie
            if score <= other_score {
                if !record.filter.contains(&OVERLAP_FILTER.to_string()) {
                    record.filter.push(OVERLAP_FILTER.to_string());
                }
            } else {
                if !other_r.filter.contains(&OVERLAP_FILTER.to_string()) {
                    other_r.filter.push(OVERLAP_FILTER.to_string());
                }
            }
        }

        current_variants.push((record, classification));
    }

    // Update classification for filtered records
    for (record, classification) in records.iter_mut() {
        if record.filter.contains(&OVERLAP_FILTER.to_string()) {
            classification.clone_from(&classifier.classify(record));
        }
    }
}

/// Wrap the results from processing support vcf
#[derive(Debug, Default)]
struct VCFResults {
    output_lines: Vec<(String, u32, String)>, // chrom, pos, vcf line
    processed_positions: HashMapSet<usize>,
    insertions: Vec<(VariantRecord, Classification)>,
    minors: Vec<VariantRecord>,
}

// Assumes rows are split for snps vs indels
fn process_main_vcf(
    vcf_file: &str,
    consensus: &mut HashMap<String, Vec<char>>,
    classifier: &Classifier,
    skip_indels: bool,
    caller_name: &Option<String>,
    verbose: bool,
) -> Result<VCFResults> {
    let vcf_reader = VCFReader::from_path(vcf_file)?;

    let mut results = VCFResults::default();

    let mut potential_output_records: Vec<(VariantRecord, Classification)> = Vec::new();

    for chrom in consensus.keys() {
        results
            .processed_positions
            .insert(chrom.clone(), HashSet::new());
    }

    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {count} records");
        }
        let mut record = record?;

        if !consensus.contains_key(&record.chrom) {
            return Err(format!("Chrom {} not found in consensus", record.chrom).into());
        }

        if skip_indels && record.is_indel() {
            continue;
        }

        // Add caller name to all records
        if let Some(caller) = caller_name {
            record
                .info
                .insert(CALLER.to_owned(), RecordValue::String(caller.to_owned()));
        }

        // Need 0-based position
        let pos: usize = record.pos_idx();

        // Update processed sites
        for i in 0..record.ref_bases.len() {
            results
                .processed_positions
                .insert_loc(&record.chrom, pos + i);
        }

        let classification = classifier.classify(&mut record);

        // simplify record to reduce ram
        simplify_record(&mut record);

        potential_output_records.push((record, classification));
    }

    mark_overlaps(&mut potential_output_records, classifier);

    // Any simple ref call which is overlapping with another variant can be ignored
    // The remaining variants are all informative
    potential_output_records.retain(|(r, c)| {
        if r.genotype().unwrap().is_hom_ref()
            && !c.has_minor_population
            && r.filter.contains(&OVERLAP_FILTER.to_string())
        {
            return false;
        }
        return true;
    });

    // We now apply all filtered records as we can write over the unfiltered variants afterwards
    for (r, c) in potential_output_records.iter() {
        if !c.is_filtered {
            continue;
        }
        apply_variant_simple(&r.chrom, c, consensus);
    }
    // Apply all non filtered records (which we've gauranteed don't overlap) apart from insertions
    for (r, c) in potential_output_records.iter() {
        if c.is_filtered {
            continue;
        }

        // Insertions only get applied at the very end, due to changing consensus length
        if c.change == Change::Ins || c.change == Change::ComplexIns {
            results.insertions.push((r.clone(), c.clone()));
            continue;
        }

        apply_variant_simple(&r.chrom, c, consensus);
    }

    results.minors = potential_output_records
        .iter()
        .filter(|(_, c)| c.has_minor_population)
        .map(|(r, _)| r.clone())
        .collect();

    results.output_lines = potential_output_records
        .into_iter()
        .map(|(r, _)| (r.chrom.clone(), r.pos, r.to_string()))
        .collect();

    return Ok(results);
}

/// Reads support vcf file to provide data on sites not in main vcf
///
/// Skips any indels, and will not apply any snps, only masking and nulls
fn process_support_vcf(
    vcf_file: &str,
    consensus: &mut HashMap<String, Vec<char>>,
    classifier: &Classifier,
    positions_already_processed: &HashMapSet<usize>,
    caller_name: &Option<String>,
    verbose: bool,
) -> Result<VCFResults> {
    let vcf_reader = VCFReader::from_path(vcf_file)?;

    let mut results = VCFResults::default();

    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {count} records of support vcf");
        }
        let mut record = record?;

        // Skips indels in support vcf
        if record.is_indel() {
            continue;
        }

        // Need to convert to 0-based
        let pos: usize = record.pos_idx();

        if let Some(processed_sites) = positions_already_processed.get(&record.chrom)
            && processed_sites.contains(&pos)
        {
            continue;
        }

        if !consensus.contains_key(&record.chrom) {
            return Err(format!("Chrom {} not found in consensus", record.chrom).into());
        }
        let chrom_seq: &mut Vec<char> = consensus
            .get_mut(&record.chrom)
            .expect("Chrom not found in consensus, but just checked that it is!");

        // By this point we know it is a single nucleotide ref or change
        if results
            .processed_positions
            .contains_loc(&record.chrom, &pos)
        {
            return Err(format!(
                "Site {}:{} already processed. Multiple snps on same site not allowed",
                record.chrom, record.chrom
            )
            .into());
        }
        results.processed_positions.insert_loc(&record.chrom, &pos);

        let mut classification = classifier.classify(&mut record);

        // mark snps if not calling them
        if !classifier.params.call_snps_in_support && !record.genotype().unwrap().is_hom_ref() {
            classification.change = Change::Null;
            classification.is_filtered = true;
            classification.new_bases = NULL.to_string();
            record.filter.push(SNP_IN_SUPPORT_VCF.to_string());
        }

        // Apply the changes, which we know are single base changes
        chrom_seq[pos] = classification
            .new_bases
            .chars()
            .next()
            .expect("New bases empty");

        // Add the caller name to the record
        if let Some(caller) = caller_name {
            record
                .info
                .insert(CALLER.to_owned(), RecordValue::String(caller.to_owned()));
        }

        // simplify record to reduce ram
        simplify_record(&mut record);

        if classification.has_minor_population {
            results.minors.push(record.clone());
        }

        results
            .output_lines
            .push((record.chrom.clone(), record.pos, record.to_string()));
    }

    return Ok(results);
}

fn write_creation_report(
    consensus: &HashMap<String, Vec<char>>,
    minors: &[VariantRecord],
    cluster_window: usize,
    output_file: &str,
) -> Result<()> {
    let mut letter_counts: HashMap<char, i32> = HashMap::new();
    let mut genome_length: i32 = 0;
    for (_, seq) in consensus.iter() {
        genome_length += seq.len() as i32;
        for base in seq.iter() {
            *letter_counts.entry(*base).or_insert(0) += 1;
        }
    }

    let null_calls = *letter_counts.get(&NULL).unwrap_or(&0);
    let fixed_coverage: f32 = 100.0 * (genome_length - null_calls) as f32 / genome_length as f32;

    // Count mixtures. For time being we mix hets and minor populations together
    let mixed_calls = minors.len() as i32;
    let mut mixed_snp_clusters = 0;

    let mut snps: Vec<(String, usize)> = minors
        .iter()
        .filter(|r| !r.is_indel())
        .map(|r| (r.chrom.clone(), r.pos_idx()))
        .collect();

    snps.sort_unstable();
    let mixed_snps = snps.len() as i32;

    if mixed_snps > 0 {
        let (mut last_chrom, mut last_site) = snps[0].clone();
        mixed_snp_clusters += 1; // first site is always a new cluster
        for (chrom, site) in snps.into_iter().skip(1) {
            if chrom != last_chrom || site - last_site > cluster_window {
                // new cluster
                mixed_snp_clusters += 1;
            }
            last_site = site;
            last_chrom = chrom;
        }
    }

    let quality_stats = SequencingQuality {
        genome_length,
        null_calls,
        deleted_calls: *letter_counts.get(&DELETED).unwrap_or(&0),
        fixed_coverage,
        mixed_snps,
        mixed_snp_clusters,
        mixed_calls,
    };
    let report = GenomeCreationReport {
        sequencing_quality: quality_stats,
    };

    let json_report = serde_json::to_string_pretty(&report)?;
    std::fs::write(output_file, json_report)?;

    Ok(())
}

fn simplify_record(record: &mut VariantRecord) {
    let mut new_info = IndexMap::new();
    for key in [CALLER, MIXED, REMOVED_MINORS] {
        if let Some(value) = record.info.get(key) {
            new_info.insert(key.to_owned(), value.clone());
        }
    }
    record.info = new_info;

    let mut format = IndexMap::new();
    format.insert(
        "GT".to_owned(),
        RecordValue::String(record.genotype().unwrap().to_string()),
    );

    format.insert(
        "DP".to_owned(),
        if let Some(dp) = record.depth() {
            RecordValue::Integer(*dp)
        } else {
            RecordValue::Integer(0)
        },
    );

    if let Some((adf, adr)) = record.strand_depths() {
        format.insert("ADF".to_owned(), RecordValue::IntegerArray(adf.clone()));
        format.insert("ADR".to_owned(), RecordValue::IntegerArray(adr.clone()));
    }
    if let Some(ad) = record.allele_depths() {
        format.insert("COV".to_owned(), RecordValue::IntegerArray(ad.clone()));
    }
    record.format = format;
}

/// Makes an empty record for a site which is missing from the vcf
///
/// expects 1-based position input
fn make_empty_record(chrom: &str, pos: u32, ref_bases: &str) -> VariantRecord {
    let mut record = VariantRecord::empty_record();
    record.chrom = chrom.to_owned();
    record.pos = pos;
    record.ref_bases = ref_bases.to_owned();
    record.set_genotype(Genotype::new());

    record
        .format
        .insert("DP".to_string(), RecordValue::Integer(0));
    record
        .format
        .insert("ADF".to_string(), RecordValue::IntegerArray(vec![0]));
    record
        .format
        .insert("ADR".to_string(), RecordValue::IntegerArray(vec![0]));
    record
        .format
        .insert("COV".to_string(), RecordValue::IntegerArray(vec![0]));

    return record;
}

fn write_vcf(
    record_lines: &[(String, u32, String)],
    output_file: &str,
    main_vcf: &str,
    support_vcf: Option<&str>,
) -> Result<()> {
    let mut header = VCFHeader::new_std_spec();
    header.add_misc_line(format!(
        "##rundial_consensus_version={}",
        env!("CARGO_PKG_VERSION")
    ));

    // Copy filters from main and support vcf
    for file in [Some(main_vcf), support_vcf].iter().flatten() {
        let reader = VCFReader::from_path(file)?;
        for line in reader.header().lines.iter() {
            if matches!(line, HeaderLine::Filter(_)) {
                header.add_header_line(line.clone());
            }
        }
    }

    header.add_filter_line(OVERLAP_FILTER.to_string(), OVERLAP_FILTER_DESC.to_string());
    header.add_filter_line(
        SNP_IN_SUPPORT_VCF.to_string(),
        SNP_IN_SUPPORT_VCF_DESC.to_string(),
    );
    header.add_filter_line(HET_GT.to_string(), HET_GT_DESC.to_string());

    header.add_info_line(
        CALLER.to_owned(),
        HeaderNumber::One,
        HeaderType::String,
        CALLER_DESC.to_owned(),
    );
    header.add_info_line(
        MIXED.to_owned(),
        HeaderNumber::Flag,
        HeaderType::Flag,
        MIXED_DESC.to_owned(),
    );
    header.add_info_line(
        REMOVED_MINORS.to_owned(),
        HeaderNumber::One,
        HeaderType::String,
        REMOVED_MINORS_DESC.to_owned(),
    );

    header.add_format_line(
        "GT".to_owned(),
        HeaderNumber::One,
        HeaderType::String,
        "Genotype".to_owned(),
    );
    header.add_format_line(
        "DP".to_owned(),
        HeaderNumber::One,
        HeaderType::Integer,
        "Basic depth".to_owned(),
    );
    header.add_format_line(
        "ADF".to_owned(),
        HeaderNumber::R,
        HeaderType::Integer,
        "Allelic depths for the ref and alt alleles on forward strand".to_owned(),
    );
    header.add_format_line(
        "ADR".to_owned(),
        HeaderNumber::R,
        HeaderType::Integer,
        "Allelic depths for the ref and alt alleles on reverse strand".to_owned(),
    );
    header.add_format_line(
        "COV".to_owned(),
        HeaderNumber::R,
        HeaderType::Integer,
        "Allelic depths for the ref and alt alleles in the order listed".to_owned(),
    );

    // Now write records
    let mut writer = VCFWriter::to_path(output_file, header)?;
    for (_, _, line) in record_lines.iter() {
        writer.write_line(line)?;
    }
    Ok(())
}

pub fn make_consensus(
    main_vcf: &str,
    support_vcf: Option<&str>,
    ref_fasta: &str,
    output_root: &str,
    params: &str,
    verbose: bool,
) -> Result<()> {
    let params: ConsensusParams = serde_yaml::from_reader(
        File::open(params).map_err(|e| format!("Failed to read params file. Error: {e}"))?,
    )?;

    let mut consensus = read_fasta(ref_fasta)?;
    let classifier = Classifier::new(&params, false);

    let mut main_results = process_main_vcf(
        main_vcf,
        &mut consensus,
        &classifier,
        params.skip_indels,
        &params.main_caller,
        verbose,
    )?;
    if verbose {
        let num_process_positions = main_results
            .processed_positions
            .values()
            .map(|s| s.len())
            .sum::<usize>();
        println!(
            "From main VCF: Processed {} sites. Found {} records to output",
            num_process_positions,
            main_results.output_lines.len(),
        );
    }

    let classifier = Classifier::new(&params, true);
    let support_results = match support_vcf {
        Some(file) => process_support_vcf(
            file,
            &mut consensus,
            &classifier,
            &main_results.processed_positions,
            &params.support_caller,
            verbose,
        )?,
        None => VCFResults::default(),
    };
    if verbose {
        let num_process_positions = support_results
            .processed_positions
            .values()
            .map(|s| s.len())
            .sum::<usize>();
        println!(
            "From support VCF: Processed {} sites. Found {} records to output",
            num_process_positions,
            support_results.output_lines.len(),
        );
    }

    main_results
        .output_lines
        .extend(support_results.output_lines);
    main_results
        .processed_positions
        .extend_all_chroms(support_results.processed_positions);
    main_results.minors.extend(support_results.minors);

    // Mask missing sites
    if params.mask_missing_sites {
        for (chrom, seq) in consensus.iter_mut() {
            let all_sites: HashSet<usize> = (0..seq.len()).collect();

            if let Some(processed_sites) = main_results.processed_positions.get(chrom) {
                for i in all_sites.difference(processed_sites) {
                    main_results.output_lines.push((
                        chrom.to_string(),
                        (i + 1) as u32,
                        make_empty_record(chrom, (i + 1) as u32, &seq[*i].to_string()).to_string(),
                    ));
                    seq[*i] = NULL;
                }
            } else {
                for i in all_sites {
                    main_results.output_lines.push((
                        chrom.to_string(),
                        (i + 1) as u32,
                        make_empty_record(chrom, (i + 1) as u32, &seq[i].to_string()).to_string(),
                    ));
                    seq[i] = NULL;
                }
            }
        }
    }

    write_creation_report(
        &consensus,
        &main_results.minors,
        12,
        &(output_root.to_owned() + ".report.json"),
    )?;

    save_fasta(&consensus, &(output_root.to_owned() + ".fasta"))?;

    // Apply insertion changes in reverse order!
    let mut variable_len_consensus = consensus;
    main_results.insertions.sort_by_key(|(_r, c)| c.pos);
    main_results.insertions.reverse();
    for (record, classification) in main_results.insertions.into_iter() {
        if classification.is_filtered {
            panic!("Filtered insertions should have been removed earlier");
        }

        apply_variant_simple(&record.chrom, &classification, &mut variable_len_consensus);
    }
    save_fasta(
        &variable_len_consensus,
        &(output_root.to_owned() + ".variable_length.fasta"),
    )?;

    // Output records to vcf
    main_results.output_lines.sort();

    write_vcf(
        &main_results.output_lines,
        &(output_root.to_owned() + ".vcf"),
        main_vcf,
        support_vcf,
    )?;

    Ok(())
}

#[cfg(test)]
mod tests;

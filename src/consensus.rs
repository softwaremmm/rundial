pub mod parameter_struct;
use core::panic;
use std::io::Write;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufReader,
};

use niffler;
use noodles::fasta::record::Sequence;

use crate::vcf::vcf_header::{HeaderLine, HeaderNumber, HeaderType};
use crate::vcf::{Genotype, RecordValue, VCFHeader, VCFReader, VCFWriter, VariantRecord};

use indexmap::IndexMap;
use noodles::fasta::{
    self as noodles_fasta,
    record::{Definition, Record},
};
use ordered_float::OrderedFloat;

pub use parameter_struct::{ConsensusParams, GenomeCreationReport, HetOption, SequencingQuality};

pub mod bed;
pub use bed::Bed;
mod classifier;
use classifier::{repeat_char, Change, Classification, Classifier};
mod contig_set;
use contig_set::ContigSet;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
type HashMapSet<T> = HashMap<String, HashSet<T>>;

const NULL: char = 'N';
const FILTERED: char = 'F';
const HET: char = 'Z';
const MASKED: char = 'M';
const DELETED: char = '-';

const OVERLAP_FILTER: &str = "OVERLAP_BETTER_VARIANT";
const OVERLAP_FILTER_DESC: &str = "This variant overlaps with a higher priority variant.";
const SNP_IN_SUPPORT_VCF: &str = "SNP_IN_SUPPORT_VCF";
const SNP_IN_SUPPORT_VCF_DESC: &str =
    "This variant is a SNP in the support VCF, which should only be used for information of nulls.";
const HET_IN_SUPPORT_VCF: &str = "HET_IN_SUPPORT_VCF";
const HET_IN_SUPPORT_VCF_DESC: &str =
    "This variant is a HET in the support VCF, which may be undesired.";
const CALLER: &str = "CALLER";
const CALLER_DESC: &str = "The variant caller that made the call.";

fn apply_variant(
    chrom: &str,
    classification: &Classification,
    chrom_seq: &mut HashMap<String, Vec<char>>,
    processed_sites: &HashMapSet<usize>,
) -> HashSet<usize> {
    let mut sites_set: HashSet<usize> = HashSet::new();
    let pos = classification.pos;

    // select the correct chromosome sequence
    let chrom_seq: &mut Vec<char> = chrom_seq
        .get_mut(chrom)
        .expect("Chrom not found in consensus");

    match classification.change {
        Change::Null => {
            assert!(classification.ref_bases.len() == classification.new_bases.len());
            for (i, base) in classification.new_bases.chars().enumerate() {
                if processed_sites.contains_loc(chrom, &(pos + i)) {
                    continue;
                }
                chrom_seq[pos + i] = base;
                sites_set.insert(pos + i);
            }
        }
        Change::HetMask => {
            assert!(classification.ref_bases.len() == classification.new_bases.len());
            for (i, base) in classification.new_bases.chars().enumerate() {
                chrom_seq[pos + i] = base;
                sites_set.insert(pos + i);
            }
        }
        Change::Ref => {
            for i in 0..classification.ref_bases.len() {
                if processed_sites.contains_loc(chrom, &(pos + i)) {
                    continue;
                }
                sites_set.insert(pos + i);
            }
        }
        Change::Snp => {
            assert!(classification.ref_bases.len() == 1 && classification.new_bases.len() == 1);
            sites_set.insert(pos);
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
                sites_set.insert(pos + i);
                chrom_seq[pos + i] = base;
            }
        }
        // sites_set does not make sense for indels which change consensus length
        Change::Ins | Change::ComplexIns => {
            // We assume that ref and alt have already been simplified
            let new_bases = classification.new_bases.chars();
            let ref_length = classification.ref_bases.len();
            chrom_seq.splice(pos..pos + ref_length, new_bases);
        }
    }
    return sites_set;
}

fn _check_indel_ref_matches_seq(
    chrom_seq: &[char],
    pos: usize, // 0-based pos
    expected_ref_bases: &str,
    skip_first: bool,
) -> bool {
    let null_sites = [NULL, FILTERED, HET, MASKED];
    let skip_num = if skip_first { 1 } else { 0 };
    for (i, base) in expected_ref_bases.chars().skip(skip_num).enumerate() {
        if null_sites.contains(&chrom_seq[pos + i]) {
            continue;
        }
        if chrom_seq[pos + i] != base {
            return false;
        }
    }
    true
}

fn score_variant(
    record: &VariantRecord,
    classification: &Classification,
) -> (i32, i32, OrderedFloat<f32>, i32) {
    let filter_score = if !classification.is_filtered { 2 } else { 0 };

    // If filtered, then snps > indels > null gt calls
    // otherwise: snps (and het snps) > indels > hom ref snp > hom ref indel > null gt calls
    let type_score = match (
        classification.is_filtered,
        record.genotype(),
        classification.has_indel_alleles,
    ) {
        (_, None, _) => 0,
        (_, Some(gt), _) if gt.is_null() => 0,
        (true, _, true) => 1,
        (true, _, false) => 2,
        (false, Some(gt), true) if gt.is_hom_ref() => 1,
        (false, Some(gt), false) if gt.is_hom_ref() => 2,
        (false, _, true) => 3,
        (false, _, false) => 4,
    };

    let qual_score = match record.qual {
        Some(f) => OrderedFloat(f),
        _ => OrderedFloat(0.0),
    };
    let depth_score = match record.depth() {
        Some(d) => *d,
        _ => 0,
    };
    return (filter_score, type_score, qual_score, depth_score);
}

fn overlaps(this: &Classification, other: &Classification) -> bool {
    fn range(c: &Classification) -> (f32, f32) {
        // Position changes on float scale where x.0 is the middle of a base x
        // so snp at base 2 would be 2.0-2.0
        // insertion just after base 2 would be 2.5-2.5
        // deletion at base 2-3 would be 2.0-3.0
        return match c.change {
            Change::Null
            | Change::HetMask
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

/// Adds "OverlapsBetterVariant" to Filter column for indels that overlap with other indels
/// that have a higher score
fn mark_overlaps(records: &mut [(VariantRecord, Classification)], classifier: &Classifier) {
    records.sort_by_key(|(r, c)| (r.chrom.clone(), c.pos));

    // let mut processed_variants: Vec<(VariantRecord, Classification)> = Vec::new();
    let mut current_variants: Vec<(&mut VariantRecord, &mut Classification)> = Vec::new();

    for (record, classification) in records.iter_mut() {
        // Discard variants which are too far back to possibly overlap with current variant
        let (mut in_range, _out_of_range): (Vec<_>, Vec<_>) =
            current_variants.into_iter().partition(|(other_r, _)| {
                if other_r.chrom != record.chrom {
                    return false;
                }
                if other_r.pos_idx() + other_r.ref_bases.len() < record.pos_idx() - 1 {
                    return false;
                }

                return true;
            });

        // Check for clashes and filter accordingly
        let score = score_variant(record, classification);
        for (other_r, other_c) in in_range.iter_mut() {
            if !overlaps(classification, other_c) {
                continue;
            }

            let other_score = score_variant(other_r, other_c);
            if score <= other_score && !record.filter.contains(&OVERLAP_FILTER.to_string()) {
                record.filter.push(OVERLAP_FILTER.to_string());
            }
            if score >= other_score && !other_r.filter.contains(&OVERLAP_FILTER.to_string()) {
                other_r.filter.push(OVERLAP_FILTER.to_string());
            }
        }

        current_variants = in_range;
        current_variants.push((record, classification));
    }

    // Update classification for filtered records
    for (record, classification) in records.iter_mut() {
        if record.filter.contains(&OVERLAP_FILTER.to_string()) {
            classification.clone_from(&classifier.classify(record));
        }
    }
}

// Assumes rows are split for snps vs indels
#[allow(clippy::type_complexity)]
fn process_main_vcf(
    vcf_file: &str,
    consensus: &mut HashMap<String, Vec<char>>,
    classifier: &Classifier,
    skip_indels: bool,
    verbose: bool,
) -> Result<(
    Vec<VariantRecord>,
    Vec<(VariantRecord, Classification)>,
    HashMapSet<usize>,
    HashMapSet<usize>,
)> {
    let vcf_reader = VCFReader::from_path(vcf_file)?;

    let mut output_records: Vec<VariantRecord> = Vec::new();
    let mut insertions: Vec<(VariantRecord, Classification)> = Vec::new();
    let mut potential_output_records: Vec<(VariantRecord, Classification)> = Vec::new();
    // All sites seen in the vcf. Used to mark sites that are not in the vcf
    let mut seen_sites: HashMapSet<usize> = HashMap::new();
    let mut processed_positions: HashMapSet<usize> = HashMap::new();
    let mut het_sites: HashMapSet<usize> = HashMap::new();
    for chrom in consensus.keys() {
        processed_positions.insert(chrom.clone(), HashSet::new());
        seen_sites.insert(chrom.clone(), HashSet::new());
        het_sites.insert(chrom.clone(), HashSet::new());
    }

    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {} records", count);
        }
        let mut record = record?;

        if !consensus.contains_key(&record.chrom) {
            return Err(format!("Chrom {} not found in consensus", record.chrom).into());
        }
        let chrom_seq: &mut Vec<char> = consensus
            .get_mut(&record.chrom)
            .expect("Chrom not found in consensus");

        // Need to convert to 0-based
        let pos: usize = record.pos_idx();
        if classifier.is_masked(&record) {
            for i in 0..record.ref_bases.len() {
                chrom_seq[pos + i] = MASKED;
            }
            continue;
        }

        if skip_indels && record.is_indel() {
            continue;
        }

        // Update seen sites
        for i in 0..record.ref_bases.len() {
            seen_sites.insert_loc(&record.chrom, pos + i);
        }

        let classification = classifier.classify(&mut record);
        if classification.is_simple_ref() {
            processed_positions.insert_loc(&record.chrom, &pos);
            continue;
        }

        potential_output_records.push((record, classification));
    }

    mark_overlaps(&mut potential_output_records, classifier);

    // TODO: remove all the cloning in the following
    // Apply unfiltered records
    for (r, c) in potential_output_records.iter() {
        if c.is_filtered {
            continue;
        }
        if c.change == Change::Ins || c.change == Change::ComplexIns {
            insertions.push((r.clone(), c.clone()));
            continue;
        }

        let affected_sites = apply_variant(&r.chrom, c, consensus, &processed_positions);
        processed_positions.extend_chrom(&r.chrom, &affected_sites);

        if c.is_het {
            het_sites.extend_chrom(&r.chrom, &affected_sites);
        }

        output_records.push(r.clone());
    }

    // Apply filtered non-indel changes records
    // We don't look for hets in this section due to filters
    for (r, c) in potential_output_records.iter() {
        if !c.is_filtered || r.is_indel() {
            continue;
        }

        let affected_sites = apply_variant(&r.chrom, c, consensus, &processed_positions);
        // Don't output filtered records that don't change the consensus
        if affected_sites.is_empty() {
            continue;
        }

        processed_positions.extend_chrom(&r.chrom, affected_sites);
        output_records.push(r.clone());
    }

    // Apply filtered indels
    for (r, c) in potential_output_records.iter() {
        if !c.is_filtered || !r.is_indel() {
            continue;
        }

        let affected_sites = apply_variant(&r.chrom, c, consensus, &processed_positions);

        // Don't output filtered records that don't change the consensus
        if affected_sites.is_empty() {
            // Unless the only filters are MIN_FRS and overlaps
            let allowed_filters: HashSet<String> = HashSet::from([
                "MIN_FRS".to_owned(),
                "OverlapWithSnp".to_owned(),
                "OverlapWithIndel".to_owned(),
                OVERLAP_FILTER.to_owned(),
            ]);
            if r.filter.iter().all(|f| allowed_filters.contains(f)) {
                output_records.push(r.clone());
            }
            continue;
        }

        processed_positions.extend_chrom(&r.chrom, affected_sites);
        output_records.push(r.clone());
    }

    // Now add simple_ref sites to processed sites
    for (chrom, sites) in seen_sites.iter() {
        processed_positions.extend_chrom(chrom, sites);
    }

    return Ok((output_records, insertions, processed_positions, het_sites));
}

/// Reads support vcf file to provide data on sites not in main vcf
///
/// Skips any indels, and will not apply any snps, only masking and nulls
fn process_support_vcf(
    vcf_file: &str,
    consensus: &mut HashMap<String, Vec<char>>,
    classifier: &Classifier,
    positions_already_processed: &HashMapSet<usize>,
    verbose: bool,
) -> Result<(Vec<VariantRecord>, HashMapSet<usize>, HashMapSet<usize>)> {
    let vcf_reader = VCFReader::from_path(vcf_file)?;

    let mut output_records: Vec<VariantRecord> = Vec::new();
    let mut processed_positions: HashMapSet<usize> = HashMap::new();
    let mut het_sites: HashMapSet<usize> = HashMap::new();

    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {} records of support vcf", count);
        }
        let mut record = record?;
        // Need to convert to 0-based
        let pos: usize = record.pos_idx();

        if let Some(processed_sites) = positions_already_processed.get(&record.chrom) {
            if processed_sites.contains(&pos) {
                continue;
            }
        }

        if !consensus.contains_key(&record.chrom) {
            return Err(format!("Chrom {} not found in consensus", record.chrom).into());
        }
        let chrom_seq: &mut Vec<char> = consensus
            .get_mut(&record.chrom)
            .expect("Chrom not found in consensus, but just checked that it is!");

        if classifier.is_masked(&record) {
            for i in 0..record.ref_bases.len() {
                chrom_seq[pos + i] = MASKED;
            }
        }

        // Skips indels in support vcf
        if record.is_indel() {
            continue;
        }

        // By this point we know it is a single nucleotide ref or change
        if processed_positions.contains_loc(&record.chrom, &pos) {
            return Err(format!(
                "Site {}:{} already processed. Multiple snps on same site not allowed",
                record.chrom, record.chrom
            )
            .into());
        }
        processed_positions.insert_loc(&record.chrom, &pos);

        let mut classification = classifier.classify(&mut record);

        // See if it is a simple ref that can be skipped
        if classification.is_simple_ref() {
            continue;
        }

        // marks snps and hets
        if classification.change == Change::Snp {
            classification.change = Change::Null;
            classification.is_filtered = true;
            classification.new_bases = repeat_char(FILTERED, record.ref_bases.len());
            record.filter.push(SNP_IN_SUPPORT_VCF.to_string());
        } else if classification.is_het {
            record
                .info
                .insert(HET_IN_SUPPORT_VCF.to_string(), RecordValue::Flag);
            het_sites.insert_loc(&record.chrom, &pos);
        }

        // Apply the changes, which we know are single base changes
        chrom_seq[pos] = classification
            .new_bases
            .chars()
            .next()
            .expect("New bases empty");
        output_records.push(record);
    }

    return Ok((output_records, processed_positions, het_sites));
}

fn read_fasta(fasta_file: &str) -> Result<HashMap<String, Vec<char>>> {
    let mut consensus: HashMap<String, Vec<char>> = HashMap::new();

    let (reader, _format) = niffler::from_path(fasta_file).map_err(|e| {
        format!(
            "Failed to open fasta input file {}. Error: {}",
            fasta_file, e
        )
    })?;
    let buf_reader = BufReader::new(reader);
    let mut fasta_reader = noodles_fasta::Reader::new(buf_reader);
    for record in fasta_reader.records() {
        let record = record?;
        let chrom = String::from_utf8_lossy(record.definition().name()).to_string();
        let seq: Vec<char> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|c| *c as char)
            .collect();
        consensus.insert(chrom, seq);
    }

    return Ok(consensus);
}

fn potentially_gzipped_writer(file: &str) -> Result<Box<dyn Write>> {
    let (nif_format, level) = if file.ends_with(".gz") {
        (
            niffler::compression::Format::Gzip,
            niffler::compression::Level::One,
        )
    } else {
        (
            niffler::compression::Format::No,
            niffler::compression::Level::Zero,
        )
    };

    let niffler_writer = niffler::to_path(file, nif_format, level)
        .map_err(|e| format!("Failed to open fasta output file {}. Error: {}", file, e))?;

    return Ok(niffler_writer);
}

fn save_fasta(consensus: &HashMap<String, Vec<char>>, output_file: &str) -> Result<()> {
    let mut writer = noodles_fasta::Writer::new(potentially_gzipped_writer(output_file)?);

    for (chrom, seq) in consensus.iter() {
        let definition = Definition::new(chrom.clone(), None);
        let sequence = Sequence::from(seq.iter().collect::<String>().as_bytes().to_vec());
        let record = Record::new(definition, sequence);
        writer.write_record(&record)?;
    }
    return Ok(());
}

fn clean_fasta_characters(consensus: &mut HashMap<String, Vec<char>>) {
    for (_, seq) in consensus.iter_mut() {
        for base in seq.iter_mut() {
            if [FILTERED, HET, MASKED].contains(base) {
                *base = NULL;
            }
        }
    }
}

fn write_creation_report(
    consensus: &HashMap<String, Vec<char>>,
    het_count: Option<i32>,
    output_file: &str,
) -> Result<()> {
    let mut letter_counts: HashMap<char, i32> = HashMap::new();
    let mut total_length: i32 = 0;
    for (_, seq) in consensus.iter() {
        total_length += seq.len() as i32;
        for base in seq.iter() {
            *letter_counts.entry(*base).or_insert(0) += 1;
        }
    }

    let mut all_null_counts = 0;
    for base in [NULL, FILTERED, HET, MASKED].iter() {
        all_null_counts += letter_counts.get(base).unwrap_or(&0);
    }

    let fixed_cov: f32 = 100.0 * (total_length - all_null_counts) as f32 / total_length as f32;

    let mixed_count = match het_count {
        Some(c) => c,
        None => *letter_counts.get(&HET).unwrap_or(&0),
    };

    let quality_stats = SequencingQuality {
        genome_length: total_length,
        null_calls: all_null_counts,
        mixed_calls: mixed_count,
        fixed_coverage: fixed_cov,
        null_genotype_calls: *letter_counts.get(&NULL).unwrap_or(&0),
        filtered_calls: *letter_counts.get(&FILTERED).unwrap_or(&0),
        masked_calls: *letter_counts.get(&MASKED).unwrap_or(&0),
        deleted_calls: *letter_counts.get(&DELETED).unwrap_or(&0),
    };
    let report = GenomeCreationReport {
        sequencing_quality: quality_stats,
    };

    let json_report = serde_json::to_string_pretty(&report)?;
    std::fs::write(output_file, json_report)?;

    Ok(())
}

fn write_vcf(
    records: &[VariantRecord],
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

    header.add_info_line(
        CALLER.to_owned(),
        HeaderNumber::One,
        HeaderType::String,
        CALLER_DESC.to_owned(),
    );
    header.add_info_line(
        HET_IN_SUPPORT_VCF.to_owned(),
        HeaderNumber::Flag,
        HeaderType::Flag,
        HET_IN_SUPPORT_VCF_DESC.to_owned(),
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
    for record in records.iter() {
        let mut output_record = record.clone();
        output_record.info = IndexMap::new();
        if let Some(caller) = record.info.get(CALLER) {
            output_record.info.insert(CALLER.to_owned(), caller.clone());
        }

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
        output_record.format = format;
        writer.write_record(&output_record)?;
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
) -> Result<HashMap<String, Vec<char>>> {
    let params: ConsensusParams = serde_yaml::from_reader(
        File::open(params).map_err(|e| format!("Failed to read params file. Error: {}", e))?,
    )?;

    let mut consensus = read_fasta(ref_fasta)?;
    let classifier = Classifier::new(&params);

    let skip_indels: bool = params.skip_indels;

    let (mut output_records, mut insertions, mut processed_positions, mut het_sites) =
        process_main_vcf(main_vcf, &mut consensus, &classifier, skip_indels, verbose)?;
    if let Some(caller) = params.main_caller {
        for r in output_records.iter_mut() {
            r.info
                .insert(CALLER.to_owned(), RecordValue::String(caller.to_owned()));
        }
    }
    if verbose {
        let num_process_positions = processed_positions.values().map(|s| s.len()).sum::<usize>();
        println!(
            "From main VCF: Processed {} sites. Found {} records to output, and {} het sites",
            num_process_positions,
            output_records.len(),
            het_sites.values().map(|s| s.len()).sum::<usize>()
        );
    }

    let (mut support_records, support_processed_positions, support_het_sites) = match support_vcf {
        Some(file) => process_support_vcf(
            file,
            &mut consensus,
            &classifier,
            &processed_positions,
            verbose,
        )?,
        None => (Vec::new(), HashMap::new(), HashMap::new()),
    };
    if let Some(caller) = params.support_caller {
        for r in support_records.iter_mut() {
            r.info
                .insert(CALLER.to_owned(), RecordValue::String(caller.to_owned()));
        }
    }
    if verbose {
        let num_process_positions = support_processed_positions
            .values()
            .map(|s| s.len())
            .sum::<usize>();
        println!(
            "From support VCF: Processed {} sites. Found {} records to output, and {} het sites",
            num_process_positions,
            support_records.len(),
            support_het_sites.values().map(|s| s.len()).sum::<usize>(),
        );
    }

    output_records.extend(support_records);

    for (chrom, sites) in support_processed_positions.into_iter() {
        processed_positions.extend_chrom(&chrom, sites);
    }
    for (chrom, sites) in support_het_sites.into_iter() {
        het_sites.extend_chrom(&chrom, sites);
    }

    fn make_empty_record(chrom: &str, pos: &usize, ref_bases: &str) -> VariantRecord {
        let mut record = VariantRecord::empty_record();
        record.chrom = chrom.to_owned();
        record.pos = *pos as u32 + 1;
        record.ref_bases = ref_bases.to_owned();
        record.set_genotype(Genotype::new());
        record
            .info
            .insert("DP".to_string(), RecordValue::Integer(0));
        record
            .info
            .insert("ADF".to_string(), RecordValue::IntegerArray(vec![0]));
        record
            .info
            .insert("ADR".to_string(), RecordValue::IntegerArray(vec![0]));
        record
            .format
            .insert("AD".to_string(), RecordValue::IntegerArray(vec![0]));
        record.update_depths();
        return record;
    }

    // Mask missing sites
    if params.mask_missing_sites {
        for (chrom, seq) in consensus.iter_mut() {
            let all_sites: HashSet<usize> = (0..seq.len()).collect();

            if let Some(processed_sites) = processed_positions.get(chrom) {
                for i in all_sites.difference(processed_sites) {
                    output_records.push(make_empty_record(chrom, i, &seq[*i].to_string()));
                    seq[*i] = NULL;
                }
            } else {
                for i in all_sites {
                    output_records.push(make_empty_record(chrom, &i, &seq[i].to_string()));
                    seq[i] = NULL;
                }
            }
        }
    }

    save_fasta(&consensus, &(output_root.to_owned() + ".full.fasta"))?;

    write_creation_report(
        &consensus,
        Some(het_sites.values().map(|s| s.len()).sum::<usize>() as i32),
        &(output_root.to_owned() + ".report.json"),
    )?;

    let mut clean_consensus = consensus.clone();
    clean_fasta_characters(&mut clean_consensus);
    save_fasta(&clean_consensus, &(output_root.to_owned() + ".fasta"))?;

    // Apply insertion changes in reverse order!
    let mut variable_len_consensus = clean_consensus;
    insertions.sort_by_key(|(_r, c)| c.pos);
    insertions.reverse();
    for (record, classification) in insertions.into_iter() {
        if classification.is_filtered {
            panic!("Filtered insertions should have been removed earlier");
        }

        apply_variant(
            &record.chrom,
            &classification,
            &mut variable_len_consensus,
            &processed_positions,
        );
        output_records.push(record);
    }
    clean_fasta_characters(&mut variable_len_consensus);
    save_fasta(
        &variable_len_consensus,
        &(output_root.to_owned() + ".variable_length.fasta"),
    )?;

    output_records.sort_by_key(|r| (r.chrom.clone(), r.pos, r.is_indel()));
    write_vcf(
        &output_records,
        &(output_root.to_owned() + ".vcf"),
        main_vcf,
        support_vcf,
    )?;

    Ok(consensus)
}

#[cfg(test)]
mod tests;

pub mod parameter_struct;
use core::panic;
use std::cell::RefCell;
use std::io::BufWriter;
use std::rc::Rc;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufReader,
};

use crate::vcf::vcf_header::{HeaderLine, HeaderNumber, HeaderType};
use crate::vcf::{RecordValue, VCFHeader, VCFReader, VCFWriter, VariantRecord};
use bio::io::fasta;
use indexmap::IndexMap;
use ordered_float::OrderedFloat;

pub use parameter_struct::{ConsensusParams, GenomeCreationReport, HetOption, SequencingQuality};

pub mod bed;
pub use bed::Bed;
mod classifier;
use classifier::{repeat_char, Change, Classification, Classifier};

const NULL: char = 'N';
const FILTERED: char = 'F';
const HET: char = 'Z';
const MASKED: char = 'M';
const DELETED: char = '-';

type Loc = (String, usize);

fn apply_variant(
    classification: &Classification,
    chrom_seq: &mut Vec<char>,
    processed_sites: &HashSet<usize>,
) -> HashSet<usize> {
    let mut sites_set: HashSet<usize> = HashSet::new();
    let pos = classification.pos;

    match classification.change {
        Change::Null => {
            for i in 0..classification.ref_bases.len() {
                if processed_sites.contains(&(pos + i)) {
                    continue;
                }
                chrom_seq[pos + i] = classification
                    .new_bases
                    .chars()
                    .nth(i)
                    .expect("New bases should have same length as ref bases for Nulls");
                sites_set.insert(pos + i);
            }
        }
        Change::Ref => {
            for i in 0..classification.ref_bases.len() {
                if processed_sites.contains(&(pos + i)) {
                    continue;
                }
                sites_set.insert(pos + i);
            }
        }
        Change::Snp => {
            sites_set.insert(pos);
            chrom_seq[pos] = classification
                .new_bases
                .chars()
                .next()
                .expect("Snp new bases empty");
        }
        Change::Del => {
            let replacement_str: String = classification.new_bases.clone()
                + &repeat_char(
                    DELETED,
                    classification.ref_bases.len() - classification.new_bases.len(),
                );
            for (i, base) in replacement_str.chars().enumerate() {
                if i == 0 {
                    // First base in a deletion row should be same in ref and alt
                    // so skip it if already processed by another row
                    if processed_sites.contains(&pos) {
                        continue;
                    }
                }

                sites_set.insert(pos + i);
                chrom_seq[pos + i] = base;
            }
        }
        // sites_set does not make sense for indels which change consensus length
        Change::Ins => {
            let new_bases = classification.new_bases.chars().skip(1);
            chrom_seq.splice(pos + 1..pos + 1, new_bases);
        }
        Change::ComplexIndel => {
            let new_bases = classification.new_bases.chars();
            chrom_seq.splice(pos..pos + classification.ref_bases.len(), new_bases);
        }
    }
    return sites_set;
}

fn check_indel_ref_matches_seq(
    chrom_seq: &Vec<char>,
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

/// Scores an indel based on the record. Higher is better
fn score_indel(record: &VariantRecord) -> (i32, OrderedFloat<f32>, i32) {
    let filter_score = if record.filter.is_empty() { 1 } else { 0 };
    let qual_score = match record.qual {
        Some(f) => OrderedFloat(f),
        _ => OrderedFloat(0.0),
    };
    let depth_score = match record.depth() {
        Some(d) => *d,
        _ => 0,
    };
    return (filter_score, qual_score, depth_score);
}

/// Adds "OverlapWithIndel" to Filter column for indels that overlap with other indels
/// that have a higher score
fn filter_overlapping_indels(indels: &mut Vec<VariantRecord>) -> () {
    let indels: Vec<Rc<RefCell<&mut VariantRecord>>> = indels
        .into_iter()
        .map(|r| Rc::new(RefCell::new(r)))
        .collect();

    // Then check for overlaps with other indels
    let mut indel_positions: HashMap<Loc, Vec<Rc<RefCell<&mut VariantRecord>>>> = HashMap::new();
    for record_ptr in indels.iter() {
        let record: &VariantRecord = &record_ptr.borrow();
        let pos: usize = record.pos as usize - 1;
        let chrom: String = record.chrom.clone();
        for i in pos..pos + record.ref_bases.len() {
            let loc: Loc = (chrom.clone(), i);
            indel_positions
                .entry(loc)
                .or_insert(Vec::new())
                .push(Rc::clone(record_ptr));
        }
    }

    // For each site with multiple indels, mark all but best as filtered
    for (_, records) in indel_positions.into_iter() {
        if records.len() <= 1 {
            continue;
        }

        let best_record: &Rc<RefCell<&mut VariantRecord>> = records
            .iter()
            .max_by_key(|r| score_indel(&r.borrow()))
            .unwrap();

        for record in records.iter() {
            if Rc::ptr_eq(record, best_record) {
                continue;
            }
            record
                .borrow_mut()
                .filter
                .push("OverlapWithIndel".to_string());
        }
    }
}

// Assumes rows are split for snps vs indels
fn process_main_vcf(
    vcf_file: &str,
    consensus: &mut HashMap<String, Vec<char>>,
    classifier: &Classifier,
    skip_indels: bool,
    verbose: bool,
) -> Result<
    (
        Vec<VariantRecord>,
        Vec<VariantRecord>,
        HashMap<String, HashSet<usize>>,
        HashMap<String, HashSet<usize>>,
    ),
    Box<dyn std::error::Error>,
> {
    let vcf_reader = VCFReader::new(BufReader::new(
        File::open(vcf_file).map_err(|e| format!("Failed to read input vcf file. Error: {}", e))?,
    ))?;

    let mut indels: Vec<VariantRecord> = Vec::new();
    let mut output_records: Vec<VariantRecord> = Vec::new();
    let mut processed_positions: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut het_sites: HashMap<String, HashSet<usize>> = HashMap::new();

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
            .expect("Chrom not found in consensus, but just checked that it is!");

        // Need to convert to 0-based
        let pos: usize = record.pos as usize - 1;
        if classifier.is_masked(&record) {
            for i in 0..record.ref_bases.len() {
                chrom_seq[pos + i] = MASKED;
            }
            continue;
        }

        if record.is_indel() {
            // Process indels later
            indels.push(record);
            continue;
        }

        // Check that the position has not been processed already
        let positions = processed_positions
            .entry(record.chrom.clone())
            .or_insert(HashSet::new());
        if positions.contains(&pos) {
            return Err(format!(
                "Site {}:{} already processed. Multiple snps on same site not allowed",
                record.chrom, record.chrom
            )
            .into());
        }
        positions.insert(pos);

        let classification = classifier.classify(&mut record);

        // See if it is a simple ref that can be skipped
        if classification.change == Change::Ref
            && !classification.has_minor_population
            && !classification.is_het
        {
            continue;
        }

        if classification.is_het {
            het_sites
                .entry(record.chrom.clone())
                .or_insert(HashSet::new())
                .insert(pos);
        }

        // Apply the changes, which we know are single base changes
        chrom_seq[pos] = classification
            .new_bases
            .chars()
            .next()
            .expect("New bases empty");

        output_records.push(record);
    }

    if skip_indels {
        return Ok((output_records, Vec::new(), processed_positions, het_sites));
    }

    // Now process the indels
    // First check for overlaps withs snps

    for mut record in indels.iter_mut() {
        // TODO: Could consider normalising indels by removing alts not in GT or at minor indel threshold

        let classification = classifier.classify(&mut record);

        let chrom_seq: &mut Vec<char> = consensus
            .get_mut(&record.chrom)
            .expect("Chrom not found in consensus, but just checked that it is!");
        let pos: usize = record.pos as usize - 1;

        let skip_first_base = matches!(classification.change, Change::Del | Change::Ins);
        if !check_indel_ref_matches_seq(chrom_seq, pos, &classification.ref_bases, skip_first_base)
        {
            record.filter.push("OverlapWithSnp".to_string());
        }
    }

    // Check for overlaps with other indels
    filter_overlapping_indels(&mut indels);

    let mut insertions: Vec<VariantRecord> = Vec::new();
    // Now apply the indels
    for mut indel in indels.into_iter() {
        let pos: usize = indel.pos as usize - 1;
        let classification = classifier.classify(&mut indel);
        let chrom_seq: &mut Vec<char> = consensus
            .get_mut(&indel.chrom)
            .expect("Chrom not found in consensus, but checked that it is!");
        let processed_sites = processed_positions
            .entry(indel.chrom.clone())
            .or_insert(HashSet::new());
        let chrom_het_sites = het_sites
            .entry(indel.chrom.clone())
            .or_insert(HashSet::new());

        if matches!(classification.change, Change::Ins | Change::ComplexIndel) {
            // These change consensus length, so dealth with separately
            insertions.push(indel);
            continue;
        }

        let sites_set = apply_variant(&classification, chrom_seq, &processed_sites);
        processed_sites.extend(sites_set.iter());

        if classification.is_het {
            chrom_het_sites.extend(sites_set.iter());
        } else {
            for site in sites_set.iter() {
                chrom_het_sites.remove(&site);
            }
        }

        match classification.change {
            Change::Null => {
                if !sites_set.is_empty() {
                    output_records.push(indel);
                } else {
                    // If only filters were MIN_FRS ands overlaps then output anyway
                    let allowed_filters: HashSet<String> = HashSet::from([
                        "MIN_FRS".to_owned(),
                        "OverlapWithSnp".to_owned(),
                        "OverlapWithIndel".to_owned(),
                    ]);
                    if indel.filter.iter().any(|f| !allowed_filters.contains(f)) {
                        output_records.push(indel);
                    }
                }
            }
            Change::Ref => {
                if !sites_set.is_empty() {
                    output_records.push(indel);
                }
            }
            Change::Snp | Change::Del => {
                output_records.push(indel);
            }
            Change::Ins | Change::ComplexIndel => {
                panic!("Should not have insertions here");
            }
        }
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
    positions_already_processed: &HashMap<String, HashSet<usize>>,
    verbose: bool,
) -> Result<(Vec<VariantRecord>, HashMap<String, HashSet<usize>>), Box<dyn std::error::Error>> {
    let vcf_reader = VCFReader::new(BufReader::new(
        File::open(vcf_file).map_err(|e| format!("Failed to read input vcf file. Error: {}", e))?,
    ))?;

    let mut output_records: Vec<VariantRecord> = Vec::new();
    let mut processed_positions: HashMap<String, HashSet<usize>> = HashMap::new();

    for (count, record) in vcf_reader.enumerate() {
        if verbose && count % 100000 == 0 && count != 0 {
            println!("Processed {} records of support vcf", count);
        }
        let mut record = record?;
        // Need to convert to 0-based
        let pos: usize = record.pos as usize - 1;

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

        if record.is_indel() {
            continue;
        }

        // By this point we know it is a single nucleotide ref or change
        let positions = processed_positions
            .entry(record.chrom.clone())
            .or_insert(HashSet::new());
        if positions.contains(&pos) {
            return Err(format!(
                "Site {}:{} already processed. Multiple snps on same site not allowed",
                record.chrom, record.chrom
            )
            .into());
        }
        positions.insert(pos);

        let mut classification = classifier.classify(&mut record);

        // See if it is a simple ref that can be skipped
        if classification.change == Change::Ref
            && !classification.has_minor_population
            && !classification.is_het
        {
            continue;
        }

        todo!("Need to check for hets here");

        if classification.change == Change::Snp {
            classification.change = Change::Null;
            classification.is_filtered = true;
            classification.new_bases = repeat_char(FILTERED, record.ref_bases.len());
            record.filter.push("SnpInSupportVCF".to_string());
        }

        output_records.push(record);

        // Apply the changes, which we know are single base changes
        chrom_seq[pos] = classification
            .new_bases
            .chars()
            .next()
            .expect("New bases empty");
    }

    return Ok((output_records, processed_positions));
}

fn read_fasta(fasta_file: &str) -> Result<HashMap<String, Vec<char>>, Box<dyn std::error::Error>> {
    let reader = fasta::Reader::from_file(fasta_file)?;
    let mut consensus: HashMap<String, Vec<char>> = HashMap::new();
    for record in reader.records() {
        let record = record.expect("Error during fasta record reading");
        let chrom = record.id().to_string();
        let seq: Vec<char> = record.seq().iter().map(|c| *c as char).collect();
        consensus.insert(chrom, seq);
    }
    return Ok(consensus);
}

fn save_fasta(
    consensus: &HashMap<String, Vec<char>>,
    output_file: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = fasta::Writer::to_file(output_file)?;
    for (chrom, seq) in consensus.iter() {
        writer.write(chrom, None, &seq.iter().collect::<String>().as_bytes())?;
    }
    return Ok(());
}

fn clean_fasta_characters(consensus: &mut HashMap<String, Vec<char>>) {
    for (_, seq) in consensus.iter_mut() {
        for base in seq.iter_mut() {
            if [FILTERED, HET, MASKED, DELETED].contains(base) {
                *base = NULL;
            }
        }
    }
}

fn write_creation_report(
    consensus: &HashMap<String, Vec<char>>,
    het_count: Option<i32>,
    output_file: &str,
) -> Result<(), Box<dyn std::error::Error>> {
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
        None => letter_counts.get(&HET).unwrap_or(&0).clone(),
    };

    let quality_stats = SequencingQuality {
        genome_length: total_length,
        null_calls: all_null_counts,
        mixed_calls: mixed_count,
        fixed_coverage: fixed_cov,
        null_genotype_calls: letter_counts.get(&NULL).unwrap_or(&0).clone(),
        filtered_calls: letter_counts.get(&FILTERED).unwrap_or(&0).clone(),
        masked_calls: letter_counts.get(&MASKED).unwrap_or(&0).clone(),
        deleted_calls: letter_counts.get(&DELETED).unwrap_or(&0).clone(),
    };
    let report = GenomeCreationReport {
        sequencing_quality: quality_stats,
    };

    let json_report = serde_json::to_string_pretty(&report)?;
    std::fs::write(output_file, json_report)?;

    Ok(())
}

fn write_vcf(
    records: &Vec<VariantRecord>,
    output_file: &str,
    main_vcf: &str,
    support_vcf: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut header = VCFHeader::new_std_spec();
    header.add_misc_line(format!(
        "##rundial consensus version{}",
        env!("CARGO_PKG_VERSION")
    ));

    // Copy filters from main and support vcf
    for file in [Some(main_vcf), support_vcf].iter() {
        if let Some(file) = file {
            let reader = VCFReader::new(BufReader::new(
                File::open(file)
                    .map_err(|e| format!("Failed to read main vcf file. Error: {}", e))?,
            ))?;
            for line in reader.header().lines.iter() {
                if matches!(line, HeaderLine::Filter(_)) {
                    header.add_header_line(line.clone());
                }
            }
        }
    }

    header.add_filter_line(
        "OverlapWithIndel".to_owned(),
        "This variant overlaps with an indel of higher quality.".to_owned(),
    );
    header.add_filter_line(
        "OverlapWithSnp".to_owned(),
        "This variant overlaps with an snp. Snps are given higher priority.".to_owned(),
    );

    header.add_info_line(
        "CALLER".to_owned(),
        HeaderNumber::One,
        HeaderType::String,
        "The variant caller that made the call.".to_owned(),
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
    let mut writer = VCFWriter::new(
        BufWriter::new(
            File::open(output_file)
                .map_err(|e| format!("Failed to open output vcf file. Error: {}", e))?,
        ),
        header,
    )?;
    for record in records.iter() {
        let mut output_record = record.clone();
        output_record.info = IndexMap::from([(
            "Caller".to_owned(),
            RecordValue::String("consensus".to_owned()),
        )]);

        let mut format = IndexMap::new();
        format.insert("GT".to_owned(), record.format.get("GT").unwrap().clone());
        if let Some(dp) = record.depth() {
            format.insert("DP".to_owned(), RecordValue::Integer(*dp));
        }
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
) -> Result<HashMap<String, Vec<char>>, Box<dyn std::error::Error>> {
    let params: ConsensusParams = serde_yaml::from_reader(
        File::open(params).map_err(|e| format!("Failed to read params file. Error: {}", e))?,
    )?;

    let mut consensus = read_fasta(ref_fasta)?;
    let classifier = Classifier::new(&params);

    let skip_indels: bool = params.skip_indels;

    let (mut main_records, mut insertions, mut processed_positions, het_sites) =
        process_main_vcf(main_vcf, &mut consensus, &classifier, skip_indels, verbose)?;
    if verbose {
        let num_process_positions = processed_positions.values().map(|s| s.len()).sum::<usize>();
        println!(
            "From main VCF: Processed {} sites. Found {} records to output, and {} het sites",
            num_process_positions,
            main_records.len(),
            het_sites.len()
        );
    }
    let (support_records, support_processed_positions) = match support_vcf {
        Some(file) => process_support_vcf(
            file,
            &mut consensus,
            &classifier,
            &processed_positions,
            verbose,
        )?,
        None => (Vec::new(), HashMap::new()),
    };
    if verbose {
        let num_process_positions = support_processed_positions
            .values()
            .map(|s| s.len())
            .sum::<usize>();
        println!(
            "From support VCF: Processed {} sites. Found {} records to output",
            num_process_positions,
            main_records.len(),
        );
    }

    main_records.extend(support_records);

    for (chrom, sites) in support_processed_positions.into_iter() {
        processed_positions
            .entry(chrom)
            .or_insert(HashSet::new())
            .extend(sites);
    }

    if params.mask_missing_sites {
        for (chrom, seq) in consensus.iter_mut() {
            let all_sites: HashSet<usize> = (0..seq.len()).collect();

            if let Some(processed_sites) = processed_positions.get(chrom) {
                for i in all_sites.difference(processed_sites) {
                    seq[*i] = NULL;
                }
            } else {
                for i in all_sites {
                    seq[i] = NULL;
                }
            }
        }
    }

    save_fasta(&consensus, &(output_root.to_owned() + ".full.fasta"))?;

    write_creation_report(
        &consensus,
        Some(het_sites.len() as i32),
        "genome_creation_report.json",
    )?;

    let mut clean_consensus = consensus.clone();
    clean_fasta_characters(&mut clean_consensus);
    save_fasta(&clean_consensus, &(output_root.to_owned() + ".fasta"))?;

    write_vcf(
        &main_records,
        &(output_root.to_owned() + ".vcf"),
        main_vcf,
        support_vcf,
    )?;

    Ok(consensus)
}

#[cfg(test)]
mod tests;

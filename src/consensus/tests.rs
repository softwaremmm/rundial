use super::*;
use crate::vcf::variant_record::tests::standard_header;

use pretty_assertions::assert_eq;

fn make_simple_chrom(seq: &str) -> HashMap<String, Vec<char>> {
    let mut chrom_seq: HashMap<String, Vec<char>> = HashMap::new();
    chrom_seq.insert("chrom".to_string(), seq.chars().collect());
    chrom_seq
}

fn test_apply(
    ref_seq: &str,
    _processed_sites: Vec<usize>,
    pos: usize,
    change: Change,
    ref_bases: &str,
    new_bases: &str,
) -> (String, HashSet<usize>) {
    let mut chrom_seq = make_simple_chrom(ref_seq);
    let mut processed_sites: HashMap<String, HashSet<usize>> = HashMap::new();
    processed_sites.extend_chrom("chrom", _processed_sites);

    let variant = Classification {
        pos,
        ref_bases: ref_bases.to_string(),
        new_bases: new_bases.to_string(),
        change,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_alleles: false,
    };
    let set_sites = apply_variant("chrom", &variant, &mut chrom_seq, &processed_sites);
    return (chrom_seq["chrom"].iter().collect(), set_sites);
}

#[test]
fn test_apply_null_variant() {
    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 0, Change::Null, "A", "N");
    assert_eq!(chrom_seq, "AAAAA");
    assert_eq!(set_sites.len(), 0);

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 3, Change::Null, "A", "N");
    assert_eq!(chrom_seq, "AAANA");
    assert_eq!(set_sites, HashSet::from([3]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 2, Change::Null, "AAA", "ZZZ");
    assert_eq!(chrom_seq, "AAAZZ");
    assert_eq!(set_sites, HashSet::from([3, 4]));
}

#[test]
fn test_apply_ref_variant() {
    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 0, Change::Ref, "A", "A");
    assert_eq!(chrom_seq, "AAAAA");
    assert_eq!(set_sites.len(), 0);

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 3, Change::Ref, "A", "A");
    assert_eq!(chrom_seq, "AAAAA");
    assert_eq!(set_sites, HashSet::from([3]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 2, Change::Ref, "AAA", "AAA");
    assert_eq!(chrom_seq, "AAAAA");
    assert_eq!(set_sites, HashSet::from([3, 4]));
}

#[test]
fn test_apply_snp_variant() {
    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1, 2], 0, Change::Snp, "A", "T");
    assert_eq!(chrom_seq, "TAAAA");
    assert_eq!(set_sites, HashSet::from([0]));
}

#[test]
fn test_apply_del_variant() {
    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 2, Change::Del, "AAA", "A");
    assert_eq!(chrom_seq, "AAA--");
    assert_eq!(set_sites, HashSet::from([2, 3, 4]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 2, Change::Del, "AA", "");
    assert_eq!(chrom_seq, "AA--A");
    assert_eq!(set_sites, HashSet::from([2, 3]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 0, Change::ComplexDel, "AA", "G");
    assert_eq!(chrom_seq, "G-AAA");
    assert_eq!(set_sites, HashSet::from([0, 1]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 0, Change::Mnp, "AA", "GT");
    assert_eq!(chrom_seq, "GTAAA");
    assert_eq!(set_sites, HashSet::from([0, 1]));
}

#[test]
fn test_apply_ins_variant() {
    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 0, Change::Ins, "A", "ATT");
    assert_eq!(chrom_seq, "ATTAAAA");
    assert_eq!(set_sites, HashSet::from([]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 1, Change::Ins, "", "TT");
    assert_eq!(chrom_seq, "ATTAAAA");
    assert_eq!(set_sites, HashSet::from([]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 0, Change::Ins, "", "TT");
    assert_eq!(chrom_seq, "TTAAAAA");
    assert_eq!(set_sites, HashSet::from([]));

    let (chrom_seq, set_sites) = test_apply("AAAAA", vec![0, 1], 1, Change::ComplexIns, "A", "GT");
    assert_eq!(chrom_seq, "AGTAAA");
    assert_eq!(set_sites, HashSet::from([]));
}

#[test]
fn test_check_indel_ref_matches_seq() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().collect();

    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "A", false),
        true
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "AAA", false),
        true
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "T", false),
        false
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "ATA", false),
        false
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "TAA", false),
        false
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "TAA", true),
        true
    );

    // Nulls are allowed
    let chrom_seq: Vec<char> = "ANNAA".chars().collect();
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "ATA", false),
        true
    );
    assert_eq!(
        _check_indel_ref_matches_seq(&chrom_seq, 0, "ATAT", false),
        false
    );
}

#[test]
fn test_overlap() {
    fn make_classification(pos: usize, ref_bases: &str, change: Change) -> Classification {
        Classification {
            pos,
            ref_bases: ref_bases.to_string(),
            new_bases: ref_bases.to_string(),
            change,
            is_het: false,
            has_minor_population: false,
            is_filtered: false,
            has_indel_alleles: false,
        }
    }
    let simple_changes = [Change::Ref, Change::Snp, Change::Null, Change::Mnp];
    let ins_changes = [Change::Ins, Change::ComplexIns];
    let del_changes = [Change::Del, Change::ComplexDel];
    for c1 in simple_changes.iter().chain(del_changes.iter()) {
        for c2 in simple_changes.iter().chain(del_changes.iter()) {
            println!("{:?} {:?}", c1, c2);
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(0, "A", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(0, "AA", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(1, "A", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(2, "A", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "AA", c1.clone()),
                    &make_classification(2, "A", c2.clone())
                ),
                true
            );
        }
    }
    for c1 in ins_changes.iter() {
        for c2 in simple_changes.iter() {
            assert_eq!(
                overlaps(
                    &make_classification(0, "A", c1.clone()),
                    &make_classification(0, "", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(0, "A", c1.clone()),
                    &make_classification(0, "A", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(0, "A", c1.clone()),
                    &make_classification(1, "", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(0, "A", c1.clone()),
                    &make_classification(1, "A", c2.clone())
                ),
                false
            );

            assert_eq!(
                overlaps(
                    &make_classification(0, "AC", c1.clone()),
                    &make_classification(1, "T", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(0, "AC", c1.clone()),
                    &make_classification(1, "", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(0, "AC", c1.clone()),
                    &make_classification(2, "", c2.clone())
                ),
                false
            );
        }
    }
    for c1 in ins_changes.iter() {
        for c2 in ins_changes.iter() {
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(0, "", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(0, "A", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(1, "", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(2, "", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(2, "", c2.clone())
                ),
                true
            );
        }
    }
    for c1 in ins_changes.iter() {
        for c2 in del_changes.iter() {
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(0, "A", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(1, "A", c2.clone())
                ),
                true
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "", c1.clone()),
                    &make_classification(2, "A", c2.clone())
                ),
                false
            );
            assert_eq!(
                overlaps(
                    &make_classification(1, "A", c1.clone()),
                    &make_classification(2, "A", c2.clone())
                ),
                true
            );
        }
    }
}

fn make_classifier() -> Classifier {
    let params = ConsensusParams {
        skip_indels: false,
        mask: None,
        mask_missing_sites: true,
        het_snp_option: HetOption::Mask,
        het_indel_option: HetOption::Mask,
        use_filters: true,
        filter_ignore_list: None,
        het_pc_threshold: None,
        minor_pop_threshold: Some(5),
        main_caller: None,
        support_caller: None,
    };
    return Classifier::new(&params);
}

#[test]
fn test_score_variants() {
    let header = standard_header();
    let classifier = make_classifier();
    let record_to_score = |r: &str| {
        let mut record = VariantRecord::from_string(&header, r).unwrap();
        let classification = classifier.classify(&mut record);
        return score_variant(&record, &classification);
    };
    let filter_score =
        record_to_score("ref\t1\tid\tTCG\tTAC\t244.589\tFILTER\tDP=28\tGT:AD\t./.:0,28");
    let null_score = record_to_score("ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t./.:0,28");
    let ref_indel_score =
        record_to_score("ref\t1\tid\tTCG\tT\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,0");
    let ref_snp_score = record_to_score("ref\t1\tid\tC\tT\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,0");
    let indel_score = record_to_score("ref\t1\tid\tT\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let snp_score = record_to_score("ref\t1\tid\tT\tC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let high_qual_score =
        record_to_score("ref\t1\tid\tT\tC\t2440.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let high_dp_score =
        record_to_score("ref\t1\tid\tT\tC\t2440.589\tPASS\tDP=280\tGT:AD\t1/1:28,0");
    assert!(filter_score < null_score);
    assert!(null_score < ref_indel_score);
    assert!(ref_indel_score < ref_snp_score);
    assert!(ref_snp_score < indel_score);
    assert!(indel_score < snp_score);
    assert!(snp_score < high_qual_score);
    assert!(high_qual_score < high_dp_score);

    let het_snp_score = record_to_score("ref\t1\tid\tT\tC\t244.589\tPASS\tDP=28\tGT:AD\t0/1:14,14");
    assert!(het_snp_score == snp_score);
}

#[test]
fn test_mark_overlaps() {
    let header = standard_header();
    let classifier = make_classifier();
    let record_to_pair = |r: &str| {
        let mut record = VariantRecord::from_string(&header, r).unwrap();
        let classification = classifier.classify(&mut record);
        return (record, classification);
    };

    let mut records: Vec<(VariantRecord, Classification)> = vec![
        record_to_pair("ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:0,28"),
        record_to_pair("ref\t2\tid\tCG\tC\t244.589\tF\tDP=28\tGT:AD\t1/1:0,28"),
        record_to_pair("ref\t2\tid\tCG\tA\t300.589\tPASS\tDP=28\tGT:AD\t1/1:0,28"),
        record_to_pair("ref\t5\tid\tCG\tC\t244.589\tF\tDP=28\tGT:AD\t1/1:0,28"),
    ];

    mark_overlaps(&mut records, &classifier);
    assert!(records[0].0.filter.contains(&OVERLAP_FILTER.to_owned()));
    assert!(records[0].1.is_filtered);
    assert_eq!(
        records[0].1,
        Classification {
            pos: 0,
            ref_bases: "TCG".to_owned(),
            new_bases: "FFF".to_owned(),
            change: Change::Null,
            is_het: false,
            has_minor_population: false,
            is_filtered: true,
            has_indel_alleles: true,
        }
    );

    assert!(records[1].0.filter.contains(&OVERLAP_FILTER.to_owned()));
    assert!(records[1].1.is_filtered);
    assert!(records[1].1.new_bases == "FF");

    assert!(!records[2].0.filter.contains(&OVERLAP_FILTER.to_owned()));
    assert!(!records[2].0.filter.contains(&OVERLAP_FILTER.to_owned()));
}

#[test]
fn test_read_write_fasta() {
    let seq = read_fasta("test_data/simple.fasta").unwrap();
    assert_eq!(
        seq,
        HashMap::from([
            ("chrom_1".to_string(), "AAAAANFFZZZMMMM-X".chars().collect()),
            ("chrom_2".to_string(), "CCCCC".chars().collect())
        ])
    );

    save_fasta(&seq, "tests/test_outputs/saved.fasta").unwrap();
    let seq2 = read_fasta("tests/test_outputs/saved.fasta").unwrap();
    assert_eq!(seq, seq2);
}

#[test]
fn test_clean_fasta_characters() {
    let mut seq = read_fasta("test_data/simple.fasta").unwrap();
    clean_fasta_characters(&mut seq);
    assert_eq!(
        seq,
        HashMap::from([
            ("chrom_1".to_string(), "AAAAANNNNNNNNNN-X".chars().collect()),
            ("chrom_2".to_string(), "CCCCC".chars().collect())
        ])
    );
}

#[test]
fn test_write_creation_report() {
    let seq = read_fasta("test_data/simple.fasta").unwrap();
    let het_count: Option<i32> = None;

    write_creation_report(&seq, het_count, "tests/test_outputs/creation_report_1.json").unwrap();
    let report = std::fs::read_to_string("tests/test_outputs/creation_report_1.json").unwrap();
    let expected_report = std::fs::read_to_string("test_data/creation_report.json").unwrap();
    let report_json: serde_json::Value = serde_json::from_str(&report).unwrap();
    let expected_report_json: serde_json::Value = serde_json::from_str(&expected_report).unwrap();
    assert_eq!(report_json, expected_report_json);

    let het_count: Option<i32> = Some(5);
    write_creation_report(&seq, het_count, "tests/test_outputs/creation_report_2.json").unwrap();
    let report = std::fs::read_to_string("tests/test_outputs/creation_report_2.json").unwrap();
    let report_json: serde_json::Value = serde_json::from_str(&report).unwrap();
    println!("{:?}", report_json);
    assert_eq!(
        report_json
            .get("Sequencing Quality")
            .unwrap()
            .get("Mixed calls")
            .unwrap(),
        5
    );
}

#[test]
fn test_write_vcf() {
    let header = standard_header();
    let records: Vec<VariantRecord> = [
        "ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:0,28",
        "ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28;CALLER=bcftools\tGT:AD\t1/1:0,28",
    ]
    .iter()
    .map(|r| VariantRecord::from_string(&header, r).unwrap())
    .collect();

    write_vcf(
        &records,
        "tests/test_outputs/write_vcf.vcf",
        "test_data/example.vcf",
        None,
    )
    .unwrap();
    let output = std::fs::read_to_string("tests/test_outputs/write_vcf.vcf").unwrap();
    let expected_output = std::fs::read_to_string("test_data/write_vcf.vcf").unwrap();
    assert_eq!(output, expected_output);
}

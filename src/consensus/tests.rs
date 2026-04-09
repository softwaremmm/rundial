use std::{fs, io::Read};

use super::*;
use crate::{
    consensus::parameter_struct::MinorPopParams, vcf::variant_record::tests::standard_header,
};

use pretty_assertions::assert_eq;

fn make_simple_chrom(seq: &str) -> HashMap<String, Vec<char>> {
    let mut chrom_seq: HashMap<String, Vec<char>> = HashMap::new();
    chrom_seq.insert("chrom".to_string(), seq.chars().collect());
    chrom_seq
}

fn test_apply(
    ref_seq: &str,
    pos: usize,
    change: Change,
    ref_bases: &str,
    new_bases: &str,
) -> String {
    let mut chrom_seq = make_simple_chrom(ref_seq);

    let variant = Classification {
        pos,
        ref_bases: ref_bases.to_string(),
        new_bases: new_bases.to_string(),
        change,
        ..Default::default()
    };
    apply_variant_simple("chrom", &variant, &mut chrom_seq);
    return chrom_seq["chrom"].iter().collect();
}

#[test]
fn test_apply_null_variant() {
    let chrom_seq = test_apply("AAAAA", 0, Change::Null, "A", "N");
    assert_eq!(chrom_seq, "NAAAA");

    let chrom_seq = test_apply("AAAAA", 3, Change::Null, "A", "N");
    assert_eq!(chrom_seq, "AAANA");

    let chrom_seq = test_apply("AAAAA", 2, Change::Null, "AAA", "FFF");
    assert_eq!(chrom_seq, "AAFFF");
}

#[test]
fn test_apply_ref_variant() {
    let chrom_seq = test_apply("AAAAA", 0, Change::Ref, "A", "A");
    assert_eq!(chrom_seq, "AAAAA");

    let chrom_seq = test_apply("AAAAA", 3, Change::Ref, "A", "A");
    assert_eq!(chrom_seq, "AAAAA");

    let chrom_seq = test_apply("AAAAA", 2, Change::Ref, "AAA", "AAA");
    assert_eq!(chrom_seq, "AAAAA");

    // empty change
    let chrom_seq = test_apply("AAAAA", 2, Change::Ref, "", "");
    assert_eq!(chrom_seq, "AAAAA");
}

#[test]
fn test_apply_snp_variant() {
    let chrom_seq = test_apply("AAAAA", 0, Change::Snp, "A", "T");
    assert_eq!(chrom_seq, "TAAAA");
}

#[test]
fn test_apply_del_variant() {
    let chrom_seq = test_apply("AAAAA", 2, Change::Del, "AAA", "A");
    assert_eq!(chrom_seq, "AAA--");

    let chrom_seq = test_apply("AAAAA", 2, Change::Del, "AA", "");
    assert_eq!(chrom_seq, "AA--A");

    let chrom_seq = test_apply("AAAAA", 0, Change::ComplexDel, "AA", "G");
    assert_eq!(chrom_seq, "G-AAA");

    let chrom_seq = test_apply("AAAAA", 0, Change::Mnp, "AA", "GT");
    assert_eq!(chrom_seq, "GTAAA");
}

#[test]
fn test_apply_ins_variant() {
    let chrom_seq = test_apply("AAAAA", 0, Change::Ins, "A", "ATT");
    assert_eq!(chrom_seq, "ATTAAAA");

    let chrom_seq = test_apply("AAAAA", 1, Change::Ins, "", "TT");
    assert_eq!(chrom_seq, "ATTAAAA");

    let chrom_seq = test_apply("AAAAA", 0, Change::Ins, "", "TT");
    assert_eq!(chrom_seq, "TTAAAAA");

    let chrom_seq = test_apply("AAAAA", 1, Change::ComplexIns, "A", "GT");
    assert_eq!(chrom_seq, "AGTAAA");
}

#[test]
fn test_overlap() {
    fn make_classification(pos: usize, ref_bases: &str, change: Change) -> Classification {
        Classification {
            pos,
            ref_bases: ref_bases.to_string(),
            new_bases: ref_bases.to_string(),
            change,
            ..Default::default()
        }
    }
    let simple_changes = [Change::Ref, Change::Snp, Change::Null, Change::Mnp];
    let ins_changes = [Change::Ins, Change::ComplexIns];
    let del_changes = [Change::Del, Change::ComplexDel];
    for c1 in simple_changes.iter().chain(del_changes.iter()) {
        for c2 in simple_changes.iter().chain(del_changes.iter()) {
            println!("{c1:?} {c2:?}");
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
        mask_missing_sites: true,
        use_filters: true,
        filter_ignore_list: None,
        minor_pop_thresholds: Some(MinorPopParams {
            threshold: 5,
            strand_bias: None,
            min_frs: None,
        }),
        remove_sub_minor_pops: false,
        call_snps_in_support: false,
        remove_minor_pops_in_support: false,
        main_caller: None,
        support_caller: None,
    };
    return Classifier::new(&params, false);
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

    let indel = record_to_score("ref\t1\tid\tT\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let filter_indel =
        record_to_score("ref\t1\tid\tTCG\tTAC\t244.589\tFILTER\tDP=28\tGT:AD\t1/1:0,28");

    let snp = record_to_score("ref\t1\tid\tT\tC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let ref_snp = record_to_score("ref\t1\tid\tC\tT\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,0");

    let ref_indel = record_to_score("ref\t1\tid\tTCG\tT\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,0");
    let minor_indel = record_to_score("ref\t1\tid\tTCG\tT\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,5");
    let filter_minor_indel =
        record_to_score("ref\t1\tid\tTCG\tT\t244.589\tFILTER\tDP=28\tGT:AD\t0/0:28,5");

    //het snp will get filter flag
    let het_snp = record_to_score("ref\t1\tid\tT\tC\t244.589\tPASS\tDP=28\tGT:AD\t0/1:14,14");
    let filter_snp = record_to_score("ref\t1\tid\tT\tC\t244.589\tFILTER\tDP=28\tGT:AD\t0/0:0,28");
    let filter_ref = record_to_score("ref\t1\tid\tT\t.\t244.589\tFILTER\tDP=28\tGT:AD\t0/0:1");

    let null = record_to_score("ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t./.:0,28");

    assert!(indel > filter_indel);
    assert!(filter_indel > snp);
    assert!(snp == ref_snp);
    assert!(ref_snp > ref_indel);
    assert!(ref_indel == minor_indel);
    assert!(minor_indel > filter_minor_indel);
    assert!(filter_minor_indel > filter_snp);
    assert!(filter_snp == het_snp);
    assert!(filter_snp == filter_ref);
    assert!(filter_ref > null);

    let std = record_to_score("ref\t1\tid\tT\tC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let high_qual = record_to_score("ref\t1\tid\tT\tC\t2440.589\tPASS\tDP=28\tGT:AD\t1/1:28,0");
    let high_dp = record_to_score("ref\t1\tid\tT\tC\t2440.589\tPASS\tDP=280\tGT:AD\t1/1:28,0");

    assert!(std < high_qual);
    assert!(high_qual < high_dp);
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
            new_bases: "NNN".to_owned(),
            change: Change::Null,
            is_filtered: true,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    assert!(records[1].0.filter.contains(&OVERLAP_FILTER.to_owned()));
    assert!(records[1].1.is_filtered);
    assert!(records[1].1.new_bases == "NN");

    assert!(!records[2].0.filter.contains(&OVERLAP_FILTER.to_owned()));
    assert!(!records[2].0.filter.contains(&OVERLAP_FILTER.to_owned()));
}

#[test]
fn test_read_write_fasta() {
    let seq = read_fasta("test_data/simple.fasta").unwrap();
    assert_eq!(
        seq,
        HashMap::from([
            ("chrom_1".to_string(), "AAAAANNNNNNNNNN-N".chars().collect()),
            ("chrom_2".to_string(), "CCCCC".chars().collect())
        ])
    );

    save_fasta(&seq, "tests/test_outputs/saved.fasta").unwrap();
    let seq2 = read_fasta("tests/test_outputs/saved.fasta").unwrap();
    assert_eq!(seq, seq2);
}

#[test]
fn test_read_write_fasta_multiline() {
    let seq = read_fasta("test_data/multiline.fasta").unwrap();
    assert_eq!(
        seq,
        HashMap::from([
            ("chrom_1".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".chars().collect()),
        ])
    );

    save_fasta(&seq, "tests/test_outputs/saved_multiline.fasta").unwrap();
    // Check that the saved file is the same as the original from actual lines
    let original = fs::read_to_string("test_data/multiline.fasta").unwrap();
    let saved = fs::read_to_string("tests/test_outputs/saved_multiline.fasta").unwrap();
    assert_eq!(original, saved);
}

#[test]
fn test_read_write_fasta_gzipped() {
    let seq = read_fasta("test_data/simple.fasta.gz").unwrap();
    assert_eq!(
        seq,
        HashMap::from([
            ("chrom_1".to_string(), "AAAAANNNNNNNNNN-N".chars().collect()),
            ("chrom_2".to_string(), "CCCCC".chars().collect())
        ])
    );

    save_fasta(&seq, "tests/test_outputs/saved.fasta.gz").unwrap();
    let seq2 = read_fasta("tests/test_outputs/saved.fasta.gz").unwrap();
    assert_eq!(seq, seq2);

    let mut file = File::open("tests/test_outputs/saved.fasta.gz").unwrap();
    let mut buffer = [0; 3];
    file.read_exact(&mut buffer).unwrap();
    assert!(buffer == [0x1f, 0x8b, 0x08]);
}

#[test]
fn test_write_creation_report() {
    let seq = read_fasta("test_data/simple.fasta").unwrap();

    fn simple_minor(chrom: &str, pos: u32, is_indel: bool) -> VariantRecord {
        let mut record = VariantRecord::empty_record();
        record.chrom = chrom.to_string();
        record.pos = pos;
        if is_indel {
            record.ref_bases = "A".to_string();
            record.alt = vec!["AT".to_string()];
        } else {
            record.ref_bases = "A".to_string();
            record.alt = vec!["C".to_string()];
        }
        record
    }

    // mock data does not actually match fasta
    let minors = Vec::from([
        simple_minor("chrom1", 1, false),
        simple_minor("chrom1", 5, false),
        simple_minor("chrom1", 6, true),
        simple_minor("chrom2", 1, false),
    ]);

    write_creation_report(
        &seq,
        &minors,
        12,
        "tests/test_outputs/creation_report_1.json",
    )
    .unwrap();
    let report = std::fs::read_to_string("tests/test_outputs/creation_report_1.json").unwrap();
    let expected_report = std::fs::read_to_string("test_data/creation_report.json").unwrap();
    let report_json: serde_json::Value = serde_json::from_str(&report).unwrap();
    let expected_report_json: serde_json::Value = serde_json::from_str(&expected_report).unwrap();
    assert_eq!(report_json, expected_report_json);
}

#[test]
fn test_simplify_and_write_vcf() {
    let header = standard_header();
    let records: Vec<VariantRecord> = [
        "ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:0,28",
        "ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28;CALLER=bcftools\tGT:AD\t1/1:0,28",
    ]
    .iter()
    .map(|r| VariantRecord::from_string(&header, r).unwrap())
    .map(|mut r| {
        simplify_record(&mut r);
        r
    })
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

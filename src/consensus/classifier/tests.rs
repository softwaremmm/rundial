use super::*;
use crate::vcf::variant_record::tests::standard_header;

#[test]
fn test_repeat_char() {
    assert_eq!(repeat_char('a', 3), "aaa");
    assert_eq!(repeat_char('a', 0), "");
}

#[test]
fn test_simplify_ref_alt() {
    assert_eq!(
        simplify_ref_alt("A", "N"),
        (0, "A".to_string(), "N".to_string())
    );
    assert_eq!(
        simplify_ref_alt("A", "T"),
        (0, "A".to_string(), "T".to_string())
    );
    assert_eq!(
        simplify_ref_alt("A", "TCG"),
        (0, "A".to_string(), "TCG".to_string())
    );
    assert_eq!(
        simplify_ref_alt("A", "ATCG"),
        (1, "".to_string(), "TCG".to_string())
    );
    assert_eq!(
        simplify_ref_alt("AC", "ACTCG"),
        (2, "".to_string(), "TCG".to_string())
    );
    assert_eq!(
        simplify_ref_alt("AC", "ACTCGC"),
        (1, "".to_string(), "CTCG".to_string())
    );
    assert_eq!(
        simplify_ref_alt("TTAC", "TTG"),
        (2, "AC".to_string(), "G".to_string())
    );
    assert_eq!(
        simplify_ref_alt("TTAC", "TTA"),
        (3, "C".to_string(), "".to_string())
    );
    assert_eq!(
        simplify_ref_alt("ATTTT", "ATT"),
        (1, "TT".to_string(), "".to_string())
    );
    assert_eq!(
        simplify_ref_alt("TTTTA", "TTA"),
        (1, "TT".to_string(), "".to_string())
    );
    assert_eq!(
        simplify_ref_alt("TTTTAC", "TTA"),
        (2, "TTAC".to_string(), "A".to_string())
    );
}

#[test]
fn test_change_from_ref_alt() {
    assert_eq!(Change::from_ref_alt("A", "N"), Change::Null);
    assert_eq!(Change::from_ref_alt("A", "AN"), Change::Null);
    assert_eq!(Change::from_ref_alt("A", "AF"), Change::Null);
    assert_eq!(Change::from_ref_alt("A", "ZZ"), Change::HetMask);
    assert_eq!(Change::from_ref_alt("A", "MA"), Change::Null);
    assert_eq!(Change::from_ref_alt("N", "N"), Change::Null);

    assert_eq!(Change::from_ref_alt("A", "A"), Change::Ref);
    assert_eq!(Change::from_ref_alt("AT", "AT"), Change::Ref);

    assert_eq!(Change::from_ref_alt("A", "T"), Change::Snp);

    assert_eq!(Change::from_ref_alt("A", "AT"), Change::Ins);
    assert_eq!(Change::from_ref_alt("A", "ATGC"), Change::Ins);
    assert_eq!(Change::from_ref_alt("AT", "ATGC"), Change::Ins);

    assert_eq!(Change::from_ref_alt("AT", "A"), Change::Del);
    assert_eq!(Change::from_ref_alt("AGCT", "A"), Change::Del);
    assert_eq!(Change::from_ref_alt("AGCT", "AG"), Change::Del);

    // We don't simplify indels here
    assert_eq!(Change::from_ref_alt("AT", "TT"), Change::Mnp);
    assert_eq!(Change::from_ref_alt("AT", "CG"), Change::Mnp);
    assert_eq!(Change::from_ref_alt("AT", "AG"), Change::Mnp);

    assert_eq!(Change::from_ref_alt("AT", "CGA"), Change::ComplexIns);
    assert_eq!(Change::from_ref_alt("ATA", "AG"), Change::ComplexDel);
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
        overriding_filters: None,
        minor_pop_threshold: Some(5),
        support_minor_pop_threshold: None,
        minor_pop_strand_bias: Some(0.05),
        main_caller: None,
        support_caller: None,
    };
    return Classifier::new(&params, false);
}

#[test]
fn test_get_flags() {
    let mut c = make_classifier();
    let mut record = VariantRecord::empty_record();

    record.filter = vec!["PASS".to_string()];
    assert_eq!(c.get_flags(&record), Vec::<String>::new());

    record.filter = vec!["A", "B", ".", "RefCall"]
        .into_iter()
        .map(|x| x.to_string())
        .collect();
    assert_eq!(c.get_flags(&record), vec!["A".to_string(), "B".to_string()]);

    c.params.filter_ignore_list = Some(vec!["A".to_string()]);
    assert_eq!(c.get_flags(&record), vec!["B".to_string()]);

    c.params.use_filters = false;
    assert_eq!(c.get_flags(&record), Vec::<String>::new());
}

#[test]
fn test_is_masked() {
    let mut c = make_classifier();
    c.mask = Some(Bed::from_file("test_data/mask.bed"));
    let mut record = VariantRecord::empty_record();
    record.chrom = "other_chrom".to_string();
    record.pos = 1;
    assert_eq!(c.is_masked(&record), false);

    // remember that bed files are 0-based
    record.pos = 10;
    assert_eq!(c.is_masked(&record), false);
    record.pos = 11;
    assert_eq!(c.is_masked(&record), true);
}

#[test]
fn test_classify_simple() {
    let header = standard_header();
    let c = make_classifier();

    // ref case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\t.\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "T".to_string(),
            change: Change::Ref,
            ..Default::default()
        }
    );

    // Ref indel case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tTA\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,1",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize + 1,
            ref_bases: "".to_string(),
            new_bases: "".to_string(),
            change: Change::Ref,
            has_indel_alleles: true,
            ..Default::default()
        }
    );
    // Ref indel non standard case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tCA\t244.589\tPASS\tDP=28\tGT:AD\t0/0:28,1",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "T".to_string(),
            change: Change::Ref,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    // null case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\t.\t244.589\tPASS\tDP=1\tGT:AD\t./.:1",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "N".to_string(),
            change: Change::Null,
            ..Default::default()
        }
    );
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tTTT\tA\t244.589\tPASS\tDP=1\tGT:AD\t./.:1,1",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "TTT".to_string(),
            new_bases: "NNN".to_string(),
            change: Change::Null,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    // Masking not applied here
    let mut record = VariantRecord::from_string(
        &header,
        "other_chrom\t15\tid\tT\t.\t244.589\tF\tDP=28\tGT:AD\t0/0:28",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "F".to_string(),
            change: Change::Null,
            is_filtered: true,
            ..Default::default()
        }
    );

    // SNP case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA\t244.589\t.\tDP=28\tGT:AD\t1/1:1,27",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "A".to_string(),
            change: Change::Snp,
            ..Default::default()
        }
    );

    // indel case, note that variations on this are handled in tests for Change::from_ref_alt
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tTA\t244.589\t.\tDP=28\tGT:AD\t1/1:1,27",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: 1,
            ref_bases: "".to_string(),
            new_bases: "A".to_string(),
            change: Change::Ins,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tTCGA\tTCC\t244.589\t.\tDP=28\tGT:AD\t1/1:1,27",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: record.pos_idx() + 2,
            ref_bases: "GA".to_string(),
            new_bases: "C".to_string(),
            change: Change::ComplexDel,
            has_indel_alleles: true,
            ..Default::default()
        }
    );
}

#[test]
fn test_classify_het() {
    let header = standard_header();
    let mut c = make_classifier();

    // Masking
    let mut snp_record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA\t244.589\t.\tDP=28\tGT:AD\t0/1:1,27",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut snp_record),
        Classification {
            pos: (snp_record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "Z".to_string(),
            change: Change::HetMask,
            is_het: true,
            ..Default::default()
        }
    );

    let mut indel_record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tTAA\tT\t244.589\t.\tDP=28\tGT:AD\t0/1:1,27",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut indel_record),
        Classification {
            pos: (indel_record.pos - 1) as usize + 1,
            ref_bases: "AA".to_string(),
            new_bases: "ZZ".to_string(),
            change: Change::HetMask,
            is_het: true,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    c.params.het_snp_option = HetOption::Ref;
    c.params.het_indel_option = HetOption::Ref;
    assert_eq!(c.classify(&mut snp_record).new_bases, "T".to_string());
    assert_eq!(c.classify(&mut snp_record).is_het, true);
    assert_eq!(c.classify(&mut indel_record).new_bases, "TAA".to_string());
    assert_eq!(c.classify(&mut indel_record).is_het, true);

    c.params.het_snp_option = HetOption::Alt;
    c.params.het_indel_option = HetOption::Alt;
    assert_eq!(c.classify(&mut snp_record).new_bases, "A".to_string());
    assert_eq!(c.classify(&mut snp_record).is_het, true);
    assert_eq!(c.classify(&mut indel_record).ref_bases, "AA".to_string());
    assert_eq!(c.classify(&mut indel_record).new_bases, "".to_string());
    assert_eq!(c.classify(&mut indel_record).is_het, true);

    c.params.het_snp_option = HetOption::Best;
    c.params.het_indel_option = HetOption::Best;
    assert_eq!(c.classify(&mut snp_record).new_bases, "A".to_string());
    assert_eq!(c.classify(&mut snp_record).is_het, true);
    assert_eq!(c.classify(&mut indel_record).ref_bases, "AA".to_string());
    assert_eq!(c.classify(&mut indel_record).new_bases, "".to_string());
    assert_eq!(c.classify(&mut indel_record).is_het, true);
}

#[test]
fn test_classify_filtered() {
    let header = standard_header();
    let mut c = make_classifier();

    // Standard filter
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\t.\t244.589\tF\tDP=28\tGT:AD\t0/0:28",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "F".to_string(),
            change: Change::Null,
            is_filtered: true,
            ..Default::default()
        }
    );

    // indel filter
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tTAA\tTA\t244.589\tF\tDP=28\tGT:AD\t0/0:28,1",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "TAA".to_string(),
            new_bases: "FFF".to_string(),
            change: Change::Null,
            is_filtered: true,
            has_indel_alleles: true,
            ..Default::default()
        }
    );

    c.params.use_filters = false;
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\t.\t244.589\tF\tDP=28\tGT:AD\t0/0:28",
    )
    .unwrap();
    assert_eq!(
        c.classify(&mut record),
        Classification {
            pos: (record.pos - 1) as usize,
            ref_bases: "T".to_string(),
            new_bases: "T".to_string(),
            change: Change::Ref,
            ..Default::default()
        }
    );
}

#[test]
#[should_panic]
fn test_classify_half_null_genotype() {
    let header = standard_header();
    let c = make_classifier();
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\t.\t244.589\tPASS\tDP=1\tGT:AD\t./1:1,0",
    )
    .unwrap();
    c.classify(&mut record);
}

#[test]
fn test_classify_minor_population() {
    let header = standard_header();
    let mut c = make_classifier();

    // standard case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t1/1:1,27,5",
    )
    .unwrap();
    assert!(c.classify(&mut record).has_minor_population);

    // not affected by filter
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tF\tDP=28\tGT:AD\t1/1:1,27,5",
    )
    .unwrap();
    assert!(c.classify(&mut record).has_minor_population);

    // Het case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t0/1:1,27,5",
    )
    .unwrap();
    assert!(c.classify(&mut record).has_minor_population);
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t1/2:1,27,5",
    )
    .unwrap();
    assert!(!c.classify(&mut record).has_minor_population);

    // null gt case
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t./.:1,27,5",
    )
    .unwrap();
    assert!(c.classify(&mut record).has_minor_population);

    // too low
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t1/1:1,27,4",
    )
    .unwrap();
    assert!(!c.classify(&mut record).has_minor_population);

    // change threshold
    c.minor_pop_threshold = Some(4);
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t1/1:1,27,4",
    )
    .unwrap();
    assert!(c.classify(&mut record).has_minor_population);

    // if no threshold
    c.minor_pop_threshold = None;
    let mut record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA,C\t244.589\tPASS\tDP=28\tGT:AD\t1/1:1,27,4",
    )
    .unwrap();
    assert!(!c.classify(&mut record).has_minor_population);

    // minor allele strand bias
    c.minor_pop_threshold = Some(2);
    let record = VariantRecord::from_string(
        &header,
        "ref\t1\tid\tT\tA\t244.589\tPASS\tDP=101;ADF=50,1;ADR=45,4\tGT:AD\t0/0:95,5",
    )
    .unwrap();

    // Pass strand bias test
    c.params.minor_pop_strand_bias = Some(0.20);
    let mut tmp_record = record.clone();
    assert!(c.classify(&mut tmp_record).has_minor_population);
    assert_eq!(
        tmp_record.to_string(),
        "ref\t1\tid\tT\tA\t244.589\tPASS\tDP=101;ADF=50,1;ADR=45,4\tGT:AD\t0/0:95,5"
    );

    // Fail strand bias test
    c.params.minor_pop_strand_bias = Some(0.21);
    let mut tmp_record = record.clone();
    assert!(!c.classify(&mut tmp_record).has_minor_population);
    assert_eq!(tmp_record.to_string(),
        "ref\t1\tid\tT\t.\t244.589\tPASS\tDP=101;ADF=50;ADR=45;FILTERED_MINOR_ALLELES=1:4\tGT:AD\t0/0:95"
    );
}

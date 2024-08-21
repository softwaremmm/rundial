use super::*;
use crate::vcf::variant_record::tests::standard_header;

#[test]
fn test_apply_null_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1, 2]);

    // Null cases
    let mut null_change = Classification {
        pos: 0,
        ref_bases: "A".to_string(),
        new_bases: "N".to_string(),
        change: Change::Null,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&null_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 0);

    null_change.pos = 3;
    let set_sites = apply_variant(&null_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 1);
    assert_eq!(tmp_chrom_seq[3], 'N');

    null_change.pos = 2;
    null_change.ref_bases = "AAA".to_string();
    null_change.new_bases = "ZZZ".to_string();
    let set_sites = apply_variant(&null_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 2);
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAAZZ");
}

#[test]
fn test_apply_ref_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1, 2]);

    let mut ref_change = Classification {
        pos: 0,
        ref_bases: "A".to_string(),
        new_bases: "A".to_string(),
        change: Change::Ref,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&ref_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 0);

    ref_change.pos = 3;
    let set_sites = apply_variant(&ref_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 1);
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAAAA");

    ref_change.pos = 2;
    ref_change.ref_bases = "AAA".to_string();
    ref_change.new_bases = "AAA".to_string();
    let set_sites = apply_variant(&ref_change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites.len(), 2);
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAAAA");
}

#[test]
fn test_apply_snp_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1, 2]);

    let change = Classification {
        pos: 0,
        ref_bases: "A".to_string(),
        new_bases: "T".to_string(),
        change: Change::Snp,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites, HashSet::from([0]));
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "TAAAA");
}

#[test]
fn test_apply_del_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1]);

    let mut change = Classification {
        pos: 2,
        ref_bases: "AAA".to_string(),
        new_bases: "A".to_string(),
        change: Change::Del,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites, HashSet::from([2, 3, 4]));
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAA--");

    // A deletion will not "set" the first site if it is already set
    // This is useful if there has a snp say at that site.
    change.pos = 1;
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites, HashSet::from([2, 3]));
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AA--A");
}

#[test]
fn test_apply_ins_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1]);

    let change = Classification {
        pos: 2,
        ref_bases: "A".to_string(),
        new_bases: "ATT".to_string(),
        change: Change::Ins,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites, HashSet::from([]));
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAATTAA");
}

#[test]
fn test_apply_cmplx_variant() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();
    let processed_sites = HashSet::from([0, 1]);

    let change = Classification {
        pos: 2,
        ref_bases: "AAA".to_string(),
        new_bases: "AT".to_string(),
        change: Change::ComplexIndel,
        is_het: false,
        has_minor_population: false,
        is_filtered: false,
        has_indel_form: false,
    };
    let mut tmp_chrom_seq = chrom_seq.clone();
    let set_sites = apply_variant(&change, &mut tmp_chrom_seq, &processed_sites);
    assert_eq!(set_sites, HashSet::from([]));
    assert_eq!(&tmp_chrom_seq.iter().collect::<String>(), "AAAT");
}

#[test]
fn test_check_indel_ref_matches_seq() {
    let chrom_seq: Vec<char> = repeat_char('A', 5).chars().into_iter().collect();

    assert_eq!(check_indel_ref_matches_seq(&chrom_seq, 0, "A", false), true);
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "AAA", false),
        true
    );
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "T", false),
        false
    );
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "ATA", false),
        false
    );
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "TAA", false),
        false
    );
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "TAA", true),
        true
    );

    // Nulls are allowed
    let chrom_seq: Vec<char> = "ANNAA".chars().into_iter().collect();
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "ATA", false),
        true
    );
    assert_eq!(
        check_indel_ref_matches_seq(&chrom_seq, 0, "ATAT", false),
        false
    );
}

#[test]
fn test_filter_overlapping_indels() {
    let header = standard_header();
    let mut records: Vec<VariantRecord> = Vec::new();

    records.push(
        VariantRecord::from_string(
            &header,
            "ref\t1\tid\tTCG\tTAC\t244.589\tPASS\tDP=28\tGT:AD\t1/1:0,28",
        )
        .unwrap(),
    );
    records.push(
        VariantRecord::from_string(
            &header,
            "ref\t2\tid\tCG\tC\t244.589\tF\tDP=28\tGT:AD\t1/1:0,28",
        )
        .unwrap(),
    );
    records.push(
        VariantRecord::from_string(
            &header,
            "ref\t2\tid\tCG\tA\t300.589\tPASS\tDP=28\tGT:AD\t1/1:0,28",
        )
        .unwrap(),
    );
    records.push(
        VariantRecord::from_string(
            &header,
            "ref\t5\tid\tCG\tC\t244.589\tF\tDP=28\tGT:AD\t1/1:0,28",
        )
        .unwrap(),
    );

    filter_overlapping_indels(&mut records);
    assert!(records[0].filter.contains(&"OverlapWithIndel".to_owned()));
    assert!(records[1].filter.contains(&"OverlapWithIndel".to_owned()));
    assert!(!records[2].filter.contains(&"OverlapWithIndel".to_owned()));
    assert!(!records[2].filter.contains(&"OverlapWithIndel".to_owned()));
}

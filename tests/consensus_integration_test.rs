use pretty_assertions::assert_eq;
use rundial::consensus::make_consensus;
use std::fs::{create_dir_all, read_to_string};

#[test]
fn test_make_consensus_single() {
    let main_vcf = "test_data/single_consensus/consensus_input.vcf";
    let support_vcf = None;
    let ref_fasta = "test_data/single_consensus/consensus_ref.fasta";
    let output_root = "tests/test_outputs/single_consensus/output";
    let params = "test_data/single_consensus/consensus_params.yml";

    create_dir_all("tests/test_outputs/single_consensus").unwrap();

    make_consensus(main_vcf, support_vcf, ref_fasta, output_root, params, true).unwrap();

    let expected_root = "test_data/single_consensus/consensus_expected";
    for ending in &[".full.fasta", ".fasta", ".variable_length.fasta", ".vcf"] {
        let output = output_root.to_string() + ending;
        let expected = expected_root.to_string() + ending;
        println!("Comparing {} to {}", output, expected);
        assert_eq!(
            read_to_string(expected).unwrap(),
            read_to_string(output).unwrap()
        );
    }
}

#[test]
fn test_make_consensus_support() {
    let folder = "test_data/support_consensus/consensus".to_string();
    let main_vcf = folder.clone() + "_main.vcf";
    let support_vcf = folder.clone() + "_support.vcf";
    let ref_fasta = folder.clone() + "_ref.fasta";
    let params = folder.clone() + "_params.yml";
    let output_root = "tests/test_outputs/support_consensus/output";

    create_dir_all("tests/test_outputs/support_consensus").unwrap();

    make_consensus(
        &main_vcf,
        Some(&support_vcf),
        &ref_fasta,
        output_root,
        &params,
        true,
    )
    .unwrap();

    let expected_root = "test_data/support_consensus/consensus_expected";
    for ending in &[".full.fasta", ".fasta", ".variable_length.fasta", ".vcf"] {
        let output = output_root.to_string() + ending;
        let expected = expected_root.to_string() + ending;
        println!("Comparing {} to {}", output, expected);
        assert_eq!(
            read_to_string(expected).unwrap(),
            read_to_string(output).unwrap()
        );
    }
}

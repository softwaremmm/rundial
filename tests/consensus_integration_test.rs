use pretty_assertions::assert_eq;
use rundial::consensus::make_consensus;
use std::fs::{create_dir_all, read_to_string};

fn update_expectations() -> bool {
    std::env::var("UPDATE_EXPECTATIONS")
        .map(|val| val == "1" || val.to_lowercase() == "true")
        .unwrap_or(false)
}

fn compare_files(expected: &str, result: &str) {
    let expected_content = match read_to_string(expected) {
        Ok(content) => content,
        Err(_) => {
            println!("Failed to read expected file '{expected}'");
            "".to_string()
        }
    };
    let result_content = read_to_string(result).unwrap();
    let equal = expected_content == result_content;

    if !equal {
        println!("Files differ:\nExpected: {expected}\nResult: {result}");

        if update_expectations() {
            std::fs::write(expected, &result_content).expect("Failed to update expected file");
            println!("Updated expected file: {expected}");
        } else {
            assert_eq!(expected_content, result_content);
        }
    }
}

#[test]
fn test_make_consensus_single() {
    let main_vcf = "test_data/single_consensus/consensus_input.vcf";
    let support_vcf = None;
    let ref_fasta = "test_data/single_consensus/consensus_ref.fasta";
    let output_root = "tests/test_outputs/single_consensus/output";
    let params = "test_data/single_consensus/consensus_params.yml";

    create_dir_all("tests/test_outputs/single_consensus").unwrap();

    make_consensus(main_vcf, support_vcf, ref_fasta, output_root, params, true).unwrap();

    let expected_root = "test_data/single_consensus/expected";
    for ending in &[".fasta", ".variable_length.fasta", ".vcf", ".report.json"] {
        let output = output_root.to_string() + ending;
        let expected = expected_root.to_string() + ending;
        compare_files(&expected, &output);
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

    let expected_root = "test_data/support_consensus/expected";
    for ending in &[".fasta", ".variable_length.fasta", ".vcf"] {
        let output = output_root.to_string() + ending;
        let expected = expected_root.to_string() + ending;
        compare_files(&expected, &output);
    }
}

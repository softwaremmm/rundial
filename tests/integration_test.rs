use rundial::filter_vcf::filter_vcf;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

fn read_lines<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    reader.lines().collect()
}

#[test]
fn test_filter_vcf() {
    let vcf_path = "test_data/unfiltered.vcf";
    let params_path = "test_data/filter_params.yml";
    let output_path = "tests/test_outputs/test_filter_vcf.vcf";
    let expected_output_path = "test_data/filtered.vcf";

    let _ = PathBuf::from(vcf_path);
    let _ = PathBuf::from(output_path);
    let _ = PathBuf::from(params_path);

    filter_vcf(
        PathBuf::from(vcf_path),
        PathBuf::from(output_path),
        PathBuf::from(params_path),
        false,
        true,
    )
    .unwrap();

    let output_lines = read_lines(output_path).unwrap();
    let expected_output_lines = read_lines(expected_output_path).unwrap();

    if output_lines != expected_output_lines {
        println!("Output does not match expectation for filtered vcf");
        println!("expectation: {}, result: {}", expected_output_path, output_path);
        panic!();
    }
}

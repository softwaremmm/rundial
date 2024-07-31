use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::io::{self, BufWriter};
use std::path::PathBuf;

use noodles::vcf;
use noodles::vcf::header::record::value::{map::Filter, Map};

use phf::phf_map;

use crate::parameter_structs::FilterParams;

static DESCRIPTIONS: phf::Map<&'static str, &'static str> = phf_map! {
    "low_qual" => "Quality is less than ?",
    "high_qual" => "Quality is greater than ?",
};

struct Filterer {
    // a list of functions with the signature fn(&vcf::Record, f32) -> bool
    filters: Vec<Box<dyn Fn(&vcf::Record, f32) -> bool>>,

}

fn is_low_qual(record: &vcf::Record, threshold: f32) -> bool {
    if let Some(Ok(quality_score)) = record.quality_score() {
        return quality_score < threshold;
    }
    return false;
}

fn is_low_depth(record: &vcf::Record, threshold: f32) -> bool {
    return false;
}

fn add_filter_to_header(header: &mut vcf::Header, key: &str, desc: &str) {
    let filter = Map::<Filter>::new(desc);
    header.filters_mut().insert(key.to_string(), filter);
}

fn add_filters_to_header(header: &mut vcf::Header, params: &FilterParams) {
    let all_params: Vec<HashMap<String, f32>> = vec![
        params.parameters.clone(),
        params.ref_parameters.clone(),
        params.snp_parameters.clone(),
        params.indel_parameters.clone(),
    ]
    .into_iter()
    .filter_map(|x| x)
    .collect();

    let all_keys = all_params
        .iter()
        .flat_map(|x| x.keys())
        .collect::<HashSet<&String>>();

    for key in all_keys {
        let thresholds: String = all_params
            .iter()
            .filter(|x| x.contains_key(key))
            .map(|x| x[key].to_string())
            .collect::<Vec<String>>()
            .join(", ");

        let desc: String = DESCRIPTIONS.get(key)
            .map(|&desc| desc.to_string())
            .unwrap_or_else(|| format!("{} - thresholds: ?", key))
            .replace("?", &thresholds);

        add_filter_to_header(header, key, &desc);
    }
}

pub fn filter_vcf(
    in_vcf: PathBuf,
    out_vcf: PathBuf,
    params: PathBuf,
    overwrite: bool,
    verbose: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let params: FilterParams = serde_yaml::from_reader(File::open(params)?)?;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(in_vcf)?;
    let mut header = reader.read_header()?;
    add_filters_to_header(&mut header, &params);

    let mut writer = vcf::io::Writer::new(BufWriter::new(File::create(out_vcf)?));
    writer.write_header(&header)?;

    let low_qual_filter = Map::<Filter>::new("Quality below threshold");
    header
        .filters_mut()
        .insert("low_qual".to_string(), low_qual_filter);

    Ok(())
}

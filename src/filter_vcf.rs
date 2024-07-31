use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::Read;
use std::io::{self, BufWriter};
use std::path::PathBuf;

use noodles::vcf;
use noodles::vcf::header::record::value::{map::Filter, Map};
use noodles::vcf::variant::record::info::field::Value;

use phf::phf_map;

use crate::parameter_structs::FilterParams;

const LOW_DP: &'static str = "low_depth";
const LOW_QUAL: &'static str = "low_qual";

static DESCRIPTIONS: phf::Map<&'static str, &'static str> = phf_map! {
    "low_qual" => "Quality is less than ?",
    "low_dp" => "Depth is less than ?",
};

type FilterFunction = Box<dyn Fn(&vcf::Header, &vcf::Record) -> String>;

struct Filterer<'a> {
    // The list of filter functions to apply
    filters: Vec<FilterFunction>,
    header: &'a vcf::Header,
}

impl<'a> Filterer<'a> {
    fn add_filter(&mut self, filter: FilterFunction) {
        self.filters.push(filter);
    }

    fn create_from_params(header: &'a vcf::Header, params: &Option<HashMap<String, f32>>) -> Self {
        let mut filterer = Self {
            filters: Vec::new(),
            header,
        };
        let Some(params) = params else {
            return filterer;
        };

        let new_header = header.clone();
        if let Some(&threshold) = params.get("low_qual") {
            filterer.add_filter(Box::new(move |header, record| {
                is_low_qual(header, record, threshold)
            }));
        }
        if let Some(&threshold) = params.get("low_depth") {
            filterer.add_filter(Box::new(move |header, record| {
                is_low_depth(header, record, threshold)
            }));
        }

        filterer
    }

    fn filter(&self, record: &vcf::Record) -> String {
        self.filters
            .iter()
            .map(|f| f(self.header, record))
            .filter(|s| s != "")
            .collect::<Vec<String>>()
            .join(";")
    }
}

fn is_low_qual(_: &vcf::Header, record: &vcf::Record, threshold: f32) -> String {
    if let Some(Ok(quality_score)) = record.quality_score() {
        if quality_score < threshold {
            return "low_qual".to_string();
        }
    }
    return "".to_string();
}

fn is_low_depth(header: &vcf::Header, record: &vcf::Record, threshold: f32) -> String {
    if let Some(Ok(Some(value))) = record.info().get(header, "DP") {
        if let Value::Integer(dp) = value {
            if dp < threshold as i32 {
                return "low_depth".to_string();
            }
        }
    }
    return "".to_string();
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

        let desc: String = DESCRIPTIONS
            .get(key)
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

    let std_filterer = Filterer::create_from_params(&header, &params.parameters);

    for result in reader.records() {
        let record: vcf::Record = result?;
        let filter_str = std_filterer.filter(&record);
        println!("{}", filter_str);
        let filters = record.filters();
        filters.

        writer.write_record(&header, &record)?;

    }

    Ok(())
}

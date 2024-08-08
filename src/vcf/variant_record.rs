use core::fmt;

use indexmap::IndexMap;
use std::error::Error;

use crate::vcf::{VCFError, VCFHeader};

pub enum RecordValue {
    Flag,
    Integer(i32),
    Float(f32),
    String(String),
    IntegerArray(Vec<i32>),
    FloatArray(Vec<f32>),
    StringArray(Vec<String>),
    Missing, // This is a . in the VCF file
}

impl fmt::Display for RecordValue {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                RecordValue::Integer(i) => i.to_string(),
                RecordValue::Float(fl) => format!("{:?}", fl),
                RecordValue::Flag => String::from(""),
                RecordValue::IntegerArray(arr) => arr
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                RecordValue::FloatArray(arr) => arr
                    .iter()
                    .map(|f| f.to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                RecordValue::String(s) => s.to_string(),
                RecordValue::StringArray(arr) => arr.join(","),
                RecordValue::Missing => String::from("."),
            }
        )
    }
}

pub struct Genotype {
    pub allele1: i32,
    pub allele2: i32,
}
impl Default for Genotype {
    fn default() -> Self {
        Self::new()
    }
}
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match (self.allele1, self.allele2) {
            (-1, -1) => return write!(f, "./."),
            (a1, -1) => return write!(f, "{}/.", a1),
            (-1, a2) => return write!(f, "./{}", a2),
            (a1, a2) => return write!(f, "{a1}/{a2}"),
        }
    }
}
impl Genotype {
    pub fn new() -> Self {
        Self {
            allele1: -1,
            allele2: -1,
        }
    }

    pub fn is_het(&self) -> bool {
        if self.allele1 == -1 || self.allele2 == -1 {
            return false;
        }
        return self.allele1 != self.allele2;
    }

    pub fn is_hom(&self) -> bool {
        if self.allele1 == -1 || self.allele2 == -1 {
            return false;
        }
        return self.allele1 == self.allele2;
    }

    pub fn is_hom_ref(&self) -> bool {
        return self.allele1 == 0 && self.allele2 == 0;
    }

    pub fn from_string(s: &str) -> Result<Self, Box<dyn Error>> {
        let alleles: Vec<&str> = s.split('/').collect();
        if alleles.len() != 2 {
            return Err(
                VCFError::InvalidField(format!("Genotype field {s} could not be parsed")).into(),
            );
        }
        let allele1: i32 = match (alleles[0], alleles[0].parse::<i32>()) {
            (".", _) => -1,
            (_, Ok(i)) if i >= -1 => i,
            _ => {
                return Err(VCFError::InvalidField(format!(
                    "Genotype field {s} could not be parsed"
                ))
                .into())
            }
        };
        let allele2: i32 = match (alleles[1], alleles[1].parse::<i32>()) {
            (".", _) => -1,
            (_, Ok(i)) if i >= -1 => i,
            _ => {
                return Err(VCFError::InvalidField(format!(
                    "Genotype field {s} could not be parsed"
                ))
                .into())
            }
        };
        return Ok(Self { allele1, allele2 });
    }
}

pub struct VariantRecord {
    // Core values
    pub chrom: String,
    pub pos: u32,
    pub id: Option<String>,
    pub ref_bases: String,
    pub alt: Vec<String>,
    pub qual: Option<f32>,
    pub filter: Vec<String>,
    pub info: IndexMap<String, RecordValue>,
    pub format: IndexMap<String, RecordValue>,

    // Derived Values
    genotype: Option<Genotype>,
    pub depth: Option<i32>,
    pub allele_depths: Option<Vec<i32>>,
    pub strand_depths: Option<(Vec<i32>, Vec<i32>)>,
}

impl fmt::Display for VariantRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dot = String::from(".");

        let pos_str: String = self.pos.to_string();
        let id_str: &str = if let Some(id) = &self.id { id } else { &dot };

        let alt_str: String = if self.alt.is_empty() {
            String::from(".")
        } else {
            self.alt.join(",")
        };

        let qual_str: String = if let Some(qual) = self.qual {
            qual.to_string()
        } else {
            String::from(".")
        };

        let filter_str: String = if self.filter.is_empty() {
            String::from("PASS")
        } else {
            self.filter.join(";")
        };

        let info_str: String = self
            .info
            .iter()
            .map(|(k, v)| {
                if let RecordValue::Flag = v {
                    return k.to_string();
                }
                format!("{}={}", k, v)
            })
            .collect::<Vec<String>>()
            .join(";");

        let (keys, values): (Vec<&str>, Vec<String>) = self
            .format
            .iter()
            .map(|(k, v)| (k.as_str(), v.to_string()))
            .unzip();
        let key_str: String = keys.join(":");
        let fmt_value_str: String = values.join(":");

        write!(
            f,
            "{}",
            [
                &self.chrom,
                &pos_str,
                id_str,
                &self.ref_bases,
                &alt_str,
                &qual_str,
                &filter_str,
                &info_str,
                &key_str,
                &fmt_value_str,
            ]
            .join("\t")
        )
    }
}

impl fmt::Debug for VariantRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Record: {}", self)
    }
}

impl VariantRecord {
    pub fn empty_record() -> Self {
        Self {
            chrom: String::new(),
            pos: 0,
            id: None,
            ref_bases: String::new(),
            alt: Vec::new(),
            qual: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: IndexMap::new(),
            genotype: None,
            depth: None,
            allele_depths: None,
            strand_depths: None,
        }
    }

    pub fn from_string(header: &VCFHeader, line: &str) -> Result<Self, Box<dyn Error>> {
        // split line by \t
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 10 {
            return Err(VCFError::InvalidRecord(format!(
                "VCF row should have 10 fields. Row: {line}"
            ))
            .into());
        }

        let mut record: Self = Self {
            chrom: fields[0].to_string(),
            pos: fields[1].parse()?,
            id: if fields[2].is_empty() {
                None
            } else {
                Some(fields[2].to_string())
            },
            ref_bases: fields[3].to_string(),
            alt: if fields[4] == "." || fields[4].is_empty() {
                Vec::new()
            } else {
                fields[4].split(',').map(|s| s.to_string()).collect()
            },
            qual: if fields[5] == "." || fields[5].is_empty() {
                None
            } else {
                match fields[5].parse::<f32>() {
                    Err(e) => {
                        return Err(VCFError::InvalidField(format!(
                            "Could not parse QUAL={}.\nError: {e:?}",
                            fields[5]
                        ))
                        .into())
                    }
                    Ok(f) => Some(f),
                }
            },
            filter: match fields[6] {
                "." | "" | "PASS" => Vec::new(),
                _ => fields[6].split(';').map(|s| s.to_string()).collect(),
            },
            info: str_to_info(header, fields[7])?,
            format: str_to_format(header, fields[8], fields[9])?,
            genotype: None,
            depth: None,
            allele_depths: None,
            strand_depths: None,
        };

        // find genotypes
        if let Some(RecordValue::String(genotype_str)) = record.format.get("GT") {
            record.genotype = Some(Genotype::from_string(genotype_str)?);
        } else {
            return Err(VCFError::InvalidRecord(format!(
                "Genotype field not found in VCF row: {line}"
            ))
            .into());
        }

        // Calculate depths
        record.update_depths();

        return Ok(record);
    }

    pub fn genotype(&self) -> Option<&Genotype> {
        return self.genotype.as_ref();
    }

    pub fn set_genotype(&mut self, genotype: Genotype) {
        self.format["GT"] = RecordValue::String(genotype.to_string());
        self.genotype = Some(genotype);
    }

    pub fn main_allele(&self) -> i32 {
        if let Some(genotype) = &self.genotype {
            if genotype.is_hom() {
                return genotype.allele1;
            }
            if genotype.allele1 == -1 || genotype.allele2 == -1 {
                return std::cmp::max(genotype.allele1, genotype.allele2);
            }
            if let Some(depths) = &self.allele_depths {
                let dp1 = depths[genotype.allele1 as usize];
                let dp2 = depths[genotype.allele2 as usize];
                return if dp1 > dp2 {
                    genotype.allele1
                } else {
                    genotype.allele2
                };
            }
            return genotype.allele1;
        }
        return -1;
    }

    pub fn is_indel(&self) -> bool {
        // Indel if any allele is longer than 1. Even if GT = 0/0
        return self.alt.iter().any(|a| a.len() > 1) || self.ref_bases.len() > 1;
    }

    pub fn is_snp(&self) -> bool {
        return !self.alt.is_empty() && !self.is_indel();
    }

    pub fn update_depths(&mut self) {
        // Calculate depths
        if let Some(RecordValue::Integer(dp)) = self.info.get("DP") {
            self.depth = Some(*dp);
        }

        if let Some(forward) = self.info.get("ADF") {
            if let Some(reverse) = self.info.get("ADR") {
                if let (RecordValue::IntegerArray(f), RecordValue::IntegerArray(r)) =
                    (forward, reverse)
                {
                    self.strand_depths = Some((f.clone(), r.clone()));
                }
            }
        }

        if let Some(depths) = self.format.get("AD") {
            if let RecordValue::IntegerArray(arr) = depths {
                self.allele_depths = Some(arr.clone())
            }
        } else if let Some((forward, reverse)) = &self.strand_depths {
            let depths: Vec<i32> = forward
                .iter()
                .zip(reverse.iter())
                .map(|(f, r)| f + r)
                .collect();
            self.allele_depths = Some(depths);
        }
    }
}

fn str_to_info(
    header: &VCFHeader,
    line: &str,
) -> Result<IndexMap<String, RecordValue>, Box<dyn Error>> {
    let mut info: IndexMap<String, RecordValue> = IndexMap::new();
    if line == "." || line.is_empty() {
        return Ok(info);
    }
    for field in line.split(';') {
        if !field.contains('=') {
            info.insert(field.to_string(), header.parse_info_value(field, field)?);
            continue;
        }
        let mut key_value = field.splitn(2, '=');
        let (Some(key), Some(value)) = (key_value.next(), key_value.next()) else {
            return Err(VCFError::InvalidField(format!(
                "field {field} could not be parsed into a key value pair"
            ))
            .into());
        };
        info.insert(key.to_string(), header.parse_info_value(key, value)?);
    }
    return Ok(info);
}

fn str_to_format(
    header: &VCFHeader,
    keys_str: &str,
    values_str: &str,
) -> Result<IndexMap<String, RecordValue>, Box<dyn Error>> {
    let keys: Vec<String> = keys_str.split(':').map(|s| s.to_string()).collect();
    let values: Vec<&str> = values_str.split(':').collect();
    if keys.len() != values.len() {
        return Err(VCFError::InvalidRecord(format!(
            "Mismatching number of keys and values: {keys:?}, {values:?}"
        ))
        .into());
    }
    let parsed_values: Vec<RecordValue> = std::iter::zip(keys.iter(), values)
        .map(|(k, v)| header.parse_format_value(k, v))
        .collect::<Result<_, _>>()?;

    let format: IndexMap<String, RecordValue> =
        IndexMap::from_iter(std::iter::zip(keys, parsed_values));
    return Ok(format);
}

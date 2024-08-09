use core::fmt;

use indexmap::IndexMap;
use std::error::Error;

use crate::vcf::{Genotype, RecordValue, VCFError, VCFHeader};

/// Represents a VCF record
///
/// This struct is used to store the values of a VCF record.
/// The info and format fields are stored as IndexMap<String, RecordValue>
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
    depth: Option<i32>,
    allele_depths: Option<Vec<i32>>,
    strand_depths: Option<(Vec<i32>, Vec<i32>)>,
}

impl fmt::Display for VariantRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dot = String::from(".");

        let pos_str: String = self.pos.to_string();
        let id_str: &str = if let Some(id) = &self.id { id } else { &dot };

        let ref_str: &str = if self.ref_bases.is_empty() {
            &dot
        } else {
            &self.ref_bases
        };

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

        let mut info_str: String = self
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
        if info_str.is_empty() {
            info_str = dot.clone();
        }

        if self.format.is_empty() {
            write!(
                f,
                "{}",
                [
                    &self.chrom,
                    &pos_str,
                    id_str,
                    ref_str,
                    &alt_str,
                    &qual_str,
                    &filter_str,
                    &info_str,
                ]
                .join("\t")
            )
        } else {
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
                    ref_str,
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

    /// Builds record from line of vcf
    /// 
    /// The header is needed to parse the info and filter values.
    /// Filter values are not checked with the header
    pub fn from_string(header: &VCFHeader, line: &str) -> Result<Self, Box<dyn Error>> {
        // split line by \t
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() < 8 {
            return Err(VCFError::InvalidRecord(format!(
                "VCF row should have at least 8 fields. Row: {line}"
            ))
            .into());
        }
        if fields.len() == 9 {
            return Err(VCFError::InvalidRecord(format!(
                "VCF row seems to have Format keys but no vlaues. Row: {line}"
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
            format: IndexMap::new(),
            genotype: None,
            depth: None,
            allele_depths: None,
            strand_depths: None,
        };

        // Add format is available
        if fields.len() > 8 {
            record.format = str_to_format(header, fields[8], fields[9])?;
        
            // find genotypes
            if let Some(RecordValue::String(genotype_str)) = record.format.get("GT") {
                record.genotype = Some(Genotype::from_string(genotype_str)?);
            } else {
                return Err(VCFError::InvalidRecord(format!(
                    "Genotype field not found in VCF row: {line}"
                ))
                .into());
            }

            if fields.len() > 10 {
                println!("More than 10 columns provided. One one sample currently supported.")
            }
        }


        // Calculate depths
        record.update_depths();

        return Ok(record);
    }

    /// Get reference to genotype
    /// 
    /// To change use [VariantRecord::set_genotype]
    pub fn genotype(&self) -> Option<&Genotype> {
        return self.genotype.as_ref();
    }

    /// Set variants genotype
    /// 
    /// This also updates the format field
    pub fn set_genotype(&mut self, genotype: Genotype) {
        self.format["GT"] = RecordValue::String(genotype.to_string());
        self.genotype = Some(genotype);
    }

    /// Returns allele with greater depth
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

    /// Returns true if variant represents a potential indel
    /// 
    /// This is not affected by genotype
    pub fn is_indel(&self) -> bool {
        // Indel if any allele is longer than 1. Even if GT = 0/0
        return self.alt.iter().any(|a| a.len() > 1) || self.ref_bases.len() > 1;
    }

    /// Returns true if variant represents a potential snp, and not an indel
    /// 
    /// This is not affected by genotype
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

    pub fn depth(&self) -> &Option<i32> {
        return &self.depth;
    }

    pub fn allele_depths(&self) -> &Option<Vec<i32>> {
        return &self.allele_depths;
    }

    pub fn strand_depths(&self) -> &Option<(Vec<i32>, Vec<i32>)> {
        return &self.strand_depths;
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
            info.insert(field.to_string(), header.parse_info_value(field, "")?);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::vcf_header::*;

    fn standard_header() -> VCFHeader {
        let mut header = VCFHeader::new();
        header.samples.push(String::from("sample"));
        header.add_header_line(HeaderLine::Filter(FilterHeader {
            id: String::from("PASS"),
            desc: String::from("All filters passed"),
        }));

        // Add DP, ADF, ADR, DP4, and MQ info headers
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("DP"),
            number: HeaderNumber::One,
            header_type: HeaderType::Integer,
            desc: String::from("Total Depth"),
        }));
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("ADF"),
            number: HeaderNumber::R,
            header_type: HeaderType::Integer,
            desc: String::from("Depth on forward strand"),
        }));
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("ADR"),
            number: HeaderNumber::R,
            header_type: HeaderType::Integer,
            desc: String::from("Depth on reverse strand"),
        }));
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("DP4"),
            number: HeaderNumber::Multiple(4),
            header_type: HeaderType::Integer,
            desc: String::from("Depth on strands for ref and alts"),
        }));
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("MQ"),
            number: HeaderNumber::One,
            header_type: HeaderType::Float,
            desc: String::from("Mapping Quality"),
        }));
        header.add_header_line(HeaderLine::Info(InfoHeader {
            id: String::from("INDEL"),
            number: HeaderNumber::Flag,
            header_type: HeaderType::Flag,
            desc: String::from("Is indel"),
        }));

        // Add GT and AD format headers
        header.add_header_line(HeaderLine::Format(FormatHeader {
            id: String::from("GT"),
            number: HeaderNumber::One,
            header_type: HeaderType::String,
            desc: String::from("Genotype"),
        }));
        header.add_header_line(HeaderLine::Format(FormatHeader {
            id: String::from("AD"),
            number: HeaderNumber::R,
            header_type: HeaderType::Integer,
            desc: String::from("Allelic depths"),
        }));

        return header;
    }

    #[test]
    fn test_from_string_to_string() {
        let std_header = standard_header();
        let record_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7".to_string();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &record_string,
        ).unwrap();

        assert_eq!(record.to_string(), record_string);

        assert_eq!(record.chrom, "ref");
        assert_eq!(record.pos, 1);
        assert_eq!(record.id, Some("id".to_string()));
        assert_eq!(record.ref_bases, "T");
        assert_eq!(record.alt, vec!["G".to_string(), "C".to_string()]);
        assert_eq!(record.qual, Some(244.589));
        assert_eq!(record.filter, vec!["F1".to_string(), "F2".to_string()]);
        assert_eq!(record.info["DP"], RecordValue::Integer(28));
        assert_eq!(record.info["ADF"], RecordValue::IntegerArray(vec![1, 2, 3]));
        assert_eq!(record.info["ADR"], RecordValue::IntegerArray(vec![2, 3, 4]));
        assert_eq!(record.info["DP4"], RecordValue::IntegerArray(vec![10, 8, 1, 5]));
        assert_eq!(record.info["MQ"], RecordValue::Float(53.0));
        assert_eq!(record.format["GT"], RecordValue::String("0/1".to_string()));
        assert_eq!(record.format["AD"], RecordValue::IntegerArray(vec![5, 6, 7]));

        assert_eq!(record.genotype().unwrap(), &Genotype{allele1: 0, allele2: 1});
        assert_eq!(record.main_allele(), 1);
        assert_eq!(record.depth().unwrap(), 28);
        assert_eq!(record.allele_depths(), &Some(vec![5, 6, 7]));
        assert_eq!(record.strand_depths(), &Some((vec![1, 2, 3], vec![2, 3, 4])));

        // If no format values, then should stop earlier
        let record_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0".to_string();
        let record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &record_string,
        ).unwrap();

        assert_eq!(record.to_string(), record_string);

    }

    #[test]
    fn test_empty_value() {
        let record = VariantRecord::empty_record();
        let record_string: String = "\t0\t.\t.\t.\t.\tPASS\t.".to_string();
        assert_eq!(record.to_string(), record_string);
    }

    #[test]
    fn test_set_genotype() {
        let std_header = standard_header();
        let record_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t0/1:5,6,7".to_string();
        let mut record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &record_string,
        ).unwrap();
        
        assert_eq!(record.genotype().unwrap(), &Genotype{allele1: 0, allele2: 1});
        
        record.set_genotype(Genotype { allele1: 2, allele2: 5 });
        assert_eq!(record.genotype().unwrap(), &Genotype{allele1: 2, allele2: 5});


        let new_record_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;MQ=53.0\tGT:AD\t2/5:5,6,7".to_string();
        assert_eq!(record.to_string(), new_record_string);
    }

    #[test]
    fn test_indel_snp_flags() {
        let std_header = standard_header();

        let snp_string: String = "ref\t1\tid\tT\tG,C\t244.589\tF1;F2\tDP=28\tGT\t0/0".to_string();
        let snp_record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &snp_string,
        ).unwrap();

        assert!(snp_record.is_snp());
        assert!(!snp_record.is_indel());

        let indel_string: String = "ref\t1\tid\tT\tGC\t244.589\tF1;F2\tDP=28\tGT\t0/0".to_string();
        let indel_record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &indel_string,
        ).unwrap();

        assert!(!indel_record.is_snp());
        assert!(indel_record.is_indel());

        let ref_string: String = "ref\t1\tid\tT\t.\t244.589\tF1;F2\tDP=28\tGT\t0/0".to_string();
        let ref_record: VariantRecord = VariantRecord::from_string(
            &std_header, 
            &ref_string,
        ).unwrap();

        assert!(!ref_record.is_snp());
        assert!(!ref_record.is_indel());
    }

    #[test]
    fn test_str_to_info() {
        let std_header = standard_header();
        let info = str_to_info(&std_header, "DP=28;ADF=1,2,3;ADR=2,3,4;DP4=10,8,1,5;INDEL;MQ=53.0").unwrap();
        assert_eq!(info["MQ"], RecordValue::Float(53.0));
        assert_eq!(info["INDEL"], RecordValue::Flag);
        assert!(info.get("T").is_none());

        assert!(str_to_info(&std_header, "T=1").is_err());
        assert!(str_to_info(&std_header, "DP=1;T=1=2").is_err());
        assert!(str_to_info(&std_header, "DP= 1").is_err());
    }

    #[test]
    fn test_str_to_format() {
        let std_header = standard_header();
        let format = str_to_format(&std_header,"GT:AD", "0/1:1,2,3").unwrap();
        assert_eq!(format["AD"], RecordValue::IntegerArray(vec![1,2,3]));
        assert_eq!(format["GT"], RecordValue::String("0/1".to_string()));

        // Missing values allowed
        assert!(str_to_format(&std_header, "GT:AD", "0/1:.").is_ok());

        // Check some errors
        assert!(str_to_format(&std_header, "GT:AD", "0/1").is_err_and(|e| e.to_string().contains("Mismatching number of keys and values")));
        assert!(str_to_format(&std_header, "GT:ADR", "0/1:1,2,3").is_err_and(|e| e.to_string().contains("No header found for key")));
    }

}

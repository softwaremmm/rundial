use crate::vcf::{RecordValue, VCFError};
use regex::Regex;
use std::fmt;
use std::{collections::HashMap, error::Error};

/// Possible Number values for INFO and FORMAT fields
/// 
/// A: One for each alternative allele
/// G: one for each genotype
/// R: One for each allele including ref
/// One: One value always
/// Flag: so 0 values
/// Multiple: Multiple values
/// Unknown: When there is a "." in the header
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum HeaderNumber {
    A,
    G,
    R,
    One,
    Flag,
    Multiple(i32),
    Unknown,
}

impl HeaderNumber {
    pub fn from_string(s: &str) -> Self {
        match s {
            "A" => HeaderNumber::A,
            "G" => HeaderNumber::G,
            "R" => HeaderNumber::R,
            "1" => HeaderNumber::One,
            "0" => HeaderNumber::Flag,
            _ => {
                if let Ok(n) = s.parse::<i32>() {
                    return HeaderNumber::Multiple(n);
                } else {
                    return HeaderNumber::Unknown;
                }
            }
        }
    }
}

impl fmt::Display for HeaderNumber {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            HeaderNumber::A => return write!(f, "A"),
            HeaderNumber::G => return write!(f, "G"),
            HeaderNumber::R => return write!(f, "R"),
            HeaderNumber::One => return write!(f, "1"),
            HeaderNumber::Flag => return write!(f, "0"),
            HeaderNumber::Multiple(n) => return write!(f, "{}", n),
            HeaderNumber::Unknown => return write!(f, "."),
        }
    }
}

/// Possible types for INFO and FORMAT fields
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum HeaderType {
    Flag,
    Integer,
    Float,
    String,
}

impl HeaderType {
    pub fn from_string(s: &str) -> Self {
        match s {
            "Flag" => HeaderType::Flag,
            "Integer" => HeaderType::Integer,
            "Float" => HeaderType::Float,
            "String" => HeaderType::String,
            _ => HeaderType::String,
        }
    }
}

impl fmt::Display for HeaderType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(
            f,
            "{}",
            match self {
                HeaderType::Flag => "Flag",
                HeaderType::Integer => "Integer",
                HeaderType::Float => "Float",
                HeaderType::String => "String",
            }
        );
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FilterHeader {
    pub id: String,
    pub desc: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct InfoHeader {
    pub id: String,
    pub number: HeaderNumber,
    pub header_type: HeaderType,
    pub desc: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FormatHeader {
    pub id: String,
    pub number: HeaderNumber,
    pub header_type: HeaderType,
    pub desc: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct MiscHeader {
    pub line: String,
}

/// Possible types of header lines in a VCF file
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum HeaderLine {
    Misc(MiscHeader),
    Info(InfoHeader),
    Format(FormatHeader),
    Filter(FilterHeader),
}

impl fmt::Display for HeaderLine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            HeaderLine::Info(h) => {
                write!(
                    f,
                    "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                    h.id, h.number, h.header_type, h.desc
                )
            }
            HeaderLine::Format(h) => {
                write!(
                    f,
                    "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
                    h.id, h.number, h.header_type, h.desc
                )
            }
            HeaderLine::Filter(h) => {
                write!(f, "##FILTER=<ID={},Description=\"{}\">", h.id, h.desc)
            }
            HeaderLine::Misc(h) => {
                write!(f, "{}", h.line)
            }
        }
    }
}

/// VCFHeader struct
///
/// Contains the header lines of a VCF file as a vector of [HeaderLine]'s.
/// The vcf specification line, and the column headers are stored separately.
#[derive(Debug, Clone)]
pub struct VCFHeader {
    pub lines: Vec<HeaderLine>,
    pub samples: Vec<String>,
    filters: HashMap<String, FilterHeader>,
    infos: HashMap<String, InfoHeader>,
    formats: HashMap<String, FormatHeader>,
}

impl fmt::Display for VCFHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut header_strings: Vec<String> = Vec::new();

        header_strings.extend(
            self.lines
                .iter()
                .map(|l| l.to_string())
                .collect::<Vec<String>>(),
        );

        if self.samples.is_empty() {
            header_strings.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string());
        } else {
            let mut columns_str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_string();
            for s in self.samples.iter() {
                columns_str.push('\t');
                columns_str.push_str(s);
            }
            header_strings.push(columns_str);
        }

        let mut header_string = header_strings.join("\n");
        header_string.push('\n');
        write!(f, "{}", header_string)
    }
}

fn get_capture<'b>(re: &Regex, s: &'b str) -> Option<&'b str> {
    if let Some(cap) = re.captures(s) {
        return Some(cap.get(1).unwrap().as_str());
    } else {
        return None;
    }
}

impl Default for VCFHeader {
    fn default() -> Self {
        Self::new()
    }
}

impl VCFHeader {
    pub fn new() -> Self {
        return VCFHeader {
            lines: Vec::new(),
            samples: Vec::new(),
            filters: HashMap::new(),
            infos: HashMap::new(),
            formats: HashMap::new(),
        };
    }

    /// Returns a map from ID to FilterHeader
    pub fn filters(&self) -> &HashMap<String, FilterHeader> {
        return &self.filters;
    }

    /// Returns a map from ID to InfoHeader
    pub fn infos(&self) -> &HashMap<String, InfoHeader> {
        return &self.infos;
    }

    /// Returns a map from ID to FormatHeader
    pub fn formats(&self) -> &HashMap<String, FormatHeader> {
        return &self.formats;
    }

    /// When changes are made to the VCFHeader, the hashmaps should be updated
    fn update_hashmaps(&mut self) {
        self.filters.clear();
        self.infos.clear();
        self.formats.clear();
        for h in self.lines.iter() {
            match h {
                HeaderLine::Filter(h) => {
                    self.filters.insert(h.id.clone(), h.clone());
                }
                HeaderLine::Info(h) => {
                    self.infos.insert(h.id.clone(), h.clone());
                }
                HeaderLine::Format(h) => {
                    self.formats.insert(h.id.clone(), h.clone());
                }
                HeaderLine::Misc(_) => {}
            }
        }
    }

    /// Creates a VCFHeader from a vector of strings
    pub fn from_lines(lines: Vec<String>) -> Self {
        let id_re = Regex::new(r"ID=([^,]+)").unwrap();
        let number_re = Regex::new(r"Number=([^,]+)").unwrap();
        let type_re = Regex::new(r"Type=([^,]+)").unwrap();
        let desc_re = Regex::new(r##"Description="(.*)""##).unwrap();

        let mut header = VCFHeader::new();

        for line in lines.into_iter() {
            if line.starts_with("##FILTER") {
                if let (Some(id), Some(desc)) =
                    (get_capture(&id_re, &line), get_capture(&desc_re, &line))
                {
                    let h = FilterHeader {
                        id: id.to_string(),
                        desc: desc.to_string(),
                    };
                    header.lines.push(HeaderLine::Filter(h));
                }
            } else if line.starts_with("##INFO") || line.starts_with("##FORMAT") {
                if let (Some(id), Some(number), Some(t), Some(desc)) = (
                    get_capture(&id_re, &line),
                    get_capture(&number_re, &line),
                    get_capture(&type_re, &line),
                    get_capture(&desc_re, &line),
                ) {
                    if line.starts_with("##INFO") {
                        let h = InfoHeader {
                            id: id.to_string(),
                            number: HeaderNumber::from_string(number),
                            header_type: HeaderType::from_string(t),
                            desc: desc.to_string(),
                        };
                        header.lines.push(HeaderLine::Info(h));
                    } else {
                        let h = FormatHeader {
                            id: id.to_string(),
                            number: HeaderNumber::from_string(number),
                            header_type: HeaderType::from_string(t),
                            desc: desc.to_string(),
                        };
                        header.lines.push(HeaderLine::Format(h));
                    }
                }
            } else if line.starts_with("#CHROM") {
                if line.contains("FORMAT") {
                    header.samples = line
                        .split('\t')
                        .skip(9)
                        .map(|s| s.to_string())
                        .collect();
                }
            } else if line.starts_with("##") {
                let h = HeaderLine::Misc(MiscHeader {
                    line: line.trim().to_string(),
                });
                header.lines.push(h);
            } else {
                panic!("Header line lacks leading ##: {line}");
            }
        }

        header.update_hashmaps();

        return header;
    }

    /// Add a new header line to the VCFHeader
    ///
    /// Will return true if an existing Filter, Info, or Format header was replaced
    pub fn add_header_line(&mut self, new_line: HeaderLine) -> bool {
        match &new_line {
            HeaderLine::Misc(h) => {
                if h.line.starts_with("##FILTER")
                    || h.line.starts_with("##INFO")
                    || h.line.starts_with("##FORMAT")
                {
                    panic!("Misc header line must not start with ##FILTER, ##INFO, ##FORMAT");
                }
                else if h.line.starts_with("#CHROM") {
                    panic!("Misc header line must not start with #CHROM. Instead set header.samples");
                }
                self.lines.push(new_line);
                return false;
            }
            HeaderLine::Filter(h) => {
                let is_present = self.filters.contains_key(&h.id);
                self.filters.insert(h.id.clone(), h.clone());
                self.lines.retain(|l| match l {
                    HeaderLine::Filter(f) => f.id != h.id,
                    _ => true,
                });
                self.lines.push(new_line);
                return is_present;
            }
            HeaderLine::Info(h) => {
                let is_present = self.infos.contains_key(&h.id);
                self.infos.insert(h.id.clone(), h.clone());
                self.lines.retain(|l| match l {
                    HeaderLine::Filter(f) => f.id != h.id,
                    _ => true,
                });
                self.lines.push(new_line);
                return is_present;
            }
            HeaderLine::Format(h) => {
                let is_present = self.formats.contains_key(&h.id);
                self.formats.insert(h.id.clone(), h.clone());
                self.lines.retain(|l| match l {
                    HeaderLine::Filter(f) => f.id != h.id,
                    _ => true,
                });
                self.lines.push(new_line);
                return is_present;
            }
        }
    }


    /// Set all the lines in the VCFHeader
    ///
    /// Will recalculate the hashmaps
    pub fn set_all_lines(&mut self, new_lines: Vec<HeaderLine>) {
        self.lines = new_lines;
        self.update_hashmaps();
    }

    /// Sorts the header lines
    ///
    /// Will keep misc headers in order, and then sort the rest (infos, formats, filters)
    pub fn sort_header(&mut self) {
        // Will keep misc headers in order (sort in stable), and then sort the rest
        let mut misc_lines: Vec<HeaderLine> = self
            .lines
            .iter()
            .filter(|l| matches!(l, HeaderLine::Misc(_))).cloned()
            .collect();
        let mut other_lines: Vec<HeaderLine> = self
            .lines
            .iter()
            .filter(|l| !matches!(l, HeaderLine::Misc(_))).cloned()
            .collect();
        other_lines.sort();
        misc_lines.extend(other_lines);
        self.lines = misc_lines;
    }

    /// Parse a value based on the header number and type
    fn parse_value(
        number: &HeaderNumber,
        value_type: &HeaderType,
        value: &str,
    ) -> Result<RecordValue, Box<dyn Error>> {
        if value == "." {
            return Ok(RecordValue::Missing);
        }

        let is_list = matches!(
            number,
            HeaderNumber::A | HeaderNumber::G | HeaderNumber::R | HeaderNumber::Multiple(_)
        ) || (matches!(number, HeaderNumber::Unknown) && value.contains(','));
        if is_list {
            match value_type {
                HeaderType::Flag => {
                    return Err(VCFError::InvalidHeader(
                        "Flag type but number suggests a list.".to_string(),
                    )
                    .into())
                }
                HeaderType::Integer => {
                    let values: Vec<_> = value
                        .split(',')
                        .map(|v| v.parse::<i32>())
                        .collect::<Result<_, _>>()?;
                    return Ok(RecordValue::IntegerArray(values));
                }
                HeaderType::Float => {
                    let values: Vec<_> = value
                        .split(',')
                        .map(|v| v.parse::<f32>())
                        .collect::<Result<_, _>>()?;
                    return Ok(RecordValue::FloatArray(values));
                }
                HeaderType::String => {
                    let values: Vec<String> = value.split(',').map(|s| s.to_string()).collect();
                    return Ok(RecordValue::StringArray(values));
                }
            }
        } else {
            match value_type {
                HeaderType::Flag => return Ok(RecordValue::Flag),
                HeaderType::Integer => return Ok(RecordValue::Integer(value.parse::<i32>()?)),
                HeaderType::Float => return Ok(RecordValue::Float(value.parse::<f32>()?)),
                HeaderType::String => return Ok(RecordValue::String(value.to_string())),
            }
        }
    }

    /// Parse an INFO value given the key
    pub fn parse_info_value(&self, key: &str, value: &str) -> Result<RecordValue, Box<dyn Error>> {
        let header_line: &InfoHeader = match self.infos.get(key) {
            None => {
                return Err(Box::new(VCFError::InvalidField(format!(
                    "No header found for key {key}. Value given: {value}"
                ))))
            }
            Some(h) => h,
        };

        Self::parse_value(&header_line.number, &header_line.header_type, value).map_err(|_| {
            VCFError::InvalidField(format!("Failed to parse {value}.\nIt has key {key}")).into()
        })
    }

    /// Parse a FORMAT value given the key
    pub fn parse_format_value(
        &self,
        key: &str,
        value: &str,
    ) -> Result<RecordValue, Box<dyn Error>> {
        let header_line: &FormatHeader = match self.formats.get(key) {
            None => {
                return Err(Box::new(VCFError::InvalidField(format!(
                    "No header found for key {key}. Value given: {value}"
                ))))
            }
            Some(h) => h,
        };

        Self::parse_value(&header_line.number, &header_line.header_type, value).map_err(|_| {
            VCFError::InvalidField(format!("Failed to parse {value}.\nIt has key {key}")).into()
        })
    }
}

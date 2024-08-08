use crate::vcf::{RecordValue, VCFError};
use regex::Regex;
use std::{collections::HashMap, error::Error};
use std::{fmt, vec};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum HeaderNumber {
    A, // One for each alternative allele
    G, // one for each genotype
    R, // One for each allele including ref
    One,
    Flag, // so 0
    Multiple(i32),
    Unknown, // Not seen this case before. Should have a "." in the header
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

#[derive(Debug, Clone)]
pub struct MiscHeader {
    pub line: String,
}
impl PartialEq for MiscHeader {
    fn eq(&self, _other: &Self) -> bool {
        return true;
    }
}
impl Eq for MiscHeader {}
impl PartialOrd for MiscHeader {
    fn partial_cmp(&self, _other: &Self) -> Option<std::cmp::Ordering> {
        return Some(std::cmp::Ordering::Equal);
    }
}
impl Ord for MiscHeader {
    fn cmp(&self, _other: &Self) -> std::cmp::Ordering {
        return std::cmp::Ordering::Equal;
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum HeaderLine {
    MiscHeader(MiscHeader),
    InfoHeader(InfoHeader),
    FormatHeader(FormatHeader),
    FilterHeader(FilterHeader),
}

impl HeaderLine {
    pub fn to_string(&self) -> String {
        match self {
            HeaderLine::InfoHeader(h) => {
                return format!(
                    "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                    h.id, h.number, h.header_type, h.desc
                );
            }
            HeaderLine::FormatHeader(h) => {
                return format!(
                    "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
                    h.id, h.number, h.header_type, h.desc
                );
            }
            HeaderLine::FilterHeader(h) => {
                return format!("##FILTER=<ID={},Description=\"{}\">", h.id, h.desc);
            }
            HeaderLine::MiscHeader(h) => {
                return format!("{}", h.line);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct VCFHeader {
    pub lines: Vec<HeaderLine>,
    pub vcf_spec: Option<HeaderLine>,
    pub column_headers: Option<HeaderLine>,
    filters: HashMap<String, FilterHeader>,
    infos: HashMap<String, InfoHeader>,
    formats: HashMap<String, FormatHeader>,
}

impl fmt::Display for VCFHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut header_strings: Vec<String> = Vec::new();

        if let Some(vcf_spec) = &self.vcf_spec {
            header_strings.push(vcf_spec.to_string());
        }

        header_strings.extend(
            self.lines
                .iter()
                .map(|l| l.to_string())
                .collect::<Vec<String>>(),
        );

        if let Some(column_headers) = &self.column_headers {
            header_strings.push(column_headers.to_string());
        }

        let mut header_string = header_strings.join("\n");
        header_string.push('\n');
        write!(f, "{}", header_string)
    }
}

fn get_capture<'a, 'b>(re: &'a Regex, s: &'b str) -> Option<&'b str> {
    if let Some(cap) = re.captures(s) {
        return Some(cap.get(1).unwrap().as_str());
    } else {
        return None;
    }
}

impl VCFHeader {
    pub fn new() -> Self {
        return VCFHeader {
            lines: Vec::new(),
            vcf_spec: None,
            column_headers: None,
            filters: HashMap::new(),
            infos: HashMap::new(),
            formats: HashMap::new(),
        };
    }

    pub fn filters(&self) -> &HashMap<String, FilterHeader> {
        return &self.filters;
    }

    pub fn infos(&self) -> &HashMap<String, InfoHeader> {
        return &self.infos;
    }

    pub fn formats(&self) -> &HashMap<String, FormatHeader> {
        return &self.formats;
    }

    fn update_hashmaps(&mut self) {
        self.filters.clear();
        self.infos.clear();
        self.formats.clear();
        for h in self.lines.iter() {
            match h {
                HeaderLine::FilterHeader(h) => {
                    self.filters.insert(h.id.clone(), h.clone());
                }
                HeaderLine::InfoHeader(h) => {
                    self.infos.insert(h.id.clone(), h.clone());
                }
                HeaderLine::FormatHeader(h) => {
                    self.formats.insert(h.id.clone(), h.clone());
                }
                HeaderLine::MiscHeader(_) => {}
            }
        }
    }

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
                    header.lines.push(HeaderLine::FilterHeader(h));
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
                            number: HeaderNumber::from_string(&number),
                            header_type: HeaderType::from_string(t),
                            desc: desc.to_string(),
                        };
                        header.lines.push(HeaderLine::InfoHeader(h));
                    } else {
                        let h = FormatHeader {
                            id: id.to_string(),
                            number: HeaderNumber::from_string(&number),
                            header_type: HeaderType::from_string(t),
                            desc: desc.to_string(),
                        };
                        header.lines.push(HeaderLine::FormatHeader(h));
                    }
                }
            } else if line.starts_with("##fileformat") {
                let h = HeaderLine::MiscHeader(MiscHeader {
                    line: line.trim().to_string(),
                });
                header.vcf_spec = Some(h);
            } else if line.starts_with("#CHROM") {
                let h = HeaderLine::MiscHeader(MiscHeader {
                    line: line.trim().to_string(),
                });
                header.column_headers = Some(h);
            } else if line.starts_with("##") {
                let h = HeaderLine::MiscHeader(MiscHeader {
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

    pub fn add_header_line(&mut self, new_line: HeaderLine) -> bool {
        // Adds new line to header, if key not already present.
        // Returns false if key is already present
        match &new_line {
            HeaderLine::FilterHeader(h) => {
                if self.filters.contains_key(&h.id) {
                    return false;
                }
            }
            HeaderLine::InfoHeader(h) => {
                if self.infos.contains_key(&h.id) {
                    return false;
                }
            }
            HeaderLine::FormatHeader(h) => {
                if self.formats.contains_key(&h.id) {
                    return false;
                }
            }
            HeaderLine::MiscHeader(h) => {
                if h.line.starts_with("##FILTER")
                    || h.line.starts_with("##INFO")
                    || h.line.starts_with("##FORMAT")
                {
                    panic!("Misc header line must not start with ##FILTER, ##INFO or ##FORMAT");
                }
                if h.line.starts_with("##fileformat") {
                    if self.vcf_spec.is_some() {
                        return false;
                    }
                    self.vcf_spec = Some(new_line);
                    return true;
                }
                if h.line.starts_with("#CHROM") {
                    if self.column_headers.is_some() {
                        return false;
                    }
                    self.column_headers = Some(new_line);
                    return true;
                }
            }
        }
        self.lines.push(new_line);
        return true;
    }

    pub fn set_all_lines(&mut self, new_lines: Vec<HeaderLine>) {
        self.lines = new_lines;
        self.update_hashmaps();
    }

    pub fn sort_header(&mut self) {
        // Will keep misc headers in order (sort in stable), and then sort the rest
        self.lines.sort();
    }

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

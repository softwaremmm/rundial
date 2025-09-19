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
            HeaderNumber::Multiple(n) => return write!(f, "{n}"),
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
impl fmt::Display for FilterHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "##FILTER=<ID={},Description=\"{}\">", self.id, self.desc);
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct InfoHeader {
    pub id: String,
    pub number: HeaderNumber,
    pub header_type: HeaderType,
    pub desc: String,
}
impl fmt::Display for InfoHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(
            f,
            "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
            self.id, self.number, self.header_type, self.desc
        );
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FormatHeader {
    pub id: String,
    pub number: HeaderNumber,
    pub header_type: HeaderType,
    pub desc: String,
}
impl fmt::Display for FormatHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(
            f,
            "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
            self.id, self.number, self.header_type, self.desc
        );
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct MiscHeader {
    pub line: String,
}
impl fmt::Display for MiscHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.line);
    }
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
                write!(f, "{h}")
            }
            HeaderLine::Format(h) => {
                write!(f, "{h}")
            }
            HeaderLine::Filter(h) => {
                write!(f, "{h}")
            }
            HeaderLine::Misc(h) => {
                write!(f, "{h}")
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
            let mut columns_str =
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_string();
            for s in self.samples.iter() {
                columns_str.push('\t');
                columns_str.push_str(s);
            }
            header_strings.push(columns_str);
        }

        let mut header_string = header_strings.join("\n");
        header_string.push('\n');
        write!(f, "{header_string}")
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

    /// Returns a VCFHeader with the standard VCFv4.2 specification
    pub fn new_std_spec() -> Self {
        return VCFHeader {
            lines: vec![HeaderLine::Misc(MiscHeader {
                line: "##fileformat=VCFv4.2".to_string(),
            })],
            samples: vec!["sample".to_string()],
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
                    header.samples = line.split('\t').skip(9).map(|s| s.to_string()).collect();
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
                } else if h.line.starts_with("#CHROM") {
                    panic!(
                        "Misc header line must not start with #CHROM. Instead set header.samples"
                    );
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
                    HeaderLine::Info(f) => f.id != h.id,
                    _ => true,
                });
                self.lines.push(new_line);
                return is_present;
            }
            HeaderLine::Format(h) => {
                let is_present = self.formats.contains_key(&h.id);
                self.formats.insert(h.id.clone(), h.clone());
                self.lines.retain(|l| match l {
                    HeaderLine::Format(f) => f.id != h.id,
                    _ => true,
                });
                self.lines.push(new_line);
                return is_present;
            }
        }
    }

    /// Add a new FILTER line to the VCFHeader
    ///
    /// Will return true if an existing FILTER header was replaced
    pub fn add_filter_line(&mut self, id: String, desc: String) -> bool {
        let h = FilterHeader { id, desc };
        return self.add_header_line(HeaderLine::Filter(h));
    }

    /// Add a new INFO line to the VCFHeader
    ///
    /// Will return true if an existing INFO header was replaced
    pub fn add_info_line(
        &mut self,
        id: String,
        number: HeaderNumber,
        header_type: HeaderType,
        desc: String,
    ) -> bool {
        let h = InfoHeader {
            id,
            number,
            header_type,
            desc,
        };
        return self.add_header_line(HeaderLine::Info(h));
    }

    /// Add a new FORMAT line to the VCFHeader
    ///
    /// Will return true if an existing FORMAT header was replaced
    pub fn add_format_line(
        &mut self,
        id: String,
        number: HeaderNumber,
        header_type: HeaderType,
        desc: String,
    ) -> bool {
        let h = FormatHeader {
            id,
            number,
            header_type,
            desc,
        };
        return self.add_header_line(HeaderLine::Format(h));
    }

    /// Add a new Misc line to the VCFHeader
    pub fn add_misc_line(&mut self, line: String) {
        let h = MiscHeader { line };
        self.add_header_line(HeaderLine::Misc(h));
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
    pub fn sort(&mut self) {
        // Will keep misc headers in order (sort in stable), and then sort the rest
        let mut misc_lines: Vec<HeaderLine> = self
            .lines
            .iter()
            .filter(|l| matches!(l, HeaderLine::Misc(_)))
            .cloned()
            .collect();
        let mut other_lines: Vec<HeaderLine> = self
            .lines
            .iter()
            .filter(|l| !matches!(l, HeaderLine::Misc(_)))
            .cloned()
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
            match (value_type, value) {
                (HeaderType::Flag, "") => return Ok(RecordValue::Flag),
                (HeaderType::Flag, _) => {
                    return Err(VCFError::InvalidHeader(
                        "Flag type but value is not empty.".to_string(),
                    )
                    .into())
                }
                (HeaderType::Integer, _) => return Ok(RecordValue::Integer(value.parse::<i32>()?)),
                (HeaderType::Float, _) => return Ok(RecordValue::Float(value.parse::<f32>()?)),
                (HeaderType::String, _) => return Ok(RecordValue::String(value.to_string())),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_number_from_string() {
        assert_eq!(HeaderNumber::from_string("A"), HeaderNumber::A);
        assert_eq!(HeaderNumber::from_string("G"), HeaderNumber::G);
        assert_eq!(HeaderNumber::from_string("R"), HeaderNumber::R);
        assert_eq!(HeaderNumber::from_string("1"), HeaderNumber::One);
        assert_eq!(HeaderNumber::from_string("0"), HeaderNumber::Flag);
        assert_eq!(HeaderNumber::from_string("42"), HeaderNumber::Multiple(42));
        assert_eq!(HeaderNumber::from_string("42.0"), HeaderNumber::Unknown);
        assert_eq!(HeaderNumber::from_string("K"), HeaderNumber::Unknown);
    }

    #[test]
    fn test_header_number_display() {
        assert_eq!(HeaderNumber::A.to_string(), "A");
        assert_eq!(HeaderNumber::G.to_string(), "G");
        assert_eq!(HeaderNumber::R.to_string(), "R");
        assert_eq!(HeaderNumber::One.to_string(), "1");
        assert_eq!(HeaderNumber::Flag.to_string(), "0");
        assert_eq!(HeaderNumber::Multiple(42).to_string(), "42");
        assert_eq!(HeaderNumber::Unknown.to_string(), ".");
    }

    #[test]
    fn test_header_type_from_string() {
        assert_eq!(HeaderType::from_string("Flag"), HeaderType::Flag);
        assert_eq!(HeaderType::from_string("Integer"), HeaderType::Integer);
        assert_eq!(HeaderType::from_string("Float"), HeaderType::Float);
        assert_eq!(HeaderType::from_string("String"), HeaderType::String);
        assert_eq!(HeaderType::from_string("42.0"), HeaderType::String);
        assert_eq!(HeaderType::from_string("Other"), HeaderType::String);
    }

    #[test]
    fn test_header_type_display() {
        assert_eq!(HeaderType::Flag.to_string(), "Flag");
        assert_eq!(HeaderType::Integer.to_string(), "Integer");
        assert_eq!(HeaderType::Float.to_string(), "Float");
        assert_eq!(HeaderType::String.to_string(), "String");
    }

    fn example_filter_header() -> FilterHeader {
        FilterHeader {
            id: "PASS".to_string(),
            desc: "All filters passed".to_string(),
        }
    }

    fn example_info_header() -> InfoHeader {
        InfoHeader {
            id: "DP".to_string(),
            number: HeaderNumber::One,
            header_type: HeaderType::Integer,
            desc: "Total Depth".to_string(),
        }
    }

    fn example_format_header() -> FormatHeader {
        FormatHeader {
            id: "GT".to_string(),
            number: HeaderNumber::One,
            header_type: HeaderType::String,
            desc: "Genotype".to_string(),
        }
    }

    fn example_misc_header() -> MiscHeader {
        MiscHeader {
            line: "##fileformat=VCFv4.2".to_string(),
        }
    }

    fn example_column_header() -> String {
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample".to_string()
    }

    #[test]
    fn test_header_display() {
        assert_eq!(
            example_filter_header().to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        );

        assert_eq!(
            example_info_header().to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
        );

        assert_eq!(
            example_format_header().to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        );

        assert_eq!(example_misc_header().to_string(), "##fileformat=VCFv4.2");
    }

    #[test]
    fn test_header_from_lines() {
        let lines = vec![
            example_misc_header().to_string(),
            example_filter_header().to_string(),
            example_info_header().to_string(),
            example_format_header().to_string(),
            example_column_header(),
        ];

        let header = VCFHeader::from_lines(lines.clone());

        assert_eq!(header.to_string(), lines.join("\n") + "\n");

        assert_eq!(header.lines.len(), 4);
        assert!(header.filters.contains_key("PASS"));
        assert!(header.infos.contains_key("DP"));
        assert!(header.formats.contains_key("GT"));
        assert_eq!(header.samples, vec!["sample".to_string()]);
    }

    #[test]
    fn test_set_all_lines() {
        let mut header = VCFHeader::new();
        header.set_all_lines(
            vec![
                HeaderLine::Misc(example_misc_header()),
                HeaderLine::Filter(example_filter_header()),
                HeaderLine::Info(example_info_header()),
                HeaderLine::Format(example_format_header()),
            ]
            .clone(),
        );
        assert_eq!(header.lines.len(), 4);
        assert!(header.filters.contains_key("PASS"));
        assert!(header.infos.contains_key("DP"));
    }

    #[test]
    fn test_sort() {
        let lines = vec![
            example_format_header().to_string(),
            example_filter_header().to_string(),
            example_info_header().to_string(),
            example_column_header(),
            example_misc_header().to_string(),
        ];

        let mut header = VCFHeader::from_lines(lines.clone());
        header.sort();

        let sorted_lines = [
            example_misc_header().to_string(),
            example_info_header().to_string(),
            example_format_header().to_string(),
            example_filter_header().to_string(),
            example_column_header(),
        ];
        assert_eq!(header.to_string(), sorted_lines.join("\n") + "\n");
    }

    #[test]
    fn test_add_header_line() {
        let mut header = VCFHeader::new();

        assert!(!header.add_header_line(HeaderLine::Misc(example_misc_header())));
        assert!(!header.add_header_line(HeaderLine::Misc(example_misc_header())));
        assert!(!header.add_header_line(HeaderLine::Misc(example_misc_header())));

        assert!(!header.add_header_line(HeaderLine::Filter(example_filter_header())));
        assert!(header.add_header_line(HeaderLine::Filter(example_filter_header())));

        assert!(!header.add_header_line(HeaderLine::Info(example_info_header())));
        assert!(header.add_header_line(HeaderLine::Info(example_info_header())));

        assert!(!header.add_header_line(HeaderLine::Format(example_format_header())));
        assert!(header.add_header_line(HeaderLine::Format(example_format_header())));

        let new_info = InfoHeader {
            id: "INDEL".to_string(),
            number: HeaderNumber::Flag,
            header_type: HeaderType::Flag,
            desc: "Is indel".to_string(),
        };

        let replace_info = InfoHeader {
            id: "INDEL".to_string(),
            number: HeaderNumber::A,
            header_type: HeaderType::Float,
            desc: "different".to_string(),
        };

        assert!(!header.add_header_line(HeaderLine::Info(new_info.clone())));
        assert!(header.add_header_line(HeaderLine::Info(replace_info.clone())));
        assert!(header.infos.get("INDEL").unwrap() == &replace_info);
    }

    #[test]
    fn test_parse_value() {
        // test individuals
        assert_eq!(
            VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Flag, "").unwrap(),
            RecordValue::Flag
        );
        assert_eq!(
            VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Integer, ".").unwrap(),
            RecordValue::Missing
        );
        assert_eq!(
            VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Integer, "42").unwrap(),
            RecordValue::Integer(42)
        );
        assert_eq!(
            VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Float, "42.0").unwrap(),
            RecordValue::Float(42.0)
        );
        assert_eq!(
            VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::String, "hello").unwrap(),
            RecordValue::String("hello".to_string())
        );
        // test lists
        for number in [
            HeaderNumber::A,
            HeaderNumber::G,
            HeaderNumber::R,
            HeaderNumber::Multiple(3),
        ] {
            assert_eq!(
                VCFHeader::parse_value(&number, &HeaderType::Integer, "1,2,3").unwrap(),
                RecordValue::IntegerArray(vec![1, 2, 3])
            );
            assert_eq!(
                VCFHeader::parse_value(&number, &HeaderType::Float, "1.0,2.0,3.0").unwrap(),
                RecordValue::FloatArray(vec![1.0, 2.0, 3.0])
            );
            assert_eq!(
                VCFHeader::parse_value(&number, &HeaderType::String, "hello,world").unwrap(),
                RecordValue::StringArray(vec!["hello".to_string(), "world".to_string()])
            );
            // Test that a flag type with multiple values is an error
            assert!(VCFHeader::parse_value(&number, &HeaderType::Flag, "1,2,3").is_err());
        }

        // Check for some errors if type does not match value
        assert!(VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Flag, "42").is_err());
        assert!(VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Integer, "42.0").is_err());
        assert!(VCFHeader::parse_value(&HeaderNumber::One, &HeaderType::Float, "hello").is_err());
    }

    #[test]
    fn test_parse_filter_and_info() {
        let mut header = VCFHeader::new();
        header.add_header_line(HeaderLine::Info(example_info_header()));
        header.add_header_line(HeaderLine::Format(example_format_header()));

        assert_eq!(
            header.parse_info_value("DP", "42").unwrap(),
            RecordValue::Integer(42)
        );
        assert_eq!(
            header.parse_info_value("DP", ".").unwrap(),
            RecordValue::Missing
        );
        assert!(header.parse_info_value("DP", "42.0").is_err());

        assert_eq!(
            header.parse_format_value("GT", "0/1").unwrap(),
            RecordValue::String("0/1".to_string())
        );
        assert_eq!(
            header.parse_format_value("GT", ".").unwrap(),
            RecordValue::Missing
        );

        assert!(header
            .parse_info_value("MISSING", "42.0")
            .is_err_and(|e| e.to_string().contains("No header found for key MISSING")));
    }
}

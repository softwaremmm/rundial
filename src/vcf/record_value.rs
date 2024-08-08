use core::fmt;

/// Represents a value for INFO or FORMAT in a VCF record
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

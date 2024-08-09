use core::fmt;

/// Represents a value for INFO or FORMAT in a VCF record
#[derive(Debug, Clone, PartialEq)]
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
                    .map(|f| format!("{:?}", f))
                    .collect::<Vec<String>>()
                    .join(","),
                RecordValue::String(s) => s.to_string(),
                RecordValue::StringArray(arr) => arr.join(","),
                RecordValue::Missing => String::from("."),
            }
        )
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_value_display() {
        let int = RecordValue::Integer(42);
        assert_eq!(int.to_string(), "42");

        let fl = RecordValue::Float(42.0);
        assert_eq!(fl.to_string(), "42.0");

        let flag = RecordValue::Flag;
        assert_eq!(flag.to_string(), "");

        let int_arr = RecordValue::IntegerArray(vec![1, 2, 3]);
        assert_eq!(int_arr.to_string(), "1,2,3");

        let fl_arr = RecordValue::FloatArray(vec![1.0, 2.0, 3.0]);
        assert_eq!(fl_arr.to_string(), "1.0,2.0,3.0");

        let s = RecordValue::String(String::from("hello"));
        assert_eq!(s.to_string(), "hello");

        let s_arr = RecordValue::StringArray(vec![String::from("hello"), String::from("world")]);
        assert_eq!(s_arr.to_string(), "hello,world");

        let missing = RecordValue::Missing;
        assert_eq!(missing.to_string(), ".");
    }
}
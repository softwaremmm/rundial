use core::fmt;

/// Represents a value for INFO or FORMAT in a VCF record
#[derive(Debug, Clone)]
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

impl PartialEq for RecordValue {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (RecordValue::Flag, RecordValue::Flag) => true,
            (RecordValue::Integer(i1), RecordValue::Integer(i2)) => i1 == i2,
            (RecordValue::Float(f1), RecordValue::Float(f2)) => {
                if f1.is_nan() && f2.is_nan() {
                    true
                } else {
                    f1 == f2
                }
            }
            (RecordValue::String(s1), RecordValue::String(s2)) => s1 == s2,
            (RecordValue::IntegerArray(arr1), RecordValue::IntegerArray(arr2)) => arr1 == arr2,
            (RecordValue::FloatArray(arr1), RecordValue::FloatArray(arr2)) => {
                if arr1.len() != arr2.len() {
                    return false;
                }
                for (f1, f2) in arr1.iter().zip(arr2.iter()) {
                    if f1.is_nan() && f2.is_nan() {
                        continue;
                    }
                    if f1 != f2 {
                        return false;
                    }
                }
                return true;
            }
            (RecordValue::StringArray(arr1), RecordValue::StringArray(arr2)) => arr1 == arr2,
            (RecordValue::Missing, RecordValue::Missing) => true,
            _ => false,
        }
    }
}
impl Eq for RecordValue {}

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

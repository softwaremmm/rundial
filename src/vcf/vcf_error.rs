use core::fmt;
use std::error::Error;

/// Common Error types for VCF parsing
#[derive(Debug, Clone)]
pub enum VCFError {
    InvalidRecord(String),
    InvalidField(String),
    InvalidHeader(String),
}
impl fmt::Display for VCFError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            VCFError::InvalidRecord(s) => write!(f, "Invalid VCF Record, {}", s),
            VCFError::InvalidField(s) => write!(f, "Invalid VCF Field, {}", s),
            VCFError::InvalidHeader(s) => write!(f, "Invalid VCF Header. {}", s),
        }
    }
}
impl Error for VCFError {}
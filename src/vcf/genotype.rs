use core::fmt;
use std::error::Error;
use crate::vcf::VCFError;

/// Represents a genotype in a VCF record
/// 
/// assumes haploid genotype, and does not support phased genotypes
#[derive(Debug, Clone, PartialEq, Eq)]
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotype_display() {
        let gt = Genotype { allele1: 0, allele2: 1 };
        assert_eq!(gt.to_string(), "0/1");

        let gt = Genotype { allele1: 1, allele2: 2 };
        assert_eq!(gt.to_string(), "1/2");

        let gt = Genotype { allele1: 0, allele2: 0 };
        assert_eq!(gt.to_string(), "0/0");

        let gt = Genotype { allele1: -1, allele2: 0 };
        assert_eq!(gt.to_string(), "./0");

        let gt = Genotype { allele1: 0, allele2: -1 };
        assert_eq!(gt.to_string(), "0/.");
    }

    #[test]
    fn test_genotype_from_string() {
        let gt = Genotype::from_string("0/1").unwrap();
        assert_eq!(gt.allele1, 0);
        assert_eq!(gt.allele2, 1);

        let gt = Genotype::from_string("1/1").unwrap();
        assert_eq!(gt.allele1, 1);
        assert_eq!(gt.allele2, 1);

        let gt = Genotype::from_string("0/0").unwrap();
        assert_eq!(gt.allele1, 0);
        assert_eq!(gt.allele2, 0);

        let gt = Genotype::from_string("./0").unwrap();
        assert_eq!(gt.allele1, -1);
        assert_eq!(gt.allele2, 0);

        let gt = Genotype::from_string("./.").unwrap();
        assert_eq!(gt.allele1, -1);
        assert_eq!(gt.allele2, -1);

        let gt = Genotype::new();
        assert_eq!(gt.allele1, -1);
        assert_eq!(gt.allele2, -1);

        let gt = Genotype::from_string("0/.").unwrap();
        assert_eq!(gt.allele1, 0);
        assert_eq!(gt.allele2, -1);

        let gt = Genotype::from_string("0/1/2");
        assert!(gt.is_err());

        let gt = Genotype::from_string("0/a");
        assert!(gt.is_err());

        let gt = Genotype::from_string("-2/0");
        assert!(gt.is_err());
    }

    #[test]
    fn test_flag_functions() {
        let gt = Genotype { allele1: 0, allele2: 1 };
        assert!(gt.is_het());
        assert!(!gt.is_hom());
        assert!(!gt.is_hom_ref());

        let gt = Genotype { allele1: 0, allele2: 0 };
        assert!(!gt.is_het());
        assert!(gt.is_hom());
        assert!(gt.is_hom_ref());

        let gt = Genotype { allele1: 1, allele2: 1 };
        assert!(!gt.is_het());
        assert!(gt.is_hom());
        assert!(!gt.is_hom_ref());

        let gt = Genotype { allele1: -1, allele2: 0 };
        assert!(!gt.is_het());
        assert!(!gt.is_hom());
        assert!(!gt.is_hom_ref());
    }
}
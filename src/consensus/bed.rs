//! Bed file format
//! 
//! Used for specifying genomes regions, particularly for masking

use core::panic;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::fs::File;

/// Bed file format
/// 
/// Browser Extensible Data (BED) format is a simple format for specifying genomic regions.
/// The start column is inclusive and the end column is exclusive.
/// Bed files are 0-based normally, meaning the first base is 0. UNLIKE VCF FILES!!
#[derive(Debug, Clone)]
pub struct Bed {
    regions: HashMap<String, HashSet<u32>>,
}

impl Bed {
    /// Make a new Bed struct from a file
    pub fn from_file(file: &str) -> Self {
        let mut regions = HashMap::new();

        let reader = BufReader::new(File::open(file).unwrap());
        for line in reader.lines() {
            let line: String = line.unwrap();
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 {
                panic!("Bed file must have at least 3 columns");
            }
            let chrom = parts[0].to_string();
            let start = parts[1].parse::<u32>().unwrap();
            let end = parts[2].parse::<u32>().unwrap();

            let region = regions.entry(chrom).or_insert(HashSet::new());
            for i in start..end {
                region.insert(i);
            }
        }
        Bed {
            regions,
        }
    }

    /// Check if a position is in the bed file
    pub fn contains(&self, chrom: &str, pos: &u32) -> bool {
        if let Some(region) = self.regions.get(chrom) {
            if region.contains(pos) {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed() {
        let bed = Bed::from_file("test_data/mask.bed");
        assert!(bed.contains("NC_000962.3", &0));
        assert!(bed.contains("NC_000962.3", &9));
        assert!(!bed.contains("NC_000962.3", &10));
        assert!(!bed.contains("NC_000962.3", &19));
        assert!(bed.contains("NC_000962.3", &20));
        assert!(bed.contains("other_chrom", &19));
        assert!(!bed.contains("missing", &19));
    }
}
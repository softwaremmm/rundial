use crate::vcf::VariantRecord;

use super::{Bed, ConsensusParams, HetOption};
use super::{FILTERED, HET, MASKED, NULL};

pub fn repeat_char(c: char, n: usize) -> String {
    std::iter::repeat(c).take(n).collect()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Change {
    Null,
    Ref,
    Snp,
    Del,
    Ins,
    ComplexIndel,
}
impl Change {
    pub fn from_ref_alt(ref_bases: &str, alt_bases: &str) -> Self {
        if alt_bases
            .chars()
            .any(|c| [NULL, FILTERED, HET, MASKED].contains(&c))
        {
            return Change::Null;
        }
        if ref_bases == alt_bases {
            return Change::Ref;
        }
        if ref_bases.len() == alt_bases.len() {
            if ref_bases.len() == 1 {
                return Change::Snp;
            } else {
                return Change::ComplexIndel;
            }
        } else if ref_bases.len() > alt_bases.len() && ref_bases.starts_with(alt_bases) {
            return Change::Del;
        } else if ref_bases.len() < alt_bases.len() && alt_bases.starts_with(ref_bases) {
            return Change::Ins;
        } else {
            return Change::ComplexIndel;
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Classification {
    pub pos: usize, // 0-based
    pub ref_bases: String,
    pub new_bases: String,
    pub change: Change,
    pub is_het: bool,
    pub has_minor_population: bool,
    pub is_filtered: bool,
    pub has_indel_form: bool, // Either ref or alt has multiple bases, gets processed later
}
impl Default for Classification {
    fn default() -> Self {
        Classification {
            pos: 0,
            ref_bases: "".to_string(),
            new_bases: "".to_string(),
            change: Change::Null,
            is_het: false,
            has_minor_population: false,
            is_filtered: false,
            has_indel_form: false,
        }
    }
}

pub struct Classifier {
    pub params: ConsensusParams,
    pub mask: Option<Bed>,
}

impl Classifier {
    pub fn new(params: &ConsensusParams) -> Self {
        let mask: Option<Bed> = params.mask.as_ref().map(|s| Bed::from_file(s));

        Classifier {
            params: params.clone(),
            mask,
        }
    }

    pub fn get_flags(&self, record: &VariantRecord) -> Vec<String> {
        if !self.params.use_filters {
            return Vec::new();
        }

        let mut flags: Vec<String> = record
            .filter
            .clone()
            .into_iter()
            .filter(|f| {
                if f == "PASS" || f == "." || f == "RefCall" {
                    return false;
                }
                return true;
            })
            .collect();

        if let Some(ignore_list) = &self.params.filter_ignore_list {
            flags.retain(|f| !ignore_list.contains(f));
        }

        return flags;
    }

    pub fn is_masked(&self, record: &VariantRecord) -> bool {
        if let Some(mask) = &self.mask {
            return mask.contains(&record.chrom, &(record.pos - 1));
        }
        return false;
    }

    /// Determines how to apply a variant based on parameters set
    ///
    /// Only considers the variant in isolation, so will not check for clashes etc.
    /// A Classification will contain the new bases to apply, the change type,
    /// and some flags to aid in the consensus process
    pub fn classify(&self, record: &mut VariantRecord) -> Classification {
        let mut classification = Classification::default();
        classification.pos = record.pos_idx();
        classification.ref_bases = record.ref_bases.clone();
        // True if any of the alleles are multiple bases
        classification.has_indel_form = record.is_indel();

        // updates the flags in the record
        record.filter = self.get_flags(record);
        classification.is_filtered = !record.filter.is_empty();

        let gt = record.genotype().expect("Genotype not found");
        classification.is_het = gt.is_het();

        if let (Some(allelic_depths), Some(minor_threshold)) =
            (record.allele_depths(), self.params.minor_pop_threshold)
        {
            // look for alleles not in GT which have depth > minor_pop_threshold
            classification.has_minor_population =
                allelic_depths.iter().enumerate().any(|(i, depth)| {
                    let i = i as i32;
                    if i == gt.allele1 || i == gt.allele2 {
                        return false;
                    }
                    return *depth >= minor_threshold;
                });
        }

        // If filtered, just mask the site
        if classification.is_filtered {
            classification.change = Change::Null;
            classification.new_bases = repeat_char(FILTERED, record.ref_bases.len());
            return classification;
        }

        // Case match based on genotype.
        // Hets for indels vs snps may be handled different based on parameters
        match (gt.allele1, gt.allele2) {
            (-1, -1) => {
                classification.change = Change::Null;
                classification.new_bases = repeat_char(NULL, record.ref_bases.len());
            }
            (-1, _) | (_, -1) => {
                panic!("Does not support partial null GT like ./1")
            }
            (0, 0) => {
                classification.change = Change::Ref;
                classification.new_bases = record.ref_bases.clone();
            }
            (i, j) if i == j => {
                classification.new_bases = record.alt[(i - 1) as usize].clone();
                classification.change =
                    Change::from_ref_alt(&classification.ref_bases, &classification.new_bases);
            }
            (i, j) => {
                // Difficult het case
                classification.is_het = true;
                let het_option = if classification.has_indel_form {
                    &self.params.het_indel_option
                } else {
                    &self.params.het_snp_option
                };

                match het_option {
                    HetOption::Mask => {
                        classification.change = Change::Null;
                        classification.new_bases = repeat_char(HET, record.ref_bases.len());
                    }
                    HetOption::Ref if i == 0 || j == 0 => {
                        classification.change = Change::Ref;
                        classification.new_bases = record.ref_bases.clone();
                    }
                    HetOption::Alt if i == 0 || j == 0 => {
                        let alt_allele = std::cmp::max(i, j);
                        classification.new_bases = record.alt[(alt_allele - 1) as usize].clone();
                        classification.change = Change::from_ref_alt(
                            &classification.ref_bases,
                            &classification.new_bases,
                        );
                    }
                    _ => {
                        if let Some(main_allele) = record.main_allele() {
                            if main_allele == 0 {
                                classification.new_bases = record.ref_bases.clone();
                            } else {
                                classification.new_bases =
                                    record.alt[(main_allele - 1) as usize].clone();
                            }
                            classification.change = Change::from_ref_alt(
                                &classification.ref_bases,
                                &classification.new_bases,
                            );
                        } else {
                            println!("Main allele not found, masking het site");
                            classification.change = Change::Null;
                            classification.new_bases = repeat_char(HET, record.ref_bases.len());
                        }
                    }
                }
            }
        }

        classification
    }
}

#[cfg(test)]
mod tests;

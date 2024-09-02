use std::iter;

use crate::vcf::VariantRecord;

use super::{Bed, ConsensusParams, HetOption};
use super::{FILTERED, HET, MASKED, NULL};

pub fn repeat_char(c: char, n: usize) -> String {
    std::iter::repeat(c).take(n).collect()
}

/// Simplifies ref-alt pair by removing matching trailing and leading bases.
/// 1. Removes shared first base if the same
/// 2. Removes all shared trailing bases
/// 3. Removes all shared leading bases
///
/// Returns the number of leading bases removed, the new ref and alt strings
fn simplify_ref_alt(ref_bases: &str, alt_bases: &str) -> (usize, String, String) {
    let mut new_ref = "".to_string();
    let mut new_alt = "".to_string();
    let mut ref_iter = ref_bases.chars();
    let mut alt_iter = alt_bases.chars();
    let mut final_char_ref: Option<char> = None;
    let mut final_char_alt: Option<char> = None;
    let mut counter = 0;

    // Remove shared first base as a special check first
    if let (Some(r), Some(c)) = (ref_bases.chars().next(), alt_bases.chars().next()) {
        if r == c {
            counter += 1;
            ref_iter.next();
            alt_iter.next();
        }
    } else {
        // If one is empty then no change to make
        return (0, ref_bases.to_string(), alt_bases.to_string());
    }

    // Remove matching trailing bases
    let mut ref_iter = ref_iter.rev();
    let mut alt_iter = alt_iter.rev();
    for r in ref_iter.by_ref() {
        if let Some(c) = alt_iter.next() {
            if r == c {
                continue;
            }
            final_char_ref = Some(r);
            final_char_alt = Some(c);
            break;
        }
        final_char_ref = Some(r);
        break;
    }

    // Reverse iterators to original order and add final chars
    let mut ref_iter: Box<dyn Iterator<Item = char>> = if let Some(c) = final_char_ref {
        Box::new(ref_iter.rev().chain(iter::once(c)))
    } else {
        Box::new(ref_iter.rev())
    };

    let mut alt_iter: Box<dyn Iterator<Item = char>> = if let Some(c) = final_char_alt {
        Box::new(alt_iter.rev().chain(iter::once(c)))
    } else {
        Box::new(alt_iter.rev())
    };

    // Remove matching leading bases
    for r in ref_iter.by_ref() {
        if let Some(c) = alt_iter.next() {
            if r == c {
                counter += 1;
                continue;
            }
            new_ref.push(r);
            new_alt.push(c);
            break;
        }
        new_ref.push(r);
        break;
    }
    new_ref.extend(ref_iter);
    new_alt.extend(alt_iter);

    (counter, new_ref, new_alt)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Change {
    Null,
    HetMask,
    Ref,
    Snp,
    Del,
    Ins,
    Mnp,        // multiple nucleotide polymorphism
    ComplexDel, // change which decreases length but is not simple deletion
    ComplexIns, // change which increases length but is not simple insertion
}
impl Change {
    pub fn from_ref_alt(ref_bases: &str, alt_bases: &str) -> Self {
        if alt_bases.chars().any(|c| c == HET) {
            return Change::HetMask;
        }
        if alt_bases
            .chars()
            .any(|c| [NULL, FILTERED, MASKED].contains(&c))
        {
            return Change::Null;
        }
        if ref_bases == alt_bases {
            return Change::Ref;
        }
        match ref_bases.len().cmp(&alt_bases.len()) {
            std::cmp::Ordering::Less => {
                if alt_bases.starts_with(ref_bases) {
                    return Change::Ins;
                }
                return Change::ComplexIns;
            }
            std::cmp::Ordering::Greater => {
                if ref_bases.starts_with(alt_bases) {
                    return Change::Del;
                }
                return Change::ComplexDel;
            }
            std::cmp::Ordering::Equal => {
                if ref_bases.len() == 1 {
                    return Change::Snp;
                } else {
                    return Change::Mnp;
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Classification {
    pub pos: usize, // 0-based
    pub change: Change,
    pub ref_bases: String,
    pub new_bases: String,
    pub is_het: bool,
    pub has_minor_population: bool,
    pub is_filtered: bool,
    pub has_indel_alleles: bool, // Either ref or alt has multiple bases, gets processed later
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
            has_indel_alleles: false,
        }
    }
}
impl Classification {
    pub fn is_simple_ref(&self) -> bool {
        self.change == Change::Ref && !self.is_het && !self.has_minor_population
    }

    fn set_ref_alt_and_simplify(&mut self, ref_bases: &str, alt_bases: &str) {
        let (leading_bases, ref_bases, alt_bases) = simplify_ref_alt(ref_bases, alt_bases);
        self.ref_bases = ref_bases;
        self.new_bases = alt_bases;
        self.pos += leading_bases;
        self.change = Change::from_ref_alt(&self.ref_bases, &self.new_bases);
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
        // updates the flags in the record
        record.filter = self.get_flags(record);

        let mut classification = Classification {
            pos: record.pos_idx(),
            ref_bases: record.ref_bases.clone(),
            has_indel_alleles: record.is_indel(), // True if any of the alleles are multiple bases
            is_filtered: !record.filter.is_empty(),
            ..Default::default()
        };

        let gt = record.genotype().expect("Genotype not found");
        classification.is_het = gt.is_het();

        if let (Some(allelic_depths), Some(minor_threshold)) =
            (record.allele_depths(), self.params.minor_pop_threshold)
        {
            // look for alleles not in GT which have depth >= minor_pop_threshold
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
                classification
                    .set_ref_alt_and_simplify(&record.ref_bases, &record.alt[(i - 1) as usize]);
            }
            (i, j) => {
                // Difficult het case
                classification.is_het = true;
                let het_option = if record.is_indel() {
                    &self.params.het_indel_option
                } else {
                    &self.params.het_snp_option
                };

                match het_option {
                    HetOption::Mask => {
                        classification.change = Change::HetMask;
                        classification.new_bases = repeat_char(HET, record.ref_bases.len());
                    }
                    HetOption::Ref if i == 0 || j == 0 => {
                        classification.change = Change::Ref;
                        classification.new_bases = record.ref_bases.clone();
                    }
                    HetOption::Alt if i == 0 || j == 0 => {
                        let alt_allele = std::cmp::max(i, j);
                        classification.set_ref_alt_and_simplify(
                            &record.ref_bases,
                            &record.alt[(alt_allele - 1) as usize],
                        );
                    }
                    _ => {
                        if let Some(main_allele) = record.main_allele() {
                            if main_allele == 0 {
                                classification.new_bases = record.ref_bases.clone();
                            } else {
                                classification.set_ref_alt_and_simplify(
                                    &record.ref_bases,
                                    &record.alt[(main_allele - 1) as usize],
                                );
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

use std::iter;

use crate::vcf::VariantRecord;

use super::MIXED;
use super::{Bed, ConsensusParams, HetOption};
use super::{FILTERED, HET, MASKED, NULL};

pub fn repeat_char(c: char, n: usize) -> String {
    std::iter::repeat_n(c, n).collect()
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

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Change {
    #[default]
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
    pub is_het: bool, // Check if GT is heterozygous OR has minor population and only filter is MIN_FRS
    pub has_minor_population: bool, // Alleles which are not in GT but have depth >= minor_pop_threshold
    pub is_filtered: bool,          // Has any (non-ignored) filter. MIN_FRS alone is considered het
    pub has_indel_alleles: bool,    // Either ref or alt has multiple bases, gets processed later
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
    pub minor_pop_threshold: Option<i32>,
}

impl Classifier {
    pub fn new(params: &ConsensusParams, is_support: bool) -> Self {
        let mask: Option<Bed> = params.mask.as_ref().map(|s| Bed::from_file(s));
        let minor_pop_threshold = if is_support {
            params.support_minor_pop_threshold
        } else {
            params.minor_pop_threshold
        };

        Classifier {
            params: params.clone(),
            mask,
            minor_pop_threshold,
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

    /// Returns a list of overriding flags for the record
    /// These are flags in the support vcf which should be applied to the consensus
    /// even if the main vcf lacks them
    pub fn get_overriding_flags(&self, record: &VariantRecord) -> Option<Vec<String>> {
        if let Some(overriding_flags) = &self.params.overriding_filters {
            let flags: Vec<String> = self
                .get_flags(record)
                .into_iter()
                .filter(|filter| overriding_flags.contains(filter))
                .collect();
            if flags.is_empty() {
                return None;
            }

            return Some(flags);
        }

        return None;
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

        let gt = record.genotype().expect("Genotype not found").clone();

        if let Some(minor_threshold) = self.minor_pop_threshold {
            if let Some(allelic_depths) = record.allele_depths() {
                for (i, depth) in allelic_depths.iter().enumerate() {
                    let i = i as i32;
                    if i == gt.allele1 || i == gt.allele2 {
                        continue; // Skip alleles in GT
                    }
                    if *depth >= minor_threshold {
                        classification.has_minor_population = true;
                        break;
                    }
                }
            }
        }

        classification.is_het = gt.is_het();

        // special check for het via minor population
        if classification.has_minor_population
            && record.filter.contains(&"MIN_FRS".to_string())
            && record.filter.len() == 1
        {
            classification.is_het = true;
            classification.is_filtered = false; // treat as het rather than filtered
        }

        if classification.is_het {
            // Then add this to the info field of the record
            record
                .info
                .insert(MIXED.to_string(), crate::vcf::RecordValue::Flag);
        }

        /// Indel is standard if ref and all alts start with the same base
        /// idea being that first base is not really part of the change
        fn is_standard_indel(ref_bases: &str, alt_bases: &[String]) -> bool {
            if let Some(first_base) = ref_bases.chars().next() {
                if alt_bases.iter().all(|alt| alt.starts_with(first_base)) {
                    return true;
                }
            }
            return false;
        }

        // If filtered, just mask the site
        if classification.is_filtered {
            classification.change = Change::Null;
            classification.new_bases = repeat_char(FILTERED, record.ref_bases.len());
            return classification;
        }

        // Case match based on genotype.
        // Hets for indels vs snps may be handled different based on parameters
        match (classification.is_het, gt.allele1, gt.allele2) {
            (_, -1, -1) => {
                classification.change = Change::Null;
                classification.new_bases = repeat_char(NULL, record.ref_bases.len());
            }
            (_, -1, _) | (_, _, -1) => {
                panic!("Does not support partial null GT like ./1")
            }
            (false, 0, 0) => {
                classification.change = Change::Ref;
                if classification.has_indel_alleles
                    && is_standard_indel(&record.ref_bases, &record.alt)
                {
                    classification.ref_bases = record
                        .ref_bases
                        .chars()
                        .skip(1)
                        .collect::<String>()
                        .to_string();
                    classification.pos += 1;
                }
                classification.new_bases = classification.ref_bases.clone();
            }
            (false, i, j) if i == j => {
                classification
                    .set_ref_alt_and_simplify(&record.ref_bases, &record.alt[(i - 1) as usize]);
            }
            (false, i, j) => {
                panic!("Not considered het yet alleles are not equal: {i} {j}");
            }
            (true, i, j) => {
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

                        if classification.has_indel_alleles
                            && is_standard_indel(&record.ref_bases, &record.alt)
                        {
                            classification.ref_bases = record
                                .ref_bases
                                .chars()
                                .skip(1)
                                .collect::<String>()
                                .to_string();
                            classification.pos += 1;
                        }

                        classification.new_bases = repeat_char(HET, classification.ref_bases.len());
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

use std::iter;

use crate::consensus::parameter_struct::MinorPopParams;
use crate::vcf::{Genotype, VariantRecord};

use super::ConsensusParams;
use super::{HET_GT, MIXED, NULL, REMOVED_MINORS};

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
        if alt_bases
            .chars()
            .any(|c| !['A', 'C', 'G', 'T'].contains(&c))
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
            has_minor_population: false,
            is_filtered: false,
            has_indel_alleles: false,
        }
    }
}
impl Classification {
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
    pub minor_pop_threshold: Option<MinorPopParams>,
    pub remove_minor_pops: bool,     // If true will remove minor alleles
    pub remove_sub_minor_pops: bool, // If true will remove minor alleles below the threshold
}

impl Classifier {
    pub fn new(params: &ConsensusParams, is_support: bool) -> Self {
        if is_support {
            return Classifier {
                params: params.clone(),
                minor_pop_threshold: params.minor_pop_thresholds.clone(),
                remove_minor_pops: params.remove_minor_pops_in_support,
                remove_sub_minor_pops: params.remove_sub_minor_pops,
            };
        }

        Classifier {
            params: params.clone(),
            minor_pop_threshold: params.minor_pop_thresholds.clone(),
            remove_minor_pops: false,
            remove_sub_minor_pops: params.remove_sub_minor_pops,
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

    /// Determines how to apply a variant based on parameters set
    ///
    /// Only considers the variant in isolation, so will not check for clashes etc.
    /// A Classification will contain the new bases to apply, the change type,
    /// and some flags to aid in the consensus process
    pub fn classify(&self, record: &mut VariantRecord) -> Classification {
        let mut gt = record.genotype().expect("Genotype not found").clone();

        // Fail fast on strange genotypes
        if (gt.allele1 == -1 || gt.allele2 == -1) && gt.allele1 != gt.allele2 {
            panic!(
                "Half null genotype not supported at {}:{}. Found GT {}/{}",
                record.chrom,
                record.pos_idx(),
                gt.allele1,
                gt.allele2
            );
        }

        // We do not handle het genotypes properly, print warning
        if gt.is_het() {
            println!(
                "Het genotype record found at {}:{} with GT {}/{}. Setting to first allele and adding {HET_GT} flag",
                record.chrom,
                record.pos_idx(),
                gt.allele1,
                gt.allele2
            );
            record.filter.push(HET_GT.to_string());
            gt = Genotype {
                allele1: gt.allele1,
                allele2: gt.allele1,
            };
            record.set_genotype(gt.clone());
        }

        // updates the flags in the record
        record.filter = self.get_flags(record);

        let mut classification = Classification {
            pos: record.pos_idx(),
            ref_bases: record.ref_bases.clone(),
            has_indel_alleles: record.is_indel(), // True if any of the alleles are multiple bases
            is_filtered: !record.filter.is_empty(),
            ..Default::default()
        };

        let mut minor_alleles = Vec::new();
        let mut sub_minor_alleles = Vec::new();

        if let Some(minor_threshold) = &self.minor_pop_threshold
            && let Some(allelic_depths) = record.allele_depths()
        {
            for (i, depth) in allelic_depths.iter().enumerate() {
                if i as i32 == gt.allele1 {
                    continue; // Skip alleles in GT
                }

                if *depth < minor_threshold.threshold {
                    sub_minor_alleles.push(i);
                    continue;
                }

                if let Some(min_frs) = minor_threshold.min_frs {
                    let total_depth: i32 = allelic_depths.iter().sum();
                    if total_depth > 0 && (*depth as f32 / total_depth as f32) < min_frs {
                        sub_minor_alleles.push(i);
                        continue;
                    }
                }

                if let Some(strand_bias) = minor_threshold.strand_bias
                    && let Some((forward, reverse)) = record.strand_depths()
                {
                    let forward_depth = forward[i];
                    let reverse_depth = reverse[i];
                    let min_strand = forward_depth.min(reverse_depth) as f32;

                    let total = (forward_depth + reverse_depth) as f32;

                    if min_strand / total < strand_bias {
                        sub_minor_alleles.push(i);
                        continue;
                    }
                }

                // Allele is minor if it passes all thresholds
                minor_alleles.push(i);
            }
        } else {
            // all non-gt alt alleles are considered sub_minor
            for i in 1..record.alt.len() + 1 {
                if i as i32 == gt.allele1 {
                    continue; // Skip alleles in GT
                }
                sub_minor_alleles.push(i);
            }
        }

        classification.has_minor_population = !minor_alleles.is_empty() && !self.remove_minor_pops;

        // Remove minor alleles if specified, otherwise just add flag
        let mut alleles_to_remove = Vec::new();
        if self.remove_minor_pops {
            alleles_to_remove.extend(minor_alleles);
        }
        if self.remove_sub_minor_pops {
            alleles_to_remove.extend(sub_minor_alleles);
        }
        // never remove allele 0
        alleles_to_remove.retain(|&i| i != 0);
        alleles_to_remove.sort();
        alleles_to_remove.reverse(); // Remove from highest index to avoid messing up indices of other alleles

        if !alleles_to_remove.is_empty() {
            let info_value = if let Some(allele_depths) = record.allele_depths() {
                alleles_to_remove
                    .iter()
                    .map(|i| {
                        format!(
                            "{}({})",
                            record.get_allele_bases(*i as i32).expect("missing allele"),
                            allele_depths[*i]
                        )
                    })
                    .collect::<Vec<String>>()
            } else {
                alleles_to_remove
                    .iter()
                    .map(|i| {
                        record
                            .get_allele_bases(*i as i32)
                            .expect("missing alt")
                            .to_string()
                    })
                    .collect::<Vec<String>>()
            };

            // Add to Info field
            record.info.insert(
                REMOVED_MINORS.to_string(),
                crate::vcf::RecordValue::StringArray(info_value),
            );

            for i in alleles_to_remove {
                record.remove_allele(i as i32);
            }

            // update gt
            gt = record.genotype().expect("Genotype not found").clone();
        }

        if classification.has_minor_population {
            record
                .info
                .insert(MIXED.to_string(), crate::vcf::RecordValue::Flag);
        }

        /// Indel is standard if ref and all alts start with the same base
        /// idea being that first base is not really part of the change
        fn is_standard_indel(ref_bases: &str, alt_bases: &[String]) -> bool {
            if let Some(first_base) = ref_bases.chars().next()
                && alt_bases.iter().all(|alt| alt.starts_with(first_base))
            {
                return true;
            }
            return false;
        }

        // If filtered, just mask the site
        if classification.is_filtered {
            classification.change = Change::Null;
            classification.new_bases = repeat_char(NULL, record.ref_bases.len());
            return classification;
        }

        // Case match based on genotype.
        // Hets for indels vs snps may be handled different based on parameters
        match gt.allele1 {
            -1 => {
                classification.change = Change::Null;
                classification.new_bases = repeat_char(NULL, record.ref_bases.len());
            }
            0 => {
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
            i => {
                classification
                    .set_ref_alt_and_simplify(&record.ref_bases, &record.alt[(i - 1) as usize]);
            }
        }

        classification
    }
}

#[cfg(test)]
mod tests;

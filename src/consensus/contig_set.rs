use std::collections::{HashMap, HashSet};

/// Helper trait to allow extending onto a set like HashMap<String, HashSet<usize>>
///
/// e.g. `sites.extend_chrom("chr1", vec![1, 2, 3]);`
pub trait ContigSet<T> {
    fn extend_chrom<I>(&mut self, chrom: &str, iter: I)
    where
        I: IntoIterator<Item = T>;
    fn contains_loc(&self, chrom: &str, pos: T) -> bool;
    fn insert_loc(&mut self, chrom: &str, pos: T) -> bool;
}
impl<T> ContigSet<T> for HashMap<String, HashSet<T>>
where
    T: std::hash::Hash + Eq,
{
    fn extend_chrom<I>(&mut self, chrom: &str, iter: I)
    where
        I: IntoIterator<Item = T>,
    {
        let sites = self.entry(chrom.to_string()).or_default();
        sites.extend(iter);
    }

    fn contains_loc(&self, chrom: &str, pos: T) -> bool {
        if let Some(sites) = self.get(chrom) {
            return sites.contains(&pos);
        }
        return false;
    }

    fn insert_loc(&mut self, chrom: &str, pos: T) -> bool {
        let sites = self.entry(chrom.to_string()).or_default();
        sites.insert(pos)
    }
}
impl<'a, T> ContigSet<&'a T> for HashMap<String, HashSet<T>>
where
    T: std::hash::Hash + Eq + Clone,
{
    fn extend_chrom<I>(&mut self, chrom: &str, iter: I)
    where
        I: IntoIterator<Item = &'a T>,
    {
        let sites = self.entry(chrom.to_string()).or_default();
        sites.extend(iter.into_iter().cloned());
    }

    fn contains_loc(&self, chrom: &str, pos: &'a T) -> bool {
        if let Some(sites) = self.get(chrom) {
            return sites.contains(pos);
        }
        return false;
    }

    fn insert_loc(&mut self, chrom: &str, pos: &'a T) -> bool {
        let sites = self.entry(chrom.to_string()).or_default();
        sites.insert(pos.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_set() {
        let mut sites: HashMap<String, HashSet<usize>> = HashMap::new();
        sites.extend_chrom("chr1", vec![1, 2, 3]);
        assert_eq!(sites.contains_loc("chr1", 1), true);
        assert_eq!(sites.contains_loc("chr1", 4), false);
        assert_eq!(sites.insert_loc("chr1", 4), true);
        assert_eq!(sites.insert_loc("chr1", 4), false);
        assert_eq!(sites.contains_loc("chr1", 4), true);
        assert_eq!(sites.contains_loc("chr2", 4), false);
        assert_eq!(sites.insert_loc("chr2", 1), true);

        let expectations: HashMap<String, HashSet<usize>> = HashMap::from([
            ("chr1".to_string(), [1, 2, 3, 4].iter().cloned().collect()),
            ("chr2".to_string(), [1].iter().cloned().collect()),
        ]);
        assert_eq!(sites, expectations);
    }
}

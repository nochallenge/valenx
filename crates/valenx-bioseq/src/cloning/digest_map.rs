//! Restriction map and virtual-gel fragment computation.
//!
//! Given a sequence and a set of enzymes, [`restriction_map`] computes
//! every fragment a complete digest would produce — its size and end
//! types — and presents them as a virtual agarose gel
//! ([`virtual_gel`]: fragments sorted large→small, the way a gel
//! ladder reads).

use crate::cloning::restriction::{self, CutSite, Enzyme, OverhangType};
use crate::error::Result;
use crate::seq::Seq;

/// One fragment produced by a digest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fragment {
    /// 0-based start coordinate on the original sequence.
    pub start: usize,
    /// 0-based end coordinate (exclusive).
    pub end: usize,
    /// Fragment length in base pairs.
    pub length: usize,
    /// Overhang type at the fragment's 5′ (start) end. `None` for the
    /// ends of a linear molecule (an uncut native end).
    pub left_end: Option<OverhangType>,
    /// Overhang type at the fragment's 3′ (end) end.
    pub right_end: Option<OverhangType>,
}

/// A complete restriction map of a sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RestrictionMap {
    /// Every cut site found, sorted by position.
    pub cut_sites: Vec<CutSite>,
    /// The fragments produced, in sequence order.
    pub fragments: Vec<Fragment>,
    /// `true` if the digested molecule was circular.
    pub circular: bool,
}

impl RestrictionMap {
    /// Fragment count.
    pub fn fragment_count(&self) -> usize {
        self.fragments.len()
    }

    /// Total length of all fragments (should equal the sequence
    /// length).
    pub fn total_length(&self) -> usize {
        self.fragments.iter().map(|f| f.length).sum()
    }
}

/// Computes the restriction map of `seq` for a complete digest with
/// every enzyme in `enzymes`.
///
/// The top-strand cut position of each [`CutSite`] is the fragment
/// boundary. A linear molecule with `k` cuts yields `k + 1`
/// fragments; a circular molecule with `k ≥ 1` cuts yields `k`
/// fragments (and `1` — the whole circle, uncut — if `k == 0`).
pub fn restriction_map(seq: &Seq, enzymes: &[&Enzyme]) -> Result<RestrictionMap> {
    let n = seq.len();
    let mut cut_sites: Vec<CutSite> = Vec::new();
    for e in enzymes {
        cut_sites.extend(restriction::digest(seq, e)?);
    }
    cut_sites.sort_by_key(|c| c.top_cut_pos);

    let circular = seq.is_circular();
    let fragments = build_fragments(&cut_sites, n, circular);
    Ok(RestrictionMap {
        cut_sites,
        fragments,
        circular,
    })
}

/// Builds fragments from sorted cut positions.
fn build_fragments(cut_sites: &[CutSite], n: usize, circular: bool) -> Vec<Fragment> {
    if n == 0 {
        return Vec::new();
    }
    // Boundary positions (top-strand cuts), deduplicated.
    let mut cuts: Vec<usize> = cut_sites.iter().map(|c| c.top_cut_pos).collect();
    cuts.sort_unstable();
    cuts.dedup();

    let mut fragments = Vec::new();
    if cuts.is_empty() {
        // Undigested — one fragment spanning the whole molecule.
        fragments.push(Fragment {
            start: 0,
            end: n,
            length: n,
            left_end: None,
            right_end: None,
        });
        return fragments;
    }

    if circular {
        // Walk consecutive cut pairs; the last fragment wraps the
        // origin back to the first cut.
        for i in 0..cuts.len() {
            let start = cuts[i];
            let end = cuts[(i + 1) % cuts.len()];
            let length = if end > start {
                end - start
            } else {
                n - start + end // wraps the origin
            };
            fragments.push(Fragment {
                start,
                end,
                length,
                // Every cut end on a circle is an enzyme-generated end.
                left_end: Some(overhang_at(cut_sites, start)),
                right_end: Some(overhang_at(cut_sites, end)),
            });
        }
    } else {
        // Linear: 0 -> cut0 -> cut1 -> ... -> n.
        let mut boundaries = vec![0usize];
        boundaries.extend(cuts.iter().copied());
        boundaries.push(n);
        boundaries.dedup();
        for w in boundaries.windows(2) {
            let (start, end) = (w[0], w[1]);
            let left = if start == 0 {
                None
            } else {
                Some(overhang_at(cut_sites, start))
            };
            let right = if end == n {
                None
            } else {
                Some(overhang_at(cut_sites, end))
            };
            fragments.push(Fragment {
                start,
                end,
                length: end - start,
                left_end: left,
                right_end: right,
            });
        }
    }
    fragments
}

/// The overhang type of whichever cut produced the boundary at `pos`.
fn overhang_at(cut_sites: &[CutSite], pos: usize) -> OverhangType {
    cut_sites
        .iter()
        .find(|c| c.top_cut_pos == pos)
        .map(|c| c.overhang)
        .unwrap_or(OverhangType::Blunt)
}

/// A virtual-gel lane: fragment sizes sorted large→small, the way they
/// migrate on an agarose gel (and the way a ladder is read).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VirtualGel {
    /// Fragment sizes in bp, descending.
    pub band_sizes: Vec<usize>,
}

impl VirtualGel {
    /// Number of distinct bands (sizes that differ; co-migrating
    /// fragments of equal size count once).
    pub fn distinct_bands(&self) -> usize {
        let mut sizes: Vec<usize> = self.band_sizes.clone();
        sizes.dedup();
        sizes.len()
    }

    /// The largest band size, or `None` for an empty gel.
    pub fn largest(&self) -> Option<usize> {
        self.band_sizes.first().copied()
    }

    /// The smallest band size, or `None` for an empty gel.
    pub fn smallest(&self) -> Option<usize> {
        self.band_sizes.last().copied()
    }
}

/// Builds a virtual gel from a restriction map — the fragment lengths
/// sorted descending.
pub fn virtual_gel(map: &RestrictionMap) -> VirtualGel {
    let mut sizes: Vec<usize> = map.fragments.iter().map(|f| f.length).collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    VirtualGel { band_sizes: sizes }
}

/// Convenience: digest `seq` with `enzymes` and return the virtual gel
/// directly.
pub fn digest_to_gel(seq: &Seq, enzymes: &[&Enzyme]) -> Result<VirtualGel> {
    Ok(virtual_gel(&restriction_map(seq, enzymes)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cloning::restriction::enzyme_by_name;
    use crate::seq::{SeqKind, Topology};

    #[test]
    fn linear_single_cut_two_fragments() {
        // EcoRI site at index 5; G^AATTC -> top cut at 6.
        let s = Seq::new(SeqKind::Dna, "AAAAAGAATTCAAAAA").unwrap(); // 16 bp
        let map = restriction_map(&s, &[enzyme_by_name("EcoRI").unwrap()]).unwrap();
        assert_eq!(map.fragment_count(), 2);
        assert_eq!(map.total_length(), 16);
        // First fragment 0..6, second 6..16.
        assert_eq!(map.fragments[0].length, 6);
        assert_eq!(map.fragments[1].length, 10);
        // Native ends are None; the cut end is a 5' overhang.
        assert_eq!(map.fragments[0].left_end, None);
        assert_eq!(map.fragments[0].right_end, Some(OverhangType::FivePrime));
    }

    #[test]
    fn linear_no_cut_one_fragment() {
        let s = Seq::new(SeqKind::Dna, "AAAAAAAAAA").unwrap();
        let map = restriction_map(&s, &[enzyme_by_name("EcoRI").unwrap()]).unwrap();
        assert_eq!(map.fragment_count(), 1);
        assert_eq!(map.fragments[0].length, 10);
    }

    #[test]
    fn circular_single_cut_one_fragment() {
        // One EcoRI site on a circular molecule -> linearized -> 1
        // fragment the size of the whole circle.
        let s = Seq::with_topology(SeqKind::Dna, "AAAAAGAATTCAAAAA", Topology::Circular).unwrap();
        let map = restriction_map(&s, &[enzyme_by_name("EcoRI").unwrap()]).unwrap();
        assert!(map.circular);
        assert_eq!(map.fragment_count(), 1);
        assert_eq!(map.fragments[0].length, 16);
    }

    #[test]
    fn circular_two_cuts_two_fragments() {
        // Two EcoRI sites on a circle -> 2 fragments.
        let s = Seq::with_topology(
            SeqKind::Dna,
            "GAATTCAAAAGAATTCAAAA", // sites at 0 and 10
            Topology::Circular,
        )
        .unwrap();
        let map = restriction_map(&s, &[enzyme_by_name("EcoRI").unwrap()]).unwrap();
        assert_eq!(map.fragment_count(), 2);
        assert_eq!(map.total_length(), 20);
    }

    #[test]
    fn double_digest_three_fragments() {
        // EcoRI then BamHI -> 2 cuts -> 3 linear fragments.
        let s = Seq::new(SeqKind::Dna, "AAAGAATTCAAAGGATCCAAA").unwrap(); // 21 bp
        let map = restriction_map(
            &s,
            &[
                enzyme_by_name("EcoRI").unwrap(),
                enzyme_by_name("BamHI").unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(map.fragment_count(), 3);
        assert_eq!(map.total_length(), 21);
    }

    #[test]
    fn virtual_gel_sorted_descending() {
        let s = Seq::new(SeqKind::Dna, "AAAGAATTCAAAGGATCCAAA").unwrap();
        let gel = digest_to_gel(
            &s,
            &[
                enzyme_by_name("EcoRI").unwrap(),
                enzyme_by_name("BamHI").unwrap(),
            ],
        )
        .unwrap();
        // Bands must be in descending order.
        for w in gel.band_sizes.windows(2) {
            assert!(w[0] >= w[1], "gel not sorted: {:?}", gel.band_sizes);
        }
        assert_eq!(gel.band_sizes.iter().sum::<usize>(), 21);
        assert!(gel.largest().unwrap() >= gel.smallest().unwrap());
    }

    #[test]
    fn distinct_bands_counts_unique_sizes() {
        let gel = VirtualGel {
            band_sizes: vec![100, 100, 50],
        };
        assert_eq!(gel.distinct_bands(), 2);
    }
}

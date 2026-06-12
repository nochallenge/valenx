//! DNA / RNA base-pair detection.
//!
//! A base pair is recognised geometrically: two nucleotides are
//! paired when their bases are coplanar-ish, close, and joined by the
//! hydrogen bonds characteristic of a Watson-Crick or wobble pair.
//!
//! The detector works on the base atoms (the purine / pyrimidine ring
//! and the H-bonding edge atoms). It checks:
//!
//! 1. the C1′–C1′ distance is in the canonical 9–11 Å window,
//! 2. at least two donor/acceptor pairs of the correct base-pair
//!    geometry sit within an H-bond distance.

use crate::structure::{Model, Residue, ResidueKind};
use nalgebra::Point3;

/// The classification of a detected base pair.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum BasePairKind {
    /// A–T / A–U Watson-Crick pair.
    WatsonCrickAT,
    /// G–C Watson-Crick pair.
    WatsonCrickGC,
    /// G–U / G–T wobble pair.
    Wobble,
    /// A geometrically detected pair that is not one of the canonical
    /// classes (a non-canonical / Hoogsteen-style pair).
    NonCanonical,
}

impl BasePairKind {
    /// Whether this pair is a canonical Watson-Crick pair.
    pub fn is_watson_crick(&self) -> bool {
        matches!(
            self,
            BasePairKind::WatsonCrickAT | BasePairKind::WatsonCrickGC
        )
    }
}

/// A detected base pair between two residues of a model.
#[derive(Clone, Debug, PartialEq)]
pub struct BasePair {
    /// `(chain_index, residue_index)` of the first base.
    pub residue_a: (usize, usize),
    /// `(chain_index, residue_index)` of the second base.
    pub residue_b: (usize, usize),
    /// The classified pair type.
    pub kind: BasePairKind,
    /// Number of detected inter-base hydrogen bonds.
    pub hbond_count: usize,
    /// C1′–C1′ distance, ångström.
    pub c1_distance: f64,
}

/// Maximum donor–acceptor distance for an inter-base hydrogen bond,
/// ångström.
const HBOND_MAX: f64 = 3.5;
/// Canonical C1′–C1′ distance window.
const C1_MIN: f64 = 8.0;
const C1_MAX: f64 = 11.5;

/// The single-letter base of a nucleotide residue (`A`, `C`, `G`,
/// `T`, `U`), or `None` if it is not a recognised nucleotide.
pub fn base_letter(residue: &Residue) -> Option<char> {
    match residue.kind() {
        ResidueKind::Dna | ResidueKind::Rna => {
            Some(crate::structure::residue_one_letter(&residue.name))
        }
        _ => None,
    }
}

/// The H-bonding edge atom names for each base, paired as
/// `(donor_or_acceptor_atoms)`. These are the standard Watson-Crick
/// edge atoms.
fn wc_edge_atoms(base: char) -> &'static [&'static str] {
    match base {
        'A' => &["N1", "N6"],
        'T' | 'U' => &["N3", "O4", "O2"],
        'G' => &["N1", "N2", "O6"],
        'C' => &["N3", "N4", "O2"],
        _ => &[],
    }
}

/// Count inter-base hydrogen bonds between two nucleotides: edge
/// atoms of `a` within [`HBOND_MAX`] of edge atoms of `b`.
fn count_hbonds(a: &Residue, b: &Residue, base_a: char, base_b: char) -> usize {
    let edge_a = wc_edge_atoms(base_a);
    let edge_b = wc_edge_atoms(base_b);
    let mut count = 0;
    for na in edge_a {
        if let Some(pa) = a.primary_atom(na) {
            for nb in edge_b {
                if let Some(pb) = b.primary_atom(nb) {
                    let d = pa.distance(pb);
                    // 2.4 A floor rejects same-atom overlaps.
                    if (2.4..=HBOND_MAX).contains(&d) {
                        count += 1;
                    }
                }
            }
        }
    }
    count
}

/// Classify a base pair from its two base letters and H-bond count.
fn classify_pair(base_a: char, base_b: char, hbonds: usize) -> BasePairKind {
    let pair = {
        let mut p = [base_a, base_b];
        p.sort_unstable();
        (p[0], p[1])
    };
    match pair {
        ('A', 'T') | ('A', 'U') if hbonds >= 2 => BasePairKind::WatsonCrickAT,
        ('C', 'G') if hbonds >= 2 => BasePairKind::WatsonCrickGC,
        ('G', 'T') | ('G', 'U') if hbonds >= 2 => BasePairKind::Wobble,
        _ => BasePairKind::NonCanonical,
    }
}

/// The C1′ (glycosidic) atom of a nucleotide, ångström. Falls back
/// to the base-ring centroid when `C1'` is absent.
pub fn c1_prime(residue: &Residue) -> Option<Point3<f64>> {
    residue
        .primary_atom("C1'")
        .or_else(|| residue.primary_atom("C1*"))
        .map(|a| a.coord)
        .or_else(|| base_ring_centroid(residue))
}

/// Centroid of a nucleotide's base-ring atoms — a stable proxy for
/// the base position when `C1'` is missing.
pub fn base_ring_centroid(residue: &Residue) -> Option<Point3<f64>> {
    const RING: &[&str] = &["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"];
    let mut acc = nalgebra::Vector3::zeros();
    let mut n = 0;
    for name in RING {
        if let Some(a) = residue.primary_atom(name) {
            acc += a.coord.coords;
            n += 1;
        }
    }
    if n >= 3 {
        Some(Point3::from(acc / n as f64))
    } else {
        None
    }
}

/// Detect every base pair in a model.
///
/// Each ordered residue pair (across all chains) is tested; a pair is
/// reported once. Detection requires the C1′–C1′ distance in the
/// canonical window and at least two inter-base hydrogen bonds.
pub fn detect_base_pairs(model: &Model) -> Vec<BasePair> {
    // Flatten nucleotide residues with their (chain, residue) index.
    struct Nuc<'a> {
        chain: usize,
        residue: usize,
        base: char,
        c1: Point3<f64>,
        data: &'a Residue,
    }
    let mut nucs: Vec<Nuc> = Vec::new();
    for (ci, chain) in model.chains.iter().enumerate() {
        for (ri, r) in chain.residues.iter().enumerate() {
            if let (Some(base), Some(c1)) = (base_letter(r), c1_prime(r)) {
                nucs.push(Nuc {
                    chain: ci,
                    residue: ri,
                    base,
                    c1,
                    data: r,
                });
            }
        }
    }

    let mut pairs = Vec::new();
    let mut used: Vec<bool> = vec![false; nucs.len()];
    for i in 0..nucs.len() {
        if used[i] {
            continue;
        }
        // Find the best partner for nucleotide i.
        let mut best: Option<(usize, usize, f64)> = None; // (j, hbonds, dist)
        for j in (i + 1)..nucs.len() {
            if used[j] {
                continue;
            }
            // Skip directly-bonded sequence neighbours in one chain.
            if nucs[i].chain == nucs[j].chain && nucs[i].residue.abs_diff(nucs[j].residue) <= 1 {
                continue;
            }
            let d = (nucs[i].c1 - nucs[j].c1).norm();
            if !(C1_MIN..=C1_MAX).contains(&d) {
                continue;
            }
            let hb = count_hbonds(nucs[i].data, nucs[j].data, nucs[i].base, nucs[j].base);
            if hb >= 2 {
                let better = match best {
                    Some((_, bh, bd)) => hb > bh || (hb == bh && d < bd),
                    None => true,
                };
                if better {
                    best = Some((j, hb, d));
                }
            }
        }
        if let Some((j, hb, d)) = best {
            used[i] = true;
            used[j] = true;
            let kind = classify_pair(nucs[i].base, nucs[j].base, hb);
            pairs.push(BasePair {
                residue_a: (nucs[i].chain, nucs[i].residue),
                residue_b: (nucs[j].chain, nucs[j].residue),
                kind,
                hbond_count: hb,
                c1_distance: d,
            });
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain};

    /// Build a minimal A nucleotide with its WC edge atoms placed.
    fn adenine(seq: i32, c1: Point3<f64>) -> Residue {
        let mut r = Residue::new("DA", seq);
        r.atoms.push(Atom::new("C1'", "C", c1));
        // N1 acceptor and N6 donor on the WC edge, offset toward +x.
        r.atoms.push(Atom::new(
            "N1",
            "N",
            c1 + nalgebra::Vector3::new(4.5, 0.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "N6",
            "N",
            c1 + nalgebra::Vector3::new(4.5, 2.8, 0.0),
        ));
        // a couple of ring atoms so the centroid fallback works too.
        r.atoms.push(Atom::new(
            "C2",
            "C",
            c1 + nalgebra::Vector3::new(3.5, -1.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "C4",
            "C",
            c1 + nalgebra::Vector3::new(3.0, 1.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "C5",
            "C",
            c1 + nalgebra::Vector3::new(3.8, 1.5, 0.0),
        ));
        r
    }

    /// Build a minimal T nucleotide positioned to pair with `adenine`.
    fn thymine(seq: i32, c1: Point3<f64>) -> Residue {
        let mut r = Residue::new("DT", seq);
        r.atoms.push(Atom::new("C1'", "C", c1));
        // Place N3 ~2.9 A from adenine N1, O4 ~2.9 A from adenine N6.
        // Adenine N1 is at x=4.5; thymine C1' will be at x=10, so its
        // edge atoms point toward -x.
        r.atoms.push(Atom::new(
            "N3",
            "N",
            c1 - nalgebra::Vector3::new(4.6, 0.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "O4",
            "O",
            c1 - nalgebra::Vector3::new(4.6, -2.8, 0.0),
        ));
        r.atoms.push(Atom::new(
            "O2",
            "O",
            c1 - nalgebra::Vector3::new(4.6, 2.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "C2",
            "C",
            c1 - nalgebra::Vector3::new(3.5, 1.0, 0.0),
        ));
        r.atoms.push(Atom::new(
            "C4",
            "C",
            c1 - nalgebra::Vector3::new(3.5, -1.0, 0.0),
        ));
        r
    }

    #[test]
    fn base_letters() {
        let a = adenine(1, Point3::origin());
        assert_eq!(base_letter(&a), Some('A'));
        let mut prot = Residue::new("ALA", 1);
        prot.atoms.push(Atom::new("CA", "C", Point3::origin()));
        assert_eq!(base_letter(&prot), None);
    }

    #[test]
    fn detects_an_at_pair() {
        // Adenine C1' at origin, thymine C1' ~10 A away: the edge
        // atoms (built above) are within H-bond distance.
        let a = adenine(1, Point3::new(0.0, 0.0, 0.0));
        let t = thymine(2, Point3::new(10.0, 0.0, 0.0));
        let mut chain_a = Chain::new("A");
        chain_a.residues.push(a);
        let mut chain_b = Chain::new("B");
        chain_b.residues.push(t);
        let mut model = Model::new(1);
        model.chains.push(chain_a);
        model.chains.push(chain_b);

        let pairs = detect_base_pairs(&model);
        assert_eq!(pairs.len(), 1, "expected one base pair");
        assert_eq!(pairs[0].kind, BasePairKind::WatsonCrickAT);
        assert!(pairs[0].kind.is_watson_crick());
        assert!(pairs[0].hbond_count >= 2);
        assert!((C1_MIN..=C1_MAX).contains(&pairs[0].c1_distance));
    }

    #[test]
    fn no_pair_when_too_far() {
        let a = adenine(1, Point3::new(0.0, 0.0, 0.0));
        let t = thymine(2, Point3::new(40.0, 0.0, 0.0));
        let mut chain_a = Chain::new("A");
        chain_a.residues.push(a);
        let mut chain_b = Chain::new("B");
        chain_b.residues.push(t);
        let mut model = Model::new(1);
        model.chains.push(chain_a);
        model.chains.push(chain_b);
        assert!(detect_base_pairs(&model).is_empty());
    }

    #[test]
    fn ring_centroid_fallback() {
        // Adenine without C1' still yields a position from its ring.
        let mut a = adenine(1, Point3::origin());
        a.atoms.retain(|at| at.name != "C1'");
        assert!(c1_prime(&a).is_some());
    }
}

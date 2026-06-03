//! Structure-based search descriptor — a FoldSeek-class structural
//! alphabet (v1).
//!
//! FoldSeek's key idea is the **3Di alphabet**: it discretises each
//! residue's *local backbone geometry* into one of a small set of
//! structural letters, turning a 3-D fold into a 1-D string that a
//! fast sequence aligner can search.
//!
//! This v1 builds a structural-alphabet string from per-residue
//! geometric features:
//!
//! - the virtual Cα pseudo-torsion of `Cα(i−1), Cα(i), Cα(i+1),
//!   Cα(i+2)`,
//! - the virtual Cα pseudo-bond-angle at `Cα(i)`,
//! - the contact direction to the spatially nearest non-local Cα.
//!
//! Those features are quantised into a 20-letter alphabet
//! (`A`–`T`). Two structures with similar folds get similar strings,
//! and the [`descriptor_identity`] of their strings is a fast
//! fold-similarity proxy.
//!
//! ## Scope of this v1
//!
//! FoldSeek's real 3Di alphabet comes from a learnt VQ-VAE over
//! interaction geometry; this is a hand-designed quantiser in the
//! same *spirit* — deterministic, no learnt weights. It captures
//! local secondary-structure context well and is a usable fast
//! pre-filter; it will not match FoldSeek's remote-homology
//! sensitivity.

use crate::error::{BiostructError, Result};
use crate::geometry::angles::dihedral;
use crate::structure::Chain;
use nalgebra::Point3;

/// The size of the structural alphabet (`A`..=`T`).
pub const ALPHABET_SIZE: usize = 20;

/// A structural-alphabet descriptor of one chain.
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralDescriptor {
    /// One structural letter per residue with a defined local
    /// geometry (chain termini are dropped).
    pub letters: String,
    /// Residue indices the letters correspond to.
    pub residue_indices: Vec<usize>,
}

impl StructuralDescriptor {
    /// Length of the descriptor string.
    pub fn len(&self) -> usize {
        self.letters.len()
    }

    /// Whether the descriptor is empty.
    pub fn is_empty(&self) -> bool {
        self.letters.is_empty()
    }
}

/// Per-residue geometric features feeding the quantiser.
struct LocalGeometry {
    /// Cα pseudo-torsion, degrees in `(-180, 180]`.
    torsion: f64,
    /// Cα pseudo-bond-angle, degrees in `[0, 180]`.
    angle: f64,
    /// Distance to the nearest non-local Cα, ångström.
    contact_dist: f64,
}

/// Compute the structural-alphabet descriptor of a protein chain.
pub fn structural_descriptor(chain: &Chain) -> Result<StructuralDescriptor> {
    // Cα trace with original residue indices.
    let trace: Vec<(usize, Point3<f64>)> = chain
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if r.is_amino_acid() {
                r.ca().map(|a| (i, a.coord))
            } else {
                None
            }
        })
        .collect();
    if trace.len() < 4 {
        return Err(BiostructError::invalid(
            "chain",
            "structural descriptor needs at least 4 Cα atoms",
        ));
    }

    let mut letters = String::new();
    let mut indices = Vec::new();

    // The torsion needs i-1..i+2; iterate the interior.
    for k in 1..trace.len() - 2 {
        let ca_prev = trace[k - 1].1;
        let ca_i = trace[k].1;
        let ca_next = trace[k + 1].1;
        let ca_next2 = trace[k + 2].1;

        let torsion = dihedral(&ca_prev, &ca_i, &ca_next, &ca_next2).unwrap_or(0.0);

        // Pseudo bond angle at ca_i.
        let v1 = (ca_prev - ca_i).normalize();
        let v2 = (ca_next - ca_i).normalize();
        let angle = v1.dot(&v2).clamp(-1.0, 1.0).acos().to_degrees();

        // Nearest non-local Cα (sequence separation >= 4).
        let mut contact_dist = f64::INFINITY;
        for (m, (_, other)) in trace.iter().enumerate() {
            if m.abs_diff(k) >= 4 {
                let d = (ca_i - other).norm();
                if d < contact_dist {
                    contact_dist = d;
                }
            }
        }
        if !contact_dist.is_finite() {
            contact_dist = 99.0;
        }

        let geom = LocalGeometry {
            torsion,
            angle,
            contact_dist,
        };
        letters.push(quantize(&geom));
        indices.push(trace[k].0);
    }

    Ok(StructuralDescriptor {
        letters,
        residue_indices: indices,
    })
}

/// Quantise a residue's local geometry into one structural letter.
///
/// The quantiser bins the three features into a `5 × 2 × 2 = 20`
/// grid, mapped onto `A`..=`T`. The torsion (the dominant
/// secondary-structure discriminator) gets five bins; the bond angle
/// and the contact distance get two each.
fn quantize(g: &LocalGeometry) -> char {
    // Torsion bins: helix-like, sheet-like and three others.
    let t_bin = if (-90.0..=-30.0).contains(&g.torsion) {
        0 // right-handed helix-ish region
    } else if g.torsion > 90.0 || g.torsion < -150.0 {
        1 // extended / sheet-ish
    } else if (-150.0..=-90.0).contains(&g.torsion) {
        2
    } else if (-30.0..=60.0).contains(&g.torsion) {
        3
    } else {
        4
    };
    // Bond-angle bin: tightly bent vs open.
    let a_bin = if g.angle < 95.0 { 0 } else { 1 };
    // Contact bin: buried (a close non-local neighbour) vs exposed.
    let c_bin = if g.contact_dist < 8.0 { 0 } else { 1 };

    let index = t_bin * 4 + a_bin * 2 + c_bin;
    // index in 0..20 -> 'A'..'T'.
    (b'A' + index.min(ALPHABET_SIZE - 1) as u8) as char
}

/// Fraction of identical letters between two equal-length structural
/// descriptors, in `[0, 1]`. Returns an error on a length mismatch.
pub fn descriptor_identity(a: &StructuralDescriptor, b: &StructuralDescriptor) -> Result<f64> {
    if a.len() != b.len() {
        return Err(BiostructError::invalid(
            "descriptor",
            "descriptors must be equal length for identity",
        ));
    }
    if a.is_empty() {
        return Ok(1.0);
    }
    let same = a
        .letters
        .chars()
        .zip(b.letters.chars())
        .filter(|(x, y)| x == y)
        .count();
    Ok(same as f64 / a.len() as f64)
}

/// A fast ungapped fold-similarity score between two structural
/// descriptors of *possibly different* lengths: the best identity
/// over every offset of the shorter string against the longer one.
///
/// This is the FoldSeek-class "fast pre-filter" — cheap, ungapped,
/// and a usable upper-bound screen before a full structural
/// alignment.
pub fn best_ungapped_similarity(
    a: &StructuralDescriptor,
    b: &StructuralDescriptor,
) -> f64 {
    let (short, long) = if a.len() <= b.len() {
        (a.letters.as_bytes(), b.letters.as_bytes())
    } else {
        (b.letters.as_bytes(), a.letters.as_bytes())
    };
    if short.is_empty() {
        return 0.0;
    }
    let mut best = 0.0_f64;
    for offset in 0..=(long.len() - short.len()) {
        let matches = short
            .iter()
            .zip(&long[offset..offset + short.len()])
            .filter(|(x, y)| x == y)
            .count();
        let frac = matches as f64 / short.len() as f64;
        if frac > best {
            best = frac;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Residue};

    fn helix_chain(n: usize) -> Chain {
        let mut c = Chain::new("A");
        for i in 0..n {
            let t = i as f64 * 100.0_f64.to_radians();
            let ca = Point3::new(2.3 * t.cos(), 2.3 * t.sin(), i as f64 * 1.5);
            let mut r = Residue::new("ALA", i as i32 + 1);
            r.atoms.push(Atom::new("CA", "C", ca));
            c.residues.push(r);
        }
        c
    }

    fn extended_chain(n: usize) -> Chain {
        // A near-straight, slightly zig-zag extended strand.
        let mut c = Chain::new("B");
        for i in 0..n {
            let ca = Point3::new(
                i as f64 * 3.5,
                if i % 2 == 0 { 0.0 } else { 0.6 },
                0.0,
            );
            let mut r = Residue::new("VAL", i as i32 + 1);
            r.atoms.push(Atom::new("CA", "C", ca));
            c.residues.push(r);
        }
        c
    }

    #[test]
    fn descriptor_has_one_letter_per_interior_residue() {
        let chain = helix_chain(20);
        let d = structural_descriptor(&chain).unwrap();
        // interior count = n - 3 (drop first and last two).
        assert_eq!(d.len(), 20 - 3);
        assert!(d.letters.chars().all(|c| ('A'..='T').contains(&c)));
    }

    #[test]
    fn identical_chains_have_identity_one() {
        let chain = helix_chain(25);
        let d = structural_descriptor(&chain).unwrap();
        assert!((descriptor_identity(&d, &d).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn helix_and_strand_differ() {
        let helix = structural_descriptor(&helix_chain(25)).unwrap();
        let strand = structural_descriptor(&extended_chain(25)).unwrap();
        let id = descriptor_identity(&helix, &strand).unwrap();
        // Different folds must not be identical.
        assert!(id < 0.9, "helix vs strand identity unexpectedly {id}");
    }

    #[test]
    fn helix_descriptor_is_internally_consistent() {
        // A regular helix should produce a low-entropy descriptor —
        // most interior residues share the dominant letter.
        let d = structural_descriptor(&helix_chain(30)).unwrap();
        let bytes = d.letters.as_bytes();
        let mut counts = [0usize; 128];
        for &b in bytes {
            counts[b as usize] += 1;
        }
        let max = counts.iter().copied().max().unwrap();
        assert!(
            max as f64 / bytes.len() as f64 > 0.5,
            "helix descriptor not dominated by one letter"
        );
    }

    #[test]
    fn ungapped_similarity_self_is_one() {
        let d = structural_descriptor(&helix_chain(20)).unwrap();
        assert!((best_ungapped_similarity(&d, &d) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ungapped_similarity_finds_subfold() {
        // A 15-residue helix is a sub-fold of a 30-residue helix:
        // the best offset alignment should be near-perfect.
        let small = structural_descriptor(&helix_chain(18)).unwrap();
        let big = structural_descriptor(&helix_chain(36)).unwrap();
        assert!(best_ungapped_similarity(&small, &big) > 0.8);
    }

    #[test]
    fn rejects_short_chains() {
        assert!(structural_descriptor(&helix_chain(3)).is_err());
    }

    #[test]
    fn identity_rejects_length_mismatch() {
        let a = structural_descriptor(&helix_chain(20)).unwrap();
        let b = structural_descriptor(&helix_chain(25)).unwrap();
        assert!(descriptor_identity(&a, &b).is_err());
    }
}

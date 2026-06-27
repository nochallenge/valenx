//! Geometric interface complementarity between two structural chains.
//!
//! When experimental or docked coordinates exist for both partners, the
//! sequence-only coevolution signal can be reinforced by *geometry*: do
//! the two chains actually pack against each other with a well-formed,
//! shape-complementary interface? This module computes a single bounded
//! `[0, 1]` complementarity proxy from the inter-chain atom contacts of
//! two [`Chain`]s.
//!
//! ## What it measures
//!
//! For every heavy-atom pair `(a in chain A, b in chain B)` within a
//! contact shell it accumulates a smooth, distance-weighted contact
//! kernel that peaks at favourable van-der-Waals packing (~3.5–4 Å) and
//! falls off to zero by the cutoff. The summed kernel is normalised by
//! the interface size (the count of interface residues) and squashed to
//! `[0, 1]`. High = a dense, well-packed interface; low = chains that
//! barely touch.
//!
//! ## Honest scope
//!
//! This is a **packing-density proxy**, not the Lawrence-Colman shape
//! complementarity statistic `Sc` (which needs a molecular surface and
//! per-point nearest-neighbour dot products). It rewards close,
//! many-atom contact and penalises steric overlap, but it does not
//! model desolvation, electrostatics, or true surface fit. Use it as a
//! geometric tie-breaker on top of coevolution, never as a standalone
//! affinity or "binds" verdict.

use valenx_biostruct::structure::Chain;

use crate::error::PpiError;

/// Heavy-atom contact-shell cutoff, ångström. Pairs beyond this
/// contribute nothing. 6 Å comfortably covers van-der-Waals contact and
/// the first hydration-bridge shell.
pub const CONTACT_CUTOFF: f64 = 6.0;

/// Distance (Å) of peak favourability for the contact kernel — typical
/// heavy-atom van-der-Waals contact.
pub const KERNEL_PEAK: f64 = 3.8;

/// Centre-to-centre distance (Å) below which a heavy-atom pair is
/// treated as a steric clash and contributes a negative term.
pub const CLASH_DISTANCE: f64 = 2.4;

/// How much more a steric clash counts than a single favourable
/// contact. A clash means the two chains physically overlap — an
/// infeasible pose — so it must dominate the few favourable contacts
/// that surround it, not be averaged away by them. Set so that a
/// clashing interface scores below a clean one of the same contact
/// count.
pub const CLASH_PENALTY_WEIGHT: f64 = 5.0;

/// A geometric interface-complementarity report for one chain pair.
#[derive(Clone, Debug, PartialEq)]
pub struct Complementarity {
    /// The squashed complementarity score in `[0, 1]` (higher = a
    /// denser, better-packed interface).
    pub value: f64,
    /// Number of interface residues (chain-A residues with at least one
    /// heavy atom within [`CONTACT_CUTOFF`] of chain B, plus the
    /// symmetric chain-B count).
    pub interface_residues: usize,
    /// Number of heavy-atom contact pairs within [`CONTACT_CUTOFF`].
    pub contact_pairs: usize,
    /// Number of heavy-atom pairs closer than [`CLASH_DISTANCE`].
    pub clashes: usize,
}

/// Smooth contact kernel: `1` at [`KERNEL_PEAK`], decaying to `0` at
/// [`CONTACT_CUTOFF`], and going negative for clashes inside
/// [`CLASH_DISTANCE`]. Half-cosine taper keeps it continuous.
fn contact_kernel(d: f64) -> f64 {
    if d <= CLASH_DISTANCE {
        // Clash penalty: linearly worse as atoms overlap (-1 at d = 0),
        // scaled by CLASH_PENALTY_WEIGHT so an infeasible overlapping
        // pose cannot be rescued by the favourable contacts around it.
        return -CLASH_PENALTY_WEIGHT * (1.0 - d / CLASH_DISTANCE);
    }
    if d >= CONTACT_CUTOFF {
        return 0.0;
    }
    if d <= KERNEL_PEAK {
        // Rising edge: clash distance -> peak.
        let t = (d - CLASH_DISTANCE) / (KERNEL_PEAK - CLASH_DISTANCE);
        // Smooth 0->1.
        0.5 * (1.0 - (std::f64::consts::PI * (1.0 - t)).cos())
    } else {
        // Falling edge: peak -> cutoff (half cosine 1->0).
        let t = (d - KERNEL_PEAK) / (CONTACT_CUTOFF - KERNEL_PEAK);
        0.5 * (1.0 + (std::f64::consts::PI * t).cos())
    }
}

/// Compute geometric interface complementarity between two chains.
///
/// Only heavy atoms (non-hydrogen) participate. Both chains must contain
/// at least one atom or the call fails loud — the complementarity term
/// has no meaning without coordinates.
///
/// # Errors
/// [`PpiError::MissingStructure`] if either chain has no atoms.
pub fn interface_complementarity(a: &Chain, b: &Chain) -> Result<Complementarity, PpiError> {
    let atoms_a = heavy_atoms(a);
    let atoms_b = heavy_atoms(b);
    if atoms_a.is_empty() {
        return Err(PpiError::MissingStructure { what: "chain_a" });
    }
    if atoms_b.is_empty() {
        return Err(PpiError::MissingStructure { what: "chain_b" });
    }

    let cutoff_sq = CONTACT_CUTOFF * CONTACT_CUTOFF;
    let mut kernel_sum = 0.0;
    let mut contact_pairs = 0usize;
    let mut clashes = 0usize;
    // Track which residue index on each side participates, so we can
    // count interface residues without double counting.
    let mut a_iface = vec![false; a.residues.len()];
    let mut b_iface = vec![false; b.residues.len()];

    for &(ra, ref pa) in &atoms_a {
        for &(rb, ref pb) in &atoms_b {
            let d2 = (pa - pb).norm_squared();
            if d2 > cutoff_sq {
                continue;
            }
            let d = d2.sqrt();
            kernel_sum += contact_kernel(d);
            contact_pairs += 1;
            if d < CLASH_DISTANCE {
                clashes += 1;
            }
            a_iface[ra] = true;
            b_iface[rb] = true;
        }
    }

    let interface_residues =
        a_iface.iter().filter(|&&x| x).count() + b_iface.iter().filter(|&&x| x).count();

    // Normalise the packing by interface size: contacts-per-interface-
    // residue. A well-packed protein interface carries on the order of
    // a handful of favourable heavy-atom contacts per interface residue,
    // so dividing by interface_residues and squashing keeps the score
    // scale-free across interface sizes. No interface -> 0.
    let value = if interface_residues == 0 {
        0.0
    } else {
        let density = kernel_sum / interface_residues as f64;
        // Squash: density of ~2 favourable contacts/residue -> ~0.86.
        // Clashes drive kernel_sum (hence density) down, lowering the
        // score, so a clashing pose never scores like a clean one.
        let d = density.max(0.0);
        1.0 - (-d).exp()
    };

    Ok(Complementarity {
        value,
        interface_residues,
        contact_pairs,
        clashes,
    })
}

/// Heavy atoms of a chain as `(residue_index, coord)` pairs.
fn heavy_atoms(chain: &Chain) -> Vec<(usize, nalgebra::Point3<f64>)> {
    let mut out = Vec::new();
    for (ri, res) in chain.residues.iter().enumerate() {
        for atom in &res.atoms {
            if atom.is_hydrogen() {
                continue;
            }
            out.push((ri, atom.coord));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;
    use valenx_biostruct::structure::{Atom, Chain, Residue};

    /// A chain of `n` CA atoms laid along +x starting at `origin`.
    fn ca_chain(id: &str, origin: Point3<f64>, n: usize, step: f64) -> Chain {
        let mut c = Chain::new(id);
        for i in 0..n {
            let mut r = Residue::new("ALA", i as i32 + 1);
            r.atoms.push(Atom::new(
                "CA",
                "C",
                origin + nalgebra::Vector3::new(i as f64 * step, 0.0, 0.0),
            ));
            c.residues.push(r);
        }
        c
    }

    #[test]
    fn missing_atoms_fail_loud() {
        let empty = Chain::new("A");
        let real = ca_chain("B", Point3::new(0.0, 0.0, 0.0), 3, 3.8);
        let err = interface_complementarity(&empty, &real).unwrap_err();
        assert_eq!(err.code(), "missing_structure");
    }

    #[test]
    fn touching_chains_score_higher_than_far_chains() {
        let a = ca_chain("A", Point3::new(0.0, 0.0, 0.0), 4, 3.8);
        // B packed ~3.8 A away in +y from A — a real interface.
        let near = ca_chain("B", Point3::new(0.0, 3.8, 0.0), 4, 3.8);
        // B pushed 40 A away — no contact.
        let far = ca_chain("B", Point3::new(0.0, 40.0, 0.0), 4, 3.8);

        let near_s = interface_complementarity(&a, &near).unwrap();
        let far_s = interface_complementarity(&a, &far).unwrap();

        assert!(near_s.value > far_s.value);
        assert!(near_s.contact_pairs > 0);
        assert_eq!(far_s.value, 0.0);
        assert_eq!(far_s.contact_pairs, 0);
    }

    #[test]
    fn clashing_interface_scores_below_clean_one() {
        let a = ca_chain("A", Point3::new(0.0, 0.0, 0.0), 4, 3.8);
        let clean = ca_chain("B", Point3::new(0.0, 3.8, 0.0), 4, 3.8);
        // Overlapping B at ~1 A — heavy clashes.
        let clash = ca_chain("B", Point3::new(0.0, 1.0, 0.0), 4, 3.8);

        let clean_s = interface_complementarity(&a, &clean).unwrap();
        let clash_s = interface_complementarity(&a, &clash).unwrap();
        assert!(clash_s.clashes > 0);
        assert!(clash_s.value < clean_s.value);
    }
}

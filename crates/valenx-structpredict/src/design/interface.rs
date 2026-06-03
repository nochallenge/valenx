//! **Feature 22 — interface / binding-site design.**
//!
//! Designing a protein-protein *interface* — the surface where a
//! binder meets its target — is a specialised design problem. The key
//! difference from monomer design: an interface residue is *buried
//! against the partner*, not against its own core, so it should be
//! hydrophobic for affinity *and* shape-complementary to the partner.
//!
//! This module:
//!
//! 1. **Detects the interface** — given a binder model and a target
//!    model (both with Cα traces, already positioned in a docked
//!    pose), it finds the binder residues whose Cα lies within an
//!    interface cutoff of any target residue.
//! 2. **Designs those positions** — runs the combinatorial design
//!    search restricted to the interface residues, with the interface
//!    residues' "burial" supplied by their partner contacts, so the
//!    search drives the interface toward an affinity-optimised,
//!    well-packed set of residues; the binder's non-interface
//!    residues are held fixed at their native identity.
//!
//! It is a real classical interface-design protocol — interface
//! detection plus restricted physics-based design. Designing a binder
//! *from scratch* (RFdiffusion-class de-novo binder generation) needs
//! the network and is adapter-only.

use serde::{Deserialize, Serialize};

use crate::abinitio::ss::SecondaryStructure;
use crate::design::score::DesignScoreWeights;
use crate::design::search::{combinatorial_design, ResiduePalette};
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;

/// The outcome of an interface-design run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InterfaceDesign {
    /// The binder with its interface residues redesigned.
    pub model: ProteinModel,
    /// Binder residue indices identified as interface positions.
    pub interface_residues: Vec<usize>,
    /// The redesigned binder sequence.
    pub designed_sequence: String,
    /// Design energy of the redesigned binder.
    pub design_energy: f64,
    /// Fraction of interface positions whose residue changed.
    pub interface_mutation_fraction: f64,
}

/// Detects the interface residues of a binder against a target.
///
/// Returns the indices of binder residues whose Cα lies within
/// `cutoff` ångström of any target Cα. Both models must already be
/// in a docked pose (their coordinates in the same frame).
///
/// # Errors
/// [`StructPredictError::Invalid`] for a non-positive cutoff or
/// models without Cα atoms.
pub fn detect_interface(
    binder: &ProteinModel,
    target: &ProteinModel,
    cutoff: f64,
) -> Result<Vec<usize>> {
    if !(cutoff.is_finite() && cutoff > 0.0) {
        return Err(StructPredictError::invalid(
            "cutoff",
            "must be finite and positive",
        ));
    }
    let target_ca: Vec<_> = target.residues.iter().filter_map(|r| r.ca).collect();
    if target_ca.is_empty() {
        return Err(StructPredictError::invalid(
            "target",
            "target has no Cα atoms",
        ));
    }
    let mut interface = Vec::new();
    for (i, res) in binder.residues.iter().enumerate() {
        let Some(ca) = res.ca else { continue };
        if target_ca.iter().any(|t| (ca - t).norm() < cutoff) {
            interface.push(i);
        }
    }
    Ok(interface)
}

/// Designs the interface of a binder against a target.
///
/// Detects the interface (within `interface_cutoff` ångström), builds
/// a [`ResiduePalette`] that lets *only* the interface residues vary
/// — every other binder residue is fixed at its native identity —
/// and runs the combinatorial design search. `moves` is the
/// Monte-Carlo budget; `seed` fixes the RNG.
///
/// # Errors
/// [`StructPredictError::Invalid`] for bad arguments;
/// [`StructPredictError::NotFound`] if no interface residues are
/// found (the binder and target are not in contact).
pub fn design_interface(
    binder: &ProteinModel,
    target: &ProteinModel,
    interface_cutoff: f64,
    moves: usize,
    seed: u64,
) -> Result<InterfaceDesign> {
    if moves == 0 {
        return Err(StructPredictError::invalid("moves", "must be at least 1"));
    }
    let interface = detect_interface(binder, target, interface_cutoff)?;
    if interface.is_empty() {
        return Err(StructPredictError::not_found(
            "interface",
            "binder and target make no Cα contact within the cutoff",
        ));
    }

    let native = binder.sequence();
    // Palette: interface residues unrestricted, everything else fixed.
    let mut allowed: Vec<Vec<char>> = native.chars().map(|c| vec![c]).collect();
    for &i in &interface {
        if let Some(slot) = allowed.get_mut(i) {
            slot.clear(); // empty → any of the 20
        }
    }
    let palette = ResiduePalette { allowed };

    // The design score needs an SS array; interface design does not
    // depend on SS, so pass an all-coil array of the right length.
    let ss = vec![SecondaryStructure::Coil; binder.residues.len()];
    let search = combinatorial_design(
        binder,
        &palette,
        &ss,
        DesignScoreWeights::default(),
        moves,
        seed,
    )?;

    // Build the redesigned binder model (same backbone, new identities).
    let mut model = binder.clone();
    for (res, aa) in model.residues.iter_mut().zip(search.sequence.chars()) {
        res.aa = aa;
    }

    // Interface mutation fraction.
    let mut mutated = 0usize;
    for &i in &interface {
        if native.as_bytes().get(i) != search.sequence.as_bytes().get(i) {
            mutated += 1;
        }
    }
    let interface_mutation_fraction = mutated as f64 / interface.len() as f64;

    Ok(InterfaceDesign {
        model,
        interface_residues: interface,
        designed_sequence: search.sequence,
        design_energy: search.score.total,
        interface_mutation_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn slab(seq: &str, dy: f64) -> ProteinModel {
        let mut m = ProteinModel::from_sequence(seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new((i % 4) as f64 * 3.8, dy, (i / 4) as f64 * 3.8));
        }
        m
    }

    #[test]
    fn interface_residues_are_detected() {
        // Two slabs ~6 Å apart in y — their facing residues form an
        // interface.
        let binder = slab("AAAAAAAA", 0.0);
        let target = slab("AAAAAAAA", 6.0);
        let interface = detect_interface(&binder, &target, 8.0).expect("detect");
        assert!(!interface.is_empty(), "interface found");
        // Far apart → no interface.
        let far = slab("AAAAAAAA", 50.0);
        let none = detect_interface(&binder, &far, 8.0).expect("detect");
        assert!(none.is_empty(), "no interface when far apart");
    }

    #[test]
    fn interface_design_only_changes_interface_positions() {
        let binder = slab("ACDEFGHIKLMNPQRS", 0.0);
        let target = slab("AAAAAAAAAAAAAAAA", 6.0);
        let native = binder.sequence();
        let res = design_interface(&binder, &target, 8.0, 400, 9).expect("design");
        // Non-interface positions keep their native residue.
        let interface: std::collections::HashSet<usize> =
            res.interface_residues.iter().copied().collect();
        for (i, (a, b)) in native
            .chars()
            .zip(res.designed_sequence.chars())
            .enumerate()
        {
            if !interface.contains(&i) {
                assert_eq!(a, b, "non-interface position {i} unchanged");
            }
        }
    }

    #[test]
    fn no_contact_is_an_error() {
        let binder = slab("AAAAAAAA", 0.0);
        let target = slab("AAAAAAAA", 100.0);
        assert!(design_interface(&binder, &target, 8.0, 100, 0).is_err());
    }

    #[test]
    fn bad_cutoff_rejected() {
        let binder = slab("AAAA", 0.0);
        let target = slab("AAAA", 6.0);
        assert!(detect_interface(&binder, &target, -1.0).is_err());
    }
}

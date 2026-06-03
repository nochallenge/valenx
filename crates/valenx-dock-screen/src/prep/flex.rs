//! Feature 4 — flexible-sidechain selection on the receptor.
//!
//! AutoDock Vina and AutoDock 4 can let a small number of receptor
//! side chains move during docking — they are split out of the rigid
//! receptor into their own torsion tree (the `flex` PDBQT input). The
//! hard part is *which* side chains to flex: too many and the search
//! space explodes; too few and the induced-fit is missed.
//!
//! This module makes that selection. Given a receptor [`Structure`]
//! and a region of interest (a [`GridBox`] or a set of ligand-atom
//! positions), [`select_flexible`] returns the residues whose side
//! chains have at least one heavy atom within a cutoff of the region.
//! Glycine and alanine are excluded (no rotatable side chain); proline
//! is excluded (its ring-locked side chain does not flex like the
//! others).

use nalgebra::Vector3;
use valenx_biostruct::structure::{Atom, Structure};

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;

/// A receptor residue selected to flex during docking.
#[derive(Clone, Debug, PartialEq)]
pub struct FlexibleSidechain {
    /// Chain identifier the residue lives on.
    pub chain_id: String,
    /// Residue name (`"ARG"`, `"TYR"`, …).
    pub residue_name: String,
    /// Residue sequence number.
    pub seq_num: i32,
    /// Number of rotatable χ angles the side chain has (informational —
    /// see [`chi_count`]).
    pub chi_count: usize,
    /// Closest distance (Å) from any side-chain heavy atom to the
    /// region of interest.
    pub min_distance: f64,
}

/// The result of a flexible-sidechain selection — the chosen residues
/// plus the cutoff that produced them.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct FlexSelection {
    /// Residues selected to flex, sorted nearest-region first.
    pub residues: Vec<FlexibleSidechain>,
    /// The distance cutoff (Å) used for the selection.
    pub cutoff: f64,
}

impl FlexSelection {
    /// Number of residues selected.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// `true` if nothing was selected.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Total rotatable χ angles across every selected side chain — the
    /// extra torsional degrees of freedom flexible docking adds.
    pub fn total_chi(&self) -> usize {
        self.residues.iter().map(|r| r.chi_count).sum()
    }
}

/// Backbone atom names — these are *never* part of a flexible side
/// chain.
const BACKBONE: &[&str] = &["N", "CA", "C", "O", "OXT", "H", "HA"];

/// Residues with no flexible side chain (or a ring-locked one).
fn is_rigid_residue(name: &str) -> bool {
    matches!(name, "GLY" | "ALA" | "PRO")
}

/// The number of rotatable χ angles a standard amino-acid side chain
/// has. Used to score how much flexibility a residue contributes.
/// Returns `0` for the rigid residues and unrecognised names.
pub fn chi_count(residue_name: &str) -> usize {
    match residue_name {
        "ARG" => 4,
        "LYS" => 4,
        "MET" => 3,
        "GLU" | "GLN" => 3,
        "ASP" | "ASN" | "HIS" | "PHE" | "TYR" | "TRP" | "LEU" => 2,
        "ILE" => 2,
        "SER" | "THR" | "CYS" | "VAL" => 1,
        // GLY / ALA / PRO and anything unrecognised.
        _ => 0,
    }
}

/// Side-chain heavy atoms of a residue — every non-hydrogen atom whose
/// name is not a backbone atom.
fn sidechain_atoms(residue: &valenx_biostruct::structure::Residue) -> Vec<&Atom> {
    residue
        .atoms
        .iter()
        .filter(|a| !a.is_hydrogen() && !BACKBONE.contains(&a.name.as_str()))
        .collect()
}

/// Select receptor side chains to flex, given a set of points that
/// define the region of interest (typically ligand-atom positions or
/// pocket-lining atoms).
///
/// A residue is selected when it is a non-rigid amino acid with at
/// least one side-chain heavy atom within `cutoff` Å of any region
/// point. Results are sorted nearest-first.
///
/// Returns [`DockScreenError::Invalid`] if `cutoff` is non-positive,
/// or [`DockScreenError::InvalidReceptor`] if the structure has no
/// amino-acid residues.
pub fn select_flexible(
    receptor: &Structure,
    region: &[Vector3<f64>],
    cutoff: f64,
) -> Result<FlexSelection> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(DockScreenError::invalid(
            "cutoff",
            format!("flexible-sidechain cutoff must be positive, got {cutoff}"),
        ));
    }
    if region.is_empty() {
        return Err(DockScreenError::invalid(
            "region",
            "region of interest has no points",
        ));
    }
    let model = receptor.first_model();
    let mut any_amino = false;
    let mut selected: Vec<FlexibleSidechain> = Vec::new();
    let cutoff_sq = cutoff * cutoff;

    for chain in &model.chains {
        for residue in &chain.residues {
            if !residue.is_amino_acid() {
                continue;
            }
            any_amino = true;
            if is_rigid_residue(&residue.name) {
                continue;
            }
            let sc = sidechain_atoms(residue);
            if sc.is_empty() {
                continue;
            }
            // Closest distance from any side-chain heavy atom to the
            // region.
            let mut min_d2 = f64::INFINITY;
            for atom in &sc {
                let p = Vector3::new(atom.coord.x, atom.coord.y, atom.coord.z);
                for r in region {
                    let d2 = (p - r).norm_squared();
                    if d2 < min_d2 {
                        min_d2 = d2;
                    }
                }
            }
            if min_d2 <= cutoff_sq {
                selected.push(FlexibleSidechain {
                    chain_id: chain.id.clone(),
                    residue_name: residue.name.clone(),
                    seq_num: residue.seq_num,
                    chi_count: chi_count(&residue.name),
                    min_distance: min_d2.sqrt(),
                });
            }
        }
    }
    if !any_amino {
        return Err(DockScreenError::invalid_receptor(
            "structure has no amino-acid residues to flex",
        ));
    }
    selected.sort_by(|a, b| {
        a.min_distance
            .partial_cmp(&b.min_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(FlexSelection {
        residues: selected,
        cutoff,
    })
}

/// Convenience: select flexible side chains using a [`GridBox`]'s
/// eight corners plus its centre as the region of interest. A coarse
/// but cheap "flex everything lining the box" heuristic.
pub fn select_flexible_in_box(
    receptor: &Structure,
    grid: &GridBox,
    cutoff: f64,
) -> Result<FlexSelection> {
    let lo = grid.origin();
    let hi = grid.max_corner();
    let region = vec![
        grid.center,
        Vector3::new(lo.x, lo.y, lo.z),
        Vector3::new(hi.x, lo.y, lo.z),
        Vector3::new(lo.x, hi.y, lo.z),
        Vector3::new(lo.x, lo.y, hi.z),
        Vector3::new(hi.x, hi.y, lo.z),
        Vector3::new(hi.x, lo.y, hi.z),
        Vector3::new(lo.x, hi.y, hi.z),
        Vector3::new(hi.x, hi.y, hi.z),
    ];
    select_flexible(receptor, &region, cutoff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_biostruct::load_structure;

    /// A tiny two-residue receptor: an ALA (rigid) and an ARG (whose
    /// side chain reaches toward the origin).
    const RECEPTOR: &str = "\
ATOM      1  N   ALA A   1      20.000  20.000  20.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      21.000  20.000  20.000  1.00  0.00           C
ATOM      3  C   ALA A   1      22.000  20.000  20.000  1.00  0.00           C
ATOM      4  O   ALA A   1      23.000  20.000  20.000  1.00  0.00           O
ATOM      5  CB  ALA A   1      21.000  21.000  20.000  1.00  0.00           C
ATOM      6  N   ARG A   2       3.000   0.000   0.000  1.00  0.00           N
ATOM      7  CA  ARG A   2       4.000   0.000   0.000  1.00  0.00           C
ATOM      8  C   ARG A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      9  O   ARG A   2       6.000   0.000   0.000  1.00  0.00           O
ATOM     10  CB  ARG A   2       3.500   1.000   0.000  1.00  0.00           C
ATOM     11  CG  ARG A   2       2.500   1.500   0.000  1.00  0.00           C
ATOM     12  CD  ARG A   2       1.500   1.000   0.000  1.00  0.00           C
ATOM     13  NE  ARG A   2       0.800   0.500   0.000  1.00  0.00           N
ATOM     14  CZ  ARG A   2       0.300   0.200   0.000  1.00  0.00           C
END
";

    #[test]
    fn chi_count_known_values() {
        assert_eq!(chi_count("ARG"), 4);
        assert_eq!(chi_count("SER"), 1);
        assert_eq!(chi_count("GLY"), 0);
        assert_eq!(chi_count("ALA"), 0);
        assert_eq!(chi_count("PRO"), 0);
        assert_eq!(chi_count("UNKNOWN"), 0);
    }

    #[test]
    fn selects_arg_near_origin_but_not_distant_ala() {
        let s = load_structure(RECEPTOR, "x").unwrap();
        // Region = a single point at the origin. ARG's side chain
        // (CZ at 0.3,0.2,0) is within 5 Å; ALA is 35 Å away.
        let region = vec![Vector3::new(0.0, 0.0, 0.0)];
        let sel = select_flexible(&s, &region, 5.0).unwrap();
        assert_eq!(sel.len(), 1);
        assert_eq!(sel.residues[0].residue_name, "ARG");
        assert_eq!(sel.residues[0].chi_count, 4);
        assert_eq!(sel.total_chi(), 4);
    }

    #[test]
    fn ala_is_never_selected_even_when_close() {
        // Put the region right on the ALA — it must still be excluded
        // (alanine has no rotatable side chain).
        let s = load_structure(RECEPTOR, "x").unwrap();
        let region = vec![Vector3::new(21.0, 21.0, 20.0)];
        let sel = select_flexible(&s, &region, 5.0).unwrap();
        assert!(sel.is_empty(), "ALA must not be selected");
    }

    #[test]
    fn rejects_bad_cutoff_and_empty_region() {
        let s = load_structure(RECEPTOR, "x").unwrap();
        assert!(select_flexible(&s, &[Vector3::zeros()], 0.0).is_err());
        assert!(select_flexible(&s, &[Vector3::zeros()], -1.0).is_err());
        assert!(select_flexible(&s, &[], 5.0).is_err());
    }

    #[test]
    fn select_in_box_uses_corner_region() {
        let s = load_structure(RECEPTOR, "x").unwrap();
        // A box centred near the origin reaches the ARG side chain.
        let gb = GridBox::cubic([0.0, 0.0, 0.0], 8.0).unwrap();
        let sel = select_flexible_in_box(&s, &gb, 3.0).unwrap();
        assert_eq!(sel.len(), 1);
        assert_eq!(sel.residues[0].residue_name, "ARG");
    }

    #[test]
    fn distant_box_selects_nothing() {
        let s = load_structure(RECEPTOR, "x").unwrap();
        // A box far from both residues.
        let gb = GridBox::cubic([100.0, 100.0, 100.0], 8.0).unwrap();
        let sel = select_flexible_in_box(&s, &gb, 3.0).unwrap();
        assert!(sel.is_empty());
    }
}

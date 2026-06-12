//! Structure editing and analysis: B-factor statistics, interface /
//! contact-residue detection, and point-mutation modelling.
//!
//! - [`bfactor_analysis`] — per-residue mean B-factors and the
//!   normalised (Z-scored) B-factors a flexibility map wants.
//! - [`interface_residues`] — the residues of two chains that sit at
//!   their mutual interface.
//! - [`mutate_residue`] — replace a residue's sidechain with that of
//!   another amino acid (point-mutation modelling v1).

use crate::error::{BiostructError, Result};
use crate::geometry::distance::NeighborGrid;
use crate::structure::{Atom, Chain, Model, Residue};
use nalgebra::{Point3, Vector3};

// --- B-factor analysis ----------------------------------------------

/// Per-residue B-factor statistics.
#[derive(Clone, Debug, PartialEq)]
pub struct BFactorAnalysis {
    /// Mean atomic B-factor of every residue, in model residue order.
    pub per_residue_mean: Vec<f64>,
    /// Z-score-normalised per-residue B-factor: `(B − μ) / σ` over the
    /// whole model. Parallel to `per_residue_mean`.
    pub per_residue_normalized: Vec<f64>,
    /// Model-wide mean atomic B-factor.
    pub mean: f64,
    /// Model-wide B-factor standard deviation.
    pub std_dev: f64,
}

impl BFactorAnalysis {
    /// Residue indices whose normalised B-factor exceeds `z` — the
    /// most flexible / least-ordered residues.
    pub fn flexible_residues(&self, z: f64) -> Vec<usize> {
        self.per_residue_normalized
            .iter()
            .enumerate()
            .filter(|(_, b)| **b > z)
            .map(|(i, _)| i)
            .collect()
    }
}

/// Compute per-residue B-factor statistics for a model.
pub fn bfactor_analysis(model: &Model) -> Result<BFactorAnalysis> {
    let mut all: Vec<f64> = Vec::new();
    let mut per_residue_mean: Vec<f64> = Vec::new();

    for chain in &model.chains {
        for residue in &chain.residues {
            if residue.atoms.is_empty() {
                per_residue_mean.push(0.0);
                continue;
            }
            let mut sum = 0.0;
            for a in &residue.atoms {
                sum += a.b_factor;
                all.push(a.b_factor);
            }
            per_residue_mean.push(sum / residue.atoms.len() as f64);
        }
    }

    if all.is_empty() {
        return Err(BiostructError::invalid("model", "model has no atoms"));
    }
    let mean = all.iter().sum::<f64>() / all.len() as f64;
    let variance = all.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / all.len() as f64;
    let std_dev = variance.sqrt();

    let per_residue_normalized: Vec<f64> = per_residue_mean
        .iter()
        .map(|b| {
            if std_dev > 1e-9 {
                (b - mean) / std_dev
            } else {
                0.0
            }
        })
        .collect();

    Ok(BFactorAnalysis {
        per_residue_mean,
        per_residue_normalized,
        mean,
        std_dev,
    })
}

// --- interface / contact residues -----------------------------------

/// Detect the residues at the interface between two chains.
///
/// A residue of `chain_a` is an interface residue when any of its
/// atoms comes within `cutoff` ångström of any atom of `chain_b`,
/// and vice versa. Returns `(residues_of_a, residues_of_b)` as
/// residue indices into each chain.
pub fn interface_residues(
    chain_a: &Chain,
    chain_b: &Chain,
    cutoff: f64,
) -> Result<(Vec<usize>, Vec<usize>)> {
    if cutoff <= 0.0 || cutoff.is_nan() {
        return Err(BiostructError::invalid("cutoff", "must be positive"));
    }

    // Build a grid over chain B's atoms, remembering which residue of
    // B each grid point belongs to.
    let mut b_points: Vec<Point3<f64>> = Vec::new();
    let mut b_residue_of_point: Vec<usize> = Vec::new();
    for (ri, residue) in chain_b.residues.iter().enumerate() {
        for atom in &residue.atoms {
            b_points.push(atom.coord);
            b_residue_of_point.push(ri);
        }
    }
    if b_points.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let grid = NeighborGrid::new(&b_points, cutoff.max(0.01))?;

    let mut a_iface = Vec::new();
    let mut b_iface_set: Vec<bool> = vec![false; chain_b.residues.len()];

    for (ai, residue) in chain_a.residues.iter().enumerate() {
        let mut is_interface = false;
        for atom in &residue.atoms {
            let hits = grid.within(&atom.coord, cutoff);
            if !hits.is_empty() {
                is_interface = true;
                for h in hits {
                    b_iface_set[b_residue_of_point[h]] = true;
                }
            }
        }
        if is_interface {
            a_iface.push(ai);
        }
    }

    let b_iface: Vec<usize> = b_iface_set
        .iter()
        .enumerate()
        .filter(|(_, v)| **v)
        .map(|(i, _)| i)
        .collect();
    Ok((a_iface, b_iface))
}

// --- point-mutation modelling ---------------------------------------

/// The internal-coordinate sidechain template of an amino acid: the
/// atoms beyond Cβ, each placed relative to the local backbone.
///
/// Each entry is `(atom_name, element, offset_from_cbeta)`, where the
/// offset is in a local frame whose x-axis is `Cβ→Cα`. v1 templates
/// are coarse — geometrically reasonable, ideal-rotamer placements.
fn sidechain_template(resname: &str) -> Option<&'static [(&'static str, &'static str, [f64; 3])]> {
    match resname {
        "ALA" => Some(&[]), // no atoms beyond CB
        "SER" => Some(&[("OG", "O", [1.4, 0.3, 0.0])]),
        "CYS" => Some(&[("SG", "S", [1.8, 0.4, 0.0])]),
        "THR" => Some(&[
            ("OG1", "O", [1.4, 0.3, 0.8]),
            ("CG2", "C", [1.4, 0.3, -0.8]),
        ]),
        "VAL" => Some(&[
            ("CG1", "C", [1.5, 0.4, 0.9]),
            ("CG2", "C", [1.5, 0.4, -0.9]),
        ]),
        "LEU" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD1", "C", [2.8, 1.0, 0.9]),
            ("CD2", "C", [2.8, 1.0, -0.9]),
        ]),
        "ILE" => Some(&[
            ("CG1", "C", [1.5, 0.5, 0.9]),
            ("CG2", "C", [1.5, 0.5, -0.9]),
            ("CD1", "C", [2.9, 1.0, 0.9]),
        ]),
        "ASN" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("OD1", "O", [2.6, 1.1, 0.6]),
            ("ND2", "N", [2.6, 1.1, -0.6]),
        ]),
        "ASP" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("OD1", "O", [2.6, 1.1, 0.6]),
            ("OD2", "O", [2.6, 1.1, -0.6]),
        ]),
        "GLN" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD", "C", [2.9, 1.0, 0.0]),
            ("OE1", "O", [4.0, 1.6, 0.6]),
            ("NE2", "N", [4.0, 1.6, -0.6]),
        ]),
        "GLU" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD", "C", [2.9, 1.0, 0.0]),
            ("OE1", "O", [4.0, 1.6, 0.6]),
            ("OE2", "O", [4.0, 1.6, -0.6]),
        ]),
        "PHE" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD1", "C", [2.8, 1.0, 0.7]),
            ("CD2", "C", [2.8, 1.0, -0.7]),
            ("CE1", "C", [4.1, 1.6, 0.7]),
            ("CE2", "C", [4.1, 1.6, -0.7]),
            ("CZ", "C", [4.8, 1.9, 0.0]),
        ]),
        "TYR" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD1", "C", [2.8, 1.0, 0.7]),
            ("CD2", "C", [2.8, 1.0, -0.7]),
            ("CE1", "C", [4.1, 1.6, 0.7]),
            ("CE2", "C", [4.1, 1.6, -0.7]),
            ("CZ", "C", [4.8, 1.9, 0.0]),
            ("OH", "O", [6.2, 2.5, 0.0]),
        ]),
        "TRP" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD1", "C", [2.6, 1.1, 0.8]),
            ("CD2", "C", [2.6, 1.1, -0.8]),
            ("NE1", "N", [3.9, 1.6, 0.8]),
            ("CE2", "C", [3.9, 1.6, -0.8]),
            ("CE3", "C", [3.0, 1.4, -2.1]),
            ("CZ2", "C", [5.2, 2.1, -0.8]),
        ]),
        "HIS" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("ND1", "N", [2.6, 1.1, 0.7]),
            ("CD2", "C", [2.6, 1.1, -0.7]),
            ("CE1", "C", [3.9, 1.6, 0.4]),
            ("NE2", "N", [3.9, 1.6, -0.5]),
        ]),
        "LYS" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD", "C", [2.9, 1.0, 0.0]),
            ("CE", "C", [4.3, 1.6, 0.0]),
            ("NZ", "N", [5.7, 2.1, 0.0]),
        ]),
        "ARG" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("CD", "C", [2.9, 1.0, 0.0]),
            ("NE", "N", [4.3, 1.6, 0.0]),
            ("CZ", "C", [5.5, 2.1, 0.0]),
            ("NH1", "N", [6.7, 2.6, 0.7]),
            ("NH2", "N", [6.7, 2.6, -0.7]),
        ]),
        "MET" => Some(&[
            ("CG", "C", [1.5, 0.5, 0.0]),
            ("SD", "S", [2.9, 1.1, 0.0]),
            ("CE", "C", [4.4, 1.7, 0.0]),
        ]),
        "PRO" => Some(&[("CG", "C", [1.5, 0.5, 0.0]), ("CD", "C", [2.4, 1.3, 0.0])]),
        "GLY" => Some(&[]), // GLY has no CB; handled specially
        _ => None,
    }
}

/// Point-mutate the residue at `(chain_index, residue_index)` of a
/// model to the amino acid `target` (a three-letter code).
///
/// The backbone atoms (`N`, `CA`, `C`, `O`) are kept; the existing
/// sidechain is removed and the target's sidechain template is placed
/// relative to the kept backbone. The residue is renamed.
///
/// This is a real v1 sidechain-replacement modeller: it produces a
/// plausible all-atom mutant geometry. It is **not** a rotamer-search
/// / energy-minimisation modeller — the placed sidechain uses one
/// idealised template conformation and is not relaxed against the
/// surrounding structure.
pub fn mutate_residue(
    model: &mut Model,
    chain_index: usize,
    residue_index: usize,
    target: &str,
) -> Result<()> {
    let target = target.trim().to_ascii_uppercase();
    let template = sidechain_template(&target).ok_or_else(|| {
        BiostructError::invalid(
            "target",
            format!("`{target}` is not a known amino acid for mutation"),
        )
    })?;

    let chain = model
        .chains
        .get_mut(chain_index)
        .ok_or_else(|| BiostructError::invalid("chain_index", "out of range"))?;
    let residue = chain
        .residues
        .get_mut(residue_index)
        .ok_or_else(|| BiostructError::invalid("residue_index", "out of range"))?;
    if !residue.is_amino_acid() {
        return Err(BiostructError::invalid(
            "residue",
            "can only mutate an amino-acid residue",
        ));
    }

    // Backbone atoms we keep.
    let n = backbone_coord(residue, "N")?;
    let ca = backbone_coord(residue, "CA")?;
    let c = backbone_coord(residue, "C")?;

    // A local frame at CB: x toward CA→CB, y in the N-CA-C plane.
    // CB itself is rebuilt from the backbone for consistency.
    let cb = build_cbeta(&n, &ca, &c);
    let x_axis = (cb - ca).normalize();
    // y_axis: component of (C - CA) perpendicular to x.
    let mut y_axis: Vector3<f64> = c - ca;
    y_axis -= x_axis * y_axis.dot(&x_axis);
    let y_axis = if y_axis.norm() > 1e-6 {
        y_axis.normalize()
    } else {
        x_axis.cross(&Vector3::z()).normalize()
    };
    let z_axis = x_axis.cross(&y_axis).normalize();

    // Keep only backbone atoms; preserve their occupancy / B-factor.
    let kept: Vec<Atom> = residue
        .atoms
        .iter()
        .filter(|a| matches!(a.name.as_str(), "N" | "CA" | "C" | "O"))
        .cloned()
        .collect();
    // A reference B-factor / occupancy for the new sidechain atoms.
    let ref_b = residue.ca().map(|a| a.b_factor).unwrap_or(0.0);
    let ref_occ = residue.ca().map(|a| a.occupancy).unwrap_or(1.0);

    let mut new_atoms = kept;
    // Glycine has no CB; every other residue gets a CB.
    if target != "GLY" {
        new_atoms.push(Atom {
            serial: 0,
            name: "CB".to_string(),
            alt_loc: ' ',
            element: "C".to_string(),
            coord: cb,
            occupancy: ref_occ,
            b_factor: ref_b,
            charge: 0,
        });
    }
    // Place each template atom: offset is (along x, along y, along z)
    // from CB.
    for (name, element, off) in template {
        let pos = cb.coords + x_axis * off[0] + y_axis * off[1] + z_axis * off[2];
        new_atoms.push(Atom {
            serial: 0,
            name: name.to_string(),
            alt_loc: ' ',
            element: element.to_string(),
            coord: Point3::from(pos),
            occupancy: ref_occ,
            b_factor: ref_b,
            charge: 0,
        });
    }

    residue.atoms = new_atoms;
    residue.name = target;
    Ok(())
}

/// Highest-occupancy backbone-atom coordinate, or an error.
fn backbone_coord(residue: &Residue, name: &str) -> Result<Point3<f64>> {
    residue.primary_atom(name).map(|a| a.coord).ok_or_else(|| {
        BiostructError::invalid(
            "residue",
            format!("residue lacks the backbone atom `{name}` needed for mutation"),
        )
    })
}

/// Build an idealised Cβ position from the backbone `N`, `CA`, `C`.
///
/// Cβ sits ~1.53 Å from Cα, tetrahedral to the `N–CA` and `C–CA`
/// bonds. The standard construction reflects the average of the two
/// backbone-bond directions through Cα and adds the out-of-plane
/// component.
pub fn build_cbeta(n: &Point3<f64>, ca: &Point3<f64>, c: &Point3<f64>) -> Point3<f64> {
    let to_n = (n - ca).normalize();
    let to_c = (c - ca).normalize();
    // In-plane bisector pointing away from both backbone bonds.
    let bisector = -(to_n + to_c);
    let bis = if bisector.norm() > 1e-6 {
        bisector.normalize()
    } else {
        Vector3::z()
    };
    // Out-of-plane direction.
    let perp = to_c.cross(&to_n);
    let perp = if perp.norm() > 1e-6 {
        perp.normalize()
    } else {
        Vector3::x()
    };
    // Mix in-plane and out-of-plane for a tetrahedral CB.
    let dir = (bis * 0.5 + perp * 0.866).normalize();
    Point3::from(ca.coords + dir * 1.53)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protein_model() -> Model {
        let mut chain = Chain::new("A");
        for seq in 1..=4 {
            let mut r = Residue::new("ALA", seq);
            let base = (seq as f64) * 3.8;
            let mut n = Atom::new("N", "N", Point3::new(base, 0.0, 0.0));
            n.b_factor = 10.0 + seq as f64;
            let mut ca = Atom::new("CA", "C", Point3::new(base + 1.5, 0.0, 0.0));
            ca.b_factor = 12.0 + seq as f64;
            let mut c = Atom::new("C", "C", Point3::new(base + 2.5, 1.4, 0.0));
            c.b_factor = 11.0 + seq as f64;
            let mut o = Atom::new("O", "O", Point3::new(base + 2.5, 2.6, 0.0));
            o.b_factor = 13.0 + seq as f64;
            r.atoms.extend([n, ca, c, o]);
            chain.residues.push(r);
        }
        let mut model = Model::new(1);
        model.chains.push(chain);
        model
    }

    #[test]
    fn bfactor_means_and_normalisation() {
        let model = protein_model();
        let bf = bfactor_analysis(&model).unwrap();
        assert_eq!(bf.per_residue_mean.len(), 4);
        assert_eq!(bf.per_residue_normalized.len(), 4);
        // Normalised values are a Z-score: they straddle zero.
        let min = bf
            .per_residue_normalized
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max = bf
            .per_residue_normalized
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(min < 0.0 && max > 0.0);
    }

    #[test]
    fn flexible_residue_detection() {
        let model = protein_model();
        let bf = bfactor_analysis(&model).unwrap();
        // The last residue has the highest B-factors.
        let flex = bf.flexible_residues(0.5);
        assert!(flex.contains(&3));
    }

    #[test]
    fn interface_between_two_close_chains() {
        // chain A near the origin, chain B 4 A away in x: their first
        // residues should form the interface.
        let mut a = Chain::new("A");
        let mut ra = Residue::new("ALA", 1);
        ra.atoms
            .push(Atom::new("CA", "C", Point3::new(0.0, 0.0, 0.0)));
        a.residues.push(ra);
        let mut ra2 = Residue::new("ALA", 2);
        ra2.atoms
            .push(Atom::new("CA", "C", Point3::new(0.0, 0.0, 30.0)));
        a.residues.push(ra2);

        let mut b = Chain::new("B");
        let mut rb = Residue::new("ALA", 1);
        rb.atoms
            .push(Atom::new("CA", "C", Point3::new(4.0, 0.0, 0.0)));
        b.residues.push(rb);

        let (ai, bi) = interface_residues(&a, &b, 5.0).unwrap();
        assert_eq!(ai, vec![0]);
        assert_eq!(bi, vec![0]);
    }

    #[test]
    fn no_interface_when_chains_far_apart() {
        let mut a = Chain::new("A");
        let mut ra = Residue::new("ALA", 1);
        ra.atoms.push(Atom::new("CA", "C", Point3::origin()));
        a.residues.push(ra);
        let mut b = Chain::new("B");
        let mut rb = Residue::new("ALA", 1);
        rb.atoms
            .push(Atom::new("CA", "C", Point3::new(100.0, 0.0, 0.0)));
        b.residues.push(rb);
        let (ai, bi) = interface_residues(&a, &b, 5.0).unwrap();
        assert!(ai.is_empty() && bi.is_empty());
    }

    #[test]
    fn cbeta_is_near_idealised_distance() {
        let n = Point3::new(0.0, 0.0, 0.0);
        let ca = Point3::new(1.46, 0.0, 0.0);
        let c = Point3::new(2.0, 1.4, 0.0);
        let cb = build_cbeta(&n, &ca, &c);
        let d = (cb - ca).norm();
        assert!((d - 1.53).abs() < 1e-6, "CB-CA distance {d}");
    }

    #[test]
    fn mutate_ala_to_trp_adds_sidechain() {
        let mut model = protein_model();
        let before = model.chains[0].residues[1].atoms.len();
        assert_eq!(before, 4); // ALA backbone only in our test build
        mutate_residue(&mut model, 0, 1, "TRP").unwrap();
        let res = &model.chains[0].residues[1];
        assert_eq!(res.name, "TRP");
        // backbone (4) + CB + 7 template atoms.
        assert!(res.atoms.len() > before);
        assert!(res.atom("CB").is_some());
        assert!(res.atom("CG").is_some());
        assert!(res.atom("N").is_some()); // backbone preserved
    }

    #[test]
    fn mutate_to_glycine_removes_cbeta() {
        let mut model = protein_model();
        mutate_residue(&mut model, 0, 0, "GLY").unwrap();
        let res = &model.chains[0].residues[0];
        assert_eq!(res.name, "GLY");
        assert!(res.atom("CB").is_none());
        // backbone still there.
        assert_eq!(res.atoms.len(), 4);
    }

    #[test]
    fn mutate_rejects_unknown_target() {
        let mut model = protein_model();
        assert!(mutate_residue(&mut model, 0, 0, "XYZ").is_err());
    }

    #[test]
    fn mutate_rejects_bad_index() {
        let mut model = protein_model();
        assert!(mutate_residue(&mut model, 9, 0, "SER").is_err());
        assert!(mutate_residue(&mut model, 0, 99, "SER").is_err());
    }
}

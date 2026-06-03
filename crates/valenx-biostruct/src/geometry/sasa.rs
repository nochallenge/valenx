//! Solvent-accessible surface area via the Shrake-Rupley algorithm.
//!
//! Each atom is given an expanded radius `r_vdw + r_probe`. A set of
//! points is spread evenly over that sphere (a golden-spiral / Fibonacci
//! lattice). A test point is *accessible* if it lies outside the
//! expanded sphere of every other atom. The atom's SASA is the
//! accessible fraction times the sphere area `4π r²`.
//!
//! This is the same algorithm Biopython's `ShrakeRupley` and the
//! original 1973 Shrake-Rupley paper use. The number of sample points
//! trades accuracy against speed; 92-200 points is the usual range.

use crate::error::{BiostructError, Result};
use crate::geometry::distance::NeighborGrid;
use crate::geometry::vdw::vdw_radius;
use crate::structure::Model;
use nalgebra::{Point3, Vector3};

/// Default probe radius — a water molecule, 1.4 Å.
pub const WATER_PROBE: f64 = 1.4;

/// Per-atom solvent-accessible surface area.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomSasa {
    /// Index of the atom in [`Model::atoms`] order.
    pub atom_index: usize,
    /// Accessible surface area of this atom, ångström².
    pub area: f64,
}

/// A full SASA result for one model.
#[derive(Clone, Debug, PartialEq)]
pub struct SasaResult {
    /// Per-atom accessible areas, in atom order.
    pub atoms: Vec<AtomSasa>,
    /// Probe radius used.
    pub probe_radius: f64,
    /// Number of sphere sample points per atom.
    pub n_points: usize,
}

impl SasaResult {
    /// Total solvent-accessible surface area, ångström².
    pub fn total(&self) -> f64 {
        self.atoms.iter().map(|a| a.area).sum()
    }
}

/// Generate `n` near-uniformly-spread unit-sphere points using a
/// Fibonacci (golden-spiral) lattice.
pub fn sphere_points(n: usize) -> Vec<Vector3<f64>> {
    if n == 0 {
        return Vec::new();
    }
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt()); // ~2.399963
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        // z runs from near +1 to near -1.
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
        let r = (1.0 - z * z).max(0.0).sqrt();
        let theta = golden * i as f64;
        pts.push(Vector3::new(r * theta.cos(), r * theta.sin(), z));
    }
    pts
}

/// Compute per-atom SASA for a model.
///
/// `probe_radius` is the solvent probe (use [`WATER_PROBE`]).
/// `n_points` is the sphere sampling density (>= 1).
pub fn shrake_rupley(model: &Model, probe_radius: f64, n_points: usize) -> Result<SasaResult> {
    if probe_radius < 0.0 {
        return Err(BiostructError::invalid(
            "probe_radius",
            "probe radius must be non-negative",
        ));
    }
    if n_points == 0 {
        return Err(BiostructError::invalid(
            "n_points",
            "need at least one sphere sample point",
        ));
    }

    // Flatten atoms; expanded radius = vdw + probe.
    struct A {
        coord: Point3<f64>,
        radius: f64,
    }
    let atoms: Vec<A> = model
        .atoms()
        .map(|a| A {
            coord: a.coord,
            radius: vdw_radius(&a.element) + probe_radius,
        })
        .collect();
    let n = atoms.len();
    if n == 0 {
        return Ok(SasaResult {
            atoms: Vec::new(),
            probe_radius,
            n_points,
        });
    }

    let unit = sphere_points(n_points);
    // The largest reach of any two expanded spheres bounds the grid
    // cell: two atoms can only occlude each other within r_i + r_j.
    let max_r = atoms.iter().map(|a| a.radius).fold(0.0_f64, f64::max);
    let cell = 2.0 * max_r;
    let coords: Vec<Point3<f64>> = atoms.iter().map(|a| a.coord).collect();
    let grid = NeighborGrid::new(&coords, cell)?;

    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let ai = &atoms[i];
        // Candidate occluders: atoms whose expanded sphere can reach
        // atom i's expanded sphere.
        let neighbors: Vec<usize> = grid
            .within_excluding(&ai.coord, cell, i)
            .into_iter()
            .filter(|&j| {
                let d = (ai.coord - atoms[j].coord).norm();
                d < ai.radius + atoms[j].radius
            })
            .collect();

        let mut accessible = 0usize;
        'point: for u in &unit {
            // Sample point on atom i's expanded sphere.
            let test = ai.coord + u * ai.radius;
            for &j in &neighbors {
                let aj = &atoms[j];
                if (test - aj.coord).norm_squared() < aj.radius * aj.radius {
                    continue 'point; // buried by neighbour j
                }
            }
            accessible += 1;
        }
        let sphere_area = 4.0 * std::f64::consts::PI * ai.radius * ai.radius;
        let area = sphere_area * accessible as f64 / n_points as f64;
        out.push(AtomSasa {
            atom_index: i,
            area,
        });
    }

    Ok(SasaResult {
        atoms: out,
        probe_radius,
        n_points,
    })
}

/// Per-residue SASA totals, in `Model::residues` order, derived by
/// summing the per-atom areas of [`shrake_rupley`].
pub fn residue_sasa(model: &Model, probe_radius: f64, n_points: usize) -> Result<Vec<f64>> {
    let per_atom = shrake_rupley(model, probe_radius, n_points)?;
    let mut residue_totals = Vec::new();
    let mut atom_cursor = 0usize;
    for chain in &model.chains {
        for residue in &chain.residues {
            let mut sum = 0.0;
            for _ in &residue.atoms {
                sum += per_atom.atoms[atom_cursor].area;
                atom_cursor += 1;
            }
            residue_totals.push(sum);
        }
    }
    Ok(residue_totals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain, Residue};

    #[test]
    fn sphere_points_lie_on_unit_sphere() {
        for p in sphere_points(64) {
            assert!((p.norm() - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn isolated_atom_is_fully_exposed() {
        // One carbon, no neighbours: SASA ~= 4*pi*(1.7+1.4)^2.
        let mut chain = Chain::new("A");
        let mut r = Residue::new("LIG", 1);
        r.hetatm = true;
        r.atoms
            .push(Atom::new("C1", "C", Point3::origin()));
        chain.residues.push(r);
        let mut model = Model::new(1);
        model.chains.push(chain);

        let sasa = shrake_rupley(&model, WATER_PROBE, 200).unwrap();
        let expected = 4.0 * std::f64::consts::PI * (1.7 + 1.4) * (1.7 + 1.4);
        // Fibonacci sampling at 200 points is within a few percent.
        let rel = (sasa.total() - expected).abs() / expected;
        assert!(rel < 0.05, "isolated-atom SASA off by {rel}");
    }

    #[test]
    fn buried_atom_has_less_area_than_exposed() {
        // A central atom surrounded by 6 close neighbours should have
        // markedly less SASA than any of the outer atoms.
        let mut chain = Chain::new("A");
        let mut r = Residue::new("LIG", 1);
        r.hetatm = true;
        r.atoms.push(Atom::new("C0", "C", Point3::origin()));
        for (k, d) in [
            Vector3::new(2.5, 0.0, 0.0),
            Vector3::new(-2.5, 0.0, 0.0),
            Vector3::new(0.0, 2.5, 0.0),
            Vector3::new(0.0, -2.5, 0.0),
            Vector3::new(0.0, 0.0, 2.5),
            Vector3::new(0.0, 0.0, -2.5),
        ]
        .into_iter()
        .enumerate()
        {
            r.atoms
                .push(Atom::new(format!("C{}", k + 1), "C", Point3::from(d)));
        }
        chain.residues.push(r);
        let mut model = Model::new(1);
        model.chains.push(chain);

        let sasa = shrake_rupley(&model, WATER_PROBE, 200).unwrap();
        let central = sasa.atoms[0].area;
        let outer = sasa.atoms[1].area;
        assert!(
            central < outer,
            "central {central} should be less exposed than outer {outer}"
        );
    }

    #[test]
    fn residue_sasa_sums_to_total() {
        let mut chain = Chain::new("A");
        for s in 1..=3 {
            let mut r = Residue::new("ALA", s);
            r.atoms.push(Atom::new(
                "CA",
                "C",
                Point3::new(s as f64 * 4.0, 0.0, 0.0),
            ));
            chain.residues.push(r);
        }
        let mut model = Model::new(1);
        model.chains.push(chain);

        let per_res = residue_sasa(&model, WATER_PROBE, 96).unwrap();
        let total = shrake_rupley(&model, WATER_PROBE, 96).unwrap().total();
        let sum: f64 = per_res.iter().sum();
        assert!((sum - total).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_arguments() {
        let model = Model::new(1);
        assert!(shrake_rupley(&model, -1.0, 100).is_err());
        assert!(shrake_rupley(&model, 1.4, 0).is_err());
    }
}

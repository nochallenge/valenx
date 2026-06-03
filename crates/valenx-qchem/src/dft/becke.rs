//! Becke fuzzy-cell partitioning of a molecular grid.
//!
//! A molecular DFT grid is the union of atom-centred grids. Where two
//! atomic grids overlap, the integration weight must be shared so the
//! whole-molecule integral counts every region once. The **Becke 1988
//! scheme** (A. D. Becke, *J. Chem. Phys.* 88, 2547 (1988)) does this
//! with a smooth nuclear-weight function.
//!
//! ## The scheme
//!
//! For a grid point `r` belonging to atom `A`, Becke defines a cell
//! function `P_A(r)` built from the *elliptical coordinate*
//!
//! ```text
//! μ_AB = (|r − A| − |r − B|) / |A − B|
//! ```
//!
//! for every other atom `B`. A smooth step `s(μ)` — three iterations of
//! the polynomial `p(μ) = (3μ − μ³)/2` followed by `(1 − p)/2` —
//! switches `s` from 1 deep in cell `A` to 0 deep in cell `B`. The
//! unnormalised cell weight is the product `P_A = Π_{B≠A} s(μ_AB)` and
//! the **partition weight** of the point is
//!
//! ```text
//! w_A(r) = P_A(r) / Σ_C P_C(r).
//! ```
//!
//! ## Atomic-size adjustment
//!
//! Becke's appendix adds an atomic-size correction so a grid point near
//! a small atom is not unfairly stolen by a large neighbour: the plain
//! `μ_AB` is shifted to `ν_AB = μ_AB + a_AB(1 − μ_AB²)` where `a_AB` is
//! built from the ratio of Bragg-Slater radii. This module includes
//! that correction, which markedly improves the integration of
//! hetero-atomic molecules.

/// Bragg-Slater atomic radii in bohr, indexed by `Z = 1..18`.
///
/// Used only for the Becke atomic-size adjustment. `Z` outside the
/// range falls back to the carbon radius.
fn bragg_slater_radius(z: u8) -> f64 {
    // Bragg-Slater radii (Å) → bohr; H specially set to 0.35 Å per
    // Becke's recommendation.
    const R_ANG: [f64; 18] = [
        0.35, 1.40, // H (Becke value), He
        1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 1.50, // Li..Ne
        1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.80, // Na..Ar
    ];
    let ang = if z >= 1 && (z as usize) <= R_ANG.len() {
        R_ANG[(z - 1) as usize]
    } else {
        0.70
    };
    ang * crate::geometry::BOHR_PER_ANGSTROM
}

/// The Becke smooth-step polynomial `p(μ) = (3μ − μ³) / 2`.
#[inline]
fn becke_p(mu: f64) -> f64 {
    0.5 * (3.0 * mu - mu * mu * mu)
}

/// The Becke cell-boundary switch `s(μ)` — three iterations of
/// [`becke_p`] then `(1 − p)/2`. Goes smoothly from 1 at `μ = −1` to 0
/// at `μ = +1`.
fn becke_switch(mu: f64) -> f64 {
    let f = becke_p(becke_p(becke_p(mu)));
    0.5 * (1.0 - f)
}

/// Becke's atomic-size adjustment: shift the elliptical coordinate so
/// the cell boundary sits in proportion to the two Bragg-Slater radii.
///
/// `chi = R_A / R_B`; the shift parameter `a_AB` is clamped to
/// `[−0.5, 0.5]` as Becke's appendix prescribes.
fn size_adjusted_mu(mu: f64, r_a: f64, r_b: f64) -> f64 {
    let chi = r_a / r_b;
    let u = (chi - 1.0) / (chi + 1.0);
    let mut a = u / (u * u - 1.0);
    a = a.clamp(-0.5, 0.5);
    let nu = mu + a * (1.0 - mu * mu);
    nu.clamp(-1.0, 1.0)
}

/// Compute the Becke partition weight of a grid point.
///
/// `point` is the Cartesian position of the grid point; `owner` is the
/// index (into `atom_positions` / `atom_z`) of the atom whose
/// atom-centred grid the point came from. The returned weight is in
/// `[0, 1]`; for a single atom it is exactly 1.
///
/// `atom_positions` and `atom_z` must have the same length — one entry
/// per atom.
pub fn becke_weight(
    point: [f64; 3],
    owner: usize,
    atom_positions: &[[f64; 3]],
    atom_z: &[u8],
) -> f64 {
    let n = atom_positions.len();
    if n <= 1 {
        return 1.0;
    }
    // Distance from the point to every atom.
    let dist: Vec<f64> = atom_positions
        .iter()
        .map(|a| {
            let dx = point[0] - a[0];
            let dy = point[1] - a[1];
            let dz = point[2] - a[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        })
        .collect();

    // Unnormalised cell weight P_C for every atom C.
    let mut cell = vec![1.0f64; n];
    for a in 0..n {
        for b in 0..n {
            if a == b {
                continue;
            }
            let dx = atom_positions[a][0] - atom_positions[b][0];
            let dy = atom_positions[a][1] - atom_positions[b][1];
            let dz = atom_positions[a][2] - atom_positions[b][2];
            let r_ab = (dx * dx + dy * dy + dz * dz).sqrt();
            if r_ab < 1.0e-12 {
                continue;
            }
            let mu = (dist[a] - dist[b]) / r_ab;
            let nu = size_adjusted_mu(
                mu,
                bragg_slater_radius(atom_z[a]),
                bragg_slater_radius(atom_z[b]),
            );
            cell[a] *= becke_switch(nu);
        }
    }

    let total: f64 = cell.iter().sum();
    if total <= 0.0 {
        // Numerically degenerate (point exactly on a cell boundary far
        // from every nucleus): share equally.
        1.0 / n as f64
    } else {
        cell[owner] / total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_atom_weight_is_one() {
        let w = becke_weight([1.0, 2.0, 3.0], 0, &[[0.0, 0.0, 0.0]], &[6]);
        assert!((w - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn switch_endpoints() {
        // s(-1) = 1, s(+1) = 0, s(0) = 1/2.
        assert!((becke_switch(-1.0) - 1.0).abs() < 1.0e-12);
        assert!(becke_switch(1.0).abs() < 1.0e-12);
        assert!((becke_switch(0.0) - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn switch_is_monotone_decreasing() {
        let mut prev = becke_switch(-1.0);
        let mut mu = -0.9;
        while mu <= 1.0 {
            let s = becke_switch(mu);
            assert!(s <= prev + 1.0e-12, "not monotone at μ={mu}");
            prev = s;
            mu += 0.1;
        }
    }

    #[test]
    fn partition_weights_sum_to_one() {
        // For any point, the partition weights over all atoms must
        // sum to exactly 1 (the partition of unity property).
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
            [1.5, 0.0, 1.0],
        ];
        let z = vec![8u8, 1, 1];
        for &probe in &[
            [0.3, 0.1, 0.5],
            [0.0, 0.0, 1.0],
            [-1.0, 2.0, 0.7],
            [0.7, 0.0, 1.6],
        ] {
            let sum: f64 = (0..positions.len())
                .map(|a| becke_weight(probe, a, &positions, &z))
                .sum();
            assert!(
                (sum - 1.0).abs() < 1.0e-12,
                "Σ weights = {sum} at {probe:?}"
            );
        }
    }

    #[test]
    fn point_at_a_nucleus_belongs_to_that_atom() {
        // A point sitting on atom 0 gets weight ≈ 1 for atom 0 and ≈ 0
        // for the others.
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]];
        let z = vec![6u8, 6];
        let w0 = becke_weight([0.0, 0.0, 0.001], 0, &positions, &z);
        let w1 = becke_weight([0.0, 0.0, 0.001], 1, &positions, &z);
        assert!(w0 > 0.99, "w0 = {w0}");
        assert!(w1 < 0.01, "w1 = {w1}");
    }

    #[test]
    fn midpoint_of_homonuclear_pair_is_shared_equally() {
        // Exactly between two identical atoms the weight is 1/2 each.
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let z = vec![6u8, 6];
        let w0 = becke_weight([0.0, 0.0, 1.0], 0, &positions, &z);
        let w1 = becke_weight([0.0, 0.0, 1.0], 1, &positions, &z);
        assert!((w0 - 0.5).abs() < 1.0e-12, "w0 = {w0}");
        assert!((w1 - 0.5).abs() < 1.0e-12, "w1 = {w1}");
    }

    #[test]
    fn size_adjustment_shifts_boundary_toward_larger_atom() {
        // For a hetero pair the cell boundary moves toward the larger
        // atom: at the geometric midpoint the smaller atom keeps less
        // than half. Use H (small) and a large fictitious neighbour.
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let z = vec![1u8, 11]; // H and Na (large Bragg-Slater radius)
        let w_h = becke_weight([0.0, 0.0, 1.0], 0, &positions, &z);
        // The small H atom should own less than half the midpoint.
        assert!(w_h < 0.5, "H midpoint weight = {w_h}");
    }
}

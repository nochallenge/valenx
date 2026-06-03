//! Cell list + Verlet neighbour list — **roadmap feature 10**.
//!
//! A naive nonbonded loop is `O(N²)`. This module brings it down to
//! `O(N)` with the two classic tricks:
//!
//! 1. **Cell list** — the simulation box is partitioned into a grid of
//!    cells, each at least `cutoff + skin` wide. Two atoms can only
//!    interact if they sit in the same cell or in one of its 26
//!    neighbours, so each atom is checked against only a constant
//!    number of candidates. [`cell_pairs`] builds the candidate list.
//!
//! 2. **Verlet list** — from the cell-list candidates a
//!    [`NeighborList`] keeps every pair within `cutoff + skin`. The
//!    extra `skin` means the list stays valid for several steps: it
//!    only needs rebuilding once an atom has drifted more than
//!    `skin/2`. [`NeighborList::needs_rebuild`] reports that.
//!
//! For a **non-periodic** box (or a box smaller than `3·(cutoff+skin)`
//! where a 3×3×3 cell grid would double-count) the code falls back to
//! the direct `O(N²)` pair enumeration — still correct, just not
//! accelerated. The fallback is transparent to the caller.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::pbc::SimBox;

/// A Verlet neighbour list: every atom pair within `cutoff + skin`.
#[derive(Clone, Debug, PartialEq)]
pub struct NeighborList {
    /// The `(i, j)` pairs, `i < j`.
    pairs: Vec<(usize, usize)>,
    /// Interaction cutoff the list was built for (nm).
    cutoff: f64,
    /// Verlet skin buffer (nm).
    skin: f64,
    /// Atom positions captured when the list was built — used by
    /// [`needs_rebuild`](Self::needs_rebuild).
    reference: Vec<Vector3<f64>>,
}

impl NeighborList {
    /// Builds the neighbour list for `positions` in `cell`.
    ///
    /// `cutoff` is the physical interaction range; `skin` is the extra
    /// buffer (a few tenths of `cutoff` is typical).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `cutoff` or `skin` is negative / not
    /// finite.
    pub fn build(
        positions: &[Vector3<f64>],
        cell: &SimBox,
        cutoff: f64,
        skin: f64,
    ) -> Result<Self> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(MdError::invalid("cutoff", "must be finite and positive"));
        }
        if !(skin.is_finite() && skin >= 0.0) {
            return Err(MdError::invalid("skin", "must be finite and non-negative"));
        }
        let range = cutoff + skin;
        let range_sq = range * range;
        let candidates = cell_pairs(positions, cell, range);
        let mut pairs = Vec::new();
        for (i, j) in candidates {
            let d = cell.min_image(positions[i] - positions[j]);
            if d.norm_squared() <= range_sq {
                pairs.push((i, j));
            }
        }
        Ok(NeighborList {
            pairs,
            cutoff,
            skin,
            reference: positions.to_vec(),
        })
    }

    /// The neighbour pairs, each `(i, j)` with `i < j`.
    pub fn pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }

    /// Number of neighbour pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// The cutoff the list was built for.
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// The skin buffer the list was built for.
    pub fn skin(&self) -> f64 {
        self.skin
    }

    /// Whether any atom has drifted more than `skin/2` from its
    /// reference position — the standard "rebuild the Verlet list"
    /// trigger. A drift of `skin/2` for two atoms moving toward each
    /// other closes the full skin.
    ///
    /// Returns `true` (rebuild) if the position count changed.
    pub fn needs_rebuild(&self, positions: &[Vector3<f64>], cell: &SimBox) -> bool {
        if positions.len() != self.reference.len() {
            return true;
        }
        let half = 0.5 * self.skin;
        let half_sq = half * half;
        positions
            .iter()
            .zip(&self.reference)
            .any(|(p, r)| cell.min_image(p - r).norm_squared() > half_sq)
    }
}

/// Enumerates candidate pairs within `range` using a cell list,
/// falling back to the direct `O(N²)` enumeration when a grid is not
/// applicable.
///
/// Returned pairs are `(i, j)` with `i < j` and are *candidates* only
/// — the caller still distance-tests them.
pub fn cell_pairs(
    positions: &[Vector3<f64>],
    cell: &SimBox,
    range: f64,
) -> Vec<(usize, usize)> {
    let n = positions.len();
    if n < 2 {
        return Vec::new();
    }
    // Decide whether a periodic cell grid is usable. We need at least
    // 3 cells per dimension so a 3x3x3 stencil does not wrap onto
    // itself; otherwise fall back to all-pairs.
    let use_grid = cell.is_periodic()
        && cell.is_orthorhombic()
        && range > 0.0
        && cell.edge_lengths().iter().all(|&l| l >= 3.0 * range);

    if !use_grid {
        return all_pairs(n);
    }

    let edges = cell.edge_lengths();
    let ncells: [usize; 3] = [
        ((edges[0] / range).floor() as usize).max(3),
        ((edges[1] / range).floor() as usize).max(3),
        ((edges[2] / range).floor() as usize).max(3),
    ];
    let total = ncells[0] * ncells[1] * ncells[2];
    let cell_of = |p: &Vector3<f64>| -> [usize; 3] {
        let wrapped = cell.wrap(*p);
        [
            (((wrapped.x / edges[0]) * ncells[0] as f64) as isize)
                .rem_euclid(ncells[0] as isize) as usize,
            (((wrapped.y / edges[1]) * ncells[1] as f64) as isize)
                .rem_euclid(ncells[1] as isize) as usize,
            (((wrapped.z / edges[2]) * ncells[2] as f64) as isize)
                .rem_euclid(ncells[2] as isize) as usize,
        ]
    };
    let flat = |c: [usize; 3]| c[0] + ncells[0] * (c[1] + ncells[1] * c[2]);

    // Bucket atoms.
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); total];
    for (idx, p) in positions.iter().enumerate() {
        buckets[flat(cell_of(p))].push(idx);
    }

    let mut pairs = Vec::new();
    // To avoid visiting each cell-pair twice, only pair a cell with a
    // canonical half of its neighbour stencil (the 13 "forward"
    // offsets plus the cell itself).
    let forward: [[isize; 3]; 13] = [
        [1, 0, 0],
        [-1, 1, 0],
        [0, 1, 0],
        [1, 1, 0],
        [-1, -1, 1],
        [0, -1, 1],
        [1, -1, 1],
        [-1, 0, 1],
        [0, 0, 1],
        [1, 0, 1],
        [-1, 1, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];
    for cx in 0..ncells[0] {
        for cy in 0..ncells[1] {
            for cz in 0..ncells[2] {
                let here = &buckets[flat([cx, cy, cz])];
                // Within-cell pairs.
                for a in 0..here.len() {
                    for b in (a + 1)..here.len() {
                        let (i, j) = (here[a], here[b]);
                        pairs.push((i.min(j), i.max(j)));
                    }
                }
                // Forward-neighbour-cell pairs.
                for off in &forward {
                    let nx =
                        (cx as isize + off[0]).rem_euclid(ncells[0] as isize) as usize;
                    let ny =
                        (cy as isize + off[1]).rem_euclid(ncells[1] as isize) as usize;
                    let nz =
                        (cz as isize + off[2]).rem_euclid(ncells[2] as isize) as usize;
                    let there = &buckets[flat([nx, ny, nz])];
                    for &i in here {
                        for &j in there {
                            if i != j {
                                pairs.push((i.min(j), i.max(j)));
                            }
                        }
                    }
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// The direct `O(N²)` pair enumeration.
fn all_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::with_capacity(n * (n.saturating_sub(1)) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            pairs.push((i, j));
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A simple cubic lattice of `m³` atoms, spacing `s`.
    fn cubic_lattice(m: usize, s: f64) -> Vec<Vector3<f64>> {
        let mut v = Vec::new();
        for i in 0..m {
            for j in 0..m {
                for k in 0..m {
                    v.push(Vector3::new(i as f64 * s, j as f64 * s, k as f64 * s));
                }
            }
        }
        v
    }

    #[test]
    fn all_pairs_count_is_n_choose_2() {
        assert_eq!(all_pairs(5).len(), 10);
        assert_eq!(all_pairs(1).len(), 0);
    }

    #[test]
    fn cell_list_finds_same_pairs_as_brute_force() {
        // 4x4x4 lattice, spacing 0.5 nm, in a 2 nm periodic box.
        let pos = cubic_lattice(4, 0.5);
        let cell = SimBox::cubic(2.0).unwrap();
        let cutoff = 0.6;
        let skin = 0.0;

        // Brute force reference.
        let mut brute = Vec::new();
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let d = cell.min_image(pos[i] - pos[j]);
                if d.norm() <= cutoff {
                    brute.push((i, j));
                }
            }
        }
        brute.sort_unstable();

        let nl = NeighborList::build(&pos, &cell, cutoff, skin).unwrap();
        let mut got = nl.pairs().to_vec();
        got.sort_unstable();
        assert_eq!(got, brute);
    }

    #[test]
    fn small_box_falls_back_to_all_pairs_correctly() {
        // Box too small for a 3-cell grid -> fallback path.
        let pos = cubic_lattice(2, 0.4);
        let cell = SimBox::cubic(0.9).unwrap();
        let cutoff = 0.5;
        let nl = NeighborList::build(&pos, &cell, cutoff, 0.0).unwrap();
        // Reference via min-image brute force.
        let mut brute = Vec::new();
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                if cell.min_image(pos[i] - pos[j]).norm() <= cutoff {
                    brute.push((i, j));
                }
            }
        }
        let mut got = nl.pairs().to_vec();
        got.sort_unstable();
        brute.sort_unstable();
        assert_eq!(got, brute);
    }

    #[test]
    fn needs_rebuild_triggers_on_drift() {
        let pos = cubic_lattice(3, 0.6);
        let cell = SimBox::cubic(3.0).unwrap();
        let nl = NeighborList::build(&pos, &cell, 0.7, 0.2).unwrap();
        // No movement -> no rebuild.
        assert!(!nl.needs_rebuild(&pos, &cell));
        // Move one atom past skin/2 = 0.1 nm.
        let mut moved = pos.clone();
        moved[0].x += 0.15;
        assert!(nl.needs_rebuild(&moved, &cell));
    }

    #[test]
    fn rejects_bad_parameters() {
        let pos = vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)];
        let cell = SimBox::cubic(5.0).unwrap();
        assert!(NeighborList::build(&pos, &cell, -1.0, 0.0).is_err());
        assert!(NeighborList::build(&pos, &cell, 1.0, -0.1).is_err());
    }

    #[test]
    fn non_periodic_box_uses_all_pairs() {
        let pos = cubic_lattice(3, 0.3);
        let cell = SimBox::none();
        let nl = NeighborList::build(&pos, &cell, 0.4, 0.0).unwrap();
        let mut brute = Vec::new();
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                if (pos[i] - pos[j]).norm() <= 0.4 {
                    brute.push((i, j));
                }
            }
        }
        let mut got = nl.pairs().to_vec();
        got.sort_unstable();
        brute.sort_unstable();
        assert_eq!(got, brute);
    }
}

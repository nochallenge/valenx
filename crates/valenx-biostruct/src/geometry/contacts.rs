//! Residue contact maps and steric-clash detection.
//!
//! - A **contact map** records which residue pairs come within a
//!   distance cutoff (the standard 8 Å Cβ–Cβ map, or any-atom
//!   minimum-distance contacts).
//! - **Clash detection** flags atom pairs whose centres are closer
//!   than the sum of their van der Waals radii minus a tolerance —
//!   the overlap a validation report wants.

use crate::error::{BiostructError, Result};
use crate::geometry::distance::NeighborGrid;
use crate::geometry::vdw::{is_bonded, vdw_radius};
use crate::structure::Model;

/// A residue-residue contact map for one model.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactMap {
    /// Number of residues the map covers (the polymer residues of the
    /// model, in order).
    pub n: usize,
    /// Contacting `(i, j)` residue-index pairs with `i < j`, plus the
    /// distance that satisfied the cutoff.
    pub contacts: Vec<(usize, usize, f64)>,
}

impl ContactMap {
    /// Whether residues `i` and `j` are in contact (order-independent).
    pub fn in_contact(&self, i: usize, j: usize) -> bool {
        let (lo, hi) = if i <= j { (i, j) } else { (j, i) };
        self.contacts.iter().any(|(a, b, _)| *a == lo && *b == hi)
    }

    /// The contact map as a dense symmetric boolean matrix
    /// (row-major, `n × n`). The diagonal is `false`.
    pub fn to_matrix(&self) -> Vec<Vec<bool>> {
        let mut m = vec![vec![false; self.n]; self.n];
        for (i, j, _) in &self.contacts {
            m[*i][*j] = true;
            m[*j][*i] = true;
        }
        m
    }

    /// Total number of contacting pairs.
    pub fn contact_count(&self) -> usize {
        self.contacts.len()
    }
}

/// Which atom of each residue the contact map measures from.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ContactMode {
    /// Cβ–Cβ distance (Cα for glycine). The classic contact-map
    /// definition.
    CbetaCbeta,
    /// Cα–Cα distance.
    CalphaCalpha,
    /// The minimum distance over *any* atom pair of the two residues
    /// — the most sensitive but slowest mode.
    MinAtom,
}

/// Build a residue contact map for a model's polymer residues.
///
/// `cutoff` is the contact distance in ångström; sequence-adjacent
/// residues closer than `min_seq_sep` apart in the chain are skipped
/// (set `min_seq_sep = 1` to include `i, i+1`, the usual choice is
/// `2` or `3` so trivially-close backbone neighbours do not dominate).
pub fn contact_map(
    model: &Model,
    cutoff: f64,
    mode: ContactMode,
    min_seq_sep: usize,
) -> Result<ContactMap> {
    if cutoff <= 0.0 || cutoff.is_nan() {
        return Err(BiostructError::invalid("cutoff", "must be positive"));
    }
    // Flatten polymer residues across all chains, keeping a chain tag
    // so the sequence-separation rule only applies within a chain.
    struct ResInfo<'a> {
        chain: usize,
        residue: &'a crate::structure::Residue,
    }
    let mut residues: Vec<ResInfo> = Vec::new();
    for (ci, chain) in model.chains.iter().enumerate() {
        for r in &chain.residues {
            if r.is_amino_acid() || r.is_nucleotide() {
                residues.push(ResInfo {
                    chain: ci,
                    residue: r,
                });
            }
        }
    }
    let n = residues.len();
    let mut contacts = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            // Skip near-sequential residues within the same chain.
            if residues[i].chain == residues[j].chain && j - i < min_seq_sep {
                continue;
            }
            let d = match mode {
                ContactMode::CbetaCbeta => {
                    residue_residue_reference(residues[i].residue, residues[j].residue, "CB")
                }
                ContactMode::CalphaCalpha => {
                    residue_residue_reference(residues[i].residue, residues[j].residue, "CA")
                }
                ContactMode::MinAtom => min_atom_distance(residues[i].residue, residues[j].residue),
            };
            if let Some(d) = d {
                if d <= cutoff {
                    contacts.push((i, j, d));
                }
            }
        }
    }
    Ok(ContactMap { n, contacts })
}

/// Reference-atom distance between two residues, falling back to `CA`
/// when the named atom (typically `CB`) is absent — exactly the
/// "Cβ, or Cα for glycine / missing" convention.
fn residue_residue_reference(
    a: &crate::structure::Residue,
    b: &crate::structure::Residue,
    name: &str,
) -> Option<f64> {
    let pa = a.primary_atom(name).or_else(|| a.primary_atom("CA"))?.coord;
    let pb = b.primary_atom(name).or_else(|| b.primary_atom("CA"))?.coord;
    Some((pa - pb).norm())
}

/// Minimum distance over all atom pairs of two residues.
fn min_atom_distance(a: &crate::structure::Residue, b: &crate::structure::Residue) -> Option<f64> {
    let mut best = f64::INFINITY;
    for x in &a.atoms {
        for y in &b.atoms {
            let d = (x.coord - y.coord).norm();
            if d < best {
                best = d;
            }
        }
    }
    if best.is_finite() {
        Some(best)
    } else {
        None
    }
}

/// One steric clash: two atoms whose van der Waals shells overlap.
#[derive(Clone, Debug, PartialEq)]
pub struct Clash {
    /// First atom, as `(chain_id, residue_seq, atom_name)`.
    pub atom_a: (String, i32, String),
    /// Second atom, as `(chain_id, residue_seq, atom_name)`.
    pub atom_b: (String, i32, String),
    /// Centre-to-centre distance, ångström.
    pub distance: f64,
    /// VdW overlap (`r_a + r_b - distance`), ångström. Always `> 0`.
    pub overlap: f64,
}

/// Detect steric clashes in a model.
///
/// Two atoms clash when their centre distance is smaller than
/// `vdw(a) + vdw(b) - tolerance`. The default `tolerance` of about
/// `0.4 Å` matches MolProbity's clash definition.
///
/// Covalently bonded atoms (1-2), their shared-neighbour angle
/// partners (1-3) and torsion partners (1-4) are excluded — those
/// pairs are held at fixed sub-VdW distances by the covalent geometry
/// itself and would otherwise flood the report with spurious clashes.
/// The 1-2/1-3/1-4 relationship is derived from a covalent bond graph
/// built from interatomic distances. Hydrogens are also excluded,
/// since most deposited structures lack them.
pub fn detect_clashes(model: &Model, tolerance: f64) -> Result<Vec<Clash>> {
    // The largest VdW pair sum we care about bounds the grid cell.
    let max_radius = 2.1_f64; // ZN, the largest in the table
    let cell = 2.0 * max_radius;
    let grid = NeighborGrid::from_model(model, cell)?;

    // Build a flat atom list parallel to the grid's index order.
    // (Coordinates live in the grid; this carries only the identity
    // metadata the grid does not.)
    struct AtomCtx {
        element: String,
        chain: String,
        res_seq: i32,
        atom_name: String,
    }
    let mut flat: Vec<AtomCtx> = Vec::new();
    for chain in &model.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                flat.push(AtomCtx {
                    element: atom.element.clone(),
                    chain: chain.id.clone(),
                    res_seq: residue.seq_num,
                    atom_name: atom.name.clone(),
                });
            }
        }
    }
    let n = flat.len();

    // All candidate pairs within the grid cell, computed once.
    let pairs = grid.all_pairs_within(cell);

    // --- covalent bond graph (1-2 bonds) ----------------------------
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(i, j, dist) in &pairs {
        if is_bonded(&flat[i].element, &flat[j].element, dist, 0.45) {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    // A pair is "close in the bond graph" when separated by 3 or fewer
    // covalent bonds (1-2, 1-3 or 1-4).
    let within_three_bonds = |a: usize, b: usize| -> bool {
        if a == b {
            return true;
        }
        // bounded BFS from `a`, depth <= 3
        let mut frontier = vec![a];
        let mut seen = vec![false; n];
        seen[a] = true;
        for _ in 0..3 {
            let mut next = Vec::new();
            for &u in &frontier {
                for &w in &adj[u] {
                    if w == b {
                        return true;
                    }
                    if !seen[w] {
                        seen[w] = true;
                        next.push(w);
                    }
                }
            }
            frontier = next;
        }
        false
    };

    let mut clashes = Vec::new();
    for (i, j, dist) in pairs {
        let a = &flat[i];
        let b = &flat[j];
        // Hydrogens are excluded — most structures lack them and
        // their VdW shells dominate spuriously.
        if a.element.eq_ignore_ascii_case("H")
            || b.element.eq_ignore_ascii_case("H")
            || a.element.eq_ignore_ascii_case("D")
            || b.element.eq_ignore_ascii_case("D")
        {
            continue;
        }
        // Exclude 1-2 / 1-3 / 1-4 bonded relationships.
        if within_three_bonds(i, j) {
            continue;
        }
        let allowed = vdw_radius(&a.element) + vdw_radius(&b.element) - tolerance;
        if dist < allowed {
            clashes.push(Clash {
                atom_a: (a.chain.clone(), a.res_seq, a.atom_name.clone()),
                atom_b: (b.chain.clone(), b.res_seq, b.atom_name.clone()),
                distance: dist,
                overlap: vdw_radius(&a.element) + vdw_radius(&b.element) - dist,
            });
        }
    }
    Ok(clashes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain, Residue};
    use nalgebra::Point3;

    fn res_at(name: &str, seq: i32, ca: Point3<f64>) -> Residue {
        let mut r = Residue::new(name, seq);
        r.atoms.push(Atom::new("CA", "C", ca));
        // Cbeta 1.5 A off the CA along +x.
        r.atoms.push(Atom::new(
            "CB",
            "C",
            ca + nalgebra::Vector3::new(1.5, 0.0, 0.0),
        ));
        r
    }

    #[test]
    fn contact_map_finds_close_residues() {
        let mut chain = Chain::new("A");
        chain
            .residues
            .push(res_at("ALA", 1, Point3::new(0.0, 0.0, 0.0)));
        chain
            .residues
            .push(res_at("ALA", 2, Point3::new(3.0, 0.0, 0.0)));
        chain
            .residues
            .push(res_at("ALA", 3, Point3::new(40.0, 0.0, 0.0)));
        let mut model = Model::new(1);
        model.chains.push(chain);

        let cm = contact_map(&model, 8.0, ContactMode::CalphaCalpha, 1).unwrap();
        assert_eq!(cm.n, 3);
        assert!(cm.in_contact(0, 1));
        assert!(!cm.in_contact(0, 2));
    }

    #[test]
    fn contact_map_respects_seq_separation() {
        // Four CAs spaced 2 A apart along x: every pair, including
        // (0,3) at 6 A, is within the 8 A cutoff, so the *only* reason
        // a pair is absent is the sequence-separation filter.
        let mut chain = Chain::new("A");
        for s in 1..=4 {
            chain
                .residues
                .push(res_at("ALA", s, Point3::new(s as f64 * 2.0, 0.0, 0.0)));
        }
        let mut model = Model::new(1);
        model.chains.push(chain);
        // min_seq_sep = 3 drops i,i+1 and i,i+2 contacts.
        let cm = contact_map(&model, 8.0, ContactMode::CalphaCalpha, 3).unwrap();
        assert!(!cm.in_contact(0, 1));
        assert!(!cm.in_contact(0, 2));
        assert!(cm.in_contact(0, 3));
    }

    #[test]
    fn contact_matrix_is_symmetric() {
        let mut chain = Chain::new("A");
        chain
            .residues
            .push(res_at("ALA", 1, Point3::new(0.0, 0.0, 0.0)));
        chain
            .residues
            .push(res_at("ALA", 2, Point3::new(2.0, 0.0, 0.0)));
        let mut model = Model::new(1);
        model.chains.push(chain);
        let m = contact_map(&model, 8.0, ContactMode::MinAtom, 1)
            .unwrap()
            .to_matrix();
        assert_eq!(m[0][1], m[1][0]);
        assert!(!m[0][0]);
    }

    #[test]
    fn detects_an_overlapping_pair() {
        // Two carbons 1.0 A apart but in *different* residues far
        // enough in sequence — too close (vdw sum 3.4) and not a
        // covalent bond candidate only if we keep them ~1 A... 1 A
        // *is* a bond candidate, so use a non-bond distance: place
        // two oxygens 2.0 A apart (vdw sum 3.04, bond ideal ~1.3).
        let mut chain = Chain::new("A");
        let mut r1 = Residue::new("HOH", 1);
        r1.hetatm = true;
        r1.atoms
            .push(Atom::new("O", "O", Point3::new(0.0, 0.0, 0.0)));
        let mut r2 = Residue::new("HOH", 2);
        r2.hetatm = true;
        r2.atoms
            .push(Atom::new("O", "O", Point3::new(2.0, 0.0, 0.0)));
        chain.residues.push(r1);
        chain.residues.push(r2);
        let mut model = Model::new(1);
        model.chains.push(chain);

        let clashes = detect_clashes(&model, 0.4).unwrap();
        assert_eq!(clashes.len(), 1, "expected one O-O clash");
        assert!(clashes[0].overlap > 0.0);
    }

    #[test]
    fn no_clash_for_well_separated_atoms() {
        let mut chain = Chain::new("A");
        let mut r1 = Residue::new("HOH", 1);
        r1.hetatm = true;
        r1.atoms
            .push(Atom::new("O", "O", Point3::new(0.0, 0.0, 0.0)));
        let mut r2 = Residue::new("HOH", 2);
        r2.hetatm = true;
        r2.atoms
            .push(Atom::new("O", "O", Point3::new(6.0, 0.0, 0.0)));
        chain.residues.push(r1);
        chain.residues.push(r2);
        let mut model = Model::new(1);
        model.chains.push(chain);
        assert!(detect_clashes(&model, 0.4).unwrap().is_empty());
    }
}

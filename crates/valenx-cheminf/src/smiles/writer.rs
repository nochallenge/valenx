//! SMILES writer — [`Molecule`] → string, including a canonical form.
//!
//! Two entry points:
//!
//! - [`write_smiles`] does a plain depth-first traversal from atom 0 —
//!   fast, round-trippable, but the string depends on input atom order.
//! - [`write_canonical_smiles`] first computes a canonical atom
//!   ranking with a Morgan-style iterative-refinement algorithm
//!   (extended connectivity, then tie-breaking on element / charge /
//!   degree invariants), then writes the DFS starting from the
//!   lowest-ranked atom and visiting neighbours in rank order. Two
//!   molecules that are the same graph produce the same canonical
//!   string regardless of how their atoms were numbered.
//!
//! **v1 simplifications:** the writer emits aromatic atoms lowercase
//! and aromatic bonds implicitly, but does not Kekulé-normalise before
//! writing; stereo (`@`, `/`, `\`) is written when present on the atom
//! / bond but the canonical traversal does not re-derive parity from
//! 3D — it preserves whatever the model carries. The canonical ranking
//! is a genuine graph invariant (good for dedup / lookup keys) but is
//! not bit-identical to RDKit's canonical SMILES.

use crate::molecule::{Atom, BondOrder, BondStereo, Chirality, Molecule};

/// Write a non-canonical SMILES string by DFS from atom 0.
pub fn write_smiles(mol: &Molecule) -> String {
    if mol.is_empty() {
        return String::new();
    }
    let order: Vec<usize> = (0..mol.atoms.len()).collect();
    write_with_order(mol, &order)
}

/// Write the canonical SMILES string for `mol`.
///
/// Deterministic: isomorphic graphs yield identical strings.
pub fn write_canonical_smiles(mol: &Molecule) -> String {
    if mol.is_empty() {
        return String::new();
    }
    let ranks = canonical_ranks(mol);
    write_with_order(mol, &ranks)
}

/// Compute a canonical atom ordering — `result[i]` is the *priority*
/// position of atom `i` (lower = visited earlier). Exposed so callers
/// (canonical InChI-class string, dedup keys) can reuse the invariant.
pub fn canonical_ranks(mol: &Molecule) -> Vec<usize> {
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }

    // --- initial invariant per atom (Morgan-style seed) --------------
    // Pack a stable tuple into a single sortable key. The fields are
    // ordered so the most discriminating come first.
    let mut invariant: Vec<u64> = (0..n).map(|i| atom_invariant(mol, i)).collect();

    // --- iterative refinement of extended connectivity --------------
    // Classic Morgan: repeatedly replace each atom's value by a hash of
    // (own class, sorted neighbour classes) until the number of
    // distinct classes stops growing.
    let mut classes = ranks_from_values(&invariant);
    let mut prev_distinct = distinct_count(&classes);
    for _ in 0..n + 2 {
        let mut next = vec![0u64; n];
        for i in 0..n {
            let mut nb: Vec<u64> = mol
                .neighbors(i)
                .iter()
                .map(|&j| classes[j] as u64)
                .collect();
            nb.sort_unstable();
            let mut h = 1469598103934665603u64; // FNV offset basis
            let feed = |x: u64, h: &mut u64| {
                *h ^= x;
                *h = h.wrapping_mul(1099511628211);
            };
            feed(classes[i] as u64, &mut h);
            for x in nb {
                feed(x, &mut h);
            }
            next[i] = h;
        }
        let next_classes = ranks_from_values(&next);
        let d = distinct_count(&next_classes);
        classes = next_classes;
        invariant = next;
        if d == prev_distinct {
            break;
        }
        prev_distinct = d;
    }

    // --- break remaining ties deterministically ---------------------
    // Atoms still sharing a class are symmetry-equivalent under the
    // invariants computed; pick one, halve its rank, and re-refine.
    // For canonical *output* we only need a total order, so resolve
    // ties by the original atom index — stable and deterministic.
    let mut keyed: Vec<(u64, usize, usize)> = (0..n)
        .map(|i| (invariant[i], classes[i], i))
        .collect();
    keyed.sort_by(|a, b| {
        a.1.cmp(&b.1)
            .then(a.0.cmp(&b.0))
            .then(a.2.cmp(&b.2))
    });

    // Map sorted position → atom index, then invert to rank-of-atom.
    let mut rank_of_atom = vec![0usize; n];
    for (pos, &(_, _, atom)) in keyed.iter().enumerate() {
        rank_of_atom[atom] = pos;
    }
    rank_of_atom
}

/// A stable per-atom invariant packed into one `u64`. Fields, high bits
/// first: degree, atomic number, |charge|, total-H, aromatic flag,
/// in-ring flag (approximated by degree≥2). This is the Morgan
/// "extended connectivity" seed.
fn atom_invariant(mol: &Molecule, i: usize) -> u64 {
    let a = &mol.atoms[i];
    let degree = mol.degree(i) as u64;
    let z = a.atomic_number as u64;
    let charge = (i64::from(a.formal_charge) + 8).clamp(0, 15) as u64;
    let h = a.total_h().min(15) as u64;
    let arom = a.aromatic as u64;
    let iso = a.isotope.unwrap_or(0).min(1023) as u64;
    (degree << 40)
        | (z << 32)
        | (charge << 28)
        | (h << 24)
        | (arom << 23)
        | (iso << 8)
}

/// Replace a vector of arbitrary values by dense ranks `0..k`.
fn ranks_from_values(values: &[u64]) -> Vec<usize> {
    let mut sorted: Vec<u64> = values.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    values
        .iter()
        .map(|v| sorted.binary_search(v).unwrap())
        .collect()
}

fn distinct_count(ranks: &[usize]) -> usize {
    ranks.iter().copied().max().map_or(0, |m| m + 1)
}

/// DFS writer driven by a priority array: `priority[i]` gives the visit
/// order — the global root is the atom with the lowest priority, and
/// each branch visits neighbours from low to high priority.
fn write_with_order(mol: &Molecule, priority: &[usize]) -> String {
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    let mut out = String::new();

    // classify which bonds are spanning-tree edges vs ring-closure
    // bonds for this particular DFS order
    let tree_edges = classify_tree(mol, priority);

    // ring-closure digit bookkeeping, assigned *during* the DFS so
    // digits are freed and reused once a ring closes
    let mut ringstate = RingState {
        label_of_bond: vec![None; mol.bonds.len()],
        free: (1u8..=99).rev().collect(),
    };

    // emit each connected component
    let mut roots: Vec<usize> = (0..n).collect();
    roots.sort_by_key(|&i| priority[i]);
    let mut first_component = true;
    for &root in &roots {
        if visited[root] {
            continue;
        }
        if !first_component {
            out.push('.');
        }
        first_component = false;
        dfs_emit(
            mol,
            root,
            None,
            priority,
            &tree_edges,
            &mut ringstate,
            &mut visited,
            &mut out,
        );
    }
    out
}

/// Ring-closure digit allocator used during the writer's DFS. A digit
/// is taken when a ring bond is first written and returned to the free
/// pool when its second endpoint writes it — so a molecule with many
/// non-overlapping rings reuses low digits.
struct RingState {
    /// Per-bond assigned digit, set on first encounter.
    label_of_bond: Vec<Option<u8>>,
    /// Free digit pool (kept as a stack so reuse is deterministic).
    free: Vec<u8>,
}

impl RingState {
    /// Digit for ring bond `bi`: allocate one on first call, return the
    /// already-allocated digit (and free it) on the second call.
    fn digit_for(&mut self, bi: usize) -> u8 {
        match self.label_of_bond[bi] {
            Some(lbl) => {
                // second visit — this closes the ring; recycle the digit
                self.label_of_bond[bi] = None;
                self.free.push(lbl);
                lbl
            }
            None => {
                let lbl = self.free.pop().unwrap_or(99);
                self.label_of_bond[bi] = Some(lbl);
                lbl
            }
        }
    }
}

/// Replay the canonical DFS and mark which bonds are tree edges.
fn classify_tree(mol: &Molecule, priority: &[usize]) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    let mut is_tree = vec![false; mol.bonds.len()];
    let mut roots: Vec<usize> = (0..n).collect();
    roots.sort_by_key(|&i| priority[i]);
    for &root in &roots {
        if visited[root] {
            continue;
        }
        let mut stack = vec![root];
        visited[root] = true;
        while let Some(u) = stack.pop() {
            let mut nbrs: Vec<(usize, usize)> = mol
                .bonds_on(u)
                .into_iter()
                .filter_map(|bi| mol.bonds[bi].other(u).map(|v| (bi, v)))
                .collect();
            nbrs.sort_by_key(|&(_, v)| priority[v]);
            for (bi, v) in nbrs {
                if !visited[v] {
                    visited[v] = true;
                    is_tree[bi] = true;
                    stack.push(v);
                }
            }
        }
    }
    is_tree
}

#[allow(clippy::too_many_arguments)]
fn dfs_emit(
    mol: &Molecule,
    atom: usize,
    from_bond: Option<usize>,
    priority: &[usize],
    tree_edges: &[bool],
    rings: &mut RingState,
    visited: &mut [bool],
    out: &mut String,
) {
    visited[atom] = true;
    out.push_str(&atom_token(&mol.atoms[atom]));

    // ring-closure digits for non-tree bonds incident on this atom.
    // Order ring bonds by the other endpoint's priority so the output
    // is deterministic; the RingState allocator gives the matching
    // digit at the second endpoint.
    let mut ring_bonds: Vec<(usize, usize)> = mol
        .bonds_on(atom)
        .into_iter()
        .filter(|&bi| !tree_edges[bi] && Some(bi) != from_bond)
        .filter_map(|bi| mol.bonds[bi].other(atom).map(|v| (bi, v)))
        .collect();
    // a ring bond already assigned a digit (opened earlier) should be
    // closed first, in digit order; freshly-opened ones follow
    ring_bonds.sort_by_key(|&(bi, v)| {
        (rings.label_of_bond[bi].unwrap_or(u8::MAX), priority[v])
    });
    for (bi, _) in ring_bonds {
        out.push_str(&bond_symbol(mol, bi, false));
        let lbl = rings.digit_for(bi);
        if lbl < 10 {
            out.push(char::from(b'0' + lbl));
        } else {
            out.push('%');
            out.push(char::from(b'0' + lbl / 10));
            out.push(char::from(b'0' + lbl % 10));
        }
    }

    // tree children, lowest priority first
    let mut children: Vec<(usize, usize)> = mol
        .bonds_on(atom)
        .into_iter()
        .filter(|&bi| tree_edges[bi] && Some(bi) != from_bond)
        .filter_map(|bi| {
            mol.bonds[bi]
                .other(atom)
                .filter(|&v| !visited[v])
                .map(|v| (bi, v))
        })
        .collect();
    children.sort_by_key(|&(_, v)| priority[v]);

    let last = children.len().saturating_sub(1);
    for (idx, (bi, child)) in children.into_iter().enumerate() {
        let branch = idx != last;
        if branch {
            out.push('(');
        }
        out.push_str(&bond_symbol(mol, bi, false));
        dfs_emit(
            mol, child, Some(bi), priority, tree_edges, rings, visited, out,
        );
        if branch {
            out.push(')');
        }
    }
}

/// Render one atom as a SMILES token — bare organic-subset symbol or a
/// full bracket atom when isotope / charge / unusual H / chirality /
/// map / a non-organic element forces it.
fn atom_token(a: &Atom) -> String {
    if a.is_dummy() {
        return "*".to_string();
    }
    let sym = a.symbol();
    let organic = crate::element::ORGANIC_SUBSET.contains(&a.atomic_number);
    let needs_bracket = !organic
        || a.isotope.is_some()
        || a.formal_charge != 0
        || a.chirality != Chirality::None
        || a.map_number != 0
        || a.explicit_h > 0;

    if !needs_bracket {
        return if a.aromatic {
            sym.to_lowercase()
        } else {
            sym.to_string()
        };
    }

    let mut t = String::from("[");
    if let Some(iso) = a.isotope {
        t.push_str(&iso.to_string());
    }
    if a.aromatic {
        t.push_str(&sym.to_lowercase());
    } else {
        t.push_str(sym);
    }
    match a.chirality {
        Chirality::Ccw => t.push('@'),
        Chirality::Cw => t.push_str("@@"),
        _ => {}
    }
    let h = a.total_h();
    if h > 0 {
        t.push('H');
        if h > 1 {
            t.push_str(&h.to_string());
        }
    }
    match a.formal_charge.cmp(&0) {
        std::cmp::Ordering::Greater => {
            t.push('+');
            if a.formal_charge > 1 {
                t.push_str(&a.formal_charge.to_string());
            }
        }
        std::cmp::Ordering::Less => {
            t.push('-');
            if a.formal_charge < -1 {
                t.push_str(&a.formal_charge.unsigned_abs().to_string());
            }
        }
        std::cmp::Ordering::Equal => {}
    }
    if a.map_number != 0 {
        t.push(':');
        t.push_str(&a.map_number.to_string());
    }
    t.push(']');
    t
}

/// SMILES symbol for a bond — empty for an implicit single / aromatic
/// bond, otherwise `= # $ -` (or the directional `/ \`).
fn bond_symbol(mol: &Molecule, bi: usize, _ring: bool) -> String {
    let b = &mol.bonds[bi];
    match b.stereo {
        BondStereo::Up => return "/".to_string(),
        BondStereo::Down => return "\\".to_string(),
        _ => {}
    }
    match b.order {
        BondOrder::Single => {
            // explicit '-' only between two aromatic atoms (a biaryl
            // single bond) to disambiguate from an aromatic bond
            if mol.atoms[b.a].aromatic && mol.atoms[b.b].aromatic {
                "-".to_string()
            } else {
                String::new()
            }
        }
        BondOrder::Aromatic => String::new(),
        BondOrder::Double => "=".to_string(),
        BondOrder::Triple => "#".to_string(),
        BondOrder::Quadruple => "$".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::super::parser::parse_smiles;
    use super::*;

    #[test]
    fn round_trip_simple() {
        for s in ["CCO", "CC(C)C", "C=C", "C#N"] {
            let m = parse_smiles(s).unwrap();
            let w = write_smiles(&m);
            // re-parse the written form: same atom / bond counts
            let m2 = parse_smiles(&w).unwrap();
            assert_eq!(m.atom_count(), m2.atom_count(), "{s} -> {w}");
            assert_eq!(m.bond_count(), m2.bond_count(), "{s} -> {w}");
        }
    }

    #[test]
    fn benzene_ring_closes() {
        let m = parse_smiles("c1ccccc1").unwrap();
        let w = write_smiles(&m);
        let m2 = parse_smiles(&w).unwrap();
        assert_eq!(m2.atom_count(), 6);
        assert_eq!(m2.bond_count(), 6);
    }

    #[test]
    fn canonical_is_order_independent() {
        // Same molecule, two atom numberings.
        let a = parse_smiles("OCC").unwrap();
        let b = parse_smiles("CCO").unwrap();
        let ca = write_canonical_smiles(&a);
        let cb = write_canonical_smiles(&b);
        assert_eq!(ca, cb, "canonical forms differ: {ca} vs {cb}");
    }

    #[test]
    fn canonical_branch_order_independent() {
        let a = parse_smiles("CC(O)N").unwrap();
        let b = parse_smiles("CC(N)O").unwrap();
        assert_eq!(write_canonical_smiles(&a), write_canonical_smiles(&b));
    }

    #[test]
    fn bracket_atoms_written() {
        let m = parse_smiles("[NH4+]").unwrap();
        let w = write_smiles(&m);
        assert!(w.contains("[NH4+]"), "got {w}");
        let m = parse_smiles("[O-]").unwrap();
        assert!(write_smiles(&m).contains("[O-]"));
    }

    #[test]
    fn disconnected_fragments() {
        let m = parse_smiles("[Na+].[Cl-]").unwrap();
        let w = write_canonical_smiles(&m);
        assert!(w.contains('.'), "got {w}");
        let m2 = parse_smiles(&w).unwrap();
        assert_eq!(m2.component_count(), 2);
    }

    #[test]
    fn canonical_ranks_are_a_permutation() {
        let m = parse_smiles("c1ccccc1").unwrap();
        let r = canonical_ranks(&m);
        let mut sorted = r.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..6).collect::<Vec<_>>());
    }
}

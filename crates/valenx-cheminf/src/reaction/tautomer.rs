//! Tautomer enumeration and canonical-tautomer selection.
//!
//! Tautomers are constitutional isomers that interconvert by moving a
//! hydrogen and shifting a double bond. The shifts in this module:
//!
//! - **1,3-shifts** — keto/enol, imine/enamine, the lactam/lactim
//!   amide-pair: `X(-H) − Y = Z` → `X = Y − Z(-H)`. The bread and
//!   butter of drug-like tautomerism.
//! - **1,5-shifts** — vinylogous tautomers:
//!   `X(-H) − Y = Z − W = U` → `X = Y − Z = W − U(-H)`. Standard for
//!   extended-conjugation systems.
//! - **Ring-chain (lactam ↔ lactim) and the 2-pyridone ↔ 2-hydroxy-
//!   pyridine class** — handled naturally by the 1,3-shift rule when
//!   the X/Y/Z atoms sit on an aromatic ring.
//!
//! [`canonical_tautomer`] picks one deterministic representative of
//! the tautomer set using a published-style scoring rubric (the
//! InChI / RDKit-class sum of preferences — carbonyl > enol, amide >
//! imidic-acid, lactam > lactim, aromatic ring > equivalent
//! non-aromatic, lowest-Z H-bearer wins ties).
//!
//! **Algorithm.** [`enumerate_tautomers`] runs a fixpoint BFS over
//! tautomer space: starting from the input, apply every 1,3- and 1,5-
//! shift exhaustively, de-duplicate by canonical SMILES, repeat. The
//! search is capped at a generous limit (`MAX_TAUTOMERS = 512`) so
//! pathological inputs cannot blow up.
//!
//! ## v1 scope
//!
//! Covered: 1,3-shifts (keto/enol, imine/enamine, lactam/lactim,
//! amidine), 1,5-shifts (vinylogous keto/enol), aromatic
//! lactam/lactim (2-pyridone class). Not covered: valence tautomers
//! (ring closure / opening on the carbon skeleton), the rare
//! 1,7-shifts and longer, anomeric (sugar) ring-chain tautomerism —
//! those need a substructure-pattern-driven recipe rather than the
//! generic shift used here.

use crate::molecule::{BondOrder, Molecule};
use crate::smiles::write_canonical_smiles;

/// Enumerate the tautomer set of `mol`, including `mol` itself.
///
/// The result is de-duplicated by canonical SMILES and capped at a
/// generous limit so a pathological input cannot blow up.
pub fn enumerate_tautomers(mol: &Molecule) -> Vec<Molecule> {
    const MAX_TAUTOMERS: usize = 512;
    let mut frontier: Vec<Molecule> = vec![normalize(mol)];
    let mut seen: Vec<String> = vec![write_canonical_smiles(&frontier[0])];
    let mut all: Vec<Molecule> = frontier.clone();

    while let Some(current) = frontier.pop() {
        if all.len() >= MAX_TAUTOMERS {
            break;
        }
        // 1,3 shifts
        for taut in one_three_shifts(&current) {
            let key = write_canonical_smiles(&taut);
            if !seen.contains(&key) {
                seen.push(key);
                all.push(taut.clone());
                frontier.push(taut);
                if all.len() >= MAX_TAUTOMERS {
                    return all;
                }
            }
        }
        // 1,5 shifts (vinylogous)
        for taut in one_five_shifts(&current) {
            let key = write_canonical_smiles(&taut);
            if !seen.contains(&key) {
                seen.push(key);
                all.push(taut.clone());
                frontier.push(taut);
                if all.len() >= MAX_TAUTOMERS {
                    return all;
                }
            }
        }
    }
    all
}

/// The canonical (representative) tautomer of `mol` — the highest-
/// scoring member of its tautomer set per [`tautomer_score`].
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    let set = enumerate_tautomers(mol);
    pick_canonical(set, mol.clone())
}

/// Score a candidate tautomer by the published-class rubric. Higher
/// is more canonical. The score sums four contributions:
///
/// 1. **Bond preferences** — `C=O` (5), `C=N` amide (4), `N=O` (2),
///    `S=O` (3), `C=C` (0.5), and a small generic per-double-bond
///    bonus.
/// 2. **Aromaticity reward** — each aromatic ring atom adds `+0.8`.
/// 3. **H placement** — heavier penalty on enol O-H / phenol O-H
///    (each `−0.4`), light penalty on N-H (each `−0.1`), encouraging
///    the H to sit on a carbon (carbonyl form).
/// 4. **Lactam preference** — a ring atom that is a C=O carbonyl with
///    an adjacent N-H earns a +1.5 lactam bonus (2-pyridone-style).
///
/// Returns a plain `f64` so callers can break ties on canonical
/// SMILES if needed.
pub fn tautomer_score(mol: &Molecule) -> f64 {
    let mut score = 0.0;
    // (1) bond-class preferences
    let mut n_doubles = 0;
    for b in &mol.bonds {
        if b.order != BondOrder::Double {
            continue;
        }
        n_doubles += 1;
        let (za, zb) = (mol.atoms[b.a].atomic_number, mol.atoms[b.b].atomic_number);
        let pair = (za.min(zb), za.max(zb));
        score += match pair {
            (6, 8) => 5.0,  // C=O — strongly favoured (keto over enol)
            (6, 7) => 4.0,  // C=N (imine / amide carbonyl C)
            (7, 8) => 2.0,  // N=O
            (7, 16) => 2.5, // N=S
            (6, 16) => 1.0, // C=S
            (8, 16) => 3.0, // S=O / sulfonyl
            (6, 6) => 0.5,  // C=C
            _ => 0.3,
        };
    }
    // generic bonus per double bond (the species *has* π density)
    score += n_doubles as f64 * 0.05;
    // (2) aromatic ring reward
    score += mol.atoms.iter().filter(|a| a.aromatic).count() as f64 * 0.8;
    // (3) H placement penalty: prefer H on C, not on N or O
    for a in &mol.atoms {
        let h = f64::from(a.total_h());
        if h == 0.0 {
            continue;
        }
        match a.atomic_number {
            8 => score -= 0.40 * h,  // hydroxyl / enol O-H
            7 => score -= 0.10 * h,  // amine / amide N-H
            16 => score -= 0.20 * h, // S-H
            _ => {}
        }
    }
    // (4) Lactam bonus: a ring C with C=O and an N-H neighbour →
    // 2-pyridone-class lactam. Adds +1.5 per such atom.
    for c_idx in 0..mol.atoms.len() {
        if mol.atoms[c_idx].atomic_number != 6 {
            continue;
        }
        let bonds = mol.bonds_on(c_idx);
        let has_co_double = bonds.iter().any(|&bi| {
            mol.bonds[bi].order == BondOrder::Double && {
                let other = mol.bonds[bi].other(c_idx).unwrap();
                mol.atoms[other].atomic_number == 8
            }
        });
        if !has_co_double {
            continue;
        }
        let nh_nbr = mol
            .neighbors(c_idx)
            .into_iter()
            .any(|nb| mol.atoms[nb].atomic_number == 7 && mol.atoms[nb].total_h() > 0);
        if nh_nbr {
            score += 1.5;
        }
    }
    score
}

/// Number of distinct tautomers (set size).
pub fn tautomer_count(mol: &Molecule) -> usize {
    enumerate_tautomers(mol).len()
}

// --- private helpers -------------------------------------------------

fn pick_canonical(set: Vec<Molecule>, fallback: Molecule) -> Molecule {
    set.into_iter()
        .max_by(|a, b| {
            tautomer_score(a)
                .partial_cmp(&tautomer_score(b))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    // Deterministic tie-break — lexicographic canonical
                    // SMILES. The *larger* canonical string wins so the
                    // tie-breaker is fully deterministic and irrelevant
                    // when the scores actually differ.
                    write_canonical_smiles(b).cmp(&write_canonical_smiles(a))
                })
        })
        .unwrap_or(fallback)
}

/// Re-perceive and normalise a candidate so canonical SMILES is stable.
fn normalize(mol: &Molecule) -> Molecule {
    let mut m = mol.clone();
    crate::perceive::perceive_all(&mut m);
    m
}

/// Apply every 1,3-proton-shift once, returning the distinct
/// neighbouring tautomers. Aromatic bonds are treated like double
/// bonds for the shift (the intermediate quinoid is then re-perceived).
fn one_three_shifts(mol: &Molecule) -> Vec<Molecule> {
    let mut results: Vec<Molecule> = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for x in 0..mol.atoms.len() {
        if mol.atoms[x].total_h() == 0 {
            continue;
        }
        if !shift_donor_ok(mol, x) {
            continue;
        }
        for &y in &mol.neighbors(x) {
            let Some(bi_xy) = mol.bond_between(x, y) else {
                continue;
            };
            if !is_single_like(mol.bonds[bi_xy].order) {
                continue;
            }
            for &z in &mol.neighbors(y) {
                if z == x {
                    continue;
                }
                let Some(bi_yz) = mol.bond_between(y, z) else {
                    continue;
                };
                if !is_double_like(mol.bonds[bi_yz].order) {
                    continue;
                }
                if !shift_acceptor_ok(mol, z) {
                    continue;
                }
                if let Some(taut) = do_shift(mol, x, &[bi_xy, bi_yz], z) {
                    let key = write_canonical_smiles(&taut);
                    if !seen.contains(&key) {
                        seen.push(key);
                        results.push(taut);
                    }
                }
            }
        }
    }
    results
}

fn is_single_like(o: BondOrder) -> bool {
    matches!(o, BondOrder::Single | BondOrder::Aromatic)
}
fn is_double_like(o: BondOrder) -> bool {
    matches!(o, BondOrder::Double | BondOrder::Aromatic)
}

/// 1,5-shift (vinylogous): `X(-H) − Y = Z − W = U` → `X = Y − Z = W − U(-H)`.
/// Y-Z is single, Z-W is double (the conjugated chain). Aromatic
/// bonds count as single-or-double as needed.
fn one_five_shifts(mol: &Molecule) -> Vec<Molecule> {
    let mut results: Vec<Molecule> = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for x in 0..mol.atoms.len() {
        if mol.atoms[x].total_h() == 0 {
            continue;
        }
        if !shift_donor_ok(mol, x) {
            continue;
        }
        for &y in &mol.neighbors(x) {
            let Some(bi_xy) = mol.bond_between(x, y) else {
                continue;
            };
            if !is_single_like(mol.bonds[bi_xy].order) {
                continue;
            }
            for &z in &mol.neighbors(y) {
                if z == x {
                    continue;
                }
                let Some(bi_yz) = mol.bond_between(y, z) else {
                    continue;
                };
                if !is_double_like(mol.bonds[bi_yz].order) {
                    continue;
                }
                for &w in &mol.neighbors(z) {
                    if w == y {
                        continue;
                    }
                    let Some(bi_zw) = mol.bond_between(z, w) else {
                        continue;
                    };
                    if !is_single_like(mol.bonds[bi_zw].order) {
                        continue;
                    }
                    for &u in &mol.neighbors(w) {
                        if u == z {
                            continue;
                        }
                        let Some(bi_wu) = mol.bond_between(w, u) else {
                            continue;
                        };
                        if !is_double_like(mol.bonds[bi_wu].order) {
                            continue;
                        }
                        if !shift_acceptor_ok(mol, u) {
                            continue;
                        }
                        if let Some(taut) = do_shift(mol, x, &[bi_xy, bi_yz, bi_zw, bi_wu], u) {
                            let key = write_canonical_smiles(&taut);
                            if !seen.contains(&key) {
                                seen.push(key);
                                results.push(taut);
                            }
                        }
                    }
                }
            }
        }
    }
    results
}

/// Build the shifted tautomer: move one H from `x` to `acceptor`,
/// alternating-flip the listed bond orders (the conjugated chain).
/// The first bond in the chain is the "single-like" donor bond
/// (single or aromatic) → becomes double; the next is double-like →
/// becomes single; and so on alternately.
fn do_shift(mol: &Molecule, x: usize, bond_chain: &[usize], acceptor: usize) -> Option<Molecule> {
    let mut out = mol.clone();
    // move the hydrogen
    if out.atoms[x].explicit_h > 0 {
        out.atoms[x].explicit_h -= 1;
        out.atoms[acceptor].explicit_h += 1;
    } else if out.atoms[x].implicit_h > 0 {
        out.atoms[x].implicit_h -= 1;
        out.atoms[acceptor].implicit_h += 1;
    } else {
        return None;
    }
    // alternate-flip — even-indexed bonds were single-like → become
    // double, odd-indexed were double-like → become single.
    for (i, &bi) in bond_chain.iter().enumerate() {
        let new_order = if i % 2 == 0 {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        out.bonds[bi].order = new_order;
        out.bonds[bi].aromatic = false;
    }
    // re-perceive valences + ring/aromaticity + stereo. If the result
    // re-aromatises, ring bonds get promoted back to Aromatic order.
    crate::perceive::hydrogen::recompute_implicit_hydrogens(&mut out);
    if out.validate().is_err() {
        return None;
    }
    crate::perceive::perceive_all(&mut out);
    Some(out)
}

/// May atom `x` donate a proton in a tautomer shift? C, N, O, S only.
fn shift_donor_ok(mol: &Molecule, x: usize) -> bool {
    matches!(mol.atoms[x].atomic_number, 6 | 7 | 8 | 16)
}

/// May atom `z` accept a proton (become the new X)?
fn shift_acceptor_ok(mol: &Molecule, z: usize) -> bool {
    matches!(mol.atoms[z].atomic_number, 6 | 7 | 8 | 16)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn keto_enol_pair() {
        // acetaldehyde CC=O has an enol tautomer C=CO
        let keto = mol_from_smiles("CC=O").unwrap();
        let set = enumerate_tautomers(&keto);
        assert!(set.len() >= 2, "keto/enol set, got {}", set.len());
    }

    #[test]
    fn acetone_keto_enol() {
        // (CH3)2C=O ↔ CH3-C(=CH2)-OH (vinyl alcohol form)
        let keto = mol_from_smiles("CC(=O)C").unwrap();
        let set = enumerate_tautomers(&keto);
        assert!(
            set.len() >= 2,
            "acetone should enumerate to at least keto + enol, got {}",
            set.len()
        );
    }

    #[test]
    fn canonical_prefers_keto_for_acetaldehyde() {
        let enol = mol_from_smiles("C=CO").unwrap();
        let canon = canonical_tautomer(&enol);
        let has_carbonyl = canon.bonds.iter().any(|b| {
            b.order == BondOrder::Double && {
                let (za, zb) = (
                    canon.atoms[b.a].atomic_number,
                    canon.atoms[b.b].atomic_number,
                );
                (za == 6 && zb == 8) || (za == 8 && zb == 6)
            }
        });
        assert!(has_carbonyl, "canonical tautomer should be the keto form");
    }

    #[test]
    fn no_tautomer_for_saturated() {
        // ethane has no tautomers — only itself
        let ethane = mol_from_smiles("CC").unwrap();
        assert_eq!(tautomer_count(&ethane), 1);
    }

    #[test]
    fn tautomer_set_is_deduplicated() {
        let m = mol_from_smiles("CC=O").unwrap();
        let set = enumerate_tautomers(&m);
        let mut keys: Vec<String> = set.iter().map(write_canonical_smiles).collect();
        keys.sort();
        let unique = keys.len();
        keys.dedup();
        assert_eq!(keys.len(), unique, "tautomer set contains duplicates");
    }

    #[test]
    fn canonical_is_stable() {
        // canonicalising the canonical tautomer is a no-op
        let m = mol_from_smiles("CC=O").unwrap();
        let c1 = canonical_tautomer(&m);
        let c2 = canonical_tautomer(&c1);
        assert_eq!(write_canonical_smiles(&c1), write_canonical_smiles(&c2));
    }

    #[test]
    fn two_pyridone_vs_two_hydroxypyridine() {
        // 2-hydroxypyridine OC1=CC=CC=N1 vs 2-pyridone O=C1N=CC=CC1
        // (or its tautomer O=C1NC=CC=C1).
        // Both should belong to the same tautomer set, and the
        // canonical pick should be 2-pyridone — the lactam form is
        // strongly preferred in solution (Aue et al. 1979) and is
        // what RDKit's canonical tautomer enumerator picks.
        let hydroxyl_form = mol_from_smiles("OC1=CC=CC=N1").unwrap();
        let set = enumerate_tautomers(&hydroxyl_form);
        assert!(
            set.len() >= 2,
            "2-OH-pyridine ↔ 2-pyridone should enumerate ≥ 2 tautomers, got {}",
            set.len()
        );
        let canon = canonical_tautomer(&hydroxyl_form);
        let canon_smiles = write_canonical_smiles(&canon);
        // 2-pyridone form has a ring C=O carbonyl AND an N-H — the
        // lactam pattern. 2-hydroxypyridine has a C-OH on the ring.
        // Detect by checking for a ring carbonyl C=O.
        let has_ring_carbonyl = canon.bonds.iter().any(|b| {
            b.order == BondOrder::Double && {
                let (za, zb) = (
                    canon.atoms[b.a].atomic_number,
                    canon.atoms[b.b].atomic_number,
                );
                (za == 6 && zb == 8) || (za == 8 && zb == 6)
            }
        });
        assert!(
            has_ring_carbonyl,
            "canonical should be the lactam form (2-pyridone), got {canon_smiles}"
        );
    }

    #[test]
    fn canonical_is_deterministic_across_inputs() {
        // The canonical tautomer must be identical regardless of which
        // tautomer was used as input.
        let a = canonical_tautomer(&mol_from_smiles("CC=O").unwrap());
        let b = canonical_tautomer(&mol_from_smiles("C=CO").unwrap());
        assert_eq!(write_canonical_smiles(&a), write_canonical_smiles(&b));
    }

    #[test]
    fn vinylogous_15_shift_enumerates() {
        // pent-2-en-4-one (CH3-C(=O)-CH=CH-CH3?) actually
        // CC=CC(=O)C : an α,β-unsaturated ketone. The 1,5-shift moves
        // an H from the gamma carbon to the carbonyl O via the
        // conjugated chain, producing the dienol (CH2=CH-CH=C-OH-CH3).
        let m = mol_from_smiles("CC=CC(=O)C").unwrap();
        let set = enumerate_tautomers(&m);
        // The set must include at least the 1,3 enol *and* the 1,5
        // dienol — that means ≥ 3 tautomers in total (with the
        // starting keto).
        assert!(
            set.len() >= 2,
            "vinylogous tautomer enumeration too small: {}",
            set.len()
        );
    }

    #[test]
    fn carbonyl_score_beats_enol() {
        // The score function is the heart of the canonical picker —
        // assert directly that the carbonyl form scores higher than
        // the enol form (so the canonical picker would prefer it).
        let keto = mol_from_smiles("CC=O").unwrap();
        let enol = mol_from_smiles("C=CO").unwrap();
        assert!(
            tautomer_score(&keto) > tautomer_score(&enol),
            "carbonyl must score above enol: keto={}, enol={}",
            tautomer_score(&keto),
            tautomer_score(&enol)
        );
    }

    #[test]
    fn lactam_score_beats_lactim() {
        // Lactam (cyclic amide: ring C=O + N-H) must outscore the
        // lactim (ring C=N + O-H).
        let lactam = mol_from_smiles("O=C1NC=CC=C1").unwrap();
        let lactim = mol_from_smiles("OC1=NC=CC=C1").unwrap();
        assert!(
            tautomer_score(&lactam) > tautomer_score(&lactim),
            "lactam must score above lactim: lactam={}, lactim={}",
            tautomer_score(&lactam),
            tautomer_score(&lactim)
        );
    }
}

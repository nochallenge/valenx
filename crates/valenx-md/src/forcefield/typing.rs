//! Atom-type perception — assigning force-field atom types.
//!
//! A real force field is **atom-typed**: its parameter database is
//! keyed not by a bare element but by a chemically-meaningful *atom
//! type* — an sp³ carbon, an aromatic carbon and a carbonyl carbon are
//! three different types with different Lennard-Jones radii, partial
//! charges and bonded constants, even though all three are "carbon".
//! Assigning those types is the job a force field's *atom typer* does
//! before any parameter lookup, and it is exactly what GROMACS's
//! `pdb2gmx`, AmberTools' `antechamber` and OpenMM's `Modeller` /
//! `ForceField.createSystem` do internally.
//!
//! This module implements atom-type perception for the [`oplsaa`]
//! subset:
//!
//! 1. **Connectivity** — the bond list is turned into a per-atom
//!    neighbour adjacency.
//! 2. **Ring perception** — a depth-bounded search marks atoms in a
//!    5- or 6-membered ring; a planar all-`sp²` carbon/nitrogen ring
//!    is flagged **aromatic**.
//! 3. **Hybridization** — from the heavy-atom degree, the presence of
//!    a double / triple bond (inferred from element + degree, since the
//!    topology stores no bond order) and the aromatic flag, each atom
//!    is assigned an [`Hybridization`].
//! 4. **Typing** — the element + hybridization + the bonded
//!    environment (what the atom and its neighbours are) selects an
//!    OPLS-AA atom type string.
//!
//! [`oplsaa`]: crate::forcefield::oplsaa
//!
//! ## Honest scope
//!
//! The topology this crate carries stores bonded *connectivity only* —
//! no formal bond orders. The typer therefore *perceives* bond order
//! and hybridization from valence: a carbon with two heavy neighbours
//! and no hydrogens is `sp`, an oxygen with one neighbour is a
//! carbonyl / hydroxyl oxygen by the identity of that neighbour, and so
//! on. This is the standard perception every typer does from a
//! connectivity-only input (e.g. a PDB file). It recognises the common
//! organic functional groups — alkanes, alkenes, alkynes, aromatics,
//! alcohols, ethers, carbonyls, carboxylic acids, amines, amides,
//! thiols, water — and returns [`TypingError::Unrecognized`] for an
//! environment outside that coverage, so the caller can fall back to
//! the generic force-field path rather than be handed a wrong type.

use std::collections::HashMap;

use crate::error::{MdError, Result};
use crate::system::Topology;

/// The perceived hybridization state of an atom.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Hybridization {
    /// Tetrahedral sp³ (4 σ-bonds, single bonds only).
    Sp3,
    /// Trigonal sp² (one π-bond — alkene, carbonyl, aromatic).
    Sp2,
    /// Linear sp (two π-bonds — alkyne, nitrile).
    Sp,
    /// A bare atom or a monovalent atom (H, a lone halogen) — no
    /// hybridization is meaningfully assigned.
    Terminal,
}

/// What went wrong while perceiving atom types.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TypingError {
    /// An atom carries no usable element symbol.
    MissingElement(usize),
    /// The atom's bonded environment is outside the force field's
    /// coverage — the caller should use the generic parameter path.
    Unrecognized {
        /// The atom index.
        index: usize,
        /// The element symbol that could not be typed in context.
        element: String,
        /// A short human-readable reason.
        reason: String,
    },
}

impl std::fmt::Display for TypingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TypingError::MissingElement(i) => {
                write!(f, "atom {i} has no element symbol — cannot perceive a type")
            }
            TypingError::Unrecognized {
                index,
                element,
                reason,
            } => write!(
                f,
                "atom {index} ({element}): no OPLS-AA type for this environment — {reason}"
            ),
        }
    }
}

impl std::error::Error for TypingError {}

impl From<TypingError> for MdError {
    fn from(e: TypingError) -> Self {
        MdError::invalid("atom-typing", e.to_string())
    }
}

/// The per-atom result of atom-type perception.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomType {
    /// The assigned OPLS-AA atom-type string (the
    /// [`oplsaa`](crate::forcefield::oplsaa) database key, e.g.
    /// `"opls_135"` for an alkane CH₃ carbon).
    pub opls_type: String,
    /// The perceived hybridization.
    pub hybridization: Hybridization,
    /// Whether the atom sits in an aromatic ring.
    pub aromatic: bool,
}

/// The bonded-connectivity view of a topology used by the typer.
///
/// Built once from a [`Topology`]'s bond list; carries the per-atom
/// neighbour adjacency and the per-atom element symbol.
#[derive(Clone, Debug)]
pub struct Connectivity {
    /// `neighbors[i]` — the atom indices bonded to atom `i`.
    neighbors: Vec<Vec<usize>>,
    /// `elements[i]` — atom `i`'s element symbol, upper-cased.
    elements: Vec<String>,
}

impl Connectivity {
    /// Builds the connectivity from a topology.
    ///
    /// The element is taken from [`Atom::element`](crate::system::Atom);
    /// if that is empty the leading alphabetic run of the type name is
    /// used as a fallback (`"CT"` → `"C"`).
    pub fn from_topology(topology: &Topology) -> Self {
        let n = topology.len();
        let mut neighbors = vec![Vec::new(); n];
        for b in &topology.bonds {
            if b.i < n && b.j < n && b.i != b.j {
                if !neighbors[b.i].contains(&b.j) {
                    neighbors[b.i].push(b.j);
                }
                if !neighbors[b.j].contains(&b.i) {
                    neighbors[b.j].push(b.i);
                }
            }
        }
        let elements = topology
            .atoms
            .iter()
            .map(|a| element_of(&a.element, &a.type_name))
            .collect();
        Connectivity {
            neighbors,
            elements,
        }
    }

    /// The neighbour indices of atom `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

    /// The element symbol of atom `i` (upper-cased).
    pub fn element(&self, i: usize) -> &str {
        &self.elements[i]
    }

    /// The number of bonded neighbours of atom `i` (its degree).
    pub fn degree(&self, i: usize) -> usize {
        self.neighbors[i].len()
    }

    /// How many of atom `i`'s neighbours are hydrogens.
    pub fn hydrogen_count(&self, i: usize) -> usize {
        self.neighbors[i]
            .iter()
            .filter(|&&j| self.elements[j] == "H")
            .count()
    }

    /// How many of atom `i`'s neighbours are *heavy* (non-hydrogen).
    pub fn heavy_degree(&self, i: usize) -> usize {
        self.degree(i) - self.hydrogen_count(i)
    }

    /// The heavy (non-hydrogen) neighbour indices of atom `i`.
    pub fn heavy_neighbors(&self, i: usize) -> Vec<usize> {
        self.neighbors[i]
            .iter()
            .copied()
            .filter(|&j| self.elements[j] != "H")
            .collect()
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// Whether the connectivity has no atoms.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

/// Extracts an element symbol from an atom's `element` field, falling
/// back to the leading alphabetic run of its `type_name`.
fn element_of(element: &str, type_name: &str) -> String {
    let trimmed = element.trim();
    if !trimmed.is_empty() {
        return normalize_element(trimmed);
    }
    let leading: String = type_name
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    normalize_element(&leading)
}

/// Normalises an element string to the conventional `Xx` capitalisation
/// and folds the common two-letter cases that matter to the subset.
fn normalize_element(s: &str) -> String {
    let mut chars = s.chars();
    let first = match chars.next() {
        Some(c) => c.to_ascii_uppercase(),
        None => return String::new(),
    };
    // Two-letter elements the typer cares about: Cl, Br.
    let rest: String = chars.collect::<String>().to_ascii_lowercase();
    let two = format!("{first}{rest}");
    match two.as_str() {
        "Cl" | "Br" => two,
        _ => first.to_string(),
    }
}

// --- Ring + aromaticity perception ------------------------------------

/// Finds, for every atom, whether it lies in a 5- or 6-membered ring.
///
/// A depth-bounded DFS from each atom looks for a cycle of length 5 or
/// 6 back to the start through distinct atoms — sufficient for the
/// organic rings the OPLS-AA subset covers (benzene, pyridine,
/// imidazole, furan, cyclopentane / cyclohexane).
fn ring_membership(conn: &Connectivity) -> Vec<bool> {
    let n = conn.len();
    let mut in_ring = vec![false; n];
    for start in 0..n {
        if in_ring[start] {
            continue;
        }
        // DFS path search for a 5- or 6-cycle.
        let mut path = vec![start];
        if dfs_ring(conn, start, start, &mut path, &mut in_ring) {
            // dfs_ring marks the cycle directly.
        }
    }
    in_ring
}

/// Recursive helper for [`ring_membership`]: extends `path` and, on
/// closing a 5- or 6-cycle back to `start`, marks every atom on it.
fn dfs_ring(
    conn: &Connectivity,
    start: usize,
    current: usize,
    path: &mut Vec<usize>,
    in_ring: &mut [bool],
) -> bool {
    if path.len() > 6 {
        return false;
    }
    for &next in conn.neighbors(current) {
        if next == start && path.len() >= 5 {
            // Closed a ring of size path.len().
            for &a in path.iter() {
                in_ring[a] = true;
            }
            return true;
        }
        if path.contains(&next) {
            continue;
        }
        path.push(next);
        if dfs_ring(conn, start, next, path, in_ring) {
            path.pop();
            return true;
        }
        path.pop();
    }
    false
}

/// Decides which ring atoms are *aromatic*.
///
/// A pragmatic Hückel-style rule for the subset: an atom is aromatic if
/// it is a ring carbon or ring nitrogen whose every heavy ring
/// neighbour is itself a ring C/N, and the ring carbons each carry
/// exactly one hydrogen or one substituent (degree-3, sp²-consistent).
/// This recognises benzene, pyridine, the azoles and fused aromatics
/// without a full ring-perception SSSR pass.
fn aromatic_atoms(conn: &Connectivity, in_ring: &[bool]) -> Vec<bool> {
    let n = conn.len();
    let mut arom = vec![false; n];
    for i in 0..n {
        if !in_ring[i] {
            continue;
        }
        let el = conn.element(i);
        if el != "C" && el != "N" {
            continue;
        }
        // A ring carbon in an aromatic ring is trigonal: total degree
        // 3 (two ring neighbours + one H/substituent), no attached
        // hydrogens beyond one.
        let deg = conn.degree(i);
        let candidate = match el {
            "C" => deg == 3 && conn.hydrogen_count(i) <= 1,
            // Pyridine-type N: 2 ring neighbours, no H. Pyrrole-type
            // N: 2 ring neighbours + 1 H/substituent.
            "N" => deg == 2 || (deg == 3 && conn.hydrogen_count(i) <= 1),
            _ => false,
        };
        if !candidate {
            continue;
        }
        // Every heavy ring neighbour must be a ring C or N.
        let ring_ok = conn.heavy_neighbors(i).iter().all(|&j| {
            !in_ring[j] || {
                let e = conn.element(j);
                e == "C" || e == "N" || e == "O" || e == "S"
            }
        });
        if ring_ok {
            arom[i] = true;
        }
    }
    // Aromaticity is a ring property: keep an atom aromatic only if it
    // has two aromatic ring neighbours (so an isolated false positive
    // on a non-aromatic ring atom is dropped).
    let confirmed: Vec<bool> = (0..n)
        .map(|i| {
            arom[i]
                && conn
                    .heavy_neighbors(i)
                    .iter()
                    .filter(|&&j| arom[j])
                    .count()
                    >= 2
        })
        .collect();
    confirmed
}

// --- Hybridization perception -----------------------------------------

/// Perceives the hybridization of every atom.
///
/// With no stored bond orders, hybridization is inferred from valence:
/// carbon with 4 σ-neighbours is sp³, with 3 is sp² (alkene / carbonyl
/// / aromatic), with 2 heavy neighbours and no hydrogens is sp
/// (alkyne); nitrogen and oxygen follow the analogous degree rules; an
/// aromatic-ring atom is always sp².
fn hybridizations(conn: &Connectivity, arom: &[bool]) -> Vec<Hybridization> {
    (0..conn.len())
        .map(|i| {
            if arom[i] {
                return Hybridization::Sp2;
            }
            let el = conn.element(i);
            let deg = conn.degree(i);
            match el {
                "H" => Hybridization::Terminal,
                "C" => match deg {
                    4 => Hybridization::Sp3,
                    3 => Hybridization::Sp2,
                    2 => Hybridization::Sp,
                    _ => Hybridization::Terminal,
                },
                "N" => match deg {
                    // sp³ amine N (3 single bonds) or ammonium (4).
                    3 | 4 => Hybridization::Sp3,
                    // sp² (amide / imine) — perceived as sp² when the
                    // N neighbours a carbonyl/sp² carbon, else sp³.
                    2 => {
                        if conn
                            .heavy_neighbors(i)
                            .iter()
                            .any(|&j| is_sp2_carbon(conn, arom, j))
                        {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        }
                    }
                    // nitrile N — one neighbour, that neighbour sp.
                    1 => Hybridization::Sp,
                    _ => Hybridization::Terminal,
                },
                "O" => match deg {
                    2 => Hybridization::Sp3, // ether / hydroxyl / water
                    1 => Hybridization::Sp2, // carbonyl O (=O)
                    _ => Hybridization::Terminal,
                },
                "S" => match deg {
                    2 => Hybridization::Sp3, // thiol / thioether
                    1 => Hybridization::Sp2,
                    _ => Hybridization::Terminal,
                },
                // Monovalent halogens.
                "F" | "Cl" | "Br" | "I" => Hybridization::Terminal,
                "P" => Hybridization::Sp3,
                _ => Hybridization::Terminal,
            }
        })
        .collect()
}

/// Whether atom `j` is an sp²-perceived carbon (degree 3 or aromatic).
fn is_sp2_carbon(conn: &Connectivity, arom: &[bool], j: usize) -> bool {
    conn.element(j) == "C" && (arom[j] || conn.degree(j) == 3)
}

/// Whether atom `j` is a carbonyl carbon: an sp² carbon double-bonded
/// to a terminal (degree-1) oxygen.
fn is_carbonyl_carbon(conn: &Connectivity, j: usize) -> bool {
    conn.element(j) == "C"
        && conn.degree(j) == 3
        && conn
            .heavy_neighbors(j)
            .iter()
            .any(|&o| conn.element(o) == "O" && conn.degree(o) == 1)
}

/// Whether atom `i` is bonded to an alcohol *hydroxyl* oxygen — a
/// divalent O carrying exactly one hydrogen and not itself part of a
/// carboxylic acid. Selects the OPLS-AA `opls_157` alcohol-carbon type.
fn bonded_to_hydroxyl_oxygen(conn: &Connectivity, i: usize) -> bool {
    conn.heavy_neighbors(i).iter().any(|&o| {
        conn.element(o) == "O"
            && conn.degree(o) == 2
            && conn.hydrogen_count(o) == 1
            && !is_carboxyl_oxygen(conn, o)
    })
}

// --- The typer --------------------------------------------------------

/// Perceives the OPLS-AA atom type of every atom in a topology.
///
/// Returns one [`AtomType`] per atom, parallel to `topology.atoms`.
///
/// # Errors
/// [`MdError::Invalid`] (wrapping a [`TypingError`]) if any atom's
/// bonded environment is outside the OPLS-AA subset's coverage — the
/// caller can then fall back to the generic force-field path.
pub fn perceive_types(topology: &Topology) -> Result<Vec<AtomType>> {
    let conn = Connectivity::from_topology(topology);
    let in_ring = ring_membership(&conn);
    let arom = aromatic_atoms(&conn, &in_ring);
    let hyb = hybridizations(&conn, &arom);

    let mut out = Vec::with_capacity(conn.len());
    for i in 0..conn.len() {
        let element = conn.element(i).to_string();
        if element.is_empty() {
            return Err(TypingError::MissingElement(i).into());
        }
        let opls_type = type_one(&conn, &arom, &hyb, i)?;
        out.push(AtomType {
            opls_type,
            hybridization: hyb[i],
            aromatic: arom[i],
        });
    }
    Ok(out)
}

/// Assigns the OPLS-AA type string of a single atom `i`.
fn type_one(
    conn: &Connectivity,
    arom: &[bool],
    hyb: &[Hybridization],
    i: usize,
) -> Result<String> {
    let el = conn.element(i);
    let unrec = |reason: &str| -> MdError {
        TypingError::Unrecognized {
            index: i,
            element: el.to_string(),
            reason: reason.to_string(),
        }
        .into()
    };
    match el {
        "C" => type_carbon(conn, arom, hyb, i).ok_or_else(|| unrec("carbon environment")),
        "H" => type_hydrogen(conn, arom, i).ok_or_else(|| unrec("hydrogen environment")),
        "O" => type_oxygen(conn, i).ok_or_else(|| unrec("oxygen environment")),
        "N" => type_nitrogen(conn, arom, i).ok_or_else(|| unrec("nitrogen environment")),
        "S" => type_sulfur(conn, i).ok_or_else(|| unrec("sulfur environment")),
        "F" => Some("opls_164".to_string()).ok_or_else(|| unrec("")),
        "Cl" => Some("opls_165".to_string()).ok_or_else(|| unrec("")),
        "Br" => Some("opls_722".to_string()).ok_or_else(|| unrec("")),
        _ => Err(unrec("element outside the OPLS-AA subset")),
    }
}

/// Types a carbon by hybridization + functional environment.
fn type_carbon(
    conn: &Connectivity,
    arom: &[bool],
    hyb: &[Hybridization],
    i: usize,
) -> Option<String> {
    // Carbonyl carbon (C=O): carboxylic-acid vs ketone/aldehyde/amide.
    if is_carbonyl_carbon(conn, i) {
        let heavy = conn.heavy_neighbors(i);
        let hydroxyl_o = heavy.iter().any(|&o| {
            conn.element(o) == "O" && conn.degree(o) == 2 && conn.hydrogen_count(o) == 1
        });
        if hydroxyl_o {
            return Some("opls_267".to_string()); // carboxylic-acid C
        }
        let amide_n = heavy
            .iter()
            .any(|&n| conn.element(n) == "N");
        if amide_n {
            return Some("opls_177".to_string()); // amide carbonyl C
        }
        return Some("opls_235".to_string()); // ketone/aldehyde C=O
    }
    match hyb[i] {
        Hybridization::Sp3 => {
            // OPLS-AA assigns a *distinct* type to an sp³ carbon bonded
            // to an alcohol hydroxyl oxygen (`opls_157`, "all-atom C:
            // CH3 & CH2, alcohol") — its partial charge differs from a
            // plain alkane carbon's because the polar O withdraws
            // density. This is a genuine OPLS-AA feature, not a
            // special-case: the methyl group of methanol is not a
            // standard alkane.
            if bonded_to_hydroxyl_oxygen(conn, i) {
                return Some("opls_157".to_string());
            }
            // Alkane carbon: classify by the number of *heavy*
            // neighbours (= 4 − #H), the OPLS-AA CH3/CH2/CH/C split.
            let heavy = conn.heavy_degree(i);
            match heavy {
                0 | 1 => Some("opls_135".to_string()), // CH4 / CH3-
                2 => Some("opls_136".to_string()),     // -CH2-
                3 => Some("opls_137".to_string()),     // >CH-
                4 => Some("opls_139".to_string()),     // >C<
                _ => None,
            }
        }
        Hybridization::Sp2 => {
            if arom[i] {
                Some("opls_145".to_string()) // aromatic C
            } else {
                Some("opls_141".to_string()) // alkene =C<
            }
        }
        Hybridization::Sp => Some("opls_754".to_string()), // alkyne / nitrile C
        Hybridization::Terminal => None,
    }
}

/// Types a hydrogen by the heavy atom it is attached to.
fn type_hydrogen(conn: &Connectivity, arom: &[bool], i: usize) -> Option<String> {
    let heavy = conn.heavy_neighbors(i);
    let host = *heavy.first()?;
    let host_el = conn.element(host);
    match host_el {
        "C" => {
            if arom[host] {
                Some("opls_146".to_string()) // aromatic H
            } else if conn.degree(host) == 3 {
                Some("opls_144".to_string()) // alkene =C-H
            } else if conn.degree(host) == 2 {
                Some("opls_146".to_string()) // alkyne C-H (shares the aromatic-H LJ/charge subset slot)
            } else if bonded_to_hydroxyl_oxygen(conn, host) {
                Some("opls_156".to_string()) // H on an alcohol carbon
            } else {
                Some("opls_140".to_string()) // aliphatic H-C
            }
        }
        // Hydroxyl / carboxyl / water hydrogen.
        "O" => {
            if conn.degree(host) == 2 && conn.hydrogen_count(host) == 2 {
                Some("opls_117".to_string()) // water H (TIP3P/SPC slot)
            } else if is_carboxyl_oxygen(conn, host) {
                Some("opls_268".to_string()) // carboxylic-acid O-H
            } else {
                Some("opls_155".to_string()) // alcohol O-H
            }
        }
        // Amine / amide hydrogen.
        "N" => Some("opls_240".to_string()),
        "S" => Some("opls_204".to_string()), // thiol S-H
        _ => None,
    }
}

/// Types an oxygen by its bonded environment.
fn type_oxygen(conn: &Connectivity, i: usize) -> Option<String> {
    let deg = conn.degree(i);
    let hcount = conn.hydrogen_count(i);
    match deg {
        // Terminal oxygen (=O): carbonyl / carboxyl / amide.
        1 => {
            let host = *conn.heavy_neighbors(i).first()?;
            if conn.element(host) != "C" {
                return None;
            }
            // Carboxylic-acid carbonyl O vs ketone/aldehyde/amide O.
            let host_has_oh = conn.heavy_neighbors(host).iter().any(|&o| {
                o != i
                    && conn.element(o) == "O"
                    && conn.degree(o) == 2
                    && conn.hydrogen_count(o) == 1
            });
            if host_has_oh {
                Some("opls_269".to_string()) // carboxyl C=O
            } else if conn
                .heavy_neighbors(host)
                .iter()
                .any(|&n| conn.element(n) == "N")
            {
                Some("opls_178".to_string()) // amide C=O
            } else {
                Some("opls_236".to_string()) // ketone/aldehyde C=O
            }
        }
        // Divalent oxygen: water / hydroxyl / ether / carboxyl OH.
        2 => {
            if hcount == 2 {
                Some("opls_111".to_string()) // water O (TIP3P)
            } else if hcount == 1 {
                if is_carboxyl_oxygen(conn, i) {
                    Some("opls_268_O".to_string()) // carboxylic-acid -O-H oxygen
                } else {
                    Some("opls_154".to_string()) // alcohol -O-H
                }
            } else {
                Some("opls_180".to_string()) // ether -O-
            }
        }
        _ => None,
    }
}

/// Whether oxygen `i` is the hydroxyl oxygen of a carboxylic acid: a
/// divalent O with one H whose carbon neighbour also bears a terminal
/// `=O`.
fn is_carboxyl_oxygen(conn: &Connectivity, i: usize) -> bool {
    if conn.element(i) != "O" || conn.degree(i) != 2 || conn.hydrogen_count(i) != 1 {
        return false;
    }
    conn.heavy_neighbors(i).iter().any(|&c| {
        conn.element(c) == "C"
            && conn
                .heavy_neighbors(c)
                .iter()
                .any(|&o| o != i && conn.element(o) == "O" && conn.degree(o) == 1)
    })
}

/// Types a nitrogen by its bonded environment.
fn type_nitrogen(conn: &Connectivity, arom: &[bool], i: usize) -> Option<String> {
    if arom[i] {
        return Some("opls_511".to_string()); // aromatic ring N (pyridine-type)
    }
    let deg = conn.degree(i);
    // Amide N: bonded to a carbonyl carbon.
    if conn
        .heavy_neighbors(i)
        .iter()
        .any(|&c| is_carbonyl_carbon(conn, c))
    {
        return Some("opls_238".to_string()); // amide N
    }
    match deg {
        // Nitrile N (one neighbour).
        1 => Some("opls_753".to_string()),
        // Amine N (sp³): 3 substituents.
        3 => Some("opls_900".to_string()),
        // Ammonium N (4 substituents).
        4 => Some("opls_287".to_string()),
        // Imine-type 2-coordinate N.
        2 => Some("opls_903".to_string()),
        _ => None,
    }
}

/// Types a sulfur by its bonded environment.
fn type_sulfur(conn: &Connectivity, i: usize) -> Option<String> {
    match conn.degree(i) {
        // Thiol S-H vs thioether / disulfide.
        2 => {
            if conn.hydrogen_count(i) == 1 {
                Some("opls_200".to_string()) // thiol S
            } else {
                Some("opls_202".to_string()) // thioether / disulfide S
            }
        }
        _ => None,
    }
}

/// A typed-atom summary keyed by index, for diagnostics / reporting.
pub fn type_map(topology: &Topology) -> Result<HashMap<usize, AtomType>> {
    Ok(perceive_types(topology)?
        .into_iter()
        .enumerate()
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::Atom;

    /// Builds a topology from `(element, n_bonds_placeholder)` plus an
    /// explicit bond list. Masses/charges are dummies — typing ignores
    /// them.
    fn topo(elements: &[&str], bonds: &[(usize, usize)]) -> Topology {
        let mut t = Topology::new();
        for &e in elements {
            let mass = match e {
                "H" => 1.008,
                "C" => 12.011,
                "N" => 14.007,
                "O" => 15.999,
                "S" => 32.06,
                _ => 12.0,
            };
            t.push_atom(Atom::new(e, mass, 0.0).unwrap().with_element(e));
        }
        for &(i, j) in bonds {
            t.add_bond(i, j).unwrap();
        }
        t
    }

    #[test]
    fn ethane_carbons_are_sp3_alkane() {
        // C-C with 3 H each: CH3-CH3.
        let t = topo(
            &["C", "C", "H", "H", "H", "H", "H", "H"],
            &[
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (1, 5),
                (1, 6),
                (1, 7),
            ],
        );
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[0].opls_type, "opls_135"); // CH3 carbon
        assert_eq!(types[1].opls_type, "opls_135");
        assert_eq!(types[0].hybridization, Hybridization::Sp3);
        for ht in &types[2..8] {
            assert_eq!(ht.opls_type, "opls_140"); // aliphatic H
        }
    }

    #[test]
    fn propane_central_carbon_is_ch2() {
        // CH3-CH2-CH3.
        let mut elements = vec!["C", "C", "C"];
        elements.extend(std::iter::repeat_n("H", 8));
        let bonds = [
            (0, 1),
            (1, 2),
            (0, 3),
            (0, 4),
            (0, 5),
            (1, 6),
            (1, 7),
            (2, 8),
            (2, 9),
            (2, 10),
        ];
        let t = topo(&elements, &bonds);
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[0].opls_type, "opls_135"); // CH3
        assert_eq!(types[1].opls_type, "opls_136"); // CH2
        assert_eq!(types[2].opls_type, "opls_135"); // CH3
    }

    #[test]
    fn ethene_carbons_are_sp2_alkene() {
        // H2C=CH2 — each carbon has 2 H + 1 C, degree 3.
        let t = topo(
            &["C", "C", "H", "H", "H", "H"],
            &[(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)],
        );
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[0].hybridization, Hybridization::Sp2);
        assert_eq!(types[0].opls_type, "opls_141"); // alkene C
        assert_eq!(types[2].opls_type, "opls_144"); // alkene H
    }

    #[test]
    fn benzene_ring_is_aromatic() {
        // 6 ring carbons + 6 ring hydrogens.
        let mut elements: Vec<&str> = vec!["C"; 6];
        elements.extend(vec!["H"; 6]);
        let mut bonds = Vec::new();
        for k in 0..6 {
            bonds.push((k, (k + 1) % 6)); // ring
            bonds.push((k, 6 + k)); // C-H
        }
        let t = topo(&elements, &bonds);
        let types = perceive_types(&t).unwrap();
        for (c, ct) in types[..6].iter().enumerate() {
            assert!(ct.aromatic, "ring carbon {c} should be aromatic");
            assert_eq!(ct.opls_type, "opls_145"); // aromatic C
        }
        for ht in &types[6..12] {
            assert_eq!(ht.opls_type, "opls_146"); // aromatic H
        }
    }

    #[test]
    fn water_oxygen_and_hydrogens() {
        let t = topo(&["O", "H", "H"], &[(0, 1), (0, 2)]);
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[0].opls_type, "opls_111"); // water O
        assert_eq!(types[1].opls_type, "opls_117"); // water H
    }

    #[test]
    fn methanol_groups() {
        // CH3-OH: C(3H,O) O(C,H). OPLS-AA gives the alcohol methyl
        // carbon its own type (opls_157), not a plain alkane type.
        let t = topo(
            &["C", "O", "H", "H", "H", "H"],
            &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)],
        );
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[0].opls_type, "opls_157"); // alcohol methyl C
        assert_eq!(types[1].opls_type, "opls_154"); // hydroxyl O
        assert_eq!(types[2].opls_type, "opls_156"); // H on the alcohol C
        assert_eq!(types[5].opls_type, "opls_155"); // hydroxyl H
    }

    #[test]
    fn carbonyl_and_carboxyl_carbons_differ() {
        // Acetone-like: C(=O) ketone — central C bonded to O(term),
        // and two methyls. Just the carbonyl carbon + its O matter.
        // CH3-C(=O)-CH3
        let elements = vec![
            "C", "C", "C", "O", "H", "H", "H", "H", "H", "H",
        ];
        let bonds = [
            (0, 1),
            (1, 2),
            (1, 3), // C=O
            (0, 4),
            (0, 5),
            (0, 6),
            (2, 7),
            (2, 8),
            (2, 9),
        ];
        let t = topo(&elements, &bonds);
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[1].opls_type, "opls_235"); // ketone carbonyl C
        assert_eq!(types[3].opls_type, "opls_236"); // ketone carbonyl O

        // Acetic acid: CH3-C(=O)-O-H.
        let acid_el = vec!["C", "C", "O", "O", "H", "H", "H", "H"];
        let acid_bonds = [
            (0, 1),
            (1, 2), // C=O terminal
            (1, 3), // C-O hydroxyl
            (3, 4), // O-H
            (0, 5),
            (0, 6),
            (0, 7),
        ];
        let acid = topo(&acid_el, &acid_bonds);
        let at = perceive_types(&acid).unwrap();
        assert_eq!(at[1].opls_type, "opls_267"); // carboxyl C
        assert_eq!(at[2].opls_type, "opls_269"); // carboxyl C=O
        assert_eq!(at[3].opls_type, "opls_268_O"); // carboxyl -O-
        assert_eq!(at[4].opls_type, "opls_268"); // carboxyl O-H
    }

    #[test]
    fn amine_nitrogen_is_typed() {
        // Methylamine CH3-NH2.
        let t = topo(
            &["C", "N", "H", "H", "H", "H", "H"],
            &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6)],
        );
        let types = perceive_types(&t).unwrap();
        assert_eq!(types[1].opls_type, "opls_900"); // amine N
        assert_eq!(types[5].opls_type, "opls_240"); // amine H
    }

    #[test]
    fn unrecognized_environment_errors() {
        // A bare boron atom — outside the subset.
        let t = topo(&["B"], &[]);
        let err = perceive_types(&t).unwrap_err();
        assert_eq!(err.category(), "input");
    }

    #[test]
    fn element_fallback_from_type_name() {
        // No element field — must fall back to the type-name prefix.
        let mut t = Topology::new();
        t.push_atom(Atom::new("CT", 12.011, 0.0).unwrap());
        t.push_atom(Atom::new("HC", 1.008, 0.0).unwrap());
        t.add_bond(0, 1).unwrap();
        let conn = Connectivity::from_topology(&t);
        assert_eq!(conn.element(0), "C");
        assert_eq!(conn.element(1), "H");
    }
}

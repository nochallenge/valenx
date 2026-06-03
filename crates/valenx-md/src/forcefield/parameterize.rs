//! Applying a real force field — atom typing → parameter lookup.
//!
//! [`parameterize`] is the path that turns a bare molecular system —
//! atoms + bonds, no numeric constants — into a fully parameterised
//! one ready to simulate, exactly as GROMACS's `pdb2gmx` or OpenMM's
//! `ForceField.createSystem` do:
//!
//! 1. **Atom typing** — [`typing::perceive_types`] assigns an OPLS-AA
//!    atom type to every atom from its element + bonded connectivity +
//!    perceived hybridization.
//! 2. **Bonded-term generation** — angles are enumerated from every
//!    bonded `i`-`j`-`k` path and proper dihedrals from every `i`-`j`-
//!    `k`-`l` path, because a topology that carries only *bonds* still
//!    has a fully-determined angle / torsion set. Improper dihedrals
//!    are added on every sp²-perceived centre to hold its planarity.
//! 3. **Parameter lookup** — each atom's LJ σ/ε and partial charge,
//!    and each generated bonded term's constants, are pulled from the
//!    [`oplsaa`] database keyed by atom type.
//!
//! The result is a [`Parameterized`]: a new [`System`] whose topology
//! carries the generated angles / dihedrals / impropers and whose
//! atoms carry the OPLS-AA partial charges, plus the matching
//! [`ForceField`] (LJ table + bonded parameters positionally parallel
//! to the new topology). It plugs straight into
//! [`Simulation::new`](crate::sim::Simulation::new).
//!
//! [`typing`]: crate::forcefield::typing
//! [`oplsaa`]: crate::forcefield::oplsaa
//!
//! ## Honest scope
//!
//! `parameterize` is the **real-force-field default** for any molecule
//! the OPLS-AA subset can type. A molecule it cannot type — an element
//! or functional group outside the encoded coverage — yields an
//! [`MdError`] from the typer or a "no parameter" error from the
//! lookup, and the caller is expected to fall back to the **generic
//! path** ([`ForceField::new`] + manual `set_lj` / `push_bond`), which
//! stays fully supported.
//!
//! The generated angle / dihedral set is the *complete graph-implied*
//! set; OPLS-AA's own topology files sometimes drop a redundant
//! torsion, but including the full set is correct (a torsion the force
//! field has no parameter for is simply skipped — see
//! [`ParameterizeOptions::strict`]).

use crate::error::{MdError, Result};
use crate::forcefield::typing::{self, AtomType, Connectivity};
use crate::forcefield::{oplsaa, CombiningRule, ForceField};
use crate::system::{System, Topology};

/// Options controlling how [`parameterize`] handles gaps.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ParameterizeOptions {
    /// If `true`, a missing bonded parameter (a bond / angle / torsion
    /// the database does not cover) is a hard [`MdError`]. If `false`
    /// (the default), a missing *torsion* is silently dropped — the
    /// force field has no term for it — while a missing *bond* or
    /// *angle* is always an error (those are structural).
    pub strict: bool,
    /// Whether to add improper-dihedral planarity restraints on
    /// sp²-perceived centres (default `true`).
    pub impropers: bool,
}

impl Default for ParameterizeOptions {
    fn default() -> Self {
        ParameterizeOptions {
            strict: false,
            impropers: true,
        }
    }
}

/// The output of [`parameterize`]: a fully-parameterised system + its
/// force field.
#[derive(Clone, Debug, PartialEq)]
pub struct Parameterized {
    /// The system with generated angles / dihedrals / impropers and
    /// OPLS-AA partial charges written onto its atoms.
    pub system: System,
    /// The force field — LJ table + bonded parameter lists parallel to
    /// `system.topology`.
    pub force_field: ForceField,
    /// The perceived atom type of every atom, parallel to
    /// `system.topology.atoms` — kept for reporting / inspection.
    pub atom_types: Vec<AtomType>,
}

impl Parameterized {
    /// The OPLS-AA atom-type string of atom `i`.
    pub fn type_of(&self, i: usize) -> Option<&str> {
        self.atom_types.get(i).map(|t| t.opls_type.as_str())
    }

    /// The total (net) charge assigned by the force field (e).
    pub fn net_charge(&self) -> f64 {
        self.system.topology.net_charge()
    }
}

/// Parameterises a system with the OPLS-AA subset force field.
///
/// Runs atom typing, generates the angle / dihedral / improper terms
/// implied by the bond graph, and looks up every nonbonded + bonded
/// parameter from the [`oplsaa`] database. Returns a [`Parameterized`]
/// bundle ready for [`Simulation::new`](crate::sim::Simulation::new).
///
/// The input system's *positions* and *box* are preserved verbatim;
/// only the topology's bonded-term lists and the atoms' charges are
/// (re)written.
///
/// # Errors
/// - [`MdError::Invalid`] if atom typing fails (an atom outside the
///   OPLS-AA subset's coverage).
/// - [`MdError::Invalid`] if a structural bond or angle parameter is
///   missing, or — under [`ParameterizeOptions::strict`] — if any
///   torsion parameter is missing.
pub fn parameterize(system: &System) -> Result<Parameterized> {
    parameterize_with(system, ParameterizeOptions::default())
}

/// [`parameterize`] with explicit [`ParameterizeOptions`].
pub fn parameterize_with(
    system: &System,
    options: ParameterizeOptions,
) -> Result<Parameterized> {
    let topo = &system.topology;
    // 1. Atom typing.
    let atom_types = typing::perceive_types(topo)?;
    let conn = Connectivity::from_topology(topo);

    // 2. Build the new topology: the same atoms (with charges
    //    rewritten) + the same bonds + generated angles / dihedrals /
    //    impropers.
    let mut new_topo = Topology::new();
    for (i, atom) in topo.atoms.iter().enumerate() {
        let t = &atom_types[i].opls_type;
        let nb = oplsaa::atom(t).ok_or_else(|| {
            MdError::invalid(
                "forcefield",
                format!("OPLS-AA type `{t}` for atom {i} has no nonbonded parameters"),
            )
        })?;
        let mut new_atom = atom.clone();
        new_atom.type_name = t.clone();
        new_atom.charge = nb.charge;
        new_topo.push_atom(new_atom);
    }
    // Bonds: copy verbatim.
    for b in &topo.bonds {
        new_topo.add_bond(b.i, b.j)?;
    }
    // Angles: every bonded i-j-k path (vertex j), each once.
    for (i, j, k) in enumerate_angles(&conn) {
        new_topo.add_angle(i, j, k)?;
    }
    // Proper dihedrals: every bonded i-j-k-l path, each once.
    let dihedral_paths = enumerate_dihedrals(&conn);
    for &(i, j, k, l) in &dihedral_paths {
        new_topo.add_dihedral(i, j, k, l)?;
    }
    // Impropers: one per sp2-perceived centre.
    let improper_centers = if options.impropers {
        enumerate_impropers(&conn, &atom_types)
    } else {
        Vec::new()
    };
    for &(c, a, b, d) in &improper_centers {
        new_topo.add_improper(c, a, b, d)?;
    }

    // 3. Build the force field — LJ table + bonded parameters parallel
    //    to the new topology's lists.
    let mut ff = ForceField::new(CombiningRule::Geometric); // OPLS-AA comb-rule 3
    // OPLS-AA 1-4 scaling: LJ and Coulomb both 0.5.
    ff.lj_14_scale = 0.5;
    ff.coulomb_14_scale = 0.5;

    // LJ table: one entry per distinct atom type.
    for at in &atom_types {
        let lj = oplsaa::lj(&at.opls_type).ok_or_else(|| {
            MdError::invalid(
                "forcefield",
                format!("OPLS-AA type `{}` has no LJ parameters", at.opls_type),
            )
        })?;
        ff.set_lj(at.opls_type.clone(), lj);
    }

    let type_of = |idx: usize| atom_types[idx].opls_type.as_str();

    // Bond parameters, in topology order.
    for b in &new_topo.bonds {
        let param = oplsaa::bond(type_of(b.i), type_of(b.j)).ok_or_else(|| {
            MdError::invalid(
                "forcefield",
                format!(
                    "no OPLS-AA bond parameter for {}-{} (atoms {}-{})",
                    type_of(b.i),
                    type_of(b.j),
                    b.i,
                    b.j
                ),
            )
        })?;
        ff.push_bond(param);
    }
    // Angle parameters, in topology order.
    for a in &new_topo.angles {
        let param = oplsaa::angle(type_of(a.i), type_of(a.j), type_of(a.k)).ok_or_else(|| {
            MdError::invalid(
                "forcefield",
                format!(
                    "no OPLS-AA angle parameter for {}-{}-{} (atoms {}-{}-{})",
                    type_of(a.i),
                    type_of(a.j),
                    type_of(a.k),
                    a.i,
                    a.j,
                    a.k
                ),
            )
        })?;
        ff.push_angle(param);
    }
    // Dihedral parameters: a torsion the database does not cover is
    // dropped from the topology (non-strict) — so the topology and the
    // force field stay length-parallel. Rebuild the dihedral list.
    let mut kept_dihedrals = Vec::new();
    for &(i, j, k, l) in &dihedral_paths {
        match oplsaa::proper_dihedral(type_of(i), type_of(j), type_of(k), type_of(l)) {
            Some(param) => {
                ff.push_dihedral(param);
                kept_dihedrals.push((i, j, k, l));
            }
            None => {
                if options.strict {
                    return Err(MdError::invalid(
                        "forcefield",
                        format!(
                            "no OPLS-AA torsion for {}-{}-{}-{}",
                            type_of(i),
                            type_of(j),
                            type_of(k),
                            type_of(l)
                        ),
                    ));
                }
                // else: silently drop — see ParameterizeOptions::strict.
            }
        }
    }
    // Improper parameters: one per centre.
    let mut kept_impropers = Vec::new();
    for &(c, a, b, d) in &improper_centers {
        if let Some(param) = oplsaa::improper(type_of(c)) {
            ff.push_improper(param);
            kept_impropers.push((c, a, b, d));
        }
    }

    // If non-strict dropped torsions/impropers, rebuild the topology's
    // dihedral / improper lists so they stay parallel to the FF.
    if kept_dihedrals.len() != new_topo.dihedrals.len()
        || kept_impropers.len() != new_topo.impropers.len()
    {
        new_topo.dihedrals.clear();
        new_topo.impropers.clear();
        for &(i, j, k, l) in &kept_dihedrals {
            new_topo.add_dihedral(i, j, k, l)?;
        }
        for &(c, a, b, d) in &kept_impropers {
            new_topo.add_improper(c, a, b, d)?;
        }
    }

    // Assemble the parameterised system, preserving positions + box +
    // velocities.
    let mut new_system = System::new(new_topo, system.positions.clone())?;
    new_system.cell = system.cell.clone();
    if system.velocities.len() == new_system.len() {
        new_system.velocities = system.velocities.clone();
    }

    // Final consistency check.
    ff.validate_against(&new_system.topology)?;

    Ok(Parameterized {
        system: new_system,
        force_field: ff,
        atom_types,
    })
}

/// Whether the OPLS-AA subset can fully type *and* parameterise a
/// system — a cheap pre-check so a caller can decide between the real
/// and the generic path without catching an error.
pub fn can_parameterize(system: &System) -> bool {
    parameterize_with(
        system,
        ParameterizeOptions {
            strict: false,
            impropers: true,
        },
    )
    .is_ok()
}

// --- Bonded-term enumeration ------------------------------------------

/// Enumerates every angle `i`-`j`-`k` (vertex `j`) implied by the bond
/// graph — each unordered `{i,k}` pair around a vertex once.
fn enumerate_angles(conn: &Connectivity) -> Vec<(usize, usize, usize)> {
    let mut out = Vec::new();
    for j in 0..conn.len() {
        let nbrs = conn.neighbors(j);
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let i = nbrs[a];
                let k = nbrs[b];
                // Canonical: smaller leg first.
                if i < k {
                    out.push((i, j, k));
                } else {
                    out.push((k, j, i));
                }
            }
        }
    }
    out
}

/// Enumerates every proper dihedral `i`-`j`-`k`-`l` implied by the
/// bond graph — each torsion about a distinct central `j`-`k` bond,
/// counted once (canonicalised so `j < k`, and `(i,l)` deterministic).
fn enumerate_dihedrals(conn: &Connectivity) -> Vec<(usize, usize, usize, usize)> {
    let mut out = Vec::new();
    for j in 0..conn.len() {
        for &k in conn.neighbors(j) {
            if k <= j {
                continue; // central bond once, as j<k
            }
            for &i in conn.neighbors(j) {
                if i == k {
                    continue;
                }
                for &l in conn.neighbors(k) {
                    if l == j || l == i {
                        continue;
                    }
                    out.push((i, j, k, l));
                }
            }
        }
    }
    out
}

/// Enumerates the improper-dihedral centres: every sp²-perceived atom
/// with exactly three heavy/total neighbours. Returns
/// `(center, n1, n2, n3)` — the centre plus its three neighbours, the
/// `i`-`j`-`k`-`l` convention the improper term expects with the
/// constrained atom first.
fn enumerate_impropers(
    conn: &Connectivity,
    types: &[AtomType],
) -> Vec<(usize, usize, usize, usize)> {
    let mut out = Vec::new();
    for (c, atom_type) in types.iter().enumerate() {
        // Only a trigonal (degree-3) sp2 centre gets a planarity
        // improper.
        if atom_type.hybridization != typing::Hybridization::Sp2 {
            continue;
        }
        let nbrs = conn.neighbors(c);
        if nbrs.len() != 3 {
            continue;
        }
        // Only place the improper if the database has a constant for
        // this centre type — keeps the count honest.
        if oplsaa::improper(&atom_type.opls_type).is_none() {
            continue;
        }
        out.push((c, nbrs[0], nbrs[1], nbrs[2]));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::Atom;
    use nalgebra::Vector3;

    /// Builds a small molecule topology + a non-periodic system from an
    /// element list, bond list and rough coordinates.
    fn molecule(
        elements: &[&str],
        bonds: &[(usize, usize)],
        positions: &[[f64; 3]],
    ) -> System {
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
        let pos = positions
            .iter()
            .map(|p| Vector3::new(p[0], p[1], p[2]))
            .collect();
        System::new(t, pos).unwrap()
    }

    /// An ethane molecule with a roughly correct geometry (nm).
    fn ethane() -> System {
        molecule(
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
            &[
                [0.0, 0.0, 0.0],
                [0.153, 0.0, 0.0],
                [-0.036, 0.103, 0.0],
                [-0.036, -0.051, 0.089],
                [-0.036, -0.051, -0.089],
                [0.189, -0.103, 0.0],
                [0.189, 0.051, 0.089],
                [0.189, 0.051, -0.089],
            ],
        )
    }

    #[test]
    fn ethane_parameterizes_with_real_oplsaa() {
        let sys = ethane();
        let p = parameterize(&sys).unwrap();
        // Carbons typed as alkane CT (opls_135).
        assert_eq!(p.type_of(0), Some("opls_135"));
        assert_eq!(p.type_of(1), Some("opls_135"));
        // Hydrogens typed opls_140.
        for h in 2..8 {
            assert_eq!(p.type_of(h), Some("opls_140"));
        }
        // Charges from the database: CT -0.18, HC +0.06; molecule
        // neutral.
        assert!((p.system.topology.atoms[0].charge - (-0.18)).abs() < 1e-9);
        assert!((p.system.topology.atoms[2].charge - 0.06).abs() < 1e-9);
        assert!(p.net_charge().abs() < 1e-9, "net q = {}", p.net_charge());
    }

    #[test]
    fn ethane_generates_the_right_bonded_terms() {
        let sys = ethane();
        let p = parameterize(&sys).unwrap();
        let t = &p.system.topology;
        // 7 bonds (C-C + 6 C-H).
        assert_eq!(t.bonds.len(), 7);
        // Angles: each carbon centres C(3 H)(1 C) -> the degree-4
        // vertex gives C(4,2)=6 angles, two carbons -> 12 angles.
        assert_eq!(t.angles.len(), 12);
        // Dihedrals: 3 H on C0 x 3 H on C1 = 9 H-C-C-H torsions.
        assert_eq!(t.dihedrals.len(), 9);
        // Ethane has no sp2 centre -> no impropers.
        assert_eq!(t.impropers.len(), 0);
    }

    #[test]
    fn ethane_bond_and_angle_params_are_published() {
        let sys = ethane();
        let p = parameterize(&sys).unwrap();
        // The C-C bond (bond 0) must be the published CT-CT.
        let cc = p.force_field.bonds()[0];
        assert!((cc.r0 - 0.1529).abs() < 1e-9, "C-C r0 = {}", cc.r0);
        // Every angle constant should be a real OPLS-AA value (>0).
        for a in p.force_field.angles() {
            assert!(a.k > 0.0);
        }
    }

    #[test]
    fn ethane_builds_a_runnable_simulation() {
        let mut sys = ethane();
        sys.cell = SimBox::none();
        let p = parameterize(&sys).unwrap();
        // The parameterised system + force field plug into Simulation.
        let mut sim = crate::sim::Simulation::new(p.system, p.force_field).unwrap();
        let report = sim.run(20).unwrap();
        assert!(report.final_potential_energy.is_finite());
        assert!(report.final_total_energy.is_finite());
    }

    #[test]
    fn water_parameterizes_to_tip3p() {
        let sys = molecule(
            &["O", "H", "H"],
            &[(0, 1), (0, 2)],
            &[
                [0.0, 0.0, 0.0],
                [0.0957, 0.0, 0.0],
                [-0.0240, 0.0927, 0.0],
            ],
        );
        let p = parameterize(&sys).unwrap();
        assert_eq!(p.type_of(0), Some("opls_111"));
        assert_eq!(p.type_of(1), Some("opls_117"));
        // TIP3P charges.
        assert!((p.system.topology.atoms[0].charge - (-0.834)).abs() < 1e-9);
        assert!((p.system.topology.atoms[1].charge - 0.417).abs() < 1e-9);
        assert!(p.net_charge().abs() < 1e-9);
        // One angle (H-O-H), no dihedrals.
        assert_eq!(p.system.topology.angles.len(), 1);
        assert_eq!(p.system.topology.dihedrals.len(), 0);
    }

    #[test]
    fn methanol_parameterizes() {
        // CH3-OH.
        let sys = molecule(
            &["C", "O", "H", "H", "H", "H"],
            &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)],
            &[
                [0.0, 0.0, 0.0],
                [0.141, 0.0, 0.0],
                [-0.036, 0.103, 0.0],
                [-0.036, -0.051, 0.089],
                [-0.036, -0.051, -0.089],
                [0.176, 0.090, 0.0],
            ],
        );
        let p = parameterize(&sys).unwrap();
        assert_eq!(p.type_of(0), Some("opls_157")); // alcohol methyl C
        assert_eq!(p.type_of(1), Some("opls_154")); // hydroxyl O
        assert_eq!(p.type_of(5), Some("opls_155")); // hydroxyl H
        // OPLS-AA methanol charges sum to zero: +0.145 + 3(+0.040)
        // − 0.683 + 0.418 = 0.
        assert!(p.net_charge().abs() < 1e-6, "net q = {}", p.net_charge());
    }

    #[test]
    fn benzene_gets_aromatic_types_and_impropers() {
        let mut elements: Vec<&str> = vec!["C"; 6];
        elements.extend(vec!["H"; 6]);
        let mut bonds = Vec::new();
        let mut pos = Vec::new();
        let r = 0.139;
        for k in 0..6 {
            bonds.push((k, (k + 1) % 6));
            bonds.push((k, 6 + k));
            let ang = k as f64 * std::f64::consts::FRAC_PI_3;
            pos.push([r * ang.cos(), r * ang.sin(), 0.0]);
        }
        for k in 0..6 {
            let ang = k as f64 * std::f64::consts::FRAC_PI_3;
            let rh = 0.247;
            pos.push([rh * ang.cos(), rh * ang.sin(), 0.0]);
        }
        let sys = molecule(&elements, &bonds, &pos);
        let p = parameterize(&sys).unwrap();
        for c in 0..6 {
            assert_eq!(p.type_of(c), Some("opls_145"));
            assert!(p.atom_types[c].aromatic);
        }
        // Each aromatic carbon (degree 3, sp2) gets a planarity
        // improper -> 6 impropers.
        assert_eq!(p.system.topology.impropers.len(), 6);
    }

    #[test]
    fn unparameterizable_molecule_errors_and_can_check_is_false() {
        // A bare boron atom — outside OPLS-AA.
        let sys = molecule(&["B"], &[], &[[0.0, 0.0, 0.0]]);
        assert!(parameterize(&sys).is_err());
        assert!(!can_parameterize(&sys));
        // Ethane, by contrast, is parameterizable.
        assert!(can_parameterize(&ethane()));
    }

    #[test]
    fn positions_and_box_are_preserved() {
        let mut sys = ethane();
        sys.cell = SimBox::cubic(3.0).unwrap();
        let original_pos = sys.positions.clone();
        let p = parameterize(&sys).unwrap();
        assert_eq!(p.system.positions, original_pos);
        assert!(p.system.cell.is_periodic());
    }
}

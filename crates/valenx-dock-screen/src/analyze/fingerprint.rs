//! Feature 20 — protein-ligand interaction fingerprint.
//!
//! Two docked poses can have similar scores yet make completely
//! different contacts. An *interaction fingerprint* (IFP) captures
//! *which* interactions a pose makes — a structural signature you can
//! compare between poses, between ligands, or against a reference.
//!
//! [`interaction_fingerprint`] classifies every receptor / ligand
//! atom pair into one of five interaction types:
//!
//! - **hydrogen bond** — a donor and an acceptor within an H-bond
//!   distance;
//! - **hydrophobic contact** — two apolar atoms (carbon / halogen)
//!   in van der Waals contact;
//! - **π-stacking** — two aromatic carbons close enough to stack;
//! - **salt bridge** — two oppositely-charged atoms within range;
//! - **halogen bond** — a ligand halogen near a receptor acceptor.
//!
//! The result is an [`InteractionFingerprint`]: a per-type contact
//! list plus per-type counts, with a Tanimoto similarity for comparing
//! two fingerprints.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::receptor::Receptor;

/// The kind of a single protein-ligand interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum InteractionKind {
    /// A donor / acceptor hydrogen bond.
    HydrogenBond,
    /// An apolar van der Waals (hydrophobic) contact.
    Hydrophobic,
    /// An aromatic π-stacking contact.
    PiStacking,
    /// An oppositely-charged salt bridge.
    SaltBridge,
    /// A halogen bond (ligand halogen ··· receptor acceptor).
    HalogenBond,
}

impl InteractionKind {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            InteractionKind::HydrogenBond => "hbond",
            InteractionKind::Hydrophobic => "hydrophobic",
            InteractionKind::PiStacking => "pi_stacking",
            InteractionKind::SaltBridge => "salt_bridge",
            InteractionKind::HalogenBond => "halogen_bond",
        }
    }

    /// All five interaction kinds, in a fixed order.
    pub fn all() -> [InteractionKind; 5] {
        [
            InteractionKind::HydrogenBond,
            InteractionKind::Hydrophobic,
            InteractionKind::PiStacking,
            InteractionKind::SaltBridge,
            InteractionKind::HalogenBond,
        ]
    }
}

/// One detected protein-ligand interaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interaction {
    /// The interaction type.
    pub kind: InteractionKind,
    /// Receptor-atom index involved.
    pub receptor_atom: usize,
    /// Ligand-atom index involved.
    pub ligand_atom: usize,
    /// Centre-to-centre distance of the contact (Å).
    pub distance: f64,
}

/// A protein-ligand interaction fingerprint for one pose.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InteractionFingerprint {
    /// Every detected interaction.
    pub interactions: Vec<Interaction>,
}

impl InteractionFingerprint {
    /// Number of interactions of a given kind.
    pub fn count(&self, kind: InteractionKind) -> usize {
        self.interactions.iter().filter(|i| i.kind == kind).count()
    }

    /// Total interaction count.
    pub fn total(&self) -> usize {
        self.interactions.len()
    }

    /// The five per-kind counts in [`InteractionKind::all`] order — a
    /// compact 5-element fingerprint vector.
    pub fn count_vector(&self) -> [usize; 5] {
        let mut v = [0usize; 5];
        for (i, k) in InteractionKind::all().iter().enumerate() {
            v[i] = self.count(*k);
        }
        v
    }

    /// The set of `(receptor_atom, ligand_atom, kind)` interaction
    /// keys — used for the bit-style Tanimoto similarity.
    fn key_set(&self) -> std::collections::BTreeSet<(usize, usize, InteractionKind)> {
        self.interactions
            .iter()
            .map(|i| (i.receptor_atom, i.ligand_atom, i.kind))
            .collect()
    }

    /// Tanimoto similarity to another fingerprint: `|A ∩ B| / |A ∪ B|`
    /// over the interaction-key sets. `1.0` for identical fingerprints,
    /// `0.0` for disjoint, and `1.0` when both are empty.
    pub fn tanimoto(&self, other: &InteractionFingerprint) -> f64 {
        let a = self.key_set();
        let b = other.key_set();
        if a.is_empty() && b.is_empty() {
            return 1.0;
        }
        let inter = a.intersection(&b).count();
        let union = a.union(&b).count();
        if union == 0 {
            0.0
        } else {
            inter as f64 / union as f64
        }
    }
}

/// Distance cutoffs (Å) for each interaction type.
mod cutoffs {
    /// Donor ··· acceptor hydrogen-bond cutoff.
    pub const HBOND: f64 = 3.5;
    /// Apolar van der Waals contact cutoff.
    pub const HYDROPHOBIC: f64 = 4.5;
    /// Aromatic-carbon π-stacking cutoff.
    pub const PI_STACKING: f64 = 5.5;
    /// Charged-atom salt-bridge cutoff.
    pub const SALT_BRIDGE: f64 = 5.0;
    /// Halogen ··· acceptor halogen-bond cutoff.
    pub const HALOGEN: f64 = 4.0;
}

/// `true` if the AD4 type is a halogen.
fn is_halogen(t: Ad4AtomType) -> bool {
    matches!(
        t,
        Ad4AtomType::F | Ad4AtomType::Cl | Ad4AtomType::Br | Ad4AtomType::I
    )
}

/// Feature 20 — compute the interaction fingerprint of a posed ligand
/// against a receptor.
///
/// `ligand_atoms` is the *posed* ligand: `(world position, AD4 type,
/// partial charge)` per atom. Each receptor / ligand atom pair within
/// the largest cutoff is classified; a pair may register more than one
/// interaction (e.g. an H-bond and a salt bridge).
pub fn interaction_fingerprint(
    receptor: &Receptor,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType, f64)],
) -> InteractionFingerprint {
    let mut fp = InteractionFingerprint::default();
    let max_cut = cutoffs::PI_STACKING; // the largest of the five
    let max_cut2 = max_cut * max_cut;

    for (li, &(lp, lt, lq)) in ligand_atoms.iter().enumerate() {
        let lp_props = lt.props();
        for (ri, ra) in receptor.atoms.iter().enumerate() {
            let d2 = (lp - ra.position).norm_squared();
            if d2 > max_cut2 {
                continue;
            }
            let d = d2.sqrt();
            let rp = ra.ad4_type.props();
            let rq = ra.partial_charge;

            // --- hydrogen bond -----------------------------------
            let donor_acceptor =
                (lp_props.is_donor && rp.is_acceptor) || (lp_props.is_acceptor && rp.is_donor);
            if donor_acceptor && d <= cutoffs::HBOND {
                fp.interactions.push(Interaction {
                    kind: InteractionKind::HydrogenBond,
                    receptor_atom: ri,
                    ligand_atom: li,
                    distance: d,
                });
            }
            // --- hydrophobic contact -----------------------------
            if lp_props.is_hydrophobic && rp.is_hydrophobic && d <= cutoffs::HYDROPHOBIC {
                fp.interactions.push(Interaction {
                    kind: InteractionKind::Hydrophobic,
                    receptor_atom: ri,
                    ligand_atom: li,
                    distance: d,
                });
            }
            // --- π-stacking (two aromatic carbons) ---------------
            if lt == Ad4AtomType::A && ra.ad4_type == Ad4AtomType::A && d <= cutoffs::PI_STACKING {
                fp.interactions.push(Interaction {
                    kind: InteractionKind::PiStacking,
                    receptor_atom: ri,
                    ligand_atom: li,
                    distance: d,
                });
            }
            // --- salt bridge (opposite formal-ish charges) -------
            if lq * rq < -0.04 && d <= cutoffs::SALT_BRIDGE {
                fp.interactions.push(Interaction {
                    kind: InteractionKind::SaltBridge,
                    receptor_atom: ri,
                    ligand_atom: li,
                    distance: d,
                });
            }
            // --- halogen bond (ligand halogen ··· acceptor) ------
            if is_halogen(lt) && rp.is_acceptor && d <= cutoffs::HALOGEN {
                fp.interactions.push(Interaction {
                    kind: InteractionKind::HalogenBond,
                    receptor_atom: ri,
                    ligand_atom: li,
                    distance: d,
                });
            }
        }
    }
    fp
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    #[test]
    fn detects_a_hydrogen_bond() {
        // Receptor donor HD, ligand acceptor OA at 2.8 Å.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::HD,
                partial_charge: 0.2,
            }],
        };
        let lig = [(Vector3::new(2.8, 0.0, 0.0), Ad4AtomType::OA, -0.4)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.count(InteractionKind::HydrogenBond), 1);
    }

    #[test]
    fn detects_a_hydrophobic_contact() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let lig = [(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.count(InteractionKind::Hydrophobic), 1);
    }

    #[test]
    fn detects_pi_stacking_between_aromatic_carbons() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::A,
                partial_charge: 0.0,
            }],
        };
        let lig = [(Vector3::new(4.5, 0.0, 0.0), Ad4AtomType::A, 0.0)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.count(InteractionKind::PiStacking), 1);
    }

    #[test]
    fn detects_a_salt_bridge() {
        // A +0.6 ligand atom and a -0.6 receptor atom.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::OA,
                partial_charge: -0.6,
            }],
        };
        let lig = [(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::N, 0.6)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.count(InteractionKind::SaltBridge), 1);
    }

    #[test]
    fn detects_a_halogen_bond() {
        // A ligand chlorine near a receptor acceptor oxygen.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::OA,
                partial_charge: -0.3,
            }],
        };
        let lig = [(Vector3::new(3.5, 0.0, 0.0), Ad4AtomType::Cl, 0.0)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.count(InteractionKind::HalogenBond), 1);
    }

    #[test]
    fn far_apart_atoms_make_no_interactions() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let lig = [(Vector3::new(20.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert_eq!(fp.total(), 0);
    }

    #[test]
    fn count_vector_has_five_entries() {
        let fp = InteractionFingerprint::default();
        assert_eq!(fp.count_vector(), [0, 0, 0, 0, 0]);
    }

    #[test]
    fn tanimoto_of_identical_fingerprints_is_one() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let lig = [(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let fp = interaction_fingerprint(&receptor, &lig);
        assert!((fp.tanimoto(&fp) - 1.0).abs() < 1e-12);
        // Two empty fingerprints are also "identical".
        let empty = InteractionFingerprint::default();
        assert_eq!(empty.tanimoto(&empty), 1.0);
    }

    #[test]
    fn tanimoto_of_disjoint_fingerprints_is_zero() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let near = [(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let far = [(Vector3::new(20.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let fp_near = interaction_fingerprint(&receptor, &near);
        let fp_far = interaction_fingerprint(&receptor, &far);
        // fp_far has no interactions → disjoint from fp_near.
        assert_eq!(fp_near.tanimoto(&fp_far), 0.0);
    }
}

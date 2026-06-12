//! Pharmacophore feature perception.
//!
//! A pharmacophore is the set of abstract interaction points a
//! molecule presents to a binding site: hydrogen-bond donors and
//! acceptors, aromatic ring centroids, hydrophobic atoms, and
//! positively / negatively ionisable groups. [`pharmacophore`] walks a
//! molecule and emits a [`PharmacophoreFeature`] for each, with a 3D
//! position when the molecule has a conformer (else the 2D depiction
//! coordinate, else the origin).
//!
//! The feature points are what a pharmacophore-based screen aligns and
//! scores; [`feature_distances`] gives the inter-feature distance
//! matrix that forms a rotation-invariant pharmacophore fingerprint.
//!
//! **v1 scope.** Features are detected by element / hybridisation /
//! charge rules — donors are N-H / O-H, acceptors are N / O lone-pair
//! atoms, hydrophobic atoms are carbons / halogens with no polar
//! neighbour, aromatic features are SSSR-ring centroids. This matches
//! the standard RDKit feature-factory categories; it does not model
//! directional H-bond vectors or the metal-coordination feature.

use crate::coords::centroid as mol_centroid;
use crate::molecule::Molecule;
use crate::perceive::rings::sssr;

/// The kinds of pharmacophore feature this crate perceives.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FeatureKind {
    /// Hydrogen-bond donor (an N-H / O-H).
    Donor,
    /// Hydrogen-bond acceptor (an N / O lone pair).
    Acceptor,
    /// Aromatic ring centroid.
    Aromatic,
    /// Hydrophobic atom (apolar carbon / halogen).
    Hydrophobic,
    /// Positively ionisable group.
    PositiveIonizable,
    /// Negatively ionisable group.
    NegativeIonizable,
}

impl FeatureKind {
    /// A short stable label for the feature kind.
    pub fn label(self) -> &'static str {
        match self {
            FeatureKind::Donor => "donor",
            FeatureKind::Acceptor => "acceptor",
            FeatureKind::Aromatic => "aromatic",
            FeatureKind::Hydrophobic => "hydrophobic",
            FeatureKind::PositiveIonizable => "pos_ionizable",
            FeatureKind::NegativeIonizable => "neg_ionizable",
        }
    }
}

/// One perceived pharmacophore feature.
#[derive(Clone, Debug, PartialEq)]
pub struct PharmacophoreFeature {
    /// The feature category.
    pub kind: FeatureKind,
    /// Atom indices that constitute the feature (one atom for a donor /
    /// acceptor, the whole ring for an aromatic feature).
    pub atoms: Vec<usize>,
    /// 3D position of the feature point.
    pub position: [f64; 3],
}

/// Perceive the pharmacophore features of `mol`.
///
/// Positions come from [`Molecule::coords`] if present, else they are
/// all the origin (the feature *kinds* are still meaningful for a 2D
/// molecule).
pub fn pharmacophore(mol: &Molecule) -> Vec<PharmacophoreFeature> {
    let mut features = Vec::new();
    let pos = |i: usize| -> [f64; 3] { mol.coords.get(i).copied().unwrap_or([0.0, 0.0, 0.0]) };

    for (i, a) in mol.atoms.iter().enumerate() {
        if a.is_hydrogen() || a.is_dummy() {
            continue;
        }
        // donor / acceptor
        if matches!(a.atomic_number, 7 | 8) {
            features.push(PharmacophoreFeature {
                kind: FeatureKind::Acceptor,
                atoms: vec![i],
                position: pos(i),
            });
            if a.total_h() > 0 {
                features.push(PharmacophoreFeature {
                    kind: FeatureKind::Donor,
                    atoms: vec![i],
                    position: pos(i),
                });
            }
        }
        // ionisable groups
        if a.formal_charge > 0 || (a.atomic_number == 7 && is_basic_amine(mol, i)) {
            features.push(PharmacophoreFeature {
                kind: FeatureKind::PositiveIonizable,
                atoms: vec![i],
                position: pos(i),
            });
        }
        if a.formal_charge < 0 || is_carboxyl_oxygen(mol, i) {
            features.push(PharmacophoreFeature {
                kind: FeatureKind::NegativeIonizable,
                atoms: vec![i],
                position: pos(i),
            });
        }
        // hydrophobic atom
        if is_hydrophobic(mol, i) {
            features.push(PharmacophoreFeature {
                kind: FeatureKind::Hydrophobic,
                atoms: vec![i],
                position: pos(i),
            });
        }
    }

    // aromatic ring centroids
    let rings = sssr(mol);
    for ring in &rings.rings {
        if ring.atoms.iter().all(|&a| mol.atoms[a].aromatic) {
            let center = ring_centroid(mol, &ring.atoms);
            features.push(PharmacophoreFeature {
                kind: FeatureKind::Aromatic,
                atoms: ring.atoms.clone(),
                position: center,
            });
        }
    }
    features
}

/// Count features of a given kind — a coarse pharmacophore descriptor.
pub fn feature_count(mol: &Molecule, kind: FeatureKind) -> usize {
    pharmacophore(mol).iter().filter(|f| f.kind == kind).count()
}

/// The inter-feature distance matrix — a rotation-invariant
/// pharmacophore signature. `result[i][j]` is the distance between
/// feature `i` and feature `j`.
pub fn feature_distances(features: &[PharmacophoreFeature]) -> Vec<Vec<f64>> {
    let n = features.len();
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let p = features[i].position;
            let q = features[j].position;
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            m[i][j] = d;
            m[j][i] = d;
        }
    }
    m
}

fn ring_centroid(mol: &Molecule, atoms: &[usize]) -> [f64; 3] {
    if mol.coords.is_empty() || atoms.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let mut c = [0.0, 0.0, 0.0];
    for &a in atoms {
        let p = mol.coords.get(a).copied().unwrap_or([0.0, 0.0, 0.0]);
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    let n = atoms.len() as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Heuristic: is the nitrogen at `i` a basic (protonatable) amine —
/// an `sp³` N not adjacent to a carbonyl (so not an amide)?
fn is_basic_amine(mol: &Molecule, i: usize) -> bool {
    if mol.atoms[i].atomic_number != 7 || mol.atoms[i].aromatic {
        return false;
    }
    // not basic if bonded to a carbonyl carbon (amide nitrogen)
    for &nbr in &mol.neighbors(i) {
        if mol.atoms[nbr].atomic_number == 6 {
            let carbonyl = mol.bonds_on(nbr).iter().any(|&bi| {
                let b = &mol.bonds[bi];
                b.order == crate::molecule::BondOrder::Double
                    && b.other(nbr)
                        .map(|o| mol.atoms[o].atomic_number == 8)
                        .unwrap_or(false)
            });
            if carbonyl {
                return false;
            }
        }
    }
    true
}

/// Is the oxygen at `i` part of a carboxyl / carboxylate group (so a
/// negatively-ionisable site)?
fn is_carboxyl_oxygen(mol: &Molecule, i: usize) -> bool {
    if mol.atoms[i].atomic_number != 8 {
        return false;
    }
    for &nbr in &mol.neighbors(i) {
        if mol.atoms[nbr].atomic_number == 6 {
            // the carbon must also carry a second oxygen
            let oxygens = mol
                .neighbors(nbr)
                .iter()
                .filter(|&&o| mol.atoms[o].atomic_number == 8)
                .count();
            if oxygens >= 2 {
                return true;
            }
        }
    }
    false
}

/// Is the atom at `i` hydrophobic — an apolar carbon or a halogen with
/// no directly-bonded polar (N/O) atom?
fn is_hydrophobic(mol: &Molecule, i: usize) -> bool {
    let z = mol.atoms[i].atomic_number;
    if !matches!(z, 6 | 9 | 17 | 35 | 53) {
        return false;
    }
    // halogens are always hydrophobic
    if matches!(z, 9 | 17 | 35 | 53) {
        return true;
    }
    // a carbon is hydrophobic if no neighbour is N / O
    !mol.neighbors(i)
        .iter()
        .any(|&v| matches!(mol.atoms[v].atomic_number, 7 | 8))
}

/// Geometric extent of the pharmacophore — the maximum distance
/// between any two feature points. `0.0` if there are < 2 features.
pub fn pharmacophore_radius(features: &[PharmacophoreFeature]) -> f64 {
    let d = feature_distances(features);
    d.iter()
        .flat_map(|row| row.iter().copied())
        .fold(0.0, f64::max)
}

/// Centroid of all feature points.
pub fn feature_centroid(features: &[PharmacophoreFeature]) -> [f64; 3] {
    if features.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let n = features.len() as f64;
    let mut c = [0.0, 0.0, 0.0];
    for f in features {
        c[0] += f.position[0];
        c[1] += f.position[1];
        c[2] += f.position[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// A convenience that re-exports the molecule centroid for callers
/// comparing the pharmacophore centre to the molecular centre.
pub fn molecule_centroid(mol: &Molecule) -> [f64; 3] {
    mol_centroid(mol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn ethanol_has_donor_and_acceptor() {
        let m = mol_from_smiles("CCO").unwrap();
        let feats = pharmacophore(&m);
        assert_eq!(feature_count(&m, FeatureKind::Donor), 1);
        assert_eq!(feature_count(&m, FeatureKind::Acceptor), 1);
        assert!(feats.iter().any(|f| f.kind == FeatureKind::Hydrophobic));
    }

    #[test]
    fn benzene_has_aromatic_feature() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        assert_eq!(feature_count(&m, FeatureKind::Aromatic), 1);
        // no donors / acceptors in pure benzene
        assert_eq!(feature_count(&m, FeatureKind::Donor), 0);
    }

    #[test]
    fn carboxylic_acid_is_negative_ionizable() {
        let m = mol_from_smiles("CC(=O)O").unwrap();
        assert!(feature_count(&m, FeatureKind::NegativeIonizable) >= 1);
    }

    #[test]
    fn amine_is_positive_ionizable() {
        let m = mol_from_smiles("CCN").unwrap();
        assert_eq!(feature_count(&m, FeatureKind::PositiveIonizable), 1);
    }

    #[test]
    fn amide_nitrogen_is_not_basic() {
        // acetamide N is not a positive-ionizable feature
        let m = mol_from_smiles("CC(=O)N").unwrap();
        assert_eq!(feature_count(&m, FeatureKind::PositiveIonizable), 0);
    }

    #[test]
    fn feature_positions_from_conformer() {
        let m = mol_from_smiles("CCO").unwrap();
        let conf = crate::coords::embed_3d(&m, 1).unwrap();
        let feats = pharmacophore(&conf);
        // with a conformer, feature distances should be non-trivial
        let radius = pharmacophore_radius(&feats);
        assert!(radius > 0.0, "pharmacophore radius = {radius}");
    }

    #[test]
    fn distance_matrix_is_symmetric() {
        let m = mol_from_smiles("CCO").unwrap();
        let conf = crate::coords::embed_3d(&m, 2).unwrap();
        let feats = pharmacophore(&conf);
        let d = feature_distances(&feats);
        for i in 0..feats.len() {
            for j in 0..feats.len() {
                assert!((d[i][j] - d[j][i]).abs() < 1e-12);
            }
        }
    }
}

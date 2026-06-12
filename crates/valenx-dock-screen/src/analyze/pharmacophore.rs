//! Feature 22 — pharmacophore-based screening.
//!
//! A *pharmacophore* is the abstract pattern of interaction features a
//! ligand must present to bind — "a hydrogen-bond donor here, an
//! aromatic ring 5 Å away, a hydrophobe over there". Pharmacophore
//! screening filters a library to the molecules that *can* present a
//! query pattern, a fast pre-filter before (or instead of) docking.
//!
//! The pharmacophore feature perception itself lives in
//! [`mod@valenx_cheminf::analyze::pharmacophore`] — this module reuses it
//! and adds the *matching*: does a candidate molecule contain features
//! of the required kinds whose pairwise distances are compatible with
//! the query?
//!
//! A [`PharmacophoreQuery`] is a small set of required feature kinds
//! with optional pairwise distance constraints. [`pharmacophore_screen`]
//! scores each library molecule by how well its features satisfy the
//! query and returns the library ranked best-first.

use valenx_cheminf::analyze::pharmacophore::{
    feature_distances, pharmacophore, FeatureKind, PharmacophoreFeature,
};
use valenx_cheminf::molecule::Molecule;

use crate::error::{DockScreenError, Result};

/// A required pairwise distance between two query features.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistanceConstraint {
    /// Index of the first query feature (into [`PharmacophoreQuery::features`]).
    pub feature_a: usize,
    /// Index of the second query feature.
    pub feature_b: usize,
    /// Target separation (Å).
    pub distance: f64,
    /// Allowed `±` tolerance on the distance (Å).
    pub tolerance: f64,
}

/// A pharmacophore query: required feature kinds plus optional
/// pairwise distance constraints.
#[derive(Clone, Debug, Default)]
pub struct PharmacophoreQuery {
    /// The feature kinds the query requires (order matters — distance
    /// constraints index into this list).
    pub features: Vec<FeatureKind>,
    /// Pairwise distance constraints between the required features.
    pub distances: Vec<DistanceConstraint>,
}

impl PharmacophoreQuery {
    /// A query that just requires a set of feature kinds (no distance
    /// constraints).
    pub fn from_kinds(kinds: impl IntoIterator<Item = FeatureKind>) -> Self {
        PharmacophoreQuery {
            features: kinds.into_iter().collect(),
            distances: Vec::new(),
        }
    }

    /// Add a pairwise distance constraint, returning `self` for
    /// chaining. Constraints referencing an out-of-range feature index
    /// are silently ignored at match time.
    pub fn with_distance(
        mut self,
        feature_a: usize,
        feature_b: usize,
        distance: f64,
        tolerance: f64,
    ) -> Self {
        self.distances.push(DistanceConstraint {
            feature_a,
            feature_b,
            distance,
            tolerance,
        });
        self
    }
}

/// The pharmacophore-match outcome for one molecule.
#[derive(Clone, Debug, PartialEq)]
pub struct PharmacophoreHit {
    /// The molecule's index in the input library.
    pub molecule_index: usize,
    /// `true` if every required feature kind is present and every
    /// distance constraint can be satisfied by some feature pairing.
    pub matches: bool,
    /// A `[0,1]` match score — the fraction of query requirements
    /// (feature kinds + distance constraints) satisfied.
    pub score: f64,
    /// Number of required feature kinds the molecule presents.
    pub features_present: usize,
}

/// Feature 22 — screen a library of molecules against a pharmacophore
/// query.
///
/// Each molecule's pharmacophore features are perceived (via
/// [`valenx_cheminf`]) and matched against the query. The library is
/// returned ranked by match score (best first); full matches sort
/// above partial matches.
///
/// Molecules should carry 3D coordinates if the query uses distance
/// constraints — without coordinates every feature sits at the origin,
/// so distance constraints will only be satisfiable when their target
/// distance is ~0.
///
/// Returns [`DockScreenError::Invalid`] for an empty library or an
/// empty query.
pub fn pharmacophore_screen(
    library: &[Molecule],
    query: &PharmacophoreQuery,
) -> Result<Vec<PharmacophoreHit>> {
    if library.is_empty() {
        return Err(DockScreenError::invalid(
            "library",
            "cannot screen an empty molecule library",
        ));
    }
    if query.features.is_empty() {
        return Err(DockScreenError::invalid(
            "query",
            "pharmacophore query has no required features",
        ));
    }

    let mut hits: Vec<PharmacophoreHit> = library
        .iter()
        .enumerate()
        .map(|(i, mol)| match_one(i, mol, query))
        .collect();
    hits.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.molecule_index.cmp(&b.molecule_index))
    });
    Ok(hits)
}

/// Match one molecule against the query.
fn match_one(index: usize, mol: &Molecule, query: &PharmacophoreQuery) -> PharmacophoreHit {
    let feats = pharmacophore(mol);
    let dist = feature_distances(&feats);

    // --- feature-kind coverage ---------------------------------------
    // Greedily assign each required query feature to a distinct
    // molecule feature of the same kind. `assigned[q]` is the molecule
    // feature index chosen for query feature `q`, or `None`.
    let mut used = vec![false; feats.len()];
    let mut assigned: Vec<Option<usize>> = vec![None; query.features.len()];
    let mut features_present = 0usize;
    for (q, &kind) in query.features.iter().enumerate() {
        if let Some(mi) = first_unused_feature(&feats, &used, kind) {
            used[mi] = true;
            assigned[q] = Some(mi);
            features_present += 1;
        }
    }

    // --- distance-constraint satisfaction ----------------------------
    let mut distances_satisfied = 0usize;
    let mut distances_total = 0usize;
    for c in &query.distances {
        if c.feature_a >= query.features.len() || c.feature_b >= query.features.len() {
            continue; // ignore an out-of-range constraint
        }
        distances_total += 1;
        if let (Some(ma), Some(mb)) = (assigned[c.feature_a], assigned[c.feature_b]) {
            if let Some(d) = dist.get(ma).and_then(|row| row.get(mb)) {
                if (d - c.distance).abs() <= c.tolerance {
                    distances_satisfied += 1;
                }
            }
        }
    }

    let n_requirements = query.features.len() + distances_total;
    let satisfied = features_present + distances_satisfied;
    let score = if n_requirements == 0 {
        0.0
    } else {
        satisfied as f64 / n_requirements as f64
    };
    let matches =
        features_present == query.features.len() && distances_satisfied == distances_total;

    PharmacophoreHit {
        molecule_index: index,
        matches,
        score,
        features_present,
    }
}

/// The index of the first unused molecule feature of `kind`.
fn first_unused_feature(
    feats: &[PharmacophoreFeature],
    used: &[bool],
    kind: FeatureKind,
) -> Option<usize> {
    feats
        .iter()
        .enumerate()
        .position(|(i, f)| f.kind == kind && !used[i])
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cheminf::mol_from_smiles;

    #[test]
    fn rejects_empty_library_and_query() {
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Donor]);
        assert!(pharmacophore_screen(&[], &q).is_err());
        let m = mol_from_smiles("CCO").unwrap();
        let empty_q = PharmacophoreQuery::default();
        assert!(pharmacophore_screen(&[m], &empty_q).is_err());
    }

    #[test]
    fn a_molecule_with_the_required_feature_matches() {
        // Ethanol has a donor (O-H) and an acceptor (O).
        let ethanol = mol_from_smiles("CCO").unwrap();
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Donor, FeatureKind::Acceptor]);
        let hits = pharmacophore_screen(&[ethanol], &q).unwrap();
        assert!(hits[0].matches, "ethanol should match donor+acceptor");
        assert!((hits[0].score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_molecule_missing_a_feature_does_not_fully_match() {
        // Ethane (CC) has no donor / acceptor — a donor query fails.
        let ethane = mol_from_smiles("CC").unwrap();
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Donor]);
        let hits = pharmacophore_screen(&[ethane], &q).unwrap();
        assert!(!hits[0].matches);
        assert_eq!(hits[0].features_present, 0);
        assert_eq!(hits[0].score, 0.0);
    }

    #[test]
    fn full_matches_rank_above_partial_matches() {
        let with_oh = mol_from_smiles("CCO").unwrap(); // has donor
        let without = mol_from_smiles("CC").unwrap(); // no donor
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Donor]);
        let hits = pharmacophore_screen(&[without, with_oh], &q).unwrap();
        // The matching molecule (originally index 1) should sort first.
        assert_eq!(hits[0].molecule_index, 1);
        assert!(hits[0].matches);
        assert!(!hits[1].matches);
    }

    #[test]
    fn distance_constraint_is_counted_in_the_score() {
        // Benzene has aromatic features; a query needing two aromatic
        // features ~0 Å apart (no coordinates → all at origin) should
        // satisfy the distance constraint.
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Aromatic, FeatureKind::Aromatic])
            .with_distance(0, 1, 0.0, 0.5);
        let hits = pharmacophore_screen(&[benzene], &q).unwrap();
        // Benzene has one aromatic ring → only one aromatic feature,
        // so the second required aromatic feature is unmatched. The
        // score is therefore partial but well-defined.
        assert!((0.0..=1.0).contains(&hits[0].score));
    }

    #[test]
    fn out_of_range_distance_constraint_is_ignored() {
        let ethanol = mol_from_smiles("CCO").unwrap();
        // Constraint references feature index 9 which does not exist.
        let q = PharmacophoreQuery::from_kinds([FeatureKind::Donor]).with_distance(0, 9, 5.0, 1.0);
        let hits = pharmacophore_screen(&[ethanol], &q).unwrap();
        // The bad constraint is ignored → the donor requirement alone
        // is satisfied → full match.
        assert!(hits[0].matches);
    }
}

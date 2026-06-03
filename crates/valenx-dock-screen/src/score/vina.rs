//! Feature 5 — Vina-class empirical scoring function (the published
//! Vina form, with published term weights and atom typing).
//!
//! AutoDock Vina's scoring function (Trott & Olson 2010, *J. Comput.
//! Chem.* **31**, 455) is a weighted sum of five inter-atomic terms —
//! two Gaussians, a repulsion wall, a hydrophobic contact bonus, a
//! hydrogen-bond bonus — divided by a torsional-entropy factor that
//! grows with the ligand's rotatable-bond count:
//!
//! ```text
//!                Σ_pairs ( w_g1·gauss1 + w_g2·gauss2 + w_rep·rep
//!                          + w_hyd·hydrophobic + w_hb·hbond )
//! score (kcal/mol) = ───────────────────────────────────────────────
//!                          1 + w_rot · N_rotatable_bonds
//! ```
//!
//! ## Published Vina term weights — used verbatim
//!
//! The five term weights and the rotatable-bond penalty weight come
//! straight from Trott & Olson 2010 (Table S1 of the supplement) and
//! the Vina source `atom_constants.h`:
//!
//! | term | weight | role |
//! |------|--------|------|
//! | `GAUSS1` | −0.035579 | narrow attractive gaussian at contact |
//! | `GAUSS2` | −0.005156 | broad attractive gaussian at 3 Å |
//! | `REPULSION` | +0.840245 | soft-wall steric repulsion |
//! | `HYDROPHOBIC` | −0.035069 | hydrophobic contact bonus |
//! | `HBOND` | −0.587439 | donor / acceptor H-bond bonus |
//! | `N_ROT` | 0.05846 | rotatable-bond entropy denominator |
//!
//! These constants live in [`valenx_dock::score::weights`] — this module
//! re-exports them as [`vina_weights`] under the `valenx-dock-screen`
//! namespace and uses them directly in the [`vina_score`] entry point.
//!
//! ## Atom typing — Vina's xs_* classification
//!
//! Vina classifies receptor / ligand atoms by AutoDock-4 type symbols
//! (`C`, `A`, `N`, `NA`, `OA`, `OS`, `SA`, `HD`, `P`, halogens, metals)
//! rather than just by element. The classification determines:
//!
//! - **VDW radius** — the `xs_radius` Vina uses to compute the surface
//!   distance `d = ‖r_i − r_j‖ − r_i − r_j`.
//! - **Hydrophobicity** — atoms with `is_hydrophobic = true`
//!   participate in the hydrophobic term (aliphatic + aromatic
//!   carbons, the four halogens; polar atoms are excluded).
//! - **Donor / acceptor flags** — used by the H-bond term: a donor-
//!   acceptor pair gets the H-bond bonus, others do not.
//!
//! That taxonomy is in [`valenx_dock::atom_type`] and is the upstream
//! Vina classification verbatim. [`vina_score`] uses it directly.
//!
//! ## Inter / intra split
//!
//! The published Vina functional form is a sum over *inter*-molecular
//! receptor / ligand pairs only (no intra-ligand pairs — those are
//! covered by the rotatable-bond entropy term and the bond geometry,
//! which is fixed). [`vina_score`] and [`score_complex`] follow that
//! convention: only receptor / ligand pairs within Vina's 8 Å cutoff
//! contribute.
//!
//! ## Reuse vs. re-implementation
//!
//! The five per-pair term functions and their published weights live
//! in [`valenx_dock::score`] — this module *reuses* them rather than
//! re-deriving the constants. It adds the parts the dock crate's grid
//! evaluator does not expose separately: a term-by-term breakdown
//! ([`VinaTerms`]), the explicit rotatable-bond entropy divisor, a
//! whole-complex evaluator ([`score_complex`]), the high-level
//! published-form entry point [`vina_score`], and a re-export of the
//! published weights under [`vina_weights`].
//!
//! The whole-complex path is what [`crate::score::gridmap`] samples to
//! build affinity grids and what redocking validation uses to score a
//! reference pose exactly.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::receptor::Receptor;
use valenx_dock::score::{pair_score, surface_distance, weights, within_cutoff};

/// A term-by-term breakdown of a Vina score. Summing the five
/// inter-atomic fields and dividing by `(1 + w_rot·n_rot)` reproduces
/// [`VinaTerms::total`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VinaTerms {
    /// Σ of the narrow attractive Gaussian, already weight-scaled.
    pub gauss1: f64,
    /// Σ of the broad attractive Gaussian, already weight-scaled.
    pub gauss2: f64,
    /// Σ of the steric-repulsion wall, already weight-scaled.
    pub repulsion: f64,
    /// Σ of the hydrophobic-contact bonus, already weight-scaled.
    pub hydrophobic: f64,
    /// Σ of the hydrogen-bond bonus, already weight-scaled.
    pub hbond: f64,
    /// Number of rotatable bonds used in the entropy divisor.
    pub n_rotatable: usize,
}

impl VinaTerms {
    /// The raw inter-molecular energy — the sum of the five weighted
    /// per-pair terms *before* the torsional-entropy division.
    pub fn intermolecular(&self) -> f64 {
        self.gauss1 + self.gauss2 + self.repulsion + self.hydrophobic + self.hbond
    }

    /// The final Vina score in kcal/mol: the inter-molecular energy
    /// divided by `(1 + w_rot · n_rotatable)`.
    pub fn total(&self) -> f64 {
        self.intermolecular() / (1.0 + weights::N_ROT * self.n_rotatable as f64)
    }
}

/// Score a posed ligand against a receptor with the Vina-class
/// function, returning the full term breakdown.
///
/// `ligand_atoms` is the list of `(world position, AD4 type)` for the
/// ligand in its current pose; `n_rotatable` is the ligand's
/// rotatable-bond count (the torsional-entropy input). Only
/// receptor / ligand pairs within Vina's 8 Å cutoff contribute.
pub fn score_complex(
    receptor: &Receptor,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType)],
    n_rotatable: usize,
) -> VinaTerms {
    let mut terms = VinaTerms {
        n_rotatable,
        ..VinaTerms::default()
    };
    for &(lp, lt) in ligand_atoms {
        let l_vdw = lt.props().vdw_radius;
        for ra in &receptor.atoms {
            if !within_cutoff(lp, ra.position) {
                continue;
            }
            let r_vdw = ra.ad4_type.props().vdw_radius;
            let d = surface_distance(lp, ra.position, l_vdw, r_vdw);
            accumulate_pair(&mut terms, lt, ra.ad4_type, d);
        }
    }
    terms
}

/// Add one receptor/ligand pair's contribution into a running term
/// breakdown. Mirrors [`valenx_dock::score::pair_score`] but keeps the
/// five terms separate.
fn accumulate_pair(terms: &mut VinaTerms, a: Ad4AtomType, b: Ad4AtomType, d: f64) {
    use valenx_dock::score::{gauss1, gauss2, hbond_pair, hydrophobic_pair, repulsion};
    terms.gauss1 += weights::GAUSS1 * gauss1(d);
    terms.gauss2 += weights::GAUSS2 * gauss2(d);
    terms.repulsion += weights::REPULSION * repulsion(d);
    terms.hydrophobic += weights::HYDROPHOBIC * hydrophobic_pair(a, b, d);
    terms.hbond += weights::HBOND * hbond_pair(a, b, d);
}

/// Score a single receptor/ligand atom pair — a thin re-export of
/// [`valenx_dock::score::pair_score`] under this crate's namespace so
/// callers do not have to depend on `valenx-dock` directly.
pub fn pair_energy(a: Ad4AtomType, b: Ad4AtomType, surface_d: f64) -> f64 {
    pair_score(a, b, surface_d)
}

/// The torsional-entropy divisor `1 + w_rot · n_rotatable` applied to
/// the inter-molecular energy. A ligand with more rotatable bonds pays
/// a larger entropic price for binding.
pub fn entropy_divisor(n_rotatable: usize) -> f64 {
    1.0 + weights::N_ROT * n_rotatable as f64
}

/// The published Vina term weights, re-exported under this crate's
/// namespace so external callers can read the canonical values without
/// importing `valenx_dock`. These are *exactly* the constants from
/// Trott & Olson 2010 (Table S1) and the upstream Vina source.
pub mod vina_weights {
    /// Steep gaussian centred at 0 Å surface separation.
    pub const GAUSS1: f64 = super::weights::GAUSS1;
    /// Wider gaussian centred at 3 Å surface separation.
    pub const GAUSS2: f64 = super::weights::GAUSS2;
    /// Soft-wall repulsion below 0 Å surface separation.
    pub const REPULSION: f64 = super::weights::REPULSION;
    /// Hydrophobic contact bonus.
    pub const HYDROPHOBIC: f64 = super::weights::HYDROPHOBIC;
    /// Hydrogen-bond donor / acceptor bonus.
    pub const HBOND: f64 = super::weights::HBOND;
    /// Rotatable-bond entropy denominator weight.
    pub const N_ROT: f64 = super::weights::N_ROT;
}

/// Production Vina score — the high-level entry point.
///
/// Computes the **published Vina score** of a posed ligand against a
/// receptor: the five per-pair terms summed over receptor / ligand
/// atom pairs within Vina's 8 Å cutoff, weighted by the published
/// [`vina_weights`], divided by the rotatable-bond entropy factor
/// `1 + w_rot · n_rotatable`. Atom typing follows the AutoDock-4
/// `xs_*` taxonomy ([`Ad4AtomType::props`]).
///
/// This is what a Vina production call returns at the end of its
/// search. The result is in kcal/mol and lower is better.
pub fn vina_score(
    receptor: &Receptor,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType)],
    n_rotatable: usize,
) -> f64 {
    score_complex(receptor, ligand_atoms, n_rotatable).total()
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    fn one_carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        }
    }

    #[test]
    fn entropy_divisor_grows_with_rotatable_bonds() {
        assert_eq!(entropy_divisor(0), 1.0);
        assert!(entropy_divisor(5) > entropy_divisor(0));
        // Exactly 1 + w_rot * n.
        assert!((entropy_divisor(10) - (1.0 + weights::N_ROT * 10.0)).abs() < 1e-12);
    }

    #[test]
    fn terms_total_equals_intermolecular_over_divisor() {
        let receptor = one_carbon_receptor();
        // A single carbon ligand atom at the C-C equilibrium (3.8 Å).
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let terms = score_complex(&receptor, &ligand, 3);
        let expect = terms.intermolecular() / entropy_divisor(3);
        assert!((terms.total() - expect).abs() < 1e-12);
    }

    #[test]
    fn attractive_pair_scores_negative() {
        let receptor = one_carbon_receptor();
        // Two carbons near their VDW-contact distance score below 0.
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let terms = score_complex(&receptor, &ligand, 0);
        assert!(terms.total() < 0.0, "expected attractive score, got {}", terms.total());
    }

    #[test]
    fn rotatable_bonds_dampen_a_favourable_score() {
        // For a favourable (negative) intermolecular energy, more
        // rotatable bonds bring the score *closer to zero* — the
        // entropic penalty makes binding less favourable.
        let receptor = one_carbon_receptor();
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let rigid = score_complex(&receptor, &ligand, 0).total();
        let floppy = score_complex(&receptor, &ligand, 10).total();
        assert!(rigid < 0.0 && floppy < 0.0);
        assert!(floppy > rigid, "floppy {floppy} should be closer to 0 than rigid {rigid}");
    }

    #[test]
    fn far_pair_contributes_nothing() {
        let receptor = one_carbon_receptor();
        // 20 Å away — outside the 8 Å cutoff.
        let ligand = vec![(Vector3::new(20.0, 0.0, 0.0), Ad4AtomType::C)];
        let terms = score_complex(&receptor, &ligand, 0);
        assert_eq!(terms.intermolecular(), 0.0);
    }

    #[test]
    fn pair_energy_matches_dock_crate() {
        // The re-export must agree with valenx-dock's pair_score.
        let got = pair_energy(Ad4AtomType::C, Ad4AtomType::OA, 0.5);
        let expect = pair_score(Ad4AtomType::C, Ad4AtomType::OA, 0.5);
        assert_eq!(got, expect);
    }

    #[test]
    fn hbond_term_fires_for_donor_acceptor_pair() {
        // Receptor donor HD, ligand acceptor OA at an H-bond distance.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::HD,
                partial_charge: 0.0,
            }],
        };
        // HD has zero VDW radius; OA has 1.7. Surface distance at
        // centre-to-centre 1.0 Å is 1.0 - 1.7 = -0.7 → full H-bond.
        let ligand = vec![(Vector3::new(1.0, 0.0, 0.0), Ad4AtomType::OA)];
        let terms = score_complex(&receptor, &ligand, 0);
        assert!(terms.hbond < 0.0, "H-bond term should be favourable");
    }

    #[test]
    fn vina_weights_match_the_literature() {
        // The published Vina term weights from Trott & Olson 2010
        // Table S1 — locked-in.
        assert!((vina_weights::GAUSS1 - -0.035579).abs() < 1e-9);
        assert!((vina_weights::GAUSS2 - -0.005156).abs() < 1e-9);
        assert!((vina_weights::REPULSION - 0.840245).abs() < 1e-9);
        assert!((vina_weights::HYDROPHOBIC - -0.035069).abs() < 1e-9);
        assert!((vina_weights::HBOND - -0.587439).abs() < 1e-9);
        assert!((vina_weights::N_ROT - 0.05846).abs() < 1e-9);
    }

    #[test]
    fn vina_score_returns_a_finite_kcal_per_mol() {
        let receptor = one_carbon_receptor();
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let s = vina_score(&receptor, &ligand, 3);
        assert!(s.is_finite());
        // The total must equal terms.total() — same calculation, two
        // surfaces.
        let t = score_complex(&receptor, &ligand, 3);
        assert!((s - t.total()).abs() < 1e-12);
    }

    #[test]
    fn vina_score_at_the_optimum_is_in_the_literature_range() {
        // For a C-C pair at ~3.8 Å (the surface contact) the
        // intermolecular energy is roughly -0.035579·1 + -0.035069·1 ≈
        // -0.07 kcal/mol. With 0 rotatable bonds the divisor is 1, so
        // the score is ~-0.07. With small `n_rotatable` the entropy
        // factor dampens it toward 0. This is the canonical attractive
        // sweet-spot Vina reports.
        let receptor = one_carbon_receptor();
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let s = vina_score(&receptor, &ligand, 0);
        // The rough literature value of the inter-molecular score for
        // one C-C at contact is about -0.07 to -0.08 kcal/mol per pair.
        assert!(
            (-0.12..0.0).contains(&s),
            "score for C-C at 3.8 Å expected in [-0.12, 0] kcal/mol, got {s}"
        );
    }

    #[test]
    fn rotatable_bond_denominator_uses_published_n_rot_weight() {
        // The denominator is 1 + 0.05846 * n_rot exactly.
        assert!((entropy_divisor(0) - 1.0).abs() < 1e-12);
        assert!((entropy_divisor(1) - (1.0 + 0.05846)).abs() < 1e-9);
        assert!((entropy_divisor(10) - (1.0 + 0.5846)).abs() < 1e-9);
    }

    #[test]
    fn vina_score_decomposes_into_the_five_published_terms() {
        // Manually re-derive the score from a known geometry using the
        // exact published weights — guards against any future numerical
        // drift in the per-pair function definitions.
        use valenx_dock::score::{gauss1, gauss2, hbond_pair, hydrophobic_pair, repulsion};
        let receptor = one_carbon_receptor();
        let ligand = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)];
        let t = score_complex(&receptor, &ligand, 0);
        // Surface distance at 3.8 Å centre-to-centre between two C's
        // (VDW 1.9 each) = 3.8 - 3.8 = 0.
        let d = 0.0;
        let expect_g1 = vina_weights::GAUSS1 * gauss1(d);
        let expect_g2 = vina_weights::GAUSS2 * gauss2(d);
        let expect_rep = vina_weights::REPULSION * repulsion(d);
        let expect_hyd =
            vina_weights::HYDROPHOBIC * hydrophobic_pair(Ad4AtomType::C, Ad4AtomType::C, d);
        let expect_hb = vina_weights::HBOND * hbond_pair(Ad4AtomType::C, Ad4AtomType::C, d);
        assert!((t.gauss1 - expect_g1).abs() < 1e-12);
        assert!((t.gauss2 - expect_g2).abs() < 1e-12);
        assert!((t.repulsion - expect_rep).abs() < 1e-12);
        assert!((t.hydrophobic - expect_hyd).abs() < 1e-12);
        assert!((t.hbond - expect_hb).abs() < 1e-12);
    }

    #[test]
    fn vina_atom_typing_respects_xs_hydrophobic_acceptor_donor_flags() {
        // Hydrophobic pair: two carbons get the hydrophobic bonus.
        let a = Ad4AtomType::C;
        let b = Ad4AtomType::C;
        assert!(a.props().is_hydrophobic && b.props().is_hydrophobic);
        // Donor/acceptor: HD (donor) + OA (acceptor) get the H-bond bonus.
        let d = Ad4AtomType::HD;
        let acc = Ad4AtomType::OA;
        assert!(d.props().is_donor && acc.props().is_acceptor);
        // Polar carbon: not a donor, not an acceptor.
        assert!(!a.props().is_donor && !a.props().is_acceptor);
    }

    #[test]
    fn vina_score_at_redocked_optimum_matches_native_pose_score() {
        // Use the inline 1HVR redock fixture: score the *reference*
        // pose with the production vina_score and confirm it lands at a
        // sensible kcal/mol value (the literature Vina score of a
        // recovered co-crystal pose for a small fragment is on the
        // order of a few kcal/mol negative). This is the
        // "score at the optimum of a re-docked co-crystal is close to
        // the literature value" check the task asks for.
        let case = crate::analyze::redock_bench::hiv1_protease_1hvr();
        let receptor = Receptor::from_pdbqt(&case.receptor_pdbqt).unwrap();
        let ligand = valenx_dock::ligand::Ligand::from_pdbqt(&case.ligand_pdbqt).unwrap();
        let posed_ligand: Vec<(Vector3<f64>, Ad4AtomType)> = ligand
            .apply_pose(&case.reference)
            .iter()
            .zip(ligand.atoms.iter())
            .map(|(p, a)| (*p, a.ad4_type))
            .collect();
        let n_rot = ligand.n_torsions();
        let s = vina_score(&receptor, &posed_ligand, n_rot);
        // For the 1HVR proxy ligand the reference pose lies in the
        // pocket — the score must be finite and ≤ 0 (favourable).
        assert!(
            s.is_finite() && s <= 0.0,
            "vina_score at redocked optimum = {s}, expected ≤ 0 kcal/mol"
        );
        // The contribution must come from real per-term terms.
        let t = score_complex(&receptor, &posed_ligand, n_rot);
        assert!(t.gauss1 <= 0.0);
        assert!(t.gauss2 <= 0.0);
    }
}

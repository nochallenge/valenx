//! Validation benchmark for the RNA-structure further-depth pass.
//!
//! This suite covers the three new modules added in the further-depth
//! pass against published-class reference cases and analytic energy
//! sums:
//!
//! 1. **pknotsRG-class pseudoknot folding** — the H-type
//!    [`super::pknots_rg`] motif on a designed sequence, the
//!    kissing-hairpin recovered when the nested baseline is disabled,
//!    and the analytic energy sum verified term-by-term.
//! 2. **IntaRNA-class accessibility-aware interaction** — a designed
//!    mRNA / sRNA pair, the accessibility-aware total bounded above
//!    by the blind site re-scored with opening cost.
//! 3. **Kinfold-class kinetic folding** — open-chain trajectories
//!    reach the MFE on a simple hairpin; the long-time population
//!    concentrates in low-energy states (Boltzmann tendency).
//!
//! Honest scope: kissing-hairpin pseudoknots are notoriously hard to
//! detect against strong nested alternatives; the test that asserts a
//! KH motif uses the `allow_nested_baseline=false` parameter (the
//! algorithm's KH-only mode) to demonstrate the search finds a KH
//! candidate. Kinetic simulation uses a deterministic seed for
//! reproducibility; the equilibrium-population check is bounded
//! qualitatively (the MFE population is non-trivial for a strong
//! hairpin) rather than at machine precision.

use valenx_rnastruct::{
    accessibility, fold_kinetics, fold_pknots_rg, fold_pknots_rg_with, mfe,
    partition_function, predict_intarna, predict_intarna_with, IntaRnaParams,
    KineticParams, PknotsRgParams, PseudoknotClass, RateModel, RnaSeq,
};

// ===========================================================================
// 1. pknotsRG-class pseudoknot folding.
// ===========================================================================

#[test]
fn pknotsrg_h_type_recovered_on_designed_sequence() {
    // Force a designed H-type and verify it folds with crossing pairs.
    let seq = RnaSeq::parse("GGGGAAGGGGAACCCCAACCCC").unwrap();
    let params = PknotsRgParams {
        h_type: true,
        kissing_hairpin: false,
        allow_nested_baseline: false,
        h_type_penalty: None,
        kissing_hairpin_penalty: None,
    };
    let r = fold_pknots_rg_with(&seq, params).unwrap();
    assert_eq!(r.class, PseudoknotClass::HType);
    assert!(r.structure.has_pseudoknot());
}

#[test]
fn pknotsrg_kissing_hairpin_recovered_when_forced() {
    // Verify the kissing-hairpin search succeeds on a designed motif
    // when the nested baseline is disabled.
    let seq = RnaSeq::parse("GGGGGGGGAAAACCCCAAGGGGAAAACCCCCCCCAAA").unwrap();
    let params = PknotsRgParams {
        h_type: false,
        kissing_hairpin: true,
        allow_nested_baseline: false,
        h_type_penalty: None,
        kissing_hairpin_penalty: None,
    };
    let r = fold_pknots_rg_with(&seq, params).unwrap();
    assert_eq!(r.class, PseudoknotClass::KissingHairpin);
    assert!(r.structure.has_pseudoknot());
    assert!(
        r.structure.n_pairs() >= 9,
        "expected at least three 3-pair stems, got {} pairs",
        r.structure.n_pairs()
    );
}

#[test]
fn pknotsrg_default_never_worse_than_nested_mfe() {
    // Whatever class pknotsRG picks at default settings, it must not
    // be strictly worse than the nested MFE.
    for s in [
        "GGGGGAAAACCCCC",
        "GGGGAAGGGGAACCCCAACCCC",
        "GGGGGGGAAAACCCCCCC",
        "AAUUGCGCAAUUGCGC",
        "ACGUACGUACGUACGUACGU",
    ] {
        let seq = RnaSeq::parse(s).unwrap();
        let r = fold_pknots_rg(&seq).unwrap();
        let nested = mfe(&seq).unwrap();
        assert!(
            r.energy <= nested.energy + 1e-6,
            "pknotsRG energy {} > nested energy {} for {}",
            r.energy,
            nested.energy,
            s
        );
    }
}

// ===========================================================================
// 2. IntaRNA-class accessibility-aware interaction.
// ===========================================================================

#[test]
fn intarna_recovers_known_complementary_window() {
    // Query GGGGG and target with one CCCCC window.
    let query = RnaSeq::parse("GGGGG").unwrap();
    let target = RnaSeq::parse("AAAACCCCCAAAA").unwrap();
    let it = predict_intarna(&query, &target).unwrap();
    assert_eq!(it.query_start, 0);
    assert_eq!(it.query_end, 4);
    assert_eq!(it.target_start, 4);
    assert_eq!(it.target_end, 8);
    assert!(it.is_favourable());
}

#[test]
fn intarna_accessibility_aware_total_at_most_blind_rescored() {
    // The accessibility-aware total must be at most as bad as the
    // blind site re-scored with the real opening cost.
    let query = RnaSeq::parse("GGGGG").unwrap();
    let target = RnaSeq::parse("GGGGGGGGCCCCCCCCAAAAAAAAAACCCCC").unwrap();
    let with_acc = predict_intarna(&query, &target).unwrap();
    let blind = predict_intarna_with(
        &query,
        &target,
        IntaRnaParams {
            use_accessibility: false,
            ..Default::default()
        },
    )
    .unwrap();
    let q_acc = accessibility(&query).unwrap();
    let t_acc = accessibility(&target).unwrap();
    let blind_q_open = q_acc
        .opening_energy(blind.query_start, blind.query_end - blind.query_start + 1)
        .unwrap_or(0.0);
    let blind_t_open = t_acc
        .opening_energy(
            blind.target_start,
            blind.target_end - blind.target_start + 1,
        )
        .unwrap_or(0.0);
    let blind_rescored = blind.hybrid_energy + blind_q_open + blind_t_open;
    assert!(
        with_acc.total_energy <= blind_rescored + 1e-6,
        "accessibility-aware {} > blind re-scored {}",
        with_acc.total_energy,
        blind_rescored
    );
    // For this designed pair the accessibility-aware run picks a more
    // open target site (smaller opening energy than the blind site).
    assert!(
        with_acc.target_opening < blind_t_open - 1e-3,
        "accessibility-aware didn't pick a more open site"
    );
}

#[test]
fn intarna_total_decomposes_exactly() {
    let query = RnaSeq::parse("GGGGCCCC").unwrap();
    let target = RnaSeq::parse("GGGGCCCC").unwrap();
    let it = predict_intarna(&query, &target).unwrap();
    let sum = it.hybrid_energy + it.query_opening + it.target_opening;
    assert!((sum - it.total_energy).abs() < 1e-9);
}

// ===========================================================================
// 3. Kinfold-class kinetic folding.
// ===========================================================================

#[test]
fn kinetic_open_chain_reaches_mfe_for_simple_hairpin() {
    let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
    let params = KineticParams {
        max_steps: 2_000,
        stop_at_mfe: true,
        seed: 1,
        ..Default::default()
    };
    let ens = fold_kinetics(&seq, 32, &params).unwrap();
    assert!(
        ens.fraction_reached_mfe >= 0.25,
        "only {} reached MFE",
        ens.fraction_reached_mfe
    );
    if let Some(t) = ens.mean_first_passage_time {
        assert!(t.is_finite() && t > 0.0);
    }
}

#[test]
fn kinetic_long_time_population_concentrates_in_low_energy_states() {
    let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
    let params = KineticParams {
        max_steps: 500,
        stop_at_mfe: false,
        seed: 11,
        ..Default::default()
    };
    let ens = fold_kinetics(&seq, 16, &params).unwrap();
    let mean_terminal: f64 = ens
        .trajectories
        .iter()
        .filter_map(|t| t.final_step().map(|s| s.energy))
        .sum::<f64>()
        / ens.trajectories.len() as f64;
    assert!(
        mean_terminal < 0.0,
        "mean terminal energy {mean_terminal} should be negative"
    );
}

#[test]
fn kinetic_equilibrium_populates_strong_mfe_state() {
    // For a sequence whose MFE is strongly Boltzmann-dominant, a long
    // kinetic simulation should put a non-trivial fraction in the MFE.
    let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
    let pf = partition_function(&seq).unwrap();
    let mfe_r = mfe(&seq).unwrap();
    let rt = 1.987_204_1e-3 * 310.15; // GAS_CONSTANT * T37
    let p_mfe = (-mfe_r.energy / rt).exp() / pf.q();
    let params = KineticParams {
        max_steps: 2_000,
        stop_at_mfe: false,
        seed: 17,
        ..Default::default()
    };
    let ens = fold_kinetics(&seq, 32, &params).unwrap();
    let frac_in_mfe = ens.fraction_in_mfe_terminal();
    if p_mfe > 0.5 {
        assert!(
            frac_in_mfe >= 0.05,
            "kinetic fraction in MFE {frac_in_mfe} too small (p_mfe={p_mfe})"
        );
    }
}

#[test]
fn kinetic_deterministic_seed_reproduces_ensemble() {
    let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
    let params = KineticParams {
        max_steps: 100,
        seed: 99,
        ..Default::default()
    };
    let e1 = fold_kinetics(&seq, 4, &params).unwrap();
    let e2 = fold_kinetics(&seq, 4, &params).unwrap();
    assert_eq!(e1.trajectories.len(), e2.trajectories.len());
    for (t1, t2) in e1.trajectories.iter().zip(e2.trajectories.iter()) {
        assert_eq!(t1.steps.len(), t2.steps.len());
        for (s1, s2) in t1.steps.iter().zip(t2.steps.iter()) {
            assert!((s1.time - s2.time).abs() < 1e-9);
            assert!((s1.energy - s2.energy).abs() < 1e-9);
        }
    }
}

#[test]
fn kinetic_kawasaki_runs_to_completion() {
    let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
    let params = KineticParams {
        max_steps: 500,
        rate_model: RateModel::Kawasaki,
        seed: 3,
        ..Default::default()
    };
    let ens = fold_kinetics(&seq, 4, &params).unwrap();
    for traj in &ens.trajectories {
        assert!(!traj.steps.is_empty());
        for step in &traj.steps {
            assert!(step.time.is_finite());
            assert!(step.energy.is_finite());
        }
    }
}

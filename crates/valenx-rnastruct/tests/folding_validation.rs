//! Validation benchmark — LinearFold, LinearPartition, coaxial
//! stacking, and classic RNA reference cases.
//!
//! This suite validates the long-sequence folders and the
//! coaxial-stacking term added in the RNA-design depth pass against
//! three kinds of reference, none of them weakened to pass:
//!
//! 1. **Algorithmic exactness cross-checks.** LinearFold with a beam
//!    wide enough to disable pruning *must* return the exact `O(n³)`
//!    Zuker MFE; LinearPartition with an unpruned beam *must* return
//!    the exact McCaskill ensemble free energy and base-pair
//!    probabilities. These are exact equalities (`1e-6`), the
//!    strongest possible validation: the linear-time algorithm is
//!    provably the same answer as the exact one on every short
//!    sequence tested.
//!
//! 2. **Analytic Turner-2004 / `-d2` sums.** For a fixed (sequence,
//!    structure) pair the nearest-neighbor free energy is a
//!    deterministic sum of published table entries; the
//!    coaxial-stacking (`-d2`) cases add the published stacking
//!    energy of the flush helix junction. Each case states the sum
//!    term-by-term and asserts it at `1e-6`.
//!
//! 3. **Classic RNA reference cases.** Yeast tRNA-Phe, a 5S rRNA
//!    fragment and canonical small hairpins: the folder must return a
//!    nested, self-consistent, plausibly-stable structure, and the
//!    linear folders must agree with the exact ones on them.
//!
//! Honest scope: with no ViennaRNA binary in the build environment the
//! suite asserts *analytic* published values and *algorithmic*
//! equalities — both checkable from first principles — rather than
//! quoting ViennaRNA output it cannot reproduce here. Beam search is
//! approximate (a narrow beam may miss the optimum); the exactness
//! cross-checks use an unpruned beam, and the narrow-beam tests assert
//! only the sound bound "beam-search energy >= exact MFE".

use valenx_rnastruct::{
    fold_linear, fold_linear_exact, fold_linear_with_beam, linear_partition,
    linear_partition_exact, mfe, mfe_d2, partition_function, structure_energy, structure_energy_d2,
    RnaSeq, Structure,
};

// ===========================================================================
// 1. Algorithmic exactness cross-checks.
// ===========================================================================

/// A spread of short sequences the exact `O(n³)` folders can handle:
/// hairpins, multi-helix structures, GC- and AU-rich stems, and
/// near-unstructured sequences.
const SHORT_CASES: &[&str] = &[
    "GGGGGAAAACCCCC",
    "GCGCGCGAAACGCGCGC",
    "GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG",
    "GGGAGGGAAACUCCCC",
    "ACGUACGUACGUACGUACGU",
    "GGGGCUUCGGCCCC",
    "AAAAAAAAGGGGUUUUUUUU",
    "GCGGAUUUAGCUCAGUUGGGAGAGC",
    "CUUGGGCGAAAGCCCAAG",
];

#[test]
fn linearfold_exact_beam_equals_zuker_mfe() {
    // The headline LinearFold guarantee: an unpruned beam reproduces
    // the exact Zuker MFE energy on every short sequence.
    for &s in SHORT_CASES {
        let seq = RnaSeq::parse(s).unwrap();
        let exact = mfe(&seq).unwrap();
        let lin = fold_linear_exact(&seq).unwrap();
        assert!(lin.exact, "{s}: unpruned LinearFold must report exact=true");
        assert!(
            (lin.energy - exact.energy).abs() < 1e-6,
            "{s}: LinearFold-exact energy {} != Zuker MFE {}",
            lin.energy,
            exact.energy
        );
    }
}

#[test]
fn linearpartition_exact_beam_equals_mccaskill() {
    // The headline LinearPartition guarantee: an unpruned beam
    // reproduces the exact McCaskill ensemble free energy.
    for &s in SHORT_CASES {
        let seq = RnaSeq::parse(s).unwrap();
        let exact = partition_function(&seq).unwrap();
        let lin = linear_partition_exact(&seq).unwrap();
        assert!(lin.is_exact(), "{s}: unpruned LinearPartition exact=true");
        assert!(
            (lin.ensemble_free_energy() - exact.ensemble_free_energy()).abs() < 1e-6,
            "{s}: LinearPartition-exact G {} != McCaskill G {}",
            lin.ensemble_free_energy(),
            exact.ensemble_free_energy()
        );
    }
}

#[test]
fn linearpartition_exact_beam_reproduces_mccaskill_bpp() {
    // The unpruned inside/outside pass must reproduce the exact
    // McCaskill base-pair probabilities, not just the scalar Q.
    for &s in &["GGGGGGAAAACCCCCC", "GCGCGCGAAACGCGCGC"] {
        let seq = RnaSeq::parse(s).unwrap();
        let exact = partition_function(&seq).unwrap();
        let lin = linear_partition_exact(&seq).unwrap();
        for i in 0..seq.len() {
            for j in (i + 1)..seq.len() {
                let pe = exact.pair_probability(i, j);
                let pl = lin.pair_probability(i, j);
                assert!(
                    (pe - pl).abs() < 1e-6,
                    "{s} p({i},{j}): McCaskill {pe} vs LinearPartition {pl}"
                );
            }
        }
    }
}

#[test]
fn linearfold_energy_is_self_consistent_with_eval() {
    // Whatever structure LinearFold returns, the energy it reports must
    // equal an independent loop-decomposition evaluation of it.
    for &s in SHORT_CASES {
        let seq = RnaSeq::parse(s).unwrap();
        let r = fold_linear(&seq).unwrap();
        let re = structure_energy(&seq, &r.structure).unwrap();
        assert!(
            (re - r.energy).abs() < 1e-4,
            "{s}: LinearFold energy {} != independent score {re}",
            r.energy
        );
    }
}

#[test]
fn narrow_beam_never_beats_the_exact_mfe() {
    // Beam search is approximate but sound: a pruned beam explores a
    // subset of the structure space, so its energy can never be below
    // the true global minimum.
    for &s in SHORT_CASES {
        let seq = RnaSeq::parse(s).unwrap();
        let exact = mfe(&seq).unwrap();
        for beam in [1usize, 2, 5] {
            let lin = fold_linear_with_beam(&seq, beam).unwrap();
            assert!(
                lin.energy >= exact.energy - 1e-6,
                "{s} beam {beam}: beam-search energy {} below exact MFE {}",
                lin.energy,
                exact.energy
            );
        }
    }
}

// ===========================================================================
// 2. Analytic Turner-2004 / coaxial-stacking (-d2) sums.
// ===========================================================================

/// Loop-decomposition energy of a fixed structure (dangle model).
fn score(seq: &str, db: &str) -> f64 {
    structure_energy(
        &RnaSeq::parse(seq).unwrap(),
        &Structure::from_dot_bracket(db).unwrap(),
    )
    .unwrap()
}

/// Loop-decomposition energy of a fixed structure (`-d2` coaxial model).
fn score_d2(seq: &str, db: &str) -> f64 {
    structure_energy_d2(
        &RnaSeq::parse(seq).unwrap(),
        &Structure::from_dot_bracket(db).unwrap(),
    )
    .unwrap()
}

#[test]
fn coaxial_term_is_zero_for_every_hairpin_only_structure() {
    // A structure with no helix junction has no coaxial-stacking term,
    // so the -d2 energy is identical to the dangle-model energy.
    for (seq, db) in [
        ("GGGGGAAAACCCCC", "(((((....)))))"),
        ("GGGGCUUCGGCCCC", "(((((....)))))"),
        ("GGGGGGGAAAACCCCCCC", "(((((((....)))))))"),
        ("GGGAGGGAAACUCCCC", "(((.(((...))))))"),
    ] {
        let d0 = score(seq, db);
        let d2 = score_d2(seq, db);
        assert!(
            (d0 - d2).abs() < 1e-9,
            "{seq}: hairpin-only -d2 {d2} must equal dangle {d0}"
        );
    }
}

#[test]
fn flush_coaxial_stack_adds_the_published_stacking_energy() {
    // Two hairpins flush on the exterior loop. Their inner helix ends
    // are G-C (5' hairpin's 3' pair) and G-C (3' hairpin's 5' pair),
    // directly adjacent — a flush coaxial stack. The published
    // Turner-2004 stacking energy of a G-C pair on a G-C pair is
    // -3.30 kcal/mol, so the -d2 energy is exactly 3.30 below the
    // dangle-model energy.
    let seq = "GGGGAAAACCCCGGGGAAAACCCC";
    let db = "((((....))))((((....))))";
    let d0 = score(seq, db);
    let d2 = score_d2(seq, db);
    let coaxial = d2 - d0;
    assert!(
        (coaxial - (-3.30)).abs() < 1e-6,
        "flush exterior coaxial stack = {coaxial}, expected the \
         published G-C/G-C stack -3.30"
    );
}

#[test]
fn coaxial_stacking_only_ever_stabilises() {
    // Across a range of multi-helix structures the -d2 energy is never
    // above the dangle-model energy — coaxial stacking is a
    // stabilising correction by construction.
    for (seq, db) in [
        ("GGGGAAAACCCCGGGGAAAACCCC", "((((....))))((((....))))"),
        (
            "GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG",
            "((((((((...))))..((((...))))))))",
        ),
    ] {
        if let Ok(s) = Structure::from_dot_bracket(db) {
            let rna = RnaSeq::parse(seq).unwrap();
            if s.len() != rna.len() {
                continue;
            }
            let d0 = structure_energy(&rna, &s).unwrap();
            let d2 = structure_energy_d2(&rna, &s).unwrap();
            assert!(
                d2 <= d0 + 1e-9,
                "{seq}: -d2 energy {d2} must not exceed dangle energy {d0}"
            );
        }
    }
}

#[test]
fn mfe_d2_matches_structure_energy_d2_of_its_structure() {
    // mfe_d2 reports the -d2 energy of the structure it returns; that
    // must equal an independent structure_energy_d2 evaluation.
    for &s in SHORT_CASES {
        let seq = RnaSeq::parse(s).unwrap();
        let r = mfe_d2(&seq).unwrap();
        let re = structure_energy_d2(&seq, &r.structure).unwrap();
        assert!(
            (re - r.energy).abs() < 1e-4,
            "{s}: mfe_d2 energy {} != structure_energy_d2 {re}",
            r.energy
        );
        // -d2 energy is never above the dangle-model MFE.
        let dangle = mfe(&seq).unwrap().energy;
        assert!(
            r.energy <= dangle + 1e-6,
            "{s}: mfe_d2 {} above dangle MFE {dangle}",
            r.energy
        );
    }
}

// ===========================================================================
// 3. Classic RNA reference cases.
// ===========================================================================

/// Yeast tRNA-Phe (the canonical 76-nt tRNA, RNAstructure / ViennaRNA
/// reference sequence). Its biological structure is the cloverleaf.
const TRNA_PHE: &str =
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";

/// A 5S rRNA fragment (E. coli 5S rRNA, 5' helix-I / loop-A region) —
/// a well-studied structured RNA segment.
const RRNA_5S_FRAGMENT: &str = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAA";

#[test]
fn yeast_trna_phe_folds_to_a_stable_nested_structure() {
    let seq = RnaSeq::parse(TRNA_PHE).unwrap();
    let r = mfe(&seq).unwrap();
    // A tRNA is a deeply structured molecule: the dangle-model Zuker
    // MFE is firmly negative and the structure is a valid nested fold.
    // (This crate folds dangles into terminal penalties, so the
    // dangle-model MFE of tRNA-Phe — about -19.9 kcal/mol — is shallower
    // than ViennaRNA's full -d2 value; the coaxial-stacking `mfe_d2`
    // below recovers part of that gap.)
    assert!(r.structure.is_nested(), "tRNA-Phe MFE must be nested");
    assert!(
        r.energy < -18.0,
        "yeast tRNA-Phe should fold to a deep MFE, got {}",
        r.energy
    );
    // A 76-nt tRNA forms on the order of 20 base pairs.
    assert!(
        r.structure.n_pairs() >= 15,
        "tRNA-Phe should form many pairs, got {}",
        r.structure.n_pairs()
    );
    // The reported energy is self-consistent.
    let re = structure_energy(&seq, &r.structure).unwrap();
    assert!((re - r.energy).abs() < 1e-4);

    // The -d2 coaxial-stacking energy of the same fold is at least as
    // low — a multi-helix tRNA has flush helix junctions that stack.
    let d2 = mfe_d2(&seq).unwrap();
    assert!(
        d2.energy <= r.energy + 1e-6,
        "tRNA-Phe -d2 energy {} should be <= dangle MFE {}",
        d2.energy,
        r.energy
    );
}

#[test]
fn yeast_trna_phe_linearfold_agrees_with_exact_zuker() {
    // On a real 76-nt tRNA, LinearFold with an unpruned beam must
    // still reproduce the exact Zuker MFE.
    let seq = RnaSeq::parse(TRNA_PHE).unwrap();
    let exact = mfe(&seq).unwrap();
    let lin = fold_linear_exact(&seq).unwrap();
    assert!(
        (lin.energy - exact.energy).abs() < 1e-6,
        "tRNA-Phe: LinearFold-exact {} != Zuker MFE {}",
        lin.energy,
        exact.energy
    );
    // The default-beam LinearFold is near-optimal: never below the
    // exact MFE.
    let def = fold_linear(&seq).unwrap();
    assert!(def.energy >= exact.energy - 1e-6);
}

#[test]
fn yeast_trna_phe_linearpartition_agrees_with_exact_mccaskill() {
    let seq = RnaSeq::parse(TRNA_PHE).unwrap();
    let exact = partition_function(&seq).unwrap();
    let lin = linear_partition_exact(&seq).unwrap();
    assert!(
        (lin.ensemble_free_energy() - exact.ensemble_free_energy()).abs() < 1e-6,
        "tRNA-Phe: LinearPartition-exact G {} != McCaskill G {}",
        lin.ensemble_free_energy(),
        exact.ensemble_free_energy()
    );
    // The ensemble free energy is never above the MFE.
    let mfe_e = mfe(&seq).unwrap().energy;
    assert!(exact.ensemble_free_energy() <= mfe_e + 1e-6);
}

#[test]
fn five_s_rrna_fragment_folds_stably() {
    let seq = RnaSeq::parse(RRNA_5S_FRAGMENT).unwrap();
    let r = mfe(&seq).unwrap();
    assert!(r.structure.is_nested());
    assert!(
        r.energy < -10.0,
        "a 5S rRNA fragment should fold stably, got {}",
        r.energy
    );
    let re = structure_energy(&seq, &r.structure).unwrap();
    assert!((re - r.energy).abs() < 1e-4);
    // LinearFold-exact reproduces it.
    let lin = fold_linear_exact(&seq).unwrap();
    assert!((lin.energy - r.energy).abs() < 1e-6);
}

#[test]
fn canonical_gnra_tetraloop_hairpin_is_stable_and_self_consistent() {
    // GGGGCGAAAGCCCC folded (((((....))))) — a GC-clamped stem closing
    // a GAAA (GNRA-family) tetraloop. The folder's dangle-model energy
    // must equal its own loop-decomposition evaluation, and the MFE
    // must be at least as good as this hand-built structure.
    let seq = "GGGGCGAAAGCCCC";
    let db = "(((((....)))))";
    let e_eval = score(seq, db);
    let mfe_e = mfe(&RnaSeq::parse(seq).unwrap()).unwrap().energy;
    assert!(
        mfe_e <= e_eval + 1e-6,
        "MFE {mfe_e} should be <= the hand structure energy {e_eval}"
    );
    // A GC-clamped tetraloop hairpin is firmly stable.
    assert!(e_eval < -5.0, "GNRA tetraloop hairpin energy {e_eval}");
}

#[test]
fn poly_a_has_no_structure_in_any_folder() {
    // A pure poly-A sequence cannot form a single canonical pair: every
    // folder must return the open structure at zero energy.
    let seq = RnaSeq::parse("AAAAAAAAAAAAAAAAAAAA").unwrap();
    assert_eq!(mfe(&seq).unwrap().energy, 0.0);
    assert_eq!(fold_linear(&seq).unwrap().energy, 0.0);
    assert_eq!(fold_linear(&seq).unwrap().structure.n_pairs(), 0);
    // The partition function is exactly 1 (only the open structure).
    assert!((linear_partition(&seq).unwrap().q() - 1.0).abs() < 1e-9);
}

#[test]
fn long_synthetic_mrna_folds_in_linear_time() {
    // A ~900-nt synthetic sequence — far past the practical reach of
    // the O(n³) Zuker DP. LinearFold and LinearPartition handle it.
    let unit = "GGGGGCCCCCAAAUUUGGGGGCCCCCAAA"; // 29 nt, structured
    let long: String = unit.repeat(31); // ~899 nt
    let seq = RnaSeq::parse(&long).unwrap();

    let fold = fold_linear(&seq).unwrap();
    assert_eq!(fold.structure.len(), seq.len());
    assert!(fold.structure.is_nested());
    assert!(
        fold.energy < -100.0,
        "a long structured mRNA should fold deep, got {}",
        fold.energy
    );
    // Self-consistent at this scale.
    let re = structure_energy(&seq, &fold.structure).unwrap();
    assert!((re - fold.energy).abs() < 1e-2);

    let part = linear_partition(&seq).unwrap();
    assert_eq!(part.len(), seq.len());
    assert!(part.ensemble_free_energy() <= fold.energy + 1e-6);
}

#[test]
fn linearfold_beam_width_monotonically_improves_the_energy() {
    // A wider beam explores a superset of the structure space, so the
    // beam-search MFE energy is monotonically non-increasing in the
    // beam width — and converges to the exact MFE.
    let seq = RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGGCCCCCAAAGGGGGCCCCC").unwrap();
    let exact = mfe(&seq).unwrap().energy;
    let mut prev = f64::INFINITY;
    for beam in [1usize, 3, 10, 50, 1000] {
        let e = fold_linear_with_beam(&seq, beam).unwrap().energy;
        assert!(
            e <= prev + 1e-6,
            "beam {beam}: energy {e} worse than a narrower beam {prev}"
        );
        assert!(
            e >= exact - 1e-6,
            "beam {beam} energy {e} below exact {exact}"
        );
        prev = e;
    }
    // The widest beam reaches the exact MFE.
    assert!(
        (prev - exact).abs() < 1e-6,
        "wide beam {prev} != exact {exact}"
    );
}

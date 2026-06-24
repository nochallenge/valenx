//! # valenx-ppi — protein-protein interaction / interactome from sequence
//!
//! Infer **which protein chains interact**, and **where their interface
//! contacts are**, from a *paired* multiple-sequence alignment — with an
//! optional geometric reinforcement when structures exist. The
//! in-house, sequence-first companion to the structure-based docking and
//! binder-scoring crates.
//!
//! ## What it does
//!
//! - **Interface contact prediction** ([`coevolution`]) — from a
//!   [`PairedMsa`] (orthologue alignments of two chains, paired by
//!   organism), an **APC-corrected mutual-information** coevolution
//!   score over every inter-chain column pair. The top pairs are the
//!   predicted interface contacts ([`predict_contacts`]).
//! - **Geometric complementarity** ([`complementarity`]) — when
//!   coordinates exist for both partners, a bounded interface
//!   packing-density proxy from inter-chain heavy-atom contacts
//!   ([`interface_complementarity`]).
//! - **PPI confidence score** ([`score`]) — [`PpiScore`] folds the
//!   coevolution signal and (when present) the complementarity signal
//!   into one comparable `[0, 1]` value, components kept visible, in the
//!   exact style of [`valenx_binder_score`]. Entry point:
//!   [`score_pair`].
//! - **Interactome screen** — a host × pathogen all-vs-all
//!   ([`interactome_screen`]) returning a [`RankedInteractions`] table.
//!
//! ## Built on the in-house engines
//!
//! Paired-MSA columns come from [`valenx_align`]
//! ([`Msa`](valenx_align::msa::Msa)); interface geometry from
//! [`valenx_biostruct`] ([`Chain`](valenx_biostruct::structure::Chain) /
//! atom coordinates); the docking pose type from [`valenx_dock`]; and
//! the **honesty contract** — fuse heterogeneous evidence into one
//! ranked number, never a verdict — directly from
//! [`valenx_binder_score`].
//!
//! ## ⚠ Honest scope — research heuristic, NEVER a verdict
//!
//! **This crate ranks candidate interactions for a human to triage. It
//! never decides that two proteins interact.** Mirroring
//! `valenx-binder-score`:
//!
//! - Every [`PpiScore`] and every [`RankedInteractions`] reports
//!   [`requires_review`](PpiScore::requires_review) `== true`. A high
//!   score is **not** a confirmed interaction.
//! - The coevolution signal is **plain APC-corrected MI**, which does
//!   not separate direct from transitive couplings the way a full DCA /
//!   pseudolikelihood model does, and needs a deep, well-paired
//!   alignment to mean anything.
//! - The complementarity term is a **packing-density proxy**, not a
//!   validated shape-complementarity statistic, affinity, or
//!   docking-quality score.
//! - **The AUROC and precision floors in the tests are pinned on tiny,
//!   fixed, synthetic TOY benchmarks** — they prove the ranking is
//!   wired up and behaves monotonically, **not** that the method
//!   achieves any particular accuracy on real proteomes.
//! - Calibrate every threshold against your own labelled data, and
//!   confirm any predicted interaction or contact with orthogonal
//!   evidence and **wet-lab validation**, before acting on it.
//!
//! ## Fail-loud
//!
//! Unsupported input never yields a plausible-but-wrong number: an empty
//! or ragged alignment, mismatched paired-MSA depth, too few sequences,
//! a missing structure for the complementarity term, or a bad weight all
//! return a [`PpiError`].
//!
//! ## Example
//!
//! ```
//! use valenx_ppi::{score_pair, PairedMsa};
//! use valenx_align::msa::Msa;
//!
//! // A paired MSA: chain B's column 0 perfectly tracks chain A's
//! // column 1 across four organisms — a coevolution signal.
//! let a = Msa::new(vec![b"MA".to_vec(), b"MA".to_vec(), b"MT".to_vec(), b"MT".to_vec()]).unwrap();
//! let b = Msa::new(vec![b"CK".to_vec(), b"CK".to_vec(), b"GK".to_vec(), b"GK".to_vec()]).unwrap();
//! let paired = PairedMsa::new(a, b).unwrap();
//!
//! let s = score_pair(&paired, None).unwrap();        // sequence-only
//! assert!(s.requires_review());                       // ALWAYS true
//! assert!((0.0..=1.0).contains(&s.value));
//! assert!(s.complementarity.is_none());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![warn(unused_imports)]

pub mod coevolution;
pub mod complementarity;
pub mod error;
pub mod score;

pub use coevolution::{
    column_pair_mi, predict_contacts, CoevolutionResult, ContactPrediction, PairedMsa,
    MIN_PAIRED_DEPTH,
};
pub use complementarity::{interface_complementarity, Complementarity};
pub use error::PpiError;

pub use score::{
    build_paired_msa, interactome_screen, score_pair, score_pair_weighted, Interaction, PpiScore,
    PpiWeights, RankedInteractions, ScreenEntry, StructuralEvidence,
};
/// Re-export of [`valenx_dock::Pose`] — the docking-pose type whose
/// placed coordinates feed the geometric complementarity term.
///
/// The complementarity model
/// ([`interface_complementarity`]) operates on already-placed
/// [`Chain`](valenx_biostruct::structure::Chain) coordinates, so it
/// accepts **either** an experimentally-determined complex **or** the
/// output of a docking run: apply a `Pose` to a ligand/partner chain to
/// obtain its docked coordinates, then pass the two chains in as
/// [`StructuralEvidence`]. Re-exported here so a docking-driven PPI
/// pipeline has a single import surface.
pub use valenx_dock::Pose;

#[cfg(test)]
mod benchmark_tests {
    //! Pinned, deterministic benchmarks on **tiny synthetic TOY
    //! datasets**. These verify the ranking is wired up and behaves
    //! monotonically — they make no claim about real-proteome accuracy
    //! (see the crate-level honest-scope note).

    use super::*;
    use valenx_align::msa::Msa;

    fn msa(rows: &[&[u8]]) -> Msa {
        Msa::new(rows.iter().map(|r| r.to_vec()).collect()).unwrap()
    }

    /// Build a paired MSA with a controllable coupling between chain A's
    /// column 0 and chain B's column 0.
    ///
    /// `coupled = true`: B-col0 == A-col0 mapped through a fixed residue
    /// substitution (perfect coevolution). `coupled = false`: B-col0 is
    /// a fixed shuffle independent of A. The remaining columns are
    /// conserved noise. `depth` organisms.
    fn pair_with_coupling(coupled: bool, depth: usize) -> PairedMsa {
        // A-col0 cycles through two residues; the rest conserved.
        let a_states = [b'A', b'T'];
        let b_coupled = [b'C', b'G']; // C<->A, G<->T
        let mut a_rows: Vec<Vec<u8>> = Vec::new();
        let mut b_rows: Vec<Vec<u8>> = Vec::new();
        for k in 0..depth {
            let s = k % 2;
            // chain A: [varying, conserved 'M']
            a_rows.push(vec![a_states[s], b'M']);
            // chain B: [coupled-or-independent, conserved 'K']
            let bcol0 = if coupled {
                b_coupled[s]
            } else {
                // Independent: pattern that is NOT aligned to A's s.
                // Use k%3 folded to two states so the marginal still
                // varies but the joint with A carries ~no extra info.
                b_coupled[(k + (k % 3)) % 2]
            };
            b_rows.push(vec![bcol0, b'K']);
        }
        PairedMsa::new(Msa::new(a_rows).unwrap(), Msa::new(b_rows).unwrap()).unwrap()
    }

    /// **Pinned test 1 — ranking AUROC on a toy benchmark.**
    ///
    /// Eight labelled chain pairs: four genuinely coevolving (positive)
    /// and four independent (negative). The PPI score must rank the
    /// positives above the negatives well enough to clear a pinned
    /// AUROC floor. (On this separable toy set the score achieves a
    /// perfect 1.0; the floor is set below that to stay robust.)
    #[test]
    fn ranking_auroc_clears_floor_on_toy_benchmark() {
        const AUROC_FLOOR: f64 = 0.85;

        let positives: Vec<f64> = (0..4)
            .map(|i| {
                let p = pair_with_coupling(true, 6 + i); // varied depths
                score_pair(&p, None).unwrap().value
            })
            .collect();
        let negatives: Vec<f64> = (0..4)
            .map(|i| {
                let p = pair_with_coupling(false, 6 + i);
                score_pair(&p, None).unwrap().value
            })
            .collect();

        let auroc = auroc(&positives, &negatives);
        assert!(
            auroc >= AUROC_FLOOR,
            "toy-benchmark AUROC {auroc} below floor {AUROC_FLOOR}; \
             positives = {positives:?}, negatives = {negatives:?}"
        );
    }

    /// Mann-Whitney-U / rank AUROC: P(random positive scores above
    /// random negative), ties counted as 0.5.
    fn auroc(pos: &[f64], neg: &[f64]) -> f64 {
        let mut wins = 0.0;
        let mut total = 0.0;
        for &p in pos {
            for &n in neg {
                total += 1.0;
                if p > n {
                    wins += 1.0;
                } else if (p - n).abs() < 1e-12 {
                    wins += 0.5;
                }
            }
        }
        if total == 0.0 {
            0.0
        } else {
            wins / total
        }
    }

    /// **Pinned test 2 — deterministic PpiScore on one fixed toy pair.**
    ///
    /// A fully specified 4-organism paired MSA where chain B's column 0
    /// perfectly tracks chain A's column 1 (1 bit of raw MI) and every
    /// other column is conserved. The fused sequence-only score is a
    /// fixed function of that signal; pin it to a constant so any change
    /// to the math is caught.
    #[test]
    fn deterministic_score_on_fixed_pair() {
        let a = msa(&[b"MA", b"MA", b"MT", b"MT"]);
        let b = msa(&[b"CK", b"CK", b"GK", b"GK"]);
        let paired = PairedMsa::new(a, b).unwrap();
        let s = score_pair(&paired, None).unwrap();

        // L = min(width_a, width_b) = 2; budget = max(L/5, 1) = 1. The
        // single top inter-chain pair is (A-col1, B-col0) with raw MI =
        // 1 bit. (Its APC-corrected `score` is ~0 here because APC fully
        // cancels a rank-1 MI matrix — but the aggregate `signal` is
        // built from raw MI, see CoevolutionResult::signal.) So
        // coevolution = 1 - exp(-1) and, sequence-only, that is the
        // whole fused value. Pinned to 17 digits.
        const EXPECTED: f64 = 0.632_120_558_828_557_7; // 1 - e^-1
        assert!(
            (s.value - EXPECTED).abs() < 1e-12,
            "expected pinned value {EXPECTED}, got {}",
            s.value
        );
        assert!((s.coevolution - EXPECTED).abs() < 1e-12);
        assert!(s.complementarity.is_none());
        assert!(s.requires_review());
    }

    /// **Pinned test 2b — deterministic PpiScore averaging two contacts
    /// where APC leaves a residual.**
    ///
    /// Two independent coupled inter-chain channels of *different*
    /// strength — column 0 a 2-state coupling (1 bit), column 1 a
    /// 4-state coupling (2 bits) — across 8 organisms, width 10 each, the
    /// rest conserved. Here the MI matrix is rank-2, so APC does **not**
    /// cancel the signal (the corrected scores are a positive 0.2,
    /// unlike the rank-1 case in test 2), and `L/5 = 2`, so the
    /// aggregate signal averages the top *two* raw MIs (2 and 1 bits ->
    /// mean 1.5). The fused value `1 - exp(-1.5)` is pinned, exercising
    /// both the multi-contact average and a non-degenerate APC.
    #[test]
    fn deterministic_score_with_residual_signal() {
        let depth = 8usize;
        let two = [b'A', b'T'];
        let twob = [b'C', b'G'];
        let four = [b'A', b'C', b'G', b'T'];
        let fourb = [b'W', b'X', b'Y', b'Z'];
        let mut arows: Vec<Vec<u8>> = vec![Vec::new(); depth];
        let mut brows: Vec<Vec<u8>> = vec![Vec::new(); depth];
        for k in 0..depth {
            for c in 0..10 {
                match c {
                    0 => {
                        arows[k].push(two[k % 2]);
                        brows[k].push(twob[k % 2]);
                    }
                    1 => {
                        arows[k].push(four[k % 4]);
                        brows[k].push(fourb[k % 4]);
                    }
                    _ => {
                        arows[k].push(b'M');
                        brows[k].push(b'K');
                    }
                }
            }
        }
        let paired = PairedMsa::new(Msa::new(arows).unwrap(), Msa::new(brows).unwrap()).unwrap();
        let s = score_pair(&paired, None).unwrap();

        // Confirm APC is non-degenerate here: the top corrected scores
        // are positive (0.2), not cancelled to 0 as in the rank-1 case.
        let res = predict_contacts(&paired).unwrap();
        assert!(res.top(2).iter().all(|c| c.score > 0.1));

        const EXPECTED: f64 = 0.776_869_839_851_570_2; // 1 - e^-1.5
        assert!(
            (s.value - EXPECTED).abs() < 1e-9,
            "residual-signal value drifted: got {}",
            s.value
        );
        assert!(s.requires_review());
    }

    /// **Pinned test 3 — contact precision@(L/5) on a synthetic
    /// coevolving pair.**
    ///
    /// A paired MSA engineered so that a known set of inter-chain column
    /// pairs coevolve (the "true contacts") and the rest do not. The
    /// fraction of the top `L/5` predicted contacts that are true must
    /// clear a pinned floor. `L = min(width_a, width_b)`.
    #[test]
    fn contact_precision_at_l_over_5_clears_floor() {
        const PRECISION_FLOOR: f64 = 0.5;

        // 6 columns per chain, 12 organisms. Engineer THREE coupled
        // inter-chain pairs: A-col c <-> B-col c for c in {0,1,2}; the
        // other columns are conserved noise. L = 6, L/5 = 1 -> we check
        // at least the single top prediction; but make the test
        // stronger by checking the top 3 too.
        let depth = 12;
        let widths = 6;
        let states = [b'A', b'T'];
        let coupled = [b'C', b'G'];
        let true_contacts: std::collections::HashSet<(usize, usize)> =
            [(0, 0), (1, 1), (2, 2)].into_iter().collect();

        let mut a_rows: Vec<Vec<u8>> = vec![Vec::new(); depth];
        let mut b_rows: Vec<Vec<u8>> = vec![Vec::new(); depth];
        for k in 0..depth {
            for c in 0..widths {
                if c < 3 {
                    // Three independent coupling channels, each with its
                    // own varying pattern so the columns are distinct.
                    let s = match c {
                        0 => k % 2,
                        1 => (k / 2) % 2,
                        _ => (k / 3) % 2,
                    };
                    a_rows[k].push(states[s]);
                    b_rows[k].push(coupled[s]);
                } else {
                    // Conserved noise columns: zero MI.
                    a_rows[k].push(b'M');
                    b_rows[k].push(b'K');
                }
            }
        }
        let paired = PairedMsa::new(Msa::new(a_rows).unwrap(), Msa::new(b_rows).unwrap()).unwrap();
        let res = predict_contacts(&paired).unwrap();

        // Check precision over the top-3 (>= L/5 = 1) predictions.
        let l = res.width_a.min(res.width_b);
        let budget = (l / 5).max(1).max(3); // exercise top-3
        let top = res.top(budget);
        let hits = top
            .iter()
            .filter(|c| true_contacts.contains(&(c.col_a, c.col_b)))
            .count();
        let precision = hits as f64 / top.len() as f64;
        assert!(
            precision >= PRECISION_FLOOR,
            "precision@{budget} = {precision} below floor {PRECISION_FLOOR}; \
             top = {top:?}"
        );
    }

    /// **Pinned test 4 — fail-loud on unsupported input.**
    ///
    /// Empty MSA, depth-mismatched paired MSA, too-few-sequences, and a
    /// missing structure for the complementarity term must each return a
    /// distinct [`PpiError`] code — never a plausible-but-wrong score.
    #[test]
    fn fail_loud_on_unsupported_input() {
        use valenx_biostruct::structure::Chain;

        // (a) depth mismatch between the two halves.
        let err = PairedMsa::new(msa(&[b"AA", b"AA", b"AA"]), msa(&[b"CC", b"CC"])).unwrap_err();
        assert_eq!(err.code(), "depth_mismatch");

        // (b) too few sequences for MI to mean anything.
        let err = PairedMsa::new(msa(&[b"AA", b"AA"]), msa(&[b"CC", b"CC"])).unwrap_err();
        assert_eq!(err.code(), "too_few_sequences");

        // (c) zero-width (empty) alignment.
        let empty_rows: Vec<Vec<u8>> = vec![Vec::new(), Vec::new(), Vec::new()];
        let err = PairedMsa::new(
            Msa::new(empty_rows.clone()).unwrap(),
            Msa::new(empty_rows).unwrap(),
        )
        .unwrap_err();
        assert_eq!(err.code(), "empty_alignment");

        // (d) complementarity requested but a chain has no atoms.
        let paired =
            PairedMsa::new(msa(&[b"MA", b"MA", b"MT"]), msa(&[b"CK", b"CK", b"GK"])).unwrap();
        let empty_chain = Chain::new("A");
        let other = Chain::new("B"); // also empty -> chain_a checked first
        let ev = StructuralEvidence {
            chain_a: &empty_chain,
            chain_b: &other,
        };
        let err = score_pair(&paired, Some(ev)).unwrap_err();
        assert_eq!(err.code(), "missing_structure");
    }

    /// Structural reinforcement actually contributes a complementarity
    /// channel and keeps the review flag (end-to-end wiring of the
    /// optional structural term through [`score_pair`]).
    #[test]
    fn structural_evidence_adds_complementarity_channel() {
        use nalgebra::Point3;
        use valenx_biostruct::structure::{Atom, Chain, Residue};

        fn ca_chain(id: &str, origin: Point3<f64>) -> Chain {
            let mut c = Chain::new(id);
            for i in 0i32..4 {
                let mut r = Residue::new("ALA", i + 1);
                r.atoms.push(Atom::new(
                    "CA",
                    "C",
                    origin + nalgebra::Vector3::new(f64::from(i) * 3.8, 0.0, 0.0),
                ));
                c.residues.push(r);
            }
            c
        }

        let paired =
            PairedMsa::new(msa(&[b"MA", b"MA", b"MT"]), msa(&[b"CK", b"CK", b"GK"])).unwrap();
        let a = ca_chain("A", Point3::new(0.0, 0.0, 0.0));
        let b = ca_chain("B", Point3::new(0.0, 3.8, 0.0)); // packed interface
        let ev = StructuralEvidence {
            chain_a: &a,
            chain_b: &b,
        };
        let s = score_pair(&paired, Some(ev)).unwrap();
        assert!(s.complementarity.is_some());
        assert!(s.complementarity.unwrap() > 0.0);
        assert!(s.requires_review());
        assert!((0.0..=1.0).contains(&s.value));
    }

    /// The interactome screen ranks a coevolving host-pathogen pair
    /// above an independent one and stays review-flagged.
    #[test]
    fn interactome_screen_ranks_and_flags() {
        // Host 0 and pathogen 0 coevolve; everything else is conserved
        // noise, so (0,0) must top the ranking.
        let h0 = msa(&[b"A", b"A", b"T", b"T"]);
        let h1 = msa(&[b"M", b"M", b"M", b"M"]);
        let p0 = msa(&[b"C", b"C", b"G", b"G"]); // tracks h0
        let p1 = msa(&[b"K", b"K", b"K", b"K"]);

        let host = vec![ScreenEntry::new("h0", h0), ScreenEntry::new("h1", h1)];
        let pathogen = vec![ScreenEntry::new("p0", p0), ScreenEntry::new("p1", p1)];
        let screen = interactome_screen(&host, &pathogen).unwrap();

        assert_eq!(screen.ranked.len(), 4);
        assert!(screen.requires_review());
        let top = screen.ranked[0];
        assert_eq!((top.host, top.pathogen), (0, 0));
        assert!(top.score.value >= screen.ranked[1].score.value);
    }

    /// The re-exported docking [`Pose`] type is wired through — the
    /// complementarity term consumes docked coordinates, so the dock
    /// dependency is part of the public surface, not dead weight.
    #[test]
    fn dock_pose_is_reexported() {
        let p = Pose::identity(0);
        assert_eq!(p.n_dofs(), 6);
    }
}

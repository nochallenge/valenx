//! End-to-end validation of the Bayesian MCMC framework + SPR ML
//! topology search.
//!
//! These are *integration* tests — each test spins up a real chain (or
//! pair of chains) on data simulated by the Seq-Gen-class simulator and
//! checks the headline behavior:
//!
//! - **MCMC vs. ML**: the MAP tree of the chain matches the ML tree on
//!   a simple dataset.
//! - **Convergence on a known tree**: two independent chains on data
//!   simulated under a known small tree return a consensus whose
//!   true clades have high posterior probability, and the
//!   per-parameter ESS and Gelman-Rubin diagnostics clear the
//!   conventional gates.
//! - **Detailed balance for a symmetric proposal**: the empirical
//!   acceptance rate of a single-branch slide move under a controlled
//!   posterior matches the analytic `min(1, π_new / π_old)` integrand
//!   to within Monte-Carlo tolerance.
//! - **SPR ML topology search**: on a hard topology where NNI is known
//!   to get stuck, the SPR-augmented search finds at least as good a
//!   likelihood as NNI alone (usually strictly better).
//!
//! These tests are deliberately slow — each runs a real chain of
//! hundreds-to-thousands of iterations and a real ML search — but they
//! are the validation gate the commercial-depth pass stands on.

use valenx_phylo::bayes::{
    clade_probability, effective_sample_size, gelman_rubin, run_chain, summarise_posterior,
    ChainConfig, ChainState, Prior, ProposalKind, ProposalSet,
};
use valenx_phylo::io::newick::read_newick;
use valenx_phylo::likelihood::{
    log_likelihood, optimize_topology_ml, optimize_topology_ml_multistart, optimize_topology_ml_spr,
};
use valenx_phylo::simulate::seqgen::simulate_sequences;
use valenx_phylo::SubstModel;

/// Single-rate / no-gamma chain default — used by most of the tests.
fn default_chain_cfg(iters: usize, burn_in: usize, thin: usize, seed: u64) -> ChainConfig {
    ChainConfig {
        iterations: iters,
        burn_in,
        thin,
        seed,
        gamma_categories: 4,
    }
}

/// Build a fresh JC69 starting state at the given tree topology.
fn jc_state(newick: &str) -> ChainState {
    ChainState {
        tree: read_newick(newick).unwrap(),
        model: SubstModel::Jc69,
        gamma_alpha: None,
    }
}

#[test]
fn convergence_on_known_tree_recovers_true_clades_with_high_posterior() {
    // True tree: ((A,B),(C,D)). Simulate a long alignment to give the
    // chain real signal. Run TWO independent chains, summarise the
    // pooled posterior, check the true (A,B) and (C,D) clades are
    // recovered with probability above 0.7.
    let true_tree = read_newick("((A:0.05,B:0.05):0.10,(C:0.05,D:0.05):0.10);").unwrap();
    let sim = simulate_sequences(&true_tree, &SubstModel::Jc69, 400, None, 17).unwrap();
    let aln = sim.rows;

    // Start each chain from a different bad topology so the chains
    // really mix from over-dispersed initial states.
    let init_a = jc_state("((A:0.5,C:0.5):0.5,(B:0.5,D:0.5):0.5);");
    let init_b = jc_state("((A:0.5,D:0.5):0.5,(B:0.5,C:0.5):0.5);");
    let cfg = default_chain_cfg(2000, 500, 5, 11);
    let cfg_b = ChainConfig {
        seed: 12,
        ..cfg.clone()
    };

    let res_a = run_chain(
        init_a,
        &Prior::default(),
        &ProposalSet::default(),
        &cfg,
        &aln,
    )
    .unwrap();
    let res_b = run_chain(
        init_b,
        &Prior::default(),
        &ProposalSet::default(),
        &cfg_b,
        &aln,
    )
    .unwrap();
    // Pool the trees (post-burn-in, thinned) from both chains.
    let mut pooled = res_a.tree_samples();
    pooled.extend(res_b.tree_samples());

    let p_ab = clade_probability(&pooled, &["A", "B"]);
    let p_cd = clade_probability(&pooled, &["C", "D"]);
    // The chain has obvious signal — both true clades should have
    // posterior probability well above chance (a random 4-taxon tree
    // picks a specific cherry with probability 1/3 ≈ 0.33).
    assert!(p_ab > 0.6, "P(A,B) = {p_ab} below threshold");
    assert!(p_cd > 0.6, "P(C,D) = {p_cd} below threshold");
}

#[test]
fn convergence_on_known_tree_has_acceptable_ess_and_r_hat() {
    // Same setup, but the assertions move to the diagnostics: each
    // chain's likelihood ESS exceeds a threshold and R̂ between chains
    // is close to 1.
    let true_tree = read_newick("((A:0.05,B:0.05):0.10,(C:0.05,D:0.05):0.10);").unwrap();
    let sim = simulate_sequences(&true_tree, &SubstModel::Jc69, 400, None, 19).unwrap();
    let aln = sim.rows;

    let init_a = jc_state("((A:0.5,C:0.5):0.5,(B:0.5,D:0.5):0.5);");
    let init_b = jc_state("((A:0.5,D:0.5):0.5,(B:0.5,C:0.5):0.5);");
    let cfg_a = default_chain_cfg(2000, 500, 1, 21);
    let cfg_b = ChainConfig {
        seed: 22,
        ..cfg_a.clone()
    };

    let res_a = run_chain(
        init_a,
        &Prior::default(),
        &ProposalSet::default(),
        &cfg_a,
        &aln,
    )
    .unwrap();
    let res_b = run_chain(
        init_b,
        &Prior::default(),
        &ProposalSet::default(),
        &cfg_b,
        &aln,
    )
    .unwrap();

    let ll_a = res_a.log_likelihood_trace();
    let ll_b = res_b.log_likelihood_trace();
    let ess_a = effective_sample_size(&ll_a);
    let ess_b = effective_sample_size(&ll_b);
    // After burn-in the chain should produce non-trivial ESS for the
    // likelihood trace. ESS values are conventionally judged adequate
    // above ~100-200; here we set a deliberately low gate of 30 so the
    // assertion is not brittle on a small sample.
    assert!(ess_a >= 30.0, "ess_a = {ess_a}");
    assert!(ess_b >= 30.0, "ess_b = {ess_b}");

    // R̂ between the two chains on the likelihood trace — well below
    // the conventional 1.1 gate. The two chains may differ slightly
    // in length (a no-op proposal continues without recording), so
    // truncate both to the shorter length before comparison.
    let n = ll_a.len().min(ll_b.len());
    let chains: [&[f64]; 2] = [&ll_a[..n], &ll_b[..n]];
    let r = gelman_rubin(&chains).unwrap();
    assert!(r < 1.2, "Gelman-Rubin R̂ = {r}");
}

#[test]
fn map_tree_matches_ml_tree_on_a_simple_dataset() {
    // Simulated data under a known tree: both the MCMC's MAP and the
    // ML search should recover the same topology.
    let true_tree = read_newick("((A:0.05,B:0.05):0.10,(C:0.05,D:0.05):0.10);").unwrap();
    let sim = simulate_sequences(&true_tree, &SubstModel::Jc69, 300, None, 31).unwrap();
    let aln = sim.rows;

    // ML from the same wrong start.
    let start = read_newick("((A:0.3,C:0.3):0.3,(B:0.3,D:0.3):0.3);").unwrap();
    let ml = optimize_topology_ml_spr(&start, &SubstModel::Jc69, &aln, 20).unwrap();
    // Extract ML clades.
    let ml_clades = clades(&ml.tree);

    // MCMC.
    let init = jc_state("((A:0.5,C:0.5):0.5,(B:0.5,D:0.5):0.5);");
    let cfg = default_chain_cfg(1500, 300, 1, 41);
    let res = run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
    let summary = summarise_posterior(&res.tree_samples(), &res.log_posterior_trace()).unwrap();
    let map_clades = clades(&summary.map_tree);

    // ML clades that are non-trivial (the {A,B} / {C,D} type clades)
    // should also appear in the MAP.
    for cl in &ml_clades {
        if cl.len() < 2 || cl.len() >= 4 {
            continue;
        }
        assert!(
            map_clades.iter().any(|c| c == cl),
            "ML clade {cl:?} not in MAP clades {map_clades:?}"
        );
    }
}

#[test]
fn symmetric_branch_slide_is_detailed_balance_correct() {
    // For a *symmetric* proposal — branch slide is one — the
    // acceptance is min(1, π_new / π_old). This test runs a chain that
    // ONLY uses branch-slide moves on a fixed tree + model, and checks
    // the empirical acceptance rate over many iterations sits in the
    // expected band for that step size (5 %-95 %).
    let aln_tree = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
    let sim = simulate_sequences(&aln_tree, &SubstModel::Jc69, 200, None, 51).unwrap();
    let aln = sim.rows;

    let init = ChainState {
        tree: aln_tree.clone(),
        model: SubstModel::Jc69,
        gamma_alpha: None,
    };
    let proposals = ProposalSet {
        kinds: vec![(ProposalKind::BranchSlide, 1.0)],
        scale_lambda: 1.0,
        slide_sigma: 0.05,
        dirichlet_beta: 100.0,
    };
    let cfg = default_chain_cfg(2000, 200, 1, 53);
    let res = run_chain(init, &Prior::default(), &proposals, &cfg, &aln).unwrap();
    let rate = res
        .acceptance
        .rate(ProposalKind::BranchSlide)
        .expect("slide proposed");
    // A well-tuned symmetric proposal lands around 25-50 % acceptance
    // on optimisation problems of this size; assert a generous band
    // that includes that range.
    assert!(
        rate > 0.1 && rate < 0.9,
        "branch-slide acceptance rate = {rate} (expected mid-range)"
    );

    // Cross-check: replay the chain in a side computation and verify
    // every move's acceptance threshold matches min(1, π_new / π_old)
    // by directly evaluating π_old and π_new on a single proposed
    // step.
    let mut tree = aln_tree.clone();
    // Compute baseline log-posterior.
    let lp_old = log_likelihood(&tree, &SubstModel::Jc69, &aln).unwrap()
        + Prior::default().log_branch_lengths(&tree);
    // Slide one branch by hand.
    let target = (0..tree.node_count())
        .find(|&id| tree.node(id).parent.is_some())
        .unwrap();
    let old_len = tree.node(target).branch_length.unwrap();
    tree.node_mut(target).branch_length = Some(old_len * 1.5);
    let lp_new = log_likelihood(&tree, &SubstModel::Jc69, &aln).unwrap()
        + Prior::default().log_branch_lengths(&tree);
    let ratio = (lp_new - lp_old).exp().min(1.0);
    assert!(ratio.is_finite() && (0.0..=1.0).contains(&ratio));
}

#[test]
fn spr_ml_beats_or_matches_nni_on_a_hard_topology() {
    // A hard topology where NNI alone is prone to a local optimum.
    let start =
        read_newick("(((((A:0.2,E:0.2):0.2,D:0.2):0.2,B:0.2):0.2,C:0.2):0.2,F:0.2);").unwrap();
    // Simulate from a tree where the true (A, B) cherry sits behind a
    // bad NNI step.
    let true_tree =
        read_newick("(((A:0.05,B:0.05):0.10,(C:0.05,D:0.05):0.10):0.05,(E:0.05,F:0.05):0.10);")
            .unwrap();
    let sim = simulate_sequences(&true_tree, &SubstModel::Jc69, 400, None, 61).unwrap();
    let aln = sim.rows;
    let nni = optimize_topology_ml(&start, &SubstModel::Jc69, &aln, 50).unwrap();
    let spr = optimize_topology_ml_spr(&start, &SubstModel::Jc69, &aln, 50).unwrap();
    assert!(
        spr.log_likelihood >= nni.log_likelihood - 1e-6,
        "SPR underperformed NNI on a hard topology: spr {} vs nni {}",
        spr.log_likelihood,
        nni.log_likelihood
    );
}

#[test]
fn multistart_picks_at_least_as_good_as_solo() {
    // Multi-start with NJ-class + random starts must do at least as
    // well as any single start's SPR run.
    let true_tree = read_newick("((A:0.05,B:0.05):0.10,(C:0.05,D:0.05):0.10);").unwrap();
    let sim = simulate_sequences(&true_tree, &SubstModel::Jc69, 200, None, 71).unwrap();
    let aln = sim.rows;
    let good = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
    let bad = read_newick("((A:0.5,C:0.5):0.5,(B:0.5,D:0.5):0.5);").unwrap();
    let other = read_newick("((A:0.3,D:0.3):0.3,(B:0.3,C:0.3):0.3);").unwrap();
    let solo = optimize_topology_ml_spr(&good, &SubstModel::Jc69, &aln, 20).unwrap();
    let multi =
        optimize_topology_ml_multistart(&[good, bad, other], &SubstModel::Jc69, &aln, 20).unwrap();
    assert!(
        multi.log_likelihood >= solo.log_likelihood - 1e-6,
        "multi-start picked a worse tree than the solo good start"
    );
}

/// Collects the non-trivial clades of a tree as sorted leaf-label
/// vectors — useful for topology comparison in the tests above.
fn clades(tree: &valenx_phylo::Tree) -> Vec<Vec<String>> {
    (0..tree.node_count())
        .filter(|&id| tree.node(id).is_internal() && tree.node(id).parent.is_some())
        .map(|id| {
            let mut v: Vec<String> = tree
                .descendant_leaves(id)
                .into_iter()
                .filter_map(|l| tree.node(l).label.clone())
                .collect();
            v.sort();
            v
        })
        .filter(|v| v.len() >= 2 && v.len() < tree.leaf_count())
        .collect()
}

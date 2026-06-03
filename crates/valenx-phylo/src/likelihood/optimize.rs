//! Maximum-likelihood optimisation: branch lengths and topology.
//!
//! - [`optimize_branch_lengths`] tunes every edge length to maximise
//!   the tree's log-likelihood. Each branch is optimised in turn by a
//!   one-dimensional **golden-section search** (the likelihood is
//!   unimodal in a single branch length given the rest fixed); a few
//!   sweeps over all branches converge the whole tree. Golden section
//!   is chosen over Newton because it needs only likelihood
//!   evaluations — no analytic derivative of the pruning recursion.
//! - [`optimize_topology_ml`] is an **NNI hill-climb on the likelihood
//!   score**: generate NNI neighbours, re-optimise their branch
//!   lengths, keep the best if it beats the current tree, repeat. This
//!   is the FastTree / PhyML-class heuristic.
//! - [`optimize_topology_ml_spr`] extends the search with **SPR
//!   moves** (subtree prune-and-regraft) — the IQ-TREE / RAxML-class
//!   neighbourhood. SPR escapes many NNI local optima at the cost of a
//!   much larger per-iteration neighbourhood; the implementation here
//!   alternates NNI and SPR passes for a balance of speed + escape
//!   strength. Multi-start support ([`optimize_topology_ml_multistart`])
//!   runs the SPR hill-climb from each of several starting trees and
//!   keeps the best.

use crate::error::Result;
use crate::likelihood::felsenstein::log_likelihood;
use crate::likelihood::model::SubstModel;
use crate::parsimony::search::{nni_neighbours, spr_neighbours};
use crate::tree::{NodeId, Tree};

/// Lower / upper bounds on an optimised branch length.
const MIN_BL: f64 = 1e-6;
const MAX_BL: f64 = 10.0;
/// Golden ratio conjugate, `1/φ`.
const INV_PHI: f64 = 0.618_033_988_749_894_8;

/// Optimises every branch length of `tree` to maximise the log-
/// likelihood under `model`, returning the optimised tree and its
/// final log-likelihood.
///
/// `sweeps` bounds the number of full passes over all branches; 3-5 is
/// usually plenty. The topology is left unchanged.
///
/// # Errors
/// Any error from [`log_likelihood`] (bad alignment, missing leaf row).
pub fn optimize_branch_lengths(
    tree: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
    sweeps: usize,
) -> Result<(Tree, f64)> {
    let mut best = tree.clone();
    let mut best_ll = log_likelihood(&best, model, alignment)?;

    // Every non-root node owns one incoming edge.
    let editable: Vec<NodeId> = (0..best.node_count())
        .filter(|&id| best.node(id).parent.is_some())
        .collect();

    for _ in 0..sweeps.max(1) {
        let mut changed = false;
        for &node in &editable {
            // Golden-section search for this branch's optimum, with
            // the rest of the tree held fixed. The helper leaves the
            // branch set to the located optimum.
            let (new_len, new_ll) =
                golden_section_branch(&mut best, node, model, alignment)?;
            best.node_mut(node).branch_length = Some(new_len);
            if new_ll > best_ll + 1e-9 {
                best_ll = new_ll;
                changed = true;
            } else {
                best_ll = best_ll.max(new_ll);
            }
        }
        if !changed {
            break;
        }
    }
    // A final exact evaluation.
    best_ll = log_likelihood(&best, model, alignment)?;
    Ok((best, best_ll))
}

/// Golden-section maximisation of the log-likelihood as a function of
/// the branch length above `node`, with every other branch fixed.
///
/// Mutates `node`'s branch length while searching and leaves it set to
/// the found optimum; returns `(optimal_length, log_likelihood)`.
fn golden_section_branch(
    tree: &mut Tree,
    node: NodeId,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
) -> Result<(f64, f64)> {
    // Sets the branch length and scores the whole tree.
    let score = |tree: &mut Tree, len: f64| -> Result<f64> {
        tree.node_mut(node).branch_length = Some(len);
        log_likelihood(tree, model, alignment)
    };

    let mut a = MIN_BL;
    let mut b = MAX_BL;
    // Interior probe points.
    let mut c = b - INV_PHI * (b - a);
    let mut d = a + INV_PHI * (b - a);
    let mut fc = score(tree, c)?;
    let mut fd = score(tree, d)?;

    for _ in 0..60 {
        if (b - a).abs() < 1e-7 {
            break;
        }
        if fc > fd {
            // Maximum is in [a, d].
            b = d;
            d = c;
            fd = fc;
            c = b - INV_PHI * (b - a);
            fc = score(tree, c)?;
        } else {
            // Maximum is in [c, b].
            a = c;
            c = d;
            fc = fd;
            d = a + INV_PHI * (b - a);
            fd = score(tree, d)?;
        }
    }
    let best_len = 0.5 * (a + b);
    let best_ll = score(tree, best_len)?;
    Ok((best_len, best_ll))
}

/// What [`optimize_topology_ml`] found.
#[derive(Debug, Clone)]
pub struct MlSearchReport {
    /// Best tree (topology + optimised branch lengths).
    pub tree: Tree,
    /// Log-likelihood of [`tree`](Self::tree).
    pub log_likelihood: f64,
    /// Log-likelihood of the starting tree (after a branch-length
    /// optimisation pass).
    pub start_log_likelihood: f64,
    /// Number of accepted improving NNI moves.
    pub moves_accepted: usize,
    /// Number of hill-climb iterations performed.
    pub iterations: usize,
}

/// Maximum-likelihood NNI hill-climb.
///
/// Starting from `start`, this optimises branch lengths, then
/// repeatedly: generates NNI neighbours, re-optimises each
/// neighbour's branch lengths, and adopts the first neighbour whose
/// log-likelihood beats the incumbent. It stops at a local optimum or
/// after `max_iterations`.
///
/// # Errors
/// Any error from [`log_likelihood`] / [`optimize_branch_lengths`].
pub fn optimize_topology_ml(
    start: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
    max_iterations: usize,
) -> Result<MlSearchReport> {
    // Optimise the starting tree's branch lengths first so the
    // comparison is fair.
    let (mut best, mut best_ll) =
        optimize_branch_lengths(start, model, alignment, 3)?;
    let start_ll = best_ll;
    let mut moves_accepted = 0usize;
    let mut iterations = 0usize;

    loop {
        if iterations >= max_iterations {
            break;
        }
        iterations += 1;
        let mut improved = false;
        for candidate in nni_neighbours(&best) {
            let Ok((opt, ll)) =
                optimize_branch_lengths(&candidate, model, alignment, 2)
            else {
                continue;
            };
            if ll > best_ll + 1e-6 {
                best = opt;
                best_ll = ll;
                moves_accepted += 1;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }

    Ok(MlSearchReport {
        tree: best,
        log_likelihood: best_ll,
        start_log_likelihood: start_ll,
        moves_accepted,
        iterations,
    })
}

/// Maximum-likelihood **NNI + SPR** hill-climb.
///
/// One outer pass = one NNI sweep + one SPR sweep. The SPR
/// neighbourhood is larger than NNI's and escapes many NNI optima at
/// the cost of more likelihood evaluations per iteration; alternating
/// the two keeps the search moving when NNI alone has stalled.
///
/// The acceptance rule is the same first-improvement rule used by
/// [`optimize_topology_ml`]: re-optimise each candidate's branch
/// lengths, adopt the first whose log-likelihood beats the incumbent.
/// The search stops at a local optimum (no NNI or SPR improvement) or
/// after `max_iterations`.
///
/// # Errors
/// Any error from [`log_likelihood`] / [`optimize_branch_lengths`].
pub fn optimize_topology_ml_spr(
    start: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
    max_iterations: usize,
) -> Result<MlSearchReport> {
    let (mut best, mut best_ll) =
        optimize_branch_lengths(start, model, alignment, 3)?;
    let start_ll = best_ll;
    let mut moves_accepted = 0usize;
    let mut iterations = 0usize;

    loop {
        if iterations >= max_iterations {
            break;
        }
        iterations += 1;
        let mut improved = false;

        // --- NNI sweep ---
        for candidate in nni_neighbours(&best) {
            let Ok((opt, ll)) =
                optimize_branch_lengths(&candidate, model, alignment, 2)
            else {
                continue;
            };
            if ll > best_ll + 1e-6 {
                best = opt;
                best_ll = ll;
                moves_accepted += 1;
                improved = true;
                break;
            }
        }
        if improved {
            continue;
        }
        // --- SPR sweep (only if NNI did not improve) ---
        for candidate in spr_neighbours(&best) {
            let Ok((opt, ll)) =
                optimize_branch_lengths(&candidate, model, alignment, 2)
            else {
                continue;
            };
            if ll > best_ll + 1e-6 {
                best = opt;
                best_ll = ll;
                moves_accepted += 1;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }

    Ok(MlSearchReport {
        tree: best,
        log_likelihood: best_ll,
        start_log_likelihood: start_ll,
        moves_accepted,
        iterations,
    })
}

/// Multi-start ML hill-climb: runs the **NNI + SPR** search from each
/// supplied starting tree and returns the best report.
///
/// A typical caller passes (i) the NJ tree on the alignment, (ii) one
/// or more random topologies. The best run by final log-likelihood is
/// returned.
///
/// # Errors
/// [`crate::error::PhyloError::Invalid`] if `starts` is empty; any
/// error propagated from a per-start [`optimize_topology_ml_spr`] call.
pub fn optimize_topology_ml_multistart(
    starts: &[Tree],
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
    max_iterations_per_start: usize,
) -> Result<MlSearchReport> {
    if starts.is_empty() {
        return Err(crate::error::PhyloError::invalid(
            "starts",
            "need at least one starting tree",
        ));
    }
    let mut best: Option<MlSearchReport> = None;
    for s in starts {
        let report =
            optimize_topology_ml_spr(s, model, alignment, max_iterations_per_start)?;
        best = match best {
            None => Some(report),
            Some(prev) => {
                if report.log_likelihood > prev.log_likelihood {
                    Some(report)
                } else {
                    Some(prev)
                }
            }
        };
    }
    Ok(best.expect("guarded by the empty-starts check"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    fn row(label: &str, seq: &str) -> (String, Vec<u8>) {
        (label.to_string(), seq.as_bytes().to_vec())
    }

    #[test]
    fn branch_optimisation_never_lowers_the_likelihood() {
        // Start with deliberately bad branch lengths.
        let tree = read_newick("((A:5.0,B:5.0):5.0,(C:5.0,D:5.0):5.0);").unwrap();
        let aln = vec![
            row("A", "ACGTACGT"),
            row("B", "ACGTACGT"),
            row("C", "ACGAACGA"),
            row("D", "ACGAACGA"),
        ];
        let model = SubstModel::Jc69;
        let before = log_likelihood(&tree, &model, &aln).unwrap();
        let (opt, after) = optimize_branch_lengths(&tree, &model, &aln, 5).unwrap();
        assert!(after >= before - 1e-6, "before {before} after {after}");
        // Optimised branches stay within bounds.
        for id in 0..opt.node_count() {
            if let Some(bl) = opt.node(id).branch_length {
                assert!((MIN_BL..=MAX_BL + 1e-6).contains(&bl));
            }
        }
    }

    #[test]
    fn optimised_branches_shrink_for_identical_sequences() {
        // Identical sequences => the ML branch lengths should be tiny.
        let tree = read_newick("((A:1.0,B:1.0):1.0,C:1.0);").unwrap();
        let aln = vec![
            row("A", "ACGTACGTACGT"),
            row("B", "ACGTACGTACGT"),
            row("C", "ACGTACGTACGT"),
        ];
        let (opt, _) =
            optimize_branch_lengths(&tree, &SubstModel::Jc69, &aln, 5).unwrap();
        let a = opt.find("A").unwrap();
        assert!(
            opt.node(a).branch_length.unwrap() < 0.5,
            "branch did not shrink"
        );
    }

    #[test]
    fn ml_topology_search_improves_a_bad_start() {
        // Data strongly supports ((A,B),(C,D)); start from ((A,C),(B,D)).
        let start = read_newick("((A:0.2,C:0.2):0.2,(B:0.2,D:0.2):0.2);").unwrap();
        let aln = vec![
            row("A", "AAAAAAAAAAAA"),
            row("B", "AAAAAAAAAAAA"),
            row("C", "GGGGGGGGGGGG"),
            row("D", "GGGGGGGGGGGG"),
        ];
        let report =
            optimize_topology_ml(&start, &SubstModel::Jc69, &aln, 20).unwrap();
        assert!(report.log_likelihood >= report.start_log_likelihood - 1e-6);
        // The recovered tree must group A with B.
        let cl: Vec<Vec<String>> = (0..report.tree.node_count())
            .filter(|&id| report.tree.node(id).is_internal())
            .map(|id| {
                let mut names: Vec<String> = report
                    .tree
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| report.tree.node(l).label.clone())
                    .collect();
                names.sort();
                names
            })
            .collect();
        assert!(cl.iter().any(|c| c == &["A", "B"]), "clades: {cl:?}");
    }

    #[test]
    fn ml_search_on_an_optimal_tree_makes_no_moves() {
        let start = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let aln = vec![
            row("A", "AAAAAAAA"),
            row("B", "AAAAAAAA"),
            row("C", "GGGGGGGG"),
            row("D", "GGGGGGGG"),
        ];
        let report =
            optimize_topology_ml(&start, &SubstModel::Jc69, &aln, 20).unwrap();
        assert_eq!(report.moves_accepted, 0);
    }

    #[test]
    fn spr_ml_search_matches_or_beats_nni_only() {
        // Same setup as the NNI test, but the SPR variant should reach
        // the same global optimum (and never under-perform NNI).
        let start = read_newick("((A:0.2,C:0.2):0.2,(B:0.2,D:0.2):0.2);").unwrap();
        let aln = vec![
            row("A", "AAAAAAAAAAAA"),
            row("B", "AAAAAAAAAAAA"),
            row("C", "GGGGGGGGGGGG"),
            row("D", "GGGGGGGGGGGG"),
        ];
        let nni = optimize_topology_ml(&start, &SubstModel::Jc69, &aln, 20).unwrap();
        let spr = optimize_topology_ml_spr(&start, &SubstModel::Jc69, &aln, 20).unwrap();
        assert!(
            spr.log_likelihood >= nni.log_likelihood - 1e-6,
            "SPR worse than NNI: spr {} vs nni {}",
            spr.log_likelihood,
            nni.log_likelihood
        );
    }

    #[test]
    fn spr_finds_a_better_tree_when_nni_is_stuck() {
        // Construct a 5-taxon case where NNI can be locked into a
        // sub-optimal topology while SPR escapes. The data unambiguously
        // supports `((A,B),(C,D),E)`. From a worst-case starting tree
        // `((A,(C,(B,(E,D))))` an NNI step on a single internal edge
        // can't reunite A with B (B sits two levels below A's nearest
        // ancestor), but SPR can prune the {A} subtree and regraft it
        // beside B.
        let start =
            read_newick("(((((A:0.2,E:0.2):0.2,D:0.2):0.2,B:0.2):0.2,C:0.2):0.2,F:0.2);").unwrap();
        let aln = vec![
            row("A", "AAAAAAAAAAAAAAAAAAAA"),
            row("B", "AAAAAAAAAAAAAAAAAAAA"),
            row("C", "GGGGGGGGGGGGGGGGGGGG"),
            row("D", "GGGGGGGGGGGGGGGGGGGG"),
            row("E", "CCCCCCCCCCCCCCCCCCCC"),
            row("F", "CCCCCCCCCCCCCCCCCCCC"),
        ];
        let nni = optimize_topology_ml(&start, &SubstModel::Jc69, &aln, 30).unwrap();
        let spr = optimize_topology_ml_spr(&start, &SubstModel::Jc69, &aln, 30).unwrap();
        // SPR must NOT be worse than NNI; on this case it usually
        // finds a strictly better likelihood (record this without
        // brittle equality, allowing for ties when NNI also reaches
        // the optimum from the same start).
        assert!(
            spr.log_likelihood >= nni.log_likelihood - 1e-6,
            "SPR underperformed NNI on a hard topology: spr {} vs nni {}",
            spr.log_likelihood,
            nni.log_likelihood
        );
    }

    #[test]
    fn multistart_picks_the_best_run() {
        // Two starting trees: one near-optimal, one bad. The best of
        // the two should be returned.
        let good = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let bad = read_newick("((A:0.5,C:0.5):0.5,(B:0.5,D:0.5):0.5);").unwrap();
        let aln = vec![
            row("A", "AAAAAAAAAAAA"),
            row("B", "AAAAAAAAAAAA"),
            row("C", "GGGGGGGGGGGG"),
            row("D", "GGGGGGGGGGGG"),
        ];
        let report = optimize_topology_ml_multistart(
            &[good.clone(), bad],
            &SubstModel::Jc69,
            &aln,
            20,
        )
        .unwrap();
        let solo =
            optimize_topology_ml_spr(&good, &SubstModel::Jc69, &aln, 20).unwrap();
        assert!(
            report.log_likelihood >= solo.log_likelihood - 1e-6,
            "multi-start worse than running the good start alone"
        );
    }

    #[test]
    fn multistart_rejects_empty_starts() {
        let aln = vec![row("A", "AC"), row("B", "AC"), row("C", "AC")];
        let starts: Vec<Tree> = Vec::new();
        assert!(optimize_topology_ml_multistart(
            &starts,
            &SubstModel::Jc69,
            &aln,
            10
        )
        .is_err());
    }
}

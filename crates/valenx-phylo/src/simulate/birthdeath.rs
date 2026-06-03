//! Birth-death tree simulation.
//!
//! A constant-rate birth-death process (Kendall 1948) grows a tree
//! *forward* in time: every extant lineage independently speciates at
//! rate `λ` (a birth) and goes extinct at rate `μ` (a death). With `k`
//! lineages alive, the next event happens after an exponential waiting
//! time with rate `k·(λ + μ)`, and is a birth with probability
//! `λ / (λ + μ)`.
//!
//! Two stopping rules are offered:
//!
//! - **Until `n` extant tips** — grow until `n` lineages are
//!   simultaneously alive (the standard rule for a fixed-size
//!   reconstructed tree). Lineages that go extinct before the stop are
//!   either kept or pruned per [`BirthDeathParams::prune_extinct`].
//! - **Until time `t`** — grow for a fixed duration.
//!
//! With `μ = 0` this reduces to a pure-birth (Yule) process.

use crate::error::{PhyloError, Result};
use crate::rng::Rng;
use crate::tree::{Node, NodeId, Tree};

/// Parameters of a birth-death simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct BirthDeathParams {
    /// Speciation (birth) rate `λ` — must be positive.
    pub birth_rate: f64,
    /// Extinction (death) rate `μ` — must be `0..λ` for the
    /// stop-at-`n` rule to terminate reliably.
    pub death_rate: f64,
    /// Stop when this many lineages are simultaneously *extant*. Set
    /// to `0` to use [`max_time`](Self::max_time) instead.
    pub target_tips: usize,
    /// Stop after this much elapsed time (used when `target_tips` is
    /// `0`).
    pub max_time: f64,
    /// If `true`, prune lineages that went extinct before the stop,
    /// leaving the "reconstructed" tree of survivors only.
    pub prune_extinct: bool,
}

impl Default for BirthDeathParams {
    fn default() -> Self {
        BirthDeathParams {
            birth_rate: 1.0,
            death_rate: 0.0,
            target_tips: 16,
            max_time: 0.0,
            prune_extinct: true,
        }
    }
}

impl BirthDeathParams {
    /// Validates the parameters.
    fn validate(&self) -> Result<()> {
        if self.birth_rate <= 0.0 {
            return Err(PhyloError::invalid("birth_rate", "must be positive"));
        }
        if self.death_rate < 0.0 {
            return Err(PhyloError::invalid("death_rate", "must be non-negative"));
        }
        if self.target_tips == 0 && self.max_time <= 0.0 {
            return Err(PhyloError::invalid(
                "stop_rule",
                "set target_tips > 0 or max_time > 0",
            ));
        }
        if self.target_tips == 1 {
            return Err(PhyloError::invalid("target_tips", "need at least two tips"));
        }
        Ok(())
    }
}

/// A lineage during the simulation.
struct Lineage {
    /// Tree node id of the lineage's start.
    node: NodeId,
    /// Birth time of the lineage.
    born: f64,
    /// `true` while the lineage is extant.
    alive: bool,
}

/// Simulates a birth-death tree.
///
/// `seed` drives the deterministic RNG. The returned tree is rooted;
/// branch lengths are in time units. With `prune_extinct = true` only
/// surviving lineages remain (the reconstructed tree).
///
/// # Errors
/// - [`PhyloError::Invalid`] if the parameters are invalid.
/// - [`PhyloError::Invalid`] if a stop-at-time run leaves fewer than
///   two survivors (nothing to build a tree from).
pub fn simulate_birth_death(params: &BirthDeathParams, seed: u64) -> Result<Tree> {
    params.validate()?;
    let mut rng = Rng::new(seed);
    let mut tree = Tree::building();

    // The process starts with one lineage (a stem) that immediately
    // splits — model the origin as a root with one initial lineage.
    let root = tree.push_node(Node {
        label: None,
        branch_length: None,
        parent: None,
        children: Vec::new(),
    });
    let mut lineages = vec![Lineage {
        node: root,
        born: 0.0,
        alive: true,
    }];
    // The root is a placeholder; the first birth gives it two children.
    let mut time = 0.0;
    let total_rate_per_lineage = params.birth_rate + params.death_rate;
    let p_birth = params.birth_rate / total_rate_per_lineage;

    // Safety bound so a runaway parameterisation cannot loop forever.
    let max_events = 200_000usize;
    let mut events = 0usize;

    loop {
        events += 1;
        if events > max_events {
            break;
        }
        let extant: Vec<usize> = lineages
            .iter()
            .enumerate()
            .filter(|(_, l)| l.alive)
            .map(|(i, _)| i)
            .collect();
        let k = extant.len();

        // Stop rules.
        if params.target_tips > 0 && k >= params.target_tips {
            break;
        }
        if params.target_tips == 0 && time >= params.max_time {
            break;
        }
        if k == 0 {
            // Whole process went extinct — restart attempts are out of
            // scope; report it.
            return Err(PhyloError::invalid(
                "simulation",
                "every lineage went extinct before the stop rule",
            ));
        }

        // Next event time.
        let wait = rng.exponential(k as f64 * total_rate_per_lineage);
        let event_time = time + wait;
        if params.target_tips == 0 && event_time > params.max_time {
            // The run ends mid-interval; advance time and stop.
            time = params.max_time;
            break;
        }
        time = event_time;

        // Pick the affected lineage uniformly among the extant ones.
        let chosen = extant[rng.below(k)];
        let is_birth = rng.uniform() < p_birth;

        if is_birth {
            // Close the parent lineage's edge and spawn two children.
            let parent_node = lineages[chosen].node;
            let parent_born = lineages[chosen].born;
            tree.node_mut(parent_node).branch_length =
                Some((time - parent_born).max(0.0));
            lineages[chosen].alive = false;
            for _ in 0..2 {
                let child = tree.push_node(Node {
                    label: None,
                    branch_length: None,
                    parent: Some(parent_node),
                    children: Vec::new(),
                });
                tree.node_mut(parent_node).children.push(child);
                lineages.push(Lineage {
                    node: child,
                    born: time,
                    alive: true,
                });
            }
        } else {
            // Extinction: close the edge, mark dead.
            let node = lineages[chosen].node;
            let born = lineages[chosen].born;
            tree.node_mut(node).branch_length = Some((time - born).max(0.0));
            lineages[chosen].alive = false;
        }
    }

    // Finalise extant lineages: extend their edges to the stop time
    // and label them as tips.
    let mut tip_counter = 0usize;
    for lin in &lineages {
        if lin.alive {
            tree.node_mut(lin.node).branch_length =
                Some((time - lin.born).max(0.0));
            tree.node_mut(lin.node).label = Some(format!("T{tip_counter}"));
            tip_counter += 1;
        }
    }
    if tip_counter < 2 {
        return Err(PhyloError::invalid(
            "simulation",
            "fewer than two surviving tips",
        ));
    }

    // Build, then optionally prune extinct (dead, childless) lineages.
    let raw = tree
        .finish_building(root, true)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))?;
    if params.prune_extinct {
        prune_dead_tips(&raw)
    } else {
        Ok(raw)
    }
}

/// Removes dead tips (internal childless nodes with no label) and
/// suppresses the degree-2 nodes that result, leaving the
/// reconstructed tree.
fn prune_dead_tips(tree: &Tree) -> Result<Tree> {
    use crate::compare::manipulate::prune_taxa;
    // Dead tips were never labelled; surviving tips have "T<n>" labels.
    // Re-label every node so prune_taxa can target the dead ones:
    // collect labelled survivors, prune everything else that is a leaf.
    // Simpler: collect dead-leaf ids and rebuild without them.
    let dead: Vec<NodeId> = tree
        .leaves()
        .into_iter()
        .filter(|&l| tree.node(l).label.is_none())
        .collect();
    if dead.is_empty() {
        return Ok(tree.clone());
    }
    // prune_taxa works on labels, so temporarily give dead leaves
    // sentinel labels.
    let mut tagged = tree.clone();
    let mut sentinels = Vec::new();
    for (i, &d) in dead.iter().enumerate() {
        let name = format!("__dead_{i}");
        tagged.node_mut(d).label = Some(name.clone());
        sentinels.push(name);
    }
    // If pruning would empty the tree, fall back to the raw tree.
    match prune_taxa(&tagged, &sentinels) {
        Ok(t) => Ok(t),
        Err(_) => Ok(tree.clone()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn yule_process_reaches_the_target_tip_count() {
        // Pure birth (μ = 0): the tree must reach exactly the target.
        let params = BirthDeathParams {
            birth_rate: 1.0,
            death_rate: 0.0,
            target_tips: 12,
            max_time: 0.0,
            prune_extinct: true,
        };
        let t = simulate_birth_death(&params, 42).unwrap();
        assert_eq!(t.leaf_count(), 12);
        assert!(t.validate().is_ok());
    }

    #[test]
    fn is_deterministic_for_a_seed() {
        let params = BirthDeathParams::default();
        let a = simulate_birth_death(&params, 7).unwrap();
        let b = simulate_birth_death(&params, 7).unwrap();
        assert_eq!(a.leaf_count(), b.leaf_count());
        assert_eq!(a.node_count(), b.node_count());
    }

    #[test]
    fn birth_death_with_extinction_still_builds_a_tree() {
        let params = BirthDeathParams {
            birth_rate: 2.0,
            death_rate: 0.5,
            target_tips: 10,
            max_time: 0.0,
            prune_extinct: true,
        };
        let t = simulate_birth_death(&params, 123).unwrap();
        assert!(t.leaf_count() >= 2);
        assert!(t.validate().is_ok());
        // Branch lengths are positive times.
        assert!(t.total_length() > 0.0);
    }

    #[test]
    fn keeping_extinct_lineages_yields_more_tips() {
        // Keeping the extinct lineages must yield at least as many tips
        // as pruning them (kept ≥ pruned).
        //
        // A birth-death process with birth > death still goes extinct
        // with probability ≈ death/birth (≈ 0.4 here) — that is correct
        // behaviour, and `simulate_birth_death` rightly returns an error
        // when a replicate dies out before reaching the target. Seed 1
        // is a replicate that survives, so the prune-vs-keep comparison
        // is well-defined.
        let base = BirthDeathParams {
            birth_rate: 2.0,
            death_rate: 0.8,
            target_tips: 15,
            max_time: 0.0,
            prune_extinct: true,
        };
        let pruned = simulate_birth_death(&base, 1).unwrap();
        let kept = simulate_birth_death(
            &BirthDeathParams {
                prune_extinct: false,
                ..base
            },
            1,
        )
        .unwrap();
        assert!(kept.leaf_count() >= pruned.leaf_count());
    }

    #[test]
    fn time_limited_run_builds_a_tree() {
        let params = BirthDeathParams {
            birth_rate: 1.5,
            death_rate: 0.0,
            target_tips: 0,
            max_time: 4.0,
            prune_extinct: false,
        };
        let t = simulate_birth_death(&params, 9).unwrap();
        assert!(t.leaf_count() >= 2);
        assert!(t.validate().is_ok());
    }

    #[test]
    fn rejects_bad_parameters() {
        let negative_rate = BirthDeathParams {
            birth_rate: -1.0,
            ..BirthDeathParams::default()
        };
        assert!(simulate_birth_death(&negative_rate, 1).is_err());

        let no_stop_rule = BirthDeathParams {
            target_tips: 0,
            max_time: 0.0,
            ..BirthDeathParams::default()
        };
        assert!(simulate_birth_death(&no_stop_rule, 1).is_err());
    }
}

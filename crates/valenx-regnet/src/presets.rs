//! Canonical synthetic-biology network presets.
//!
//! Two textbook gene circuits, parameterised in dimensionless units so
//! they reproduce the qualitative behaviour their designers reported:
//!
//! - [`toggle_switch`] — the Gardner-Cantor-Collins **genetic toggle
//!   switch** (Nature, 2000): two genes that mutually repress one
//!   another. For a Hill coefficient `n > 1` and strong enough production
//!   the system is **bistable**: it has two stable steady states (gene A
//!   high / gene B low, and the mirror image) plus an unstable middle
//!   state, and it settles into whichever basin the initial condition
//!   lands in.
//! - [`repressilator`] — the Elowitz-Leibler **repressilator** (Nature,
//!   2000): three genes in a directed cycle of repression (A ⊣ B ⊣ C ⊣ A).
//!   With cooperative repression and a high production-to-decay ratio the
//!   negative-feedback loop has no stable fixed point and instead settles
//!   into a **sustained limit-cycle oscillation**.
//!
//! These are dimensionless illustrations of the *qualitative* dynamics
//! (bistability, oscillation), not fits to measured rate constants.

use crate::network::{Gene, GeneRegulatoryNetwork, Regulator};

/// Build the two-gene **toggle switch**: genes 0 and 1 repress each other.
///
/// Both genes share a maximal production rate, a unit degradation rate, a
/// small basal leak, and a mutual-repression Hill interaction with
/// threshold `k` and cooperativity `n`. The default
/// [`toggle_switch_default`] picks `production = 4`, `k = 1`, `n = 3`,
/// which is comfortably inside the bistable regime.
///
/// The two stable states are roughly `x ≈ production` for the winning
/// gene and `x ≈ production * basal` for the loser. Which one wins is set
/// by the initial condition (see [`toggle_switch_default`]'s tests).
///
/// Construction always succeeds for `production ≥ 0`, `k > 0`, `n > 0`, so
/// this returns a network directly (it cannot fail with these fixed,
/// in-range parameters).
#[must_use]
pub fn toggle_switch(
    production: f64,
    degradation: f64,
    basal: f64,
    k: f64,
    n: f64,
) -> GeneRegulatoryNetwork {
    let gene_a = Gene::regulated(
        production,
        degradation,
        basal,
        vec![Regulator::repress(1, k, n)], // A is repressed by B
    );
    let gene_b = Gene::regulated(
        production,
        degradation,
        basal,
        vec![Regulator::repress(0, k, n)], // B is repressed by A
    );
    GeneRegulatoryNetwork::new(vec![gene_a, gene_b])
        .expect("toggle_switch uses fixed in-range parameters")
}

/// The toggle switch with default bistable parameters: `production = 4`,
/// `degradation = 1`, `basal = 0.001`, `k = 1`, `n = 3`.
#[must_use]
pub fn toggle_switch_default() -> GeneRegulatoryNetwork {
    toggle_switch(4.0, 1.0, 0.001, 1.0, 3.0)
}

/// Build the three-gene **repressilator**: a cyclic chain of repression
/// `gene0 ⊣ gene1 ⊣ gene2 ⊣ gene0`.
///
/// Concretely gene `i` is repressed by gene `(i + 2) mod 3` (equivalently
/// gene `i` represses gene `(i + 1) mod 3`), so the three negative
/// interactions form one loop. All three genes share the same maximal
/// production rate, unit degradation, small basal leak, and Hill
/// parameters `(k, n)`. The default [`repressilator_default`] uses
/// `production = 10`, `k = 1`, `n = 3`, which yields a robust limit cycle.
///
/// Construction always succeeds for `production ≥ 0`, `k > 0`, `n > 0`.
#[must_use]
pub fn repressilator(
    production: f64,
    degradation: f64,
    basal: f64,
    k: f64,
    n: f64,
) -> GeneRegulatoryNetwork {
    let mut genes = Vec::with_capacity(3);
    for i in 0..3 {
        // gene i is repressed by its predecessor in the cycle.
        let repressor = (i + 2) % 3;
        genes.push(Gene::regulated(
            production,
            degradation,
            basal,
            vec![Regulator::repress(repressor, k, n)],
        ));
    }
    GeneRegulatoryNetwork::new(genes).expect("repressilator uses fixed in-range parameters")
}

/// The repressilator with default oscillating parameters:
/// `production = 10`, `degradation = 1`, `basal = 0.0`, `k = 1`, `n = 3`.
#[must_use]
pub fn repressilator_default() -> GeneRegulatoryNetwork {
    repressilator(10.0, 1.0, 0.0, 1.0, 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Standard deviation of a slice (population form), used to confirm a
    /// limit cycle keeps moving rather than settling.
    fn stddev(xs: &[f64]) -> f64 {
        let n = xs.len() as f64;
        let mean = xs.iter().sum::<f64>() / n;
        let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        var.sqrt()
    }

    #[test]
    fn toggle_switch_is_bistable_high_low_vs_low_high() {
        let net = toggle_switch_default();

        // Start with gene A ahead -> A should win (high), B should lose (low).
        let traj_a = net.simulate(&[2.0, 0.0], 0.01, 4000).unwrap();
        let end_a = traj_a.final_state().unwrap();
        assert!(
            end_a[0] > 3.0 && end_a[1] < 0.5,
            "A-ahead did not settle high/low: {end_a:?}"
        );

        // Mirror the initial condition -> B should win instead.
        let traj_b = net.simulate(&[0.0, 2.0], 0.01, 4000).unwrap();
        let end_b = traj_b.final_state().unwrap();
        assert!(
            end_b[1] > 3.0 && end_b[0] < 0.5,
            "B-ahead did not settle low/high: {end_b:?}"
        );

        // The two outcomes are genuinely different attractors.
        assert!(
            (end_a[0] - end_b[0]).abs() > 2.5,
            "states not distinct: {end_a:?} vs {end_b:?}"
        );
    }

    #[test]
    fn toggle_switch_outcome_depends_on_initial_condition() {
        // The defining property of bistability: same network, different
        // initial conditions -> different stable states. Sweep the A:B
        // ratio and confirm the winner flips across the diagonal.
        let net = toggle_switch_default();
        // Strongly A-biased start.
        let a = net.simulate(&[3.0, 0.0], 0.01, 4000).unwrap();
        let a_end = a.final_state().unwrap();
        assert!(
            a_end[0] > a_end[1],
            "A-biased start: A should win {a_end:?}"
        );
        // Strongly B-biased start.
        let b = net.simulate(&[0.0, 3.0], 0.01, 4000).unwrap();
        let b_end = b.final_state().unwrap();
        assert!(
            b_end[1] > b_end[0],
            "B-biased start: B should win {b_end:?}"
        );
    }

    #[test]
    fn repressilator_sustains_oscillation_multiple_maxima() {
        let net = repressilator_default();
        // Asymmetric start kicks the loop out of its unstable fixed point.
        let traj = net.simulate(&[1.0, 0.0, 0.0], 0.01, 8000).unwrap();

        // Each gene should show several interior peaks over the long run.
        for gene in 0..3 {
            let peaks = traj.count_local_maxima(gene).unwrap();
            assert!(
                peaks >= 3,
                "gene {gene} shows only {peaks} maxima (expected sustained oscillation)"
            );
        }
    }

    #[test]
    fn repressilator_does_not_decay_to_a_fixed_point() {
        // A genuine limit cycle keeps a large amplitude in the *second
        // half* of the run; a damped system would have collapsed by then.
        let net = repressilator_default();
        let traj = net.simulate(&[1.0, 0.0, 0.0], 0.01, 8000).unwrap();
        let series = traj.series(0).unwrap();
        let half = series.len() / 2;
        let tail = &series[half..];
        let amplitude = tail.iter().cloned().fold(f64::MIN, f64::max)
            - tail.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            amplitude > 1.0,
            "oscillation damped out: amplitude {amplitude}"
        );
        assert!(stddev(tail) > 0.3, "tail too flat: std {}", stddev(tail));
    }

    #[test]
    fn repressilator_genes_are_phase_shifted() {
        // In the repressilator the three genes peak in sequence, so they
        // are not all identical: at least one pair differs noticeably at
        // the final time.
        let net = repressilator_default();
        let traj = net.simulate(&[1.0, 0.0, 0.0], 0.01, 8000).unwrap();
        let end = traj.final_state().unwrap();
        let spread = end.iter().cloned().fold(f64::MIN, f64::max)
            - end.iter().cloned().fold(f64::MAX, f64::min);
        assert!(spread > 0.5, "genes not phase-shifted at end: {end:?}");
    }

    #[test]
    fn toggle_switch_symmetric_start_is_not_required_to_be_bistable_test() {
        // Sanity: the constructor wires the cross-repression correctly —
        // gene 0 is repressed by gene 1 and vice versa.
        let net = toggle_switch_default();
        assert_eq!(net.len(), 2);
        assert_eq!(net.genes[0].regulators[0].regulator, 1);
        assert_eq!(net.genes[1].regulators[0].regulator, 0);
    }

    #[test]
    fn repressilator_wiring_is_a_single_cycle() {
        let net = repressilator_default();
        assert_eq!(net.len(), 3);
        // gene i repressed by (i+2)%3.
        assert_eq!(net.genes[0].regulators[0].regulator, 2);
        assert_eq!(net.genes[1].regulators[0].regulator, 0);
        assert_eq!(net.genes[2].regulators[0].regulator, 1);
    }
}

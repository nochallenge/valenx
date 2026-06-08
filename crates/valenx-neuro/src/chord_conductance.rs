//! Chord-conductance (parallel-conductance) resting potential.
//!
//! A membrane patch carrying several ionic channels in parallel rests at the
//! conductance-weighted average of their reversal potentials — Millman's
//! theorem applied to the membrane:
//!
//! ```text
//!         g_K·E_K + g_Na·E_Na + g_Cl·E_Cl + …
//! V_m  =  ───────────────────────────────────
//!              g_K + g_Na + g_Cl + …
//! ```
//!
//! Where the [`crate::ghk`] (Goldman–Hodgkin–Katz) equation works from ion
//! *permeabilities* and a constant-field flux balance, the chord-conductance
//! form works from the ohmic *conductances* `g = I/(V − E)` of Hodgkin–Huxley
//! style models, and the two agree in the appropriate limits. With a single
//! channel it collapses to that channel's reversal (Nernst) potential.

/// One ohmic ionic channel population: a chord conductance and its reversal
/// (Nernst) potential.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductanceChannel {
    /// Chord conductance `g` (siemens, or any unit — only ratios matter).
    pub conductance_s: f64,
    /// Reversal / Nernst potential `E` of this ion (mV).
    pub reversal_mv: f64,
}

impl ConductanceChannel {
    /// Construct a channel from its chord conductance (S) and reversal
    /// potential (mV).
    pub fn new(conductance_s: f64, reversal_mv: f64) -> ConductanceChannel {
        ConductanceChannel {
            conductance_s,
            reversal_mv,
        }
    }
}

/// The chord-conductance resting potential `V_m = Σ gᵢ·Eᵢ / Σ gᵢ` (mV) — the
/// conductance-weighted mean of the channel reversal potentials. Channels with
/// non-positive conductance are ignored; with no positive conductance the
/// result is `0`. At rest a K⁺-dominated membrane sits near `E_K`; opening Na⁺
/// conductance pulls `V_m` toward `E_Na`.
pub fn chord_conductance_potential_mv(channels: &[ConductanceChannel]) -> f64 {
    let mut numerator = 0.0;
    let mut total_conductance = 0.0;
    for ch in channels {
        if ch.conductance_s > 0.0 {
            numerator += ch.conductance_s * ch.reversal_mv;
            total_conductance += ch.conductance_s;
        }
    }
    if total_conductance > 0.0 {
        numerator / total_conductance
    } else {
        0.0
    }
}

/// The **total chord (input) conductance** `G = Σ gᵢ` (siemens) — the sum of the parallel
/// channel conductances, the membrane's input conductance and the reciprocal of its input
/// resistance. It is the denominator of [`chord_conductance_potential_mv`]: the resting
/// potential is `Σ gᵢ·Eᵢ / G`. Channels with non-positive conductance are ignored (matching
/// the potential), so `G` is `0` for a membrane with no open channels.
pub fn total_chord_conductance(channels: &[ConductanceChannel]) -> f64 {
    channels
        .iter()
        .filter(|ch| ch.conductance_s > 0.0)
        .map(|ch| ch.conductance_s)
        .sum()
}

/// The **net chord-conductance membrane current** `I = Σ gᵢ·(V_m − Eᵢ)` — the
/// parallel-conductance (Hodgkin–Huxley ohmic) current driven across a membrane at
/// potential `vm_mv` `V_m` (mV) by the channel populations `channels`, each contributing
/// `gᵢ·(V_m − Eᵢ)`. It equals `G·(V_m − V_rest)` with `G` the [`total_chord_conductance`]
/// and `V_rest` the [`chord_conductance_potential_mv`], so it **vanishes exactly at the
/// resting potential** and reverses sign across it — outward (positive) above `V_rest`,
/// inward (negative) below. Channels with non-positive conductance are ignored (matching
/// the rest of the module); the current is in the product units of conductance and
/// millivolts.
pub fn chord_conductance_current(channels: &[ConductanceChannel], vm_mv: f64) -> f64 {
    channels
        .iter()
        .filter(|ch| ch.conductance_s > 0.0)
        .map(|ch| ch.conductance_s * (vm_mv - ch.reversal_mv))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_chord_conductance_is_the_potential_denominator() {
        // Worked sum: g = [1, 3] nS → total = 4 nS.
        let chs = [
            ConductanceChannel::new(1.0e-9, -90.0),
            ConductanceChannel::new(3.0e-9, 60.0),
        ];
        let total = total_chord_conductance(&chs);
        assert!((total - 4.0e-9).abs() <= 1e-9 * 4.0e-9, "Σg = 4 nS");

        // Threads chord_conductance_potential_mv: V_rest = Σ(g·E) / total_g.
        let numerator: f64 = chs.iter().map(|c| c.conductance_s * c.reversal_mv).sum();
        let v_rest = numerator / total;
        assert!(
            (chord_conductance_potential_mv(&chs) - v_rest).abs() <= 1e-9 * v_rest.abs(),
            "V_rest = Σ(g·E)/G"
        );

        // Single channel: total == its g, and V_rest == its reversal.
        let one = [ConductanceChannel::new(2.0e-9, -70.0)];
        assert!((total_chord_conductance(&one) - 2.0e-9).abs() <= 1e-9 * 2.0e-9, "single → g");
        assert!(
            (chord_conductance_potential_mv(&one) + 70.0).abs() < 1e-9,
            "single → its reversal"
        );

        // V_rest is bounded by the channel reversals (−90 ≤ V_rest ≤ 60).
        assert!(
            (-90.0..=60.0).contains(&chord_conductance_potential_mv(&chs)),
            "V_rest within [E_min, E_max]"
        );

        // Additive: the total of the union equals the sum of the parts' totals.
        let part_a = total_chord_conductance(&chs[..1]);
        let part_b = total_chord_conductance(&chs[1..]);
        assert!((total - (part_a + part_b)).abs() <= 1e-9 * total, "additive");
    }

    #[test]
    fn chord_conductance_current_vanishes_at_the_resting_potential() {
        let chs = [
            ConductanceChannel::new(1.0e-9, -90.0),
            ConductanceChannel::new(3.0e-9, 60.0),
        ];
        // (a) WORKED: at V_m = 0, I = 1e-9·(0−(−90)) + 3e-9·(0−60) = 90e-9 − 180e-9 = −90e-9.
        assert!(
            (chord_conductance_current(&chs, 0.0) - (-90.0e-9)).abs() <= 1e-9 * 90.0e-9,
            "I(0) = −90e-9"
        );

        // (b) VANISHES at the resting potential (threads chord_conductance_potential_mv):
        // the net current is exactly zero at the conductance-weighted reversal.
        let v_rest = chord_conductance_potential_mv(&chs);
        assert!(chord_conductance_current(&chs, v_rest).abs() < 1e-15, "I(V_rest) = 0");

        // (c) THREAD total_chord_conductance + potential (non-tautological): I = G·(V − V_rest).
        for &vm in &[-100.0_f64, -65.0, 0.0, 30.0] {
            let expected = total_chord_conductance(&chs) * (vm - v_rest);
            assert!(
                (chord_conductance_current(&chs, vm) - expected).abs()
                    <= 1e-9 * expected.abs().max(1e-12),
                "I = G·(V−V_rest) at V={vm}"
            );
        }

        // (d) SIGN + LINEARITY: outward above the rest, inward below, linear in V − V_rest.
        assert!(chord_conductance_current(&chs, v_rest + 10.0) > 0.0, "above rest → outward");
        assert!(chord_conductance_current(&chs, v_rest - 10.0) < 0.0, "below rest → inward");
        assert!(
            (chord_conductance_current(&chs, v_rest + 20.0)
                - 2.0 * chord_conductance_current(&chs, v_rest + 10.0))
            .abs()
                <= 1e-9 * chord_conductance_current(&chs, v_rest + 20.0).abs(),
            "linear in V − V_rest"
        );

        // (e) Non-positive conductances are ignored (matching the family convention).
        let with_dead = [
            ConductanceChannel::new(1.0e-9, -90.0),
            ConductanceChannel::new(3.0e-9, 60.0),
            ConductanceChannel::new(0.0, 1000.0),
            ConductanceChannel::new(-2.0e-9, -1000.0),
        ];
        assert!(
            (chord_conductance_current(&with_dead, 0.0) - chord_conductance_current(&chs, 0.0)).abs()
                <= 1e-9 * 90.0e-9,
            "dead channels ignored"
        );
    }

    #[test]
    fn chord_conductance_is_the_conductance_weighted_reversal_mean() {
        // A purely K⁺ membrane rests exactly at E_K.
        let k_only = chord_conductance_potential_mv(&[ConductanceChannel::new(3.0, -90.0)]);
        assert!((k_only + 90.0).abs() < 1e-12, "K-only {k_only}");
        // Equal Na⁺/K⁺ conductance → the midpoint of the two reversals.
        let mid = chord_conductance_potential_mv(&[
            ConductanceChannel::new(1.0, 60.0),
            ConductanceChannel::new(1.0, -90.0),
        ]);
        assert!((mid + 15.0).abs() < 1e-12, "midpoint {mid}");
        // K⁺-dominated (g_K = 10·g_Na) → a realistic ~ −76 mV resting potential.
        let rest = chord_conductance_potential_mv(&[
            ConductanceChannel::new(1.0, 60.0),
            ConductanceChannel::new(10.0, -90.0),
        ]);
        assert!((rest - (-840.0 / 11.0)).abs() < 1e-12, "rest {rest}");
        assert!((-80.0..=-70.0).contains(&rest), "physiological rest {rest}");
        // Non-positive conductances are ignored; an empty / zero set → 0.
        assert_eq!(
            chord_conductance_potential_mv(&[
                ConductanceChannel::new(0.0, 60.0),
                ConductanceChannel::new(-1.0, 10.0),
            ]),
            0.0
        );
        assert_eq!(chord_conductance_potential_mv(&[]), 0.0);
    }
}

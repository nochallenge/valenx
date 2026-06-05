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

#[cfg(test)]
mod tests {
    use super::*;

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

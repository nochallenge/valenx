//! Goldman–Hodgkin–Katz (GHK) voltage equation — the steady resting membrane
//! potential set by *several* permeant ions at once.
//!
//! Where the Nernst potential gives the equilibrium of a single ion, a real
//! membrane sits at the weighted compromise of every ion that can cross it,
//! each weighted by its permeability `P`:
//!
//! ```text
//!         R·T      P_K[K]_o + P_Na[Na]_o + … + P_Cl[Cl]_i + …
//! V_m  =  ───· ln ───────────────────────────────────────────
//!          F       P_K[K]_i + P_Na[Na]_i + … + P_Cl[Cl]_o + …
//! ```
//!
//! Cations contribute their *outside* concentration to the numerator and
//! *inside* to the denominator; anions (e.g. Cl⁻) flip, reflecting their
//! opposite charge. Only monovalent ions are handled — the classic GHK form.
//! Permeabilities are relative (only ratios matter) and concentrations may be
//! in any shared unit. With a single permeant ion the equation collapses back
//! to that ion's Nernst potential.

use crate::nernst::{thermal_voltage_mv, FARADAY};

/// A monovalent permeant ion for the [`ghk_potential_mv`] equation.
///
/// `permeability` is the relative membrane permeability `P` (only ratios
/// matter); `outside` and `inside` are the extra- and intracellular
/// concentrations in any shared unit.
#[derive(Clone, Copy, Debug)]
pub struct GhkIon {
    /// Relative membrane permeability `P` (arbitrary units; only ratios matter).
    pub permeability: f64,
    /// External (extracellular) concentration.
    pub outside: f64,
    /// Internal (intracellular) concentration.
    pub inside: f64,
}

impl GhkIon {
    /// Construct a permeant ion from its permeability and its outside / inside
    /// concentrations.
    pub fn new(permeability: f64, outside: f64, inside: f64) -> GhkIon {
        GhkIon {
            permeability,
            outside,
            inside,
        }
    }
}

/// The GHK resting membrane potential in **millivolts** at absolute
/// temperature `temp_k` (K), from the permeant monovalent `cations` and
/// `anions`.
///
/// Each cation adds `P·[outside]` to the numerator and `P·[inside]` to the
/// denominator; each anion flips (`P·[inside]` up top, `P·[outside]` below).
/// At least one ion with a positive permeability and concentration must be
/// supplied — an empty mixture yields a non-finite result the caller should
/// guard.
pub fn ghk_potential_mv(temp_k: f64, cations: &[GhkIon], anions: &[GhkIon]) -> f64 {
    let numerator: f64 = cations.iter().map(|c| c.permeability * c.outside).sum::<f64>()
        + anions.iter().map(|a| a.permeability * a.inside).sum::<f64>();
    let denominator: f64 = cations.iter().map(|c| c.permeability * c.inside).sum::<f64>()
        + anions.iter().map(|a| a.permeability * a.outside).sum::<f64>();
    thermal_voltage_mv(temp_k) * (numerator / denominator).ln()
}

/// The Goldman–Hodgkin–Katz **current** (constant-field flux) equation for a
/// single permeant ion — the nonlinear, *rectifying* current–voltage relation
/// the constant-field assumption (Goldman 1943) predicts, in contrast to the
/// ohmic [`crate::ionic_current`] `g·(V−E)`:
///
/// ```text
///                        c_in − c_out·e^(−u)              z·F·V_m
///   I = P · z · F · u · ─────────────────────  ,    u =  ─────────
///                            1 − e^(−u)                     R·T
/// ```
///
/// `permeability` `P` and the concentrations set the absolute scale (supply
/// consistent units — e.g. `P` in m·s⁻¹ with concentrations in mol·m⁻³ gives a
/// current density in A·m⁻²); the *shape* — the reversal and the rectification —
/// is independent of that choice. Two signatures define it:
///
/// * the current is **exactly zero at the ion's Nernst reversal potential**
///   `E = (R·T/zF)·ln(c_out/c_in)` and reverses sign across it; and
/// * with `c_in ≠ c_out` the `I–V` curve **rectifies** — its limiting slope is
///   proportional to the concentration on the side the current flows *from*
///   (outward ∝ `c_in`, inward ∝ `c_out`) — the behaviour an ohmic channel
///   cannot reproduce.
///
/// The removable `0/0` at `V_m = 0` is resolved with its analytic limit
/// `I → P·z·F·(c_in − c_out)`. `valence` must be non-zero and the concentrations
/// positive; the caller guards degenerate inputs.
pub fn ghk_current_density(
    temp_k: f64,
    permeability: f64,
    valence: f64,
    v_m_mv: f64,
    c_out: f64,
    c_in: f64,
) -> f64 {
    let u = valence * v_m_mv / thermal_voltage_mv(temp_k);
    // u·(…)/(1 − e^(−u)) has a removable singularity at u = 0 where
    // u/(1 − e^(−u)) → 1, so the bracket → (c_in − c_out): the constant-field
    // current passes finitely through V_m = 0 rather than evaluating 0/0.
    let bracket = if u.abs() < 1.0e-9 {
        c_in - c_out
    } else {
        u * (c_in - c_out * (-u).exp()) / (1.0 - (-u).exp())
    };
    permeability * valence * FARADAY * bracket
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nernst::{nernst_potential_mv, BODY_TEMPERATURE_K, FARADAY};

    #[test]
    fn ghk_reduces_to_nernst_for_a_single_cation() {
        // One permeant cation ⇒ the GHK potential is exactly its Nernst
        // potential, independent of the (lone) permeability value.
        let k = GhkIon::new(2.7, 5.0, 140.0);
        let ghk = ghk_potential_mv(BODY_TEMPERATURE_K, &[k], &[]);
        let nernst = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 5.0, 140.0);
        assert!((ghk - nernst).abs() < 1e-9, "GHK {ghk} vs Nernst {nernst}");
    }

    #[test]
    fn ghk_reduces_to_nernst_for_a_single_anion() {
        // One permeant anion ⇒ its Nernst potential with valence −1 (the
        // numerator/denominator flip exactly mirrors the z = −1 sign).
        let cl = GhkIon::new(1.0, 110.0, 10.0);
        let ghk = ghk_potential_mv(BODY_TEMPERATURE_K, &[], &[cl]);
        let nernst = nernst_potential_mv(BODY_TEMPERATURE_K, -1.0, 110.0, 10.0);
        assert!((ghk - nernst).abs() < 1e-9, "GHK {ghk} vs Nernst {nernst}");
    }

    #[test]
    fn resting_potential_is_physiological_and_near_e_k() {
        // A mammalian-ish mix: K dominates permeability, Na leaks a little,
        // Cl is moderately permeant. P_K : P_Na : P_Cl = 1 : 0.05 : 0.45.
        let k = GhkIon::new(1.0, 5.0, 140.0);
        let na = GhkIon::new(0.05, 145.0, 15.0);
        let cl = GhkIon::new(0.45, 110.0, 10.0);
        let vm = ghk_potential_mv(BODY_TEMPERATURE_K, &[k, na], &[cl]);
        // A resting cell sits well on the hyperpolarized side.
        assert!((-85.0..=-45.0).contains(&vm), "resting Vm {vm} mV");
        // It is bracketed by the K and Na equilibrium potentials, and sits
        // much closer to E_K because potassium dominates the permeability.
        let e_k = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 5.0, 140.0);
        let e_na = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 145.0, 15.0);
        assert!(e_k < vm && vm < e_na, "Vm {vm} must lie between E_K {e_k} and E_Na {e_na}");
        assert!((vm - e_k).abs() < (vm - e_na).abs(), "resting Vm should hug E_K");
    }

    #[test]
    fn raising_sodium_permeability_depolarises() {
        // The action-potential upstroke: when Na permeability shoots up, the
        // membrane swings from rest toward E_Na (depolarizes).
        let k = GhkIon::new(1.0, 5.0, 140.0);
        let cl = GhkIon::new(0.45, 110.0, 10.0);
        let rest = ghk_potential_mv(
            BODY_TEMPERATURE_K,
            &[k, GhkIon::new(0.05, 145.0, 15.0)],
            &[cl],
        );
        let firing = ghk_potential_mv(
            BODY_TEMPERATURE_K,
            &[k, GhkIon::new(5.0, 145.0, 15.0)],
            &[cl],
        );
        assert!(firing > rest, "more Na permeability depolarizes: {rest} -> {firing}");
        assert!(firing > 0.0, "a Na-dominated membrane overshoots positive: {firing}");
    }

    #[test]
    fn ghk_current_vanishes_at_the_nernst_reversal() {
        // The defining property: the constant-field current crosses zero exactly
        // at the ion's equilibrium potential, and is large just a few mV away.
        let t = BODY_TEMPERATURE_K;
        let (c_out, c_in, z) = (4.0, 140.0, 1.0); // a K⁺ gradient
        let e_rev = nernst_potential_mv(t, z, c_out, c_in);
        let i_rev = ghk_current_density(t, 1.0, z, e_rev, c_out, c_in);
        let i_near = ghk_current_density(t, 1.0, z, e_rev + 10.0, c_out, c_in);
        assert!(i_rev.abs() < 1e-6, "current must vanish at reversal, got {i_rev}");
        assert!(i_near.abs() > 1.0, "and be clearly non-zero 10 mV away, got {i_near}");
    }

    #[test]
    fn ghk_current_is_finite_at_zero_voltage() {
        // The removable 0/0 at V_m = 0 resolves to the analytic limit P·z·F·(c_in − c_out).
        let t = BODY_TEMPERATURE_K;
        let i0 = ghk_current_density(t, 1.0, 1.0, 0.0, 4.0, 140.0);
        assert!(i0.is_finite(), "no 0/0 blow-up at V_m = 0");
        let expected = FARADAY * (140.0 - 4.0); // P = z = 1
        assert!((i0 - expected).abs() < 1e-6 * expected, "V=0 limit {i0} vs {expected}");
        // More K⁺ inside than out ⇒ a net outward (positive) diffusive current.
        assert!(i0 > 0.0, "outward diffusion at rest gives a positive current");
    }

    #[test]
    fn ghk_current_rectifies_with_an_asymmetric_gradient() {
        // The hallmark that separates GHK from an ohmic channel: with c_in ≠ c_out
        // the I–V curve is asymmetric. Its limiting (large-|V|) slope is ∝ the
        // source-side concentration — outward ∝ c_in, inward ∝ c_out — so the
        // slope ratio recovers c_in/c_out.
        let t = BODY_TEMPERATURE_K;
        let (c_out, c_in) = (10.0, 100.0);
        let slope_out = (ghk_current_density(t, 1.0, 1.0, 300.0, c_out, c_in)
            - ghk_current_density(t, 1.0, 1.0, 250.0, c_out, c_in))
            / 50.0;
        let slope_in = (ghk_current_density(t, 1.0, 1.0, -250.0, c_out, c_in)
            - ghk_current_density(t, 1.0, 1.0, -300.0, c_out, c_in))
            / 50.0;
        let ratio = slope_out / slope_in;
        let expected = c_in / c_out;
        assert!(
            (ratio - expected).abs() < 0.05 * expected,
            "rectification slope ratio {ratio} vs c_in/c_out {expected}"
        );
        // An ohmic channel would be symmetric (ratio 1); GHK clearly is not.
        assert!(ratio > 5.0, "GHK strongly rectifies here, ratio {ratio}");
    }
}

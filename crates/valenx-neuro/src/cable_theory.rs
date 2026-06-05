//! Passive-cable constants (Rall cable theory).
//!
//! Before an axon spikes, it behaves as a leaky RC cable. Two constants set
//! how a subthreshold signal spreads and settles:
//!
//! - **Length (space) constant** `λ = sqrt(d·R_m / (4·R_i))` — the distance
//!   over which a steady voltage decays to `1/e`. Larger fibres, higher
//!   membrane resistance, and lower axoplasmic resistivity all let signals
//!   reach farther.
//! - **Membrane time constant** `τ = R_m·C_m` — how fast the membrane charges.
//!
//! and the dimensionless **electrotonic length** `L = ℓ/λ` rates a real
//! segment's length against `λ` (a compartment with `L ≪ 1` is electrically
//! compact; `L ≫ 1` is "long").
//!
//! Units are the cellular-physiology CGS-ish convention: `d` in cm, specific
//! membrane resistance `R_m` in Ω·cm², axial resistivity `R_i` in Ω·cm,
//! specific membrane capacitance `C_m` in F/cm² — giving `λ` in cm and `τ` in
//! seconds. These are the exact passive results, not a fit.

/// Length (space) constant `λ = sqrt(d·R_m / (4·R_i))` in cm, from fibre
/// diameter `d_cm` (cm), specific membrane resistance `r_m_ohm_cm2` (Ω·cm²),
/// and axial resistivity `r_i_ohm_cm` (Ω·cm).
pub fn length_constant_cm(d_cm: f64, r_m_ohm_cm2: f64, r_i_ohm_cm: f64) -> f64 {
    (d_cm * r_m_ohm_cm2 / (4.0 * r_i_ohm_cm)).sqrt()
}

/// Membrane time constant `τ = R_m·C_m` in seconds, from specific membrane
/// resistance `r_m_ohm_cm2` (Ω·cm²) and capacitance `c_m_f_cm2` (F/cm²).
pub fn time_constant_s(r_m_ohm_cm2: f64, c_m_f_cm2: f64) -> f64 {
    r_m_ohm_cm2 * c_m_f_cm2
}

/// Dimensionless electrotonic length `L = ℓ / λ` for a physical length
/// `length_cm` and a length constant `lambda_cm` (both cm). `lambda_cm` must
/// be positive.
pub fn electrotonic_length(length_cm: f64, lambda_cm: f64) -> f64 {
    length_cm / lambda_cm
}

/// Semi-infinite cable input resistance `R_∞ = √(r_m·r_a)` in ohms — the
/// resistance seen looking into one end of a long cylindrical fibre, from fibre
/// diameter `d_cm` (cm), specific membrane resistance `r_m_ohm_cm2` (Ω·cm²) and
/// axial resistivity `r_i_ohm_cm` (Ω·cm). Equivalently `R_∞ = r_a·λ`, the axial
/// resistance per unit length `r_a = 4·R_i/(π·d²)` times the length constant λ;
/// it falls as `d^(−3/2)`, so a thicker fibre presents a much lower input load.
pub fn semi_infinite_input_resistance(d_cm: f64, r_m_ohm_cm2: f64, r_i_ohm_cm: f64) -> f64 {
    // Per-unit-length membrane (Ω·cm) and axial (Ω/cm) resistances.
    let r_m = r_m_ohm_cm2 / (std::f64::consts::PI * d_cm);
    let r_a = 4.0 * r_i_ohm_cm / (std::f64::consts::PI * d_cm * d_cm);
    (r_m * r_a).sqrt()
}

/// The fraction of the steady-state voltage a passive membrane has reached at
/// time `t_s` (s) after a step of current — the RC charging curve
/// `1 − e^(−t/τ)`, with `tau_s` the membrane time constant (s). It reaches
/// ~63% at one time constant and ~95% at three. Returns `0` for a non-positive
/// time constant.
pub fn charging_fraction(t_s: f64, tau_s: f64) -> f64 {
    if tau_s <= 0.0 {
        return 0.0;
    }
    1.0 - (-t_s / tau_s).exp()
}

/// The time (s) for a passive membrane to charge to `fraction` of its
/// steady-state voltage, `t = −τ·ln(1 − fraction)` — the inverse of
/// [`charging_fraction`]. `None` for a `fraction` outside `[0, 1)` (the membrane
/// approaches but never reaches 100%) or a non-positive time constant.
pub fn time_to_charge_fraction(fraction: f64, tau_s: f64) -> Option<f64> {
    if (0.0..1.0).contains(&fraction) && tau_s > 0.0 {
        Some(-tau_s * (1.0 - fraction).ln())
    } else {
        None
    }
}

/// Rall's 3/2-power rule — the **equivalent-cylinder** diameter for a set of
/// daughter branches meeting at a node. Matching the cable input conductance at
/// the branch point requires `d_eq^(3/2) = Σ dᵢ^(3/2)`, hence
/// `d_eq = (Σ dᵢ^(3/2))^(2/3)` (all diameters in the same units). When every
/// branch obeys this rule the whole dendritic tree collapses to one uniform
/// cylinder, so the passive [`length_constant_cm`] / [`time_constant_s`]
/// analysis applies to the tree as a single cable. Non-positive diameters are
/// ignored; an empty branch list gives `0`.
pub fn rall_equivalent_diameter(child_diameters: &[f64]) -> f64 {
    child_diameters
        .iter()
        .filter(|&&d| d > 0.0)
        .map(|&d| d.powf(1.5))
        .sum::<f64>()
        .powf(2.0 / 3.0)
}

/// Steady-state voltage attenuation `V/V₀ = e^(−L)` along a semi-infinite
/// passive cable, as a function of the dimensionless electrotonic distance
/// `L = x/λ` (see [`electrotonic_length`]). This is the spatial counterpart to
/// the temporal [`charging_fraction`]: a steady subthreshold voltage decays to
/// `1/e` (~37%) over one length constant and to ~5% by `L = 3`. Returns `1` at
/// the injection site (`L = 0`); a negative `L` clamps to that injection value.
pub fn steady_state_attenuation(electrotonic_length: f64) -> f64 {
    (-electrotonic_length.max(0.0)).exp()
}

/// Input resistance of a finite passive cable with a **sealed** (closed,
/// no-axial-leak) far end: `R_in = R_∞·coth(L)`, from the semi-infinite input
/// resistance `r_inf_ohms` ([`semi_infinite_input_resistance`]) and the
/// electrotonic length `electrotonic_length` `L` ([`electrotonic_length`]). A
/// sealed terminal traps current, so a short cable presents a *higher*
/// resistance than the semi-infinite case; as `L → ∞` it relaxes to `R_∞`.
/// Returns `R_∞` for non-positive `L`.
pub fn sealed_end_input_resistance(r_inf_ohms: f64, electrotonic_length: f64) -> f64 {
    if electrotonic_length > 0.0 {
        // coth(L) = 1 / tanh(L).
        r_inf_ohms / electrotonic_length.tanh()
    } else {
        r_inf_ohms
    }
}

/// Input resistance of a finite passive cable with an **open** (killed,
/// clamped-to-rest) far end: `R_in = R_∞·tanh(L)`. An open terminal sinks
/// current freely, so a short cable presents a *lower* resistance than the
/// semi-infinite case; as `L → ∞` it relaxes to `R_∞`. Returns `0` for
/// non-positive `L` (a zero-length open cable shorts straight to rest).
pub fn open_end_input_resistance(r_inf_ohms: f64, electrotonic_length: f64) -> f64 {
    if electrotonic_length > 0.0 {
        r_inf_ohms * electrotonic_length.tanh()
    } else {
        0.0
    }
}

/// The membrane low-pass cutoff (corner) frequency `f_c = 1/(2π·τ)` in hertz,
/// from the membrane time constant `tau_s` (s). The passive RC membrane behaves
/// as a first-order low-pass filter: voltage signals slower than `f_c` pass
/// through, faster ones are attenuated. With `τ ≈ 10 ms` a neuron passes only
/// ~16 Hz. The frequency-domain companion of [`time_constant_s`]. Returns `0`
/// for a non-positive time constant.
pub fn membrane_cutoff_frequency(tau_s: f64) -> f64 {
    if tau_s > 0.0 {
        1.0 / (2.0 * std::f64::consts::PI * tau_s)
    } else {
        0.0
    }
}

/// The frequency-dependent (**AC**) length constant
/// `|λ_AC(f)| = λ / (1 + (2π·f·τ)²)^(1/4)` (same units as `dc_length_constant`) —
/// how far a *sinusoidal* signal at `frequency_hz` spreads passively along a
/// cable, given the steady (DC) length constant `dc_length_constant` (`λ`) and
/// the membrane time constant `time_constant_s` (`τ`, s). Because the membrane is
/// a low-pass RC filter, faster signals are shunted by the membrane capacitance
/// and decay over a *shorter* distance than slow ones: `|λ_AC|` equals `λ` at DC
/// (`f = 0`), falls monotonically, and at the membrane cutoff `f_c = 1/(2πτ)`
/// ([`membrane_cutoff_frequency`]) has dropped to `λ/2^(1/4) ≈ 0.84·λ`. This is
/// the spatial counterpart of the temporal low-pass filtering — dendrites blur
/// fast synaptic input in space as well as time. Returns `0` for non-physical
/// input (`λ` or `τ` non-positive, `f` negative, or non-finite).
pub fn ac_length_constant(dc_length_constant: f64, time_constant_s: f64, frequency_hz: f64) -> f64 {
    if !dc_length_constant.is_finite()
        || dc_length_constant <= 0.0
        || !time_constant_s.is_finite()
        || time_constant_s <= 0.0
        || !frequency_hz.is_finite()
        || frequency_hz < 0.0
    {
        return 0.0;
    }
    let omega_tau = 2.0 * std::f64::consts::PI * frequency_hz * time_constant_s;
    dc_length_constant / (1.0 + omega_tau * omega_tau).powf(0.25)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_constant_matches_the_formula_and_scales() {
        // λ = sqrt(d·R_m/(4·R_i)); d=4e-4, R_m=1e4, R_i=100 → sqrt(0.01) = 0.1 cm.
        assert!((length_constant_cm(4e-4, 1e4, 100.0) - 0.1).abs() < 1e-12);
        // λ ∝ √d: quadrupling the diameter doubles λ.
        let l1 = length_constant_cm(4e-4, 1e4, 100.0);
        let l4 = length_constant_cm(16e-4, 1e4, 100.0);
        assert!((l4 - 2.0 * l1).abs() < 1e-12);
        // λ ∝ 1/√(R_i): quadrupling axial resistivity halves λ.
        let lr = length_constant_cm(4e-4, 1e4, 400.0);
        assert!((lr - 0.5 * l1).abs() < 1e-12);
    }

    #[test]
    fn time_constant_is_rm_times_cm() {
        // τ = R_m·C_m; 1e4 Ω·cm² × 1e-6 F/cm² = 0.01 s (10 ms).
        assert!((time_constant_s(1e4, 1e-6) - 0.01).abs() < 1e-15);
        // Linear in each factor.
        assert!((time_constant_s(2e4, 1e-6) - 2.0 * time_constant_s(1e4, 1e-6)).abs() < 1e-15);
    }

    #[test]
    fn electrotonic_length_is_length_over_lambda() {
        assert!((electrotonic_length(0.2, 0.1) - 2.0).abs() < 1e-12);
        // A compartment much shorter than λ is electrically compact (L ≪ 1).
        assert!(electrotonic_length(0.01, 0.1) < 0.2);
    }

    #[test]
    fn semi_infinite_input_resistance_relates_to_length_constant() {
        // d=4e-4 cm (4 µm), R_m=1e4 Ω·cm², R_i=100 Ω·cm.
        let r_inf = semi_infinite_input_resistance(4e-4, 1e4, 100.0);
        assert!(r_inf > 0.0 && r_inf.is_finite());
        // R_∞ = r_a·λ — cross-check against the length constant.
        let lambda = length_constant_cm(4e-4, 1e4, 100.0);
        let r_a = 4.0 * 100.0 / (std::f64::consts::PI * 4e-4 * 4e-4);
        assert!((r_inf - r_a * lambda).abs() / r_inf < 1e-9, "R_∞ = r_a·λ: {r_inf}");
        // ~80 MΩ for this thin fibre.
        assert!((r_inf / 1e6 - 79.6).abs() < 1.0, "R_∞ {} MΩ", r_inf / 1e6);
        // R_∞ ∝ d^(−3/2): a 4× thicker fibre has 1/8 the input resistance.
        let thick = semi_infinite_input_resistance(16e-4, 1e4, 100.0);
        assert!((thick - r_inf / 8.0).abs() / r_inf < 1e-9, "d^(−3/2) scaling");
    }

    #[test]
    fn charging_fraction_follows_the_rc_curve() {
        // 0 at t=0, 1 − 1/e (~63.2%) at one τ, ~95.0% at three τ.
        assert!(charging_fraction(0.0, 0.01).abs() < 1e-12);
        assert!(
            (charging_fraction(0.01, 0.01) - (1.0 - 1.0 / std::f64::consts::E)).abs() < 1e-12
        );
        assert!((charging_fraction(0.03, 0.01) - 0.9502).abs() < 1e-3);
        // Monotonically rising toward 1; a non-positive τ gives 0.
        assert!(charging_fraction(0.05, 0.01) > charging_fraction(0.03, 0.01));
        assert_eq!(charging_fraction(0.01, 0.0), 0.0);
    }

    #[test]
    fn time_to_charge_fraction_inverts_the_curve() {
        let tau = 0.01;
        // Reaching 1 − 1/e takes exactly one τ.
        let t = time_to_charge_fraction(1.0 - 1.0 / std::f64::consts::E, tau).unwrap();
        assert!((t - tau).abs() < 1e-12, "one τ to 63%, got {t}");
        // Round-trip: charge to f over time t, and the fraction at t is f again.
        let f = 0.8;
        let tf = time_to_charge_fraction(f, tau).unwrap();
        assert!((charging_fraction(tf, tau) - f).abs() < 1e-12);
        // 100% (or out of range) is never reached → None.
        assert!(time_to_charge_fraction(1.0, tau).is_none());
        assert!(time_to_charge_fraction(-0.1, tau).is_none());
    }

    #[test]
    fn rall_equivalent_diameter_follows_the_three_halves_rule() {
        // A lone branch is its own equivalent cylinder.
        assert!((rall_equivalent_diameter(&[2.0]) - 2.0).abs() < 1e-12);
        // Two equal d-branches → 2^(2/3)·d (≈ 1.587 d), NOT 2 d — the cross-
        // sectional areas do not simply add.
        let d = 1.5;
        let eq = rall_equivalent_diameter(&[d, d]);
        assert!(
            (eq - 2.0_f64.powf(2.0 / 3.0) * d).abs() < 1e-12,
            "two equal branches: {eq}"
        );
        // The defining identity holds for an uneven tree: d_eq^(3/2) = Σ dᵢ^(3/2).
        let kids = [1.0, 0.8, 0.6];
        let parent = rall_equivalent_diameter(&kids);
        let rhs: f64 = kids.iter().map(|d| d.powf(1.5)).sum();
        assert!((parent.powf(1.5) - rhs).abs() < 1e-12, "3/2 identity");
        // Negative/zero diameters are ignored; an empty list gives 0.
        assert!((rall_equivalent_diameter(&[d, -1.0, 0.0]) - d).abs() < 1e-12);
        assert_eq!(rall_equivalent_diameter(&[]), 0.0);
    }

    #[test]
    fn steady_state_attenuation_decays_exponentially_over_lambda() {
        // V/V₀ = 1 at the injection site, 1/e at one λ, ~5% by L = 3.
        assert!((steady_state_attenuation(0.0) - 1.0).abs() < 1e-12);
        assert!(
            (steady_state_attenuation(1.0) - 1.0 / std::f64::consts::E).abs() < 1e-12
        );
        assert!((steady_state_attenuation(3.0) - 0.049_787).abs() < 1e-5);
        // Composes with electrotonic_length: a 0.2 cm segment on a λ = 0.1 cm
        // cable (L = 2) attenuates to e^(−2).
        let l = electrotonic_length(0.2, 0.1);
        assert!((steady_state_attenuation(l) - (-2.0_f64).exp()).abs() < 1e-12);
        // Monotonically decreasing; a negative L clamps to the injection value.
        assert!(steady_state_attenuation(2.0) < steady_state_attenuation(1.0));
        assert_eq!(steady_state_attenuation(-1.0), 1.0);
    }

    #[test]
    fn finite_cable_input_resistances_bracket_the_semi_infinite_case() {
        let r_inf = 80.0e6;
        // As L → ∞ both boundary conditions relax to the semi-infinite R_∞.
        assert!((sealed_end_input_resistance(r_inf, 10.0) - r_inf).abs() / r_inf < 1e-6);
        assert!((open_end_input_resistance(r_inf, 10.0) - r_inf).abs() / r_inf < 1e-6);
        // A sealed short cable traps current → higher than R_∞; open sinks it → lower.
        assert!(sealed_end_input_resistance(r_inf, 0.5) > r_inf);
        assert!(open_end_input_resistance(r_inf, 0.5) < r_inf);
        // Known values at L = 1: coth(1) ≈ 1.31304, tanh(1) ≈ 0.76159.
        assert!((sealed_end_input_resistance(r_inf, 1.0) / r_inf - 1.31304).abs() < 1e-4);
        assert!((open_end_input_resistance(r_inf, 1.0) / r_inf - 0.76159).abs() < 1e-4);
        // Sealed always exceeds open at the same length (trapped vs sunk current).
        assert!(
            sealed_end_input_resistance(r_inf, 0.8) > open_end_input_resistance(r_inf, 0.8)
        );
        // Degenerate length: sealed → R_∞, open → 0 (a direct short to rest).
        assert_eq!(sealed_end_input_resistance(r_inf, 0.0), r_inf);
        assert_eq!(open_end_input_resistance(r_inf, 0.0), 0.0);
    }

    #[test]
    fn membrane_cutoff_frequency_is_the_rc_corner() {
        // f_c = 1/(2π·τ); τ = 10 ms → ≈ 15.915 Hz.
        assert!((membrane_cutoff_frequency(0.01) - 15.915).abs() < 1e-2);
        // Inversely proportional to τ: a 2× slower membrane halves the cutoff.
        assert!(
            (membrane_cutoff_frequency(0.02) - 0.5 * membrane_cutoff_frequency(0.01)).abs() < 1e-9
        );
        // Round-trip: f_c · (2π·τ) = 1.
        let tau = 0.005;
        assert!(
            (membrane_cutoff_frequency(tau) * 2.0 * std::f64::consts::PI * tau - 1.0).abs() < 1e-12
        );
        // Non-positive τ → 0.
        assert_eq!(membrane_cutoff_frequency(0.0), 0.0);
    }

    #[test]
    fn ac_length_constant_low_pass_filters_in_space() {
        let lambda = 0.1; // cm (DC length constant)
        let tau = 0.01; // s (10 ms)
        // At DC the AC length constant equals the steady length constant.
        assert!((ac_length_constant(lambda, tau, 0.0) - lambda).abs() < 1e-12);
        // It matches the closed form at an arbitrary frequency.
        let f = 25.0;
        let wt = 2.0 * std::f64::consts::PI * f * tau;
        let expected = lambda / (1.0 + wt * wt).powf(0.25);
        assert!((ac_length_constant(lambda, tau, f) - expected).abs() < 1e-12);
        // At the membrane cutoff frequency it has dropped to λ/2^(1/4) ≈ 0.841·λ.
        let fc = membrane_cutoff_frequency(tau);
        let at_cutoff = ac_length_constant(lambda, tau, fc);
        assert!(
            (at_cutoff - lambda / 2.0_f64.powf(0.25)).abs() < 1e-9,
            "at cutoff {at_cutoff}"
        );
        // Strictly decreasing in frequency, never exceeds λ, and → 0 at high f.
        assert!(ac_length_constant(lambda, tau, 10.0) > ac_length_constant(lambda, tau, 100.0));
        assert!(ac_length_constant(lambda, tau, 1.0) <= lambda);
        assert!(ac_length_constant(lambda, tau, 1.0e6) < 0.01 * lambda);
        // Non-physical inputs → 0.
        assert_eq!(ac_length_constant(0.0, tau, f), 0.0);
        assert_eq!(ac_length_constant(lambda, 0.0, f), 0.0);
        assert_eq!(ac_length_constant(lambda, tau, -1.0), 0.0);
        assert_eq!(ac_length_constant(f64::NAN, tau, f), 0.0);
    }
}

//! Electrode–tissue impedance: what the stimulating electrode "sees."
//!
//! A lumped model of a disk microelectrode in tissue: the **access
//! (spreading) resistance** `R_a = 1/(4σa)` in series with a **constant-phase
//! element (CPE)** modelling the electrode–electrolyte double layer. The
//! impedance is capacitive (large) at low frequency and falls to the resistive
//! `R_a` plateau at high frequency.

/// A constant-phase element: `Z_CPE(ω) = 1 / (Q·(jω)ⁿ)`.
#[derive(Debug, Clone, Copy)]
pub struct Cpe {
    /// Magnitude parameter `Q` (S·sⁿ).
    pub q: f64,
    /// Exponent `n` in [0, 1] (1 = ideal capacitor, 0 = pure resistor).
    pub n: f64,
}

impl Default for Cpe {
    fn default() -> Self {
        // Representative microelectrode double-layer values.
        Self { q: 1.0e-5, n: 0.9 }
    }
}

/// A disk microelectrode of radius `a` in tissue of conductivity `σ`, with a
/// double-layer CPE in series.
#[derive(Debug, Clone, Copy)]
pub struct ElectrodeImpedance {
    a_m: f64,
    sigma_s_m: f64,
    cpe: Cpe,
}

impl ElectrodeImpedance {
    /// A disk electrode of radius `a_um` (µm) in tissue of conductivity
    /// `sigma_s_m` (S/m), with double-layer `cpe`.
    pub fn disk(a_um: f64, sigma_s_m: f64, cpe: Cpe) -> Self {
        Self {
            a_m: a_um * 1.0e-6,
            sigma_s_m,
            cpe,
        }
    }

    /// Access (spreading) resistance `R_a = 1/(4σa)` (Ω) — the high-frequency
    /// resistive plateau.
    pub fn access_resistance_ohm(&self) -> f64 {
        1.0 / (4.0 * self.sigma_s_m * self.a_m)
    }

    /// Magnitude of the total impedance `|Z(ω)| = |R_a + Z_CPE|` (Ω) at
    /// frequency `freq_hz`.
    pub fn magnitude_ohm(&self, freq_hz: f64) -> f64 {
        let w = 2.0 * std::f64::consts::PI * freq_hz;
        let mag_cpe = 1.0 / (self.cpe.q * w.powf(self.cpe.n));
        let phase = self.cpe.n * std::f64::consts::FRAC_PI_2;
        let re = self.access_resistance_ohm() + mag_cpe * phase.cos();
        let im = -mag_cpe * phase.sin();
        (re * re + im * im).sqrt()
    }

    /// Phase angle of the total impedance `∠Z(ω) = ∠(R_a + Z_CPE)` (degrees) at
    /// frequency `freq_hz` — the second half of the electrode's Bode response and the
    /// companion to [`ElectrodeImpedance::magnitude_ohm`] (which forms the same
    /// `Re`/`Im` and keeps only `|Z|`). The double-layer CPE lags the current, so the
    /// phase is **negative** (capacitive): it tends to the bare-CPE constant `−n·90°`
    /// at low frequency, where the reactive double layer dominates, and rises toward
    /// `0°` (purely resistive) at high frequency, where the access-resistance plateau
    /// `R_a` takes over. Together with the magnitude this is what an
    /// impedance-spectroscopy sweep of the electrode reports.
    pub fn phase_deg(&self, freq_hz: f64) -> f64 {
        let w = 2.0 * std::f64::consts::PI * freq_hz;
        let mag_cpe = 1.0 / (self.cpe.q * w.powf(self.cpe.n));
        let phase = self.cpe.n * std::f64::consts::FRAC_PI_2;
        let re = self.access_resistance_ohm() + mag_cpe * phase.cos();
        let im = -mag_cpe * phase.sin();
        im.atan2(re).to_degrees()
    }

    /// The electrode's **corner (crossover) frequency** `f_c` (Hz) — where the
    /// double-layer CPE magnitude equals the access resistance `R_a`, i.e. the Bode
    /// "knee" separating the low-frequency **capacitive** regime (the CPE dominates,
    /// `|Z| → ∞` as `f → 0`) from the high-frequency **resistive** plateau
    /// (`|Z| → R_a`). Setting `1/(Q·ωⁿ) = R_a` gives `ω_c = (1/(Q·R_a))^(1/n)` and
    /// `f_c = ω_c/2π`. It is the natural single-number summary of the
    /// [`ElectrodeImpedance::magnitude_ohm`] / [`ElectrodeImpedance::phase_deg`]
    /// sweep — below it the impedance is reactive, above it nearly real.
    pub fn crossover_frequency_hz(&self) -> f64 {
        let w_c = (1.0 / (self.cpe.q * self.access_resistance_ohm())).powf(1.0 / self.cpe.n);
        w_c / (2.0 * std::f64::consts::PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disk_access_resistance_matches_formula() {
        // a = 50 µm, σ = 0.2 S/m → R_a = 1/(4σa).
        let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe::default());
        let r_a = z.access_resistance_ohm();
        let expect = 1.0 / (4.0 * 0.2 * 50.0e-6);
        assert!(
            (r_a / expect - 1.0).abs() < 1e-9,
            "R_a={r_a} expect={expect}"
        );
    }

    #[test]
    fn impedance_is_capacitive_low_and_resistive_high() {
        let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe::default());
        let lo = z.magnitude_ohm(1.0); // 1 Hz
        let hi = z.magnitude_ohm(1.0e5); // 100 kHz
        assert!(
            lo > hi,
            "low-f capacitive must exceed high-f resistive: lo={lo} hi={hi}"
        );
        // At high frequency the CPE shorts out and |Z| → R_a.
        let r_a = z.access_resistance_ohm();
        assert!(
            (hi - r_a).abs() / r_a < 0.2,
            "high-f magnitude should approach the R_a plateau; hi={hi} R_a={r_a}"
        );
    }

    #[test]
    fn phase_is_the_bode_companion_of_the_magnitude() {
        use std::f64::consts::{FRAC_PI_2, PI};
        let n = 0.9;
        let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe { q: 1.0e-5, n });

        // High frequency: the CPE shorts out → the impedance is the resistive R_a
        // plateau, so the phase → 0°.
        assert!(
            z.phase_deg(1.0e9).abs() < 1e-3,
            "high-f phase → 0° (resistive)"
        );

        // Low frequency: the double-layer CPE dominates → the phase → the bare-CPE
        // constant phase −n·90° (= −81° for n = 0.9).
        let lo = z.phase_deg(1.0e-9);
        assert!(
            (lo - (-n * 90.0)).abs() < 1e-3,
            "low-f phase → −n·90°, got {lo}"
        );

        // At an arbitrary frequency: recompute Re/Im of Z = R_a + Z_CPE independently
        // and check (i) phase_deg = atan2(Im, Re) and (ii) the (magnitude, phase) pair
        // reconstructs that same complex impedance — tying phase_deg to the existing
        // magnitude_ohm via |Z|·cosφ = Re and |Z|·sinφ = Im.
        let f = 1.0e3;
        let w = 2.0 * PI * f;
        let mag_cpe = 1.0 / (1.0e-5 * w.powf(n));
        let cpe_phase = n * FRAC_PI_2;
        let re = z.access_resistance_ohm() + mag_cpe * cpe_phase.cos();
        let im = -mag_cpe * cpe_phase.sin();
        let phi = z.phase_deg(f);
        assert!(
            (phi - im.atan2(re).to_degrees()).abs() < 1e-9,
            "phase = atan2(Im, Re)"
        );
        let mag = z.magnitude_ohm(f);
        let phi_rad = phi.to_radians();
        assert!((mag * phi_rad.cos() - re).abs() < 1e-9, "|Z|·cosφ = Re Z");
        assert!((mag * phi_rad.sin() - im).abs() < 1e-9, "|Z|·sinφ = Im Z");

        // The double layer lags the current: phase is strictly negative and never
        // steeper than the −n·90° CPE asymptote.
        assert!(
            phi < 0.0 && phi > -n * 90.0,
            "capacitive lag in (−n·90°, 0°): {phi}"
        );
    }

    #[test]
    fn crossover_frequency_is_where_the_cpe_equals_the_access_resistance() {
        use std::f64::consts::{FRAC_PI_2, PI};
        let (q, n) = (1.0e-5, 0.9);
        let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe { q, n });
        let f_c = z.crossover_frequency_hz();
        assert!(f_c > 0.0, "corner frequency is positive, got {f_c}");

        // Defining property: at f_c the bare-CPE magnitude 1/(Q·ωⁿ) equals R_a
        // (threads access_resistance_ohm).
        let r_a = z.access_resistance_ohm();
        let w = 2.0 * PI * f_c;
        let mag_cpe = 1.0 / (q * w.powf(n));
        assert!((mag_cpe - r_a).abs() / r_a < 1e-9, "|Z_CPE(f_c)| = R_a");

        // At the knee the two equal-magnitude legs (R_a + Z_CPE at constant phase
        // −nπ/2) sum to |Z| = R_a·√(2(1+cos(nπ/2))) — threads magnitude_ohm.
        let expected = r_a * (2.0 * (1.0 + (n * FRAC_PI_2).cos())).sqrt();
        assert!(
            (z.magnitude_ohm(f_c) - expected).abs() / expected < 1e-9,
            "|Z(f_c)| from two equal-magnitude legs"
        );

        // Capacitive below the knee, resistive above: |Z| a decade below exceeds
        // |Z| a decade above.
        assert!(
            z.magnitude_ohm(f_c / 10.0) > z.magnitude_ohm(f_c * 10.0),
            "magnitude falls through the corner"
        );
    }
}

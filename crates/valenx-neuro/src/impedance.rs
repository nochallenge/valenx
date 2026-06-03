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
        assert!((r_a / expect - 1.0).abs() < 1e-9, "R_a={r_a} expect={expect}");
    }

    #[test]
    fn impedance_is_capacitive_low_and_resistive_high() {
        let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe::default());
        let lo = z.magnitude_ohm(1.0); // 1 Hz
        let hi = z.magnitude_ohm(1.0e5); // 100 kHz
        assert!(lo > hi, "low-f capacitive must exceed high-f resistive: lo={lo} hi={hi}");
        // At high frequency the CPE shorts out and |Z| → R_a.
        let r_a = z.access_resistance_ohm();
        assert!(
            (hi - r_a).abs() / r_a < 0.2,
            "high-f magnitude should approach the R_a plateau; hi={hi} R_a={r_a}"
        );
    }
}

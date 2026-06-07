//! Ohmic ionic current — the per-ion membrane current that an
//! equilibrium potential (see [`crate::nernst`]) actually drives.
//!
//! A channel population behaves as a conductance in series with a battery set
//! to the ion's reversal (Nernst) potential `E_rev`. The **driving force**
//! `V_m − E_rev` is the voltage across that conductance, and the current is
//! Ohm's law on it:
//!
//! ```text
//! I = g · (V_m − E_rev)
//! ```
//!
//! This is the building block of every conductance-based (Hodgkin–Huxley-style)
//! membrane model: the total membrane current is the sum of one such term per
//! ionic species (plus capacitive current). The current vanishes exactly at the
//! reversal potential and reverses sign as `V_m` crosses it.
//!
//! With `V_m` and `E_rev` in millivolts and `g` in mS/cm², the current is in
//! µA/cm² (`mS·mV = µA`).

/// Electrochemical driving force on an ion in **millivolts** — the membrane
/// potential minus the ion's reversal (Nernst) potential, `V_m − E_rev`. Its
/// sign sets the current direction; it is zero at equilibrium (`V_m = E_rev`).
pub fn driving_force_mv(vm_mv: f64, e_rev_mv: f64) -> f64 {
    vm_mv - e_rev_mv
}

/// Ohmic ionic current density `I = g·(V_m − E_rev)`. With `conductance` in
/// mS/cm² and the potentials in mV, the result is in µA/cm². Returns zero at
/// the reversal potential and reverses sign as `vm_mv` crosses `e_rev_mv`.
pub fn ionic_current(conductance: f64, vm_mv: f64, e_rev_mv: f64) -> f64 {
    conductance * driving_force_mv(vm_mv, e_rev_mv)
}

/// Ohmic **power dissipated** by an ionic current `P = I·(V_m − E_rev) = g·(V_m −
/// E_rev)²` — the Joule heating of the channel conductance, and the electrical part
/// of the metabolic cost of holding current against the ion's battery. With
/// `conductance` in mS/cm² and the potentials in mV (the [`ionic_current`]
/// convention), the result is in **nW/cm²** (`µA·mV = nW`).
///
/// Unlike the signed [`ionic_current`], the dissipated power is **non-negative on
/// both sides of the reversal potential**: below `E_rev` the current flows inward
/// (negative) and above it outward (positive), but the product with the driving
/// force is a perfect square `g·(V_m − E_rev)²` either way. It is zero only at the
/// reversal potential, where no current flows.
pub fn ionic_power_density(conductance: f64, vm_mv: f64, e_rev_mv: f64) -> f64 {
    let driving = driving_force_mv(vm_mv, e_rev_mv);
    ionic_current(conductance, vm_mv, e_rev_mv) * driving
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ionic_power_density_is_the_dissipated_joule_heating() {
        let g = 0.36;
        // Threads ionic_current and driving_force; equals g·(V−E)²; and is
        // non-negative on BOTH sides of the reversal potential (unlike the current).
        for &(vm, e) in &[(-65.0_f64, -90.0_f64), (-100.0, -90.0), (60.0, -90.0), (-80.0, 60.0)] {
            let p = ionic_power_density(g, vm, e);
            assert!(
                (p - ionic_current(g, vm, e) * driving_force_mv(vm, e)).abs() < 1e-12,
                "P = I·(V−E)"
            );
            assert!((p - g * (vm - e) * (vm - e)).abs() < 1e-9, "P = g·(V−E)²");
            assert!(p >= 0.0, "dissipated power is non-negative: {p} at V={vm}, E={e}");
        }
        // Below reversal the current is inward (negative) yet the power is positive.
        assert!(ionic_current(g, -100.0, -90.0) < 0.0 && ionic_power_density(g, -100.0, -90.0) > 0.0);

        // Worked value: g = 0.36 mS/cm², ΔV = 25 mV → P = 0.36·625 = 225 nW/cm².
        assert!((ionic_power_density(0.36, -65.0, -90.0) - 225.0).abs() < 1e-9);
        // No dissipation at the reversal potential (no current flows).
        assert!(ionic_power_density(0.36, -77.0, -77.0).abs() < 1e-12, "zero at E_rev");
        // Quadratic in the driving force: doubling V_m − E quadruples the power.
        let p25 = ionic_power_density(0.36, -65.0, -90.0); // ΔV = 25
        let p50 = ionic_power_density(0.36, -40.0, -90.0); // ΔV = 50
        assert!((p50 - 4.0 * p25).abs() / (4.0 * p25) < 1e-9, "quadratic in driving force");
    }

    #[test]
    fn driving_force_is_vm_minus_reversal() {
        // Zero at the reversal potential — equilibrium, no net drive.
        assert!(driving_force_mv(-90.0, -90.0).abs() < 1e-12);
        // At rest (−65 mV): K⁺ is gently driven outward (+25 mV), Na⁺ strongly
        // inward (−125 mV) — the asymmetry behind the action potential.
        assert!((driving_force_mv(-65.0, -90.0) - 25.0).abs() < 1e-12);
        assert!((driving_force_mv(-65.0, 60.0) - (-125.0)).abs() < 1e-12);
    }

    #[test]
    fn ionic_current_is_ohmic_and_vanishes_at_reversal() {
        // No current flows when V_m sits at the reversal potential.
        assert!(ionic_current(0.36, -77.0, -77.0).abs() < 1e-12);
        // Ohm's law I = g·(V_m − E), linear in the conductance.
        assert!((ionic_current(0.36, -65.0, -90.0) - 0.36 * 25.0).abs() < 1e-9);
        assert!(
            (ionic_current(0.72, -65.0, -90.0) - 2.0 * ionic_current(0.36, -65.0, -90.0)).abs()
                < 1e-9
        );
        // The current reverses sign as V_m crosses the reversal potential.
        let below = ionic_current(0.36, -100.0, -90.0); // V_m < E → inward (−)
        let above = ionic_current(0.36, -80.0, -90.0); // V_m > E → outward (+)
        assert!(below < 0.0 && above > 0.0, "current reverses across E: {below}, {above}");
    }
}

//! Gibbs–Donnan equilibrium — the passive partition of permeant ions across a
//! membrane that traps an **impermeant** fixed charge.
//!
//! Where the Nernst potential is the equilibrium of one ion and the GHK
//! equation is the steady resting potential set by several *permeant* ions, the
//! Donnan equilibrium is a different idea: a membrane freely permeable to a
//! small 1:1 salt (a univalent cation and anion) but **not** to a large charged
//! species trapped on one side — intracellular protein, or the fixed charge of
//! a gel or cartilage. That trapped charge forces an unequal distribution of
//! the permeant ions, an osmotic imbalance, and a standing **Donnan potential**.
//!
//! With an external 1:1 salt bath at concentration `c` and an internal
//! univalent impermeant anion at concentration `A`, electroneutrality inside
//! (`[cat]_i = [an]_i + A`) together with the Donnan condition that the permeant
//! product is equal across the membrane (`[cat]_i·[an]_i = c²`) gives the
//! **Donnan ratio**
//!
//! ```text
//!         [cat]_i     [an]_o      A + √(A² + 4c²)
//!   r  =  ───────  =  ──────  =  ────────────────
//!         [cat]_o     [an]_i           2c
//! ```
//!
//! and a **Donnan potential** `V = −(R·T/F)·ln r`, the Nernst potential of either
//! permeant ion at equilibrium. The trapped anion concentrates cations inside
//! (`r > 1`) and holds the interior electronegative (`V < 0`).

use crate::nernst::thermal_voltage_mv;

/// The **Donnan ratio** `r = (A + √(A² + 4c²)) / (2c)` — the equilibrium ratio of
/// permeant-cation concentration inside to outside (equivalently, permeant-anion
/// outside to inside) for an external 1:1 salt at concentration `salt_conc` (`c`)
/// and an internal univalent impermeant anion at concentration `impermeant_charge`
/// (`A`). Concentrations may be in any shared unit. `r = 1` (no partition) when
/// `A = 0`; `r > 1` and grows with `A`. Returns `1` for non-physical input
/// (`c ≤ 0` or `A < 0`, non-finite), where no partition is defined.
pub fn donnan_ratio(salt_conc: f64, impermeant_charge: f64) -> f64 {
    if !salt_conc.is_finite()
        || salt_conc <= 0.0
        || !impermeant_charge.is_finite()
        || impermeant_charge < 0.0
    {
        return 1.0;
    }
    let a = impermeant_charge;
    let c = salt_conc;
    (a + (a * a + 4.0 * c * c).sqrt()) / (2.0 * c)
}

/// The **Donnan potential** `V = −(R·T/F)·ln r` in **millivolts** at absolute
/// temperature `temp_k` (K) — the standing membrane potential set by the
/// trapped impermeant charge, equal to the Nernst potential of either permeant
/// ion at Donnan equilibrium. `salt_conc` (`c`) and `impermeant_charge` (`A`)
/// are as in [`donnan_ratio`]. It is `0` with no impermeant charge and negative
/// (interior electronegative) once `A > 0`.
pub fn donnan_potential_mv(temp_k: f64, salt_conc: f64, impermeant_charge: f64) -> f64 {
    -thermal_voltage_mv(temp_k) * donnan_ratio(salt_conc, impermeant_charge).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nernst::BODY_TEMPERATURE_K;

    #[test]
    fn no_impermeant_charge_gives_no_partition() {
        // A = 0 → the salt distributes evenly: r = 1 and V = 0.
        let r = donnan_ratio(100.0, 0.0);
        assert!((r - 1.0).abs() < 1e-12, "ratio {r}");
        assert!(donnan_potential_mv(BODY_TEMPERATURE_K, 100.0, 0.0).abs() < 1e-12);
    }

    #[test]
    fn trapped_anion_concentrates_cations_and_makes_the_inside_negative() {
        // c = 100, A = 50: r = (50 + √(2500 + 40000))/200 = (50 + √42500)/200.
        let (c, a) = (100.0_f64, 50.0_f64);
        let r = donnan_ratio(c, a);
        let expected = (a + (a * a + 4.0 * c * c).sqrt()) / (2.0 * c);
        assert!((r - expected).abs() < 1e-12);
        assert!(r > 1.0, "trapped anion concentrates cations inside: r = {r}");

        // The interior is electronegative (V < 0) — the Donnan potential.
        let v = donnan_potential_mv(BODY_TEMPERATURE_K, c, a);
        assert!(v < 0.0, "inside should be electronegative: V = {v} mV");
        // It equals the Nernst potential of the permeant cation at equilibrium:
        // V = −V_T·ln r.
        assert!((v + thermal_voltage_mv(BODY_TEMPERATURE_K) * r.ln()).abs() < 1e-9);

        // The Donnan condition holds: with [cat]_i = r·c and [an]_i = c/r, the
        // permeant product equals c² and electroneutrality [cat]_i = [an]_i + A.
        let cat_i = r * c;
        let an_i = c / r;
        assert!((cat_i * an_i - c * c).abs() < 1e-6, "Donnan product = c²");
        assert!((cat_i - an_i - a).abs() < 1e-6, "electroneutrality [cat]_i = [an]_i + A");
    }

    #[test]
    fn ratio_grows_monotonically_with_the_trapped_charge() {
        let c = 120.0;
        let r0 = donnan_ratio(c, 0.0);
        let r1 = donnan_ratio(c, 40.0);
        let r2 = donnan_ratio(c, 120.0);
        assert!(r0 < r1 && r1 < r2, "r should grow with A: {r0} < {r1} < {r2}");
        assert!((r0 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn non_physical_inputs_return_the_neutral_ratio() {
        assert_eq!(donnan_ratio(0.0, 50.0), 1.0);
        assert_eq!(donnan_ratio(-1.0, 50.0), 1.0);
        assert_eq!(donnan_ratio(100.0, -1.0), 1.0);
        assert_eq!(donnan_ratio(f64::NAN, 50.0), 1.0);
        // A neutral ratio means a zero potential.
        assert_eq!(donnan_potential_mv(BODY_TEMPERATURE_K, 0.0, 50.0), 0.0);
    }
}

//! Flow-rate-driven duct sizing helpers.

/// Convert a target volumetric flow (CFM) and maximum allowable
/// velocity (FPM) to a *rectangular* duct (w, h) in inches that
/// satisfies the velocity constraint at the closest 1 in-on-a-side
/// square shape. Square is the v1 convention; v2 can route aspect-
/// ratio preferences.
pub fn cfm_to_duct_size(cfm: f64, max_velocity_fpm: f64) -> (f64, f64) {
    if cfm <= 0.0 || max_velocity_fpm <= 0.0 {
        return (0.0, 0.0);
    }
    // Required area in ft^2 = CFM / FPM.
    let area_ft2 = cfm / max_velocity_fpm;
    // Convert to in^2 and pick a square side.
    let area_in2 = area_ft2 * 144.0;
    let side_in = area_in2.sqrt().ceil();
    (side_in, side_in)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cfm_to_duct_size_for_500_cfm_at_700_fpm_is_at_least_10() {
        let (w, h) = cfm_to_duct_size(500.0, 700.0);
        // 500/700 = 0.714 ft² = 102.86 in² → side ≈ 10.14 → ceil to 11 in.
        assert!(w >= 10.0);
        assert!((w - h).abs() < 1e-9);
    }

    #[test]
    fn cfm_to_duct_size_rejects_zero_inputs() {
        assert_eq!(cfm_to_duct_size(0.0, 700.0), (0.0, 0.0));
        assert_eq!(cfm_to_duct_size(500.0, 0.0), (0.0, 0.0));
    }
}

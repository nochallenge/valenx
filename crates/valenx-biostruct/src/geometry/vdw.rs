//! Van der Waals radii by element.
//!
//! A compact Bondi-style radius table (ångström), used by both the
//! steric-clash detector and the Shrake-Rupley SASA sampler. Elements
//! not in the table fall back to a generic 1.70 Å (carbon).

/// Bondi van der Waals radius of an element symbol, ångström.
///
/// The element symbol is matched case-insensitively. Unknown elements
/// return the carbon radius `1.70`.
pub fn vdw_radius(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "H" | "D" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        "SE" => 1.90,
        "B" => 1.92,
        "SI" => 2.10,
        // common metal ions in structures — ionic-ish radii
        "NA" => 2.27,
        "MG" => 1.73,
        "K" => 2.75,
        "CA" => 2.31,
        "MN" => 2.05,
        "FE" => 2.05,
        "CO" => 2.00,
        "NI" => 2.00,
        "CU" => 2.00,
        "ZN" => 2.10,
        _ => 1.70,
    }
}

/// Covalent radius of an element, ångström — used to decide whether
/// two atoms are *bonded* (and therefore should be excluded from
/// clash detection). Pyykkö-style single-bond covalent radii.
pub fn covalent_radius(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "H" | "D" => 0.31,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "BR" => 1.20,
        "I" => 1.39,
        "SE" => 1.20,
        "B" => 0.84,
        "SI" => 1.11,
        "NA" => 1.66,
        "MG" => 1.41,
        "K" => 2.03,
        "CA" => 1.76,
        "MN" => 1.39,
        "FE" => 1.32,
        "CO" => 1.26,
        "NI" => 1.24,
        "CU" => 1.32,
        "ZN" => 1.22,
        _ => 0.77,
    }
}

/// Whether two atoms `a` and `b` separated by `dist` ångström are
/// plausibly covalently bonded: `dist` within `tol` of the sum of
/// their covalent radii.
pub fn is_bonded(elem_a: &str, elem_b: &str, dist: f64, tol: f64) -> bool {
    let ideal = covalent_radius(elem_a) + covalent_radius(elem_b);
    dist > 0.1 && dist <= ideal + tol
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_elements() {
        assert!((vdw_radius("C") - 1.70).abs() < 1e-12);
        assert!((vdw_radius("o") - 1.52).abs() < 1e-12);
        assert!((vdw_radius("ZN") - 2.10).abs() < 1e-12);
    }

    #[test]
    fn unknown_falls_back_to_carbon() {
        assert!((vdw_radius("XX") - 1.70).abs() < 1e-12);
    }

    #[test]
    fn covalent_bond_test() {
        // a typical C-C single bond is ~1.54 A; 0.76+0.76 = 1.52.
        assert!(is_bonded("C", "C", 1.54, 0.45));
        // two carbons 3.5 A apart are not bonded.
        assert!(!is_bonded("C", "C", 3.5, 0.45));
    }
}

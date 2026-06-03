//! Multi-contact electrode: field superposition + current steering.
//!
//! A multi-contact lead (a DBS or high-density neural electrode) drives several
//! contacts with independently programmed current fractions. Because
//! `−∇·(σ∇φ)=I` is **linear**, the field is the superposition of the
//! per-contact fields, and shifting current between contacts **steers** the
//! stimulation focus without physically moving the lead. This module wraps the
//! anisotropic FEM solver to express both.

use crate::aniso_field::{AnisoTissue, SolvedField};

/// A multi-contact electrode embedded in a tissue block: contact node indices
/// driven by programmable current fractions.
pub struct ContactArray<'a> {
    tissue: &'a AnisoTissue,
    contacts: Vec<usize>,
}

impl<'a> ContactArray<'a> {
    /// Build an array from contact node indices (see
    /// [`AnisoTissue::node_at_steps`]).
    pub fn new(tissue: &'a AnisoTissue, contacts: Vec<usize>) -> Self {
        Self { tissue, contacts }
    }

    /// Number of contacts.
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    /// Whether the array has no contacts.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// Solve the grounded-boundary field with total current `total_ua` split
    /// across the contacts by `fractions` (one per contact; typically summing
    /// to one). The result is the superposition of the per-contact fields.
    pub fn solve(&self, total_ua: f64, fractions: &[f64]) -> SolvedField {
        let loads: Vec<(usize, f64)> = self
            .contacts
            .iter()
            .zip(fractions)
            .map(|(&node, &f)| (node, total_ua * f))
            .collect();
        self.tissue.solve_currents(&loads, |_| 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aniso_field::Conductivity;
    use nalgebra::{Matrix3, Vector3};

    fn iso(s: f64) -> Conductivity {
        Matrix3::from_diagonal(&Vector3::new(s, s, s))
    }

    #[test]
    fn fields_superpose_linearly() {
        // φ from two contacts together must equal the sum of the two single-
        // contact fields, to solver tolerance — the linearity steering needs.
        let tissue = AnisoTissue::homogeneous(40.0, 15, iso(0.2));
        let c0 = tissue.node_at_steps(-3, 0, 0);
        let c1 = tissue.node_at_steps(3, 0, 0);
        let both = tissue.solve_currents(&[(c0, 60.0), (c1, 40.0)], |_| 0.0);
        let only0 = tissue.solve_currents(&[(c0, 60.0)], |_| 0.0);
        let only1 = tissue.solve_currents(&[(c1, 40.0)], |_| 0.0);
        let peak = both.phi_mv().iter().cloned().fold(0.0_f64, f64::max);
        let max_dev = both
            .phi_mv()
            .iter()
            .zip(only0.phi_mv())
            .zip(only1.phi_mv())
            .map(|((&b, &a0), &a1)| (b - (a0 + a1)).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_dev < 1.0e-6 * peak.max(1.0),
            "FEM field must superpose (linearity): max deviation {max_dev}, peak {peak}"
        );
    }

    #[test]
    fn current_steering_moves_the_focus() {
        // Two contacts at x = ∓4 mm. As current shifts from the left contact to
        // the right, the field focus must move left → right through centre.
        let tissue = AnisoTissue::homogeneous(40.0, 21, iso(0.2));
        let left = tissue.node_at_steps(-4, 0, 0);
        let right = tissue.node_at_steps(4, 0, 0);
        let array = ContactArray::new(&tissue, vec![left, right]);
        assert_eq!(array.len(), 2);
        let focus = |f_left: f64| array.solve(100.0, &[f_left, 1.0 - f_left]).focus_x_m();
        let left_heavy = focus(0.9);
        let balanced = focus(0.5);
        let right_heavy = focus(0.1);
        assert!(
            left_heavy < balanced && balanced < right_heavy,
            "focus steers left→right as current shifts: {left_heavy:.4} < {balanced:.4} < {right_heavy:.4}"
        );
        assert!(left_heavy < 0.0 && right_heavy > 0.0, "the heavier side sets the focus sign");
        assert!(balanced.abs() < 1.0e-3, "symmetric drive focuses near centre: {balanced}");
    }
}

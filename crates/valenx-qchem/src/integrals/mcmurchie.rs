//! McMurchie-Davidson Hermite-Gaussian recursions.
//!
//! The McMurchie-Davidson (1978) scheme expands a product of two
//! Cartesian Gaussians in *Hermite* Gaussians centred on the product
//! point. Two families of recursions do all the work:
//!
//! - the **`E` coefficients** `E_t^{ij}` — the Hermite-expansion
//!   coefficients of the 1-D product of an `i`-degree Gaussian on `A`
//!   and a `j`-degree Gaussian on `B`;
//! - the **`R` coefficients** `R_{tuv}` — Hermite Coulomb integrals
//!   built from the Boys function, used by the nuclear-attraction and
//!   electron-repulsion integrals.
//!
//! Every primitive integral in the crate — overlap, kinetic, nuclear,
//! ERI, multipole — is assembled from these two recursions, which keeps
//! the integral code small and uniform.

use super::boys::boys_array;

/// The Hermite expansion coefficient `E_t^{ij}` for one Cartesian
/// direction.
///
/// `E(i, j, t)` is the coefficient of the order-`t` Hermite Gaussian in
/// the expansion of `x_A^i x_B^j exp(-a x_A²) exp(-b x_B²)`, where
/// `x_A = x - A`, `x_B = x - B`, the total exponent is `p = a + b`,
/// `q = a b / p` and `AB = A - B`.
///
/// The recursion (Helgaker, Jørgensen, Olsen, *Molecular
/// Electronic-Structure Theory*, eq. 9.5.6) is
///
/// ```text
/// E_t^{i+1,j} = (1/2p) E_{t-1}^{ij} - (q AB / a) E_t^{ij} + (t+1) E_{t+1}^{ij}
/// E_t^{i,j+1} = (1/2p) E_{t-1}^{ij} + (q AB / b) E_t^{ij} + (t+1) E_{t+1}^{ij}
/// ```
///
/// with `E_0^{00} = exp(-q AB²)` and `E_t^{ij} = 0` for `t < 0` or
/// `t > i + j`.
pub fn hermite_e(i: i32, j: i32, t: i32, ab: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    let q = a * b / p;
    if t < 0 || t > i + j {
        0.0
    } else if i == 0 && j == 0 {
        (-q * ab * ab).exp()
    } else if j == 0 {
        // decrement i
        (1.0 / (2.0 * p)) * hermite_e(i - 1, j, t - 1, ab, a, b)
            - (q * ab / a) * hermite_e(i - 1, j, t, ab, a, b)
            + (t as f64 + 1.0) * hermite_e(i - 1, j, t + 1, ab, a, b)
    } else {
        // decrement j
        (1.0 / (2.0 * p)) * hermite_e(i, j - 1, t - 1, ab, a, b)
            + (q * ab / b) * hermite_e(i, j - 1, t, ab, a, b)
            + (t as f64 + 1.0) * hermite_e(i, j - 1, t + 1, ab, a, b)
    }
}

/// The Hermite Coulomb integral `R_{tuv}^{(n)}`.
///
/// `R(t, u, v, n, p, pc)` is the order-`(t,u,v)` Hermite Coulomb
/// integral with auxiliary index `n`, total exponent `p` and
/// centre-separation vector `pc = P - C`. The recursion (Helgaker
/// et al. eq. 9.9.18-20) reduces every index toward
/// `R_{000}^{(n)} = (-2p)^n F_n(p |PC|²)`.
pub fn hermite_r(t: i32, u: i32, v: i32, n: i32, p: f64, pc: [f64; 3], boys: &[f64]) -> f64 {
    if t < 0 || u < 0 || v < 0 {
        return 0.0;
    }
    if t == 0 && u == 0 && v == 0 {
        return (-2.0 * p).powi(n) * boys[n as usize];
    }
    if t > 0 {
        let mut val = 0.0;
        if t > 1 {
            val += (t as f64 - 1.0) * hermite_r(t - 2, u, v, n + 1, p, pc, boys);
        }
        val += pc[0] * hermite_r(t - 1, u, v, n + 1, p, pc, boys);
        val
    } else if u > 0 {
        let mut val = 0.0;
        if u > 1 {
            val += (u as f64 - 1.0) * hermite_r(t, u - 2, v, n + 1, p, pc, boys);
        }
        val += pc[1] * hermite_r(t, u - 1, v, n + 1, p, pc, boys);
        val
    } else {
        let mut val = 0.0;
        if v > 1 {
            val += (v as f64 - 1.0) * hermite_r(t, u, v - 2, n + 1, p, pc, boys);
        }
        val += pc[2] * hermite_r(t, u, v - 1, n + 1, p, pc, boys);
        val
    }
}

/// Build the Boys-function table `F_0 … F_{n_max}` for argument
/// `p · |PC|²` — the input every [`hermite_r`] call tree shares.
pub fn coulomb_boys_table(n_max: i32, p: f64, pc: [f64; 3]) -> Vec<f64> {
    let r2 = pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2];
    boys_array(n_max.max(0) as usize, p * r2)
}

/// The Gaussian product centre `P = (a A + b B) / (a + b)` for two
/// primitives.
#[inline]
pub fn product_centre(a: f64, centre_a: [f64; 3], b: f64, centre_b: [f64; 3]) -> [f64; 3] {
    let p = a + b;
    [
        (a * centre_a[0] + b * centre_b[0]) / p,
        (a * centre_a[1] + b * centre_b[1]) / p,
        (a * centre_a[2] + b * centre_b[2]) / p,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e000_is_the_gaussian_prefactor() {
        // E_0^{00} = exp(-q AB²) with q = ab/(a+b).
        let (a, b, ab) = (1.1, 0.9, 0.6);
        let q = a * b / (a + b);
        let e = hermite_e(0, 0, 0, ab, a, b);
        assert!((e - (-q * ab * ab).exp()).abs() < 1.0e-14);
    }

    #[test]
    fn e_vanishes_above_total_degree() {
        // E_t^{ij} = 0 for t > i+j.
        assert_eq!(hermite_e(1, 1, 3, 0.4, 1.0, 1.0), 0.0);
        assert_eq!(hermite_e(2, 0, 3, 0.4, 1.0, 1.0), 0.0);
    }

    #[test]
    fn e_symmetric_at_zero_separation() {
        // With AB = 0, E_t^{ij} should be symmetric under i<->j when
        // a == b.
        let eij = hermite_e(2, 1, 1, 0.0, 1.3, 1.3);
        let eji = hermite_e(1, 2, 1, 0.0, 1.3, 1.3);
        assert!((eij - eji).abs() < 1.0e-12);
    }

    #[test]
    fn r000_matches_boys_definition() {
        // R_{000}^{(n)} = (-2p)^n F_n(p R²).
        let p = 1.4;
        let pc = [0.3, -0.2, 0.5];
        let boys = coulomb_boys_table(3, p, pc);
        for n in 0..=3 {
            let r = hermite_r(0, 0, 0, n, p, pc, &boys);
            let expect = (-2.0 * p).powi(n) * boys[n as usize];
            assert!((r - expect).abs() < 1.0e-12, "n={n}");
        }
    }

    #[test]
    fn product_centre_midpoint_for_equal_exponents() {
        let p = product_centre(1.0, [0.0, 0.0, 0.0], 1.0, [2.0, 0.0, 0.0]);
        assert!((p[0] - 1.0).abs() < 1.0e-14);
    }
}

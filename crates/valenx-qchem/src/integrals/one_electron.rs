//! One-electron integrals — overlap, kinetic energy, nuclear attraction
//! and Cartesian multipoles.
//!
//! Every routine here builds the requested integral between two
//! [`BasisFunction`]s by looping over their contracted primitives and
//! summing the primitive integral, which is assembled from the
//! McMurchie-Davidson [`hermite_e`] / [`hermite_r`] recursions.
//!
//! ## The primitive formulae
//!
//! For two primitives with exponents `a`, `b`, total exponent
//! `p = a + b`, product centre `P` and Cartesian exponents
//! `(la, ma, na)`, `(lb, mb, nb)`:
//!
//! - **Overlap** factorises into three 1-D overlaps; each is
//!   `E_0` for that direction times `√(π/p)`.
//! - **Kinetic energy** uses the standard `b`-dependent combination of
//!   overlaps with `lb` shifted by `0, ±2`.
//! - **Nuclear attraction** is a triple Hermite sum over the `R`
//!   coefficients, weighted by `-Z` and `2π/p`.
//! - **Multipole** integrals (here the dipole) shift one extra
//!   Cartesian power onto the operator centre.

use super::mcmurchie::{coulomb_boys_table, hermite_e, hermite_r, product_centre};
use crate::basis::BasisFunction;

/// `√π` reused throughout the overlap / kinetic prefactors.
const SQRT_PI: f64 = 1.772_453_850_905_516;

/// The overlap integral `S_{μν} = ⟨μ|ν⟩` between two basis functions.
pub fn overlap(mu: &BasisFunction, nu: &BasisFunction) -> f64 {
    let mut s = 0.0;
    for pa in &mu.primitives {
        for pb in &nu.primitives {
            s += pa.coefficient
                * pb.coefficient
                * primitive_overlap(
                    pa.exponent,
                    mu.centre,
                    mu.cart,
                    pb.exponent,
                    nu.centre,
                    nu.cart,
                );
        }
    }
    s
}

/// The kinetic-energy integral `T_{μν} = ⟨μ|-½∇²|ν⟩`.
pub fn kinetic(mu: &BasisFunction, nu: &BasisFunction) -> f64 {
    let mut t = 0.0;
    for pa in &mu.primitives {
        for pb in &nu.primitives {
            t += pa.coefficient
                * pb.coefficient
                * primitive_kinetic(
                    pa.exponent,
                    mu.centre,
                    mu.cart,
                    pb.exponent,
                    nu.centre,
                    nu.cart,
                );
        }
    }
    t
}

/// The nuclear-attraction integral `V_{μν} = ⟨μ|Σ_C -Z_C/|r-C| |ν⟩`,
/// summed over the point charges `charges` (each `(Z, position)`).
pub fn nuclear_attraction(
    mu: &BasisFunction,
    nu: &BasisFunction,
    charges: &[(f64, [f64; 3])],
) -> f64 {
    let mut v = 0.0;
    for pa in &mu.primitives {
        for pb in &nu.primitives {
            let cc = pa.coefficient * pb.coefficient;
            for &(z, centre_c) in charges {
                v += cc
                    * (-z)
                    * primitive_nuclear(
                        pa.exponent,
                        mu.centre,
                        mu.cart,
                        pb.exponent,
                        nu.centre,
                        nu.cart,
                        centre_c,
                    );
            }
        }
    }
    v
}

/// A Cartesian multipole integral `⟨μ| (x-Ox)^ex (y-Oy)^ey (z-Oz)^ez |ν⟩`
/// about the operator origin `origin`. `(ex, ey, ez)` are the operator
/// powers — `(1,0,0)` is the x-dipole, `(0,0,0)` the overlap.
pub fn multipole(
    mu: &BasisFunction,
    nu: &BasisFunction,
    origin: [f64; 3],
    powers: (u32, u32, u32),
) -> f64 {
    let mut m = 0.0;
    for pa in &mu.primitives {
        for pb in &nu.primitives {
            m += pa.coefficient
                * pb.coefficient
                * primitive_multipole(
                    pa.exponent,
                    mu.centre,
                    mu.cart,
                    pb.exponent,
                    nu.centre,
                    nu.cart,
                    origin,
                    powers,
                );
        }
    }
    m
}

// --- primitive integrals ---------------------------------------------

/// 1-D overlap of two primitive Gaussians of degree `i`, `j`.
fn overlap_1d(i: i32, j: i32, ab: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    hermite_e(i, j, 0, ab, a, b) * (SQRT_PI / p.sqrt())
}

/// Primitive Cartesian overlap — a product of three 1-D overlaps.
#[allow(clippy::too_many_arguments)]
fn primitive_overlap(
    a: f64,
    centre_a: [f64; 3],
    ca: (u32, u32, u32),
    b: f64,
    centre_b: [f64; 3],
    cb: (u32, u32, u32),
) -> f64 {
    let sx = overlap_1d(ca.0 as i32, cb.0 as i32, centre_a[0] - centre_b[0], a, b);
    let sy = overlap_1d(ca.1 as i32, cb.1 as i32, centre_a[1] - centre_b[1], a, b);
    let sz = overlap_1d(ca.2 as i32, cb.2 as i32, centre_a[2] - centre_b[2], a, b);
    sx * sy * sz
}

/// Primitive kinetic-energy integral via the `b`-dependent overlap
/// combination (Helgaker et al. eq. 9.3.40).
#[allow(clippy::too_many_arguments)]
fn primitive_kinetic(
    a: f64,
    centre_a: [f64; 3],
    ca: (u32, u32, u32),
    b: f64,
    centre_b: [f64; 3],
    cb: (u32, u32, u32),
) -> f64 {
    let (la, ma, na) = (ca.0 as i32, ca.1 as i32, ca.2 as i32);
    let (lb, mb, nb) = (cb.0 as i32, cb.1 as i32, cb.2 as i32);
    let abx = centre_a[0] - centre_b[0];
    let aby = centre_a[1] - centre_b[1];
    let abz = centre_a[2] - centre_b[2];

    let sx = overlap_1d(la, lb, abx, a, b);
    let sy = overlap_1d(ma, mb, aby, a, b);
    let sz = overlap_1d(na, nb, abz, a, b);

    // One-dimensional kinetic factors:
    // T_j = b(2j+1) S_j - 2b² S_{j+2} - ½ j(j-1) S_{j-2}.
    let tx = kinetic_1d(la, lb, abx, a, b);
    let ty = kinetic_1d(ma, mb, aby, a, b);
    let tz = kinetic_1d(na, nb, abz, a, b);

    tx * sy * sz + sx * ty * sz + sx * sy * tz
}

/// 1-D kinetic factor `T_j`.
fn kinetic_1d(i: i32, j: i32, ab: f64, a: f64, b: f64) -> f64 {
    let term0 = b * (2.0 * j as f64 + 1.0) * overlap_1d(i, j, ab, a, b);
    let term_up = -2.0 * b * b * overlap_1d(i, j + 2, ab, a, b);
    let term_dn = if j >= 2 {
        -0.5 * (j as f64) * (j as f64 - 1.0) * overlap_1d(i, j - 2, ab, a, b)
    } else {
        0.0
    };
    term0 + term_up + term_dn
}

/// Primitive nuclear-attraction integral for one point charge at `C`
/// (the `-Z` weight is applied by the caller).
#[allow(clippy::too_many_arguments)]
fn primitive_nuclear(
    a: f64,
    centre_a: [f64; 3],
    ca: (u32, u32, u32),
    b: f64,
    centre_b: [f64; 3],
    cb: (u32, u32, u32),
    centre_c: [f64; 3],
) -> f64 {
    let p = a + b;
    let centre_p = product_centre(a, centre_a, b, centre_b);
    let pc = [
        centre_p[0] - centre_c[0],
        centre_p[1] - centre_c[1],
        centre_p[2] - centre_c[2],
    ];
    let (la, ma, na) = (ca.0 as i32, ca.1 as i32, ca.2 as i32);
    let (lb, mb, nb) = (cb.0 as i32, cb.1 as i32, cb.2 as i32);
    let n_max = la + lb + ma + mb + na + nb;
    let boys = coulomb_boys_table(n_max, p, pc);

    let abx = centre_a[0] - centre_b[0];
    let aby = centre_a[1] - centre_b[1];
    let abz = centre_a[2] - centre_b[2];

    let mut val = 0.0;
    for t in 0..=(la + lb) {
        let ex = hermite_e(la, lb, t, abx, a, b);
        if ex == 0.0 {
            continue;
        }
        for u in 0..=(ma + mb) {
            let ey = hermite_e(ma, mb, u, aby, a, b);
            if ey == 0.0 {
                continue;
            }
            for v in 0..=(na + nb) {
                let ez = hermite_e(na, nb, v, abz, a, b);
                if ez == 0.0 {
                    continue;
                }
                val += ex * ey * ez * hermite_r(t, u, v, 0, p, pc, &boys);
            }
        }
    }
    val * 2.0 * std::f64::consts::PI / p
}

/// Primitive Cartesian multipole integral about `origin`.
///
/// A power `(x - Ox)^e` is folded onto basis function `μ`: it shifts
/// the centre to the operator origin and raises the degree, expanded
/// via the binomial `(x-Ox)^e = Σ_k C(e,k) (x-A)^k (A-Ox)^{e-k}`.
#[allow(clippy::too_many_arguments)]
fn primitive_multipole(
    a: f64,
    centre_a: [f64; 3],
    ca: (u32, u32, u32),
    b: f64,
    centre_b: [f64; 3],
    cb: (u32, u32, u32),
    origin: [f64; 3],
    powers: (u32, u32, u32),
) -> f64 {
    let mx = multipole_1d(
        ca.0 as i32,
        cb.0 as i32,
        centre_a[0],
        centre_b[0],
        origin[0],
        powers.0,
        a,
        b,
    );
    let my = multipole_1d(
        ca.1 as i32,
        cb.1 as i32,
        centre_a[1],
        centre_b[1],
        origin[1],
        powers.1,
        a,
        b,
    );
    let mz = multipole_1d(
        ca.2 as i32,
        cb.2 as i32,
        centre_a[2],
        centre_b[2],
        origin[2],
        powers.2,
        a,
        b,
    );
    mx * my * mz
}

/// 1-D multipole integral with operator power `e` about `o`.
#[allow(clippy::too_many_arguments)]
fn multipole_1d(i: i32, j: i32, ax: f64, bx: f64, o: f64, e: u32, a: f64, b: f64) -> f64 {
    // (x-o)^e = Σ_k C(e,k) (x-A)^k (A-o)^{e-k}.
    let ab = ax - bx;
    let a_minus_o = ax - o;
    let mut sum = 0.0;
    for k in 0..=e {
        let coeff = binomial(e, k) * a_minus_o.powi((e - k) as i32);
        sum += coeff * overlap_1d(i + k as i32, j, ab, a, b);
    }
    sum
}

/// Integer binomial coefficient `C(n, k)`.
fn binomial(n: u32, k: u32) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut num = 1.0;
    let mut den = 1.0;
    for i in 0..k {
        num *= (n - i) as f64;
        den *= (i + 1) as f64;
    }
    num / den
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::Primitive;

    /// A single normalised s primitive at `centre` with exponent `a`.
    fn s_primitive(a: f64, centre: [f64; 3]) -> BasisFunction {
        // Normalise the s primitive: N = (2a/π)^{3/4}.
        let n = (2.0 * a / std::f64::consts::PI).powf(0.75);
        BasisFunction {
            atom_index: 0,
            centre,
            cart: (0, 0, 0),
            primitives: vec![Primitive {
                exponent: a,
                coefficient: n,
            }],
        }
    }

    #[test]
    fn s_self_overlap_is_one() {
        let f = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let s = overlap(&f, &f);
        assert!((s - 1.0).abs() < 1.0e-12, "self-overlap = {s}");
    }

    #[test]
    fn overlap_is_symmetric() {
        let f = s_primitive(0.8, [0.0, 0.0, 0.0]);
        let g = s_primitive(1.3, [0.0, 0.0, 1.1]);
        assert!((overlap(&f, &g) - overlap(&g, &f)).abs() < 1.0e-13);
    }

    #[test]
    fn overlap_decays_with_distance() {
        let f = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let near = s_primitive(1.0, [0.0, 0.0, 0.5]);
        let far = s_primitive(1.0, [0.0, 0.0, 3.0]);
        assert!(overlap(&f, &near) > overlap(&f, &far));
        assert!(overlap(&f, &far) > 0.0);
    }

    #[test]
    fn kinetic_self_energy_of_s_gaussian() {
        // ⟨g|-½∇²|g⟩ for a normalised s Gaussian of exponent a equals
        // 3a/2.
        let a = 1.0;
        let f = s_primitive(a, [0.0, 0.0, 0.0]);
        let t = kinetic(&f, &f);
        assert!((t - 1.5 * a).abs() < 1.0e-10, "T = {t}");
    }

    #[test]
    fn kinetic_is_symmetric() {
        let f = s_primitive(0.7, [0.0, 0.0, 0.0]);
        let g = s_primitive(1.1, [0.4, 0.0, 0.9]);
        assert!((kinetic(&f, &g) - kinetic(&g, &f)).abs() < 1.0e-12);
    }

    #[test]
    fn nuclear_attraction_is_negative() {
        // A positive nuclear charge attracts: V < 0.
        let f = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let v = nuclear_attraction(&f, &f, &[(1.0, [0.0, 0.0, 0.0])]);
        assert!(v < 0.0, "V = {v}");
    }

    #[test]
    fn nuclear_attraction_self_value_at_origin() {
        // For a normalised s Gaussian centred on a unit point charge at
        // the origin, V = -⟨g|1/r|g⟩. With g² = (2a/π)^{3/2} e^{-2ar²}
        // and ∫ e^{-2ar²}/r d³r = π/a, this gives
        //   ⟨g|1/r|g⟩ = (2a/π)^{3/2}·π/a = 2√2·√(a/π),
        // so the analytic reference is V = -2√2·√(a/π).
        let a = 1.0;
        let f = s_primitive(a, [0.0, 0.0, 0.0]);
        let v = nuclear_attraction(&f, &f, &[(1.0, [0.0, 0.0, 0.0])]);
        let expect =
            -2.0 * std::f64::consts::SQRT_2 * (a / std::f64::consts::PI).sqrt();
        assert!((v - expect).abs() < 1.0e-9, "V = {v}, expect {expect}");
    }

    #[test]
    fn multipole_with_zero_power_is_overlap() {
        let f = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let g = s_primitive(1.2, [0.0, 0.0, 0.6]);
        let m = multipole(&f, &g, [0.0, 0.0, 0.0], (0, 0, 0));
        assert!((m - overlap(&f, &g)).abs() < 1.0e-12);
    }

    #[test]
    fn dipole_of_centred_gaussian_vanishes() {
        // A Gaussian centred at the origin has zero dipole about the
        // origin by symmetry.
        let f = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let dz = multipole(&f, &f, [0.0, 0.0, 0.0], (0, 0, 1));
        assert!(dz.abs() < 1.0e-12, "dz = {dz}");
    }

    #[test]
    fn dipole_tracks_displacement() {
        // ⟨g|z|g⟩ about the origin for a normalised s Gaussian centred
        // at z0 equals z0.
        let f = s_primitive(1.0, [0.0, 0.0, 1.3]);
        let dz = multipole(&f, &f, [0.0, 0.0, 0.0], (0, 0, 1));
        assert!((dz - 1.3).abs() < 1.0e-10, "dz = {dz}");
    }
}

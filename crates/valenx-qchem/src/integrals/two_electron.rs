//! Two-electron repulsion integrals (ERIs).
//!
//! The electron-repulsion integral in the chemists' notation
//!
//! ```text
//! (μν|λσ) = ∫∫ μ(r₁) ν(r₁) (1/r₁₂) λ(r₂) σ(r₂) dr₁ dr₂
//! ```
//!
//! is the most expensive object in a Hartree-Fock calculation. This
//! module assembles it with the McMurchie-Davidson scheme: the bra pair
//! `μν` and the ket pair `λσ` are each Hermite-expanded, and the two
//! Hermite sets are contracted through the Hermite Coulomb integrals
//! `R_{tuv}`.
//!
//! ## The tensor
//!
//! [`EriTensor`] stores the full four-index tensor `(μν|λσ)` for an
//! `n`-function basis as a dense `n⁴` array. With the small basis sets
//! this crate targets (a few dozen functions) the dense tensor is a few
//! megabytes and the simplest correct choice.
//!
//! ## v1 performance caveat
//!
//! The build loops every `(μ,ν,λ,σ)` quartet and only uses the trivial
//! `(μν|λσ) = (νμ|λσ) = (μν|σλ) = (λσ|μν)` 8-fold permutational
//! symmetry to fill equivalent entries from one computed value. There
//! is no Schwarz-inequality screening and no shell-pair blocking, so
//! the cost is `O(n⁴)` primitive-integral evaluations — fine for the
//! small-molecule regime, not competitive with a production integral
//! engine.

use super::mcmurchie::{coulomb_boys_table, hermite_e, hermite_r, product_centre};
use crate::basis::{BasisFunction, BasisSet};

/// `2 π^{5/2}` — the universal prefactor of the McMurchie-Davidson ERI.
const TWO_PI_POW_5_2: f64 = 34.986_836_655_249_725;

/// The dense four-index electron-repulsion tensor for a basis set.
#[derive(Clone, Debug)]
pub struct EriTensor {
    /// Number of basis functions `n`.
    pub n: usize,
    /// Flat `n⁴` storage in row-major `[μ][ν][λ][σ]` order.
    data: Vec<f64>,
}

impl EriTensor {
    /// Allocate a zeroed tensor for an `n`-function basis.
    pub fn zeros(n: usize) -> Self {
        EriTensor {
            n,
            data: vec![0.0; n * n * n * n],
        }
    }

    /// Flat index of `(μ, ν, λ, σ)`.
    #[inline]
    fn idx(&self, mu: usize, nu: usize, la: usize, si: usize) -> usize {
        ((mu * self.n + nu) * self.n + la) * self.n + si
    }

    /// The integral `(μν|λσ)`.
    #[inline]
    pub fn get(&self, mu: usize, nu: usize, la: usize, si: usize) -> f64 {
        self.data[self.idx(mu, nu, la, si)]
    }

    /// Store `value` into `(μν|λσ)`.
    #[inline]
    pub fn set(&mut self, mu: usize, nu: usize, la: usize, si: usize, value: f64) {
        let i = self.idx(mu, nu, la, si);
        self.data[i] = value;
    }

    /// Compute the full ERI tensor for `basis`.
    ///
    /// The unique quartets `μ ≥ ν`, `λ ≥ σ`, `(μν) ≥ (λσ)` are
    /// evaluated once and the 8-fold permutational symmetry mirrors
    /// each value into its equivalents.
    pub fn build(basis: &BasisSet) -> Self {
        let n = basis.n_functions();
        let mut tensor = EriTensor::zeros(n);
        let f = &basis.functions;
        for mu in 0..n {
            for nu in 0..=mu {
                let bra = mu * (mu + 1) / 2 + nu;
                for la in 0..n {
                    for si in 0..=la {
                        let ket = la * (la + 1) / 2 + si;
                        if ket > bra {
                            continue;
                        }
                        let value = eri(&f[mu], &f[nu], &f[la], &f[si]);
                        // 8-fold permutational symmetry.
                        for &(i, j, k, l) in &[
                            (mu, nu, la, si),
                            (nu, mu, la, si),
                            (mu, nu, si, la),
                            (nu, mu, si, la),
                            (la, si, mu, nu),
                            (si, la, mu, nu),
                            (la, si, nu, mu),
                            (si, la, nu, mu),
                        ] {
                            tensor.set(i, j, k, l, value);
                        }
                    }
                }
            }
        }
        tensor
    }
}

/// The contracted electron-repulsion integral `(μν|λσ)`.
pub fn eri(
    mu: &BasisFunction,
    nu: &BasisFunction,
    la: &BasisFunction,
    si: &BasisFunction,
) -> f64 {
    let mut total = 0.0;
    for pa in &mu.primitives {
        for pb in &nu.primitives {
            for pc in &la.primitives {
                for pd in &si.primitives {
                    total += pa.coefficient
                        * pb.coefficient
                        * pc.coefficient
                        * pd.coefficient
                        * primitive_eri(
                            (pa.exponent, mu.centre, mu.cart),
                            (pb.exponent, nu.centre, nu.cart),
                            (pc.exponent, la.centre, la.cart),
                            (pd.exponent, si.centre, si.cart),
                        );
                }
            }
        }
    }
    total
}

/// One primitive ERI via the McMurchie-Davidson Hermite contraction.
///
/// Each tuple is `(exponent, centre, cartesian-exponents)`.
fn primitive_eri(
    bra_a: (f64, [f64; 3], (u32, u32, u32)),
    bra_b: (f64, [f64; 3], (u32, u32, u32)),
    ket_c: (f64, [f64; 3], (u32, u32, u32)),
    ket_d: (f64, [f64; 3], (u32, u32, u32)),
) -> f64 {
    let (a, centre_a, ca) = bra_a;
    let (b, centre_b, cb) = bra_b;
    let (c, centre_c, cc) = ket_c;
    let (d, centre_d, cd) = ket_d;

    let p = a + b;
    let q = c + d;
    let centre_p = product_centre(a, centre_a, b, centre_b);
    let centre_q = product_centre(c, centre_c, d, centre_d);
    let alpha = p * q / (p + q);
    let pq = [
        centre_p[0] - centre_q[0],
        centre_p[1] - centre_q[1],
        centre_p[2] - centre_q[2],
    ];

    // Bra / ket Cartesian exponents.
    let (la, ma, na) = (ca.0 as i32, ca.1 as i32, ca.2 as i32);
    let (lb, mb, nb) = (cb.0 as i32, cb.1 as i32, cb.2 as i32);
    let (lc, mc, nc) = (cc.0 as i32, cc.1 as i32, cc.2 as i32);
    let (ld, md, nd) = (cd.0 as i32, cd.1 as i32, cd.2 as i32);

    let ab = [
        centre_a[0] - centre_b[0],
        centre_a[1] - centre_b[1],
        centre_a[2] - centre_b[2],
    ];
    let cd_sep = [
        centre_c[0] - centre_d[0],
        centre_c[1] - centre_d[1],
        centre_c[2] - centre_d[2],
    ];

    let n_max = la + lb + ma + mb + na + nb + lc + ld + mc + md + nc + nd;
    let boys = coulomb_boys_table(n_max, alpha, pq);

    let mut total = 0.0;
    // Bra Hermite indices t1,u1,v1 ; ket Hermite indices t2,u2,v2.
    for t1 in 0..=(la + lb) {
        let ex1 = hermite_e(la, lb, t1, ab[0], a, b);
        if ex1 == 0.0 {
            continue;
        }
        for u1 in 0..=(ma + mb) {
            let ey1 = hermite_e(ma, mb, u1, ab[1], a, b);
            if ey1 == 0.0 {
                continue;
            }
            for v1 in 0..=(na + nb) {
                let ez1 = hermite_e(na, nb, v1, ab[2], a, b);
                if ez1 == 0.0 {
                    continue;
                }
                let bra = ex1 * ey1 * ez1;
                for t2 in 0..=(lc + ld) {
                    let ex2 = hermite_e(lc, ld, t2, cd_sep[0], c, d);
                    if ex2 == 0.0 {
                        continue;
                    }
                    for u2 in 0..=(mc + md) {
                        let ey2 = hermite_e(mc, md, u2, cd_sep[1], c, d);
                        if ey2 == 0.0 {
                            continue;
                        }
                        for v2 in 0..=(nc + nd) {
                            let ez2 = hermite_e(nc, nd, v2, cd_sep[2], c, d);
                            if ez2 == 0.0 {
                                continue;
                            }
                            // Ket Hermite functions enter with sign
                            // (-1)^{t2+u2+v2}.
                            let sign = if (t2 + u2 + v2) % 2 == 0 {
                                1.0
                            } else {
                                -1.0
                            };
                            let r = hermite_r(
                                t1 + t2,
                                u1 + u2,
                                v1 + v2,
                                0,
                                alpha,
                                pq,
                                &boys,
                            );
                            total += bra * ex2 * ey2 * ez2 * sign * r;
                        }
                    }
                }
            }
        }
    }

    total * TWO_PI_POW_5_2 / (p * q * (p + q).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::Primitive;

    fn s_primitive(a: f64, centre: [f64; 3]) -> BasisFunction {
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
    fn eri_of_one_centre_s_gaussian() {
        // The self-repulsion (ss|ss) of a normalised s Gaussian of
        // exponent a, all four on the same centre. The McMurchie-Davidson
        // closed form is
        //   (ss|ss) = (2a/π)³ · 2π^{5/2} / [(2a)(2a)√(4a)]
        //           = 2·√(a/π).
        // (See Szabo & Ostlund, Appendix A — the two-electron integral
        //  over identical 1s Gaussians.) For a = 1, (ss|ss) = 2/√π.
        let a = 1.0;
        let g = s_primitive(a, [0.0, 0.0, 0.0]);
        let v = eri(&g, &g, &g, &g);
        let expect = 2.0 * (a / std::f64::consts::PI).sqrt();
        assert!((v - expect).abs() < 1.0e-9, "(ss|ss) = {v}, expect {expect}");
    }

    #[test]
    fn eri_is_positive() {
        let g = s_primitive(0.9, [0.0, 0.0, 0.0]);
        let h = s_primitive(1.4, [0.0, 0.0, 1.2]);
        assert!(eri(&g, &g, &h, &h) > 0.0);
    }

    #[test]
    fn eri_permutational_symmetry() {
        let f0 = s_primitive(0.8, [0.0, 0.0, 0.0]);
        let f1 = s_primitive(1.1, [0.0, 0.0, 0.7]);
        let f2 = s_primitive(0.6, [0.5, 0.0, 0.0]);
        let f3 = s_primitive(1.3, [0.0, 0.4, 0.0]);
        let base = eri(&f0, &f1, &f2, &f3);
        // (μν|λσ) = (νμ|λσ) = (μν|σλ) = (λσ|μν).
        assert!((base - eri(&f1, &f0, &f2, &f3)).abs() < 1.0e-12);
        assert!((base - eri(&f0, &f1, &f3, &f2)).abs() < 1.0e-12);
        assert!((base - eri(&f2, &f3, &f0, &f1)).abs() < 1.0e-12);
    }

    #[test]
    fn eri_decays_with_pair_separation() {
        let g = s_primitive(1.0, [0.0, 0.0, 0.0]);
        let near = s_primitive(1.0, [0.0, 0.0, 1.0]);
        let far = s_primitive(1.0, [0.0, 0.0, 6.0]);
        assert!(eri(&g, &g, &near, &near) > eri(&g, &g, &far, &far));
    }

    #[test]
    fn tensor_build_respects_symmetry() {
        use crate::basis::{AngularMomentum, BasisSet, Shell};
        let shell = |a: f64, c: [f64; 3]| Shell {
            atom_index: 0,
            centre: c,
            angular: AngularMomentum::S,
            primitives: vec![Primitive {
                exponent: a,
                coefficient: 1.0,
            }],
        }
        .normalised();
        let shells = vec![
            shell(1.0, [0.0, 0.0, 0.0]),
            shell(1.2, [0.0, 0.0, 1.4]),
        ];
        let functions: Vec<BasisFunction> = shells
            .iter()
            .map(|s| BasisFunction {
                atom_index: s.atom_index,
                centre: s.centre,
                cart: (0, 0, 0),
                primitives: s.primitives.clone(),
            })
            .collect();
        let basis = BasisSet {
            name: "test",
            shells,
            functions,
        };
        let t = EriTensor::build(&basis);
        for mu in 0..2 {
            for nu in 0..2 {
                for la in 0..2 {
                    for si in 0..2 {
                        let v = t.get(mu, nu, la, si);
                        assert!((v - t.get(nu, mu, la, si)).abs() < 1.0e-13);
                        assert!((v - t.get(la, si, mu, nu)).abs() < 1.0e-13);
                    }
                }
            }
        }
    }
}

//! Lebedev angular quadrature.
//!
//! A Lebedev quadrature integrates a function over the surface of the
//! unit sphere — `∫ f dΩ ≈ 4π Σ_k w_k f(Ω_k)` — exactly for every
//! spherical harmonic up to a given order. It is the natural angular
//! half of an atom-centred molecular grid: the radial quadrature
//! supplies the `r` points, the Lebedev quadrature the `(θ, φ)`
//! directions, and the two tensor-multiply.
//!
//! ## The grids shipped
//!
//! Lebedev grids exist only for a fixed set of point counts. This
//! module ships four octahedral grids — **6, 26, 50 and 110 points**,
//! exact for spherical harmonics through degree 3, 7, 11 and 17 —
//! generated from the octahedral *generator* orbits with the published
//! Lebedev weights (V. I. Lebedev, *Zh. Vychisl. Mat. Mat. Fiz.*; the
//! values here match the canonical `Lebedev-Laikov` tables). Each grid
//! is built once from its generators, so the point list is exact to
//! machine precision and the weights sum to 1.
//!
//! Still-larger Lebedev grids (194, 302, 590 …) follow the identical
//! generator pattern with more orbits; the four here span the
//! coarse / medium / fine / converged spectrum a small-molecule DFT
//! grid needs and keep the table compact.

/// A single Lebedev quadrature point: a unit direction and its weight.
///
/// The weights of a grid sum to `1`; the spherical integral is
/// `4π Σ_k w_k f(p_k)`.
#[derive(Copy, Clone, Debug)]
pub struct LebedevPoint {
    /// Unit direction `(x, y, z)` on the sphere.
    pub dir: [f64; 3],
    /// Quadrature weight (the grid's weights sum to 1).
    pub weight: f64,
}

/// Push the six `(±1,0,0)`-type octahedral vertices, each with weight
/// `w`.
fn add_a1(out: &mut Vec<LebedevPoint>, w: f64) {
    for &s in &[1.0, -1.0] {
        out.push(LebedevPoint {
            dir: [s, 0.0, 0.0],
            weight: w,
        });
        out.push(LebedevPoint {
            dir: [0.0, s, 0.0],
            weight: w,
        });
        out.push(LebedevPoint {
            dir: [0.0, 0.0, s],
            weight: w,
        });
    }
}

/// Push the twelve `(±1,±1,0)/√2`-type octahedral edge points, each
/// with weight `w`.
fn add_a2(out: &mut Vec<LebedevPoint>, w: f64) {
    let c = std::f64::consts::FRAC_1_SQRT_2;
    for &sx in &[1.0, -1.0] {
        for &sy in &[1.0, -1.0] {
            out.push(LebedevPoint {
                dir: [sx * c, sy * c, 0.0],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [sx * c, 0.0, sy * c],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [0.0, sx * c, sy * c],
                weight: w,
            });
        }
    }
}

/// Push the eight `(±1,±1,±1)/√3`-type octahedral face points, each
/// with weight `w`.
fn add_a3(out: &mut Vec<LebedevPoint>, w: f64) {
    let c = 1.0 / 3.0_f64.sqrt();
    for &sx in &[1.0, -1.0] {
        for &sy in &[1.0, -1.0] {
            for &sz in &[1.0, -1.0] {
                out.push(LebedevPoint {
                    dir: [sx * c, sy * c, sz * c],
                    weight: w,
                });
            }
        }
    }
}

/// Push the 24 points of a `(l, l, m)` generator orbit (`m² = 1−2l²`),
/// each with weight `w`.
fn add_bk(out: &mut Vec<LebedevPoint>, l: f64, w: f64) {
    let m = (1.0 - 2.0 * l * l).max(0.0).sqrt();
    // The 24 sign / axis permutations of (l, l, m) with the large
    // component m on each of the three axes.
    for &sl1 in &[1.0, -1.0] {
        for &sl2 in &[1.0, -1.0] {
            for &sm in &[1.0, -1.0] {
                out.push(LebedevPoint {
                    dir: [sl1 * l, sl2 * l, sm * m],
                    weight: w,
                });
                out.push(LebedevPoint {
                    dir: [sl1 * l, sm * m, sl2 * l],
                    weight: w,
                });
                out.push(LebedevPoint {
                    dir: [sm * m, sl1 * l, sl2 * l],
                    weight: w,
                });
            }
        }
    }
}

/// Push the 24 points of a `(p, q, 0)` generator orbit (`p² + q² = 1`,
/// `p ≠ q`), each with weight `w`.
fn add_ck(out: &mut Vec<LebedevPoint>, p: f64, w: f64) {
    let q = (1.0 - p * p).max(0.0).sqrt();
    // 24 points: (p, q, 0) and (q, p, 0) with the zero on each axis.
    for &sp in &[1.0, -1.0] {
        for &sq in &[1.0, -1.0] {
            out.push(LebedevPoint {
                dir: [sp * p, sq * q, 0.0],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [sp * q, sq * p, 0.0],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [sp * p, 0.0, sq * q],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [sp * q, 0.0, sq * p],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [0.0, sp * p, sq * q],
                weight: w,
            });
            out.push(LebedevPoint {
                dir: [0.0, sp * q, sq * p],
                weight: w,
            });
        }
    }
}

/// The 6-point Lebedev grid — the octahedron vertices. Exact for
/// spherical harmonics through degree 3.
pub fn lebedev_6() -> Vec<LebedevPoint> {
    let mut v = Vec::with_capacity(6);
    add_a1(&mut v, 1.0 / 6.0);
    v
}

/// The 26-point Lebedev grid (`a1 + a2 + a3` orbits). Exact through
/// degree 7.
pub fn lebedev_26() -> Vec<LebedevPoint> {
    let mut v = Vec::with_capacity(26);
    add_a1(&mut v, 1.0 / 21.0);
    add_a2(&mut v, 4.0 / 105.0);
    add_a3(&mut v, 27.0 / 840.0);
    v
}

/// The 50-point Lebedev grid (`a1 + a2 + a3` plus two `bk` orbits).
/// Exact through degree 11 — the Lebedev-Laikov order-11 grid.
pub fn lebedev_50() -> Vec<LebedevPoint> {
    let mut v = Vec::with_capacity(50);
    add_a1(&mut v, 0.012_698_412_698_412_698);
    add_a2(&mut v, 0.022_574_955_908_289_24);
    add_a3(&mut v, 0.021_093_75);
    add_bk(&mut v, 0.301_511_344_577_763_74, 0.020_173_335_537_918_87);
    v
}

/// The 110-point Lebedev grid (`a1 + a3`, three `bk` orbits and one
/// `ck` orbit). Exact through degree 17 — the Lebedev-Laikov order-17
/// grid `LD0110`.
pub fn lebedev_110() -> Vec<LebedevPoint> {
    let mut v = Vec::with_capacity(110);
    add_a1(&mut v, 0.003_828_270_494_937_162);
    add_a3(&mut v, 0.009_793_737_512_487_513);
    add_bk(&mut v, 0.185_115_635_334_736_2, 0.008_211_737_283_191_111);
    add_bk(&mut v, 0.690_421_048_382_292_2, 0.009_942_814_891_178_103);
    add_bk(&mut v, 0.395_689_473_055_941_9, 0.009_595_471_336_070_962);
    add_ck(&mut v, 0.478_369_028_812_150_2, 0.009_694_996_361_663_029);
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Lebedev grid's weights must sum to 1 (so `4π Σ w = 4π`, the
    /// surface area of the unit sphere).
    fn weights_sum_to_one(grid: &[LebedevPoint]) {
        let s: f64 = grid.iter().map(|p| p.weight).sum();
        assert!((s - 1.0).abs() < 1.0e-12, "Σw = {s}");
    }

    #[test]
    fn all_points_are_unit_vectors() {
        for grid in [lebedev_6(), lebedev_26(), lebedev_50(), lebedev_110()] {
            for p in &grid {
                let r2 = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1] + p.dir[2] * p.dir[2];
                assert!((r2 - 1.0).abs() < 1.0e-12, "|dir|² = {r2}");
            }
        }
    }

    #[test]
    fn grid_sizes_are_as_named() {
        assert_eq!(lebedev_6().len(), 6);
        assert_eq!(lebedev_26().len(), 26);
        assert_eq!(lebedev_50().len(), 50);
        assert_eq!(lebedev_110().len(), 110);
    }

    #[test]
    fn weights_normalised() {
        weights_sum_to_one(&lebedev_6());
        weights_sum_to_one(&lebedev_26());
        weights_sum_to_one(&lebedev_50());
        weights_sum_to_one(&lebedev_110());
    }

    /// Integrate a polynomial in the direction cosines and compare with
    /// the exact average over the sphere. `⟨x²⟩ = 1/3`, `⟨x⁴⟩ = 1/5`,
    /// `⟨x²y²⟩ = 1/15` — a 26-point grid (degree 7) reproduces all
    /// three exactly.
    #[test]
    fn integrates_low_degree_polynomials_exactly() {
        let grid = lebedev_26();
        let avg = |f: &dyn Fn([f64; 3]) -> f64| -> f64 {
            grid.iter().map(|p| p.weight * f(p.dir)).sum::<f64>()
        };
        assert!((avg(&|d| d[0] * d[0]) - 1.0 / 3.0).abs() < 1.0e-12);
        assert!((avg(&|d| d[0].powi(4)) - 1.0 / 5.0).abs() < 1.0e-12);
        assert!((avg(&|d| d[0] * d[0] * d[1] * d[1]) - 1.0 / 15.0).abs() < 1.0e-12);
        // An odd polynomial integrates to zero by symmetry.
        assert!(avg(&|d| d[0].powi(3)).abs() < 1.0e-12);
    }

    /// The 50-point grid is exact through degree 11, so it must also
    /// reproduce `⟨x⁶⟩ = 1/7` and `⟨x⁸⟩ = 1/9`.
    #[test]
    fn fifty_point_grid_is_exact_through_degree_eleven() {
        let grid = lebedev_50();
        let avg = |f: &dyn Fn([f64; 3]) -> f64| -> f64 {
            grid.iter().map(|p| p.weight * f(p.dir)).sum::<f64>()
        };
        assert!((avg(&|d| d[0].powi(6)) - 1.0 / 7.0).abs() < 1.0e-10);
        assert!((avg(&|d| d[0].powi(8)) - 1.0 / 9.0).abs() < 1.0e-10);
        // Mixed degree-8 monomial: ⟨x⁴y⁴⟩ = 9/(7·9·11·... ) — use the
        // closed form 3·3/(945) = 1/105.
        assert!((avg(&|d| d[0].powi(4) * d[1].powi(4)) - 1.0 / 105.0).abs() < 1.0e-10);
    }

    /// The 110-point grid is exact through degree 17 — it must
    /// reproduce `⟨x¹⁴⟩ = 1/15` and `⟨x¹⁶⟩ = 1/17`, which the smaller
    /// grids cannot.
    #[test]
    fn hundred_ten_point_grid_is_exact_through_degree_seventeen() {
        let grid = lebedev_110();
        let avg = |f: &dyn Fn([f64; 3]) -> f64| -> f64 {
            grid.iter().map(|p| p.weight * f(p.dir)).sum::<f64>()
        };
        for k in (2..=16).step_by(2) {
            let exact = 1.0 / (k as f64 + 1.0);
            let got = avg(&|d| d[0].powi(k));
            assert!(
                (got - exact).abs() < 1.0e-10,
                "⟨x^{k}⟩ = {got}, expected {exact}"
            );
        }
        // An odd monomial still vanishes by symmetry.
        assert!(avg(&|d| d[0].powi(15)).abs() < 1.0e-10);
    }
}

//! The curvature engine: from a metric to its full Riemann curvature.
//!
//! Given a [`Spacetime`] and a point, [`curvature_at`] returns the Christoffel
//! symbols, the Riemann/Ricci/Einstein tensors and the Ricci and Kretschmann
//! scalars — all in the coordinate basis.
//!
//! ## Method
//!
//! 1. Evaluate the metric `g_μν` (plain `f64`) and invert it for `g^{μν}`.
//! 2. Obtain the first partials `∂_c g_{ab}` and second partials
//!    `∂_c∂_d g_{ab}` *exactly* via one [`HyperDual`] pass per coordinate pair.
//! 3. Assemble the Christoffel symbols and their gradients algebraically from
//!    those partials — we never differentiate through the matrix inverse, which
//!    is the standard, numerically robust route.
//! 4. Build the Riemann tensor, then contract for Ricci, the Ricci scalar, the
//!    Einstein tensor and the Kretschmann invariant.
//!
//! Index conventions: `christoffel[a][b][c] = Γ^a_{bc}`,
//! `riemann[a][b][c][d] = R^a_{bcd}`, `ricci[a][b] = R_{ab}`.

// These are dense rank-3/4 tensor contractions over the 4 spacetime indices.
// Explicit `for i in 0..4` index loops mirror the Einstein-summation formulas
// one-to-one and index several tensors per iteration; rewriting them with
// `enumerate()` (as clippy::needless_range_loop suggests) would obscure the
// math and is error-prone, so the lint is disabled for this module.
#![allow(clippy::needless_range_loop)]

use crate::autodiff::{Dual, HyperDual};
use crate::metric::Spacetime;
use crate::tensor::{self, Mat4};
use crate::{RelativityError, Result};

type T3 = [[[f64; 4]; 4]; 4];
type T4 = [[[[f64; 4]; 4]; 4]; 4];

/// The local curvature of a spacetime at one coordinate point.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Curvature {
    /// Christoffel symbols of the second kind, `Γ^a_{bc}` (indices `[a][b][c]`).
    pub christoffel: T3,
    /// Riemann curvature tensor `R^a_{bcd}` (mixed; indices `[a][b][c][d]`).
    pub riemann: T4,
    /// Ricci tensor `R_{ab}`.
    pub ricci: Mat4,
    /// Ricci scalar `R = g^{ab} R_{ab}`.
    pub ricci_scalar: f64,
    /// Einstein tensor `G_{ab} = R_{ab} − ½ R g_{ab}`.
    pub einstein: Mat4,
    /// Kretschmann scalar `K = R_{abcd} R^{abcd}`.
    pub kretschmann: f64,
}

/// Compute the full curvature of `spacetime` at coordinate point `x`.
///
/// # Errors
/// Returns [`RelativityError::CoordinateSingularity`] if the metric is
/// non-finite or non-invertible at `x` (e.g. on a horizon or coordinate axis),
/// rather than producing silent NaNs.
pub fn curvature_at(spacetime: &impl Spacetime, x: [f64; 4]) -> Result<Curvature> {
    let g = spacetime.metric::<f64>(x);
    if !tensor::all_finite(&g) {
        return Err(RelativityError::CoordinateSingularity(format!(
            "metric is non-finite at {x:?}"
        )));
    }
    let ginv = tensor::inverse(&g).ok_or_else(|| {
        RelativityError::CoordinateSingularity(format!("metric is singular at {x:?}"))
    })?;

    let (dg, d2g) = metric_partials(spacetime, x);

    // ∂_c g^{ab} = − g^{ae} (∂_c g_{ef}) g^{fb}.
    let mut dginv: T3 = [[[0.0; 4]; 4]; 4];
    for (c, dgc) in dg.iter().enumerate() {
        dginv[c] = neg_inv_deriv(&ginv, dgc);
    }

    let gamma = christoffel(&ginv, &dg);
    let dgamma = christoffel_grad(&ginv, &dginv, &dg, &d2g);
    let riemann = riemann(&gamma, &dgamma);
    let ricci = ricci_from_riemann(&riemann);
    let ricci_scalar = contract2(&ginv, &ricci);
    let einstein = einstein_tensor(&ricci, ricci_scalar, &g);
    let kretschmann = kretschmann(&riemann, &g, &ginv);

    if !ricci_scalar.is_finite() || !kretschmann.is_finite() {
        return Err(RelativityError::CoordinateSingularity(format!(
            "curvature is non-finite at {x:?}"
        )));
    }

    Ok(Curvature {
        christoffel: gamma,
        riemann,
        ricci,
        ricci_scalar,
        einstein,
        kretschmann,
    })
}

/// First partials `dg[c][a][b] = ∂_c g_{ab}` and second partials
/// `d2g[c][d][a][b] = ∂_c∂_d g_{ab}`, via one hyper-dual evaluation per
/// coordinate pair `(c, d)` with `c ≤ d` (the rest follow by symmetry).
fn metric_partials(s: &impl Spacetime, x: [f64; 4]) -> (T3, T4) {
    let mut dg: T3 = [[[0.0; 4]; 4]; 4];
    let mut d2g: T4 = [[[[0.0; 4]; 4]; 4]; 4];
    for c in 0..4 {
        for d in c..4 {
            let xhd = [
                HyperDual::seed(x[0], c == 0, d == 0),
                HyperDual::seed(x[1], c == 1, d == 1),
                HyperDual::seed(x[2], c == 2, d == 2),
                HyperDual::seed(x[3], c == 3, d == 3),
            ];
            let ghd = s.metric::<HyperDual>(xhd);
            for a in 0..4 {
                for b in 0..4 {
                    let comp = ghd[a][b];
                    dg[c][a][b] = comp.d1; // ∂_c g_ab
                    dg[d][a][b] = comp.d2; // ∂_d g_ab
                    d2g[c][d][a][b] = comp.d12; // ∂_c∂_d g_ab
                    d2g[d][c][a][b] = comp.d12;
                }
            }
        }
    }
    (dg, d2g)
}

/// Christoffel symbols `Γ^a_{bc}` at a point.
///
/// Cheaper than [`curvature_at`]: only the metric's first derivatives are
/// needed (4 dual-number passes), so this is what the geodesic integrator calls
/// on every step.
///
/// # Errors
/// [`RelativityError::CoordinateSingularity`] if the metric is non-finite or
/// non-invertible at `x`.
pub fn christoffel_at(spacetime: &impl Spacetime, x: [f64; 4]) -> Result<T3> {
    let g = spacetime.metric::<f64>(x);
    if !tensor::all_finite(&g) {
        return Err(RelativityError::CoordinateSingularity(format!(
            "metric is non-finite at {x:?}"
        )));
    }
    let ginv = tensor::inverse(&g).ok_or_else(|| {
        RelativityError::CoordinateSingularity(format!("metric is singular at {x:?}"))
    })?;
    Ok(christoffel(&ginv, &first_partials(spacetime, x)))
}

/// First partials `dg[c][a][b] = ∂_c g_{ab}` via 4 dual-number passes (one per
/// coordinate). Used where second derivatives are not required.
fn first_partials(s: &impl Spacetime, x: [f64; 4]) -> T3 {
    let mut dg: T3 = [[[0.0; 4]; 4]; 4];
    for c in 0..4 {
        let xd: [Dual; 4] = std::array::from_fn(|i| {
            if i == c {
                Dual::variable(x[i])
            } else {
                Dual::constant(x[i])
            }
        });
        let gd = s.metric::<Dual>(xd);
        for a in 0..4 {
            for b in 0..4 {
                dg[c][a][b] = gd[a][b].d;
            }
        }
    }
    dg
}

/// `−g⁻¹ (∂_c g) g⁻¹`, i.e. `∂_c g^{ab} = −g^{ae}(∂_c g_{ef})g^{fb}`.
fn neg_inv_deriv(ginv: &Mat4, dgc: &Mat4) -> Mat4 {
    let mut out = [[0.0; 4]; 4];
    for (a, outa) in out.iter_mut().enumerate() {
        for (b, outab) in outa.iter_mut().enumerate() {
            let mut s = 0.0;
            for e in 0..4 {
                for f in 0..4 {
                    s += ginv[a][e] * dgc[e][f] * ginv[f][b];
                }
            }
            *outab = -s;
        }
    }
    out
}

/// Christoffel symbols `Γ^a_{bc} = ½ g^{ad}(∂_b g_{dc} + ∂_c g_{db} − ∂_d g_{bc})`.
fn christoffel(ginv: &Mat4, dg: &T3) -> T3 {
    let mut gamma: T3 = [[[0.0; 4]; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                let mut s = 0.0;
                for d in 0..4 {
                    s += ginv[a][d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]);
                }
                gamma[a][b][c] = 0.5 * s;
            }
        }
    }
    gamma
}

/// Coordinate gradient of the Christoffel symbols, `dgamma[e][a][b][c] = ∂_e Γ^a_{bc}`.
fn christoffel_grad(ginv: &Mat4, dginv: &T3, dg: &T3, d2g: &T4) -> T4 {
    let mut out: T4 = [[[[0.0; 4]; 4]; 4]; 4];
    for e in 0..4 {
        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    let mut s = 0.0;
                    for d in 0..4 {
                        let inner = dg[b][d][c] + dg[c][d][b] - dg[d][b][c];
                        let inner_grad = d2g[e][b][d][c] + d2g[e][c][d][b] - d2g[e][d][b][c];
                        s += dginv[e][a][d] * inner + ginv[a][d] * inner_grad;
                    }
                    out[e][a][b][c] = 0.5 * s;
                }
            }
        }
    }
    out
}

/// Riemann tensor
/// `R^a_{bcd} = ∂_c Γ^a_{bd} − ∂_d Γ^a_{bc} + Γ^a_{ce}Γ^e_{bd} − Γ^a_{de}Γ^e_{bc}`.
fn riemann(gamma: &T3, dgamma: &T4) -> T4 {
    let mut r: T4 = [[[[0.0; 4]; 4]; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let mut val = dgamma[c][a][b][d] - dgamma[d][a][b][c];
                    for e in 0..4 {
                        val += gamma[a][c][e] * gamma[e][b][d] - gamma[a][d][e] * gamma[e][b][c];
                    }
                    r[a][b][c][d] = val;
                }
            }
        }
    }
    r
}

/// Ricci tensor `R_{bd} = R^a_{bad}` (contract the first and third indices).
fn ricci_from_riemann(r: &T4) -> Mat4 {
    let mut ric = [[0.0; 4]; 4];
    for b in 0..4 {
        for d in 0..4 {
            let mut s = 0.0;
            for a in 0..4 {
                s += r[a][b][a][d];
            }
            ric[b][d] = s;
        }
    }
    ric
}

/// Contract a symmetric rank-2 tensor with the inverse metric: `g^{ab} M_{ab}`.
fn contract2(ginv: &Mat4, m: &Mat4) -> f64 {
    let mut s = 0.0;
    for a in 0..4 {
        for b in 0..4 {
            s += ginv[a][b] * m[a][b];
        }
    }
    s
}

/// Einstein tensor `G_{ab} = R_{ab} − ½ R g_{ab}`.
fn einstein_tensor(ricci: &Mat4, ricci_scalar: f64, g: &Mat4) -> Mat4 {
    let mut e = [[0.0; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            e[a][b] = ricci[a][b] - 0.5 * ricci_scalar * g[a][b];
        }
    }
    e
}

/// Kretschmann scalar `K = R_{abcd} R^{abcd}`, computed by fully lowering and
/// fully raising the (mixed) Riemann tensor and contracting.
fn kretschmann(riem: &T4, g: &Mat4, ginv: &Mat4) -> f64 {
    // Fully lowered: R_{abcd} = g_{ae} R^e_{bcd}.
    let mut rd: T4 = [[[[0.0; 4]; 4]; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let mut s = 0.0;
                    for e in 0..4 {
                        s += g[a][e] * riem[e][b][c][d];
                    }
                    rd[a][b][c][d] = s;
                }
            }
        }
    }
    // Fully raised: R^{abcd} = g^{bp} g^{cq} g^{dw} R^a_{pqw}.
    let mut ru: T4 = [[[[0.0; 4]; 4]; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let mut s = 0.0;
                    for p in 0..4 {
                        for q in 0..4 {
                            for w in 0..4 {
                                s += ginv[b][p] * ginv[c][q] * ginv[d][w] * riem[a][p][q][w];
                            }
                        }
                    }
                    ru[a][b][c][d] = s;
                }
            }
        }
    }
    let mut k = 0.0;
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    k += rd[a][b][c][d] * ru[a][b][c][d];
                }
            }
        }
    }
    k
}

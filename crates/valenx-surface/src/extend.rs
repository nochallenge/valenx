//! Surface extrapolation — extend a surface beyond one of its
//! parametric edges by linearly extrapolating the boundary control
//! points along the surface's tangent there.
//!
//! v1 algorithm: for each control point on the chosen edge `e`,
//! compute the difference vector from the *adjacent* interior CP
//! and add a new CP at `edge_cp + length * unit_diff`. Continuity
//! at the original edge is preserved (the original CPs become the
//! interior of the extended surface).

use nalgebra::Vector3;

use crate::nurbs_surface::NurbsSurface;
use crate::sew::Edge;

/// Extend `surface` beyond `edge` by `length` along the surface's
/// tangent at that edge.
///
/// The returned NurbsSurface has one extra row (or column) of CPs
/// in the relevant direction; its u (or v) knot vector is a
/// clamped open-uniform vector of the new length.
pub fn extrapolate(surface: &NurbsSurface, edge: Edge, length: f64) -> NurbsSurface {
    let nu = surface.nu();
    let nv = surface.nv();
    match edge {
        Edge::UMax => {
            let mut cps = surface.control_points.clone();
            let mut weights = surface.weights.clone();
            let mut new_row = Vec::with_capacity(nv);
            let mut new_w = Vec::with_capacity(nv);
            for j in 0..nv {
                let last = surface.control_points[nu - 1][j];
                let prev = surface.control_points[nu - 2][j];
                let diff = last - prev;
                let unit = if diff.norm() < 1e-30 {
                    Vector3::zeros()
                } else {
                    diff.normalize()
                };
                new_row.push(last + length * unit);
                new_w.push(surface.weights[nu - 1][j]);
            }
            cps.push(new_row);
            weights.push(new_w);
            let nu_new = cps.len();
            let u_knots = open_uniform_knots(nu_new, surface.u_degree);
            NurbsSurface::new_unchecked(
                surface.u_degree,
                surface.v_degree,
                u_knots,
                surface.v_knots.clone(),
                cps,
                weights,
            )
        }
        Edge::UMin => {
            let mut cps = vec![Vec::new()];
            let mut weights = vec![Vec::new()];
            for j in 0..nv {
                let first = surface.control_points[0][j];
                let second = surface.control_points[1][j];
                let diff = first - second;
                let unit = if diff.norm() < 1e-30 {
                    Vector3::zeros()
                } else {
                    diff.normalize()
                };
                cps[0].push(first + length * unit);
                weights[0].push(surface.weights[0][j]);
            }
            for i in 0..nu {
                cps.push(surface.control_points[i].clone());
                weights.push(surface.weights[i].clone());
            }
            let nu_new = cps.len();
            let u_knots = open_uniform_knots(nu_new, surface.u_degree);
            NurbsSurface::new_unchecked(
                surface.u_degree,
                surface.v_degree,
                u_knots,
                surface.v_knots.clone(),
                cps,
                weights,
            )
        }
        Edge::VMax => {
            let mut cps = surface.control_points.clone();
            let mut weights = surface.weights.clone();
            for i in 0..nu {
                let last = surface.control_points[i][nv - 1];
                let prev = surface.control_points[i][nv - 2];
                let diff = last - prev;
                let unit = if diff.norm() < 1e-30 {
                    Vector3::zeros()
                } else {
                    diff.normalize()
                };
                cps[i].push(last + length * unit);
                weights[i].push(surface.weights[i][nv - 1]);
            }
            let nv_new = cps[0].len();
            let v_knots = open_uniform_knots(nv_new, surface.v_degree);
            NurbsSurface::new_unchecked(
                surface.u_degree,
                surface.v_degree,
                surface.u_knots.clone(),
                v_knots,
                cps,
                weights,
            )
        }
        Edge::VMin => {
            let mut cps = Vec::with_capacity(nu);
            let mut weights = Vec::with_capacity(nu);
            for i in 0..nu {
                let first = surface.control_points[i][0];
                let second = surface.control_points[i][1];
                let diff = first - second;
                let unit = if diff.norm() < 1e-30 {
                    Vector3::zeros()
                } else {
                    diff.normalize()
                };
                let mut row = Vec::with_capacity(nv + 1);
                let mut wrow = Vec::with_capacity(nv + 1);
                row.push(first + length * unit);
                wrow.push(surface.weights[i][0]);
                for j in 0..nv {
                    row.push(surface.control_points[i][j]);
                    wrow.push(surface.weights[i][j]);
                }
                cps.push(row);
                weights.push(wrow);
            }
            let nv_new = cps[0].len();
            let v_knots = open_uniform_knots(nv_new, surface.v_degree);
            NurbsSurface::new_unchecked(
                surface.u_degree,
                surface.v_degree,
                surface.u_knots.clone(),
                v_knots,
                cps,
                weights,
            )
        }
    }
}

fn open_uniform_knots(n_cp: usize, degree: usize) -> Vec<f64> {
    let p = degree;
    let m = n_cp + p + 1;
    let mut k = vec![0.0; m];
    if n_cp <= p + 1 {
        // Pure Bezier — `vec![0.0; m]` already zeroes the first
        // `p + 1` entries; set the trailing `p + 1` entries to 1.0.
        for kv in k.iter_mut().skip(m - p - 1) {
            *kv = 1.0;
        }
        return k;
    }
    let n_internal = n_cp - p - 1;
    for (i, kv) in k.iter_mut().enumerate().take(m) {
        if i <= p {
            *kv = 0.0;
        } else if i >= n_cp {
            *kv = 1.0;
        } else {
            let idx = i - p;
            *kv = idx as f64 / (n_internal + 1) as f64;
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    fn planar_unit_square() -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(u, v, 0.0)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    #[test]
    fn extend_umax_grows_in_x() {
        let s = planar_unit_square();
        let ext = extrapolate(&s, Edge::UMax, 1.0);
        // 5 rows of CPs now in u.
        assert_eq!(ext.nu(), 5);
        assert_eq!(ext.nv(), 4);
        // The new row (last) should have x = 1 + length * unit_diff.
        // Original last - second-last CP delta in x is 1/3; unit
        // length is 1; so new x = 1 + 1*1 = 2.
        for j in 0..4 {
            let cp = ext.control_points[4][j];
            assert!((cp.x - 2.0).abs() < 1e-10, "j={j}: x={}", cp.x);
            let v_expected = j as f64 / 3.0;
            assert!((cp.y - v_expected).abs() < 1e-10);
            assert!(cp.z.abs() < 1e-10);
        }
    }

    #[test]
    fn extend_vmin_grows_in_negative_y() {
        let s = planar_unit_square();
        let ext = extrapolate(&s, Edge::VMin, 0.5);
        // 5 columns of CPs now.
        assert_eq!(ext.nv(), 5);
        assert_eq!(ext.nu(), 4);
        // The new first column should have y = 0 + 0.5 * (-1) = -0.5.
        for i in 0..4 {
            let cp = ext.control_points[i][0];
            assert!((cp.y - (-0.5)).abs() < 1e-10, "i={i}: y={}", cp.y);
        }
    }

    #[test]
    fn extend_umax_preserves_original_interior() {
        let s = planar_unit_square();
        let ext = extrapolate(&s, Edge::UMax, 1.0);
        // First 4 rows of the extended surface must equal the
        // original CPs.
        for i in 0..4 {
            for j in 0..4 {
                assert!((ext.control_points[i][j] - s.control_points[i][j]).norm() < 1e-10);
            }
        }
    }
}

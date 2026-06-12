//! Phase 188 — `Prs3d_Drawer::SetUIsoAspect` /
//! `SetVIsoAspect` — display isoparametric lines on surfaces.
//!
//! ## What OCCT does
//!
//! Each NURBS / B-Rep face carries a U-parameter range and a V-parameter
//! range. OCCT's isoline display walks N equally-spaced U values,
//! evaluates the surface at `(u_i, v)` for `v ∈ [v_min, v_max]` (and
//! vice versa for V isolines), then renders the resulting curves as
//! thin grey overlays on the face. Standard count is 10 isolines per
//! direction (`Prs3d_Drawer::SetUIsoAspect(Prs3d_LineAspect, 10)`);
//! the spacing reveals surface curvature visually.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 188.5). The function takes a
//! [`valenx_surface::NurbsSurface`] and samples its isoparametric
//! curves directly: `u_count` lines of constant `u` (each sampled at
//! `samples_per_line` points across the surface's `v_range`) and
//! `v_count` lines of constant `v`. Every sample point comes from the
//! real `NurbsSurface::evaluate` tensor-product evaluator, so the
//! isolines follow the true surface curvature — exactly the overlay
//! geometry OCCT's `SetUIsoAspect` draws. The renderer consumes the
//! returned polylines as a thin-line wire overlay.

use valenx_surface::NurbsSurface;

use crate::error::OcctVizError;

/// Minimum isoline count per direction.
pub const MIN_ISOS: u32 = 1;
/// Maximum isoline count per direction.
pub const MAX_ISOS: u32 = 50;
/// Minimum sample points along each isoline.
pub const MIN_SAMPLES: u32 = 2;
/// Maximum sample points along each isoline.
pub const MAX_SAMPLES: u32 = 1000;

/// Isoparametric line overlay sampled off a surface.
#[derive(Clone, Debug, Default)]
pub struct Isolines {
    /// Lines of constant `u` — each inner `Vec` is one polyline
    /// running along the `v` direction.
    pub u_isolines: Vec<Vec<[f64; 3]>>,
    /// Lines of constant `v` — each inner `Vec` is one polyline
    /// running along the `u` direction.
    pub v_isolines: Vec<Vec<[f64; 3]>>,
}

/// Sample the isoparametric line overlay for `surface`.
///
/// `u_count` / `v_count` are the number of isolines in each
/// direction (`[MIN_ISOS, MAX_ISOS]`). `samples_per_line` is the
/// number of points each isoline is sampled into
/// (`[MIN_SAMPLES, MAX_SAMPLES]`).
///
/// Isolines are placed at the *interior* of the parameter range when
/// `count == 1` (the single line sits at the mid-parameter); for
/// `count >= 2` the lines span the range endpoints inclusive.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if any count or `samples_per_line`
///   falls outside its valid range.
///
/// # Example
///
/// ```
/// # use valenx_surface::NurbsSurface;
/// # use nalgebra::Vector3;
/// use valenx_occt_viz::prs3d_drawer_isolines::prs3d_drawer_isolines;
/// // A flat bilinear (degree-1) patch.
/// let cps = vec![
///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
///     vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
/// ];
/// let w = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
/// let knots = vec![0.0, 0.0, 1.0, 1.0];
/// let surf = NurbsSurface::new(1, 1, knots.clone(), knots, cps, w).unwrap();
/// let iso = prs3d_drawer_isolines(&surf, 5, 5, 12).unwrap();
/// assert_eq!(iso.u_isolines.len(), 5);
/// assert_eq!(iso.u_isolines[0].len(), 12);
/// ```
pub fn prs3d_drawer_isolines(
    surface: &NurbsSurface,
    u_count: u32,
    v_count: u32,
    samples_per_line: u32,
) -> Result<Isolines, OcctVizError> {
    if !(MIN_ISOS..=MAX_ISOS).contains(&u_count) {
        return Err(OcctVizError::bad_input(
            "u_count",
            format!("must be in [{MIN_ISOS}, {MAX_ISOS}] (got {u_count})"),
        ));
    }
    if !(MIN_ISOS..=MAX_ISOS).contains(&v_count) {
        return Err(OcctVizError::bad_input(
            "v_count",
            format!("must be in [{MIN_ISOS}, {MAX_ISOS}] (got {v_count})"),
        ));
    }
    if !(MIN_SAMPLES..=MAX_SAMPLES).contains(&samples_per_line) {
        return Err(OcctVizError::bad_input(
            "samples_per_line",
            format!("must be in [{MIN_SAMPLES}, {MAX_SAMPLES}] (got {samples_per_line})"),
        ));
    }

    let (u_lo, u_hi) = surface.u_range();
    let (v_lo, v_hi) = surface.v_range();
    let n = samples_per_line as usize;

    // Lines of constant u: walk u positions, sample along v.
    let mut u_isolines = Vec::with_capacity(u_count as usize);
    for ui in 0..u_count {
        let u = param_at(ui, u_count, u_lo, u_hi);
        let mut line = Vec::with_capacity(n);
        for s in 0..n {
            let t = s as f64 / (n - 1) as f64;
            let v = v_lo + (v_hi - v_lo) * t;
            let p = surface.evaluate(u, v);
            line.push([p.x, p.y, p.z]);
        }
        u_isolines.push(line);
    }

    // Lines of constant v: walk v positions, sample along u.
    let mut v_isolines = Vec::with_capacity(v_count as usize);
    for vi in 0..v_count {
        let v = param_at(vi, v_count, v_lo, v_hi);
        let mut line = Vec::with_capacity(n);
        for s in 0..n {
            let t = s as f64 / (n - 1) as f64;
            let u = u_lo + (u_hi - u_lo) * t;
            let p = surface.evaluate(u, v);
            line.push([p.x, p.y, p.z]);
        }
        v_isolines.push(line);
    }

    Ok(Isolines {
        u_isolines,
        v_isolines,
    })
}

/// Parameter value for isoline `idx` of `count` across `[lo, hi]`.
/// A single isoline sits at the midpoint; multiple isolines span the
/// range endpoints inclusive.
fn param_at(idx: u32, count: u32, lo: f64, hi: f64) -> f64 {
    if count == 1 {
        0.5 * (lo + hi)
    } else {
        let t = idx as f64 / (count - 1) as f64;
        lo + (hi - lo) * t
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// A flat unit bilinear patch on the XY plane.
    fn flat_patch() -> NurbsSurface {
        let cps = vec![
            vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
            vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
        ];
        let w = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        NurbsSurface::new(1, 1, knots.clone(), knots, cps, w).unwrap()
    }

    #[test]
    fn rejects_zero_u() {
        let err = prs3d_drawer_isolines(&flat_patch(), 0, 10, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_high_v() {
        let err = prs3d_drawer_isolines(&flat_patch(), 10, 100, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_too_few_samples() {
        let err = prs3d_drawer_isolines(&flat_patch(), 10, 10, 1).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn samples_the_requested_counts() {
        let iso = prs3d_drawer_isolines(&flat_patch(), 7, 4, 20).unwrap();
        assert_eq!(iso.u_isolines.len(), 7);
        assert_eq!(iso.v_isolines.len(), 4);
        for line in &iso.u_isolines {
            assert_eq!(line.len(), 20);
        }
        for line in &iso.v_isolines {
            assert_eq!(line.len(), 20);
        }
    }

    #[test]
    fn isolines_lie_on_the_surface() {
        // On a flat z=0 patch every isoline sample must have z == 0
        // and stay within the unit square.
        let iso = prs3d_drawer_isolines(&flat_patch(), 10, 10, 10).unwrap();
        for line in iso.u_isolines.iter().chain(iso.v_isolines.iter()) {
            for p in line {
                assert!(p[2].abs() < 1e-9, "flat patch isoline left z=0");
                assert!((-1e-9..=1.0 + 1e-9).contains(&p[0]));
                assert!((-1e-9..=1.0 + 1e-9).contains(&p[1]));
            }
        }
    }

    #[test]
    fn single_isoline_sits_at_midparameter() {
        // u_count == 1 → the lone u-isoline is at u = 0.5, so on the
        // flat patch all its points have x == 0.5.
        let iso = prs3d_drawer_isolines(&flat_patch(), 1, 1, 8).unwrap();
        assert_eq!(iso.u_isolines.len(), 1);
        for p in &iso.u_isolines[0] {
            assert!((p[0] - 0.5).abs() < 1e-9, "single u-isoline not at u=0.5");
        }
    }
}

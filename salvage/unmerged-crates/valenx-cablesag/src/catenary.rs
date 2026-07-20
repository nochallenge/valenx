//! Catenary cable under its own weight (load uniform along the arc).
//!
//! When the load `w` is uniform per unit *arc* length — a bare cable or
//! chain hanging under gravity — the equilibrium shape is the catenary,
//! not the parabola. With the origin at the low point and the catenary
//! parameter `a = H / w` (`H` the horizontal tension):
//!
//! `y(x) = a (cosh(x / a) - 1)`.
//!
//! For level supports a span `L` apart the standard results are
//! (see Irvine 1981, or any statics text):
//!
//! `d = a (cosh(L / 2a) - 1)`        (mid-span sag)
//! `s = 2 a sinh(L / 2a)`            (total arc length)
//! `T(x) = H cosh(x / a) = w (a + y)` (tension; max at the supports)
//! `T_max = H cosh(L / 2a) = w (a + d)`.
//!
//! The catenary and the [`crate::parabola`] agree in the shallow limit
//! (`L / 2a -> 0`) and the catenary sags slightly more for the same
//! `w, L, H` because its weight is concentrated nearer the supports.

use crate::error::{positive, CableError};

/// Catenary parameter `a = H / w` (a length).
///
/// Sets the scale of the curve: large `a` (high tension or light cable)
/// gives a shallow, nearly straight cable; small `a` gives a deep sag.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::catenary::parameter_a;
/// // H = 1000 N, w = 10 N/m  ->  a = 100 m.
/// let a = parameter_a(10.0, 1000.0).unwrap();
/// assert!((a - 100.0).abs() < 1e-12);
/// ```
pub fn parameter_a(w: f64, h_tension: f64) -> Result<f64, CableError> {
    let w = positive("w", w)?;
    let h = positive("h_tension", h_tension)?;
    Ok(h / w)
}

/// Mid-span sag of a catenary, `d = a (cosh(L / 2a) - 1)`.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::catenary::sag;
/// // a = 100, L = 100  ->  d = 100 (cosh(0.5) - 1) ~ 12.7626 m.
/// let d = sag(10.0, 100.0, 1000.0).unwrap();
/// assert!((d - 12.762_596_5).abs() < 1e-6);
/// ```
pub fn sag(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let l = positive("span_l", span_l)?;
    let a = parameter_a(w, h_tension)?;
    Ok(a * ((l / (2.0 * a)).cosh() - 1.0))
}

/// Vertical offset `y` below the low point at horizontal station `x`.
///
/// `y(x) = a (cosh(x / a) - 1)`. Unlike the parabola this form is valid
/// for any finite `x`; the model leaves the caller to keep `x` on the
/// physical span.
///
/// # Errors
///
/// Returns [`CableError`] if `w`, `h_tension` are non-positive or `x` is
/// non-finite.
pub fn profile_y(w: f64, h_tension: f64, x: f64) -> Result<f64, CableError> {
    let a = parameter_a(w, h_tension)?;
    let x = crate::error::finite("x", x)?;
    Ok(a * ((x / a).cosh() - 1.0))
}

/// Total arc length between level supports, `s = 2 a sinh(L / 2a)`.
///
/// Always exceeds the span `L` (since `sinh t > t` for `t > 0`).
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::catenary::arc_length;
/// let s = arc_length(10.0, 100.0, 1000.0).unwrap();
/// assert!(s > 100.0);
/// ```
pub fn arc_length(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let l = positive("span_l", span_l)?;
    let a = parameter_a(w, h_tension)?;
    Ok(2.0 * a * (l / (2.0 * a)).sinh())
}

/// Maximum cable tension, at the supports: `T = H cosh(L / 2a)`.
///
/// Equivalently `T = w (a + d)` with `d` the catenary [`sag`]: the
/// tension at height `y` above the low point is `H + w y`, and the
/// supports sit at `y = d`. Always `>= H`, since `cosh >= 1`.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::catenary::max_tension;
/// // T = 1000 * cosh(0.5) ~ 1127.6 N, and T >= H = 1000.
/// let t = max_tension(10.0, 100.0, 1000.0).unwrap();
/// assert!((t - 1_127.625_965).abs() < 1e-3);
/// assert!(t >= 1000.0);
/// ```
pub fn max_tension(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let l = positive("span_l", span_l)?;
    let h = positive("h_tension", h_tension)?;
    let a = parameter_a(w, h_tension)?;
    Ok(h * (l / (2.0 * a)).cosh())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::CableError;

    const EPS: f64 = 1e-9;

    // Ground-truth worked example:
    //   w = 10 N/m, L = 100 m, H = 1000 N  ->  a = 100 m, L/2a = 0.5.
    //   cosh(0.5) = 1.1276259652063807
    //   sinh(0.5) = 0.5210953054937474
    //   d = 100 (cosh 0.5 - 1) = 12.762596520638076 m
    //   s = 2*100*sinh(0.5)    = 104.21906109874948 m
    //   T = 1000 * cosh(0.5)   = 1127.6259652063807 N
    const W: f64 = 10.0;
    const L: f64 = 100.0;
    const H: f64 = 1000.0;

    #[test]
    fn parameter_a_is_h_over_w() {
        let a = parameter_a(W, H).unwrap();
        assert!((a - 100.0).abs() < EPS, "got {a}");
    }

    #[test]
    fn sag_matches_cosh_ground_truth() {
        let d = sag(W, L, H).unwrap();
        assert!((d - 12.762_596_520_638_076).abs() < 1e-9, "got {d}");
    }

    #[test]
    fn arc_length_matches_sinh_ground_truth_and_exceeds_span() {
        let s = arc_length(W, L, H).unwrap();
        assert!((s - 104.219_061_098_749_48).abs() < 1e-9, "got {s}");
        assert!(s > L);
    }

    #[test]
    fn max_tension_matches_cosh_and_exceeds_h() {
        let t = max_tension(W, L, H).unwrap();
        assert!((t - 1_127.625_965_206_380_7).abs() < 1e-9, "got {t}");
        assert!(t >= H);
    }

    #[test]
    fn max_tension_equals_w_times_a_plus_sag() {
        // Independent identity: T = w (a + d).
        let t = max_tension(W, L, H).unwrap();
        let a = parameter_a(W, H).unwrap();
        let d = sag(W, L, H).unwrap();
        assert!((t - W * (a + d)).abs() < 1e-6, "got {t}");
    }

    #[test]
    fn profile_endpoints_and_symmetry() {
        // y(0) = 0, y(+/-L/2) = sag, even in x.
        let d = sag(W, L, H).unwrap();
        assert!(profile_y(W, H, 0.0).unwrap().abs() < EPS);
        let y_end = profile_y(W, H, L / 2.0).unwrap();
        assert!((y_end - d).abs() < 1e-9, "got {y_end}");
        let yp = profile_y(W, H, 23.0).unwrap();
        let yn = profile_y(W, H, -23.0).unwrap();
        assert!((yp - yn).abs() < EPS);
    }

    #[test]
    fn tension_ge_h_invariant_over_a_sweep() {
        for &w in &[0.5, 5.0, 50.0] {
            for &l in &[10.0, 200.0, 2000.0] {
                for &h in &[100.0, 1000.0, 1.0e5] {
                    let t = max_tension(w, l, h).unwrap();
                    assert!(t >= h, "T={t} < H={h}");
                }
            }
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(matches!(
            sag(0.0, L, H),
            Err(CableError::NonPositive { name: "w", .. })
        ));
        assert!(matches!(
            parameter_a(W, -5.0),
            Err(CableError::NonPositive {
                name: "h_tension",
                ..
            })
        ));
        assert!(matches!(
            profile_y(W, H, f64::INFINITY),
            Err(CableError::NotFinite { name: "x", .. })
        ));
    }
}

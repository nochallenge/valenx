//! Parabolic cable under a uniform load distributed along the horizontal.
//!
//! When the load `w` (force per unit *horizontal* length) is constant —
//! the standard idealisation for a suspension-bridge deck hung from a
//! cable, or a light cable with a heavy uniform deck — the equilibrium
//! shape is an exact parabola. With level supports a horizontal span `L`
//! apart and horizontal tension component `H`, taking the origin at the
//! low point:
//!
//! `y(x) = w x^2 / (2 H)`,  for  `x` in `[-L/2, L/2]`.
//!
//! The headline results (all derived from cable statics, see Irvine,
//! *Cable Structures*, MIT Press 1981, ch. 1):
//!
//! `d   = w L^2 / (8 H)`          (mid-span sag)
//! `T   = sqrt(H^2 + (w L / 2)^2)`  (tension at a support, the maximum)
//! `V   = w L / 2`                (vertical reaction at each support)
//!
//! Because the support tension adds the vertical reaction in quadrature
//! with `H`, it always satisfies `T >= H`, with equality only in the
//! degenerate `w -> 0` limit.

use crate::error::{in_range, positive, CableError};

/// Mid-span sag of a parabolic cable, `d = w L^2 / (8 H)`.
///
/// # Parameters
///
/// `w` — uniform load per unit horizontal length (force/length, > 0).
/// `span_l` — horizontal distance between the level supports (> 0).
/// `h_tension` — horizontal tension component (force, > 0).
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::parabola::sag;
/// // w = 2 N/m, L = 100 m, H = 2500 N  ->  d = 2*100^2/(8*2500) = 1.0 m
/// let d = sag(2.0, 100.0, 2500.0).unwrap();
/// assert!((d - 1.0).abs() < 1e-12);
/// ```
pub fn sag(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let w = positive("w", w)?;
    let l = positive("span_l", span_l)?;
    let h = positive("h_tension", h_tension)?;
    Ok(w * l * l / (8.0 * h))
}

/// Horizontal tension required to hold a target mid-span `sag`.
///
/// Inverse of [`sag`]: `H = w L^2 / (8 d)`. Handy for design, where the
/// allowable sag is the given and the cable pretension is the unknown.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::parabola::horizontal_tension_for_sag;
/// let h = horizontal_tension_for_sag(2.0, 100.0, 1.0).unwrap();
/// assert!((h - 2500.0).abs() < 1e-9);
/// ```
pub fn horizontal_tension_for_sag(w: f64, span_l: f64, target_sag: f64) -> Result<f64, CableError> {
    let w = positive("w", w)?;
    let l = positive("span_l", span_l)?;
    let d = positive("target_sag", target_sag)?;
    Ok(w * l * l / (8.0 * d))
}

/// Vertical reaction at each support, `V = w L / 2`.
///
/// For level supports the load splits evenly, so each end carries half
/// the total weight `w L`.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
pub fn support_vertical_reaction(w: f64, span_l: f64) -> Result<f64, CableError> {
    let w = positive("w", w)?;
    let l = positive("span_l", span_l)?;
    Ok(w * l / 2.0)
}

/// Maximum cable tension, which occurs at the supports.
///
/// `T = sqrt(H^2 + (w L / 2)^2)`. The horizontal component is constant
/// along the cable at `H`; the vertical component grows from zero at the
/// low point to `w L / 2` at each support, so the magnitude is greatest
/// there. The result satisfies `T >= H` for all admissible inputs.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::parabola::max_tension;
/// // H = 2500, w L / 2 = 100  ->  T = sqrt(2500^2 + 100^2)
/// let t = max_tension(2.0, 100.0, 2500.0).unwrap();
/// assert!((t - (2500.0f64.powi(2) + 100.0f64.powi(2)).sqrt()).abs() < 1e-9);
/// assert!(t >= 2500.0);
/// ```
pub fn max_tension(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let h = positive("h_tension", h_tension)?;
    let v = support_vertical_reaction(w, span_l)?;
    Ok((h * h + v * v).sqrt())
}

/// Vertical offset `y` below the low point at horizontal station `x`.
///
/// `y(x) = w x^2 / (2 H)`, with `x` measured from the mid-span low point
/// and constrained to the half-span `[-L/2, L/2]`. At `x = L/2` this
/// equals the mid-span [`sag`].
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive, or
/// if `x` lies outside `[-L/2, L/2]`.
///
/// # Example
///
/// ```
/// use valenx_cablesag::parabola::{profile_y, sag};
/// let d = sag(2.0, 100.0, 2500.0).unwrap();
/// // At the support x = L/2 the profile returns the full sag.
/// let y_end = profile_y(2.0, 100.0, 2500.0, 50.0).unwrap();
/// assert!((y_end - d).abs() < 1e-12);
/// // At the low point the profile is zero.
/// assert!(profile_y(2.0, 100.0, 2500.0, 0.0).unwrap().abs() < 1e-12);
/// ```
pub fn profile_y(w: f64, span_l: f64, h_tension: f64, x: f64) -> Result<f64, CableError> {
    let w = positive("w", w)?;
    let l = positive("span_l", span_l)?;
    let h = positive("h_tension", h_tension)?;
    let half = l / 2.0;
    let x = in_range("x", x, -half, half, "x must lie in [-L/2, L/2]")?;
    Ok(w * x * x / (2.0 * h))
}

/// Arc length of the parabolic cable between the two level supports.
///
/// Uses the exact closed-form integral of `sqrt(1 + (y')^2)` for a
/// parabola. With `c = 4 d / L` (so `y' = c` at the support) and
/// `u = c` the half-span slope:
///
/// `s = (L / 2) * [ sqrt(1 + u^2) + asinh(u) / u ]`.
///
/// In the shallow-cable limit this reduces to the familiar series
/// `s ~ L (1 + 8 d^2 / (3 L^2))`, which [`arc_length_shallow`] returns.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::parabola::{arc_length, sag};
/// // Arc length always exceeds the straight-line span.
/// let s = arc_length(2.0, 100.0, 2500.0).unwrap();
/// assert!(s > 100.0);
/// ```
pub fn arc_length(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let l = positive("span_l", span_l)?;
    let d = sag(w, span_l, h_tension)?;
    // Slope at the support: y'(L/2) = w (L/2) / H = 4 d / L.
    let u = 4.0 * d / l;
    // s = (L/2) [ sqrt(1+u^2) + asinh(u)/u ].
    Ok((l / 2.0) * ((1.0 + u * u).sqrt() + u.asinh() / u))
}

/// Shallow-cable (small-sag) series approximation of the arc length.
///
/// `s ~ L (1 + 8 d^2 / (3 L^2))`, the leading correction to the span.
/// Valid when the sag-to-span ratio `d / L` is small (the usual
/// engineering regime, say below ~1/8). The dropped next term
/// `-32 d^4 / (5 L^3)` is negative, so this two-term series slightly
/// overestimates the exact [`arc_length`], and the gap grows with
/// `d / L`.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
pub fn arc_length_shallow(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let l = positive("span_l", span_l)?;
    let d = sag(w, span_l, h_tension)?;
    Ok(l * (1.0 + 8.0 * d * d / (3.0 * l * l)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::CableError;

    const EPS: f64 = 1e-9;

    // Ground-truth worked example used across the suite:
    //   w = 2 N/m, L = 100 m, H = 2500 N.
    //   d = w L^2 / (8 H) = 2 * 10000 / 20000 = 1.0 m   (exact)
    //   V = w L / 2 = 100 N
    //   T = sqrt(2500^2 + 100^2) = sqrt(6260000) = 2502.0 N (approx)
    const W: f64 = 2.0;
    const L: f64 = 100.0;
    const H: f64 = 2500.0;

    #[test]
    fn sag_matches_hand_worked_value() {
        let d = sag(W, L, H).unwrap();
        assert!((d - 1.0).abs() < EPS, "got {d}");
    }

    #[test]
    fn horizontal_tension_inverts_sag() {
        let d = sag(W, L, H).unwrap();
        let h = horizontal_tension_for_sag(W, L, d).unwrap();
        assert!((h - H).abs() < 1e-6, "got {h}");
    }

    #[test]
    fn vertical_reaction_is_half_total_load() {
        let v = support_vertical_reaction(W, L).unwrap();
        assert!((v - 100.0).abs() < EPS, "got {v}");
        // Two reactions carry the full weight w*L.
        assert!((2.0 * v - W * L).abs() < EPS);
    }

    #[test]
    fn max_tension_matches_quadrature_and_exceeds_h() {
        let t = max_tension(W, L, H).unwrap();
        let expected = (H * H + 100.0_f64 * 100.0).sqrt();
        assert!((t - expected).abs() < 1e-6, "got {t}");
        // Ground-truth invariant: T >= H always.
        assert!(t >= H);
    }

    #[test]
    fn max_tension_decomposes_into_h_and_v() {
        // T^2 should equal H^2 + V^2 exactly (Pythagoras on the
        // constant horizontal and the support vertical component).
        let t = max_tension(W, L, H).unwrap();
        let v = support_vertical_reaction(W, L).unwrap();
        assert!((t * t - (H * H + v * v)).abs() < 1e-3);
    }

    #[test]
    fn tension_ge_h_invariant_over_a_sweep() {
        // T >= H must hold for every admissible (w, L, H).
        for &w in &[0.1, 1.0, 5.0, 50.0] {
            for &l in &[10.0, 100.0, 1000.0] {
                for &h in &[100.0, 2500.0, 1.0e6] {
                    let t = max_tension(w, l, h).unwrap();
                    assert!(t >= h, "T={t} < H={h} for w={w} L={l}");
                }
            }
        }
    }

    #[test]
    fn profile_endpoints_and_symmetry() {
        // y(0) = 0, y(+/-L/2) = sag, and y is even in x.
        let d = sag(W, L, H).unwrap();
        assert!(profile_y(W, L, H, 0.0).unwrap().abs() < EPS);
        let y_end = profile_y(W, L, H, L / 2.0).unwrap();
        assert!((y_end - d).abs() < EPS, "got {y_end}");
        let yp = profile_y(W, L, H, 17.0).unwrap();
        let yn = profile_y(W, L, H, -17.0).unwrap();
        assert!((yp - yn).abs() < EPS);
    }

    #[test]
    fn profile_quarter_point_is_a_quarter_of_sag() {
        // y(L/4) / y(L/2) = (L/4)^2 / (L/2)^2 = 1/4 for a parabola.
        let y_q = profile_y(W, L, H, L / 4.0).unwrap();
        let y_e = profile_y(W, L, H, L / 2.0).unwrap();
        assert!((y_q - 0.25 * y_e).abs() < EPS, "got {y_q}");
    }

    #[test]
    fn arc_length_exceeds_span_and_beats_chord_bound() {
        // The arc is longer than the chord but, for this shallow case,
        // less than going straight down and across (L + 2d).
        let s = arc_length(W, L, H).unwrap();
        assert!(s > L, "arc {s} <= span {L}");
        assert!(s < L + 2.0 * 1.0);
    }

    #[test]
    fn arc_length_known_series_value() {
        // For d/L = 1/100, the shallow series gives
        //   s = L (1 + 8 (1/100)^2 / 3) = 100 * (1 + 8/30000)
        //     = 100.0266666...  m.
        let s = arc_length_shallow(W, L, H).unwrap();
        let expected = L * (1.0 + 8.0 * (1.0 / 100.0_f64).powi(2) / 3.0);
        assert!((s - expected).abs() < EPS, "got {s}");
        assert!((s - 100.026_666_666_666_67).abs() < 1e-6, "got {s}");
    }

    #[test]
    fn shallow_series_approximates_exact_for_small_sag() {
        // With d/L = 1/100 the two arc lengths agree to ~1e-6 relative.
        let s_exact = arc_length(W, L, H).unwrap();
        let s_approx = arc_length_shallow(W, L, H).unwrap();
        assert!((s_exact - s_approx).abs() < 1e-4, "{s_exact} vs {s_approx}");
        // The two-term series drops a negative quartic term, so it
        // slightly overestimates the exact arc length.
        assert!(s_approx >= s_exact - EPS, "{s_approx} vs {s_exact}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(matches!(
            sag(0.0, L, H),
            Err(CableError::NonPositive { name: "w", .. })
        ));
        assert!(matches!(
            sag(W, -1.0, H),
            Err(CableError::NonPositive { name: "span_l", .. })
        ));
        assert!(matches!(
            max_tension(W, L, f64::NAN),
            Err(CableError::NotFinite {
                name: "h_tension",
                ..
            })
        ));
        // x outside the half-span is rejected.
        assert!(matches!(
            profile_y(W, L, H, 60.0),
            Err(CableError::OutOfDomain { name: "x", .. })
        ));
    }
}

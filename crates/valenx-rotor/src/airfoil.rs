//! Airfoil section polars: lift `Cl(alpha)` and drag `Cd(alpha)`.
//!
//! Blade-element-momentum theory needs, at each radial station, the 2-D
//! section lift and drag coefficients as a function of the local angle of
//! attack `alpha` (radians). This module provides two polars:
//!
//! - [`Polar::Analytic`] — a closed-form thin-airfoil approximation,
//!   `Cl = 2*pi*sin(alpha)` clamped to a stall ceiling with a gentle linear
//!   post-stall slope, and a parabolic drag bucket `Cd = cd0 + k*Cl^2`. It
//!   needs no data and is handy for sanity checks and teaching, but it is
//!   only a rough stand-in for a real section.
//! - [`Polar::Table`] — a caller-supplied lookup table of `(alpha, Cl, Cd)`
//!   samples, linearly interpolated and clamped (held flat) outside the
//!   tabulated range. This is how measured / XFOIL airfoil data is fed in.
//!
//! ## Honest scope
//!
//! Both polars are 2-D, steady, incompressible and Reynolds-independent.
//! The analytic polar in particular is a thin-airfoil idealisation: real
//! airfoils have a finite lift-curve slope below `2*pi`, a camber offset,
//! a Reynolds- and Mach-dependent drag bucket, and a stall behaviour that
//! a single clamp does not capture. For anything beyond a rough estimate,
//! supply measured data via [`Polar::Table`].

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::error::RotorError;

/// Default thin-airfoil stall ceiling on `|Cl|` for the analytic polar.
///
/// Thin-airfoil theory (`Cl = 2*pi*sin(alpha)`) has no stall; real
/// sections lose lift past `~10-15 deg`. We clamp the linear region at
/// this `Cl` (`~ +/-1.5`, a typical max for a moderate-Reynolds section)
/// and apply a gentle post-stall slope beyond it.
pub const DEFAULT_CL_MAX: f64 = 1.5;

/// Default parabolic-drag offset `cd0` (minimum profile drag) for the
/// analytic polar.
pub const DEFAULT_CD0: f64 = 0.01;

/// Default induced-drag factor `k` in `Cd = cd0 + k*Cl^2` for the
/// analytic polar.
pub const DEFAULT_K: f64 = 0.02;

/// Post-stall lift slope, `dCl/dalpha` (per radian), applied beyond the
/// stall angle in the analytic polar. Small and negative-feeling in
/// practice (lift keeps rising only weakly); kept positive and shallow so
/// the polar stays single-valued and the BEMT residual stays well-behaved.
const POST_STALL_SLOPE: f64 = 0.5;

/// One sample of a tabulated airfoil polar.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PolarSample {
    /// Angle of attack `alpha` (radians).
    pub alpha_rad: f64,
    /// Section lift coefficient `Cl` at this `alpha`.
    pub cl: f64,
    /// Section drag coefficient `Cd` at this `alpha` (should be > 0).
    pub cd: f64,
}

/// A section lift/drag polar.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Polar {
    /// Closed-form thin-airfoil polar with a stall clamp and parabolic
    /// drag. See the module docs for the exact form.
    Analytic {
        /// Stall ceiling on `|Cl|` (the linear thin-airfoil region is
        /// clamped here). Must be finite and positive.
        cl_max: f64,
        /// Minimum profile drag `cd0`. Must be finite and positive.
        cd0: f64,
        /// Induced-drag factor `k` in `Cd = cd0 + k*Cl^2`. Must be finite
        /// and non-negative.
        k: f64,
    },
    /// Tabulated polar, linearly interpolated in `alpha`, held flat at the
    /// endpoints. Samples must be strictly increasing in `alpha`.
    Table(Vec<PolarSample>),
}

impl Default for Polar {
    /// The default polar is the analytic thin-airfoil polar with the
    /// crate's default `cl_max`, `cd0` and `k`.
    fn default() -> Self {
        Polar::Analytic {
            cl_max: DEFAULT_CL_MAX,
            cd0: DEFAULT_CD0,
            k: DEFAULT_K,
        }
    }
}

impl Polar {
    /// Build the default analytic thin-airfoil polar.
    pub fn analytic_default() -> Self {
        Polar::default()
    }

    /// Build a validated analytic polar with explicit coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`RotorError::NonPositive`] if `cl_max` or `cd0` is not
    /// finite and positive, and [`RotorError::NotFinite`] if `k` is not
    /// finite or is negative.
    pub fn analytic(cl_max: f64, cd0: f64, k: f64) -> Result<Self, RotorError> {
        let cl_max = crate::error::require_positive("cl_max", cl_max)?;
        let cd0 = crate::error::require_positive("cd0", cd0)?;
        if !(k.is_finite() && k >= 0.0) {
            return Err(RotorError::NotFinite {
                name: "k",
                value: k,
            });
        }
        Ok(Polar::Analytic { cl_max, cd0, k })
    }

    /// Build a validated tabulated polar from `(alpha, Cl, Cd)` samples.
    ///
    /// The samples must contain at least two entries, be strictly
    /// increasing in `alpha`, and all be finite. (Negative `Cd` is
    /// allowed structurally but discouraged; the BEMT integration treats
    /// `Cd` as supplied.)
    ///
    /// # Errors
    ///
    /// Returns [`RotorError::BadPolar`] if the table is too short, has a
    /// non-finite entry, or is not strictly increasing in `alpha`.
    pub fn table(samples: Vec<PolarSample>) -> Result<Self, RotorError> {
        if samples.len() < 2 {
            return Err(RotorError::BadPolar {
                reason: format!("need at least 2 samples, got {}", samples.len()),
            });
        }
        for (i, s) in samples.iter().enumerate() {
            if !(s.alpha_rad.is_finite() && s.cl.is_finite() && s.cd.is_finite()) {
                return Err(RotorError::BadPolar {
                    reason: format!("sample {i} has a non-finite value"),
                });
            }
            if i > 0 && s.alpha_rad <= samples[i - 1].alpha_rad {
                return Err(RotorError::BadPolar {
                    reason: format!(
                        "alpha not strictly increasing at sample {i} \
                         ({} <= {})",
                        s.alpha_rad,
                        samples[i - 1].alpha_rad
                    ),
                });
            }
        }
        Ok(Polar::Table(samples))
    }

    /// Lift coefficient `Cl` at angle of attack `alpha` (radians).
    ///
    /// Returns a finite value for any finite `alpha`.
    pub fn cl(&self, alpha: f64) -> f64 {
        match self {
            Polar::Analytic { cl_max, .. } => analytic_cl(alpha, *cl_max),
            Polar::Table(samples) => interp(samples, alpha, |s| s.cl),
        }
    }

    /// Drag coefficient `Cd` at angle of attack `alpha` (radians).
    ///
    /// Returns a finite value for any finite `alpha`.
    pub fn cd(&self, alpha: f64) -> f64 {
        match self {
            Polar::Analytic { cl_max, cd0, k } => {
                let cl = analytic_cl(alpha, *cl_max);
                cd0 + k * cl * cl
            }
            Polar::Table(samples) => interp(samples, alpha, |s| s.cd),
        }
    }

    /// Both coefficients at once: `(Cl, Cd)`.
    pub fn coefficients(&self, alpha: f64) -> (f64, f64) {
        (self.cl(alpha), self.cd(alpha))
    }
}

/// Analytic thin-airfoil lift: `2*pi*sin(alpha)` in the linear region,
/// clamped at `+/-cl_max` with a shallow linear post-stall slope so the
/// polar stays single-valued and bounded for any finite `alpha`.
fn analytic_cl(alpha: f64, cl_max: f64) -> f64 {
    // Linear (attached-flow) thin-airfoil lift.
    let linear = 2.0 * PI * alpha.sin();
    if linear.abs() <= cl_max {
        return linear;
    }
    // Past stall: the angle at which the linear branch first reached
    // +/-cl_max, then a shallow linear decay/rise beyond it. Using sin
    // keeps this well-defined for all alpha.
    let sign = if linear >= 0.0 { 1.0 } else { -1.0 };
    // Angle (in the principal branch) where |2*pi*sin(a)| == cl_max.
    let stall_sin = (cl_max / (2.0 * PI)).min(1.0);
    let stall_alpha = stall_sin.asin(); // in [0, pi/2]
    let excess = alpha.abs() - stall_alpha;
    sign * (cl_max + POST_STALL_SLOPE * excess).max(0.0)
}

/// Linear interpolation of a tabulated field, clamped (held flat) outside
/// the tabulated `alpha` range. `samples` is assumed validated: length
/// >= 2 and strictly increasing in `alpha`.
fn interp(samples: &[PolarSample], alpha: f64, field: impl Fn(&PolarSample) -> f64) -> f64 {
    let first = &samples[0];
    let last = &samples[samples.len() - 1];
    if alpha <= first.alpha_rad {
        return field(first);
    }
    if alpha >= last.alpha_rad {
        return field(last);
    }
    // Find the bracketing pair. Linear scan is fine for the small tables
    // typical of airfoil polars (tens of points).
    for w in samples.windows(2) {
        let (lo, hi) = (&w[0], &w[1]);
        if alpha >= lo.alpha_rad && alpha <= hi.alpha_rad {
            let span = hi.alpha_rad - lo.alpha_rad;
            // span > 0 because the table is strictly increasing.
            let t = (alpha - lo.alpha_rad) / span;
            return field(lo) + t * (field(hi) - field(lo));
        }
    }
    // Unreachable for a validated, in-range alpha; fall back to the last
    // value rather than panic.
    field(last)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * b.abs().max(1.0)
    }

    #[test]
    fn analytic_is_linear_near_zero() {
        let p = Polar::analytic_default();
        // Small alpha: Cl ~ 2*pi*alpha (sin ~ alpha).
        let a = 0.05_f64;
        assert!(close(p.cl(a), 2.0 * PI * a.sin()));
        // Symmetric.
        assert!(close(p.cl(-a), -p.cl(a)));
    }

    #[test]
    fn analytic_cl_is_clamped_and_finite_past_stall() {
        let p = Polar::analytic(DEFAULT_CL_MAX, DEFAULT_CD0, DEFAULT_K).unwrap();
        // Deep stall and even past 90 deg stays finite and bounded near
        // the ceiling (the post-stall slope is shallow).
        for &a in &[0.5, 1.0, 1.4, 2.0, 3.0, -2.0] {
            let cl = p.cl(a);
            assert!(cl.is_finite());
            assert!(cl.abs() <= DEFAULT_CL_MAX + 1.0, "cl = {cl} at alpha = {a}");
        }
    }

    #[test]
    fn analytic_drag_bucket() {
        let p = Polar::analytic(DEFAULT_CL_MAX, 0.012, 0.03).unwrap();
        // At alpha = 0, Cl = 0 so Cd = cd0.
        assert!(close(p.cd(0.0), 0.012));
        // Drag rises with |alpha| (more lift -> more induced drag).
        assert!(p.cd(0.1) > p.cd(0.0));
    }

    #[test]
    fn table_interpolates_and_clamps() {
        let p = Polar::table(vec![
            PolarSample {
                alpha_rad: -0.1,
                cl: -0.6,
                cd: 0.02,
            },
            PolarSample {
                alpha_rad: 0.0,
                cl: 0.0,
                cd: 0.01,
            },
            PolarSample {
                alpha_rad: 0.1,
                cl: 0.6,
                cd: 0.02,
            },
        ])
        .unwrap();
        // Midpoint interpolation.
        assert!(close(p.cl(0.05), 0.3));
        assert!(close(p.cd(0.05), 0.015));
        // Clamp below / above the table.
        assert!(close(p.cl(-1.0), -0.6));
        assert!(close(p.cl(1.0), 0.6));
        assert!(close(p.cd(1.0), 0.02));
    }

    #[test]
    fn table_rejects_bad_input() {
        // Too few samples.
        assert!(Polar::table(vec![PolarSample {
            alpha_rad: 0.0,
            cl: 0.0,
            cd: 0.01
        }])
        .is_err());
        // Non-increasing alpha.
        assert!(Polar::table(vec![
            PolarSample {
                alpha_rad: 0.1,
                cl: 0.0,
                cd: 0.01
            },
            PolarSample {
                alpha_rad: 0.1,
                cl: 0.1,
                cd: 0.01
            },
        ])
        .is_err());
        // Non-finite entry.
        assert!(Polar::table(vec![
            PolarSample {
                alpha_rad: 0.0,
                cl: f64::NAN,
                cd: 0.01
            },
            PolarSample {
                alpha_rad: 0.1,
                cl: 0.1,
                cd: 0.01
            },
        ])
        .is_err());
    }

    #[test]
    fn analytic_rejects_bad_coeffs() {
        assert!(Polar::analytic(-1.0, 0.01, 0.02).is_err());
        assert!(Polar::analytic(1.5, 0.0, 0.02).is_err());
        assert!(Polar::analytic(1.5, 0.01, -0.1).is_err());
        assert!(Polar::analytic(1.5, 0.01, f64::NAN).is_err());
    }
}

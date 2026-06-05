//! The angle-of-attack sweep — building a lift / drag polar.
//!
//! A single wind-tunnel run gives the coefficients at one angle of
//! attack. For a wing or an aircraft the engineering question is how
//! the coefficients vary *across* the angle-of-attack range — the
//! **drag polar** (`Cd` vs `Cl`) and the lift curve (`Cl` vs `α`).
//! That tells you the lift-curve slope, the zero-lift angle, the
//! stall angle and the maximum lift-to-drag ratio.
//!
//! This module runs the steady solver at a sequence of angles of
//! attack and assembles the [`PolarCurve`].
//!
//! # Honest scope
//!
//! A real sweep — each point is a genuine steady RANS solve. The two
//! v1 caveats are inherited from the solver: the immersed-boundary
//! Cartesian grid limits the absolute accuracy of each point, and a
//! steady RANS solver does not capture post-stall behaviour well (deep
//! stall is massively unsteady — a transient or scale-resolving run is
//! needed there). The pre-stall polar — the part an aircraft actually
//! cruises on — is the reliable range.

use crate::api::{run_on_tunnel, AeroRequest};
use crate::domain::WindTunnel;
use crate::error::AeroError;
use crate::geometry::TriMesh;
use crate::wind::Wind;

/// One point of a lift / drag polar — the coefficients at one angle of
/// attack.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarPoint {
    /// The angle of attack (radians).
    pub alpha: f64,
    /// The drag coefficient at this angle.
    pub cd: f64,
    /// The lift coefficient at this angle.
    pub cl: f64,
    /// The pitch-moment coefficient at this angle.
    pub cm: f64,
    /// `true` if the solve at this angle converged.
    pub converged: bool,
}

impl PolarPoint {
    /// The lift-to-drag ratio `L/D` — the efficiency of the body at
    /// this angle. Returns `0` if the drag is non-positive.
    pub fn lift_to_drag(&self) -> f64 {
        if self.cd > 1e-9 {
            self.cl / self.cd
        } else {
            0.0
        }
    }
}

/// A lift / drag polar — the coefficients across an angle-of-attack
/// sweep.
#[derive(Clone, Debug)]
pub struct PolarCurve {
    /// The polar points, ordered by ascending angle of attack.
    pub points: Vec<PolarPoint>,
}

impl PolarCurve {
    /// The maximum lift-to-drag ratio over the sweep — the best
    /// cruise efficiency.
    pub fn max_lift_to_drag(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.lift_to_drag())
            .fold(f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }

    /// The angle of attack of maximum lift — the stall angle estimate
    /// (radians). Returns `None` for an empty curve.
    pub fn stall_angle(&self) -> Option<f64> {
        self.points
            .iter()
            .max_by(|a, b| a.cl.partial_cmp(&b.cl).unwrap_or(std::cmp::Ordering::Equal))
            .map(|p| p.alpha)
    }

    /// The maximum lift coefficient reached over the sweep.
    pub fn max_lift(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.cl)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// The lift-curve slope `dCl/dα` (per radian), estimated by a
    /// least-squares fit over the points whose angle is below the
    /// stall angle (the linear pre-stall range). Returns `0` if there
    /// are too few pre-stall points.
    pub fn lift_curve_slope(&self) -> f64 {
        let stall = self.stall_angle().unwrap_or(f64::INFINITY);
        let pre: Vec<&PolarPoint> = self
            .points
            .iter()
            .filter(|p| p.alpha < stall - 1e-9)
            .collect();
        if pre.len() < 2 {
            return 0.0;
        }
        // Least-squares slope of Cl vs alpha.
        let n = pre.len() as f64;
        let sx: f64 = pre.iter().map(|p| p.alpha).sum();
        let sy: f64 = pre.iter().map(|p| p.cl).sum();
        let sxx: f64 = pre.iter().map(|p| p.alpha * p.alpha).sum();
        let sxy: f64 = pre.iter().map(|p| p.alpha * p.cl).sum();
        let denom = n * sxx - sx * sx;
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (n * sxy - sx * sy) / denom
        }
    }

    /// The polar point of best lift-to-drag ratio — the most efficient
    /// (best-glide / cruise) angle of attack. Returns `None` for an empty
    /// curve. Drag-free points contribute `L/D = 0` (see
    /// [`PolarPoint::lift_to_drag`]), so they only win when every point is
    /// drag-free.
    pub fn best_lift_to_drag_point(&self) -> Option<PolarPoint> {
        self.points.iter().copied().max_by(|a, b| {
            a.lift_to_drag()
                .partial_cmp(&b.lift_to_drag())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// The polar point of least drag — the bottom of the drag bucket (the
    /// parasitic, near-zero-lift drag floor). Returns `None` for an empty curve.
    pub fn min_drag_point(&self) -> Option<PolarPoint> {
        self.points
            .iter()
            .copied()
            .min_by(|a, b| a.cd.partial_cmp(&b.cd).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// The zero-lift angle of attack (radians) — where the lift curve crosses
    /// `Cl = 0`, found by linear interpolation between the first pair of polar
    /// points that bracket zero lift (negative for a cambered section, ~0 for a
    /// symmetric one). `None` if the sweep never crosses zero lift.
    pub fn zero_lift_angle(&self) -> Option<f64> {
        for w in self.points.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            let brackets = (lo.cl <= 0.0 && hi.cl >= 0.0) || (lo.cl >= 0.0 && hi.cl <= 0.0);
            if brackets {
                let dcl = hi.cl - lo.cl;
                if dcl.abs() < 1e-12 {
                    return Some(lo.alpha);
                }
                let t = -lo.cl / dcl;
                return Some(lo.alpha + t * (hi.alpha - lo.alpha));
            }
        }
        None
    }

    /// The induced-drag factor `k` in the parabolic drag polar
    /// `Cd = Cd₀ + k·Cl²` — the least-squares slope of `Cd` against `Cl²` over
    /// the sweep. A larger `k` means drag climbs faster with lift (lower span
    /// efficiency). Returns `0` for fewer than two points or a degenerate
    /// (constant-`Cl²`) fit.
    pub fn induced_drag_factor(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        let n = self.points.len() as f64;
        let sx: f64 = self.points.iter().map(|p| p.cl * p.cl).sum();
        let sy: f64 = self.points.iter().map(|p| p.cd).sum();
        let sxx: f64 = self.points.iter().map(|p| (p.cl * p.cl).powi(2)).sum();
        let sxy: f64 = self.points.iter().map(|p| p.cl * p.cl * p.cd).sum();
        let denom = n * sxx - sx * sx;
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (n * sxy - sx * sy) / denom
        }
    }

    /// The parasitic (zero-lift) drag coefficient `Cd₀` — the intercept of the
    /// parabolic drag polar `Cd = Cd₀ + k·Cl²`, i.e. the modeled drag at zero
    /// lift. With [`PolarCurve::induced_drag_factor`] (`k`) it fully defines the
    /// fitted polar. Computed from the same least-squares fit as
    /// `mean(Cd) − k·mean(Cl²)`. Returns `0` for fewer than two points.
    pub fn parasitic_drag_coefficient(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        let k = self.induced_drag_factor();
        let n = self.points.len() as f64;
        let mean_cd: f64 = self.points.iter().map(|p| p.cd).sum::<f64>() / n;
        let mean_cl2: f64 = self.points.iter().map(|p| p.cl * p.cl).sum::<f64>() / n;
        mean_cd - k * mean_cl2
    }

    /// The longitudinal static-stability slope `dCm/dCl` — the least-squares
    /// slope of the pitching-moment coefficient against the lift coefficient
    /// over the sweep. A **negative** slope is statically stable (a lift
    /// increase yields a nose-down restoring moment); its magnitude about the
    /// centre of gravity is the static margin. Returns `0` for fewer than two
    /// points or a degenerate (constant-`Cl`) fit.
    pub fn pitch_stability_slope(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        let n = self.points.len() as f64;
        let sx: f64 = self.points.iter().map(|p| p.cl).sum();
        let sy: f64 = self.points.iter().map(|p| p.cm).sum();
        let sxx: f64 = self.points.iter().map(|p| p.cl * p.cl).sum();
        let sxy: f64 = self.points.iter().map(|p| p.cl * p.cm).sum();
        let denom = n * sxx - sx * sx;
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (n * sxy - sx * sy) / denom
        }
    }

    /// The polar point of maximum **endurance** — the angle of attack that
    /// maximises `Cl^1.5 / Cd`, i.e. minimises the power required for steady
    /// level flight (power ∝ `Cd / Cl^1.5`). This is the best-loiter /
    /// minimum-sink condition; for a parabolic polar it sits at `√3` times the
    /// lift of the best lift-to-drag (best-range) point, so it is a distinct,
    /// higher-`Cl` optimum. Only positive-lift, finite-drag points qualify;
    /// returns `None` if none do.
    pub fn best_endurance_point(&self) -> Option<PolarPoint> {
        self.points
            .iter()
            .copied()
            .filter(|p| p.cl > 0.0 && p.cd > 1e-9)
            .max_by(|a, b| {
                let ea = a.cl.powf(1.5) / a.cd;
                let eb = b.cl.powf(1.5) / b.cd;
                ea.partial_cmp(&eb).unwrap_or(std::cmp::Ordering::Equal)
            })
    }
}

/// Run an angle-of-attack sweep and assemble the lift / drag polar.
///
/// `angles` is the list of angles of attack (radians) to solve at;
/// `base` is the request template (its `pitch` is overridden per
/// angle). The tunnel is built **once** at the first angle's
/// orientation and reused for every angle by re-orienting the wind —
/// the body voxelization is identical across the sweep (it is the wind
/// that rotates, not the body), which is both correct and far cheaper
/// than rebuilding the grid per angle.
///
/// Returns an [`AeroError`] only if the case setup is invalid; an
/// individual angle that fails to converge is recorded with
/// `converged == false`, not aborted.
pub fn aoa_sweep(
    body: &TriMesh,
    base: &AeroRequest,
    angles: &[f64],
) -> Result<PolarCurve, AeroError> {
    if angles.is_empty() {
        return Err(AeroError::BadParameter {
            name: "angles",
            reason: "the sweep needs at least one angle of attack".into(),
        });
    }

    // Build the tunnel once. The wind direction nominally points
    // along +x; the angle of attack is applied as the wind's pitch,
    // so the same Cartesian grid + body voxelization serves every
    // angle.
    let first_wind = Wind::new(
        base.speed,
        base.yaw,
        angles[0],
        base.air,
        base.turbulence_intensity,
    )?;
    let tunnel =
        WindTunnel::build_with(body, first_wind, base.boundary, base.sizing)?;

    let mut points = Vec::with_capacity(angles.len());
    for &alpha in angles {
        // Re-orient the wind in a fresh tunnel struct (cheap clone —
        // the grid + voxelization are shared, only the wind changes).
        let mut t = tunnel.clone();
        let wind = Wind::new(
            base.speed,
            base.yaw,
            alpha,
            base.air,
            base.turbulence_intensity,
        )?;
        t.wind = wind;
        // The reference area follows the new wind direction.
        t.reference_area = body.frontal_area(wind.direction()).max(1e-9);

        let mut req = *base;
        req.pitch = alpha;
        let result = run_on_tunnel(&t, &req)?;
        points.push(PolarPoint {
            alpha,
            cd: result.coefficients.cd,
            cl: result.coefficients.cl,
            cm: result.coefficients.cmy,
            converged: result.converged,
        });
    }

    // Order by angle.
    points.sort_by(|a, b| a.alpha.partial_cmp(&b.alpha).unwrap_or(std::cmp::Ordering::Equal));
    Ok(PolarCurve { points })
}

/// Build an evenly-spaced list of angle-of-attack values (radians)
/// from `start` to `end` inclusive, with `count` points — a
/// convenience for [`aoa_sweep`].
pub fn linspace_degrees(start_deg: f64, end_deg: f64, count: usize) -> Vec<f64> {
    let count = count.max(1);
    let d2r = std::f64::consts::PI / 180.0;
    if count == 1 {
        return vec![start_deg * d2r];
    }
    (0..count)
        .map(|i| {
            let t = i as f64 / (count - 1) as f64;
            (start_deg + t * (end_deg - start_deg)) * d2r
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::box_body;
    use crate::turbulence::TurbulenceModel;
    use nalgebra::Vector3;

    #[test]
    fn linspace_degrees_spans_the_range() {
        let a = linspace_degrees(-4.0, 8.0, 4);
        assert_eq!(a.len(), 4);
        let d2r = std::f64::consts::PI / 180.0;
        assert!((a[0] - (-4.0 * d2r)).abs() < 1e-12);
        assert!((a[3] - (8.0 * d2r)).abs() < 1e-12);
    }

    #[test]
    fn polar_point_lift_to_drag() {
        let p = PolarPoint {
            alpha: 0.1,
            cd: 0.05,
            cl: 0.5,
            cm: 0.0,
            converged: true,
        };
        // L/D = 0.5/0.05 = 10.
        assert!((p.lift_to_drag() - 10.0).abs() < 1e-9);
        // Zero drag → zero L/D (guarded).
        let p0 = PolarPoint {
            alpha: 0.0,
            cd: 0.0,
            cl: 0.5,
            cm: 0.0,
            converged: true,
        };
        assert_eq!(p0.lift_to_drag(), 0.0);
    }

    #[test]
    fn polar_curve_diagnostics_on_a_synthetic_polar() {
        // A hand-built polar: lift rises linearly, peaks, then falls
        // (a stall). The diagnostics must find the peak.
        let pts = vec![
            PolarPoint { alpha: 0.0, cd: 0.02, cl: 0.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.1, cd: 0.03, cl: 0.5, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.2, cd: 0.05, cl: 1.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.3, cd: 0.12, cl: 1.2, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.4, cd: 0.25, cl: 0.8, cm: 0.0, converged: true },
        ];
        let curve = PolarCurve { points: pts };
        // Max lift is 1.2 at alpha = 0.3.
        assert!((curve.max_lift() - 1.2).abs() < 1e-9);
        assert!((curve.stall_angle().unwrap() - 0.3).abs() < 1e-9);
        // The lift-curve slope (pre-stall) is positive.
        assert!(curve.lift_curve_slope() > 0.0);
        // Max L/D is positive.
        assert!(curve.max_lift_to_drag() > 0.0);
    }

    #[test]
    fn best_lift_to_drag_point_is_the_polar_peak() {
        // L/D by point: 0, 16.7, 20.0 (peak), 10.0, 3.2.
        let pts = vec![
            PolarPoint { alpha: 0.0, cd: 0.02, cl: 0.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.1, cd: 0.03, cl: 0.5, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.2, cd: 0.05, cl: 1.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.3, cd: 0.12, cl: 1.2, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.4, cd: 0.25, cl: 0.8, cm: 0.0, converged: true },
        ];
        let curve = PolarCurve { points: pts };
        let best = curve.best_lift_to_drag_point().expect("non-empty curve");
        // Peak efficiency is at α = 0.2 rad (L/D = 1.0 / 0.05 = 20).
        assert!((best.alpha - 0.2).abs() < 1e-12, "best-glide α");
        assert!((best.lift_to_drag() - 20.0).abs() < 1e-9);
        // It agrees with the curve's reported maximum L/D value.
        assert!((best.lift_to_drag() - curve.max_lift_to_drag()).abs() < 1e-9);
        // An empty curve has no best point.
        assert!(PolarCurve { points: vec![] }.best_lift_to_drag_point().is_none());
    }

    #[test]
    fn min_drag_point_is_the_bottom_of_the_drag_bucket() {
        let pts = vec![
            PolarPoint { alpha: 0.0, cd: 0.02, cl: 0.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.1, cd: 0.03, cl: 0.5, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.2, cd: 0.05, cl: 1.0, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.3, cd: 0.12, cl: 1.2, cm: 0.0, converged: true },
            PolarPoint { alpha: 0.4, cd: 0.25, cl: 0.8, cm: 0.0, converged: true },
        ];
        let curve = PolarCurve { points: pts };
        let min = curve.min_drag_point().expect("non-empty curve");
        // Lowest drag is Cd = 0.02 at α = 0.
        assert!((min.cd - 0.02).abs() < 1e-12, "min drag cd {}", min.cd);
        assert!((min.alpha - 0.0).abs() < 1e-12, "min drag α");
        // No point has lower drag.
        assert!(curve.points.iter().all(|p| p.cd >= min.cd - 1e-12));
        // An empty curve has no minimum-drag point.
        assert!(PolarCurve { points: vec![] }.min_drag_point().is_none());
    }

    #[test]
    fn zero_lift_angle_interpolates_the_cl_crossing() {
        // A cambered lift curve crossing Cl = 0: Cl=−0.1 at α=−0.05, Cl=+0.3 at
        // α=+0.05. Linear interp t = 0.1/0.4 = 0.25 → α = −0.05 + 0.25·0.1.
        let curve = PolarCurve {
            points: vec![
                PolarPoint { alpha: -0.05, cd: 0.02, cl: -0.1, cm: 0.0, converged: true },
                PolarPoint { alpha: 0.05, cd: 0.03, cl: 0.3, cm: 0.0, converged: true },
            ],
        };
        let a0 = curve.zero_lift_angle().expect("crosses zero lift");
        assert!((a0 - (-0.025)).abs() < 1e-9, "zero-lift α {a0}");
        // A purely positive-lift sweep never crosses zero → None.
        let positive = PolarCurve {
            points: vec![
                PolarPoint { alpha: 0.0, cd: 0.02, cl: 0.5, cm: 0.0, converged: true },
                PolarPoint { alpha: 0.1, cd: 0.03, cl: 1.0, cm: 0.0, converged: true },
            ],
        };
        assert!(positive.zero_lift_angle().is_none());
    }

    #[test]
    fn induced_drag_factor_recovers_a_parabolic_polar() {
        // Cd = 0.02 + 0.05·Cl² exactly → the least-squares fit recovers k = 0.05.
        let pts: Vec<PolarPoint> = [0.0, 0.5, 1.0, 1.5]
            .iter()
            .map(|&cl| PolarPoint {
                alpha: 0.0,
                cd: 0.02 + 0.05 * cl * cl,
                cl,
                cm: 0.0,
                converged: true,
            })
            .collect();
        let curve = PolarCurve { points: pts };
        assert!(
            (curve.induced_drag_factor() - 0.05).abs() < 1e-9,
            "k = {}",
            curve.induced_drag_factor()
        );
        // Too few points → 0.
        assert_eq!(PolarCurve { points: vec![] }.induced_drag_factor(), 0.0);
    }

    #[test]
    fn parasitic_drag_coefficient_recovers_the_polar_intercept() {
        // The same exact polar Cd = 0.02 + 0.05·Cl² → intercept Cd₀ = 0.02.
        let pts: Vec<PolarPoint> = [-0.3, 0.0, 0.5, 1.0, 1.5]
            .iter()
            .map(|&cl| PolarPoint {
                alpha: 0.0,
                cd: 0.02 + 0.05 * cl * cl,
                cl,
                cm: 0.0,
                converged: true,
            })
            .collect();
        let curve = PolarCurve { points: pts };
        assert!(
            (curve.induced_drag_factor() - 0.05).abs() < 1e-9,
            "k = {}",
            curve.induced_drag_factor()
        );
        assert!(
            (curve.parasitic_drag_coefficient() - 0.02).abs() < 1e-9,
            "Cd0 = {}",
            curve.parasitic_drag_coefficient()
        );
        // Too few points → 0.
        assert_eq!(
            PolarCurve { points: vec![] }.parasitic_drag_coefficient(),
            0.0
        );
    }

    #[test]
    fn pitch_stability_slope_recovers_a_linear_moment_curve() {
        // A statically stable section: Cm = 0.05 − 0.12·Cl → slope dCm/dCl = −0.12.
        let pts: Vec<PolarPoint> = [-0.2, 0.2, 0.6, 1.0]
            .iter()
            .map(|&cl| PolarPoint {
                alpha: 0.0,
                cd: 0.02,
                cl,
                cm: 0.05 - 0.12 * cl,
                converged: true,
            })
            .collect();
        let curve = PolarCurve { points: pts };
        assert!(
            (curve.pitch_stability_slope() - (-0.12)).abs() < 1e-9,
            "dCm/dCl = {}",
            curve.pitch_stability_slope()
        );
        // A negative slope flags static stability.
        assert!(curve.pitch_stability_slope() < 0.0);
        // Too few points → 0.
        assert_eq!(PolarCurve { points: vec![] }.pitch_stability_slope(), 0.0);
    }

    #[test]
    fn best_endurance_point_maximises_cl_to_the_three_halves_over_cd() {
        // A parabolic polar Cd = 0.02 + 0.05·Cl² sampled over a lift sweep.
        let pts: Vec<PolarPoint> = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
            .iter()
            .map(|&cl| PolarPoint {
                alpha: cl, // a monotone stand-in for the angle of attack
                cd: 0.02 + 0.05 * cl * cl,
                cl,
                cm: 0.0,
                converged: true,
            })
            .collect();
        let curve = PolarCurve { points: pts.clone() };
        let endurance = curve.best_endurance_point().unwrap();
        // It is exactly the sample that maximises Cl^1.5 / Cd.
        let best = pts
            .iter()
            .max_by(|a, b| {
                (a.cl.powf(1.5) / a.cd)
                    .partial_cmp(&(b.cl.powf(1.5) / b.cd))
                    .unwrap()
            })
            .unwrap();
        assert!((endurance.cl - best.cl).abs() < 1e-12, "endurance cl {}", endurance.cl);
        // Max endurance sits at a higher lift than max range (best L/D).
        let range = curve.best_lift_to_drag_point().unwrap();
        assert!(
            endurance.cl > range.cl,
            "endurance cl {} should exceed best-range cl {}",
            endurance.cl,
            range.cl
        );
        // With no positive-lift points there is no endurance optimum.
        let neg = PolarCurve {
            points: vec![PolarPoint {
                alpha: 0.0,
                cd: 0.1,
                cl: -0.5,
                cm: 0.0,
                converged: true,
            }],
        };
        assert!(neg.best_endurance_point().is_none());
    }

    #[test]
    fn aoa_sweep_rejects_an_empty_angle_list() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 0.3, 1.0));
        let base = AeroRequest::new(20.0);
        assert!(aoa_sweep(&body, &base, &[]).is_err());
    }

    #[test]
    fn aoa_sweep_produces_one_point_per_angle() {
        // A short three-angle sweep over a flat-plate-ish body. The
        // tunnel is kept compact (small clearances) with a modest cell
        // cap so the three steady SIMPLE solves stay fast — but the cap
        // must still leave the cells fine enough to voxelize the thin
        // 0.2 m plate (a razor-thin body in a too-coarse grid voxelizes
        // to zero solid cells). This is a structural check (one point
        // per angle, ordered, finite).
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 0.2, 1.0));
        let base = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(crate::domain::TunnelSizing {
                upstream: 1.0,
                downstream: 2.0,
                lateral: 1.0,
                cells_across_body: 4,
                max_cells: 150_000,
            })
            .with_max_iterations(60);
        let angles = linspace_degrees(0.0, 6.0, 3);
        let curve = aoa_sweep(&body, &base, &angles).unwrap();
        assert_eq!(curve.points.len(), 3);
        // Points are ordered by angle and all coefficients finite.
        for w in curve.points.windows(2) {
            assert!(w[1].alpha >= w[0].alpha);
        }
        assert!(curve.points.iter().all(|p| p.cd.is_finite() && p.cl.is_finite()));
    }
}

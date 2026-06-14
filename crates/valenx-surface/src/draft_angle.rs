//! **Draft-angle (mould-release) analysis** for manufacturability — the
//! pull-direction check a CAD tool runs before a part can be injection-moulded
//! or die-cast.
//!
//! For a chosen **pull direction** `d` (the axis along which the mould halves
//! separate), the **draft angle** of a surface point is
//!
//! ```text
//!   draft = asin( n̂ · d̂ )      (degrees)
//! ```
//!
//! where `n̂` is the surface unit normal and `d̂` the unit pull direction. The
//! sign convention follows the surface's own normal orientation
//! ([`NurbsSurface::normal`](crate::nurbs_surface::NurbsSurface::normal) returns
//! `∂u × ∂v` normalised):
//!
//! - **0°** — the face is *parallel* to the pull (a vertical wall): zero draft,
//!   it will drag against the mould and needs draft added.
//! - **+90°** — the face squarely faces the pull (`n̂` along `d̂`): fully
//!   drafted, releases freely.
//! - **negative** — the face points *away* from the pull (`n̂ · d̂ < 0`): an
//!   **undercut**, trapped in this mould half.
//!
//! A part is *mouldable* in direction `d` when every face has draft at least a
//! required minimum (commonly 1–3°). [`draft_report`] samples the surface and
//! summarises min / max / mean draft, the undercut count, and the verdict.
//!
//! Validated analytically: a plane facing the pull → 90°; a plane tilted by `θ`
//! → `90° − θ`; a cylinder pulled along its **axis** → ~0° on the walls
//! (vertical, undrafted).
//!
//! Honest scope: a pointwise draft diagnostic from the surface normal —
//! research-grade. It is the per-point data and a summary, not a coloured
//! viewport map or a parting-line/parting-surface generator, and a step toward,
//! not an equal of, CATIA-class draft analysis.

use nalgebra::Vector3;

use crate::nurbs_surface::NurbsSurface;

/// The **draft angle** (degrees) of `surface` at `(u, v)` for unit-ised pull
/// direction `pull`: `asin(n̂ · d̂)`. `0` if the normal or pull is degenerate.
///
/// See the [module docs](self) for the sign convention (0° = wall parallel to
/// pull, +90° = faces the pull, negative = undercut).
pub fn draft_angle(surface: &NurbsSurface, u: f64, v: f64, pull: Vector3<f64>) -> f64 {
    let n = surface.normal(u, v);
    let pull_norm = pull.norm();
    if n.norm() < 1e-12 || pull_norm < 1e-12 {
        return 0.0;
    }
    let d = pull / pull_norm;
    n.dot(&d).clamp(-1.0, 1.0).asin().to_degrees()
}

/// Summary of a draft-angle analysis over a sampled `(u, v)` grid.
#[derive(Clone, Debug, PartialEq)]
pub struct DraftReport {
    /// Smallest (most negative / most undercut) sampled draft angle (degrees).
    pub min_deg: f64,
    /// Largest sampled draft angle (degrees).
    pub max_deg: f64,
    /// Mean sampled draft angle (degrees).
    pub mean_deg: f64,
    /// Number of `(u, v)` samples taken.
    pub samples: usize,
    /// Samples whose draft is below the required minimum (need more draft or
    /// are undercuts) — the faces that would fail to release.
    pub undercut_count: usize,
    /// Whether every sample meets the required minimum draft (mouldable in this
    /// pull direction).
    pub mouldable: bool,
}

/// Sample `surface` on an `n × n` `(u, v)` grid (`n = samples.max(2)`) and
/// summarise its draft angle for pull direction `pull`, flagging any sample
/// below `min_draft_deg` as failing to release.
pub fn draft_report(
    surface: &NurbsSurface,
    pull: Vector3<f64>,
    samples: usize,
    min_draft_deg: f64,
) -> DraftReport {
    let n = samples.max(2);
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let mut min_deg = f64::INFINITY;
    let mut max_deg = f64::NEG_INFINITY;
    let mut sum = 0.0;
    let mut count = 0usize;
    let mut undercut = 0usize;
    for i in 0..n {
        let s = i as f64 / (n - 1) as f64;
        let u = u0 + s * (u1 - u0);
        for j in 0..n {
            let t = j as f64 / (n - 1) as f64;
            let v = v0 + t * (v1 - v0);
            let d = draft_angle(surface, u, v, pull);
            min_deg = min_deg.min(d);
            max_deg = max_deg.max(d);
            sum += d;
            count += 1;
            if d < min_draft_deg {
                undercut += 1;
            }
        }
    }
    DraftReport {
        min_deg,
        max_deg,
        mean_deg: sum / count.max(1) as f64,
        samples: count,
        undercut_count: undercut,
        mouldable: undercut == 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The flat `z = 0` plane (normal `+ẑ`), spanning x in u and y in v.
    fn plane() -> NurbsSurface {
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
                vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .expect("valid plane")
    }

    /// A flat plane tilted by `theta` about the y-axis: `∂u = (cosθ, 0, sinθ)`,
    /// `∂v = (0, 1, 0)`, so the normal `n = ∂u × ∂v = (−sinθ, 0, cosθ)` makes
    /// angle `theta` with `+ẑ` (`n · ẑ = cosθ`).
    fn tilted_plane(theta: f64) -> NurbsSurface {
        let (c, s) = (theta.cos(), theta.sin());
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
                vec![Vector3::new(c, 0.0, s), Vector3::new(c, 1.0, s)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .expect("valid tilted plane")
    }

    /// A rational-quadratic quarter cylinder of radius `r`, arc in u, extruded
    /// along **+z** (the axis) in v.
    fn quarter_cylinder(r: f64, h: f64) -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            2,
            1,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::new(r, 0.0, 0.0), Vector3::new(r, 0.0, h)],
                vec![Vector3::new(r, r, 0.0), Vector3::new(r, r, h)],
                vec![Vector3::new(0.0, r, 0.0), Vector3::new(0.0, r, h)],
            ],
            vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]],
        )
        .expect("valid quarter cylinder")
    }

    const Z: Vector3<f64> = Vector3::new(0.0, 0.0, 1.0);

    #[test]
    fn plane_facing_the_pull_has_90_degree_draft() {
        // Normal is ±ẑ, pull +ẑ → |draft| = 90° (face squarely faces the pull).
        let s = plane();
        let d = draft_angle(&s, 0.5, 0.5, Z);
        assert!((d.abs() - 90.0).abs() < 1e-4, "plane draft {d}, want ±90");
    }

    #[test]
    fn tilted_plane_has_ninety_minus_theta_draft() {
        // n makes angle θ with the pull → draft = asin(cosθ) = 90° − θ.
        for theta_deg in [10.0_f64, 30.0, 45.0, 60.0] {
            let theta = theta_deg.to_radians();
            let s = tilted_plane(theta);
            let d = draft_angle(&s, 0.5, 0.5, Z);
            let expected = 90.0 - theta_deg;
            assert!(
                (d - expected).abs() < 1e-3,
                "tilt {theta_deg}°: draft {d}, want {expected}"
            );
        }
    }

    #[test]
    fn cylinder_pulled_along_its_axis_has_zero_draft_walls() {
        // The walls are parallel to the axis (radial normals ⟂ pull) → ~0°
        // draft everywhere: a vertical, undrafted wall.
        let s = quarter_cylinder(2.0, 1.0);
        let report = draft_report(&s, Z, 10, 3.0);
        assert!(
            report.max_deg.abs() < 0.5 && report.min_deg.abs() < 0.5,
            "cylinder wall draft spread [{}, {}] should be ~0",
            report.min_deg,
            report.max_deg
        );
        // Zero-draft walls are NOT mouldable against a 3° requirement — they
        // need draft added (the whole point of the diagnostic).
        assert!(
            !report.mouldable,
            "vertical walls should fail a 3° draft check"
        );
        assert_eq!(report.undercut_count, report.samples);
    }

    #[test]
    fn report_flags_undercuts_and_passes_well_drafted_faces() {
        // A plane facing the pull is comfortably mouldable at a 3° minimum.
        let good = draft_report(&plane(), Z, 6, 3.0);
        assert!(good.mouldable && good.undercut_count == 0, "{good:?}");
        assert!((good.mean_deg.abs() - 90.0).abs() < 1e-3);

        // A face tilted 120° from the pull faces away (n · d = cos120° < 0) →
        // negative draft, an undercut: not mouldable.
        let under = draft_report(&tilted_plane(120.0_f64.to_radians()), Z, 6, 3.0);
        assert!(
            under.min_deg < 0.0,
            "should be an undercut, min {}",
            under.min_deg
        );
        assert!(
            !under.mouldable && under.undercut_count == under.samples,
            "{under:?}"
        );
    }

    #[test]
    fn degenerate_pull_is_zero() {
        assert_eq!(draft_angle(&plane(), 0.5, 0.5, Vector3::zeros()), 0.0);
    }
}

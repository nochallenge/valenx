//! [`Requirement`] — predicates over a [`Trace`] with a quantitative margin.
//!
//! A [`Requirement`] is a single safety/performance property the autonomous run
//! must satisfy. Evaluating it against a [`Trace`] yields a [`RequirementOutcome`]
//! carrying both a boolean **pass/fail** and a signed scalar **margin**.
//!
//! ## Margin sign convention (uniform across all requirements)
//!
//! - **Positive margin** ⇒ the requirement is *satisfied*, with that much slack
//!   (in the requirement's natural units — metres, seconds, …).
//! - **Negative margin** ⇒ the requirement is *violated*, by that much.
//! - The boolean pass is exactly `margin >= 0.0`.
//!
//! This single convention lets [`crate::sweep`] take a meaningful *minimum*
//! margin across a whole suite: the worst case is the smallest margin (the
//! nearest miss, or the deepest violation).
//!
//! The built-in requirements are checked against the harness **ground-truth**
//! [`crate::trace`]`::VehicleState` (the true pose at each tick), which is the
//! right reference for a safety property; building requirements over a *noisy
//! estimate* (e.g. "the estimator's position error stays < ε") is a documented
//! extension that would read the sensor frames instead.

use nalgebra::Vector3;
use valenx_sensors::{Scene, Surface};

use crate::error::VnvError;
use crate::trace::Trace;

/// An axis-aligned bounding box (min/max corner), used by
/// [`Requirement::StayInBounds`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// The minimum corner (m).
    pub min: Vector3<f64>,
    /// The maximum corner (m).
    pub max: Vector3<f64>,
}

impl Aabb {
    /// Build an AABB from two corners, ordering each axis so `min <= max`.
    ///
    /// # Errors
    /// [`VnvError::NonFinite`] if any component is non-finite.
    pub fn new(a: Vector3<f64>, b: Vector3<f64>) -> Result<Self, VnvError> {
        if !(a.iter().all(|x| x.is_finite()) && b.iter().all(|x| x.is_finite())) {
            return Err(VnvError::NonFinite("AABB corner".into()));
        }
        Ok(Self {
            min: Vector3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z)),
            max: Vector3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z)),
        })
    }

    /// Signed clearance of a point to the box boundary: positive (and equal to
    /// the distance to the nearest face) when the point is *inside*, negative
    /// (the overshoot past the nearest violated face) when *outside*.
    ///
    /// For a point inside, the slack is `min over axes of (p−min, max−p)`. For a
    /// point outside, the violation is the largest per-axis overshoot, returned
    /// negative. (The interior measure is exact; the exterior measure is a
    /// per-axis penetration depth, which is what a "stay in bounds" margin
    /// wants — how far out of bounds the worst axis went.)
    #[must_use]
    pub fn signed_clearance(&self, p: Vector3<f64>) -> f64 {
        let mut inside = true;
        let mut worst_overshoot = 0.0_f64;
        let mut min_slack = f64::INFINITY;
        for i in 0..3 {
            let lo = self.min[i];
            let hi = self.max[i];
            let pi = p[i];
            if pi < lo {
                inside = false;
                worst_overshoot = worst_overshoot.max(lo - pi);
            } else if pi > hi {
                inside = false;
                worst_overshoot = worst_overshoot.max(pi - hi);
            } else {
                min_slack = min_slack.min((pi - lo).min(hi - pi));
            }
        }
        if inside {
            min_slack
        } else {
            -worst_overshoot
        }
    }
}

/// A single safety/performance property to verify over a run.
#[derive(Debug, Clone, PartialEq)]
pub enum Requirement {
    /// The vehicle's ground-truth position must stay at least `d` metres clear
    /// of every surface in the run's [`Scene`] at every tick.
    ///
    /// Margin = `min over ticks (distance_to_scene(pos) − d)`. Pass iff the
    /// vehicle never came within `d` of a surface.
    MinClearance {
        /// Minimum required clearance distance (m, ≥ 0).
        d: f64,
    },

    /// A target placed at `target` must be detected by the LiDAR no later than
    /// `t_max` seconds — i.e. some beam returns a hit whose endpoint lands within
    /// `tol` metres of `target` at or before `t_max`.
    ///
    /// Margin = `t_max − t_detect` (seconds of slack), or `−(final_time − t_max)`
    /// if never detected within the run (a negative shortfall measured against
    /// the deadline). Pass iff detected by `t_max`.
    DetectByTime {
        /// The world-frame target point to detect (m).
        target: Vector3<f64>,
        /// The detection deadline (s, finite, ≥ 0).
        t_max: f64,
        /// How close a beam endpoint must land to `target` to count as a
        /// detection (m, > 0).
        tol: f64,
    },

    /// The vehicle's ground-truth position must stay inside the box `bounds` at
    /// every tick.
    ///
    /// Margin = `min over ticks (signed_clearance_to_box(pos))` — positive slack
    /// when always inside, negative penetration when it left. Pass iff never out
    /// of bounds.
    StayInBounds {
        /// The allowed region.
        bounds: Aabb,
    },

    /// The vehicle's ground-truth position must never collide with (enter) any
    /// surface in the scene, treating the vehicle as a sphere of `radius` metres.
    ///
    /// This is [`Requirement::MinClearance`] with `d = radius` framed as a
    /// collision check: margin = `min over ticks (distance_to_scene(pos) −
    /// radius)`. Pass iff the vehicle's hull never touched a surface.
    NoCollision {
        /// The vehicle's collision radius (m, ≥ 0).
        radius: f64,
    },
}

/// The result of evaluating one [`Requirement`] against one [`Trace`].
#[derive(Debug, Clone, PartialEq)]
pub struct RequirementOutcome {
    /// A short label describing the requirement (for reports).
    pub label: String,
    /// Whether the requirement passed (`margin >= 0`).
    pub pass: bool,
    /// The signed margin: positive slack if satisfied, negative violation if
    /// not. Units depend on the requirement (m for clearance/bounds/collision,
    /// s for detection deadlines).
    pub margin: f64,
    /// Whether this requirement was actually *exercised* by the trace — i.e. it
    /// had the data it needs to make a real pass/fail decision. (A
    /// `DetectByTime` over a trace with a LiDAR is exercised; coverage counts
    /// these.) Always `true` for the built-ins once they evaluate successfully;
    /// it exists so requirement-coverage can distinguish "ran and decided" from
    /// "skipped".
    pub exercised: bool,
}

impl Requirement {
    /// A short stable label for this requirement (used in reports and as a
    /// requirement-coverage key).
    #[must_use]
    pub fn label(&self) -> String {
        match self {
            Requirement::MinClearance { d } => format!("MinClearance(d={d})"),
            Requirement::DetectByTime { target, t_max, tol } => {
                format!(
                    "DetectByTime(target=[{:.3},{:.3},{:.3}], t_max={t_max}, tol={tol})",
                    target.x, target.y, target.z
                )
            }
            Requirement::StayInBounds { bounds } => format!(
                "StayInBounds(min=[{:.3},{:.3},{:.3}], max=[{:.3},{:.3},{:.3}])",
                bounds.min.x, bounds.min.y, bounds.min.z, bounds.max.x, bounds.max.y, bounds.max.z
            ),
            Requirement::NoCollision { radius } => format!("NoCollision(radius={radius})"),
        }
    }

    /// Validate the requirement's own parameters (finite, in range).
    ///
    /// # Errors
    /// [`VnvError::InvalidConfig`] / [`VnvError::NonFinite`].
    pub fn validate(&self) -> Result<(), VnvError> {
        match self {
            Requirement::MinClearance { d } | Requirement::NoCollision { radius: d } => {
                if !(d.is_finite() && *d >= 0.0) {
                    return Err(VnvError::InvalidConfig(format!(
                        "clearance/radius must be finite and ≥ 0, got {d}"
                    )));
                }
            }
            Requirement::DetectByTime { target, t_max, tol } => {
                if !target.iter().all(|x| x.is_finite()) {
                    return Err(VnvError::NonFinite("DetectByTime target".into()));
                }
                if !(t_max.is_finite() && *t_max >= 0.0) {
                    return Err(VnvError::InvalidConfig(format!(
                        "t_max must be finite and ≥ 0, got {t_max}"
                    )));
                }
                if !(tol.is_finite() && *tol > 0.0) {
                    return Err(VnvError::InvalidConfig(format!(
                        "tol must be finite and > 0, got {tol}"
                    )));
                }
            }
            Requirement::StayInBounds { bounds } => {
                if !(bounds.min.iter().all(|x| x.is_finite())
                    && bounds.max.iter().all(|x| x.is_finite()))
                {
                    return Err(VnvError::NonFinite("StayInBounds bounds".into()));
                }
            }
        }
        Ok(())
    }

    /// Evaluate this requirement against `trace`, returning a
    /// [`RequirementOutcome`].
    ///
    /// # Errors
    /// - The requirement's own [`Requirement::validate`] error, or
    /// - [`VnvError::RequirementMismatch`] if the requirement cannot be applied
    ///   to this trace: an empty trace for any requirement, or a
    ///   [`Requirement::DetectByTime`] over a trace that carries no LiDAR data
    ///   (so "detect" is unanswerable rather than trivially failing).
    pub fn evaluate(&self, trace: &Trace) -> Result<RequirementOutcome, VnvError> {
        self.validate()?;
        if trace.is_empty() {
            return Err(VnvError::RequirementMismatch(format!(
                "requirement '{}' cannot evaluate an empty trace '{}'",
                self.label(),
                trace.scenario
            )));
        }

        let margin = match self {
            Requirement::MinClearance { d } | Requirement::NoCollision { radius: d } => {
                min_clearance_margin(trace, *d)
            }
            Requirement::StayInBounds { bounds } => stay_in_bounds_margin(trace, bounds),
            Requirement::DetectByTime { target, t_max, tol } => {
                if !trace.has_lidar() {
                    return Err(VnvError::RequirementMismatch(format!(
                        "DetectByTime requires LiDAR but trace '{}' has none",
                        trace.scenario
                    )));
                }
                detect_by_time_margin(trace, *target, *t_max, *tol)
            }
        };

        Ok(RequirementOutcome {
            label: self.label(),
            pass: margin >= 0.0,
            margin,
            exercised: true,
        })
    }
}

/// Minimum over ticks of (distance from the vehicle to the nearest scene
/// surface − `d`).
fn min_clearance_margin(trace: &Trace, d: f64) -> f64 {
    trace
        .frames
        .iter()
        .map(|f| distance_to_scene(&trace.scene, f.state.position) - d)
        .fold(f64::INFINITY, f64::min)
}

/// Minimum over ticks of the signed clearance of the vehicle to the box.
fn stay_in_bounds_margin(trace: &Trace, bounds: &Aabb) -> f64 {
    trace
        .frames
        .iter()
        .map(|f| bounds.signed_clearance(f.state.position))
        .fold(f64::INFINITY, f64::min)
}

/// `t_max − t_detect` if the target was detected by the deadline; otherwise a
/// negative shortfall `−(final_time − t_max)` (or, if never detected at all and
/// the run ended before the deadline, the same negative-distance-to-deadline so
/// the margin is consistently negative on failure).
fn detect_by_time_margin(trace: &Trace, target: Vector3<f64>, t_max: f64, tol: f64) -> f64 {
    for f in &trace.frames {
        if let Some(scan) = &f.lidar {
            // A beam endpoint = sensor (vehicle) position + range · world dir.
            // The harness mounts the LiDAR at the body origin aligned to the
            // body axes, so the world ray direction is orientation · body_dir.
            for beam in &scan.beams {
                if let Some(range) = beam.range {
                    let world_dir = f.state.orientation * beam.direction;
                    let endpoint = f.state.position + world_dir * range;
                    if (endpoint - target).norm() <= tol {
                        return t_max - f.time;
                    }
                }
            }
        }
    }
    // Never detected within the run: a negative margin measured against the
    // deadline. If the run extended past t_max, this is how late we are; if it
    // ended early, it is how much sooner than the deadline the run stopped
    // (still a failure, still negative).
    -(trace.final_time() - t_max).abs().max(f64::MIN_POSITIVE)
}

/// Distance from a point to the nearest surface in a scene (`+∞` for an empty
/// scene). Planes use the exact point-plane distance; spheres the |center−p| −
/// r distance (negative inside, clamped at the surface for "touching"); and
/// triangles the exact point-triangle distance.
fn distance_to_scene(scene: &Scene, p: Vector3<f64>) -> f64 {
    scene
        .surfaces
        .iter()
        .map(|s| distance_to_surface(s, p))
        .fold(f64::INFINITY, f64::min)
}

/// Distance from a point to a single surface.
fn distance_to_surface(surface: &Surface, p: Vector3<f64>) -> f64 {
    match surface {
        Surface::Plane(plane) => plane.normal.dot(&(p - plane.point)).abs(),
        Surface::Sphere(sphere) => ((p - sphere.center).norm() - sphere.radius).max(0.0),
        Surface::Triangle(tri) => point_triangle_distance(p, tri.v0, tri.v1, tri.v2),
    }
}

/// Exact unsigned distance from point `p` to triangle `(a, b, c)` (Ericson,
/// *Real-Time Collision Detection*, closest-point-on-triangle).
fn point_triangle_distance(
    p: Vector3<f64>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (p - a).norm(); // vertex region A
    }
    let bp = p - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (p - b).norm(); // vertex region B
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (p - (a + ab * v)).norm(); // edge AB
    }
    let cp = p - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (p - c).norm(); // vertex region C
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (p - (a + ac * w)).norm(); // edge AC
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).norm(); // edge BC
    }
    // Face region: project onto the triangle plane via barycentric coords.
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    (p - (a + ab * v + ac * w)).norm()
}

/// A named, validated collection of [`Requirement`]s evaluated together as the
/// acceptance criteria for a run.
#[derive(Debug, Clone, PartialEq)]
pub struct RequirementSet {
    /// The requirements (may be empty — an empty set passes vacuously, which the
    /// report records honestly).
    pub requirements: Vec<Requirement>,
}

impl RequirementSet {
    /// Build a requirement set.
    #[must_use]
    pub fn new(requirements: Vec<Requirement>) -> Self {
        Self { requirements }
    }

    /// Number of requirements.
    #[must_use]
    pub fn len(&self) -> usize {
        self.requirements.len()
    }

    /// Whether the set is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.requirements.is_empty()
    }

    /// Validate every requirement's parameters.
    ///
    /// # Errors
    /// The first requirement validation error.
    pub fn validate(&self) -> Result<(), VnvError> {
        for r in &self.requirements {
            r.validate()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn aabb_orders_corners_and_measures_interior_slack() {
        // Corners given out of order are normalised.
        let b = Aabb::new(v(10.0, 1.0, 1.0), v(0.0, -1.0, -1.0)).unwrap();
        assert_eq!(b.min, v(0.0, -1.0, -1.0));
        assert_eq!(b.max, v(10.0, 1.0, 1.0));
        // Centre point: slack to nearest face = min half-extent = 1.0 (y/z).
        assert!((b.signed_clearance(v(5.0, 0.0, 0.0)) - 1.0).abs() < 1e-12);
        // On a face ⇒ zero slack.
        assert!(b.signed_clearance(v(0.0, 0.0, 0.0)).abs() < 1e-12);
    }

    #[test]
    fn aabb_exterior_clearance_is_negative_worst_overshoot() {
        let b = Aabb::new(v(0.0, 0.0, 0.0), v(10.0, 10.0, 10.0)).unwrap();
        // Out on +x by 2 and +y by 5 ⇒ worst overshoot 5 ⇒ −5.
        assert!((b.signed_clearance(v(12.0, 15.0, 5.0)) - (-5.0)).abs() < 1e-12);
    }

    #[test]
    fn aabb_rejects_nonfinite() {
        assert!(Aabb::new(v(f64::NAN, 0.0, 0.0), v(1.0, 1.0, 1.0)).is_err());
    }

    #[test]
    fn point_plane_distance_is_perpendicular() {
        // Plane x=5, point at x=2 ⇒ distance 3.
        let plane = valenx_sensors::Plane::new(v(5.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let d = distance_to_surface(&Surface::Plane(plane), v(2.0, 100.0, -50.0));
        assert!((d - 3.0).abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn point_sphere_distance_clamps_inside_to_zero() {
        let s = valenx_sensors::Sphere::new(v(0.0, 0.0, 0.0), 2.0).unwrap();
        // Outside at x=5 ⇒ 5−2 = 3.
        assert!((distance_to_surface(&Surface::Sphere(s), v(5.0, 0.0, 0.0)) - 3.0).abs() < 1e-12);
        // Inside ⇒ clamped to 0 (a hull "inside the obstacle" reads as touching).
        assert!(distance_to_surface(&Surface::Sphere(s), v(0.5, 0.0, 0.0)).abs() < 1e-12);
    }

    #[test]
    fn point_triangle_distance_matches_known_geometry() {
        // Triangle in the plane z=0; a point straight above the centroid by 3 ⇒
        // distance 3 (face region).
        let a = v(-1.0, -1.0, 0.0);
        let b = v(1.0, -1.0, 0.0);
        let c = v(0.0, 1.0, 0.0);
        let centroid = (a + b + c) / 3.0;
        let p = centroid + v(0.0, 0.0, 3.0);
        let d = point_triangle_distance(p, a, b, c);
        assert!((d - 3.0).abs() < 1e-12, "face distance = {d}");
        // A point at vertex A ⇒ distance 0 (vertex region).
        assert!(point_triangle_distance(a, a, b, c).abs() < 1e-12);
    }
}

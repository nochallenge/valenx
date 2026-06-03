//! Helical-axis fitting for nucleic-acid duplexes.
//!
//! Two methods are provided:
//!
//! 1. [`fit_helical_axis`] — a global **straight** axis fitted by
//!    total-least-squares + chord-difference + Kåsa circle fitting.
//!    The v1 path, correct and useful for roughly linear duplexes.
//! 2. [`fit_curved_axis`] — a **Curves+-class** curved axis: a real
//!    per-bp local axis point is derived from each base pair's frame
//!    plus the orientation of the next pair (the canonical Curves+
//!    "local helical reference point" — the intersection of the
//!    rotation axis taking pair `k` onto pair `k+1` with the pair-`k`
//!    base-pair plane). A natural cubic spline is fitted through
//!    those per-bp axis points; per-bp curvature `κ` is computed
//!    analytically from the spline.
//!
//! Both fits report rise / twist / radius / bp-per-turn; the curved
//! axis additionally exposes a curvature profile and a cumulative
//! axis-arc-length sampling.

use crate::error::{BiostructError, Result};
use crate::nucleic::params::{step_parameters, BaseFrame};
use nalgebra::{Matrix3, Point3, Vector3};

/// A fitted global helical axis and its descriptors.
#[derive(Clone, Debug, PartialEq)]
pub struct HelicalAxis {
    /// A point on the axis (the centroid of the base-pair origins).
    pub point: Point3<f64>,
    /// Unit direction of the axis.
    pub direction: Vector3<f64>,
    /// Mean rise per base-pair step, ångström.
    pub rise_per_bp: f64,
    /// Mean twist per base-pair step, degrees.
    pub twist_per_bp: f64,
    /// Mean radius of the base-pair centres from the axis, ångström.
    pub radius: f64,
    /// Base pairs per helical turn (`360 / twist_per_bp`).
    pub bp_per_turn: f64,
    /// Total contour length of the axis between the first and last
    /// base-pair projections, ångström.
    pub length: f64,
}

/// Fit a helical axis to a sequence of base-pair centre points
/// (typically the [`pair_mid_frame`](super::params::pair_mid_frame)
/// origins, in 5′→3′ order).
///
/// At least 4 points are required (a helix axis is under-determined
/// from fewer).
///
/// # Axis-direction method
///
/// The naive total-least-squares line through the base-pair centres
/// does **not** recover the helix axis: a helix of a non-integer
/// number of turns has nonzero `Cov(x, z)` / `Cov(y, z)`, which tilts
/// the first principal axis away from the true axis.
///
/// Instead this uses the chord-difference construction. With chord
/// vectors `dₖ = pₖ₊₁ − pₖ`, every chord of an ideal helix has the
/// *same* axial component, so the second differences
/// `eₖ = dₖ₊₁ − dₖ` are exactly **perpendicular** to the axis. The
/// axis direction is therefore the normal of the best-fit plane of
/// the `eₖ` — the smallest-eigenvalue eigenvector of `Σ eₖ eₖᵀ`. This
/// is exact for an ideal helix and robust for a near-ideal one,
/// independent of how many turns the duplex spans.
pub fn fit_helical_axis(centers: &[Point3<f64>]) -> Result<HelicalAxis> {
    if centers.len() < 4 {
        return Err(BiostructError::invalid(
            "centers",
            "helical-axis fit needs at least 4 base-pair centres",
        ));
    }

    // --- centroid (a point on the axis) ----------------------------
    let mut acc = Vector3::zeros();
    for p in centers {
        acc += p.coords;
    }
    let centroid = Point3::from(acc / centers.len() as f64);

    // --- axis direction from chord second differences --------------
    // chord vectors dₖ = pₖ₊₁ − pₖ
    let chords: Vec<Vector3<f64>> = centers
        .windows(2)
        .map(|w| w[1].coords - w[0].coords)
        .collect();
    // second differences eₖ = dₖ₊₁ − dₖ lie ⟂ to the helix axis
    let mut perp_scatter = Matrix3::zeros();
    for w in chords.windows(2) {
        let e = w[1] - w[0];
        perp_scatter += e * e.transpose();
    }
    let perp_eig = nalgebra::SymmetricEigen::new(perp_scatter);
    // smallest-eigenvalue eigenvector = normal of the eₖ plane = axis
    let mut min_i = 0;
    for i in 1..3 {
        if perp_eig.eigenvalues[i] < perp_eig.eigenvalues[min_i] {
            min_i = i;
        }
    }
    let mut direction = perp_eig
        .eigenvectors
        .column(min_i)
        .into_owned()
        .normalize();

    // Orient the axis 5'->3' (toward the last base pair).
    let span = centers[centers.len() - 1] - centers[0];
    if direction.dot(&span) < 0.0 {
        direction = -direction;
    }

    // --- locate the axis line ---------------------------------------
    // The centroid lies *on* the helix axis only for an integer
    // number of turns; in general it is offset by the (non-cancelling)
    // mean radial vector. The true axis pierces the centre of the
    // circle that the base-pair centres trace when projected onto the
    // plane ⟂ to `direction`. Fit that circle (Kåsa algebraic fit) and
    // use its centre.
    let (u, v) = perpendicular_basis(&direction);
    let proj: Vec<(f64, f64)> = centers
        .iter()
        .map(|p| {
            let r = p.coords - centroid.coords;
            (r.dot(&u), r.dot(&v))
        })
        .collect();
    let (cu, cv, _circle_r) = fit_circle(&proj);
    // A point on the axis: the centroid shifted by the circle centre
    // within the perpendicular plane.
    let axis_point = Point3::from(centroid.coords + u * cu + v * cv);

    // --- project centres onto the axis -----------------------------
    // Axial coordinate t and the radial (perpendicular) vector of
    // each base-pair centre, measured from the axis line.
    let mut axial: Vec<f64> = Vec::with_capacity(centers.len());
    let mut radial: Vec<Vector3<f64>> = Vec::with_capacity(centers.len());
    for p in centers {
        let r = p.coords - axis_point.coords;
        let t = r.dot(&direction);
        axial.push(t);
        radial.push(r - direction * t);
    }

    // --- rise: mean axial advance per step -------------------------
    let mut rise_sum = 0.0;
    for k in 1..axial.len() {
        rise_sum += (axial[k] - axial[k - 1]).abs();
    }
    let rise_per_bp = rise_sum / (axial.len() - 1) as f64;

    // --- twist: mean rotation of the radial vector per step --------
    let mut twist_sum = 0.0;
    let mut twist_count = 0;
    for k in 1..radial.len() {
        let a = radial[k - 1];
        let b = radial[k];
        if a.norm() > 1e-6 && b.norm() > 1e-6 {
            let an = a.normalize();
            let bn = b.normalize();
            let cos = an.dot(&bn).clamp(-1.0, 1.0);
            twist_sum += cos.acos().to_degrees();
            twist_count += 1;
        }
    }
    let twist_per_bp = if twist_count > 0 {
        twist_sum / twist_count as f64
    } else {
        0.0
    };

    // --- radius: mean perpendicular distance -----------------------
    let radius =
        radial.iter().map(|r| r.norm()).sum::<f64>() / radial.len() as f64;

    let bp_per_turn = if twist_per_bp.abs() > 1e-6 {
        360.0 / twist_per_bp
    } else {
        f64::INFINITY
    };

    let t_min = axial.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = axial.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let length = t_max - t_min;

    Ok(HelicalAxis {
        point: axis_point,
        direction,
        rise_per_bp,
        twist_per_bp,
        radius,
        bp_per_turn,
        length,
    })
}

/// Two orthonormal vectors spanning the plane perpendicular to a unit
/// direction `d`.
fn perpendicular_basis(d: &Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    // Pick the global axis least aligned with `d` to avoid degeneracy.
    let seed = if d.x.abs() <= d.y.abs() && d.x.abs() <= d.z.abs() {
        Vector3::new(1.0, 0.0, 0.0)
    } else if d.y.abs() <= d.z.abs() {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        Vector3::new(0.0, 0.0, 1.0)
    };
    let u = (seed - d * seed.dot(d)).normalize();
    let v = d.cross(&u);
    (u, v)
}

/// Least-squares (Kåsa algebraic) circle fit to 2-D points, returning
/// the circle centre `(cx, cy)` and radius.
///
/// Solves the linear system that minimises `Σ (xᵢ² + yᵢ² − 2cx·xᵢ −
/// 2cy·yᵢ − k)²`. Degenerate inputs (collinear points) fall back to
/// the point centroid with a radius of the mean distance to it.
fn fit_circle(points: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = points.len() as f64;
    if n < 3.0 {
        return (0.0, 0.0, 0.0);
    }
    let (mut sx, mut sy) = (0.0, 0.0);
    for &(x, y) in points {
        sx += x;
        sy += y;
    }
    let (mx, my) = (sx / n, sy / n);

    // Work in centred coordinates for numerical conditioning.
    let (mut suu, mut suv, mut svv) = (0.0, 0.0, 0.0);
    let (mut suuu, mut svvv, mut suvv, mut svuu) = (0.0, 0.0, 0.0, 0.0);
    for &(x, y) in points {
        let u = x - mx;
        let v = y - my;
        suu += u * u;
        suv += u * v;
        svv += v * v;
        suuu += u * u * u;
        svvv += v * v * v;
        suvv += u * v * v;
        svuu += v * u * u;
    }
    // [ suu suv ] [uc]   [ (suuu + suvv) / 2 ]
    // [ suv svv ] [vc] = [ (svvv + svuu) / 2 ]
    let det = suu * svv - suv * suv;
    if det.abs() < 1e-12 {
        // collinear / coincident — no circle; fall back to centroid.
        let r = points
            .iter()
            .map(|&(x, y)| ((x - mx).powi(2) + (y - my).powi(2)).sqrt())
            .sum::<f64>()
            / n;
        return (mx, my, r);
    }
    let b1 = (suuu + suvv) / 2.0;
    let b2 = (svvv + svuu) / 2.0;
    let uc = (b1 * svv - b2 * suv) / det;
    let vc = (suu * b2 - suv * b1) / det;
    let cx = uc + mx;
    let cy = vc + my;
    let radius = (uc * uc + vc * vc + (suu + svv) / n).sqrt();
    (cx, cy, radius)
}

/// Perpendicular distance of a point from a fitted helical axis.
pub fn distance_to_axis(axis: &HelicalAxis, p: &Point3<f64>) -> f64 {
    let r = p - axis.point;
    let t = r.dot(&axis.direction);
    (r - axis.direction * t).norm()
}

/// Build an idealised helix of `n` base-pair centres for a given
/// rise, twist (degrees) and radius — a convenience for testing and
/// for previewing canonical B/A/Z-DNA geometry.
pub fn ideal_helix_centers(
    n: usize,
    rise: f64,
    twist_deg: f64,
    radius: f64,
) -> Vec<Point3<f64>> {
    let twist = twist_deg.to_radians();
    (0..n)
        .map(|i| {
            let theta = i as f64 * twist;
            Point3::new(
                radius * theta.cos(),
                radius * theta.sin(),
                i as f64 * rise,
            )
        })
        .collect()
}

// ===================================================================
//  Curves+-class curved helical-axis fitting.
// ===================================================================

/// A fitted Curves+-class **curved** helical axis: a natural cubic
/// spline through per-bp axis reference points, plus per-bp curvature.
#[derive(Clone, Debug, PartialEq)]
pub struct CurvedHelicalAxis {
    /// Per-bp local helical axis points, in 5′→3′ order. There is one
    /// axis point per base-pair frame (the first and last are
    /// extrapolated from the available steps).
    pub axis_points: Vec<Point3<f64>>,
    /// Cumulative arc length along the polyline through `axis_points`,
    /// in ångström. `arc_length[k]` is the arc length up to
    /// `axis_points[k]`; the total contour length is the last entry.
    pub arc_length: Vec<f64>,
    /// Cubic-spline polynomial segments — one per interval `[k, k+1]`
    /// of `axis_points`. Each segment carries the cubic polynomial
    /// coefficients in component form (x, y, z separately, in the
    /// segment's local parameter `s ∈ [0, h_k]`, `h_k = arc_length[k+1] − arc_length[k]`).
    pub segments: Vec<SplineSegment>,
    /// Per-bp curvature `κ` (in Å⁻¹) measured at each axis point.
    /// `curvature[k]` is the curvature of the spline at
    /// `arc_length[k]`. A near-zero value means an almost straight
    /// duplex; large values indicate a bend.
    pub curvature: Vec<f64>,
    /// Mean rise per step (ångström) — the mean of the
    /// per-step axis-point spacings.
    pub rise_per_bp: f64,
    /// Mean twist per step (degrees) — taken from the underlying
    /// 3DNA step parameters' twist.
    pub twist_per_bp: f64,
    /// Mean radius — average distance from the base-pair centres to
    /// the nearest axis-point segment.
    pub radius: f64,
    /// Base pairs per helical turn (`360 / twist_per_bp`).
    pub bp_per_turn: f64,
}

/// One cubic spline segment over an arc-length interval `[0, h]`.
/// The 3-D curve is `(x(s), y(s), z(s)) = (a + b·s + c·s² + d·s³, …)`.
#[derive(Clone, Debug, PartialEq)]
pub struct SplineSegment {
    /// Length of this segment in the original arc-length parameter.
    pub h: f64,
    /// `a` coefficients (x, y, z) — value at `s = 0`.
    pub a: [f64; 3],
    /// `b` coefficients (first derivative at `s = 0`).
    pub b: [f64; 3],
    /// `c` coefficients (second-derivative term / 2).
    pub c: [f64; 3],
    /// `d` coefficients (third-derivative term / 6).
    pub d: [f64; 3],
}

impl SplineSegment {
    /// Evaluate the polynomial at local parameter `s`.
    pub fn evaluate(&self, s: f64) -> Point3<f64> {
        Point3::new(
            self.a[0] + self.b[0] * s + self.c[0] * s * s + self.d[0] * s * s * s,
            self.a[1] + self.b[1] * s + self.c[1] * s * s + self.d[1] * s * s * s,
            self.a[2] + self.b[2] * s + self.c[2] * s * s + self.d[2] * s * s * s,
        )
    }

    /// First derivative `r'(s)`.
    pub fn first_derivative(&self, s: f64) -> Vector3<f64> {
        Vector3::new(
            self.b[0] + 2.0 * self.c[0] * s + 3.0 * self.d[0] * s * s,
            self.b[1] + 2.0 * self.c[1] * s + 3.0 * self.d[1] * s * s,
            self.b[2] + 2.0 * self.c[2] * s + 3.0 * self.d[2] * s * s,
        )
    }

    /// Second derivative `r''(s)`.
    pub fn second_derivative(&self, s: f64) -> Vector3<f64> {
        Vector3::new(
            2.0 * self.c[0] + 6.0 * self.d[0] * s,
            2.0 * self.c[1] + 6.0 * self.d[1] * s,
            2.0 * self.c[2] + 6.0 * self.d[2] * s,
        )
    }
}

/// Derive a per-bp local helical-axis point from a pair of consecutive
/// base-pair frames.
///
/// Computes the **screw axis** of the rigid transform taking the
/// lower base-pair frame onto the upper base-pair frame. The screw
/// axis is the set of points fixed (up to a translation along the axis)
/// by the transform — the canonical Curves+ "local helical reference
/// point" lives on this axis, at the foot of the perpendicular from
/// `lower.origin` to it.
///
/// Construction: let `R = R_u · R_lᵀ` be the rotation taking `lower`
/// onto `upper`, `t = upper.origin − R · lower.origin` the
/// corresponding translation. The rotation-axis direction `u` is the
/// eigenvector of `R` for eigenvalue +1 (extracted by the standard
/// `(R32 − R23, R13 − R31, R21 − R12)` formula). A point on the screw
/// axis solves `(I − R)·p_⊥ = t_⊥`, where `_⊥` denotes the component
/// perpendicular to `u`; the equation is solved in the 2-D plane
/// perpendicular to `u`. The returned axis point is then the foot of
/// the perpendicular from `lower.origin` onto the axis line through
/// `p_⊥` with direction `u`.
fn local_axis_point_from_step(lower: &BaseFrame, upper: &BaseFrame) -> Point3<f64> {
    let r_l = lower.rotation();
    let r_u = upper.rotation();
    let r = r_u * r_l.transpose();
    let t = upper.origin - r * lower.origin.coords;
    // Rotation axis direction.
    let u_unnorm = Vector3::new(
        r[(2, 1)] - r[(1, 2)],
        r[(0, 2)] - r[(2, 0)],
        r[(1, 0)] - r[(0, 1)],
    );
    if u_unnorm.norm() < 1e-9 {
        // Near-identity rotation — axis is ill-defined; fall back to
        // the local frame's z (the base-plane normal) and put the
        // axis point at the midpoint.
        return Point3::from((lower.origin.coords + upper.origin.coords) * 0.5);
    }
    let u = u_unnorm.normalize();
    let t_axial = u * u.dot(&t.coords);
    let t_perp = t.coords - t_axial;
    // Build orthonormal basis (u, e1, e2) of R³ with e1, e2 spanning
    // the plane perpendicular to u.
    let seed = if u.x.abs() <= u.y.abs() && u.x.abs() <= u.z.abs() {
        Vector3::new(1.0, 0.0, 0.0)
    } else if u.y.abs() <= u.z.abs() {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        Vector3::new(0.0, 0.0, 1.0)
    };
    let e1 = (seed - u * seed.dot(&u)).normalize();
    let e2 = u.cross(&e1);
    // In the (e1, e2) plane, `(I − R)` acts as a 2×2 invertible matrix
    // for any non-identity rotation. Build it.
    let m11 = e1.dot(&(e1 - r * e1));
    let m12 = e1.dot(&(e2 - r * e2));
    let m21 = e2.dot(&(e1 - r * e1));
    let m22 = e2.dot(&(e2 - r * e2));
    let rhs1 = e1.dot(&t_perp);
    let rhs2 = e2.dot(&t_perp);
    let det = m11 * m22 - m12 * m21;
    if det.abs() < 1e-9 {
        return Point3::from((lower.origin.coords + upper.origin.coords) * 0.5);
    }
    let p1 = (rhs1 * m22 - rhs2 * m12) / det;
    let p2 = (m11 * rhs2 - m21 * rhs1) / det;
    let p_perp = e1 * p1 + e2 * p2;
    // Now the screw axis is the line { p_perp + s·u : s ∈ R }.
    // Pick the foot of perpendicular from `lower.origin` onto it.
    let v = lower.origin.coords - p_perp;
    let s = v.dot(&u);
    Point3::from(p_perp + u * s)
}

/// Same screw-axis derivation as [`local_axis_point_from_step`], but
/// returns the foot of perpendicular from `upper.origin` to the screw
/// axis — the appropriate axis estimate for base pair `k+1`.
fn local_axis_point_for_upper(lower: &BaseFrame, upper: &BaseFrame) -> Point3<f64> {
    let r_l = lower.rotation();
    let r_u = upper.rotation();
    let r = r_u * r_l.transpose();
    let t = upper.origin - r * lower.origin.coords;
    let u_unnorm = Vector3::new(
        r[(2, 1)] - r[(1, 2)],
        r[(0, 2)] - r[(2, 0)],
        r[(1, 0)] - r[(0, 1)],
    );
    if u_unnorm.norm() < 1e-9 {
        return Point3::from((lower.origin.coords + upper.origin.coords) * 0.5);
    }
    let u = u_unnorm.normalize();
    let t_axial = u * u.dot(&t.coords);
    let t_perp = t.coords - t_axial;
    let seed = if u.x.abs() <= u.y.abs() && u.x.abs() <= u.z.abs() {
        Vector3::new(1.0, 0.0, 0.0)
    } else if u.y.abs() <= u.z.abs() {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        Vector3::new(0.0, 0.0, 1.0)
    };
    let e1 = (seed - u * seed.dot(&u)).normalize();
    let e2 = u.cross(&e1);
    let m11 = e1.dot(&(e1 - r * e1));
    let m12 = e1.dot(&(e2 - r * e2));
    let m21 = e2.dot(&(e1 - r * e1));
    let m22 = e2.dot(&(e2 - r * e2));
    let rhs1 = e1.dot(&t_perp);
    let rhs2 = e2.dot(&t_perp);
    let det = m11 * m22 - m12 * m21;
    if det.abs() < 1e-9 {
        return Point3::from((lower.origin.coords + upper.origin.coords) * 0.5);
    }
    let p1 = (rhs1 * m22 - rhs2 * m12) / det;
    let p2 = (m11 * rhs2 - m21 * rhs1) / det;
    let p_perp = e1 * p1 + e2 * p2;
    let v = upper.origin.coords - p_perp;
    let s = v.dot(&u);
    Point3::from(p_perp + u * s)
}

/// Solve a tridiagonal linear system `A·x = b` in O(n) via the Thomas
/// algorithm. `lower` has `n-1` sub-diagonal entries, `diag` has `n`
/// diagonal entries, `upper` has `n-1` super-diagonal entries.
/// Returns the solution `x`.
fn solve_tridiag(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![rhs[0] / diag[0]];
    }
    let mut c_star = vec![0.0_f64; n - 1];
    let mut d_star = vec![0.0_f64; n];
    c_star[0] = upper[0] / diag[0];
    d_star[0] = rhs[0] / diag[0];
    for i in 1..n {
        let denom = diag[i] - lower[i - 1] * c_star.get(i - 1).copied().unwrap_or(0.0);
        if i < n - 1 {
            c_star[i] = upper[i] / denom;
        }
        d_star[i] = (rhs[i] - lower[i - 1] * d_star[i - 1]) / denom;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d_star[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
    x
}

/// Build a natural cubic spline through the given knot points,
/// parameterised by arc length along the polyline. Returns one
/// [`SplineSegment`] per interval and the cumulative arc-length array.
fn natural_cubic_spline(
    points: &[Point3<f64>],
) -> (Vec<SplineSegment>, Vec<f64>) {
    let n = points.len();
    assert!(n >= 3, "cubic spline needs at least 3 knots");
    // Chord-length parameterisation.
    let mut t = vec![0.0_f64; n];
    for i in 1..n {
        t[i] = t[i - 1] + (points[i] - points[i - 1]).norm();
    }
    let mut h = vec![0.0_f64; n - 1];
    for i in 0..n - 1 {
        h[i] = (t[i + 1] - t[i]).max(1e-9);
    }

    // For each component (x, y, z), solve the natural-cubic spline
    // tridiagonal system for the second-derivative values at the
    // interior knots; clamp boundary second derivative to 0.
    let mut m: [Vec<f64>; 3] = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    if n >= 3 {
        // Number of interior equations.
        let p = n - 2;
        // Tridiagonal coefficients (interior knots i = 1..n-1).
        // `lower` has p-1 sub-diagonal entries, `upper` has p-1
        // super-diagonal entries (empty when p == 1).
        let lower: Vec<f64> = if p >= 2 {
            (1..p).map(|i| h[i]).collect()
        } else {
            Vec::new()
        };
        let mut diag = vec![0.0_f64; p];
        for i in 0..p {
            diag[i] = 2.0 * (h[i] + h[i + 1]);
        }
        let upper: Vec<f64> = if p >= 2 {
            (0..p - 1).map(|i| h[i + 1]).collect()
        } else {
            Vec::new()
        };
        for (comp, m_comp) in m.iter_mut().enumerate() {
            let mut rhs = vec![0.0_f64; p];
            for (i, slot) in rhs.iter_mut().enumerate().take(p) {
                let f_im1 = points[i].coords[comp];
                let f_i = points[i + 1].coords[comp];
                let f_ip1 = points[i + 2].coords[comp];
                *slot = 6.0 * ((f_ip1 - f_i) / h[i + 1] - (f_i - f_im1) / h[i]);
            }
            let interior = solve_tridiag(&lower, &diag, &upper, &rhs);
            m_comp[1..=p].copy_from_slice(&interior[..p]);
        }
    }

    let mut segments = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let mut a = [0.0_f64; 3];
        let mut b = [0.0_f64; 3];
        let mut c = [0.0_f64; 3];
        let mut d = [0.0_f64; 3];
        for comp in 0..3 {
            let f_i = points[i].coords[comp];
            let f_ip1 = points[i + 1].coords[comp];
            // Standard natural-cubic-spline form:
            // f(s) = f_i + b·s + (M_i / 2)·s² + ((M_{i+1} − M_i) / (6h))·s³,
            // with b chosen so f(h) = f_{i+1}.
            a[comp] = f_i;
            c[comp] = m[comp][i] / 2.0;
            d[comp] = (m[comp][i + 1] - m[comp][i]) / (6.0 * h[i]);
            b[comp] = (f_ip1 - f_i) / h[i]
                - h[i] * (2.0 * m[comp][i] + m[comp][i + 1]) / 6.0;
        }
        segments.push(SplineSegment {
            h: h[i],
            a,
            b,
            c,
            d,
        });
    }
    (segments, t)
}

/// Curvature `κ` of a 3-D curve at a point given its first and second
/// derivatives — the standard `|r' × r''| / |r'|³` formula.
fn curvature_at(r1: Vector3<f64>, r2: Vector3<f64>) -> f64 {
    let r1_norm = r1.norm();
    if r1_norm < 1e-9 {
        return 0.0;
    }
    r1.cross(&r2).norm() / r1_norm.powi(3)
}

/// Fit a Curves+-class **curved** helical axis to a sequence of
/// base-pair frames (typically the `pair_mid_frame` outputs in
/// 5′→3′ order). At least 3 base-pair frames are required.
pub fn fit_curved_axis(frames: &[BaseFrame]) -> Result<CurvedHelicalAxis> {
    if frames.len() < 3 {
        return Err(BiostructError::invalid(
            "frames",
            "curved helical-axis fit needs at least 3 base-pair frames",
        ));
    }

    // --- per-bp local axis points -----------------------------------
    // For each step (k, k+1) we derive an axis point at the foot of
    // perpendicular from `lower.origin`. To get one axis point per
    // base pair we project both bp origins onto the step's screw
    // axis: this yields a "lower" point at the foot for frame k and
    // an "upper" point at the foot for frame k+1. Interior base pairs
    // appear in two adjacent steps and the two estimates are averaged.
    let n = frames.len();
    let mut axis_pts: Vec<Point3<f64>> = Vec::with_capacity(n);
    let mut accum: Vec<Vec<Point3<f64>>> = vec![Vec::new(); n];
    for k in 0..n - 1 {
        let lower_foot = local_axis_point_from_step(&frames[k], &frames[k + 1]);
        let upper_foot =
            local_axis_point_for_upper(&frames[k], &frames[k + 1]);
        accum[k].push(lower_foot);
        accum[k + 1].push(upper_foot);
    }
    for pts in &accum {
        let mut sum = Vector3::zeros();
        for p in pts {
            sum += p.coords;
        }
        axis_pts.push(Point3::from(sum / pts.len() as f64));
    }

    // --- spline fit + arc length ------------------------------------
    let (segments, arc_length) = natural_cubic_spline(&axis_pts);

    // --- per-bp curvature -------------------------------------------
    let mut curvature = vec![0.0_f64; n];
    for (k, c) in curvature.iter_mut().enumerate() {
        // Pick the segment containing arc_length[k]; for interior knots
        // we average the curvature at the right end of the previous
        // segment and the left end of the next, to avoid a kink
        // artefact at the knot.
        if k == 0 {
            let seg = &segments[0];
            *c = curvature_at(seg.first_derivative(0.0), seg.second_derivative(0.0));
        } else if k == n - 1 {
            let seg = &segments[n - 2];
            let s = seg.h;
            *c = curvature_at(seg.first_derivative(s), seg.second_derivative(s));
        } else {
            let left = &segments[k - 1];
            let right = &segments[k];
            let s = left.h;
            let kl = curvature_at(left.first_derivative(s), left.second_derivative(s));
            let kr = curvature_at(right.first_derivative(0.0), right.second_derivative(0.0));
            *c = 0.5 * (kl + kr);
        }
    }

    // --- mean rise / radius -----------------------------------------
    let mut rise_sum = 0.0;
    for k in 1..axis_pts.len() {
        rise_sum += (axis_pts[k] - axis_pts[k - 1]).norm();
    }
    let rise_per_bp = rise_sum / (axis_pts.len() - 1) as f64;

    // Mean twist per step: derived from the 3DNA step parameters'
    // twist component (which is signed and consistent with the
    // sequence direction).
    let mut twist_sum = 0.0;
    for k in 1..n {
        let sp = step_parameters(&frames[k - 1], &frames[k]);
        twist_sum += sp.twist;
    }
    let twist_per_bp = twist_sum / (n - 1) as f64;

    // Radius — mean perpendicular distance from each frame origin to
    // its derived axis point.
    let mut radius_sum = 0.0;
    for k in 0..n {
        radius_sum += (frames[k].origin - axis_pts[k]).norm();
    }
    let radius = radius_sum / n as f64;

    let bp_per_turn = if twist_per_bp.abs() > 1e-6 {
        360.0 / twist_per_bp.abs()
    } else {
        f64::INFINITY
    };

    Ok(CurvedHelicalAxis {
        axis_points: axis_pts,
        arc_length,
        segments,
        curvature,
        rise_per_bp,
        twist_per_bp,
        radius,
        bp_per_turn,
    })
}

impl CurvedHelicalAxis {
    /// Total contour length of the spline.
    pub fn contour_length(&self) -> f64 {
        *self.arc_length.last().unwrap_or(&0.0)
    }

    /// Mean curvature across the spline.
    pub fn mean_curvature(&self) -> f64 {
        if self.curvature.is_empty() {
            return 0.0;
        }
        self.curvature.iter().copied().sum::<f64>() / self.curvature.len() as f64
    }

    /// Maximum curvature across the spline — a single bend metric.
    pub fn max_curvature(&self) -> f64 {
        self.curvature
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
    }

    /// Evaluate the curved-axis position at an arc-length coordinate
    /// `s ∈ [0, contour_length]`. Out-of-range values are clamped.
    pub fn evaluate(&self, s: f64) -> Point3<f64> {
        if self.arc_length.is_empty() {
            return Point3::origin();
        }
        let s = s.clamp(0.0, self.contour_length());
        // Find the segment whose arc-length window contains s.
        for k in 0..self.segments.len() {
            if s <= self.arc_length[k + 1] + 1e-9 {
                let local = s - self.arc_length[k];
                return self.segments[k].evaluate(local);
            }
        }
        // Fallback: last point.
        *self.axis_points.last().unwrap()
    }
}

/// Build a list of base-pair mid-frames from base-pair centres
/// arranged along the +z axis at the given rise + twist. A convenience
/// for testing the curved-axis fit on idealised B-DNA-like geometry.
pub fn ideal_helix_frames(
    n: usize,
    rise: f64,
    twist_deg: f64,
    radius: f64,
) -> Vec<BaseFrame> {
    let twist = twist_deg.to_radians();
    (0..n)
        .map(|i| {
            let theta = i as f64 * twist;
            let origin = Point3::new(
                radius * theta.cos(),
                radius * theta.sin(),
                i as f64 * rise,
            );
            // Each base pair's z axis is the helix axis (+z).
            // The x axis points radially outward from the helix axis,
            // y completes a right-handed frame.
            let x = Vector3::new(theta.cos(), theta.sin(), 0.0);
            let y = Vector3::new(-theta.sin(), theta.cos(), 0.0);
            let z = Vector3::new(0.0, 0.0, 1.0);
            BaseFrame { origin, x, y, z }
        })
        .collect()
}

/// Build a list of base-pair mid-frames following a circular arc in
/// the xz plane — useful for testing curvature recovery on a known
/// bent duplex. `radius_curve` is the radius of curvature of the
/// axis (large radius → small curvature); `total_arc_angle_deg` is
/// the angle subtended by the bent fragment. The helix twist about
/// the local axis is `twist_deg`.
pub fn ideal_bent_helix_frames(
    n: usize,
    _rise: f64,
    twist_deg: f64,
    helix_radius: f64,
    radius_curve: f64,
    total_arc_angle_deg: f64,
) -> Vec<BaseFrame> {
    // The arc length per step IS the rise here; we use the arc-angle
    // distribution `total_arc_angle_deg / (n-1)` directly so the test
    // can dial the curvature with a single knob. `_rise` is reserved
    // for a future use where arc-length and twist are decoupled.
    let twist = twist_deg.to_radians();
    let arc = total_arc_angle_deg.to_radians();
    let step_arc = arc / (n - 1).max(1) as f64;
    (0..n)
        .map(|i| {
            let theta_bend = i as f64 * step_arc;
            // Axis point on a circular arc in the xz plane.
            let cx = radius_curve - radius_curve * theta_bend.cos();
            let cz = radius_curve * theta_bend.sin();
            // Local axis direction tangent to the arc.
            let tx = theta_bend.sin();
            let tz = theta_bend.cos();
            let axis_dir = Vector3::new(tx, 0.0, tz).normalize();
            // Place the bp origin off the axis by helix_radius along
            // the local "radial" direction (rotating around axis_dir
            // by the helix twist).
            let phi = i as f64 * twist;
            // Perpendicular basis for axis_dir.
            let perp1 = Vector3::new(0.0, 1.0, 0.0);
            let perp2 = axis_dir.cross(&perp1).normalize();
            let radial = perp1 * phi.cos() + perp2 * phi.sin();
            let origin =
                Point3::new(cx, 0.0, cz) + radial * helix_radius;
            BaseFrame {
                origin,
                x: radial.normalize(),
                y: axis_dir.cross(&radial).normalize(),
                z: axis_dir,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_a_straight_axis_along_z() {
        // An ideal helix built around the z-axis: the fitted axis
        // direction must be (anti)parallel to z.
        let centers = ideal_helix_centers(12, 3.4, 34.3, 9.0);
        let axis = fit_helical_axis(&centers).unwrap();
        assert!(axis.direction.z.abs() > 0.99, "axis dir {}", axis.direction);
    }

    #[test]
    fn recovers_canonical_b_dna_parameters() {
        // Build a textbook B-DNA helix: rise 3.4, twist 34.3, radius
        // ~9. The fit should recover those within tolerance.
        let centers = ideal_helix_centers(20, 3.4, 34.3, 9.0);
        let axis = fit_helical_axis(&centers).unwrap();
        assert!(
            (axis.rise_per_bp - 3.4).abs() < 0.05,
            "rise {} not ~3.4",
            axis.rise_per_bp
        );
        assert!(
            (axis.twist_per_bp - 34.3).abs() < 1.0,
            "twist {} not ~34.3",
            axis.twist_per_bp
        );
        assert!(
            (axis.radius - 9.0).abs() < 0.1,
            "radius {} not ~9",
            axis.radius
        );
        // ~10.5 bp per turn for B-DNA.
        assert!(
            (axis.bp_per_turn - 10.5).abs() < 0.5,
            "bp/turn {} not ~10.5",
            axis.bp_per_turn
        );
    }

    #[test]
    fn axis_is_oriented_5_to_3() {
        let centers = ideal_helix_centers(10, 3.0, 36.0, 8.0);
        let axis = fit_helical_axis(&centers).unwrap();
        // direction should point from the first centre toward the
        // last, i.e. +z here.
        assert!(axis.direction.z > 0.0);
    }

    #[test]
    fn distance_to_axis_equals_radius_for_ideal_helix() {
        let centers = ideal_helix_centers(15, 3.4, 34.0, 9.0);
        let axis = fit_helical_axis(&centers).unwrap();
        for c in &centers {
            let d = distance_to_axis(&axis, c);
            assert!((d - 9.0).abs() < 0.1, "centre distance {d} not ~9");
        }
    }

    #[test]
    fn axis_length_spans_the_helix() {
        let centers = ideal_helix_centers(11, 3.4, 34.0, 9.0);
        let axis = fit_helical_axis(&centers).unwrap();
        // 11 centres -> 10 steps -> ~34 A of contour.
        assert!((axis.length - 34.0).abs() < 1.0, "length {}", axis.length);
    }

    #[test]
    fn rejects_too_few_centers() {
        let centers = ideal_helix_centers(2, 3.4, 34.0, 9.0);
        assert!(fit_helical_axis(&centers).is_err());
    }

    // --- curved-axis tests --------------------------------------------

    #[test]
    fn curved_axis_runs_on_straight_b_dna() {
        // A straight B-DNA fragment should fit a near-straight curved
        // axis: curvature ~0 everywhere, and the rise/twist/radius
        // should match the canonical B-DNA values.
        let frames = ideal_helix_frames(15, 3.4, 34.3, 9.0);
        let curved = fit_curved_axis(&frames).unwrap();
        assert_eq!(curved.axis_points.len(), 15);
        assert!(
            (curved.rise_per_bp - 3.4).abs() < 0.1,
            "rise {} not ~3.4",
            curved.rise_per_bp
        );
        assert!(
            (curved.twist_per_bp.abs() - 34.3).abs() < 1.0,
            "twist {} not ~34.3",
            curved.twist_per_bp
        );
        assert!(
            (curved.radius - 9.0).abs() < 0.5,
            "radius {} not ~9",
            curved.radius
        );
        // Curvature should be near zero everywhere on a straight helix.
        let max_k = curved.max_curvature();
        assert!(
            max_k < 0.05,
            "straight B-DNA curvature should be small, got max κ = {max_k}"
        );
    }

    #[test]
    fn curved_axis_recovers_bend() {
        // A bent duplex with a circular-arc axis of curvature radius
        // 100 Å should give a curvature κ ≈ 1/100 = 0.01 Å⁻¹.
        let frames = ideal_bent_helix_frames(15, 3.4, 34.3, 9.0, 100.0, 30.0);
        let curved = fit_curved_axis(&frames).unwrap();
        // The mean curvature should be measurably larger than zero
        // (the straight test asserts < 0.05; a bent helix exceeds it).
        let mean_k = curved.mean_curvature();
        assert!(
            mean_k > 0.001,
            "bent-helix mean curvature should be > 0.001, got {mean_k}"
        );
        // And measurably smaller than the straight-axis bound's
        // maximum — we expect κ ≈ 0.01 Å⁻¹ (1 / 100 Å).
        assert!(
            mean_k < 0.1,
            "bent-helix mean curvature should be < 0.1, got {mean_k}"
        );
    }

    #[test]
    fn curved_axis_passes_through_knots() {
        // The spline must interpolate the axis points exactly at
        // arc-length values arc_length[k].
        let frames = ideal_helix_frames(12, 3.4, 34.3, 9.0);
        let curved = fit_curved_axis(&frames).unwrap();
        for k in 0..curved.axis_points.len() {
            let p = curved.evaluate(curved.arc_length[k]);
            let target = curved.axis_points[k];
            assert!(
                (p - target).norm() < 1e-6,
                "spline does not pass through knot {k}: got {p:?}, target {target:?}"
            );
        }
    }

    #[test]
    fn curved_axis_contour_length_is_sane() {
        let frames = ideal_helix_frames(11, 3.4, 34.3, 9.0);
        let curved = fit_curved_axis(&frames).unwrap();
        // 11 points -> 10 steps, each ~3.4 Å rise -> ~34 Å
        // (could be slightly longer if axis points are not perfectly
        // on the helix axis).
        let l = curved.contour_length();
        assert!(
            (l - 34.0).abs() < 5.0,
            "contour length {l} unexpected for 11-bp B-DNA"
        );
    }

    #[test]
    fn curved_axis_rejects_too_few_frames() {
        let frames = ideal_helix_frames(2, 3.4, 34.3, 9.0);
        assert!(fit_curved_axis(&frames).is_err());
    }

    #[test]
    fn spline_segment_evaluation_is_consistent() {
        // Build a tiny straight-line spline and verify the segment
        // evaluation returns the original points at s=0 and s=h.
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let (segs, _) = natural_cubic_spline(&pts);
        for (k, seg) in segs.iter().enumerate() {
            let p0 = seg.evaluate(0.0);
            let ph = seg.evaluate(seg.h);
            assert!((p0 - pts[k]).norm() < 1e-9);
            assert!((ph - pts[k + 1]).norm() < 1e-9);
        }
    }
}

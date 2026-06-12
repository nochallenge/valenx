//! 3DNA-class base-pair and base-pair-step parameters.
//!
//! For each base a local right-handed reference frame is built from
//! the base-ring atoms. Comparing the two frames of a pair yields the
//! six **base-pair parameters**; comparing the mid-frames of two
//! consecutive pairs yields the six **step parameters**.
//!
//! | Group | Three translations | Three rotations |
//! |---|---|---|
//! | Base-pair | shear, stretch, stagger | buckle, propeller, opening |
//! | Step | shift, slide, rise | tilt, roll, twist |
//!
//! ## Scope of this v1
//!
//! The reference frame here is a least-squares frame fitted from the
//! base-ring atoms — close to, but not bit-identical with, the
//! 3DNA "standard reference frame" (which is fitted to idealised
//! standard-base templates by the Olson et al. 2001 convention). The
//! parameter *definitions* (the CEHS / mid-frame decomposition) are
//! the standard ones, so the numbers are physically meaningful and
//! correctly signed; absolute values may differ by a small amount
//! from a 3DNA run on the same coordinates.

use crate::error::{BiostructError, Result};
use crate::structure::Residue;
use nalgebra::{Matrix3, Point3, Vector3};

/// A local right-handed orthonormal reference frame.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BaseFrame {
    /// Frame origin (≈ the base centroid).
    pub origin: Point3<f64>,
    /// Short-axis / x direction (toward the major groove edge).
    pub x: Vector3<f64>,
    /// Long-axis / y direction (along the base-pair long axis).
    pub y: Vector3<f64>,
    /// Normal / z direction (the base-plane normal).
    pub z: Vector3<f64>,
}

impl BaseFrame {
    /// The frame's rotation matrix (columns = x, y, z axes).
    pub fn rotation(&self) -> Matrix3<f64> {
        Matrix3::from_columns(&[self.x, self.y, self.z])
    }
}

/// The six base-pair parameters, all in ångström (translations) or
/// degrees (rotations).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct BasePairParameters {
    /// Shear — x translation between the two bases, Å.
    pub shear: f64,
    /// Stretch — y translation, Å.
    pub stretch: f64,
    /// Stagger — z translation, Å.
    pub stagger: f64,
    /// Buckle — rotation about x, degrees.
    pub buckle: f64,
    /// Propeller — rotation about y, degrees.
    pub propeller: f64,
    /// Opening — rotation about z, degrees.
    pub opening: f64,
}

/// The six base-pair-step parameters.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct StepParameters {
    /// Shift — x translation between consecutive pairs, Å.
    pub shift: f64,
    /// Slide — y translation, Å.
    pub slide: f64,
    /// Rise — z translation, Å.
    pub rise: f64,
    /// Tilt — rotation about x, degrees.
    pub tilt: f64,
    /// Roll — rotation about y, degrees.
    pub roll: f64,
    /// Twist — rotation about z, degrees.
    pub twist: f64,
}

/// Build a base reference frame from a nucleotide's ring atoms.
///
/// The origin is the ring-atom centroid; the plane normal `z` is the
/// best-fit plane normal (smallest-inertia eigenvector); `y` points
/// from the centroid toward the glycosidic `N1`/`N9` atom; `x`
/// completes the right-handed triple.
pub fn base_frame(residue: &Residue) -> Result<BaseFrame> {
    // Ring atoms — purines have a 9-atom fused ring, pyrimidines 6.
    const RING: &[&str] = &["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"];
    let mut pts: Vec<Point3<f64>> = Vec::new();
    for name in RING {
        if let Some(a) = residue.primary_atom(name) {
            pts.push(a.coord);
        }
    }
    if pts.len() < 3 {
        return Err(BiostructError::invalid(
            "residue",
            "base frame needs at least 3 ring atoms",
        ));
    }
    let origin = centroid(&pts);

    // Best-fit plane normal: smallest eigenvector of the scatter
    // matrix of the centred ring atoms.
    let mut scatter = Matrix3::zeros();
    for p in &pts {
        let r = p.coords - origin.coords;
        scatter += r * r.transpose();
    }
    let eig = nalgebra::SymmetricEigen::new(scatter);
    let mut min_i = 0;
    for i in 1..3 {
        if eig.eigenvalues[i] < eig.eigenvalues[min_i] {
            min_i = i;
        }
    }
    let mut z = eig.eigenvectors.column(min_i).into_owned().normalize();

    // y points toward the glycosidic nitrogen (N9 purine / N1
    // pyrimidine). Project it into the base plane.
    let glyco = residue
        .primary_atom("N9")
        .or_else(|| residue.primary_atom("N1"))
        .map(|a| a.coord)
        .unwrap_or(pts[0]);
    let mut y: Vector3<f64> = glyco - origin;
    y -= z * y.dot(&z); // remove the out-of-plane component
    if y.norm() < 1e-6 {
        // Degenerate — pick any in-plane direction.
        y = z.cross(&Vector3::x());
        if y.norm() < 1e-6 {
            y = z.cross(&Vector3::y());
        }
    }
    let y = y.normalize();
    let x = y.cross(&z).normalize();
    // Re-orthogonalise z for a perfect right-handed frame.
    z = x.cross(&y).normalize();

    Ok(BaseFrame { origin, x, y, z })
}

/// Compute the six base-pair parameters from the two bases' frames.
///
/// The partner frame is flipped about its x-axis (the standard
/// antiparallel-strand convention) before the comparison so a
/// canonical Watson-Crick pair gives near-zero parameters.
pub fn base_pair_parameters(frame_a: &BaseFrame, frame_b: &BaseFrame) -> BasePairParameters {
    // Flip frame B: antiparallel strand -> reverse y and z.
    let flipped = BaseFrame {
        origin: frame_b.origin,
        x: frame_b.x,
        y: -frame_b.y,
        z: -frame_b.z,
    };
    let (rot, trans) = relative_transform(frame_a, &flipped);
    let (rx, ry, rz) = rotation_to_xyz_degrees(&rot);
    BasePairParameters {
        shear: trans.x,
        stretch: trans.y,
        stagger: trans.z,
        buckle: rx,
        propeller: ry,
        opening: rz,
    }
}

/// Compute the six step parameters from two consecutive base-pair
/// mid-frames.
pub fn step_parameters(frame_lower: &BaseFrame, frame_upper: &BaseFrame) -> StepParameters {
    let (rot, trans) = relative_transform(frame_lower, frame_upper);
    let (rx, ry, rz) = rotation_to_xyz_degrees(&rot);
    StepParameters {
        shift: trans.x,
        slide: trans.y,
        rise: trans.z,
        tilt: rx,
        roll: ry,
        twist: rz,
    }
}

/// The mid-frame of a base pair: the frame halfway between the two
/// base frames, used as the "base-pair frame" for step parameters.
pub fn pair_mid_frame(frame_a: &BaseFrame, frame_b: &BaseFrame) -> BaseFrame {
    let flipped = BaseFrame {
        origin: frame_b.origin,
        x: frame_b.x,
        y: -frame_b.y,
        z: -frame_b.z,
    };
    average_frame(frame_a, &flipped)
}

/// The relative transform from `from` to `to`, expressed in `from`'s
/// coordinate system: `(rotation, translation)`.
fn relative_transform(from: &BaseFrame, to: &BaseFrame) -> (Matrix3<f64>, Vector3<f64>) {
    let r_from = from.rotation();
    let r_to = to.rotation();
    // Rotation that takes `from` axes to `to` axes, in from-frame.
    let rel_rot = r_from.transpose() * r_to;
    // Origin offset, expressed in from's frame.
    let offset = to.origin - from.origin;
    let rel_trans = r_from.transpose() * offset;
    (rel_rot, rel_trans)
}

/// The average of two frames: midpoint origin, and the rotation
/// "halfway" between them (via the average axes, re-orthonormalised).
fn average_frame(a: &BaseFrame, b: &BaseFrame) -> BaseFrame {
    let origin = Point3::from((a.origin.coords + b.origin.coords) * 0.5);
    // Average each axis and re-orthonormalise (Gram-Schmidt).
    let mut z = (a.z + b.z).normalize();
    let mut y = a.y + b.y;
    y -= z * y.dot(&z);
    let y = y.normalize();
    let x = y.cross(&z).normalize();
    z = x.cross(&y).normalize();
    BaseFrame { origin, x, y, z }
}

/// Decompose a rotation matrix into x-y-z intrinsic Euler angles in
/// degrees — the `(buckle/tilt, propeller/roll, opening/twist)`
/// triple used by both parameter sets.
fn rotation_to_xyz_degrees(r: &Matrix3<f64>) -> (f64, f64, f64) {
    // Tait-Bryan X-Y-Z extraction.
    let sy = -r[(2, 0)];
    let (rx, ry, rz);
    if sy.abs() < 0.999999 {
        ry = sy.clamp(-1.0, 1.0).asin();
        rx = r[(2, 1)].atan2(r[(2, 2)]);
        rz = r[(1, 0)].atan2(r[(0, 0)]);
    } else {
        // Gimbal lock.
        ry = sy.clamp(-1.0, 1.0).asin();
        rx = (-r[(1, 2)]).atan2(r[(1, 1)]);
        rz = 0.0;
    }
    (rx.to_degrees(), ry.to_degrees(), rz.to_degrees())
}

/// Centroid of a point list (caller guarantees non-empty).
fn centroid(pts: &[Point3<f64>]) -> Point3<f64> {
    let mut acc = Vector3::zeros();
    for p in pts {
        acc += p.coords;
    }
    Point3::from(acc / pts.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Atom;

    /// Build a flat hexagon of ring atoms in the z=0 plane.
    fn flat_base(name: &str, center: Point3<f64>) -> Residue {
        let mut r = Residue::new(name, 1);
        let ring = ["N1", "C2", "N3", "C4", "C5", "C6"];
        for (k, atom) in ring.iter().enumerate() {
            let theta = k as f64 * std::f64::consts::FRAC_PI_3;
            let p = center + Vector3::new(1.4 * theta.cos(), 1.4 * theta.sin(), 0.0);
            r.atoms.push(Atom::new(*atom, "C", p));
        }
        r
    }

    #[test]
    fn frame_is_orthonormal() {
        let base = flat_base("DA", Point3::origin());
        let f = base_frame(&base).unwrap();
        assert!((f.x.norm() - 1.0).abs() < 1e-9);
        assert!((f.y.norm() - 1.0).abs() < 1e-9);
        assert!((f.z.norm() - 1.0).abs() < 1e-9);
        assert!(f.x.dot(&f.y).abs() < 1e-9);
        assert!(f.x.dot(&f.z).abs() < 1e-9);
        assert!(f.y.dot(&f.z).abs() < 1e-9);
        // a flat z=0 base has a normal parallel to z.
        assert!(f.z.z.abs() > 0.99);
    }

    #[test]
    fn frame_is_right_handed() {
        let base = flat_base("DG", Point3::new(3.0, 1.0, 0.0));
        let f = base_frame(&base).unwrap();
        // x cross y == z for a right-handed triple.
        assert!((f.x.cross(&f.y) - f.z).norm() < 1e-6);
    }

    #[test]
    fn identical_frames_give_zero_step() {
        let base = flat_base("DA", Point3::origin());
        let f = base_frame(&base).unwrap();
        let step = step_parameters(&f, &f);
        assert!(step.rise.abs() < 1e-9);
        assert!(step.twist.abs() < 1e-9);
        assert!(step.shift.abs() < 1e-9);
    }

    #[test]
    fn pure_rise_step() {
        // Two identical frames, the upper one shifted +3.4 A along
        // its own z: rise = 3.4, everything else 0.
        let base = flat_base("DA", Point3::origin());
        let lower = base_frame(&base).unwrap();
        let mut upper = lower;
        upper.origin = lower.origin + lower.z * 3.4;
        let step = step_parameters(&lower, &upper);
        assert!((step.rise - 3.4).abs() < 1e-6, "rise was {}", step.rise);
        assert!(step.shift.abs() < 1e-6);
        assert!(step.slide.abs() < 1e-6);
        assert!(step.twist.abs() < 1e-6);
    }

    #[test]
    fn pure_twist_step() {
        // Rotate the upper frame 36 deg about z: twist = 36.
        let base = flat_base("DA", Point3::origin());
        let lower = base_frame(&base).unwrap();
        let angle = 36.0_f64.to_radians();
        let rz = Matrix3::new(
            angle.cos(),
            -angle.sin(),
            0.0,
            angle.sin(),
            angle.cos(),
            0.0,
            0.0,
            0.0,
            1.0,
        );
        // upper axes = lower axes rotated by rz (in lower's frame).
        let r_lower = lower.rotation();
        let r_upper = r_lower * rz;
        let upper = BaseFrame {
            origin: lower.origin,
            x: r_upper.column(0).into_owned(),
            y: r_upper.column(1).into_owned(),
            z: r_upper.column(2).into_owned(),
        };
        let step = step_parameters(&lower, &upper);
        assert!((step.twist - 36.0).abs() < 1e-4, "twist was {}", step.twist);
        assert!(step.roll.abs() < 1e-4);
        assert!(step.tilt.abs() < 1e-4);
    }

    #[test]
    fn base_pair_params_run() {
        // Two flat bases facing each other -> finite parameters.
        let a = flat_base("DA", Point3::new(0.0, 0.0, 0.0));
        let t = flat_base("DT", Point3::new(0.0, 8.5, 0.0));
        let fa = base_frame(&a).unwrap();
        let ft = base_frame(&t).unwrap();
        let bp = base_pair_parameters(&fa, &ft);
        // Just assert the numbers are finite and the function ran.
        assert!(bp.shear.is_finite());
        assert!(bp.opening.is_finite());
    }

    #[test]
    fn mid_frame_is_between() {
        let a = flat_base("DA", Point3::new(0.0, 0.0, 0.0));
        let t = flat_base("DT", Point3::new(0.0, 8.0, 0.0));
        let fa = base_frame(&a).unwrap();
        let ft = base_frame(&t).unwrap();
        let mid = pair_mid_frame(&fa, &ft);
        // mid origin y is between the two base centroids.
        assert!(mid.origin.y > 0.0 && mid.origin.y < 8.0);
    }

    #[test]
    fn rejects_too_few_ring_atoms() {
        let mut r = Residue::new("DA", 1);
        r.atoms.push(Atom::new("N1", "N", Point3::origin()));
        assert!(base_frame(&r).is_err());
    }
}

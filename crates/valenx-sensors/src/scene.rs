//! Analytic ray-cast scene primitives for the simulated [`crate::lidar`].
//!
//! A [`Scene`] is a small collection of closed-form surfaces — infinite
//! [`Plane`]s, [`Sphere`]s, and [`Triangle`]s — each of which can be intersected
//! with a [`Ray`] analytically. There is no acceleration structure and no mesh
//! loader: the point of this module is a *ground-truth* world whose ray–surface
//! ranges are exact (or exact to floating-point), so the LiDAR model can be
//! pinned against hand-computed distances in tests.
//!
//! Conventions:
//! - A [`Ray`] has an `origin` and a **unit** `direction`; [`Ray::new`]
//!   normalises the direction (and rejects a zero direction).
//! - An intersection returns the ray parameter `t` (the distance along the unit
//!   direction, i.e. metres) of the *nearest* hit strictly in front of the
//!   origin (`t > 0`), or `None` for a miss.

use nalgebra::Vector3;

use crate::error::SensorError;

/// A small tolerance used to reject self-intersections and treat near-parallel
/// rays as misses.
const EPS: f64 = 1e-9;

/// A semi-infinite ray with a unit direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    /// Ray origin (m).
    pub origin: Vector3<f64>,
    /// Ray direction, **unit length** (enforced by [`Ray::new`]).
    pub direction: Vector3<f64>,
}

impl Ray {
    /// Build a ray, normalising `direction`.
    ///
    /// # Errors
    /// Returns [`SensorError::DegenerateGeometry`] if `direction` has
    /// (near-)zero length and so has no defined heading.
    pub fn new(origin: Vector3<f64>, direction: Vector3<f64>) -> Result<Self, SensorError> {
        let norm = direction.norm();
        if norm < EPS {
            return Err(SensorError::DegenerateGeometry(
                "ray direction has zero length".into(),
            ));
        }
        Ok(Self {
            origin,
            direction: direction / norm,
        })
    }

    /// The point at parameter `t` along the ray: `origin + t · direction`.
    #[must_use]
    pub fn at(&self, t: f64) -> Vector3<f64> {
        self.origin + self.direction * t
    }
}

/// An infinite plane `{ x : n · (x − point) = 0 }`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// A point on the plane (m).
    pub point: Vector3<f64>,
    /// The plane normal, **unit length** (enforced by [`Plane::new`]).
    pub normal: Vector3<f64>,
}

impl Plane {
    /// Build a plane, normalising `normal`.
    ///
    /// # Errors
    /// Returns [`SensorError::DegenerateGeometry`] if `normal` is (near-)zero.
    pub fn new(point: Vector3<f64>, normal: Vector3<f64>) -> Result<Self, SensorError> {
        let n = normal.norm();
        if n < EPS {
            return Err(SensorError::DegenerateGeometry(
                "plane normal has zero length".into(),
            ));
        }
        Ok(Self {
            point,
            normal: normal / n,
        })
    }

    /// Nearest ray–plane intersection distance `t > 0`, or `None`.
    ///
    /// Solves `n · (o + t·d − p) = 0 ⇒ t = n·(p − o) / (n·d)`. A ray parallel to
    /// the plane (`|n·d| < EPS`) misses; a hit behind the origin (`t ≤ EPS`) is
    /// not returned.
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        let denom = self.normal.dot(&ray.direction);
        if denom.abs() < EPS {
            return None; // parallel
        }
        let t = self.normal.dot(&(self.point - ray.origin)) / denom;
        (t > EPS).then_some(t)
    }
}

/// A sphere of `radius` about `center`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    /// Sphere centre (m).
    pub center: Vector3<f64>,
    /// Sphere radius (m), strictly positive (enforced by [`Sphere::new`]).
    pub radius: f64,
}

impl Sphere {
    /// Build a sphere.
    ///
    /// # Errors
    /// Returns [`SensorError::DegenerateGeometry`] if `radius` is not strictly
    /// positive (or not finite).
    pub fn new(center: Vector3<f64>, radius: f64) -> Result<Self, SensorError> {
        if !(radius.is_finite() && radius > 0.0) {
            return Err(SensorError::DegenerateGeometry(format!(
                "sphere radius must be > 0, got {radius}"
            )));
        }
        Ok(Self { center, radius })
    }

    /// Nearest ray–sphere intersection distance `t > 0`, or `None`.
    ///
    /// Substituting the ray into `|x − c|² = r²` gives a quadratic in `t` with
    /// `a = 1` (unit direction). Returns the smaller positive root (the near
    /// face); if the origin is inside the sphere the far root (the exit point) is
    /// returned.
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        let oc = ray.origin - self.center;
        let b = 2.0 * oc.dot(&ray.direction);
        let c = oc.norm_squared() - self.radius * self.radius;
        let disc = b * b - 4.0 * c;
        if disc < 0.0 {
            return None; // no real intersection
        }
        let sqrt_disc = disc.sqrt();
        let t0 = 0.5 * (-b - sqrt_disc);
        let t1 = 0.5 * (-b + sqrt_disc);
        if t0 > EPS {
            Some(t0)
        } else if t1 > EPS {
            Some(t1)
        } else {
            None
        }
    }
}

/// A single triangle, intersected with the **Möller–Trumbore** algorithm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle {
    /// First vertex (m).
    pub v0: Vector3<f64>,
    /// Second vertex (m).
    pub v1: Vector3<f64>,
    /// Third vertex (m).
    pub v2: Vector3<f64>,
}

impl Triangle {
    /// Build a triangle, rejecting a degenerate (zero-area / collinear) one.
    ///
    /// # Errors
    /// Returns [`SensorError::DegenerateGeometry`] if the three vertices are
    /// collinear (the edge cross product has near-zero length).
    pub fn new(v0: Vector3<f64>, v1: Vector3<f64>, v2: Vector3<f64>) -> Result<Self, SensorError> {
        let cross = (v1 - v0).cross(&(v2 - v0));
        if cross.norm() < EPS {
            return Err(SensorError::DegenerateGeometry(
                "triangle vertices are collinear (zero area)".into(),
            ));
        }
        Ok(Self { v0, v1, v2 })
    }

    /// Nearest ray–triangle intersection distance `t > 0`, or `None`, via the
    /// Möller–Trumbore intersection test.
    ///
    /// The test is **double-sided** (a back-facing hit counts), since a LiDAR
    /// beam ranges a surface regardless of which way its winding faces.
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        let edge1 = self.v1 - self.v0;
        let edge2 = self.v2 - self.v0;
        let pvec = ray.direction.cross(&edge2);
        let det = edge1.dot(&pvec);
        if det.abs() < EPS {
            return None; // ray parallel to the triangle plane
        }
        let inv_det = 1.0 / det;
        let tvec = ray.origin - self.v0;
        let u = tvec.dot(&pvec) * inv_det;
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let qvec = tvec.cross(&edge1);
        let v = ray.direction.dot(&qvec) * inv_det;
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = edge2.dot(&qvec) * inv_det;
        (t > EPS).then_some(t)
    }
}

/// One ray-castable surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface {
    /// An infinite plane.
    Plane(Plane),
    /// A sphere.
    Sphere(Sphere),
    /// A single triangle.
    Triangle(Triangle),
}

impl Surface {
    /// Nearest intersection distance `t > 0` of this surface with `ray`.
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        match self {
            Surface::Plane(p) => p.intersect(ray),
            Surface::Sphere(s) => s.intersect(ray),
            Surface::Triangle(t) => t.intersect(ray),
        }
    }
}

/// A collection of analytic surfaces forming the world a LiDAR ranges against.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Scene {
    /// The surfaces in the scene.
    pub surfaces: Vec<Surface>,
}

impl Scene {
    /// An empty scene (every ray misses).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a surface (builder style).
    #[must_use]
    pub fn with(mut self, surface: Surface) -> Self {
        self.surfaces.push(surface);
        self
    }

    /// Add a plane.
    pub fn push_plane(&mut self, plane: Plane) {
        self.surfaces.push(Surface::Plane(plane));
    }

    /// Add a sphere.
    pub fn push_sphere(&mut self, sphere: Sphere) {
        self.surfaces.push(Surface::Sphere(sphere));
    }

    /// Add a triangle.
    pub fn push_triangle(&mut self, triangle: Triangle) {
        self.surfaces.push(Surface::Triangle(triangle));
    }

    /// Cast `ray` against every surface and return the **nearest** hit distance
    /// `t > 0`, or `None` if the ray misses the whole scene.
    #[must_use]
    pub fn cast(&self, ray: &Ray) -> Option<f64> {
        let mut nearest: Option<f64> = None;
        for s in &self.surfaces {
            if let Some(t) = s.intersect(ray) {
                nearest = Some(match nearest {
                    Some(best) if best <= t => best,
                    _ => t,
                });
            }
        }
        nearest
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn ray_rejects_zero_direction() {
        assert!(Ray::new(v(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0)).is_err());
    }

    #[test]
    fn ray_normalises_direction() {
        let r = Ray::new(v(0.0, 0.0, 0.0), v(0.0, 0.0, 5.0)).unwrap();
        assert!((r.direction.norm() - 1.0).abs() < 1e-12);
        assert!((r.at(3.0) - v(0.0, 0.0, 3.0)).norm() < 1e-12);
    }

    #[test]
    fn plane_hit_is_exact_distance() {
        // Plane x = 10 (normal +x), ray from origin along +x ⇒ t = 10.
        let plane = Plane::new(v(10.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let t = plane.intersect(&ray).unwrap();
        assert!((t - 10.0).abs() < 1e-12, "t = {t}");
    }

    #[test]
    fn plane_behind_and_parallel_miss() {
        let plane = Plane::new(v(10.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        // Looking away from the plane: behind ⇒ miss.
        let away = Ray::new(v(0.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap();
        assert!(plane.intersect(&away).is_none());
        // Parallel to the plane ⇒ miss.
        let parallel = Ray::new(v(0.0, 0.0, 0.0), v(0.0, 1.0, 0.0)).unwrap();
        assert!(plane.intersect(&parallel).is_none());
    }

    #[test]
    fn sphere_near_face_distance() {
        // Sphere radius 1 centred at x = 5; ray along +x hits the near face at
        // x = 4 ⇒ t = 4.
        let sphere = Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap();
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let t = sphere.intersect(&ray).unwrap();
        assert!((t - 4.0).abs() < 1e-12, "t = {t}");
    }

    #[test]
    fn sphere_tangent_and_miss() {
        let sphere = Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap();
        // Offset by 2 in y (> radius) ⇒ miss.
        let miss = Ray::new(v(0.0, 2.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        assert!(sphere.intersect(&miss).is_none());
    }

    #[test]
    fn sphere_from_inside_returns_exit() {
        let sphere = Sphere::new(v(0.0, 0.0, 0.0), 2.0).unwrap();
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let t = sphere.intersect(&ray).unwrap();
        assert!((t - 2.0).abs() < 1e-12, "exit t = {t}");
    }

    #[test]
    fn triangle_moller_trumbore_distance() {
        // Triangle in the plane z = 3 spanning the origin region; a +z ray
        // through the centroid hits at t = 3.
        let tri = Triangle::new(v(-1.0, -1.0, 3.0), v(1.0, -1.0, 3.0), v(0.0, 1.0, 3.0)).unwrap();
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let t = tri.intersect(&ray).unwrap();
        assert!((t - 3.0).abs() < 1e-12, "t = {t}");
    }

    #[test]
    fn triangle_outside_barycentric_misses() {
        let tri = Triangle::new(v(-1.0, -1.0, 3.0), v(1.0, -1.0, 3.0), v(0.0, 1.0, 3.0)).unwrap();
        // Aim well outside the triangle in x.
        let ray = Ray::new(v(10.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        assert!(tri.intersect(&ray).is_none());
    }

    #[test]
    fn degenerate_geometry_rejected() {
        assert!(Sphere::new(v(0.0, 0.0, 0.0), 0.0).is_err());
        assert!(Sphere::new(v(0.0, 0.0, 0.0), -1.0).is_err());
        assert!(Plane::new(v(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0)).is_err());
        // Collinear triangle.
        assert!(Triangle::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(2.0, 0.0, 0.0)).is_err());
    }

    #[test]
    fn scene_returns_nearest_hit() {
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(10.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap());
        scene.push_sphere(Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap()); // near face t=4
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let t = scene.cast(&ray).unwrap();
        assert!((t - 4.0).abs() < 1e-12, "nearest = {t}");
    }

    #[test]
    fn empty_scene_misses() {
        let scene = Scene::new();
        let ray = Ray::new(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        assert!(scene.cast(&ray).is_none());
    }
}

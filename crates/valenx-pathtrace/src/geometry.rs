//! Ray, axis-aligned bounding box, and the triangle primitive with the
//! Möller-Trumbore intersection test.
//!
//! These are the geometric atoms the [`crate::bvh`] indexes and the
//! [`crate::tracer`] casts against. Everything is `f32` and `Copy`.

use crate::math::Vec3;

/// A ray — an origin plus a unit direction, traced over the parametric
/// half-line `origin + t·direction` for `t ∈ [t_min, t_max]`.
///
/// The `t_min` / `t_max` bounds let a caller exclude self-intersection
/// (a small positive `t_min`) and limit a shadow ray to the distance of
/// the light without paying for geometry beyond it.
#[derive(Clone, Copy, Debug)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: Vec3,
    /// Unit ray direction.
    pub direction: Vec3,
    /// Reciprocal of `direction`, precomputed once so the slab AABB
    /// test multiplies instead of dividing per box.
    pub inv_dir: Vec3,
}

impl Ray {
    /// Build a ray from an origin and a direction.
    ///
    /// The direction is **not** normalised here — callers pass a unit
    /// direction (the path tracer always does). `inv_dir` is computed
    /// from whatever is supplied; a zero component yields an infinite
    /// reciprocal, which the slab test in [`Aabb::hit`] handles
    /// correctly (the IEEE `inf`/`-inf` ordering is exactly what the
    /// branchless slab method relies on).
    #[inline]
    pub fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray {
            origin,
            direction,
            inv_dir: Vec3 {
                x: 1.0 / direction.x,
                y: 1.0 / direction.y,
                z: 1.0 / direction.z,
            },
        }
    }

    /// The point at parameter `t` along the ray.
    #[inline]
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin.add(self.direction.scale(t))
    }
}

/// An axis-aligned bounding box — the BVH node volume.
#[derive(Clone, Copy, Debug)]
pub struct Aabb {
    /// Minimum corner (componentwise smallest).
    pub min: Vec3,
    /// Maximum corner (componentwise largest).
    pub max: Vec3,
}

impl Aabb {
    /// The "empty" box — `min` at `+∞`, `max` at `−∞`. Unioning any
    /// point into it yields exactly that point's degenerate box, so it
    /// is the correct identity for an accumulating bounds fold.
    pub fn empty() -> Aabb {
        Aabb {
            min: Vec3::splat(f32::INFINITY),
            max: Vec3::splat(f32::NEG_INFINITY),
        }
    }

    /// Grow the box to include the point `p`.
    #[inline]
    pub fn expand_point(&mut self, p: Vec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// Grow the box to include another box.
    #[inline]
    pub fn expand_box(&mut self, o: &Aabb) {
        self.min = self.min.min(o.min);
        self.max = self.max.max(o.max);
    }

    /// The box centre.
    #[inline]
    pub fn centroid(&self) -> Vec3 {
        self.min.add(self.max).scale(0.5)
    }

    /// The diagonal extent (`max − min`).
    #[inline]
    pub fn extent(&self) -> Vec3 {
        self.max.sub(self.min)
    }

    /// Surface area — the cost metric of the BVH's surface-area
    /// heuristic. An empty box returns 0.
    #[inline]
    pub fn surface_area(&self) -> f32 {
        let e = self.extent();
        if e.x < 0.0 || e.y < 0.0 || e.z < 0.0 {
            return 0.0;
        }
        2.0 * (e.x * e.y + e.y * e.z + e.z * e.x)
    }

    /// Slab ray-box intersection test.
    ///
    /// Returns `true` if the ray's `[t_min, t_max]` interval overlaps
    /// the box. This is the branchless slab method: for each axis the
    /// ray's entry / exit parameters into the pair of parallel planes
    /// are computed by a multiply with the precomputed `inv_dir`, the
    /// per-axis intervals are intersected, and an overlap survives iff
    /// `t_enter ≤ t_exit`. A zero `direction` component produces an
    /// `inf` reciprocal, and the `min`/`max` of `±inf` collapses the
    /// slab to "the ray is inside or outside this pair of planes for
    /// all t" — the correct degenerate behaviour.
    #[inline]
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for axis in 0..3 {
            let inv = ray.inv_dir.axis(axis);
            let o = ray.origin.axis(axis);
            let mut t0 = (self.min.axis(axis) - o) * inv;
            let mut t1 = (self.max.axis(axis) - o) * inv;
            if inv < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            tmin = tmin.max(t0);
            tmax = tmax.min(t1);
            if tmax < tmin {
                return false;
            }
        }
        true
    }
}

/// A shading triangle — three world-space vertices, three per-vertex
/// normals, and an index into the scene's material table.
///
/// Per-vertex normals let the integrator interpolate a smooth shading
/// normal across the face (so a coarse mesh still shades smoothly). The
/// geometric normal (the face's true plane normal) is recomputed on the
/// fly for the light-leak / back-face guards.
#[derive(Clone, Copy, Debug)]
pub struct Triangle {
    /// First vertex.
    pub v0: Vec3,
    /// Second vertex.
    pub v1: Vec3,
    /// Third vertex.
    pub v2: Vec3,
    /// Shading normal at `v0`.
    pub n0: Vec3,
    /// Shading normal at `v1`.
    pub n1: Vec3,
    /// Shading normal at `v2`.
    pub n2: Vec3,
    /// Index into the [`crate::scene::Scene`] material table.
    pub material: usize,
}

/// A successful ray-triangle intersection.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    /// Ray parameter at the hit point.
    pub t: f32,
    /// World-space hit position.
    pub position: Vec3,
    /// Interpolated, normalised shading normal (already flipped to
    /// oppose the incoming ray so it always faces the viewer).
    pub normal: Vec3,
    /// Geometric (true face) normal, same hemisphere as `normal`.
    pub geo_normal: Vec3,
    /// Material index of the triangle that was hit.
    pub material: usize,
    /// True if the ray struck the triangle's back side — the shading
    /// normal was flipped to compensate.
    pub back_face: bool,
}

impl Triangle {
    /// Build a triangle with explicit per-vertex normals.
    pub fn new(verts: [Vec3; 3], normals: [Vec3; 3], material: usize) -> Triangle {
        Triangle {
            v0: verts[0],
            v1: verts[1],
            v2: verts[2],
            n0: normals[0],
            n1: normals[1],
            n2: normals[2],
            material,
        }
    }

    /// Build a flat-shaded triangle — all three vertex normals set to
    /// the geometric face normal.
    pub fn flat(verts: [Vec3; 3], material: usize) -> Triangle {
        let geo = verts[1]
            .sub(verts[0])
            .cross(verts[2].sub(verts[0]))
            .normalized()
            .unwrap_or(Vec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            });
        Triangle::new(verts, [geo; 3], material)
    }

    /// The geometric (true plane) normal of the face, normalised.
    ///
    /// Falls back to `+Z` for a degenerate (zero-area) triangle so the
    /// caller never sees a NaN normal.
    #[inline]
    pub fn geometric_normal(&self) -> Vec3 {
        self.v1
            .sub(self.v0)
            .cross(self.v2.sub(self.v0))
            .normalized()
            .unwrap_or(Vec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            })
    }

    /// Twice the triangle area — `|edge1 × edge2|`.
    #[inline]
    pub fn double_area(&self) -> f32 {
        self.v1.sub(self.v0).cross(self.v2.sub(self.v0)).length()
    }

    /// The triangle's axis-aligned bounding box.
    pub fn bounds(&self) -> Aabb {
        let mut b = Aabb::empty();
        b.expand_point(self.v0);
        b.expand_point(self.v1);
        b.expand_point(self.v2);
        b
    }

    /// The triangle's centroid (used as the BVH split key).
    #[inline]
    pub fn centroid(&self) -> Vec3 {
        self.v0.add(self.v1).add(self.v2).scale(1.0 / 3.0)
    }

    /// Möller-Trumbore ray-triangle intersection.
    ///
    /// Returns the [`Hit`] if the ray strikes this triangle at a
    /// parameter inside `(t_min, t_max)`, else `None`.
    ///
    /// # Method
    ///
    /// The classic 1997 Möller-Trumbore algorithm. The ray is solved
    /// against the triangle's barycentric plane without precomputing
    /// the plane equation:
    ///
    /// ```text
    ///   pvec = direction × edge2
    ///   det  = edge1 · pvec
    ///   u    = (origin − v0) · pvec / det
    ///   qvec = (origin − v0) × edge1
    ///   v    = direction · qvec / det
    ///   t    = edge2 · qvec / det
    /// ```
    ///
    /// A hit requires `u ≥ 0`, `v ≥ 0`, `u + v ≤ 1` (the barycentric
    /// inside test) and `t ∈ (t_min, t_max)`. A `|det|` near zero means
    /// the ray is parallel to the triangle plane — no hit.
    ///
    /// The shading normal is the barycentric blend of the three vertex
    /// normals; both it and the geometric normal are flipped if the ray
    /// hit the back face so they always oppose the incoming direction.
    #[inline]
    pub fn intersect(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
        const PARALLEL_EPS: f32 = 1e-9;
        let edge1 = self.v1.sub(self.v0);
        let edge2 = self.v2.sub(self.v0);
        let pvec = ray.direction.cross(edge2);
        let det = edge1.dot(pvec);
        // Ray parallel to the triangle plane → no intersection. Note we
        // do *not* cull back faces here (a path tracer needs both sides
        // — interiors, glass, double-sided emitters).
        if det.abs() < PARALLEL_EPS {
            return None;
        }
        let inv_det = 1.0 / det;
        let tvec = ray.origin.sub(self.v0);
        let u = tvec.dot(pvec) * inv_det;
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let qvec = tvec.cross(edge1);
        let v = ray.direction.dot(qvec) * inv_det;
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = edge2.dot(qvec) * inv_det;
        if t <= t_min || t >= t_max {
            return None;
        }

        // Barycentric weights: w for v0, u for v1, v for v2.
        let w = 1.0 - u - v;
        let mut normal = self
            .n0
            .scale(w)
            .add(self.n1.scale(u))
            .add(self.n2.scale(v))
            .normalized()
            .unwrap_or_else(|| self.geometric_normal());
        let mut geo = self.geometric_normal();

        // Flip both normals to oppose the incoming ray so the shading
        // code always has a viewer-facing normal. `det < 0` means the
        // ray hit the back face.
        let back_face = ray.direction.dot(geo) > 0.0;
        if back_face {
            normal = normal.neg();
            geo = geo.neg();
        }

        Some(Hit {
            t,
            position: ray.at(t),
            normal,
            geo_normal: geo,
            material: self.material,
            back_face,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    fn unit_triangle() -> Triangle {
        // A triangle in the z = 0 plane, facing +Z.
        Triangle::flat(
            [
                vec3(0.0, 0.0, 0.0),
                vec3(1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
            ],
            0,
        )
    }

    #[test]
    fn ray_hits_triangle_centre() {
        let tri = unit_triangle();
        // Ray straight down at the triangle's interior.
        let ray = Ray::new(vec3(0.25, 0.25, 1.0), vec3(0.0, 0.0, -1.0));
        let hit = tri.intersect(&ray, 1e-4, 1e30).expect("centre ray hits");
        assert!((hit.t - 1.0).abs() < 1e-5, "t should be 1, got {}", hit.t);
        // The hit normal must point back up at the ray (+Z).
        assert!(hit.normal.z > 0.9, "normal should face the ray");
    }

    #[test]
    fn ray_misses_outside_the_triangle() {
        let tri = unit_triangle();
        // Aimed at a point outside the barycentric simplex.
        let ray = Ray::new(vec3(0.9, 0.9, 1.0), vec3(0.0, 0.0, -1.0));
        assert!(tri.intersect(&ray, 1e-4, 1e30).is_none());
    }

    #[test]
    fn parallel_ray_does_not_intersect() {
        let tri = unit_triangle();
        // A ray travelling in the triangle's own plane.
        let ray = Ray::new(vec3(-1.0, 0.25, 0.0), vec3(1.0, 0.0, 0.0));
        assert!(tri.intersect(&ray, 1e-4, 1e30).is_none());
    }

    #[test]
    fn back_face_hit_is_flagged_and_normal_flipped() {
        let tri = unit_triangle();
        // Ray coming from below (−Z), travelling up. It strikes the
        // back of a +Z-facing triangle.
        let ray = Ray::new(vec3(0.25, 0.25, -1.0), vec3(0.0, 0.0, 1.0));
        let hit = tri.intersect(&ray, 1e-4, 1e30).expect("back ray hits");
        assert!(hit.back_face, "should report a back-face hit");
        // The shading normal was flipped to oppose the ray → faces −Z.
        assert!(hit.normal.z < -0.9, "back-face normal should be flipped");
    }

    #[test]
    fn t_range_excludes_hits_outside_the_interval() {
        let tri = unit_triangle();
        let ray = Ray::new(vec3(0.25, 0.25, 1.0), vec3(0.0, 0.0, -1.0));
        // The hit is at t = 1; a [t_min, t_max] that excludes it misses.
        assert!(tri.intersect(&ray, 1e-4, 0.5).is_none(), "t_max clips");
        assert!(tri.intersect(&ray, 2.0, 1e30).is_none(), "t_min clips");
    }

    #[test]
    fn aabb_slab_test_hits_and_misses() {
        let b = Aabb {
            min: vec3(-1.0, -1.0, -1.0),
            max: vec3(1.0, 1.0, 1.0),
        };
        // A ray aimed through the box.
        let through = Ray::new(vec3(0.0, 0.0, -5.0), vec3(0.0, 0.0, 1.0));
        assert!(b.hit(&through, 0.0, 1e30));
        // A ray that passes well to the side.
        let beside = Ray::new(vec3(5.0, 5.0, -5.0), vec3(0.0, 0.0, 1.0));
        assert!(!b.hit(&beside, 0.0, 1e30));
    }

    #[test]
    fn aabb_surface_area_of_unit_cube_is_six() {
        let b = Aabb {
            min: Vec3::ZERO,
            max: Vec3::ONE,
        };
        assert!((b.surface_area() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn triangle_bounds_enclose_all_vertices() {
        let tri = unit_triangle();
        let b = tri.bounds();
        assert_eq!(b.min, vec3(0.0, 0.0, 0.0));
        assert_eq!(b.max, vec3(1.0, 1.0, 0.0));
    }
}

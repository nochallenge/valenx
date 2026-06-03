//! Small `f32` 3-vector arithmetic used throughout the path tracer.
//!
//! The renderer keeps a dependency-free `Vec3` rather than reaching for
//! `nalgebra` so the hot inner loop (millions of ray-triangle tests)
//! has a flat, inlinable, `Copy` value type with no generic machinery.
//! Every operation here is a tiny `#[inline]` free function.

/// A 3-component single-precision vector — a point, a direction, or a
/// linear-RGB colour depending on context.
///
/// `f32` is deliberate: a path tracer is bandwidth-bound on the BVH
/// traversal and the framebuffer, and `f32` halves both. The numeric
/// error of `f32` ray-triangle intersection is well below a pixel for
/// any scene that fits in a sensible world scale; the shadow-ray
/// `EPSILON` offset in [`crate::tracer`] absorbs it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3 {
    /// X component.
    pub x: f32,
    /// Y component.
    pub y: f32,
    /// Z component.
    pub z: f32,
}

/// Shorthand constructor — `vec3(x, y, z)`.
#[inline]
pub fn vec3(x: f32, y: f32, z: f32) -> Vec3 {
    Vec3 { x, y, z }
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Vec3 = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// The all-ones vector (white, when read as a colour).
    pub const ONE: Vec3 = Vec3 {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };

    /// Build a vector with every component set to `s`.
    #[inline]
    pub fn splat(s: f32) -> Vec3 {
        Vec3 { x: s, y: s, z: s }
    }

    /// Build a vector from a `[f32; 3]` array (the colour layout the
    /// rest of `valenx-render-bridge` uses).
    #[inline]
    pub fn from_array(a: [f32; 3]) -> Vec3 {
        Vec3 {
            x: a[0],
            y: a[1],
            z: a[2],
        }
    }

    /// Convert to a `[f32; 3]` array.
    #[inline]
    pub fn to_array(self) -> [f32; 3] {
        [self.x, self.y, self.z]
    }

    /// Component-wise add.
    ///
    /// `Vec3` keeps inherent `add` / `sub` / `mul` / `neg` methods
    /// rather than `impl`ing the `std::ops` traits: the path tracer's
    /// hot loop reads more clearly as an explicit chain
    /// (`a.add(b).scale(s)`) than as operator soup, and an inherent
    /// method is unambiguous about taking `self` by value. The
    /// `should_implement_trait` lint is therefore allowed here by
    /// design.
    #[inline]
    #[allow(clippy::should_implement_trait)]
    pub fn add(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + o.x,
            y: self.y + o.y,
            z: self.z + o.z,
        }
    }

    /// Component-wise subtract.
    #[inline]
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - o.x,
            y: self.y - o.y,
            z: self.z - o.z,
        }
    }

    /// Component-wise multiply (used for colour × colour throughput).
    #[inline]
    #[allow(clippy::should_implement_trait)]
    pub fn mul(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.x * o.x,
            y: self.y * o.y,
            z: self.z * o.z,
        }
    }

    /// Scale by a scalar.
    #[inline]
    pub fn scale(self, s: f32) -> Vec3 {
        Vec3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    /// Unary negation.
    #[inline]
    #[allow(clippy::should_implement_trait)]
    pub fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Vec3) -> f32 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    /// Squared length — cheaper than [`Self::length`] when only a
    /// comparison is needed.
    #[inline]
    pub fn length_sq(self) -> f32 {
        self.dot(self)
    }

    /// Euclidean length.
    #[inline]
    pub fn length(self) -> f32 {
        self.length_sq().sqrt()
    }

    /// Return the unit vector, or `None` if the vector is degenerate
    /// (length below `1e-12`). The `Option` makes the caller decide
    /// what a zero direction means rather than silently returning NaN.
    #[inline]
    pub fn normalized(self) -> Option<Vec3> {
        let len = self.length();
        if len < 1e-12 {
            None
        } else {
            Some(self.scale(1.0 / len))
        }
    }

    /// Largest of the three components — used to test colour throughput
    /// for the Russian-roulette survival probability.
    #[inline]
    pub fn max_component(self) -> f32 {
        self.x.max(self.y).max(self.z)
    }

    /// Component-wise minimum of two vectors (AABB corner reduction).
    #[inline]
    pub fn min(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.x.min(o.x),
            y: self.y.min(o.y),
            z: self.z.min(o.z),
        }
    }

    /// Component-wise maximum of two vectors (AABB corner reduction).
    #[inline]
    pub fn max(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.x.max(o.x),
            y: self.y.max(o.y),
            z: self.z.max(o.z),
        }
    }

    /// Index a component by axis (`0 = x`, `1 = y`, `2 = z`). Used by
    /// the BVH builder to sort centroids on the split axis.
    #[inline]
    pub fn axis(self, axis: usize) -> f32 {
        match axis {
            0 => self.x,
            1 => self.y,
            _ => self.z,
        }
    }

    /// True if every component is finite (no NaN, no infinity).
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Reflect `self` (treated as an incident direction) about the unit
    /// normal `n`: `r = d − 2·(d·n)·n`.
    #[inline]
    pub fn reflect(self, n: Vec3) -> Vec3 {
        self.sub(n.scale(2.0 * self.dot(n)))
    }
}

/// Build an orthonormal `(tangent, bitangent)` basis spanning the plane
/// perpendicular to the unit vector `n`.
///
/// Uses the branch-light Duff et al. (2017) construction —
/// "Building an Orthonormal Basis, Revisited" — which is numerically
/// stable for every `n`, including the `n.z ≈ −1` pole that the naive
/// cross-with-a-fixed-axis method handles badly.
#[inline]
pub fn ortho_basis(n: Vec3) -> (Vec3, Vec3) {
    let sign = if n.z >= 0.0 { 1.0 } else { -1.0 };
    let a = -1.0 / (sign + n.z);
    let b = n.x * n.y * a;
    let tangent = Vec3 {
        x: 1.0 + sign * n.x * n.x * a,
        y: sign * b,
        z: -sign * n.x,
    };
    let bitangent = Vec3 {
        x: b,
        y: sign + n.y * n.y * a,
        z: -n.y,
    };
    (tangent, bitangent)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_and_cross_are_textbook() {
        let x = vec3(1.0, 0.0, 0.0);
        let y = vec3(0.0, 1.0, 0.0);
        // x · y = 0 (orthogonal), x × y = z (right-handed).
        assert_eq!(x.dot(y), 0.0);
        let z = x.cross(y);
        assert!((z.sub(vec3(0.0, 0.0, 1.0))).length() < 1e-7);
    }

    #[test]
    fn normalized_rejects_the_zero_vector() {
        assert!(Vec3::ZERO.normalized().is_none());
        let n = vec3(3.0, 0.0, 4.0).normalized().unwrap();
        assert!((n.length() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn reflect_flips_the_normal_component() {
        // A ray heading straight down reflects straight back up.
        let r = vec3(0.0, 0.0, -1.0).reflect(vec3(0.0, 0.0, 1.0));
        assert!((r.sub(vec3(0.0, 0.0, 1.0))).length() < 1e-6);
    }

    #[test]
    fn ortho_basis_is_orthonormal_for_every_axis() {
        // The basis must be orthonormal for tricky normals too — the
        // n.z ≈ −1 pole is the case the naive method fails.
        for n in [
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, -1.0),
            vec3(1.0, 0.0, 0.0),
            vec3(0.577, 0.577, 0.577).normalized().unwrap(),
        ] {
            let (t, b) = ortho_basis(n);
            assert!((t.length() - 1.0).abs() < 1e-5, "tangent not unit");
            assert!((b.length() - 1.0).abs() < 1e-5, "bitangent not unit");
            assert!(t.dot(n).abs() < 1e-5, "tangent not ⟂ n");
            assert!(b.dot(n).abs() < 1e-5, "bitangent not ⟂ n");
            assert!(t.dot(b).abs() < 1e-5, "tangent not ⟂ bitangent");
        }
    }

    #[test]
    fn min_max_reduce_aabb_corners() {
        let a = vec3(1.0, 5.0, 3.0);
        let b = vec3(4.0, 2.0, 6.0);
        assert_eq!(a.min(b), vec3(1.0, 2.0, 3.0));
        assert_eq!(a.max(b), vec3(4.0, 5.0, 6.0));
    }
}

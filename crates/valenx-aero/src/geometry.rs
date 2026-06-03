//! Body geometry — a triangle mesh and the helpers that extract it
//! from Valenx's canonical mesh / CAD types.
//!
//! The wind-tunnel solver treats the immersed body as a soup of
//! oriented triangles ([`TriMesh`]). That is all the immersed-boundary
//! voxelizer ([`crate::immersed`]) needs: a closed-ish triangle shell
//! it can ray-cast against and integrate forces over. Splitting the
//! geometry out here keeps the solver agnostic of where the triangles
//! came from — a [`valenx_mesh::Mesh`], a [`valenx_cad::Solid`]
//! tessellation, or a hand-built list all flow in through the same
//! [`TriMesh`].

use nalgebra::Vector3;

use crate::error::AeroError;

/// One oriented triangle — three vertices, counter-clockwise when
/// viewed from outside the body so the geometric normal points
/// outward.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Triangle {
    /// Vertex 0.
    pub a: Vector3<f64>,
    /// Vertex 1.
    pub b: Vector3<f64>,
    /// Vertex 2.
    pub c: Vector3<f64>,
}

impl Triangle {
    /// Build a triangle from three corner points.
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Triangle {
        Triangle { a, b, c }
    }

    /// The (un-normalised) geometric normal — `(b−a) × (c−a)`. Its
    /// magnitude is twice the triangle area.
    #[inline]
    pub fn raw_normal(&self) -> Vector3<f64> {
        (self.b - self.a).cross(&(self.c - self.a))
    }

    /// The unit outward normal. Returns the zero vector for a
    /// degenerate (zero-area) triangle.
    pub fn normal(&self) -> Vector3<f64> {
        let n = self.raw_normal();
        let m = n.norm();
        if m > 1e-30 {
            n / m
        } else {
            Vector3::zeros()
        }
    }

    /// The triangle area.
    #[inline]
    pub fn area(&self) -> f64 {
        0.5 * self.raw_normal().norm()
    }

    /// The centroid (arithmetic mean of the three vertices).
    #[inline]
    pub fn centroid(&self) -> Vector3<f64> {
        (self.a + self.b + self.c) / 3.0
    }

    /// True if the triangle is degenerate — effectively zero area.
    pub fn is_degenerate(&self) -> bool {
        self.area() < 1e-14
    }
}

/// An axis-aligned bounding box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Vector3<f64>,
    /// Maximum corner.
    pub max: Vector3<f64>,
}

impl Aabb {
    /// The box extent `(dx, dy, dz)`.
    #[inline]
    pub fn extent(&self) -> Vector3<f64> {
        self.max - self.min
    }

    /// The box centre.
    #[inline]
    pub fn centre(&self) -> Vector3<f64> {
        0.5 * (self.min + self.max)
    }

    /// The longest axis extent.
    #[inline]
    pub fn longest(&self) -> f64 {
        let e = self.extent();
        e.x.max(e.y).max(e.z)
    }
}

/// A triangle-soup mesh — the body dropped into the virtual wind
/// tunnel.
#[derive(Clone, Debug, Default)]
pub struct TriMesh {
    /// The body's triangles, outward-oriented.
    pub triangles: Vec<Triangle>,
}

impl TriMesh {
    /// An empty mesh.
    pub fn new() -> TriMesh {
        TriMesh {
            triangles: Vec::new(),
        }
    }

    /// Build a mesh from a triangle list.
    pub fn from_triangles(triangles: Vec<Triangle>) -> TriMesh {
        TriMesh { triangles }
    }

    /// Extract a [`TriMesh`] from a [`valenx_mesh::Mesh`].
    ///
    /// Reads every `Tri3` element block; quads / higher-order elements
    /// in the mesh are ignored (an external-aero body is a surface
    /// shell, which a tessellator delivers as triangles). Returns an
    /// [`AeroError::BadGeometry`] if no triangles survive.
    pub fn from_mesh(mesh: &valenx_mesh::Mesh) -> Result<TriMesh, AeroError> {
        use valenx_mesh::element::ElementType;
        let mut tris = Vec::new();
        for block in &mesh.element_blocks {
            if block.element_type != ElementType::Tri3 {
                continue;
            }
            for tri in block.connectivity.chunks_exact(3) {
                let (ia, ib, ic) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                if ia >= mesh.nodes.len() || ib >= mesh.nodes.len() || ic >= mesh.nodes.len() {
                    return Err(AeroError::BadGeometry(format!(
                        "triangle references node out of range (n_nodes = {})",
                        mesh.nodes.len()
                    )));
                }
                tris.push(Triangle::new(mesh.nodes[ia], mesh.nodes[ib], mesh.nodes[ic]));
            }
        }
        tris.retain(|t| !t.is_degenerate());
        if tris.is_empty() {
            return Err(AeroError::BadGeometry(
                "mesh contains no non-degenerate Tri3 elements".into(),
            ));
        }
        Ok(TriMesh::from_triangles(tris))
    }

    /// Extract a [`TriMesh`] from a [`valenx_cad::Solid`] by
    /// tessellating its BRep at the given chord tolerance.
    pub fn from_solid(solid: &valenx_cad::Solid, tolerance: f64) -> Result<TriMesh, AeroError> {
        let mesh = valenx_cad::solid_to_mesh(solid, tolerance)
            .map_err(|e| AeroError::BadGeometry(format!("tessellation failed: {e}")))?;
        TriMesh::from_mesh(&mesh)
    }

    /// Triangle count.
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    /// True if the mesh holds no triangles.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// The total surface area of the body.
    pub fn surface_area(&self) -> f64 {
        self.triangles.iter().map(|t| t.area()).sum()
    }

    /// The axis-aligned bounding box of all vertices. Returns `None`
    /// for an empty mesh.
    pub fn aabb(&self) -> Option<Aabb> {
        if self.triangles.is_empty() {
            return None;
        }
        let mut min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for t in &self.triangles {
            for v in [t.a, t.b, t.c] {
                min = min.inf(&v);
                max = max.sup(&v);
            }
        }
        Some(Aabb { min, max })
    }

    /// Translate every vertex by `delta` (used to place the body in
    /// the tunnel domain).
    pub fn translate(&mut self, delta: Vector3<f64>) {
        for t in &mut self.triangles {
            t.a += delta;
            t.b += delta;
            t.c += delta;
        }
    }

    /// Rotate every vertex about the body's bounding-box centre by the
    /// given yaw / pitch / roll (radians), applied in that order
    /// (yaw about z, pitch about y, roll about x). Used to orient the
    /// body at a yaw / angle-of-attack inside the tunnel.
    pub fn rotate_about_centre(&mut self, yaw: f64, pitch: f64, roll: f64) {
        let centre = match self.aabb() {
            Some(b) => b.centre(),
            None => return,
        };
        let rot = yaw_pitch_roll_matrix(yaw, pitch, roll);
        for t in &mut self.triangles {
            t.a = centre + rot * (t.a - centre);
            t.b = centre + rot * (t.b - centre);
            t.c = centre + rot * (t.c - centre);
        }
    }

    /// The projected frontal area of the body onto the plane normal to
    /// `dir` — the reference area for drag-coefficient
    /// normalisation.
    ///
    /// Computed by summing the *positive* projected area of every
    /// triangle facet onto the plane (each facet that faces into the
    /// wind contributes `area · |n·dir|`); for a closed convex body
    /// the front and back facets each contribute the same silhouette,
    /// so the windward sum equals the true projected frontal area.
    /// This avoids needing a rasteriser and is exact for convex
    /// bodies; for a concave body it is a slight over-estimate (a
    /// documented v1 simplification).
    pub fn frontal_area(&self, dir: Vector3<f64>) -> f64 {
        let d = match dir.try_normalize(1e-12) {
            Some(d) => d,
            None => return 0.0,
        };
        let mut windward = 0.0;
        for t in &self.triangles {
            let proj = t.raw_normal().dot(&d); // = 2·area·(n·d)
            if proj < 0.0 {
                // Facet faces into the wind (normal opposes flow dir).
                windward += -0.5 * proj;
            }
        }
        windward
    }
}

/// The yaw-pitch-roll rotation matrix `Rz(yaw)·Ry(pitch)·Rx(roll)`.
pub fn yaw_pitch_roll_matrix(yaw: f64, pitch: f64, roll: f64) -> nalgebra::Matrix3<f64> {
    use nalgebra::Matrix3;
    let (cy, sy) = (yaw.cos(), yaw.sin());
    let (cp, sp) = (pitch.cos(), pitch.sin());
    let (cr, sr) = (roll.cos(), roll.sin());
    let rz = Matrix3::new(cy, -sy, 0.0, sy, cy, 0.0, 0.0, 0.0, 1.0);
    let ry = Matrix3::new(cp, 0.0, sp, 0.0, 1.0, 0.0, -sp, 0.0, cp);
    let rx = Matrix3::new(1.0, 0.0, 0.0, 0.0, cr, -sr, 0.0, sr, cr);
    rz * ry * rx
}

/// Build a simple axis-aligned box body — `nx·ny·nz` need not be a
/// fine mesh; a box is the canonical bluff-body test geometry. The 12
/// triangles are outward-oriented.
pub fn box_body(min: Vector3<f64>, max: Vector3<f64>) -> TriMesh {
    let v = |x: f64, y: f64, z: f64| Vector3::new(x, y, z);
    let (x0, y0, z0) = (min.x, min.y, min.z);
    let (x1, y1, z1) = (max.x, max.y, max.z);
    // Eight corners.
    let c000 = v(x0, y0, z0);
    let c100 = v(x1, y0, z0);
    let c110 = v(x1, y1, z0);
    let c010 = v(x0, y1, z0);
    let c001 = v(x0, y0, z1);
    let c101 = v(x1, y0, z1);
    let c111 = v(x1, y1, z1);
    let c011 = v(x0, y1, z1);
    let mut tris = Vec::with_capacity(12);
    let mut quad = |a, b, c, d| {
        // Two CCW triangles for the quad (a,b,c,d).
        tris.push(Triangle::new(a, b, c));
        tris.push(Triangle::new(a, c, d));
    };
    // -z bottom (normal -z): order so cross points down.
    quad(c000, c010, c110, c100);
    // +z top (normal +z).
    quad(c001, c101, c111, c011);
    // -y front (normal -y).
    quad(c000, c100, c101, c001);
    // +y back (normal +y).
    quad(c010, c011, c111, c110);
    // -x left (normal -x).
    quad(c000, c001, c011, c010);
    // +x right (normal +x).
    quad(c100, c110, c111, c101);
    TriMesh::from_triangles(tris)
}

/// Build a UV-sphere body of `lat × lon` facets — the canonical
/// curved-body drag test geometry. Outward-oriented.
pub fn sphere_body(centre: Vector3<f64>, radius: f64, lat: usize, lon: usize) -> TriMesh {
    let lat = lat.max(2);
    let lon = lon.max(3);
    let mut tris = Vec::new();
    let pi = std::f64::consts::PI;
    let point = |theta: f64, phi: f64| -> Vector3<f64> {
        centre
            + radius
                * Vector3::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                )
    };
    for i in 0..lat {
        let t0 = pi * i as f64 / lat as f64;
        let t1 = pi * (i + 1) as f64 / lat as f64;
        // The top (i = 0) and bottom (i = lat-1) latitude bands touch a
        // pole: the two "top" vertices (or two "bottom" vertices) of the
        // quad collapse onto the single pole point, so one of the two
        // split triangles is degenerate (zero area, zero normal) and
        // must NOT be emitted — a pole band is a triangle fan, not a
        // quad strip.
        let at_north_pole = i == 0;
        let at_south_pole = i == lat - 1;
        for j in 0..lon {
            let p0 = 2.0 * pi * j as f64 / lon as f64;
            let p1 = 2.0 * pi * (j + 1) as f64 / lon as f64;
            let a = point(t0, p0);
            let b = point(t1, p0);
            let c = point(t1, p1);
            let d = point(t0, p1);
            // Outward-oriented triangles, skipping the degenerate one at
            // each pole.
            if !at_south_pole {
                tris.push(Triangle::new(a, b, c));
            }
            if !at_north_pole {
                tris.push(Triangle::new(a, c, d));
            }
        }
    }
    TriMesh::from_triangles(tris)
}

/// The half-thickness `y_t(x)` of a NACA 4-digit *symmetric* airfoil at
/// chord fraction `x ∈ [0, 1]`, for a maximum-thickness fraction `t`.
///
/// The standard NACA 4-digit thickness distribution:
///
/// ```text
///   y_t = 5·t·(0.2969·√x − 0.1260·x − 0.3516·x²
///              + 0.2843·x³ − 0.1015·x⁴)
/// ```
///
/// NACA 0012 (the canonical validation airfoil) is `t = 0.12`. The
/// distribution is for the *open* trailing edge of the original NACA
/// definition; the small open gap is closed by the cap facets when the
/// section is extruded into a watertight wing.
pub fn naca4_half_thickness(x: f64, t: f64) -> f64 {
    let x = x.clamp(0.0, 1.0);
    5.0 * t
        * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x
            + 0.2843 * x * x * x
            - 0.1015 * x * x * x * x)
}

/// Build a 3-D wing body — a NACA 4-digit *symmetric* airfoil section
/// extruded spanwise into a closed triangle shell.
///
/// `chord` is the chord length (along `+x`), `span` the spanwise extent
/// (along `+y`), `thickness_fraction` the airfoil's maximum thickness
/// as a fraction of chord (`0.12` for a NACA 0012), and `chord_panels`
/// the number of chordwise panels per surface. The section lies in the
/// `x`–`z` plane (chord along `x`, thickness along `z`) and is extruded
/// along `y`; the leading edge sits at the origin. The result is a
/// watertight outward-oriented mesh — two cambered surfaces, two
/// end-caps and a closed trailing edge — ready for the wind tunnel.
///
/// A symmetric section at zero geometric incidence; the angle of attack
/// is applied by yawing / pitching the *wind* (or rotating the body),
/// exactly as a real wind-tunnel sweep does.
pub fn naca_wing(
    chord: f64,
    span: f64,
    thickness_fraction: f64,
    chord_panels: usize,
) -> TriMesh {
    let n = chord_panels.max(4);
    let t = thickness_fraction.max(1e-4);
    // Chordwise stations with cosine clustering — the leading edge,
    // where the curvature is highest, gets the finest spacing.
    let stations: Vec<f64> = (0..=n)
        .map(|i| {
            let beta = std::f64::consts::PI * i as f64 / n as f64;
            0.5 * (1.0 - beta.cos())
        })
        .collect();
    // The two spanwise stations of the section extrusion.
    let y0 = 0.0;
    let y1 = span.max(1e-6);

    // Section profile points: upper surface (z = +y_t) and lower
    // surface (z = −y_t) at each chord station, scaled by the chord.
    let upper = |x_frac: f64| -> (f64, f64) {
        (x_frac * chord, naca4_half_thickness(x_frac, t) * chord)
    };

    let mut tris = Vec::new();
    let v = Vector3::new;

    // The two cambered surfaces, panel by panel along the chord.
    for w in stations.windows(2) {
        let (xa, xb) = (w[0], w[1]);
        let (xau, zau) = upper(xa);
        let (xbu, zbu) = upper(xb);
        // Upper surface quad (outward normal +z-ish): the spanwise
        // strip between chord stations xa and xb.
        let u_a0 = v(xau, y0, zau);
        let u_a1 = v(xau, y1, zau);
        let u_b0 = v(xbu, y0, zbu);
        let u_b1 = v(xbu, y1, zbu);
        tris.push(Triangle::new(u_a0, u_b0, u_b1));
        tris.push(Triangle::new(u_a0, u_b1, u_a1));
        // Lower surface quad (outward normal −z-ish): z mirrored.
        let l_a0 = v(xau, y0, -zau);
        let l_a1 = v(xau, y1, -zau);
        let l_b0 = v(xbu, y0, -zbu);
        let l_b1 = v(xbu, y1, -zbu);
        tris.push(Triangle::new(l_a0, l_b1, l_b0));
        tris.push(Triangle::new(l_a0, l_a1, l_b1));
    }

    // The two end-caps — fan the section polygon from the leading edge.
    for w in stations.windows(2) {
        let (xa, xb) = (w[0], w[1]);
        let (xau, zau) = upper(xa);
        let (xbu, zbu) = upper(xb);
        // y0 cap (outward normal −y): upper and lower triangles.
        tris.push(Triangle::new(
            v(xau, y0, zau),
            v(xau, y0, -zau),
            v(xbu, y0, zbu),
        ));
        tris.push(Triangle::new(
            v(xbu, y0, zbu),
            v(xau, y0, -zau),
            v(xbu, y0, -zbu),
        ));
        // y1 cap (outward normal +y): reversed winding.
        tris.push(Triangle::new(
            v(xau, y1, zau),
            v(xbu, y1, zbu),
            v(xau, y1, -zau),
        ));
        tris.push(Triangle::new(
            v(xbu, y1, zbu),
            v(xbu, y1, -zbu),
            v(xau, y1, -zau),
        ));
    }

    // The trailing edge — close the small open gap with a flat face
    // (the NACA definition leaves the TE open; a watertight body needs
    // it sealed). The TE half-thickness at x = 1.
    let (xte, zte) = upper(1.0);
    if zte > 1e-9 {
        // Outward normal +x.
        tris.push(Triangle::new(
            v(xte, y0, zte),
            v(xte, y1, zte),
            v(xte, y0, -zte),
        ));
        tris.push(Triangle::new(
            v(xte, y1, zte),
            v(xte, y1, -zte),
            v(xte, y0, -zte),
        ));
    }

    TriMesh::from_triangles(tris)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_normal_area_centroid() {
        let t = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        );
        // Normal points +z, area = 2.
        assert!((t.normal() - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-12);
        assert!((t.area() - 2.0).abs() < 1e-12);
        let c = t.centroid();
        assert!((c - Vector3::new(2.0 / 3.0, 2.0 / 3.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn box_body_has_twelve_triangles_and_correct_area() {
        let b = box_body(Vector3::zeros(), Vector3::new(2.0, 3.0, 4.0));
        assert_eq!(b.len(), 12);
        // Surface area of a 2×3×4 box = 2(6 + 8 + 12) = 52.
        assert!((b.surface_area() - 52.0).abs() < 1e-9);
    }

    #[test]
    fn box_body_normals_all_point_outward() {
        // Each facet normal must point away from the box centre.
        let b = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        for t in &b.triangles {
            let outward = t.centroid(); // centre at origin → centroid = outward dir
            assert!(
                t.normal().dot(&outward) > 0.0,
                "facet normal should point outward"
            );
        }
    }

    #[test]
    fn sphere_body_area_approaches_analytic() {
        // A fine UV-sphere's triangulated area approaches 4πr².
        let r = 1.5;
        let s = sphere_body(Vector3::zeros(), r, 48, 96);
        let analytic = 4.0 * std::f64::consts::PI * r * r;
        let rel = (s.surface_area() - analytic).abs() / analytic;
        assert!(rel < 0.02, "sphere area off by {rel}");
    }

    #[test]
    fn sphere_normals_point_radially_outward() {
        let s = sphere_body(Vector3::zeros(), 1.0, 12, 24);
        for t in &s.triangles {
            let radial = t.centroid().normalize();
            assert!(
                t.normal().dot(&radial) > 0.5,
                "sphere facet normal should be roughly radial-outward"
            );
        }
    }

    #[test]
    fn aabb_bounds_the_body() {
        let b = box_body(Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 6.0, 9.0));
        let bb = b.aabb().unwrap();
        assert!((bb.min - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-12);
        assert!((bb.max - Vector3::new(4.0, 6.0, 9.0)).norm() < 1e-12);
        assert!((bb.extent() - Vector3::new(3.0, 4.0, 6.0)).norm() < 1e-12);
    }

    #[test]
    fn frontal_area_of_a_box_is_the_silhouette() {
        // A 2×4×6 box, wind along +x: the frontal silhouette is the
        // y-z face = 4·6 = 24.
        let b = box_body(Vector3::zeros(), Vector3::new(2.0, 4.0, 6.0));
        let fa = b.frontal_area(Vector3::new(1.0, 0.0, 0.0));
        assert!((fa - 24.0).abs() < 1e-9, "frontal area {fa} should be 24");
        // Wind along +z: silhouette is the x-y face = 2·4 = 8.
        let fz = b.frontal_area(Vector3::new(0.0, 0.0, 1.0));
        assert!((fz - 8.0).abs() < 1e-9, "frontal area {fz} should be 8");
    }

    #[test]
    fn translate_shifts_the_bounding_box() {
        let mut b = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        b.translate(Vector3::new(5.0, 0.0, -2.0));
        let bb = b.aabb().unwrap();
        assert!((bb.min - Vector3::new(5.0, 0.0, -2.0)).norm() < 1e-12);
    }

    #[test]
    fn from_mesh_reads_tri3_blocks() {
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut m = valenx_mesh::Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let tm = TriMesh::from_mesh(&m).unwrap();
        assert_eq!(tm.len(), 1);
        assert!((tm.surface_area() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn from_mesh_rejects_an_empty_mesh() {
        let m = valenx_mesh::Mesh::new("empty");
        assert!(TriMesh::from_mesh(&m).is_err());
    }

    #[test]
    fn naca4_thickness_matches_the_standard_distribution() {
        // NACA 0012: max half-thickness ≈ 0.06·chord, occurring near
        // x/c ≈ 0.30. The thickness is zero at the leading edge.
        assert!(naca4_half_thickness(0.0, 0.12).abs() < 1e-9);
        let y_max = naca4_half_thickness(0.30, 0.12);
        // The NACA 0012 half-thickness peaks at ≈ 0.0600·c.
        assert!(
            (y_max - 0.060).abs() < 0.003,
            "NACA 0012 peak half-thickness {y_max} should be ≈ 0.060"
        );
        // Monotone fall from the peak to the trailing edge.
        assert!(naca4_half_thickness(0.6, 0.12) < y_max);
        assert!(naca4_half_thickness(0.95, 0.12) < naca4_half_thickness(0.6, 0.12));
    }

    #[test]
    fn naca_wing_is_a_watertight_outward_shell() {
        // A NACA 0012 wing must be a closed mesh with outward normals
        // and a bounding box matching the chord / span / thickness.
        let wing = naca_wing(1.0, 3.0, 0.12, 24);
        assert!(!wing.is_empty(), "the wing should have triangles");
        let bb = wing.aabb().unwrap();
        // Chord runs 0→1 along x, span 0→3 along y.
        assert!((bb.min.x - 0.0).abs() < 1e-6 && (bb.max.x - 1.0).abs() < 1e-6);
        assert!((bb.min.y - 0.0).abs() < 1e-6 && (bb.max.y - 3.0).abs() < 1e-6);
        // The section is symmetric: |z| extent ≈ 2·0.06 = 0.12.
        assert!(
            (bb.extent().z - 0.12).abs() < 0.02,
            "wing thickness extent {} should be ≈ 0.12",
            bb.extent().z
        );
        // The voxelizer must classify the wing interior as solid — a
        // point at mid-chord, mid-span, on the chord line is inside.
        assert!(
            crate::immersed::point_inside(Vector3::new(0.4, 1.5, 0.0), &wing),
            "the wing mid-section should be inside the body"
        );
        // A point well above the wing is outside.
        assert!(!crate::immersed::point_inside(
            Vector3::new(0.4, 1.5, 1.0),
            &wing
        ));
    }

    #[test]
    fn naca_wing_frontal_area_is_small_for_a_thin_section() {
        // Edge-on (wind along +x) a thin airfoil presents only its
        // thin trailing-edge-to-leading-edge silhouette — the frontal
        // area is a small fraction of the chord×span planform.
        let wing = naca_wing(1.0, 2.0, 0.12, 32);
        let frontal = wing.frontal_area(Vector3::new(1.0, 0.0, 0.0));
        let planform = 1.0 * 2.0;
        assert!(
            frontal > 0.0 && frontal < 0.25 * planform,
            "wing frontal area {frontal} should be a small fraction of \
             the planform {planform}"
        );
    }
}

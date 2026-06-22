//! Shared triangle-mesh **primitive builder** for procedural parts.
//!
//! Every product that draws a machine (bearings, gear trains, pulleys,
//! fasteners, …) used to hand-roll its own triangle soup with a private
//! `push_cylinder` copy, there was **no sphere primitive at all**, and the
//! result was a single flat colour per part. [`MeshBuilder`] fixes that root
//! geometry-quality problem: it is one lightweight, well-tested triangle
//! accumulator that emits a [`valenx_mesh::Mesh`] of [`ElementType::Tri3`]
//! triangles **with a parallel per-vertex colour buffer** laid out exactly the
//! way the renderer's coloured path consumes it.
//!
//! It is deliberately *not* a B-rep / CSG kernel — there are no booleans and no
//! true holes. For solids that need real subtraction (a bolt hole drilled
//! through a flange) use `valenx-cad`. This builder is the fast path for
//! "lath­e/extrude/stamp a pile of coloured parts and hand the renderer a mesh".
//!
//! # Colour alignment with the renderer
//!
//! The wgpu renderer's plain path
//! ([`crate::wgpu_renderer::triangles_to_vertices`]) walks the surface
//! triangles and emits **three [`Vertex`](crate::wgpu_renderer::Vertex) entries
//! per triangle** (triangle-major, then the three corners of each triangle).
//! The coloured path
//! ([`crate::wgpu_renderer::triangles_to_vertices_colored`], ~`wgpu_renderer.rs:943`)
//! reads a parallel `&[[f32; 3]]` indexed by that *same* emitted-vertex stream —
//! i.e. it expects `3 × triangle_count` colours, one per corner, in the same
//! order. [`WorkspaceProduct::vertex_colors`](crate::WorkspaceProduct) (lib.rs:341)
//! carries exactly this buffer.
//!
//! [`MeshBuilder`] keeps its `colors` vector in lock-step with the triangle
//! list: every appended [`Tri3`](ElementType::Tri3) pushes its node triple to
//! the mesh **and** three identical colour entries (the primitive's flat
//! colour) to `colors`, in the same triangle-major order. So
//! [`MeshBuilder::into_mesh_and_colors`] returns a `(Mesh, Vec<[f32; 3]>)` pair
//! that drops straight into the coloured path with
//! `colors.len() == 3 * triangle_count`, no re-indexing required.
//!
//! # Returned node ranges
//!
//! Each primitive returns the [`Range<usize>`] of **node indices** it appended
//! to the underlying mesh (`start..end`, half-open, into `mesh.nodes`).
//! Producers that animate a moving part (a spinning bearing ball, a sliding
//! rod) keep that range and later transform `mesh.nodes[range]` rigidly. Ranges
//! from successive primitives are contiguous: the `end` of part *k* equals the
//! `start` of part *k + 1*.

use std::ops::Range;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

const TAU: f64 = std::f64::consts::TAU;

/// Accumulates [`Tri3`](ElementType::Tri3) triangles plus a parallel per-vertex
/// colour buffer, then bakes them into a [`valenx_mesh::Mesh`].
///
/// Build it with [`MeshBuilder::new`], stamp primitives (each returns the
/// [`Range<usize>`] of node indices it added), then call
/// [`MeshBuilder::into_mesh_and_colors`] to get the `(Mesh, colors)` pair the
/// renderer's coloured path wants. See the [module docs](self) for the colour
/// layout and the contiguous-range guarantee.
#[derive(Clone, Debug, Default)]
pub struct MeshBuilder {
    /// All vertices, in append order. Primitive node ranges index into this.
    nodes: Vec<Vector3<f64>>,
    /// Triangle list — each entry is a node-index triple, wound
    /// counter-clockwise as seen from outside so face normals point out.
    tris: Vec<[u32; 3]>,
    /// Per-triangle-vertex colours, triangle-major (three identical entries per
    /// triangle). Always kept at `3 * tris.len()`.
    colors: Vec<[f32; 3]>,
}

impl MeshBuilder {
    /// A fresh, empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Vertices accumulated so far. The next primitive's node range starts here.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Triangles accumulated so far. The colour buffer length is always
    /// `3 × triangle_count`.
    pub fn triangle_count(&self) -> usize {
        self.tris.len()
    }

    /// Consume the builder, producing a [`valenx_mesh::Mesh`] (one
    /// [`Tri3`](ElementType::Tri3) [`ElementBlock`]) and the parallel per-vertex
    /// colour buffer aligned to the renderer's coloured path
    /// ([`crate::wgpu_renderer::triangles_to_vertices_colored`]).
    ///
    /// The returned `colors.len()` is always `3 × mesh` triangle count: one
    /// `[r, g, b]` per emitted surface vertex, triangle-major. Stats are
    /// refreshed so node / element counts are populated.
    pub fn into_mesh_and_colors(self) -> (Mesh, Vec<[f32; 3]>) {
        let mut mesh = Mesh::new("mesh-prims");
        mesh.nodes = self.nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.reserve(self.tris.len() * 3);
        for t in &self.tris {
            block.connectivity.extend_from_slice(t);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        (mesh, self.colors)
    }

    // ---- internal helpers -------------------------------------------------

    /// Push one vertex, return its index.
    fn vert(&mut self, p: Vector3<f64>) -> u32 {
        let i = self.nodes.len() as u32;
        self.nodes.push(p);
        i
    }

    /// Append one triangle `(a, b, c)` (node indices, outward winding) and its
    /// three flat-colour entries, keeping `colors` triangle-major.
    fn tri(&mut self, a: u32, b: u32, c: u32, color: [f32; 3]) {
        self.tris.push([a, b, c]);
        self.colors.push(color);
        self.colors.push(color);
        self.colors.push(color);
    }

    /// Append a quad as two triangles `(a, b, c)` + `(a, c, d)`, wound so the
    /// outward face is the one with corners listed counter-clockwise.
    fn quad(&mut self, a: u32, b: u32, c: u32, d: u32, color: [f32; 3]) {
        self.tri(a, b, c, color);
        self.tri(a, c, d, color);
    }
}

/// Build a right-handed orthonormal frame `(u, v, w)` with `w` along the
/// (normalized) `axis`. `u`/`v` span the plane perpendicular to the axis, so a
/// ring at angle `a` is `center + r·(u·cos a + v·sin a)`. The seed is chosen so
/// it is never parallel to `w`, giving a numerically stable cross product for
/// any axis direction (x / y / z / oblique). A zero or tiny axis defaults to
/// `+z`.
fn axis_frame(axis: [f64; 3]) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let mut w = Vector3::new(axis[0], axis[1], axis[2]);
    let n = w.norm();
    w = if n < 1e-12 {
        Vector3::new(0.0, 0.0, 1.0)
    } else {
        w / n
    };
    let seed = if w.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = w.cross(&seed).normalize();
    let v = w.cross(&u);
    (u, v, w)
}

impl MeshBuilder {
    /// Append a **capped cylinder** of `radius` and `length`, centred at
    /// `center` with its long axis along `axis` (so it can point along x, y, z
    /// or any oblique direction). Built from two `segments`-sided end rings, a
    /// quad band between them, and a triangle-fan cap at each end (so it reads
    /// as a solid), all outward-facing.
    ///
    /// Returns the [`Range<usize>`] of node indices appended:
    /// `2·segments` ring vertices + `2` cap centres.
    #[allow(clippy::too_many_arguments)]
    pub fn cylinder(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        radius: f64,
        length: f64,
        segments: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        // A cylinder is a two-station revolve of the rectangle profile
        // (r, -h/2) → (r, +h/2). Express it as a frustum with equal radii so we
        // get the caps for free and share one code path.
        self.cone(center, axis, radius, radius, length, segments, color)
    }

    /// Append a **hollow annular cylinder** (a tube / pipe / bearing race /
    /// washer) of inner radius `r_inner`, outer radius `r_outer` and `length`,
    /// centred at `center` along `axis`. Has an outer wall, an inner wall
    /// (inward-facing so its normals point into the bore), and a flat annular
    /// ring closing each end. `r_inner` is clamped below `r_outer`.
    ///
    /// Returns the appended node [`Range<usize>`]: `4·segments` ring vertices
    /// (inner+outer, each end).
    #[allow(clippy::too_many_arguments)]
    pub fn tube(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        r_inner: f64,
        r_outer: f64,
        length: f64,
        segments: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        let seg = segments.max(3);
        let ri = r_inner.max(0.0).min(r_outer);
        let (u, v, w) = axis_frame(axis);
        let c = Vector3::new(center[0], center[1], center[2]);
        let h = length * 0.5;
        // Four rings: outer-bottom, outer-top, inner-bottom, inner-top.
        let ring = |b: &mut Self, r: f64, z: f64| -> u32 {
            let base = b.nodes.len() as u32;
            for k in 0..seg {
                let a = k as f64 / seg as f64 * TAU;
                let off = u * (r * a.cos()) + v * (r * a.sin()) + w * z;
                b.vert(c + off);
            }
            base
        };
        let ob = ring(self, r_outer, -h);
        let ot = ring(self, r_outer, h);
        let ib = ring(self, ri, -h);
        let it = ring(self, ri, h);
        let seg_u = seg as u32;
        for k in 0..seg_u {
            let k1 = (k + 1) % seg_u;
            // Outer wall — outward (normals point away from the axis).
            self.quad(ob + k, ob + k1, ot + k1, ot + k, color);
            // Inner wall — reversed so normals point into the bore.
            self.quad(ib + k1, ib + k, it + k, it + k1, color);
            // Bottom annulus (faces -w): outer→inner.
            self.quad(ob + k, ib + k, ib + k1, ob + k1, color);
            // Top annulus (faces +w): inner→outer.
            self.quad(ot + k, ot + k1, it + k1, it + k, color);
        }
        start..self.nodes.len()
    }

    /// Append a **UV sphere** of `radius` centred at `center`, with `lat_segs`
    /// latitude bands (pole to pole) and `lon_segs` longitude divisions. The
    /// two poles are single vertices joined by triangle fans; the interior
    /// bands are quads. Outward-facing. This is the primitive the old code was
    /// missing entirely (bearing balls, ball-joint centres, fillet beads).
    ///
    /// Vertex layout (the returned [`Range<usize>`]): the north pole, then
    /// `(lat_segs − 1)` interior rings of `lon_segs` vertices each, then the
    /// south pole — i.e. `2 + (lat_segs − 1)·lon_segs` vertices.
    pub fn sphere(
        &mut self,
        center: [f64; 3],
        radius: f64,
        lat_segs: usize,
        lon_segs: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        let lat = lat_segs.max(2);
        let lon = lon_segs.max(3);
        let c = Vector3::new(center[0], center[1], center[2]);

        let north = self.vert(c + Vector3::new(0.0, 0.0, radius));
        // Interior rings i = 1 ..= lat-1 at polar angle theta = i·pi/lat.
        let mut ring_base: Vec<u32> = Vec::with_capacity(lat - 1);
        for i in 1..lat {
            let theta = i as f64 / lat as f64 * std::f64::consts::PI;
            let (st, ct) = theta.sin_cos();
            let base = self.nodes.len() as u32;
            ring_base.push(base);
            for j in 0..lon {
                let phi = j as f64 / lon as f64 * TAU;
                let (sp, cp) = phi.sin_cos();
                self.vert(
                    c + Vector3::new(radius * st * cp, radius * st * sp, radius * ct),
                );
            }
        }
        let south = self.vert(c + Vector3::new(0.0, 0.0, -radius));

        // North cap fan.
        let first = ring_base[0];
        for j in 0..lon {
            let j1 = (j + 1) % lon;
            self.tri(north, first + j as u32, first + j1 as u32, color);
        }
        // Interior quad bands.
        for r in 0..ring_base.len() - 1 {
            let a = ring_base[r];
            let b = ring_base[r + 1];
            for j in 0..lon {
                let j1 = (j + 1) % lon;
                self.quad(
                    a + j as u32,
                    b + j as u32,
                    b + j1 as u32,
                    a + j1 as u32,
                    color,
                );
            }
        }
        // South cap fan.
        let last = *ring_base.last().unwrap();
        for j in 0..lon {
            let j1 = (j + 1) % lon;
            self.tri(south, last + j1 as u32, last + j as u32, color);
        }
        start..self.nodes.len()
    }

    /// Append an axis-aligned **box** of full extents `dims` (`[dx, dy, dz]`)
    /// centred at `center`, as 8 corner vertices and 12 outward triangles.
    ///
    /// Returns the appended node [`Range<usize>`] (8 vertices).
    pub fn cuboid(&mut self, center: [f64; 3], dims: [f64; 3], color: [f32; 3]) -> Range<usize> {
        let start = self.nodes.len();
        let [cx, cy, cz] = center;
        let (hx, hy, hz) = (dims[0] * 0.5, dims[1] * 0.5, dims[2] * 0.5);
        // 8 corners, indexed by (x,y,z) sign bits.
        let mut v = [0u32; 8];
        for (i, slot) in v.iter_mut().enumerate() {
            let sx = if i & 1 == 0 { -hx } else { hx };
            let sy = if i & 2 == 0 { -hy } else { hy };
            let sz = if i & 4 == 0 { -hz } else { hz };
            *slot = self.vert(Vector3::new(cx + sx, cy + sy, cz + sz));
        }
        // Each face wound CCW as seen from outside.
        self.quad(v[0], v[2], v[3], v[1], color); // -z
        self.quad(v[4], v[5], v[7], v[6], color); // +z
        self.quad(v[0], v[1], v[5], v[4], color); // -y
        self.quad(v[2], v[6], v[7], v[3], color); // +y
        self.quad(v[0], v[4], v[6], v[2], color); // -x
        self.quad(v[1], v[3], v[7], v[5], color); // +x
        start..self.nodes.len()
    }

    /// Append a **frustum** (truncated cone) of base radius `r_base`, top radius
    /// `r_top` and `length`, centred at `center` along `axis`. A side wall plus
    /// an end cap at each end with non-zero radius (a zero-radius end becomes a
    /// pointed apex, so `r_top = 0` gives a true cone and `r_base = r_top` gives
    /// a plain capped [`cylinder`](Self::cylinder)). Outward-facing.
    ///
    /// Returns the appended node [`Range<usize>`]: `2·segments` ring vertices
    /// plus one cap centre per non-zero-radius end.
    #[allow(clippy::too_many_arguments)]
    pub fn cone(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        r_base: f64,
        r_top: f64,
        length: f64,
        segments: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        let seg = segments.max(3);
        let (u, v, w) = axis_frame(axis);
        let c = Vector3::new(center[0], center[1], center[2]);
        let h = length * 0.5;
        let ring = |b: &mut Self, r: f64, z: f64| -> u32 {
            let base = b.nodes.len() as u32;
            for k in 0..seg {
                let a = k as f64 / seg as f64 * TAU;
                b.vert(c + u * (r * a.cos()) + v * (r * a.sin()) + w * z);
            }
            base
        };
        let rb = ring(self, r_base, -h);
        let rt = ring(self, r_top, h);
        let seg_u = seg as u32;
        for k in 0..seg_u {
            let k1 = (k + 1) % seg_u;
            self.quad(rb + k, rb + k1, rt + k1, rt + k, color);
        }
        // Bottom cap (faces -w): fan, wound so the normal points along -w.
        if r_base > 1e-12 {
            let cb = self.vert(c - w * h);
            for k in 0..seg_u {
                let k1 = (k + 1) % seg_u;
                self.tri(cb, rb + k1, rb + k, color);
            }
        }
        // Top cap (faces +w).
        if r_top > 1e-12 {
            let ct = self.vert(c + w * h);
            for k in 0..seg_u {
                let k1 = (k + 1) % seg_u;
                self.tri(ct, rt + k, rt + k1, color);
            }
        }
        start..self.nodes.len()
    }

    /// Append a **torus** (ring / O-ring / bearing-cage groove) with tube
    /// (minor) radius `minor_r` swept around a circle of `major_r`, centred at
    /// `center` in the plane perpendicular to `axis`. `seg_major` divisions
    /// around the big ring, `seg_minor` around the tube cross-section, all quads
    /// and outward-facing.
    ///
    /// Returns the appended node [`Range<usize>`]: `seg_major · seg_minor`
    /// vertices.
    #[allow(clippy::too_many_arguments)]
    pub fn torus(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        major_r: f64,
        minor_r: f64,
        seg_major: usize,
        seg_minor: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        let sm = seg_major.max(3);
        let sn = seg_minor.max(3);
        let (u, v, w) = axis_frame(axis);
        let c = Vector3::new(center[0], center[1], center[2]);
        // Ring i sits at angle alpha around the axis; its centre is on the major
        // circle, and the tube cross-section is spanned by the radial direction
        // `dir` and the axis `w`.
        let base = self.nodes.len() as u32;
        for i in 0..sm {
            let alpha = i as f64 / sm as f64 * TAU;
            let (sa, ca) = alpha.sin_cos();
            let dir = u * ca + v * sa;
            let ring_c = c + dir * major_r;
            for j in 0..sn {
                let beta = j as f64 / sn as f64 * TAU;
                let (sb, cb) = beta.sin_cos();
                self.vert(ring_c + dir * (minor_r * cb) + w * (minor_r * sb));
            }
        }
        let idx = |i: usize, j: usize| base + (i % sm * sn + j % sn) as u32;
        for i in 0..sm {
            for j in 0..sn {
                self.quad(
                    idx(i, j),
                    idx(i + 1, j),
                    idx(i + 1, j + 1),
                    idx(i, j + 1),
                    color,
                );
            }
        }
        start..self.nodes.len()
    }

    /// **Lathe** a 2-D `profile` of `(r, z)` pairs around `axis` (through
    /// `center`) by `angle_deg` in `segments` angular steps — a pulley
    /// V-groove, a fastener head, a simple hull of revolution. `r` is the
    /// distance from the axis, `z` the position along it. A profile point with
    /// `r ≤ 0` is an apex (single vertex), so the band there closes with a fan
    /// (the same convention as `rocket_mesh.rs:21`). A full `360°` sweep wraps
    /// closed; a partial sweep leaves the profile open (no end caps). Built
    /// outward-facing.
    ///
    /// Returns the appended node [`Range<usize>`].
    #[allow(clippy::too_many_arguments)]
    pub fn revolve(
        &mut self,
        profile: &[[f64; 2]],
        center: [f64; 3],
        axis: [f64; 3],
        angle_deg: f64,
        segments: usize,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        if profile.len() < 2 {
            return start..self.nodes.len();
        }
        let seg = segments.max(3);
        let (u, v, w) = axis_frame(axis);
        let c = Vector3::new(center[0], center[1], center[2]);
        let full = (angle_deg.abs() - 360.0).abs() < 1e-6 || angle_deg.abs() >= 360.0;
        // Number of distinct angular stations and rim count per ring.
        let span = angle_deg.to_radians();
        let stations = if full { seg } else { seg + 1 };

        // One ring of vertices per profile point (a single apex vertex when
        // r ≤ 0). `ring_start[p]` is the first node index of profile point p;
        // `ring_len[p]` is 1 for an apex, else `stations`.
        let mut ring_start = Vec::with_capacity(profile.len());
        let mut ring_len = Vec::with_capacity(profile.len());
        for &[r, z] in profile {
            ring_start.push(self.nodes.len() as u32);
            if r <= 1e-9 {
                self.vert(c + w * z);
                ring_len.push(1usize);
            } else {
                for s in 0..stations {
                    let a = if full {
                        s as f64 / seg as f64 * TAU
                    } else {
                        s as f64 / seg as f64 * span
                    };
                    let off = u * (r * a.cos()) + v * (r * a.sin()) + w * z;
                    self.vert(c + off);
                }
                ring_len.push(stations);
            }
        }

        let rim = |s: usize| if full { s % seg } else { s };
        for p in 0..profile.len() - 1 {
            let s0 = ring_start[p];
            let s1 = ring_start[p + 1];
            let next_apex = ring_len[p + 1] == 1;
            let this_apex = ring_len[p] == 1;
            if next_apex {
                let apex = s1;
                for s in 0..seg {
                    self.tri(s0 + rim(s) as u32, s0 + rim(s + 1) as u32, apex, color);
                }
            } else if this_apex {
                let apex = s0;
                for s in 0..seg {
                    self.tri(apex, s1 + rim(s + 1) as u32, s1 + rim(s) as u32, color);
                }
            } else {
                for s in 0..seg {
                    let a0 = s0 + rim(s) as u32;
                    let a1 = s0 + rim(s + 1) as u32;
                    let b0 = s1 + rim(s) as u32;
                    let b1 = s1 + rim(s + 1) as u32;
                    self.quad(a0, a1, b1, b0, color);
                }
            }
        }
        start..self.nodes.len()
    }

    /// **Extrude** a closed 2-D `profile` (a list of `[x, y]` polygon vertices,
    /// ordered counter-clockwise) straight along +z from `z0` to `z1`, with a
    /// bottom cap, a top cap, and a side wall — a prism. Caps are
    /// centroid-fanned (the canonical loft/cap pattern from
    /// `valenx-gears/src/solid.rs:157`). Built outward-facing for a
    /// CCW profile; a clockwise profile inverts the normals.
    ///
    /// Returns the appended node [`Range<usize>`]: `2·n` wall vertices + `2`
    /// cap centroids for an `n`-vertex profile.
    pub fn extrude(
        &mut self,
        profile: &[[f64; 2]],
        z0: f64,
        z1: f64,
        color: [f32; 3],
    ) -> Range<usize> {
        let start = self.nodes.len();
        let n = profile.len();
        if n < 3 {
            return start..self.nodes.len();
        }
        let bottom = self.nodes.len() as u32;
        for p in profile {
            self.vert(Vector3::new(p[0], p[1], z0));
        }
        let top = self.nodes.len() as u32;
        for p in profile {
            self.vert(Vector3::new(p[0], p[1], z1));
        }
        // Side wall.
        for i in 0..n {
            let j = (i + 1) % n;
            self.quad(
                bottom + i as u32,
                bottom + j as u32,
                top + j as u32,
                top + i as u32,
                color,
            );
        }
        // Bottom cap (faces -z): centroid fan wound CW as seen from below.
        let mut cb = Vector3::zeros();
        for p in profile {
            cb += Vector3::new(p[0], p[1], z0);
        }
        cb /= n as f64;
        let cbi = self.vert(cb);
        for i in 0..n {
            let j = (i + 1) % n;
            self.tri(cbi, bottom + j as u32, bottom + i as u32, color);
        }
        // Top cap (faces +z).
        let mut ct = Vector3::zeros();
        for p in profile {
            ct += Vector3::new(p[0], p[1], z1);
        }
        ct /= n as f64;
        let cti = self.vert(ct);
        for i in 0..n {
            let j = (i + 1) % n;
            self.tri(cti, top + i as u32, top + j as u32, color);
        }
        start..self.nodes.len()
    }

    /// Append an **already-built [`Tri3`](ElementType::Tri3) mesh** (its nodes
    /// re-based onto this builder's node array) with one flat `color`,
    /// preserving the per-triangle colour lock-step. This is the bridge for
    /// solids the primitive lathe/extrude path can't express — chiefly the true
    /// involute spur gears emitted by [`valenx_gears::to_solid_spur`] +
    /// tessellated to a [`valenx_mesh::Mesh`] — so a gearbox can place real
    /// toothed gears alongside primitive housings/shafts in one coloured buffer.
    /// Non-`Tri3` element blocks are skipped (the gear tessellation is all
    /// `Tri3`). Triangle node indices are offset by the current node count so
    /// they reference the merged array.
    ///
    /// Returns the [`Range<usize>`] of node indices appended (the whole of
    /// `mesh.nodes`), so the caller can record a [`crate::RigidPart`] over it.
    pub fn append_tri_mesh(&mut self, mesh: &Mesh, color: [f32; 3]) -> Range<usize> {
        let start = self.nodes.len();
        let offset = self.nodes.len() as u32;
        self.nodes.extend_from_slice(&mesh.nodes);
        for blk in &mesh.element_blocks {
            if blk.element_type != ElementType::Tri3 {
                continue;
            }
            for t in blk.connectivity.chunks_exact(3) {
                self.tri(t[0] + offset, t[1] + offset, t[2] + offset, color);
            }
        }
        start..self.nodes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RED: [f32; 3] = [1.0, 0.0, 0.0];

    /// Bake a single-primitive builder and return (mesh, colors, tri_count).
    fn bake(b: MeshBuilder) -> (Mesh, Vec<[f32; 3]>, usize) {
        let tris = b.triangle_count();
        let (m, c) = b.into_mesh_and_colors();
        (m, c, tris)
    }

    /// The renderer's surface-triangle count for a baked mesh.
    fn mesh_tri_count(m: &Mesh) -> usize {
        m.element_blocks.iter().map(|blk| blk.count()).sum()
    }

    /// Every primitive must yield nodes>0, triangles>0, and a colour buffer of
    /// exactly `3 × triangle_count` aligned to the renderer's coloured path.
    fn assert_well_formed(m: &Mesh, colors: &[[f32; 3]], builder_tris: usize) {
        assert!(!m.nodes.is_empty(), "primitive produced no nodes");
        let tc = mesh_tri_count(m);
        assert!(tc > 0, "primitive produced no triangles");
        assert_eq!(tc, builder_tris, "builder/mesh triangle count disagree");
        assert_eq!(
            colors.len(),
            3 * tc,
            "colors must be 3 per triangle (renderer coloured-path layout)"
        );
        // Connectivity must reference only valid node indices.
        for blk in &m.element_blocks {
            for &idx in &blk.connectivity {
                assert!((idx as usize) < m.nodes.len(), "node index out of range");
            }
        }
    }

    #[test]
    fn cylinder_well_formed_and_node_count() {
        let mut b = MeshBuilder::new();
        let seg = 16;
        let range = b.cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0, seg, RED);
        // 2 rings of `seg` + 2 cap centres.
        assert_eq!(range, 0..(2 * seg + 2));
        let (m, c, t) = bake(b);
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn cylinder_caps_present() {
        // Caps present <=> node count is exactly 2 rings + 2 centres, and the
        // triangle count is the wall band (2·seg) plus two fans (seg each).
        let seg = 24;
        let mut b = MeshBuilder::new();
        b.cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 3.0, seg, RED);
        let (m, _c, t) = bake(b);
        assert_eq!(m.nodes.len(), 2 * seg + 2, "expected 2 rings + 2 centres");
        assert_eq!(t, 2 * seg + 2 * seg, "wall band + two cap fans");
    }

    #[test]
    fn cylinder_axis_orientation() {
        // An x-axis cylinder must span x and be thin in z; a z-axis one the
        // opposite. Confirms `axis` actually steers the long direction.
        let seg = 12;
        let half_len = 4.0;
        let r = 1.0;
        let mut bx = MeshBuilder::new();
        bx.cylinder([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], r, 2.0 * half_len, seg, RED);
        let (mx, _) = bx.into_mesh_and_colors();
        let max_x = mx.nodes.iter().map(|p| p.x).fold(f64::MIN, f64::max);
        let max_z = mx.nodes.iter().map(|p| p.z).fold(f64::MIN, f64::max);
        assert!((max_x - half_len).abs() < 1e-6, "x extent = half length");
        assert!((max_z - r).abs() < 1e-6, "z extent = radius");
    }

    #[test]
    fn tube_well_formed_and_node_count() {
        let seg = 20;
        let mut b = MeshBuilder::new();
        let range = b.tube([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 2.0, 4.0, seg, RED);
        // 4 rings (outer/inner × bottom/top) of `seg`.
        assert_eq!(range, 0..(4 * seg));
        let (m, c, t) = bake(b);
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn sphere_vertex_pattern() {
        // UV sphere: 2 poles + (lat-1) interior rings of `lon` each.
        let lat = 8;
        let lon = 12;
        let mut b = MeshBuilder::new();
        let range = b.sphere([0.0, 0.0, 0.0], 3.0, lat, lon, RED);
        let expected = 2 + (lat - 1) * lon;
        assert_eq!(range, 0..expected, "sphere vertex (lat×lon) pattern");
        let (m, c, t) = bake(b);
        assert_eq!(m.nodes.len(), expected);
        assert_well_formed(&m, &c, t);
        // Every vertex sits on the radius (it's a true sphere).
        for p in &m.nodes {
            assert!((p.norm() - 3.0).abs() < 1e-9, "vertex off the sphere radius");
        }
    }

    #[test]
    fn cuboid_well_formed() {
        let mut b = MeshBuilder::new();
        let range = b.cuboid([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], RED);
        assert_eq!(range, 0..8, "8 corners");
        let (m, c, t) = bake(b);
        assert_eq!(t, 12, "12 triangles (2 per face)");
        assert_well_formed(&m, &c, t);
        // Centred: extents are ±half-dims about the centre.
        let min_x = m.nodes.iter().map(|p| p.x).fold(f64::MAX, f64::min);
        let max_x = m.nodes.iter().map(|p| p.x).fold(f64::MIN, f64::max);
        assert!((min_x - 0.0).abs() < 1e-9 && (max_x - 2.0).abs() < 1e-9);
    }

    #[test]
    fn cone_frustum_well_formed() {
        let seg = 16;
        let mut b = MeshBuilder::new();
        // True cone: top radius zero → only the base cap centre is added.
        let range = b.cone([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 0.0, 4.0, seg, RED);
        assert_eq!(range, 0..(2 * seg + 1), "2 rings + 1 base centre (apex top)");
        let (m, c, t) = bake(b);
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn torus_well_formed_and_node_count() {
        let sm = 24;
        let sn = 10;
        let mut b = MeshBuilder::new();
        let range = b.torus([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 5.0, 1.0, sm, sn, RED);
        assert_eq!(range, 0..(sm * sn), "major×minor grid");
        let (m, c, t) = bake(b);
        assert_eq!(t, 2 * sm * sn, "two triangles per grid quad");
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn revolve_rectangle_is_a_cylinder() {
        // Revolving the rectangle profile (r,z) = (R,-h/2)->(R,h/2) a full turn
        // yields ≈ a cylinder side wall: bbox is a [-R,R]² × [-h/2,h/2] box.
        let r = 2.0;
        let h = 6.0;
        let profile = [[r, -h * 0.5], [r, h * 0.5]];
        let mut b = MeshBuilder::new();
        let range = b.revolve(&profile, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 360.0, 32, RED);
        assert!(!range.is_empty());
        let (m, c, t) = bake(b);
        assert_well_formed(&m, &c, t);
        let max_xy = m
            .nodes
            .iter()
            .map(|p| (p.x * p.x + p.y * p.y).sqrt())
            .fold(f64::MIN, f64::max);
        let max_z = m.nodes.iter().map(|p| p.z).fold(f64::MIN, f64::max);
        let min_z = m.nodes.iter().map(|p| p.z).fold(f64::MAX, f64::min);
        assert!((max_xy - r).abs() < 1e-6, "radius matches profile");
        assert!((max_z - h * 0.5).abs() < 1e-9 && (min_z + h * 0.5).abs() < 1e-9);
    }

    #[test]
    fn extrude_square_prism() {
        // A unit square extruded z=0..2 → a box: 12 wall + caps triangles.
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let mut b = MeshBuilder::new();
        let range = b.extrude(&square, 0.0, 2.0, RED);
        // 2·n wall verts + 2 cap centroids, n = 4.
        assert_eq!(range, 0..(2 * 4 + 2));
        let (m, c, t) = bake(b);
        // Side wall = 2·n tris; two caps = 2·n tris.
        assert_eq!(t, 2 * 4 + 2 * 4);
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn multi_primitive_ranges_are_contiguous() {
        // Build N parts; the end of part k must equal the start of part k+1, and
        // the final range end must equal the total node count. This is the
        // "build N parts, each a colour + a returned range" contract.
        let mut b = MeshBuilder::new();
        let r0 = b.cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 2.0, 12, [1.0, 0.0, 0.0]);
        let r1 = b.sphere([3.0, 0.0, 0.0], 1.0, 6, 8, [0.0, 1.0, 0.0]);
        let r2 = b.cuboid([0.0, 3.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 1.0]);
        let r3 = b.torus([0.0, 0.0, 3.0], [0.0, 0.0, 1.0], 2.0, 0.5, 12, 8, [1.0, 1.0, 0.0]);
        assert_eq!(r0.start, 0);
        assert_eq!(r0.end, r1.start, "part0.end == part1.start");
        assert_eq!(r1.end, r2.start, "part1.end == part2.start");
        assert_eq!(r2.end, r3.start, "part2.end == part3.start");
        let total = b.node_count();
        assert_eq!(r3.end, total, "last range covers up to the node count");
        let (m, c, t) = bake(b);
        assert_eq!(m.nodes.len(), total);
        assert_well_formed(&m, &c, t);
    }

    #[test]
    fn append_tri_mesh_rebases_and_colours_per_triangle() {
        // A 2-triangle quad mesh appended after a cuboid must: re-base its node
        // indices onto the builder's array, push exactly 3 colours per appended
        // triangle, and report the appended node range. This is the bridge the
        // gear producers use to drop a tessellated involute solid into the
        // coloured buffer.
        let mut ext = Mesh::new("ext");
        ext.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 0, 2, 3];
        ext.element_blocks.push(blk);

        let mut b = MeshBuilder::new();
        let cuboid = b.cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], RED); // 8 verts, 12 tris
        let green = [0.0, 1.0, 0.0];
        let r = b.append_tri_mesh(&ext, green);
        // Appended range is the 4 external verts, right after the cuboid's 8.
        assert_eq!(r, cuboid.end..(cuboid.end + 4));
        let (m, c, t) = bake(b);
        assert_well_formed(&m, &c, t);
        assert_eq!(t, 12 + 2, "cuboid 12 tris + appended 2 tris");
        // The last 2 triangles (6 colour entries) are the appended mesh's green.
        assert!(c[36..].iter().all(|&x| x == green), "appended tris are green");
    }

    #[test]
    fn colors_track_per_part_color() {
        // Two parts in two colours: the colour buffer's two halves must hold the
        // respective flat colours (3 entries per triangle, triangle-major).
        let mut b = MeshBuilder::new();
        let red = [1.0, 0.0, 0.0];
        let blue = [0.0, 0.0, 1.0];
        b.cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], red); // 12 tris → 36 colours
        b.cuboid([2.0, 0.0, 0.0], [1.0, 1.0, 1.0], blue); // 12 tris → 36 colours
        let (_m, c, _t) = bake(b);
        assert_eq!(c.len(), 72);
        assert!(c[..36].iter().all(|&x| x == red), "first part all red");
        assert!(c[36..].iter().all(|&x| x == blue), "second part all blue");
    }
}

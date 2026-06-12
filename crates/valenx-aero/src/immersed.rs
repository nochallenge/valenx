//! The immersed-boundary method — voxelizing an arbitrary triangle
//! mesh into the Cartesian grid.
//!
//! # The idea
//!
//! Body-fitted CFD meshes a curved car / wing exactly, but generating
//! that mesh for an arbitrary watertight body is hard. The
//! **immersed-boundary method (IBM)** sidesteps it: the grid stays a
//! plain uniform Cartesian box, and the body is represented by
//! *flagging* the cells it occupies. Cells inside the body are
//! `Solid`, cells outside are `Fluid`, and the band of fluid cells
//! adjacent to a solid cell are the `Cut` cells where the no-slip wall
//! condition is enforced (direct-forcing IBM — the velocity in a cut
//! cell is driven toward the wall velocity).
//!
//! This is what lets a user *drop* a triangle mesh into the flow and
//! run, with no body-fitted mesher in the loop.
//!
//! # Inside / outside classification
//!
//! Each pressure cell centre is classified by **ray casting**: shoot a
//! ray in `+x` and count how many body triangles it crosses; an odd
//! count means the point is inside the body. The implementation uses
//! the Möller–Trumbore ray/triangle intersection and a robust
//! parity-of-crossings count. For extra robustness against a ray that
//! grazes a shared triangle edge, the classifier casts **three** rays
//! (one per axis) and takes the majority vote — a single grazing
//! artefact cannot then flip the result.
//!
//! # Cut-cell geometry
//!
//! By default the voxelizer also reconstructs **cut-cell geometry**
//! ([`crate::cutcell`]) for every cell the body surface passes through:
//! the sub-cell fluid volume fraction, the apertured (partial) areas of
//! the six Cartesian faces, and the true clipped wall face (a polygon
//! with an exact area + outward normal). The finite-volume solver and
//! the force integration consult that geometry so the wall is the real
//! body surface, not the voxel staircase. The legacy staircased path is
//! still selectable via [`crate::cutcell::WallMethod::Staircase`].
//!
//! # Honest scope
//!
//! Robust voxelization of a *watertight* mesh. A mesh with holes or
//! flipped facets degrades gracefully (the majority vote tolerates a
//! minority of bad crossings) but is not guaranteed correct — the
//! body should be a closed shell, which any CAD tessellation produces.
//! The cut-cell reconstruction closes the *staircased-wall* accuracy
//! gap; what remains is a body-fitted unstructured mesher with a
//! near-wall prism layer — a documented Tier-3 residue.

use nalgebra::Vector3;
use rayon::prelude::*;

use crate::cutcell::{compute_cut_geometry, CellBox, CellGeometry, WallMethod};
use crate::geometry::{TriMesh, Triangle};
use crate::grid::Grid3;

/// The role a pressure cell plays relative to the immersed body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellTag {
    /// A normal fluid cell, away from the body — the Navier-Stokes
    /// equations are solved here unchanged.
    Fluid,
    /// A fluid cell immediately adjacent to the body — the immersed-
    /// boundary no-slip wall forcing is applied here.
    Cut,
    /// A cell inside the body — frozen to the (solid) body velocity,
    /// excluded from the fluid solve.
    Solid,
}

/// The voxelized representation of the body inside the grid — one
/// [`CellTag`] per pressure cell plus the wall data the solver needs.
///
/// Under the default cut-cell method ([`WallMethod::CutCell`]) the body
/// also carries per-cell [`CellGeometry`] — the fluid volume fraction,
/// the six apertured Cartesian-face areas, and the true clipped wall
/// face — for every cut cell; fully-solid and fully-fluid cells get the
/// trivial geometry. Under [`WallMethod::Staircase`] the geometry array
/// holds the staircased equivalent (volume fraction 0 or 1, apertures 0
/// or 1, no sub-cell wall) so callers can read it uniformly.
#[derive(Clone, Debug)]
pub struct ImmersedBody {
    /// The grid this voxelization is for.
    pub grid: Grid3,
    /// One tag per pressure cell, x-fastest then y then z.
    pub tags: Vec<CellTag>,
    /// Signed-distance-like proximity: the (positive) distance from a
    /// cut/fluid cell centre to the nearest body triangle, used by the
    /// wall functions for y+. `0` inside solid cells.
    pub wall_distance: Vec<f64>,
    /// The wall treatment this body was voxelized for.
    pub method: WallMethod,
    /// Per-cell cut-cell geometry — volume fraction, face apertures and
    /// the clipped wall face. One entry per pressure cell.
    pub geometry: Vec<CellGeometry>,
}

impl ImmersedBody {
    /// Linear index of pressure cell `(i, j, k)`.
    #[inline]
    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.grid.nx * (j + self.grid.ny * k)
    }

    /// The tag of pressure cell `(i, j, k)`.
    #[inline]
    pub fn tag(&self, i: usize, j: usize, k: usize) -> CellTag {
        self.tags[self.idx(i, j, k)]
    }

    /// True if cell `(i, j, k)` is solid (inside the body).
    #[inline]
    pub fn is_solid(&self, i: usize, j: usize, k: usize) -> bool {
        self.tag(i, j, k) == CellTag::Solid
    }

    /// True if cell `(i, j, k)` is a fluid or cut cell — i.e. part of
    /// the flow domain the solver advances.
    #[inline]
    pub fn is_fluid_like(&self, i: usize, j: usize, k: usize) -> bool {
        self.tag(i, j, k) != CellTag::Solid
    }

    /// True if cell `(i, j, k)` is a cut cell (adjacent to the body).
    #[inline]
    pub fn is_cut(&self, i: usize, j: usize, k: usize) -> bool {
        self.tag(i, j, k) == CellTag::Cut
    }

    /// Wall distance of cell `(i, j, k)` — distance to the nearest
    /// body triangle.
    #[inline]
    pub fn wall_distance(&self, i: usize, j: usize, k: usize) -> f64 {
        self.wall_distance[self.idx(i, j, k)]
    }

    /// The cut-cell geometry of cell `(i, j, k)` — its fluid volume
    /// fraction, the six Cartesian-face apertures, and the clipped wall
    /// face.
    #[inline]
    pub fn geometry(&self, i: usize, j: usize, k: usize) -> &CellGeometry {
        &self.geometry[self.idx(i, j, k)]
    }

    /// The fluid volume fraction `θ ∈ [0, 1]` of cell `(i, j, k)` — the
    /// share of the cell box that is outside the body. `0` for a solid
    /// cell, `1` for a clear fluid cell, partial for a cut cell.
    #[inline]
    pub fn volume_fraction(&self, i: usize, j: usize, k: usize) -> f64 {
        self.geometry[self.idx(i, j, k)].volume_fraction
    }

    /// The clipped body-surface wall face of cut cell `(i, j, k)` —
    /// area zero for a non-cut cell.
    #[inline]
    pub fn cut_face(&self, i: usize, j: usize, k: usize) -> &crate::cutcell::CutFace {
        &self.geometry[self.idx(i, j, k)].cut_face
    }

    /// True if this body was voxelized with the cut-cell method.
    #[inline]
    pub fn uses_cut_cells(&self) -> bool {
        self.method.is_cut_cell()
    }

    /// The total clipped wall area summed over every cut cell — under
    /// the cut-cell method this approximates the body's true surface
    /// area (it tiles the real surface), a useful sanity diagnostic.
    pub fn cut_wall_area(&self) -> f64 {
        self.geometry.iter().map(|g| g.cut_face.area).sum()
    }

    /// Count cells of each tag — `(fluid, cut, solid)`.
    pub fn tag_counts(&self) -> (usize, usize, usize) {
        let mut f = 0;
        let mut c = 0;
        let mut s = 0;
        for &t in &self.tags {
            match t {
                CellTag::Fluid => f += 1,
                CellTag::Cut => c += 1,
                CellTag::Solid => s += 1,
            }
        }
        (f, c, s)
    }

    /// The body volume estimate.
    ///
    /// Under the cut-cell method this sums `(1 − θ)·cellVolume` over
    /// every cell — the partial solid volume of each cut cell is
    /// counted, so the estimate is far better than a voxel count.
    /// Under the staircased method it is the solid-cell count times the
    /// cell volume.
    pub fn solid_volume(&self) -> f64 {
        let cell_vol = self.grid.dx() * self.grid.dy() * self.grid.dz();
        if self.method.is_cut_cell() {
            self.geometry
                .iter()
                .map(|g| (1.0 - g.volume_fraction) * cell_vol)
                .sum()
        } else {
            let (_, _, s) = self.tag_counts();
            s as f64 * cell_vol
        }
    }
}

/// Möller–Trumbore ray/triangle intersection.
///
/// Returns the ray parameter `t` (`origin + t·dir` is the hit point)
/// if the ray crosses the triangle at `t > eps`, else `None`. `dir`
/// need not be normalised.
fn ray_triangle(origin: Vector3<f64>, dir: Vector3<f64>, tri: &Triangle) -> Option<f64> {
    const EPS: f64 = 1e-12;
    let e1 = tri.b - tri.a;
    let e2 = tri.c - tri.a;
    let pvec = dir.cross(&e2);
    let det = e1.dot(&pvec);
    if det.abs() < EPS {
        return None; // ray parallel to the triangle plane
    }
    let inv_det = 1.0 / det;
    let tvec = origin - tri.a;
    let u = tvec.dot(&pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let qvec = tvec.cross(&e1);
    let v = dir.dot(&qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(&qvec) * inv_det;
    if t > EPS {
        Some(t)
    } else {
        None
    }
}

/// Count how many of the mesh's triangles a ray from `origin` along
/// `dir` crosses. An odd count means `origin` is inside a watertight
/// body.
fn crossing_count(origin: Vector3<f64>, dir: Vector3<f64>, mesh: &TriMesh) -> usize {
    let mut n = 0;
    for tri in &mesh.triangles {
        if ray_triangle(origin, dir, tri).is_some() {
            n += 1;
        }
    }
    n
}

/// Classify a single point as inside the body by a three-ray majority
/// vote — one ray per axis. A grazing artefact on one ray cannot flip
/// the verdict.
pub fn point_inside(point: Vector3<f64>, mesh: &TriMesh) -> bool {
    // Slightly irrational ray directions reduce the chance of a ray
    // lying exactly in a triangle's plane.
    let dirs = [
        Vector3::new(1.0, 0.013, 0.007),
        Vector3::new(0.011, 1.0, 0.009),
        Vector3::new(0.008, 0.012, 1.0),
    ];
    let mut votes = 0;
    for d in dirs {
        if crossing_count(point, d, mesh) % 2 == 1 {
            votes += 1;
        }
    }
    votes >= 2
}

/// The squared distance from a point to a triangle (closest-point on
/// the triangle, clamped to its interior / edges / vertices).
fn point_tri_dist_sq(p: Vector3<f64>, tri: &Triangle) -> f64 {
    // Ericson, Real-Time Collision Detection — closest point on a
    // triangle to a point.
    let ab = tri.b - tri.a;
    let ac = tri.c - tri.a;
    let ap = p - tri.a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.norm_squared();
    }
    let bp = p - tri.b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return bp.norm_squared();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (ap - v * ab).norm_squared();
    }
    let cp = p - tri.c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return cp.norm_squared();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (ap - w * ac).norm_squared();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (tri.b + w * (tri.c - tri.b))).norm_squared();
    }
    // Inside the face region — project onto the plane.
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    (p - (tri.a + ab * v + ac * w)).norm_squared()
}

/// The distance from a point to the nearest triangle in the mesh.
pub fn distance_to_body(point: Vector3<f64>, mesh: &TriMesh) -> f64 {
    let mut best = f64::INFINITY;
    for tri in &mesh.triangles {
        let d2 = point_tri_dist_sq(point, tri);
        if d2 < best {
            best = d2;
        }
    }
    best.sqrt()
}

/// The small-cell volume-fraction floor.
///
/// A cut cell with a fluid volume fraction below this floor is a *tiny
/// cell* — its arbitrarily small control volume would destabilise the
/// segregated SIMPLE solver (the explicit pressure-correction
/// coefficients scale like `1/θ`). Such a cell is **merged into the
/// solid** (re-tagged `Solid`) — a standard small-cell remedy: the tiny
/// fluid sliver is absorbed, removing the stiff cell entirely. Cut
/// cells above the floor are kept but their effective volume is clamped
/// to at least `SMALL_CELL_FLOOR` of the cell volume in the
/// discretisation (a flux-redistribution-style floor), so the
/// time/relaxation scaling never blows up. See [`stabilized_volume`].
pub const SMALL_CELL_FLOOR: f64 = 0.05;

/// The effective fluid volume of a cut cell for the finite-volume
/// discretisation — the geometric volume fraction, floored at
/// [`SMALL_CELL_FLOOR`] so a small cut cell does not produce a stiff
/// (near-singular) momentum / pressure coefficient. This is the
/// flux-redistribution-style floor that pairs with the cell-merging of
/// genuinely tiny cells.
#[inline]
pub fn stabilized_volume(volume_fraction: f64) -> f64 {
    volume_fraction.max(SMALL_CELL_FLOOR)
}

/// Voxelize a body mesh into a grid with the default (cut-cell) wall
/// method.
///
/// See [`voxelize_with`] for the full description; this is the
/// convenience form used by [`crate::WindTunnel`].
///
/// `mesh` should be expressed in the grid's world coordinates — place
/// the body before calling.
pub fn voxelize(grid: &Grid3, mesh: &TriMesh) -> ImmersedBody {
    voxelize_with(grid, mesh, WallMethod::default())
}

/// Voxelize a body mesh into a grid with an explicit wall method.
///
/// Pass 1 classifies every pressure-cell centre inside / outside the
/// body by the three-ray majority vote. Then:
///
/// - Under [`WallMethod::CutCell`] the cut-cell geometry
///   ([`crate::cutcell`]) is reconstructed for the band of cells the
///   body *surface* clips — the fluid volume fraction, the six
///   apertured face areas, and the true clipped wall face. A cell the
///   surface passes through is tagged `Cut`; a cell wholly inside the
///   body is `Solid`; a cell wholly outside is `Fluid`. Genuinely tiny
///   cut cells (volume fraction below [`SMALL_CELL_FLOOR`]) are
///   *merged into the solid* for numerical stability.
/// - Under [`WallMethod::Staircase`] the legacy path runs — a fluid
///   cell with a solid face-neighbour is promoted to `Cut`, and the
///   geometry array is filled with the staircased equivalent (volume
///   fraction 0 or 1, apertures 0 or 1, no sub-cell wall).
///
/// Either way the wall distance is recorded for the cut and near-fluid
/// cells (used by the turbulence wall functions for `y+`).
pub fn voxelize_with(grid: &Grid3, mesh: &TriMesh, method: WallMethod) -> ImmersedBody {
    let n = grid.cell_count();
    let mut tags = vec![CellTag::Fluid; n];
    let mut wall_distance = vec![f64::INFINITY; n];
    let mut geometry = vec![CellGeometry::full_fluid(); n];

    let bb = mesh.aabb();

    if method.is_cut_cell() {
        // Cut-cell pass — reconstruct the sub-cell geometry for the band
        // of cells the body surface clips, and tag every cell from it.
        //
        // The per-cell geometry reconstruction is independent across
        // cells, so it is computed in parallel over the rayon pool
        // (z-plane chunks). A cell whose box (slightly padded so a
        // surface grazing the face is not missed) overlaps the body
        // AABB is a candidate; `compute_cut_geometry` then decides — via
        // the exact triangle clip — whether the surface genuinely passes
        // through it. A non-candidate cell is wholly fluid.
        let margin = grid.dx().max(grid.dy()).max(grid.dz());
        let plane = grid.nx * grid.ny;
        geometry
            .par_chunks_mut(plane.max(1))
            .enumerate()
            .for_each(|(k, slab)| {
                for j in 0..grid.ny {
                    for i in 0..grid.nx {
                        let local = i + grid.nx * j;
                        let cell = CellBox::of_cell(grid, i, j, k);
                        let candidate = match bb {
                            Some(b) => {
                                cell.max.x >= b.min.x - margin
                                    && cell.min.x <= b.max.x + margin
                                    && cell.max.y >= b.min.y - margin
                                    && cell.min.y <= b.max.y + margin
                                    && cell.max.z >= b.min.z - margin
                                    && cell.min.z <= b.max.z + margin
                            }
                            None => false,
                        };
                        slab[local] = if candidate {
                            compute_cut_geometry(&cell, mesh)
                        } else {
                            CellGeometry::full_fluid()
                        };
                    }
                }
            });
        // Tag every cell from the reconstructed geometry — a cell the
        // surface clips is `Cut`, an all-solid-fraction cell is `Solid`,
        // the rest `Fluid`.
        for (idx, g) in geometry.iter().enumerate() {
            if g.cut_face.has_wall() {
                tags[idx] = CellTag::Cut;
            } else if g.volume_fraction <= 0.0 {
                tags[idx] = CellTag::Solid;
                wall_distance[idx] = 0.0;
            } else {
                tags[idx] = CellTag::Fluid;
            }
        }

        // Small-cell stabilisation — cell merging with wall transfer.
        //
        // A cut cell whose fluid volume fraction is below the floor is a
        // tiny cell — its near-singular control volume would destabilise
        // the segregated SIMPLE solver. It is *merged*: re-tagged
        // `Solid` (so the flow solve sees a clean staircased wall on the
        // face toward its open side, and never tries to advance the
        // stiff sliver), but its **clipped wall face is transferred to**
        // the adjacent larger fluid / cut cell on the open side. That
        // neighbour is promoted to `Cut` and carries the merged wall, so
        // the surface area and the pressure / shear force are preserved
        // — the wall is not discarded, only re-homed onto a cell big
        // enough to solve stably. Without the transfer a body face that
        // grazes a grid plane (every cell on it tiny) would lose its
        // entire wall and its drag contribution.
        let small: Vec<usize> = (0..n)
            .filter(|&idx| {
                tags[idx] == CellTag::Cut && geometry[idx].volume_fraction < SMALL_CELL_FLOOR
            })
            .collect();
        for &idx in &small {
            let k = idx / (grid.nx * grid.ny);
            let rem = idx % (grid.nx * grid.ny);
            let j = rem / grid.nx;
            let i = rem % grid.nx;
            let wall = geometry[idx].cut_face;
            // Pick the merge partner — the face-neighbour most aligned
            // with the outward wall normal (the open / fluid side) that
            // is not solid. Fall back to any fluid-like neighbour.
            let nbrs: [(i32, i32, i32); 6] = [
                (1, 0, 0),
                (-1, 0, 0),
                (0, 1, 0),
                (0, -1, 0),
                (0, 0, 1),
                (0, 0, -1),
            ];
            let mut best: Option<(usize, f64)> = None;
            for (di, dj, dk) in nbrs {
                let (ni, nj, nk) = (i as i32 + di, j as i32 + dj, k as i32 + dk);
                if ni < 0
                    || nj < 0
                    || nk < 0
                    || ni >= grid.nx as i32
                    || nj >= grid.ny as i32
                    || nk >= grid.nz as i32
                {
                    continue;
                }
                let nidx = ni as usize + grid.nx * (nj as usize + grid.ny * nk as usize);
                // A small cell about to be merged is not a valid host.
                if tags[nidx] == CellTag::Solid || geometry[nidx].volume_fraction < SMALL_CELL_FLOOR
                {
                    continue;
                }
                let align = (di as f64) * wall.normal.x
                    + (dj as f64) * wall.normal.y
                    + (dk as f64) * wall.normal.z;
                if best.map(|(_, a)| align > a).unwrap_or(true) {
                    best = Some((nidx, align));
                }
            }
            // Re-tag the tiny cell solid; transfer its wall to the host.
            tags[idx] = CellTag::Solid;
            geometry[idx] = CellGeometry::full_solid();
            wall_distance[idx] = 0.0;
            if let Some((host, _)) = best {
                let merged = geometry[host].cut_face.merged(&wall);
                geometry[host].cut_face = merged;
                if tags[host] == CellTag::Fluid {
                    tags[host] = CellTag::Cut;
                }
            }
        }
    } else {
        // Staircased classification — every cell centre is tested
        // inside / outside the body by the three-ray majority vote; an
        // inside cell becomes `Solid`. The body AABB bounds the work.
        // The per-cell classification is independent, so it runs over
        // the rayon pool by z-plane.
        let plane = grid.nx * grid.ny;
        tags.par_chunks_mut(plane.max(1))
            .zip(wall_distance.par_chunks_mut(plane.max(1)))
            .enumerate()
            .for_each(|(k, (tag_slab, wd_slab))| {
                for j in 0..grid.ny {
                    for i in 0..grid.nx {
                        let local = i + grid.nx * j;
                        let (cx, cy, cz) = grid.cell_centre(i, j, k);
                        let c = Vector3::new(cx, cy, cz);
                        let inside = match bb {
                            Some(b) => {
                                !(c.x < b.min.x
                                    || c.x > b.max.x
                                    || c.y < b.min.y
                                    || c.y > b.max.y
                                    || c.z < b.min.z
                                    || c.z > b.max.z)
                                    && point_inside(c, mesh)
                            }
                            None => false,
                        };
                        if inside {
                            tag_slab[local] = CellTag::Solid;
                            wd_slab[local] = 0.0;
                        }
                    }
                }
            });
        // Promote fluid cells adjacent to a solid cell to Cut, exactly
        // the legacy behaviour.
        let tags_snapshot = tags.clone();
        let tag_of = |i: usize, j: usize, k: usize| -> CellTag {
            tags_snapshot[i + grid.nx * (j + grid.ny * k)]
        };
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    let idx = i + grid.nx * (j + grid.ny * k);
                    if tags_snapshot[idx] != CellTag::Fluid {
                        continue;
                    }
                    let mut adjacent_solid = false;
                    if i > 0 && tag_of(i - 1, j, k) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if i + 1 < grid.nx && tag_of(i + 1, j, k) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if j > 0 && tag_of(i, j - 1, k) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if j + 1 < grid.ny && tag_of(i, j + 1, k) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if k > 0 && tag_of(i, j, k - 1) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if k + 1 < grid.nz && tag_of(i, j, k + 1) == CellTag::Solid {
                        adjacent_solid = true;
                    }
                    if adjacent_solid {
                        tags[idx] = CellTag::Cut;
                    }
                }
            }
        }
        // Staircased geometry — 0/1 volume fractions and apertures, no
        // sub-cell wall, so callers can read `geometry` uniformly.
        for (idx, g) in geometry.iter_mut().enumerate() {
            *g = if tags[idx] == CellTag::Solid {
                CellGeometry::full_solid()
            } else {
                CellGeometry::full_fluid()
            };
        }
    }

    // Wall distance for cut cells (and their fluid neighbours, useful
    // for y+). Only computed near the body to keep it cheap; far-field
    // fluid cells keep INFINITY.
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let idx = i + grid.nx * (j + grid.ny * k);
                if tags[idx] == CellTag::Solid {
                    wall_distance[idx] = 0.0;
                    continue;
                }
                let near = match bb {
                    Some(b) => {
                        let (cx, cy, cz) = grid.cell_centre(i, j, k);
                        let margin = 3.0 * grid.dx().max(grid.dy()).max(grid.dz());
                        cx >= b.min.x - margin
                            && cx <= b.max.x + margin
                            && cy >= b.min.y - margin
                            && cy <= b.max.y + margin
                            && cz >= b.min.z - margin
                            && cz <= b.max.z + margin
                    }
                    None => false,
                };
                if near {
                    let (cx, cy, cz) = grid.cell_centre(i, j, k);
                    wall_distance[idx] = distance_to_body(Vector3::new(cx, cy, cz), mesh);
                }
            }
        }
    }

    ImmersedBody {
        grid: *grid,
        tags,
        wall_distance,
        method,
        geometry,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{box_body, sphere_body};

    #[test]
    fn ray_triangle_hits_a_facing_triangle() {
        let tri = Triangle::new(
            Vector3::new(0.0, -1.0, -1.0),
            Vector3::new(0.0, 1.0, -1.0),
            Vector3::new(0.0, 0.0, 1.0),
        );
        // Ray from -x toward +x through the triangle.
        let t = ray_triangle(
            Vector3::new(-2.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            &tri,
        );
        assert!(t.is_some());
        assert!((t.unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn ray_triangle_misses_a_triangle_outside_the_ray() {
        let tri = Triangle::new(
            Vector3::new(0.0, 10.0, 10.0),
            Vector3::new(0.0, 11.0, 10.0),
            Vector3::new(0.0, 10.0, 11.0),
        );
        let t = ray_triangle(
            Vector3::new(-2.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            &tri,
        );
        assert!(t.is_none());
    }

    #[test]
    fn point_inside_a_box_is_classified_inside() {
        // The defining IBM correctness test: a known body, a known
        // inside point and a known outside point.
        let body = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        // Centre of the box — unambiguously inside.
        assert!(point_inside(Vector3::new(0.0, 0.0, 0.0), &body));
        // Off-centre but still inside.
        assert!(point_inside(Vector3::new(0.4, -0.6, 0.3), &body));
        // Well outside the box.
        assert!(!point_inside(Vector3::new(5.0, 0.0, 0.0), &body));
        assert!(!point_inside(Vector3::new(0.0, 0.0, 3.0), &body));
    }

    #[test]
    fn point_inside_a_sphere_matches_the_radius_test() {
        // For a sphere the inside test must agree with |p−c| < r.
        let centre = Vector3::new(2.0, 3.0, 4.0);
        let r = 1.5;
        let body = sphere_body(centre, r, 32, 64);
        for p in [
            Vector3::new(2.0, 3.0, 4.0),       // centre, inside
            Vector3::new(2.0 + 1.0, 3.0, 4.0), // r=1 < 1.5, inside
            Vector3::new(2.0 + 3.0, 3.0, 4.0), // r=3 > 1.5, outside
            Vector3::new(2.0, 3.0 + 0.5, 4.0), // inside
        ] {
            let analytic = (p - centre).norm() < r;
            assert_eq!(
                point_inside(p, &body),
                analytic,
                "inside test disagreed for {p:?}"
            );
        }
    }

    #[test]
    fn distance_to_body_is_zero_on_the_surface_and_grows_away() {
        let body = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        // A point on the +x face.
        let on_face = distance_to_body(Vector3::new(1.0, 0.0, 0.0), &body);
        assert!(on_face < 1e-9, "on-surface distance should be ~0");
        // A point 2 units off the +x face.
        let off = distance_to_body(Vector3::new(3.0, 0.0, 0.0), &body);
        assert!((off - 2.0).abs() < 1e-9, "off-face distance should be 2");
    }

    #[test]
    fn voxelize_a_box_flags_a_solid_core() {
        // A box body centred in a grid: the interior cells must come
        // out Solid, the cells around them Cut, the far cells Fluid.
        let grid = Grid3::new(20, 20, 20, 10.0, 10.0, 10.0, -5.0, -5.0, -5.0);
        let body = box_body(Vector3::new(-2.0, -2.0, -2.0), Vector3::new(2.0, 2.0, 2.0));
        let vb = voxelize(&grid, &body);
        let (f, c, s) = vb.tag_counts();
        assert!(s > 0, "the box must produce solid cells");
        assert!(c > 0, "the box must produce a cut band");
        assert!(f > 0, "there must be far-field fluid");
        assert_eq!(f + c + s, grid.cell_count());
        // The voxelized solid volume should be in the ballpark of the
        // true 4×4×4 = 64 box volume.
        let vol = vb.solid_volume();
        assert!(vol > 30.0 && vol < 110.0, "solid volume {vol} unreasonable");
        // A cell at the grid corner is far-field fluid.
        assert_eq!(vb.tag(0, 0, 0), CellTag::Fluid);
        // The centre cell is solid.
        assert!(vb.is_solid(10, 10, 10));
    }

    #[test]
    fn voxelize_empty_mesh_is_all_fluid() {
        let grid = Grid3::at_origin(6, 6, 6, 1.0, 1.0, 1.0);
        let vb = voxelize(&grid, &TriMesh::new());
        let (f, c, s) = vb.tag_counts();
        assert_eq!(s, 0);
        assert_eq!(c, 0);
        assert_eq!(f, grid.cell_count());
    }

    #[test]
    fn cut_cells_have_finite_wall_distance() {
        let grid = Grid3::new(16, 16, 16, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), 1.5, 24, 48);
        let vb = voxelize(&grid, &body);
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    if vb.is_cut(i, j, k) {
                        let d = vb.wall_distance(i, j, k);
                        assert!(d.is_finite() && d >= 0.0, "cut cell wall dist bad: {d}");
                    }
                }
            }
        }
    }

    #[test]
    fn default_voxelize_uses_the_cut_cell_method() {
        // The plain `voxelize` entry point must default to cut-cell —
        // the cut-cell path is the real wind-tunnel default.
        let grid = Grid3::new(16, 16, 16, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), 1.5, 24, 48);
        let vb = voxelize(&grid, &body);
        assert!(vb.uses_cut_cells());
        assert_eq!(vb.method, WallMethod::CutCell);
    }

    #[test]
    fn cut_cell_volume_fractions_are_partial_on_the_surface() {
        // Under the cut-cell method every cut cell on a curved body must
        // carry a genuinely partial fluid volume fraction (strictly
        // between 0 and 1), a real cut face, and apertures in range.
        let grid = Grid3::new(20, 20, 20, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), 1.5, 32, 64);
        let vb = voxelize(&grid, &body);
        let mut partial_seen = false;
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    let g = vb.geometry(i, j, k);
                    // Volume fraction always in [0, 1].
                    assert!(
                        (0.0..=1.0).contains(&g.volume_fraction),
                        "volume fraction {} out of range",
                        g.volume_fraction
                    );
                    for ap in g.face_aperture {
                        assert!((0.0..=1.0).contains(&ap));
                    }
                    match vb.tag(i, j, k) {
                        CellTag::Fluid => {
                            assert_eq!(g.volume_fraction, 1.0);
                        }
                        CellTag::Solid => {
                            assert_eq!(g.volume_fraction, 0.0);
                        }
                        CellTag::Cut => {
                            assert!(g.cut_face.has_wall(), "a cut cell must have a wall");
                            if g.volume_fraction > 0.05 && g.volume_fraction < 0.95 {
                                partial_seen = true;
                            }
                        }
                    }
                }
            }
        }
        assert!(
            partial_seen,
            "a voxelized sphere must produce genuinely partial cut cells"
        );
    }

    #[test]
    fn cut_cell_solid_volume_beats_the_voxel_count_for_a_sphere() {
        // The cut-cell solid volume sums the partial (1−θ) solid
        // fraction of every cut cell, so it must estimate the analytic
        // sphere volume far better than the staircased voxel count.
        let r = 1.5;
        let grid = Grid3::new(24, 24, 24, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), r, 48, 96);
        let analytic = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;

        let cut = voxelize_with(&grid, &body, WallMethod::CutCell);
        let stair = voxelize_with(&grid, &body, WallMethod::Staircase);
        let err_cut = (cut.solid_volume() - analytic).abs() / analytic;
        let err_stair = (stair.solid_volume() - analytic).abs() / analytic;
        assert!(
            err_cut <= err_stair,
            "cut-cell volume error {err_cut} should beat staircased {err_stair}"
        );
        // The cut-cell estimate should be good to a few percent.
        assert!(err_cut < 0.05, "cut-cell sphere volume error {err_cut}");
    }

    #[test]
    fn cut_wall_area_approximates_the_true_surface_area() {
        // The summed cut-face area tiles the real body surface, so for a
        // sphere it must approximate the analytic 4πr².
        let r = 1.5;
        let grid = Grid3::new(24, 24, 24, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), r, 48, 96);
        let vb = voxelize(&grid, &body);
        let analytic = 4.0 * std::f64::consts::PI * r * r;
        let rel = (vb.cut_wall_area() - analytic).abs() / analytic;
        assert!(
            rel < 0.10,
            "cut wall area {} vs analytic {analytic} (rel {rel})",
            vb.cut_wall_area()
        );
    }

    #[test]
    fn small_cell_merge_keeps_the_wall_on_a_grid_aligned_box() {
        // A grid-aligned box can place a face just inside a cell, making
        // every cell on that face a tiny cut cell. The small-cell
        // stabilisation merges those tiny cells into the solid but must
        // *transfer their wall* to the neighbour — the box's six faces
        // must all still carry cut-face wall area, none lost.
        let grid = Grid3::new(40, 40, 40, 10.0, 10.0, 10.0, -5.0, -5.0, -5.0);
        let body = box_body(Vector3::new(-2.0, -2.0, -2.0), Vector3::new(2.0, 2.0, 2.0));
        let vb = voxelize(&grid, &body);
        // Sum the cut-face area-vector projected onto each of the six
        // axis directions — every face of the box must contribute.
        let mut proj = [0.0f64; 6];
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    if vb.tag(i, j, k) != CellTag::Cut {
                        continue;
                    }
                    let cf = vb.cut_face(i, j, k);
                    if !cf.has_wall() {
                        continue;
                    }
                    let av = cf.normal * cf.area;
                    if av.x > 0.0 {
                        proj[0] += av.x;
                    } else {
                        proj[1] += -av.x;
                    }
                    if av.y > 0.0 {
                        proj[2] += av.y;
                    } else {
                        proj[3] += -av.y;
                    }
                    if av.z > 0.0 {
                        proj[4] += av.z;
                    } else {
                        proj[5] += -av.z;
                    }
                }
            }
        }
        // Each of the box's six 4×4 = 16 m² faces must be substantially
        // represented — none wiped out by the small-cell merge.
        for (f, &p) in proj.iter().enumerate() {
            assert!(
                p > 8.0,
                "box face {f} projected wall area {p} too small — \
                 small-cell merge lost the wall"
            );
        }
    }

    #[test]
    fn staircase_method_still_available_and_has_no_sub_cell_wall() {
        // The legacy staircased method must remain selectable and behave
        // as before — 0/1 volume fractions, no sub-cell cut face.
        let grid = Grid3::new(16, 16, 16, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let body = sphere_body(Vector3::zeros(), 1.5, 24, 48);
        let vb = voxelize_with(&grid, &body, WallMethod::Staircase);
        assert!(!vb.uses_cut_cells());
        let (f, c, s) = vb.tag_counts();
        assert!(f > 0 && c > 0 && s > 0);
        for g in &vb.geometry {
            // Staircased geometry is strictly 0/1 and carries no wall.
            assert!(g.volume_fraction == 0.0 || g.volume_fraction == 1.0);
            assert!(!g.cut_face.has_wall());
        }
    }
}

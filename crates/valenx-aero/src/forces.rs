//! Surface-force integration and the aerodynamic coefficients.
//!
//! Once the flow field has converged, the engineering answer a wind-
//! tunnel run exists to produce is a *number*: the drag and lift
//! coefficients. This module integrates the pressure and the viscous
//! wall shear over the immersed body to get the total aerodynamic
//! force and moment, then non-dimensionalises them.
//!
//! # How the force is integrated
//!
//! The immersed body is a set of cut cells, each with one or more
//! faces exposed to the fluid. For each exposed face the module reads
//! the adjacent fluid-cell pressure (the **pressure force**, normal to
//! the face) and the velocity gradient at the wall (the **viscous
//! force**, tangential — the wall shear stress `τ_w = μ_eff·∂u_t/∂n`).
//! Summed over every exposed face with the face's outward normal and
//! area, this gives the total force the air exerts on the body and the
//! moment about a reference point.
//!
//! # The coefficients
//!
//! With the force `F`, the dynamic pressure `q∞ = ½ρU∞²` and the
//! reference area `A`:
//!
//! ```text
//!   C_d = F_drag  / (q∞·A)        drag    — along the wind
//!   C_l = F_lift  / (q∞·A)        lift    — perpendicular, "up"
//!   C_s = F_side  / (q∞·A)        side    — perpendicular, lateral
//!   C_m = M       / (q∞·A·L)      moment  — about the reference point
//! ```
//!
//! # Cut-cell vs. staircased integration
//!
//! Under the default **cut-cell** method ([`crate::cutcell`]) the force
//! is integrated over the **true clipped wall faces** — each cut cell
//! carries one polygon ([`crate::cutcell::CutFace`]) with an exact area
//! and outward normal, the actual body surface inside that cell, not
//! the axis-aligned voxel faces. The pressure acts on that polygon
//! along its true normal and the wall shear acts in its true tangent
//! plane, so the integrated `Cd` / `Cl` reflect the real geometry —
//! this is where the cut-cell accuracy gain shows up. Under the legacy
//! [`crate::cutcell::WallMethod::Staircase`] method the integration
//! falls back to the staircased voxel faces.
//!
//! # Honest scope
//!
//! A real force integration. Under the cut-cell method the integrated
//! area is the true clipped surface (it tiles the body, converging with
//! the mesh); the pressure force is exact on the cut face and the
//! viscous force uses a one-sided wall-gradient estimate normal to the
//! true wall. The remaining gap to engineering tolerance is the lack of
//! a body-fitted near-wall prism layer — the boundary-layer profile is
//! still resolved on the uniform Cartesian grid.

use nalgebra::Vector3;

use crate::cutcell::WallMethod;
use crate::domain::WindTunnel;
use crate::immersed::CellTag;
use crate::solver::FlowField;
use crate::turbulence::TurbulenceModel;

/// The six exposed-face directions of a cut cell.
const FACE_DIRS: [(i32, i32, i32); 6] = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
];

/// The integrated aerodynamic force and moment on the body, plus the
/// drag breakdown.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AeroForces {
    /// Total force on the body in world axes (N).
    pub force: Vector3<f64>,
    /// Total moment about the reference point in world axes (N·m).
    pub moment: Vector3<f64>,
    /// The pressure (form / normal) part of the force (N).
    pub pressure_force: Vector3<f64>,
    /// The viscous (skin-friction / tangential) part of the force (N).
    pub viscous_force: Vector3<f64>,
}

impl AeroForces {
    /// The component of the force along a given (unit) direction.
    pub fn along(&self, dir: Vector3<f64>) -> f64 {
        self.force.dot(&dir)
    }
}

/// The non-dimensional aerodynamic coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AeroCoefficients {
    /// Drag coefficient `C_d` — force component along the wind.
    pub cd: f64,
    /// Lift coefficient `C_l` — force component "up" (perpendicular to
    /// the wind, in the vertical plane).
    pub cl: f64,
    /// Side-force coefficient `C_s` — the remaining lateral component.
    pub cs: f64,
    /// The pressure-drag part of `C_d`.
    pub cd_pressure: f64,
    /// The friction-drag part of `C_d`.
    pub cd_friction: f64,
    /// Roll-moment coefficient `C_mx` (about the wind axis).
    pub cmx: f64,
    /// Pitch-moment coefficient `C_my`.
    pub cmy: f64,
    /// Yaw-moment coefficient `C_mz`.
    pub cmz: f64,
}

/// An orthonormal wind frame — the drag / lift / side directions.
#[derive(Clone, Copy, Debug)]
pub struct WindFrame {
    /// Unit drag direction (along the free-stream).
    pub drag: Vector3<f64>,
    /// Unit lift direction ("up", perpendicular to drag).
    pub lift: Vector3<f64>,
    /// Unit side direction (lateral, perpendicular to drag and lift).
    pub side: Vector3<f64>,
}

impl WindFrame {
    /// Build the wind frame from a free-stream direction. The lift
    /// direction is chosen to lie in the vertical plane containing the
    /// drag direction (so for level flight lift is "up"); if the wind
    /// is exactly vertical, an arbitrary perpendicular is picked.
    pub fn from_direction(dir: Vector3<f64>) -> WindFrame {
        let drag = dir.try_normalize(1e-12).unwrap_or_else(Vector3::x);
        // The lift direction: project +z out of the drag direction.
        let up = Vector3::new(0.0, 0.0, 1.0);
        let lift = match (up - up.dot(&drag) * drag).try_normalize(1e-9) {
            Some(l) => l,
            None => {
                // Wind is vertical — pick +x's perpendicular.
                let alt = Vector3::new(1.0, 0.0, 0.0);
                (alt - alt.dot(&drag) * drag).normalize()
            }
        };
        let side = drag.cross(&lift);
        WindFrame { drag, lift, side }
    }
}

/// The wall-normal sampling distance from a cut cell's centre to its
/// clipped wall face — the boundary-layer height the near-wall model
/// samples for the surface force / `Cf` / `y+`.
///
/// The signed distance of the cell centre off the cut-face plane along
/// the outward normal is the true sampling height; it is floored at
/// 10 % of a cell and capped at one cell so a centroid that almost
/// touches the wall cannot produce a singular shear, with `half_cell`
/// the degenerate fallback.
fn wall_normal_distance(
    grid: &crate::grid::Grid3,
    cut: &crate::cutcell::CutFace,
    i: usize,
    j: usize,
    k: usize,
    half_cell: f64,
) -> f64 {
    let (cx, cy, cz) = grid.cell_centre(i, j, k);
    let centre = Vector3::new(cx, cy, cz);
    let d = (centre - cut.centroid).dot(&cut.normal).abs();
    let cell = 2.0 * half_cell;
    if d.is_finite() && d > 0.0 {
        d.clamp(0.1 * cell, cell)
    } else {
        half_cell.max(1e-9)
    }
}

/// Integrate the aerodynamic force and moment on the immersed body.
///
/// `tunnel` supplies the voxelized body and the air properties;
/// `flow` is the converged flow field; `moment_ref` is the world-space
/// point the moment is taken about (typically the body's bounding-box
/// centre).
///
/// Under the cut-cell method the integral runs over the true clipped
/// wall faces; under the staircased method it runs over the exposed
/// voxel faces.
pub fn integrate_forces(
    tunnel: &WindTunnel,
    flow: &FlowField,
    moment_ref: Vector3<f64>,
) -> AeroForces {
    integrate_forces_with(tunnel, flow, moment_ref, true)
}

/// [`integrate_forces`] with an explicit near-wall-model toggle.
///
/// With `wall_model` set (what [`integrate_forces`] uses) the viscous
/// surface force is the law-of-the-wall wall shear ([`crate::wallmodel`])
/// — the turbulent boundary-layer profile reconstructed from the first
/// cell. With it clear the viscous force is the legacy crude linear
/// gradient `τ_w = μ_eff·u_t/y`, kept so a caller can measure the
/// accuracy delta of the near-wall model directly. The toggle only
/// affects the cut-cell path; the staircased legacy path is unchanged.
pub fn integrate_forces_with(
    tunnel: &WindTunnel,
    flow: &FlowField,
    moment_ref: Vector3<f64>,
    wall_model: bool,
) -> AeroForces {
    if tunnel.body.method == WallMethod::CutCell {
        integrate_forces_cutcell(tunnel, flow, moment_ref, wall_model)
    } else {
        integrate_forces_staircase(tunnel, flow, moment_ref)
    }
}

/// Integrate the force / moment over the **true clipped cut faces** —
/// the cut-cell path. Each cut cell carries one wall polygon with an
/// exact area, outward normal and centroid; the pressure acts along the
/// true normal and the wall shear acts in the true tangent plane.
fn integrate_forces_cutcell(
    tunnel: &WindTunnel,
    flow: &FlowField,
    moment_ref: Vector3<f64>,
    wall_model: bool,
) -> AeroForces {
    let grid = tunnel.grid;
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let mu_lam = tunnel.wind.air.dynamic_viscosity;
    let rho = tunnel.wind.air.density;
    let nu = (mu_lam / rho).max(1e-12);
    let body = &tunnel.body;
    // The fallback half-cell scale for a degenerate wall-normal
    // distance; the true cut-face centroid distance is preferred.
    let half_cell = 0.5 * (dx * dx + dy * dy + dz * dz).sqrt() / 3.0_f64.sqrt();

    let mut pressure_force = Vector3::zeros();
    let mut viscous_force = Vector3::zeros();
    let mut moment = Vector3::zeros();

    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if body.tag(i, j, k) != CellTag::Cut {
                    continue;
                }
                let cut = body.cut_face(i, j, k);
                if !cut.has_wall() {
                    continue;
                }
                // `normal` is the outward wall normal (solid → fluid);
                // `area` is the true clipped wall-patch area.
                let normal = cut.normal;
                let area = cut.area;

                // --- pressure force ---
                // The cut cell's pressure acts on the wall patch. The
                // force on the body is −p·A·n (pressure pushes against
                // the outward-facing surface).
                let p = flow.pressure.at(i, j, k);
                let fp = -p * area * normal;
                pressure_force += fp;

                // --- viscous force ---
                // With the near-wall model the wall shear is recovered
                // from the cut cell's tangential speed and its true
                // wall-normal distance via the law of the wall, τ_w =
                // ρ·u_τ² — the turbulent boundary-layer profile, not a
                // crude linear gradient μ·u/d across the first cell.
                let vel = Vector3::new(
                    flow.u_at_cell(i, j, k),
                    flow.v_at_cell(i, j, k),
                    flow.w_at_cell(i, j, k),
                );
                let v_tan = vel - vel.dot(&normal) * normal;
                let u_t = v_tan.norm();
                let y = wall_normal_distance(&grid, cut, i, j, k, half_cell);
                let tau_mag = if wall_model {
                    crate::wallmodel::wall_shear_stress(rho, u_t, y, nu)
                } else {
                    // Legacy crude linear gradient — μ_eff·u_t/y.
                    let mu_eff = mu_lam + flow.turbulence.mu_t.at(i, j, k);
                    mu_eff * u_t / y.max(1e-9)
                };
                // The shear acts along the tangential-velocity direction
                // (opposing the body's relative motion through the air).
                let tau = if u_t > 1e-12 {
                    tau_mag * v_tan / u_t
                } else {
                    Vector3::zeros()
                };
                let fv = tau * area;
                viscous_force += fv;

                // --- moment about the reference point ---
                // The true wall-patch centroid is the force application
                // point.
                let arm = cut.centroid - moment_ref;
                moment += arm.cross(&(fp + fv));
            }
        }
    }

    AeroForces {
        force: pressure_force + viscous_force,
        moment,
        pressure_force,
        viscous_force,
    }
}

/// Integrate the force / moment over the staircased voxel faces — the
/// legacy path, kept for [`WallMethod::Staircase`].
fn integrate_forces_staircase(
    tunnel: &WindTunnel,
    flow: &FlowField,
    moment_ref: Vector3<f64>,
) -> AeroForces {
    let grid = tunnel.grid;
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let mu_lam = tunnel.wind.air.dynamic_viscosity;
    let body = &tunnel.body;

    let mut pressure_force = Vector3::zeros();
    let mut viscous_force = Vector3::zeros();
    let mut moment = Vector3::zeros();

    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if body.tag(i, j, k) != CellTag::Cut {
                    continue;
                }
                // A cut (fluid) cell — examine its six faces. A face
                // toward a solid cell is an exposed body surface.
                for &(di, dj, dk) in &FACE_DIRS {
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
                    let (ni, nj, nk) = (ni as usize, nj as usize, nk as usize);
                    if body.tag(ni, nj, nk) != CellTag::Solid {
                        continue;
                    }
                    // This face is body surface. Its outward normal
                    // (pointing into the fluid) is −(di,dj,dk).
                    let normal = Vector3::new(-di as f64, -dj as f64, -dk as f64);
                    // The face area depends on which axis it spans.
                    let area = if di != 0 {
                        dy * dz
                    } else if dj != 0 {
                        dx * dz
                    } else {
                        dx * dy
                    };
                    // The face centre — for the moment arm.
                    let (cx, cy, cz) = grid.cell_centre(i, j, k);
                    let fc = Vector3::new(
                        cx + 0.5 * di as f64 * dx,
                        cy + 0.5 * dj as f64 * dy,
                        cz + 0.5 * dk as f64 * dz,
                    );

                    // --- pressure force ---
                    let p = flow.pressure.at(i, j, k);
                    let fp = -p * area * normal;
                    pressure_force += fp;

                    // --- viscous force ---
                    let mu_eff = mu_lam + flow.turbulence.mu_t.at(i, j, k);
                    let wall_n = match (di.abs(), dj.abs(), dk.abs()) {
                        (1, 0, 0) => dx,
                        (0, 1, 0) => dy,
                        _ => dz,
                    };
                    let half = 0.5 * wall_n;
                    let vel = Vector3::new(
                        flow.u_at_cell(i, j, k),
                        flow.v_at_cell(i, j, k),
                        flow.w_at_cell(i, j, k),
                    );
                    let v_tan = vel - vel.dot(&normal) * normal;
                    let tau = mu_eff * v_tan / half.max(1e-9);
                    let fv = tau * area;
                    viscous_force += fv;

                    // --- moment ---
                    let arm = fc - moment_ref;
                    moment += arm.cross(&(fp + fv));
                }
            }
        }
    }

    AeroForces {
        force: pressure_force + viscous_force,
        moment,
        pressure_force,
        viscous_force,
    }
}

/// Non-dimensionalise the integrated forces into the aerodynamic
/// coefficients.
///
/// `tunnel` supplies the reference area, reference length and the wind;
/// `forces` is the integrated force / moment. The coefficients are
/// reported in the wind frame ([`WindFrame`]).
pub fn coefficients(tunnel: &WindTunnel, forces: &AeroForces) -> AeroCoefficients {
    let q = tunnel.dynamic_pressure();
    let a = tunnel.reference_area.max(1e-12);
    let l = tunnel.reference_length.max(1e-12);
    let denom = (q * a).max(1e-30);
    let denom_m = (q * a * l).max(1e-30);

    let frame = WindFrame::from_direction(tunnel.wind.direction());

    let cd = forces.force.dot(&frame.drag) / denom;
    let cl = forces.force.dot(&frame.lift) / denom;
    let cs = forces.force.dot(&frame.side) / denom;
    let cd_pressure = forces.pressure_force.dot(&frame.drag) / denom;
    let cd_friction = forces.viscous_force.dot(&frame.drag) / denom;

    AeroCoefficients {
        cd,
        cl,
        cs,
        cd_pressure,
        cd_friction,
        cmx: forces.moment.x / denom_m,
        cmy: forces.moment.y / denom_m,
        cmz: forces.moment.z / denom_m,
    }
}

/// One sampled point of the surface pressure distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePoint {
    /// World-space position of the body-surface face centre.
    pub position: Vector3<f64>,
    /// Outward surface normal at the point.
    pub normal: Vector3<f64>,
    /// Static pressure on the surface (Pa).
    pub pressure: f64,
    /// The pressure coefficient `Cp = (p − p∞) / q∞`.
    pub cp: f64,
    /// The skin-friction coefficient `Cf = τ_w / q∞`.
    pub cf: f64,
    /// The dimensionless wall distance `y+` of the first fluid cell.
    pub y_plus: f64,
}

/// Sample the surface pressure-coefficient, skin-friction and `y+`
/// fields over the immersed body.
///
/// Returns one [`SurfacePoint`] per body-surface patch. Under the
/// cut-cell method that is one point per true clipped wall face; under
/// the staircased method it is one per exposed voxel face. This is the
/// data a pressure-tap plot, a Cp colour map, a skin-friction streak
/// plot, or a `y+` mesh-quality check is built from.
pub fn surface_field(tunnel: &WindTunnel, flow: &FlowField) -> Vec<SurfacePoint> {
    if tunnel.body.method == WallMethod::CutCell {
        surface_field_cutcell(tunnel, flow)
    } else {
        surface_field_staircase(tunnel, flow)
    }
}

/// Sample the surface field over the true clipped cut faces.
fn surface_field_cutcell(tunnel: &WindTunnel, flow: &FlowField) -> Vec<SurfacePoint> {
    let grid = tunnel.grid;
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let rho = tunnel.wind.air.density;
    let mu_lam = tunnel.wind.air.dynamic_viscosity;
    let nu = (mu_lam / rho).max(1e-12);
    let q = tunnel.dynamic_pressure().max(1e-12);
    let body = &tunnel.body;
    let half_cell = 0.5 * (dx * dx + dy * dy + dz * dz).sqrt() / 3.0_f64.sqrt();
    let mut points = Vec::new();

    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if body.tag(i, j, k) != CellTag::Cut {
                    continue;
                }
                let cut = body.cut_face(i, j, k);
                if !cut.has_wall() {
                    continue;
                }
                let normal = cut.normal;
                let p = flow.pressure.at(i, j, k);
                let cp = p / q;

                // Skin friction and y+ from the near-wall model: the
                // friction velocity is recovered from the cut cell's
                // tangential speed and its true wall-normal distance via
                // the law of the wall, giving τ_w = ρ·u_τ², Cf = τ_w/q
                // and y+ = y·u_τ/ν — the turbulent-profile values, not a
                // linear-gradient guess.
                let vel = Vector3::new(
                    flow.u_at_cell(i, j, k),
                    flow.v_at_cell(i, j, k),
                    flow.w_at_cell(i, j, k),
                );
                let v_tan = vel - vel.dot(&normal) * normal;
                let u_t = v_tan.norm();
                let y = wall_normal_distance(&grid, cut, i, j, k, half_cell);
                let u_tau = crate::wallmodel::friction_velocity(u_t, y, nu);
                let tau_w = rho * u_tau * u_tau;
                let cf = tau_w / q;
                let y_plus = y * u_tau / nu;

                points.push(SurfacePoint {
                    position: cut.centroid,
                    normal,
                    pressure: p,
                    cp,
                    cf,
                    y_plus,
                });
            }
        }
    }
    points
}

/// Sample the surface field over the staircased voxel faces (legacy).
fn surface_field_staircase(tunnel: &WindTunnel, flow: &FlowField) -> Vec<SurfacePoint> {
    let grid = tunnel.grid;
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let rho = tunnel.wind.air.density;
    let mu_lam = tunnel.wind.air.dynamic_viscosity;
    let q = tunnel.dynamic_pressure().max(1e-12);
    let body = &tunnel.body;
    let mut points = Vec::new();

    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if body.tag(i, j, k) != CellTag::Cut {
                    continue;
                }
                for &(di, dj, dk) in &FACE_DIRS {
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
                    if body.tag(ni as usize, nj as usize, nk as usize) != CellTag::Solid {
                        continue;
                    }
                    let normal = Vector3::new(-di as f64, -dj as f64, -dk as f64);
                    let (cx, cy, cz) = grid.cell_centre(i, j, k);
                    let pos = Vector3::new(
                        cx + 0.5 * di as f64 * dx,
                        cy + 0.5 * dj as f64 * dy,
                        cz + 0.5 * dk as f64 * dz,
                    );
                    let p = flow.pressure.at(i, j, k);
                    // Cp = (p − p∞)/q∞; the gauge pressure has p∞ ≈ 0
                    // because the outlet is pinned to zero.
                    let cp = p / q;

                    // Wall shear and Cf.
                    let mu_eff = mu_lam + flow.turbulence.mu_t.at(i, j, k);
                    let wall_n = match (di.abs(), dj.abs(), dk.abs()) {
                        (1, 0, 0) => dx,
                        (0, 1, 0) => dy,
                        _ => dz,
                    };
                    let half = (0.5 * wall_n).max(1e-9);
                    let vel = Vector3::new(
                        flow.u_at_cell(i, j, k),
                        flow.v_at_cell(i, j, k),
                        flow.w_at_cell(i, j, k),
                    );
                    let v_tan = vel - vel.dot(&normal) * normal;
                    let tau_w = mu_eff * v_tan.norm() / half;
                    let cf = tau_w / q;

                    // y+ = (ρ·u_τ·y) / μ, u_τ = √(τ_w/ρ).
                    let u_tau = (tau_w / rho).sqrt();
                    let y_plus = rho * u_tau * half / mu_lam.max(1e-30);

                    points.push(SurfacePoint {
                        position: pos,
                        normal,
                        pressure: p,
                        cp,
                        cf,
                        y_plus,
                    });
                }
            }
        }
    }
    points
}

/// Summary statistics of the surface field — handy for an at-a-glance
/// report without walking every [`SurfacePoint`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceStats {
    /// Number of body-surface faces sampled.
    pub face_count: usize,
    /// Minimum surface pressure coefficient (the suction peak).
    pub cp_min: f64,
    /// Maximum surface pressure coefficient (stagnation ≈ +1).
    pub cp_max: f64,
    /// Mean dimensionless wall distance `y+`.
    pub y_plus_mean: f64,
    /// Maximum `y+`.
    pub y_plus_max: f64,
}

/// Reduce a surface field to summary statistics.
pub fn surface_stats(points: &[SurfacePoint]) -> SurfaceStats {
    if points.is_empty() {
        return SurfaceStats {
            face_count: 0,
            cp_min: 0.0,
            cp_max: 0.0,
            y_plus_mean: 0.0,
            y_plus_max: 0.0,
        };
    }
    let mut cp_min = f64::INFINITY;
    let mut cp_max = f64::NEG_INFINITY;
    let mut yp_sum = 0.0;
    let mut yp_max = 0.0f64;
    for p in points {
        cp_min = cp_min.min(p.cp);
        cp_max = cp_max.max(p.cp);
        yp_sum += p.y_plus;
        yp_max = yp_max.max(p.y_plus);
    }
    SurfaceStats {
        face_count: points.len(),
        cp_min,
        cp_max,
        y_plus_mean: yp_sum / points.len() as f64,
        y_plus_max: yp_max,
    }
}

/// A note on the turbulence model the forces were computed with — used
/// by the report to qualify the accuracy.
pub fn turbulence_note(model: TurbulenceModel) -> &'static str {
    match model {
        TurbulenceModel::Laminar => "laminar — coefficients valid only at low Reynolds number",
        TurbulenceModel::KEpsilon => {
            "k-epsilon — robust, weaker in strong adverse-pressure-gradient separation"
        }
        TurbulenceModel::KOmegaSST => {
            "k-omega SST — the external-aero standard, better separation behaviour"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::{BoundaryConditions, TunnelSizing, WindTunnel};
    use crate::geometry::box_body;
    use crate::solver::{solve_steady, BodyMotion, SolverControls};
    use crate::wind::Wind;

    /// Build a wind tunnel on a deliberately coarse grid — used by the
    /// force-integration tests that overwrite the flow field anyway
    /// (they need real array shapes, not a converged solution).
    fn build_coarse_tunnel(body: &crate::geometry::TriMesh) -> WindTunnel {
        WindTunnel::build_with(
            body,
            Wind::straight(20.0).unwrap(),
            BoundaryConditions::external_aero(),
            TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..TunnelSizing::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn wind_frame_is_orthonormal() {
        let f = WindFrame::from_direction(Vector3::new(1.0, 0.3, 0.1));
        // All three axes unit length.
        assert!((f.drag.norm() - 1.0).abs() < 1e-12);
        assert!((f.lift.norm() - 1.0).abs() < 1e-12);
        assert!((f.side.norm() - 1.0).abs() < 1e-12);
        // Mutually perpendicular.
        assert!(f.drag.dot(&f.lift).abs() < 1e-12);
        assert!(f.drag.dot(&f.side).abs() < 1e-12);
        assert!(f.lift.dot(&f.side).abs() < 1e-12);
    }

    #[test]
    fn wind_frame_straight_wind_has_standard_axes() {
        // Straight +x wind: drag = +x, lift = +z, side = −y (right-
        // handed drag×lift).
        let f = WindFrame::from_direction(Vector3::new(1.0, 0.0, 0.0));
        assert!((f.drag - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-9);
        assert!((f.lift - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-9);
    }

    #[test]
    fn uniform_pressure_field_integrates_to_zero_resultant() {
        // The analytic check the task asks for: a closed body in a
        // *uniform* pressure field feels zero net pressure force —
        // the integral of a constant pressure over a closed surface
        // vanishes. Build a box tunnel, overwrite the pressure field
        // with a constant, and confirm the integrated pressure force
        // is ~0.
        let body = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        // The flow field is fully overwritten below — only its array
        // *shapes* are needed, so a coarse grid + a single iteration is
        // all this test requires.
        let tunnel = build_coarse_tunnel(&body);
        let controls = SolverControls {
            max_iterations: 1,
            turbulence: TurbulenceModel::Laminar,
            ..SolverControls::default()
        };
        let mut flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        // Overwrite with a uniform pressure and a zero velocity so the
        // viscous force vanishes too.
        flow.pressure.fill(101325.0);
        flow.u.fill(0.0);
        flow.v.fill(0.0);
        flow.w.fill(0.0);
        let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
        // A uniform pressure over a closed (voxel) surface → ~0 net
        // force. The staircased surface is closed, so this is exact.
        assert!(
            forces.pressure_force.norm() < 1e-6 * 101325.0,
            "uniform pressure should give zero resultant, got {:?}",
            forces.pressure_force
        );
        assert!(forces.viscous_force.norm() < 1e-9);
    }

    #[test]
    fn one_sided_pressure_gives_a_force_toward_low_pressure() {
        // A physical directional check: if the pressure on the +x
        // side of a box is higher than on the −x side, the net force
        // pushes the body in −x. We fake such a field and confirm the
        // sign.
        let body = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        // The flow is fully overwritten below — a coarse grid + a single
        // iteration suffices to allocate the field arrays.
        let tunnel = build_coarse_tunnel(&body);
        let controls = SolverControls {
            max_iterations: 1,
            turbulence: TurbulenceModel::Laminar,
            ..SolverControls::default()
        };
        let mut flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        flow.u.fill(0.0);
        flow.v.fill(0.0);
        flow.w.fill(0.0);
        // Set a pressure that rises with x — the upstream face of the
        // body sees higher pressure.
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    let (cx, _, _) = tunnel.grid.cell_centre(i, j, k);
                    flow.pressure.set(i, j, k, 100.0 * cx);
                }
            }
        }
        let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
        // Higher pressure on the upstream (−x... actually +x grows) —
        // pressure increasing with x pushes the body toward −x.
        assert!(
            forces.pressure_force.x < 0.0,
            "force should point toward low pressure (−x), got {}",
            forces.pressure_force.x
        );
    }

    #[test]
    fn box_drag_coefficient_is_in_the_bluff_body_ballpark() {
        // A cube broadside to the flow has a textbook drag coefficient
        // of roughly Cd ≈ 1.0–1.3 in the turbulent regime. The
        // immersed-boundary v1 won't nail it, but it must land in a
        // physically plausible bluff-body band, not at 0 or 100.
        //
        // A deliberately coarse grid keeps the steady SIMPLE solve fast
        // enough for the test suite — the bluff-body Cd band is wide
        // and the qualitative result (a converged flow, a plausible
        // drag) does not need a fine mesh.
        let body = box_body(Vector3::new(-0.5, -0.5, -0.5), Vector3::new(0.5, 0.5, 0.5));
        let tunnel = WindTunnel::build_with(
            &body,
            Wind::straight(20.0).unwrap(),
            BoundaryConditions::external_aero(),
            TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..TunnelSizing::default()
            },
        )
        .unwrap();
        let controls = SolverControls {
            max_iterations: 60,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
        let coeff = coefficients(&tunnel, &forces);
        // A bluff body: Cd is a positive O(1) number.
        assert!(
            coeff.cd > 0.1 && coeff.cd < 5.0,
            "cube Cd {} outside the plausible bluff-body band",
            coeff.cd
        );
        // Drag = pressure drag + friction drag.
        let sum = coeff.cd_pressure + coeff.cd_friction;
        assert!(
            (sum - coeff.cd).abs() < 1e-6,
            "drag breakdown must sum to total Cd"
        );
        // For a bluff body the pressure drag dominates the friction
        // drag.
        assert!(
            coeff.cd_pressure > coeff.cd_friction,
            "bluff-body drag should be pressure-dominated"
        );
    }

    #[test]
    fn surface_field_and_stats_are_consistent() {
        let body = box_body(Vector3::new(-0.5, -0.5, -0.5), Vector3::new(0.5, 0.5, 0.5));
        // A coarse grid keeps this real solve fast — the test asserts
        // surface-field *consistency* (cp_min ≤ cp_max, finite Cp/Cf/y+),
        // which holds at any resolution.
        let tunnel = WindTunnel::build_with(
            &body,
            Wind::straight(20.0).unwrap(),
            crate::domain::BoundaryConditions::external_aero(),
            crate::domain::TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..crate::domain::TunnelSizing::default()
            },
        )
        .unwrap();
        let controls = SolverControls {
            max_iterations: 40,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        let surf = surface_field(&tunnel, &flow);
        assert!(!surf.is_empty(), "a body should have surface faces");
        let stats = surface_stats(&surf);
        assert_eq!(stats.face_count, surf.len());
        assert!(stats.cp_min <= stats.cp_max);
        // y+ must be finite and non-negative.
        assert!(stats.y_plus_max.is_finite() && stats.y_plus_max >= 0.0);
        assert!(stats.y_plus_mean >= 0.0);
        // Every surface point's Cp / Cf must be finite.
        assert!(surf.iter().all(|p| p.cp.is_finite() && p.cf.is_finite()));
    }

    #[test]
    fn empty_surface_field_gives_zero_stats() {
        let stats = surface_stats(&[]);
        assert_eq!(stats.face_count, 0);
        assert_eq!(stats.cp_min, 0.0);
    }

    /// Solve a body at 20 m/s on a coarse grid with a chosen wall
    /// method and return the drag coefficient — the shared helper for
    /// the cut-cell-vs-staircased accuracy tests.
    fn drag_with_method(
        body: &crate::geometry::TriMesh,
        method: crate::cutcell::WallMethod,
        cells_across: usize,
    ) -> f64 {
        use crate::immersed::voxelize_with;
        let mut tunnel = WindTunnel::build_with(
            body,
            Wind::straight(20.0).unwrap(),
            BoundaryConditions::external_aero(),
            TunnelSizing {
                cells_across_body: cells_across,
                max_cells: 200_000,
                ..TunnelSizing::default()
            },
        )
        .unwrap();
        // Re-voxelize with the requested method (the tunnel builder
        // defaults to cut-cell).
        tunnel.body = voxelize_with(&tunnel.grid, body, method);
        let controls = SolverControls {
            max_iterations: 120,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
        coefficients(&tunnel, &forces).cd
    }

    #[test]
    fn cutcell_cube_drag_is_at_least_as_accurate_as_staircased() {
        // A cube broadside to the flow has a textbook drag coefficient
        // of about Cd ≈ 1.05 in the turbulent regime. The staircased
        // immersed boundary over-predicts it (the voxel staircase sheds
        // a blockier wake); the cut-cell method, integrating pressure
        // over the true clipped faces, must land *at least as close* to
        // the reference — the wall-accuracy upgrade this is.
        use crate::cutcell::WallMethod;
        let cube = box_body(Vector3::new(-0.5, -0.5, -0.5), Vector3::new(0.5, 0.5, 0.5));
        let reference = 1.05;
        let cd_stair = drag_with_method(&cube, WallMethod::Staircase, 8);
        let cd_cut = drag_with_method(&cube, WallMethod::CutCell, 8);
        // Both must be plausible positive bluff-body numbers.
        assert!(
            cd_stair > 0.3 && cd_stair < 3.0,
            "staircased cube Cd {cd_stair} implausible"
        );
        assert!(
            cd_cut > 0.3 && cd_cut < 3.0,
            "cut-cell cube Cd {cd_cut} implausible"
        );
        // The cut-cell error must not exceed the staircased error — the
        // upgrade is required to be no worse, and is observed closer
        // (cube: staircased ~1.36, cut-cell ~1.26 vs the ~1.05 ref).
        let err_stair = (cd_stair - reference).abs();
        let err_cut = (cd_cut - reference).abs();
        assert!(
            err_cut <= err_stair + 1e-6,
            "cut-cell cube Cd {cd_cut} (err {err_cut}) is less accurate \
             than staircased {cd_stair} (err {err_stair}) vs ref {reference}"
        );
    }

    #[test]
    fn cutcell_sphere_drag_is_at_least_as_accurate_as_staircased() {
        // A sphere at this Reynolds number has a textbook Cd in the
        // 0.4–0.5 band. The staircased immersed boundary turns the
        // smooth sphere into a blocky shape and over-predicts the drag;
        // the cut-cell method recovers the true curved surface and must
        // land at least as close to the reference.
        use crate::cutcell::WallMethod;
        use crate::geometry::sphere_body;
        let sphere = sphere_body(Vector3::zeros(), 0.5, 32, 64);
        let reference = 0.5;
        let cd_stair = drag_with_method(&sphere, WallMethod::Staircase, 8);
        let cd_cut = drag_with_method(&sphere, WallMethod::CutCell, 8);
        assert!(
            cd_stair > 0.1 && cd_stair < 3.0,
            "staircased sphere Cd {cd_stair} implausible"
        );
        assert!(
            cd_cut > 0.1 && cd_cut < 3.0,
            "cut-cell sphere Cd {cd_cut} implausible"
        );
        let err_stair = (cd_stair - reference).abs();
        let err_cut = (cd_cut - reference).abs();
        assert!(
            err_cut <= err_stair + 1e-6,
            "cut-cell sphere Cd {cd_cut} (err {err_cut}) is less accurate \
             than staircased {cd_stair} (err {err_stair}) vs ref {reference}"
        );
    }

    #[test]
    fn flat_plate_aligned_with_the_flow_has_low_drag() {
        // A thin flat plate edge-on to the wind (its broad faces
        // parallel to the flow) presents almost no frontal area — its
        // drag is friction-dominated and small. The cut-cell method
        // must reproduce that: Cd on the plate's *planform* area is a
        // small O(0.01–0.1) number, far below a bluff body's O(1).
        let plate = box_body(
            Vector3::new(-0.5, -0.5, -0.01),
            Vector3::new(0.5, 0.5, 0.01),
        );
        // Reference area for the coefficient: the planform (1×1).
        let mut tunnel = WindTunnel::build_with(
            &plate,
            Wind::straight(20.0).unwrap(),
            BoundaryConditions::external_aero(),
            TunnelSizing {
                cells_across_body: 6,
                max_cells: 200_000,
                ..TunnelSizing::default()
            },
        )
        .unwrap();
        // Normalise on the planform area so the number is a skin-
        // friction-scale coefficient.
        tunnel.reference_area = 1.0;
        let controls = SolverControls {
            max_iterations: 80,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
        let coeff = coefficients(&tunnel, &forces);
        // Edge-on plate: drag is small and positive (friction-led).
        assert!(
            coeff.cd > 0.0 && coeff.cd < 0.3,
            "flow-aligned flat plate Cd {} should be small",
            coeff.cd
        );
        // And much smaller than a broadside bluff body's ~1.
        assert!(coeff.cd < 0.5);
    }
}

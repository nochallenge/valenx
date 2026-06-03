//! Flow-field extraction and post-processing.
//!
//! The converged [`FlowField`] is a raw block of numbers; this module
//! turns it into the things an engineer actually looks at: a wake
//! survey behind the body, streamline traces through the flow, the
//! vorticity and Q-criterion fields that reveal vortex cores, and
//! axis-aligned cut-plane slices for visualization.

use nalgebra::Vector3;

use crate::domain::WindTunnel;
use crate::grid::{Field3, Grid3};
use crate::solver::FlowField;

/// A wake-survey line — the velocity deficit sampled along a line
/// behind the body.
#[derive(Clone, Debug)]
pub struct WakeSurvey {
    /// World-space sample positions along the survey line.
    pub positions: Vec<Vector3<f64>>,
    /// Velocity-magnitude at each sample (m·s⁻¹).
    pub speed: Vec<f64>,
    /// The velocity deficit `1 − speed/U∞` at each sample — `0` in
    /// undisturbed flow, positive in the wake.
    pub deficit: Vec<f64>,
}

impl WakeSurvey {
    /// The peak velocity deficit along the survey — the depth of the
    /// wake.
    pub fn peak_deficit(&self) -> f64 {
        self.deficit.iter().cloned().fold(0.0, f64::max)
    }
}

/// Sample a wake survey on a line behind the body.
///
/// The line runs `length` metres along the wind direction starting a
/// distance `behind` past the body's downstream extent, sampled at
/// `samples` points. This is the classic wake-rake measurement —
/// integrating the deficit gives the profile drag.
pub fn wake_survey(
    tunnel: &WindTunnel,
    flow: &FlowField,
    behind: f64,
    length: f64,
    samples: usize,
) -> WakeSurvey {
    let samples = samples.max(2);
    let dir = tunnel.wind.direction();
    let u_inf = tunnel.wind.speed.max(1e-9);
    // Start point: the body's downstream extent plus the offset.
    let solid_max_x = body_downstream_extent(tunnel);
    let start = Vector3::new(
        solid_max_x + behind,
        tunnel.grid.y0 + 0.5 * tunnel.grid.ly,
        tunnel.grid.z0 + 0.5 * tunnel.grid.lz,
    );
    let mut positions = Vec::with_capacity(samples);
    let mut speed = Vec::with_capacity(samples);
    let mut deficit = Vec::with_capacity(samples);
    for s in 0..samples {
        let t = s as f64 / (samples - 1) as f64;
        let p = start + dir * (t * length);
        let sp = sample_speed(flow, &tunnel.grid, p);
        positions.push(p);
        speed.push(sp);
        deficit.push((1.0 - sp / u_inf).clamp(-1.0, 1.0));
    }
    WakeSurvey {
        positions,
        speed,
        deficit,
    }
}

/// The downstream `x`-extent of the body's solid cells.
fn body_downstream_extent(tunnel: &WindTunnel) -> f64 {
    let g = tunnel.grid;
    let mut max_x = g.x0;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                if tunnel.body.is_solid(i, j, k) {
                    let (cx, _, _) = g.cell_centre(i, j, k);
                    if cx > max_x {
                        max_x = cx;
                    }
                }
            }
        }
    }
    max_x
}

/// Trilinearly sample the cell-centred velocity magnitude at a
/// world-space point. Out-of-domain points return the free-stream-ish
/// nearest-cell value.
fn sample_speed(flow: &FlowField, grid: &Grid3, p: Vector3<f64>) -> f64 {
    let vel = sample_velocity(flow, grid, p);
    vel.norm()
}

/// Trilinearly sample the cell-centred velocity vector at a
/// world-space point.
pub fn sample_velocity(flow: &FlowField, grid: &Grid3, p: Vector3<f64>) -> Vector3<f64> {
    // Continuous cell-index coordinates (cell centres at +0.5).
    let fx = ((p.x - grid.x0) / grid.dx() - 0.5).clamp(0.0, (grid.nx - 1) as f64);
    let fy = ((p.y - grid.y0) / grid.dy() - 0.5).clamp(0.0, (grid.ny - 1) as f64);
    let fz = ((p.z - grid.z0) / grid.dz() - 0.5).clamp(0.0, (grid.nz - 1) as f64);
    let (i0, j0, k0) = (fx.floor() as usize, fy.floor() as usize, fz.floor() as usize);
    let i1 = (i0 + 1).min(grid.nx - 1);
    let j1 = (j0 + 1).min(grid.ny - 1);
    let k1 = (k0 + 1).min(grid.nz - 1);
    let (tx, ty, tz) = (fx - i0 as f64, fy - j0 as f64, fz - k0 as f64);

    let comp = |get: &dyn Fn(usize, usize, usize) -> f64| -> f64 {
        let c000 = get(i0, j0, k0);
        let c100 = get(i1, j0, k0);
        let c010 = get(i0, j1, k0);
        let c110 = get(i1, j1, k0);
        let c001 = get(i0, j0, k1);
        let c101 = get(i1, j0, k1);
        let c011 = get(i0, j1, k1);
        let c111 = get(i1, j1, k1);
        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;
        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tz) + c1 * tz
    };
    Vector3::new(
        comp(&|i, j, k| flow.u_at_cell(i, j, k)),
        comp(&|i, j, k| flow.v_at_cell(i, j, k)),
        comp(&|i, j, k| flow.w_at_cell(i, j, k)),
    )
}

/// A single streamline — a polyline traced along the velocity field.
#[derive(Clone, Debug)]
pub struct Streamline {
    /// The traced points, in order from the seed.
    pub points: Vec<Vector3<f64>>,
}

impl Streamline {
    /// The streamline's arc length.
    pub fn length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| (w[1] - w[0]).norm())
            .sum()
    }
}

/// Trace a streamline (particle trace) from a seed point.
///
/// Integrates `dx/dt = u(x)` with a fixed-step RK4 scheme for up to
/// `max_steps` steps of size `step` metres. Tracing stops if the
/// particle leaves the domain or enters a solid cell (it has hit the
/// body).
pub fn trace_streamline(
    tunnel: &WindTunnel,
    flow: &FlowField,
    seed: Vector3<f64>,
    step: f64,
    max_steps: usize,
) -> Streamline {
    let grid = &tunnel.grid;
    let mut points = Vec::with_capacity(max_steps + 1);
    let mut p = seed;
    points.push(p);
    // `step` is an ARC-LENGTH increment — each RK4 step advances the
    // trace `step` units in space, independent of the local flow speed.
    // (Integrating dx/dt = u directly with `step` as a time increment
    // would make the spatial spacing proportional to |u|, so a fast
    // uniform flow would be traced in only a handful of points.) The
    // velocity is therefore normalised to a unit direction field.
    let direction = |q: Vector3<f64>| -> Vector3<f64> {
        let v = sample_velocity(flow, grid, q);
        let n = v.norm();
        if n < 1e-9 {
            Vector3::zeros()
        } else {
            v / n
        }
    };
    for _ in 0..max_steps {
        // RK4 on the unit-direction field — fixed arc-length step.
        let k1 = direction(p);
        if k1.norm() < 1e-9 {
            break;
        }
        let k2 = direction(p + 0.5 * step * k1);
        let k3 = direction(p + 0.5 * step * k2);
        let k4 = direction(p + step * k3);
        p += (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        // Stop if out of the domain.
        if p.x < grid.x0
            || p.x > grid.x0 + grid.lx
            || p.y < grid.y0
            || p.y > grid.y0 + grid.ly
            || p.z < grid.z0
            || p.z > grid.z0 + grid.lz
        {
            break;
        }
        // Stop if inside the body.
        if point_in_solid(tunnel, p) {
            break;
        }
        points.push(p);
    }
    Streamline { points }
}

/// True if a world-space point falls in a solid cell.
fn point_in_solid(tunnel: &WindTunnel, p: Vector3<f64>) -> bool {
    let g = &tunnel.grid;
    let i = ((p.x - g.x0) / g.dx()).floor();
    let j = ((p.y - g.y0) / g.dy()).floor();
    let k = ((p.z - g.z0) / g.dz()).floor();
    if i < 0.0 || j < 0.0 || k < 0.0 {
        return false;
    }
    let (i, j, k) = (i as usize, j as usize, k as usize);
    if i >= g.nx || j >= g.ny || k >= g.nz {
        return false;
    }
    tunnel.body.is_solid(i, j, k)
}

/// The vorticity field — `ω = ∇ × u`, three cell-centred components.
#[derive(Clone, Debug)]
pub struct VorticityField {
    /// x-component of vorticity.
    pub omega_x: Field3,
    /// y-component of vorticity.
    pub omega_y: Field3,
    /// z-component of vorticity.
    pub omega_z: Field3,
    /// The vorticity magnitude `|ω|`.
    pub magnitude: Field3,
}

/// Compute the vorticity field `ω = ∇ × u` from the flow.
pub fn vorticity(flow: &FlowField) -> VorticityField {
    let g = &flow.grid;
    let (uc, vc, wc) = flow.cell_centred_velocity();
    let mut omega_x = g.scalar_field();
    let mut omega_y = g.scalar_field();
    let mut omega_z = g.scalar_field();
    let mut magnitude = g.scalar_field();
    let (dx, dy, dz) = (g.dx(), g.dy(), g.dz());

    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let cdx = |f: &Field3| {
                    if i > 0 && i + 1 < g.nx {
                        (f.at(i + 1, j, k) - f.at(i - 1, j, k)) / (2.0 * dx)
                    } else {
                        0.0
                    }
                };
                let cdy = |f: &Field3| {
                    if j > 0 && j + 1 < g.ny {
                        (f.at(i, j + 1, k) - f.at(i, j - 1, k)) / (2.0 * dy)
                    } else {
                        0.0
                    }
                };
                let cdz = |f: &Field3| {
                    if k > 0 && k + 1 < g.nz {
                        (f.at(i, j, k + 1) - f.at(i, j, k - 1)) / (2.0 * dz)
                    } else {
                        0.0
                    }
                };
                // ω = (∂w/∂y − ∂v/∂z, ∂u/∂z − ∂w/∂x, ∂v/∂x − ∂u/∂y).
                let wx = cdy(&wc) - cdz(&vc);
                let wy = cdz(&uc) - cdx(&wc);
                let wz = cdx(&vc) - cdy(&uc);
                omega_x.set(i, j, k, wx);
                omega_y.set(i, j, k, wy);
                omega_z.set(i, j, k, wz);
                magnitude.set(i, j, k, (wx * wx + wy * wy + wz * wz).sqrt());
            }
        }
    }
    VorticityField {
        omega_x,
        omega_y,
        omega_z,
        magnitude,
    }
}

/// Compute the Q-criterion field — `Q = ½(|Ω|² − |S|²)`, where `Ω` is
/// the rotation-rate tensor and `S` the strain-rate tensor.
///
/// A positive `Q` marks a region where rotation dominates strain — the
/// standard objective definition of a **vortex core**. Iso-surfaces of
/// `Q > 0` are the canonical vortex-visualization data for an
/// external-aero wake.
pub fn q_criterion(flow: &FlowField) -> Field3 {
    let g = &flow.grid;
    let (uc, vc, wc) = flow.cell_centred_velocity();
    let mut q = g.scalar_field();
    let (dx, dy, dz) = (g.dx(), g.dy(), g.dz());

    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let cdx = |f: &Field3| {
                    if i > 0 && i + 1 < g.nx {
                        (f.at(i + 1, j, k) - f.at(i - 1, j, k)) / (2.0 * dx)
                    } else {
                        0.0
                    }
                };
                let cdy = |f: &Field3| {
                    if j > 0 && j + 1 < g.ny {
                        (f.at(i, j + 1, k) - f.at(i, j - 1, k)) / (2.0 * dy)
                    } else {
                        0.0
                    }
                };
                let cdz = |f: &Field3| {
                    if k > 0 && k + 1 < g.nz {
                        (f.at(i, j, k + 1) - f.at(i, j, k - 1)) / (2.0 * dz)
                    } else {
                        0.0
                    }
                };
                let dudx = cdx(&uc);
                let dudy = cdy(&uc);
                let dudz = cdz(&uc);
                let dvdx = cdx(&vc);
                let dvdy = cdy(&vc);
                let dvdz = cdz(&vc);
                let dwdx = cdx(&wc);
                let dwdy = cdy(&wc);
                let dwdz = cdz(&wc);
                // Strain-rate tensor S and its squared norm.
                let s11 = dudx;
                let s22 = dvdy;
                let s33 = dwdz;
                let s12 = 0.5 * (dudy + dvdx);
                let s13 = 0.5 * (dudz + dwdx);
                let s23 = 0.5 * (dvdz + dwdy);
                let s_norm2 = s11 * s11
                    + s22 * s22
                    + s33 * s33
                    + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                // Rotation-rate tensor Ω and its squared norm.
                let o12 = 0.5 * (dudy - dvdx);
                let o13 = 0.5 * (dudz - dwdx);
                let o23 = 0.5 * (dvdz - dwdy);
                let o_norm2 = 2.0 * (o12 * o12 + o13 * o13 + o23 * o23);
                q.set(i, j, k, 0.5 * (o_norm2 - s_norm2));
            }
        }
    }
    q
}

/// The axis a cut plane is normal to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SliceAxis {
    /// A plane of constant `x` (a y-z slice).
    X,
    /// A plane of constant `y` (an x-z slice).
    Y,
    /// A plane of constant `z` (an x-y slice).
    Z,
}

/// A 2-D field slice through the 3-D domain — an axis-aligned cut plane
/// for visualization.
#[derive(Clone, Debug)]
pub struct FieldSlice {
    /// The axis the plane is normal to.
    pub axis: SliceAxis,
    /// The world-space coordinate of the plane along its normal axis.
    pub coordinate: f64,
    /// Width of the slice (cells along the first in-plane axis).
    pub width: usize,
    /// Height of the slice (cells along the second in-plane axis).
    pub height: usize,
    /// The sampled scalar values, row-major over `width × height`.
    pub values: Vec<f64>,
}

impl FieldSlice {
    /// Read the slice value at in-plane index `(a, b)`.
    pub fn at(&self, a: usize, b: usize) -> f64 {
        self.values[a + b * self.width]
    }

    /// The minimum and maximum value over the slice.
    pub fn range(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in &self.values {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        if self.values.is_empty() {
            (0.0, 0.0)
        } else {
            (lo, hi)
        }
    }
}

/// Cut a 2-D slice of a cell-centred field at a plane.
///
/// `field` must be a cell-centred field (`nx · ny · nz`) — the
/// pressure, a velocity component, the vorticity magnitude, etc.
/// `axis` and `coordinate` place the plane; the nearest grid layer is
/// taken.
pub fn slice_field(
    grid: &Grid3,
    field: &Field3,
    axis: SliceAxis,
    coordinate: f64,
) -> FieldSlice {
    match axis {
        SliceAxis::X => {
            let i = (((coordinate - grid.x0) / grid.dx()).floor() as i64)
                .clamp(0, grid.nx as i64 - 1) as usize;
            let mut values = vec![0.0; grid.ny * grid.nz];
            for k in 0..grid.nz {
                for j in 0..grid.ny {
                    values[j + k * grid.ny] = field.at(i, j, k);
                }
            }
            FieldSlice {
                axis,
                coordinate,
                width: grid.ny,
                height: grid.nz,
                values,
            }
        }
        SliceAxis::Y => {
            let j = (((coordinate - grid.y0) / grid.dy()).floor() as i64)
                .clamp(0, grid.ny as i64 - 1) as usize;
            let mut values = vec![0.0; grid.nx * grid.nz];
            for k in 0..grid.nz {
                for i in 0..grid.nx {
                    values[i + k * grid.nx] = field.at(i, j, k);
                }
            }
            FieldSlice {
                axis,
                coordinate,
                width: grid.nx,
                height: grid.nz,
                values,
            }
        }
        SliceAxis::Z => {
            let k = (((coordinate - grid.z0) / grid.dz()).floor() as i64)
                .clamp(0, grid.nz as i64 - 1) as usize;
            let mut values = vec![0.0; grid.nx * grid.ny];
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    values[i + j * grid.nx] = field.at(i, j, k);
                }
            }
            FieldSlice {
                axis,
                coordinate,
                width: grid.nx,
                height: grid.ny,
                values,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::WindTunnel;
    use crate::geometry::{box_body, sphere_body};
    use crate::solver::{solve_steady, BodyMotion, SolverControls};
    use crate::turbulence::TurbulenceModel;
    use crate::wind::Wind;

    fn quick_flow(body: &crate::geometry::TriMesh) -> (WindTunnel, FlowField) {
        // A deliberately coarse grid: the post-processing tests need a
        // *converged* flow field, not a fine one, and a coarse grid
        // keeps the steady SIMPLE solve fast enough for the suite.
        let tunnel = WindTunnel::build_with(
            body,
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
            max_iterations: 60,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        (tunnel, flow)
    }

    #[test]
    fn sample_velocity_recovers_a_uniform_field() {
        // A flow field set to a uniform velocity must sample back to
        // that velocity anywhere in the domain.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let (tunnel, mut flow) = quick_flow(&body);
        flow.u.fill(7.0);
        flow.v.fill(-2.0);
        flow.w.fill(3.0);
        let p = Vector3::new(
            tunnel.grid.x0 + 0.4 * tunnel.grid.lx,
            tunnel.grid.y0 + 0.6 * tunnel.grid.ly,
            tunnel.grid.z0 + 0.5 * tunnel.grid.lz,
        );
        let v = sample_velocity(&flow, &tunnel.grid, p);
        assert!((v - Vector3::new(7.0, -2.0, 3.0)).norm() < 1e-9);
    }

    #[test]
    fn wake_survey_finds_a_deficit_behind_a_body() {
        // Behind a bluff body the flow is slowed — the wake survey
        // must report a positive velocity deficit somewhere.
        let body = sphere_body(Vector3::zeros(), 0.5, 20, 40);
        let (tunnel, flow) = quick_flow(&body);
        let survey = wake_survey(&tunnel, &flow, 0.2, 2.0, 20);
        assert_eq!(survey.positions.len(), 20);
        assert!(
            survey.peak_deficit() > 0.0,
            "there should be a wake deficit behind the body"
        );
        // All deficits finite.
        assert!(survey.deficit.iter().all(|d| d.is_finite()));
    }

    #[test]
    fn streamline_through_uniform_flow_runs_straight() {
        // In a uniform +x flow a streamline is a straight line along x.
        let body = box_body(Vector3::zeros(), Vector3::new(0.3, 0.3, 0.3));
        let (tunnel, mut flow) = quick_flow(&body);
        flow.u.fill(10.0);
        flow.v.fill(0.0);
        flow.w.fill(0.0);
        // Seed well OFF the body's centreline. The 0.3 m box at the
        // origin spans only y,z ∈ [-0.15, 0.15]; a seed near the tunnel
        // wall (here a quarter of the height up) traces a straight line
        // that never enters the solid, so the trace is not cut short.
        // (Zeroing the turbulent viscosity does NOT remove the solid
        // mask — `trace_streamline` stops on `tunnel.body`, not on µ_t.)
        let seed = Vector3::new(
            tunnel.grid.x0 + 1.0,
            tunnel.grid.y0 + 0.25 * tunnel.grid.ly,
            tunnel.grid.z0 + 0.25 * tunnel.grid.lz,
        );
        let sl = trace_streamline(&tunnel, &flow, seed, 0.1, 50);
        assert!(sl.points.len() > 5, "streamline should advance");
        // The y and z coordinates should not have drifted.
        for p in &sl.points {
            assert!((p.y - seed.y).abs() < 1e-6);
            assert!((p.z - seed.z).abs() < 1e-6);
        }
        // The streamline length is positive.
        assert!(sl.length() > 0.0);
    }

    #[test]
    fn vorticity_of_uniform_flow_is_zero() {
        // A uniform velocity field has zero curl.
        let body = box_body(Vector3::zeros(), Vector3::new(0.3, 0.3, 0.3));
        let (_, mut flow) = quick_flow(&body);
        flow.u.fill(5.0);
        flow.v.fill(0.0);
        flow.w.fill(0.0);
        let vort = vorticity(&flow);
        assert!(vort.magnitude.abs_max() < 1e-9, "uniform flow has zero vorticity");
    }

    #[test]
    fn vorticity_of_solid_body_rotation_is_twice_omega() {
        // A solid-body rotation u = Ω×r about z has vorticity 2Ω·ẑ.
        let body = box_body(Vector3::zeros(), Vector3::new(0.3, 0.3, 0.3));
        let (tunnel, mut flow) = quick_flow(&body);
        let omega = 0.5;
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    let (cx, cy, _) = tunnel.grid.cell_centre(i, j, k);
                    let (rx, ry) = (cx, cy);
                    // u = (−Ω·y, Ω·x, 0).
                    let uc = -omega * ry;
                    let vc = omega * rx;
                    // Stamp onto the staggered faces (approx — uniform
                    // per cell is enough for a central-difference curl).
                    flow.u.set(i, j, k, uc);
                    flow.u.set(i + 1, j, k, uc);
                    flow.v.set(i, j, k, vc);
                    flow.v.set(i, j + 1, k, vc);
                }
            }
        }
        flow.w.fill(0.0);
        let vort = vorticity(&flow);
        // Sample an interior cell — ω_z should be ≈ 2Ω = 1.0.
        let (i, j, k) = (tunnel.grid.nx / 2, tunnel.grid.ny / 2, tunnel.grid.nz / 2);
        let wz = vort.omega_z.at(i, j, k);
        assert!(
            (wz - 2.0 * omega).abs() < 0.15,
            "solid-body rotation vorticity {wz} should be ~{}",
            2.0 * omega
        );
    }

    #[test]
    fn q_criterion_is_positive_in_a_pure_vortex() {
        // A pure rotation (no strain) has Q > 0 — rotation dominates.
        let body = box_body(Vector3::zeros(), Vector3::new(0.3, 0.3, 0.3));
        let (tunnel, mut flow) = quick_flow(&body);
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    let (cx, cy, _) = tunnel.grid.cell_centre(i, j, k);
                    flow.u.set(i, j, k, -cy);
                    flow.u.set(i + 1, j, k, -cy);
                    flow.v.set(i, j, k, cx);
                    flow.v.set(i, j + 1, k, cx);
                }
            }
        }
        flow.w.fill(0.0);
        let q = q_criterion(&flow);
        let (i, j, k) = (tunnel.grid.nx / 2, tunnel.grid.ny / 2, tunnel.grid.nz / 2);
        assert!(q.at(i, j, k) > 0.0, "pure rotation should have Q > 0");
    }

    #[test]
    fn slice_field_extracts_the_right_layer() {
        // A field that varies only with z, sliced on a z-plane, must
        // come out constant.
        let body = box_body(Vector3::zeros(), Vector3::new(0.5, 0.5, 0.5));
        let (tunnel, flow) = quick_flow(&body);
        let mut f = tunnel.grid.scalar_field();
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    f.set(i, j, k, k as f64);
                }
            }
        }
        let zc = tunnel.grid.z0 + 0.5 * tunnel.grid.lz;
        let s = slice_field(&tunnel.grid, &f, SliceAxis::Z, zc);
        assert_eq!(s.axis, SliceAxis::Z);
        assert_eq!((s.width, s.height), (tunnel.grid.nx, tunnel.grid.ny));
        // Every value in the slice is the same z layer index.
        let (lo, hi) = s.range();
        assert!((lo - hi).abs() < 1e-12, "z-slice of a z-varying field is flat");
        // A pressure slice should be finite everywhere.
        let ps = slice_field(&tunnel.grid, &flow.pressure, SliceAxis::Y, zc);
        assert!(ps.values.iter().all(|v| v.is_finite()));
    }
}


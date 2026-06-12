//! The typed request / response surface — the LLM- / agent-controllable
//! wind-tunnel API.
//!
//! Everything an external caller (a UI, a script, or an LLM agent
//! running a study) needs is here: a single [`AeroRequest`] describes
//! the case in plain typed fields, [`run_windtunnel`] executes it, and
//! a single [`AeroResult`] carries the whole answer — the
//! coefficients, the surface field, the flow field, the wake, the
//! convergence history. The builder methods on [`AeroRequest`] make a
//! valid request hard to get wrong, and every failure mode is a typed
//! [`AeroError`] with a stable `code()` an agent can branch on.

use nalgebra::Vector3;

use crate::compressible::{
    correct_coefficients, mach_number, speed_of_sound, CompressibleCoefficients,
};
use crate::domain::{BoundaryConditions, TunnelSizing, WindTunnel};
use crate::error::AeroError;
use crate::forces::{
    coefficients, integrate_forces, surface_field, surface_stats, AeroCoefficients, AeroForces,
    SurfaceStats,
};
use crate::geometry::TriMesh;
use crate::postprocess::{wake_survey, WakeSurvey};
use crate::solver::{solve_steady, BodyMotion, FlowField, SolverControls};
use crate::turbulence::TurbulenceModel;
use crate::wind::{Air, Wind};

/// A fully-typed wind-tunnel run request.
///
/// Build one with [`AeroRequest::new`] (which sets straight sea-level
/// flow at the given speed) and the `with_*` methods, or fill the
/// fields directly. The defaults are a sensible external-aero case;
/// nothing is required beyond the free-stream speed.
#[derive(Clone, Copy, Debug)]
pub struct AeroRequest {
    /// Free-stream speed `U∞` (m·s⁻¹).
    pub speed: f64,
    /// Yaw angle of the wind (radians).
    pub yaw: f64,
    /// Pitch angle / angle of attack of the wind (radians).
    pub pitch: f64,
    /// The air the tunnel is filled with.
    pub air: Air,
    /// Upstream turbulence intensity (fraction).
    pub turbulence_intensity: f64,
    /// The turbulence model to run.
    pub turbulence: TurbulenceModel,
    /// The boundary-condition set.
    pub boundary: BoundaryConditions,
    /// The tunnel-sizing policy.
    pub sizing: TunnelSizing,
    /// Maximum SIMPLE outer iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the mass-imbalance residual.
    pub tolerance: f64,
    /// The absolute air temperature (K) — used for the Mach number /
    /// compressibility correction.
    pub temperature_k: f64,
    /// If `true`, apply the Prandtl-Glauert compressibility correction
    /// to the reported coefficients when the Mach number warrants it.
    pub apply_compressibility: bool,
}

impl AeroRequest {
    /// A straight (zero-yaw, zero-pitch) sea-level run at the given
    /// speed, with k-ω SST turbulence and the standard external-aero
    /// boundary set.
    pub fn new(speed: f64) -> AeroRequest {
        AeroRequest {
            speed,
            yaw: 0.0,
            pitch: 0.0,
            air: Air::sea_level(),
            turbulence_intensity: 0.01,
            turbulence: TurbulenceModel::KOmegaSST,
            boundary: BoundaryConditions::external_aero(),
            sizing: TunnelSizing::default(),
            max_iterations: SolverControls::default().max_iterations,
            tolerance: SolverControls::default().tolerance,
            temperature_k: 288.0,
            apply_compressibility: false,
        }
    }

    /// Set the angle of attack (the wind's pitch angle, radians).
    pub fn with_angle_of_attack(mut self, pitch: f64) -> AeroRequest {
        self.pitch = pitch;
        self
    }

    /// Set the yaw angle (radians).
    pub fn with_yaw(mut self, yaw: f64) -> AeroRequest {
        self.yaw = yaw;
        self
    }

    /// Choose the turbulence model.
    pub fn with_turbulence(mut self, model: TurbulenceModel) -> AeroRequest {
        self.turbulence = model;
        self
    }

    /// Set the boundary-condition set (e.g. an automotive moving
    /// ground).
    pub fn with_boundary(mut self, bc: BoundaryConditions) -> AeroRequest {
        self.boundary = bc;
        self
    }

    /// Set the tunnel-sizing policy.
    pub fn with_sizing(mut self, sizing: TunnelSizing) -> AeroRequest {
        self.sizing = sizing;
        self
    }

    /// Cap the SIMPLE outer iterations.
    pub fn with_max_iterations(mut self, iterations: usize) -> AeroRequest {
        self.max_iterations = iterations;
        self
    }

    /// Enable the Prandtl-Glauert compressibility correction.
    pub fn with_compressibility(mut self, on: bool) -> AeroRequest {
        self.apply_compressibility = on;
        self
    }

    /// Set the air explicitly (density + viscosity, e.g. altitude
    /// air).
    pub fn with_air(mut self, air: Air) -> AeroRequest {
        self.air = air;
        self
    }

    /// Build the validated [`Wind`] this request describes.
    pub fn wind(&self) -> Result<Wind, AeroError> {
        Wind::new(
            self.speed,
            self.yaw,
            self.pitch,
            self.air,
            self.turbulence_intensity,
        )
    }

    /// Build the [`SolverControls`] this request implies.
    pub fn controls(&self) -> SolverControls {
        let mut c = SolverControls {
            turbulence: self.turbulence,
            max_iterations: self.max_iterations.max(1),
            tolerance: self.tolerance,
            ..SolverControls::default()
        };
        // A laminar request runs at the laminar default relaxation.
        if self.turbulence == TurbulenceModel::Laminar {
            c.relax_u = 0.6;
            c.relax_p = 0.25;
        }
        c
    }
}

/// The complete result of a wind-tunnel run.
#[derive(Clone, Debug)]
pub struct AeroResult {
    /// The built tunnel (grid, voxelized body, reference quantities).
    pub tunnel: WindTunnel,
    /// The converged flow field.
    pub flow: FlowField,
    /// The integrated aerodynamic force / moment.
    pub forces: AeroForces,
    /// The non-dimensional aerodynamic coefficients.
    pub coefficients: AeroCoefficients,
    /// Summary statistics of the body-surface field (Cp / y+ ranges).
    pub surface: SurfaceStats,
    /// A wake survey behind the body.
    pub wake: WakeSurvey,
    /// The compressibility-corrected coefficients, if the request
    /// enabled the correction.
    pub compressible: Option<CompressibleCoefficients>,
    /// `true` if the solver reached the convergence tolerance.
    pub converged: bool,
    /// The free-stream Reynolds number of the run.
    pub reynolds_number: f64,
    /// The free-stream Mach number of the run.
    pub mach_number: f64,
}

impl AeroResult {
    /// The headline drag coefficient.
    pub fn drag_coefficient(&self) -> f64 {
        self.coefficients.cd
    }

    /// The headline lift coefficient.
    pub fn lift_coefficient(&self) -> f64 {
        self.coefficients.cl
    }

    /// The drag *area* `Cd·A` (m²) — the quantity that actually sets a
    /// vehicle's drag force, independent of the reference-area choice.
    pub fn drag_area(&self) -> f64 {
        self.coefficients.cd * self.tunnel.reference_area
    }
}

/// Run a complete virtual-wind-tunnel study on a body.
///
/// This is the one-call entry point: it validates the request, builds
/// the tunnel around the body, solves the steady 3-D flow, integrates
/// the forces, computes the coefficients and the surface / wake
/// fields, optionally applies the compressibility correction, and
/// returns the whole [`AeroResult`].
///
/// Returns an [`AeroError`] only for genuine *input* failures (a bad
/// body, an ill-posed wind, a body that will not fit the domain) — a
/// solve that hits the iteration cap before converging still returns a
/// result, with `converged == false`.
pub fn run_windtunnel(body: &TriMesh, request: &AeroRequest) -> Result<AeroResult, AeroError> {
    let wind = request.wind()?;
    let tunnel = WindTunnel::build_with(body, wind, request.boundary, request.sizing)?;
    run_on_tunnel(&tunnel, request)
}

/// Run a study on an already-built tunnel — used by the angle-of-attack
/// sweep, which builds one tunnel and re-orients the wind.
pub(crate) fn run_on_tunnel(
    tunnel: &WindTunnel,
    request: &AeroRequest,
) -> Result<AeroResult, AeroError> {
    let controls = request.controls();
    let flow = solve_steady(tunnel, &controls, &BodyMotion::static_body());

    let moment_ref = body_centre(tunnel);
    let forces = integrate_forces(tunnel, &flow, moment_ref);
    let coeff = coefficients(tunnel, &forces);
    let surf = surface_stats(&surface_field(tunnel, &flow));

    // Wake survey: 0.5 body-lengths behind, running 3 body-lengths.
    let l = tunnel.reference_length.max(1e-6);
    let wake = wake_survey(tunnel, &flow, 0.5 * l, 3.0 * l, 24);

    // Mach number / compressibility.
    let a = speed_of_sound(request.temperature_k);
    let mach = mach_number(tunnel.wind.speed, a);
    let compressible = if request.apply_compressibility {
        Some(correct_coefficients(coeff.cd, coeff.cl, mach))
    } else {
        None
    };

    Ok(AeroResult {
        converged: flow.converged,
        reynolds_number: tunnel.reynolds_number(),
        mach_number: mach,
        tunnel: tunnel.clone(),
        flow,
        forces,
        coefficients: coeff,
        surface: surf,
        wake,
        compressible,
    })
}

/// The body's bounding-box centre (the moment reference point).
fn body_centre(tunnel: &WindTunnel) -> Vector3<f64> {
    let g = tunnel.grid;
    let mut min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                if tunnel.body.is_solid(i, j, k) {
                    let (cx, cy, cz) = g.cell_centre(i, j, k);
                    let c = Vector3::new(cx, cy, cz);
                    min = min.inf(&c);
                    max = max.sup(&c);
                    any = true;
                }
            }
        }
    }
    if any {
        0.5 * (min + max)
    } else {
        Vector3::new(g.x0 + 0.5 * g.lx, g.y0 + 0.5 * g.ly, g.z0 + 0.5 * g.lz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{box_body, sphere_body};

    #[test]
    fn request_builder_sets_fields() {
        let r = AeroRequest::new(40.0)
            .with_angle_of_attack(0.1)
            .with_yaw(0.05)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_max_iterations(50)
            .with_compressibility(true);
        assert_eq!(r.speed, 40.0);
        assert!((r.pitch - 0.1).abs() < 1e-12);
        assert!((r.yaw - 0.05).abs() < 1e-12);
        assert_eq!(r.turbulence, TurbulenceModel::KEpsilon);
        assert_eq!(r.max_iterations, 50);
        assert!(r.apply_compressibility);
    }

    #[test]
    fn request_builds_a_valid_wind_and_controls() {
        let r = AeroRequest::new(30.0);
        let w = r.wind().unwrap();
        assert!((w.speed - 30.0).abs() < 1e-12);
        let c = r.controls();
        assert_eq!(c.turbulence, TurbulenceModel::KOmegaSST);
    }

    #[test]
    fn run_windtunnel_rejects_a_bad_request() {
        // A negative speed is an ill-posed wind.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let r = AeroRequest::new(-10.0);
        let err = run_windtunnel(&body, &r).unwrap_err();
        assert_eq!(err.code(), "aero.bad_parameter");
    }

    #[test]
    fn run_windtunnel_rejects_an_empty_body() {
        let r = AeroRequest::new(20.0);
        let err = run_windtunnel(&TriMesh::new(), &r).unwrap_err();
        assert_eq!(err.category(), crate::ErrorCategory::Input);
    }

    /// A coarse but real grid for the api smoke tests — they assert
    /// result *completeness* (every field populated, finite, the right
    /// sign), which holds at any resolution. A coarse grid keeps the
    /// real end-to-end solve fast.
    fn coarse_sizing() -> TunnelSizing {
        TunnelSizing {
            cells_across_body: 4,
            max_cells: 40_000,
            ..TunnelSizing::default()
        }
    }

    #[test]
    fn run_windtunnel_produces_a_complete_result() {
        // The end-to-end smoke test: a sphere, run, every field of the
        // result must be populated and finite.
        let body = sphere_body(Vector3::zeros(), 0.5, 16, 32);
        let r = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse_sizing())
            .with_max_iterations(40);
        let result = run_windtunnel(&body, &r).unwrap();
        assert!(result.coefficients.cd.is_finite());
        assert!(result.coefficients.cl.is_finite());
        assert!(result.reynolds_number > 0.0);
        assert!(result.mach_number > 0.0 && result.mach_number < 0.3);
        assert!(!result.wake.positions.is_empty());
        assert!(result.surface.face_count > 0);
        // The drag area is Cd·A.
        let da = result.drag_area();
        assert!((da - result.coefficients.cd * result.tunnel.reference_area).abs() < 1e-9);
        // No compressibility correction was requested.
        assert!(result.compressible.is_none());
    }

    #[test]
    fn compressibility_correction_appears_when_requested() {
        // A fast run with the correction enabled — the result must
        // carry the compressible coefficients.
        let body = sphere_body(Vector3::zeros(), 0.3, 16, 32);
        let r = AeroRequest::new(200.0) // ~Mach 0.6
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse_sizing())
            .with_max_iterations(20)
            .with_compressibility(true);
        let result = run_windtunnel(&body, &r).unwrap();
        assert!(result.compressible.is_some());
        let comp = result.compressible.unwrap();
        // Mach ~0.6 → the regime is subsonic.
        assert!(result.mach_number > 0.3);
        // The compressible coefficients carry the same Mach number and
        // a finite Prandtl-Glauert factor.
        assert!((comp.mach - result.mach_number).abs() < 1e-12);
        assert!(comp.beta.is_finite() && comp.beta > 0.0);
    }

    #[test]
    fn drag_coefficient_of_a_sphere_is_plausible() {
        // A sphere at a high Reynolds number has a textbook Cd around
        // 0.1–0.5 (the exact value depends on whether the boundary
        // layer has tripped). The immersed-boundary v1 won't pin it,
        // but it must be a positive O(0.1–1) number, not garbage. A
        // coarse grid keeps the steady SIMPLE solve fast — the plausible
        // band is wide and does not need a fine mesh.
        let body = sphere_body(Vector3::zeros(), 0.5, 20, 40);
        let r = AeroRequest::new(25.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..TunnelSizing::default()
            })
            .with_max_iterations(60);
        let result = run_windtunnel(&body, &r).unwrap();
        let cd = result.drag_coefficient();
        assert!(
            cd > 0.05 && cd < 3.0,
            "sphere Cd {cd} is outside any plausible range"
        );
    }
}

//! Darcy-Weisbach head loss and the end-to-end pipe-flow solve.
//!
//! The Darcy-Weisbach equation gives the head lost to friction over a
//! straight run of pipe,
//!
//! ```text
//!   h_f = f * (L / D) * v^2 / (2 g)
//! ```
//!
//! where `f` is the Darcy friction factor (from [`crate::friction`]),
//! `L` the pipe length (m), `D` the internal diameter (m), `v` the bulk
//! velocity (m/s) and `g` the gravitational acceleration. `h_f` is a
//! *head* — a length, in metres of the flowing fluid.
//!
//! The corresponding **pressure drop** is `dP = rho * g * h_f`, i.e.
//!
//! ```text
//!   dP = f * (L / D) * rho * v^2 / 2
//! ```

use crate::error::{require_positive, PipeFlowError};
use crate::friction::{friction_factor, FrictionResult};
use crate::reynolds::reynolds_number;

/// Standard gravitational acceleration at the Earth's surface (m/s^2),
/// the conventional value used to convert between head and pressure.
pub const G_STANDARD: f64 = 9.806_65;

/// Darcy-Weisbach frictional head loss `h_f = f (L/D) v^2 / (2 g)`,
/// using standard gravity [`G_STANDARD`].
///
/// Inputs are SI: `friction_factor` dimensionless, `length` and
/// `diameter` in m, `velocity` in m/s. Returns the head loss in metres
/// of fluid. `length`, `diameter` and `velocity` must be finite and
/// strictly positive; `friction_factor` must be finite and positive.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::headloss::head_loss;
///
/// // f=0.02, L=100 m, D=0.1 m, v=2 m/s.
/// let hf = head_loss(0.02, 100.0, 0.1, 2.0).unwrap();
/// // 0.02 * (100/0.1) * 4 / (2*9.80665) = 4.0789 m.
/// assert!((hf - 4.0789).abs() < 1e-3);
/// ```
pub fn head_loss(
    friction_factor: f64,
    length: f64,
    diameter: f64,
    velocity: f64,
) -> Result<f64, PipeFlowError> {
    head_loss_g(friction_factor, length, diameter, velocity, G_STANDARD)
}

/// Darcy-Weisbach head loss with a caller-supplied gravitational
/// acceleration `g` (m/s^2). See [`head_loss`].
pub fn head_loss_g(
    friction_factor: f64,
    length: f64,
    diameter: f64,
    velocity: f64,
    g: f64,
) -> Result<f64, PipeFlowError> {
    let f = require_positive("friction_factor", friction_factor)?;
    let length = require_positive("length", length)?;
    let diameter = require_positive("diameter", diameter)?;
    let velocity = require_positive("velocity", velocity)?;
    let g = require_positive("g", g)?;
    Ok(f * (length / diameter) * velocity * velocity / (2.0 * g))
}

/// Frictional **pressure drop** `dP = f (L/D) rho v^2 / 2` (Pa).
///
/// This is the Darcy-Weisbach loss expressed as a pressure rather than a
/// head; it does not depend on `g`. `rho` is the fluid density (kg/m^3).
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::headloss::pressure_drop;
///
/// // Water (rho=1000), f=0.02, L=100 m, D=0.1 m, v=2 m/s.
/// let dp = pressure_drop(0.02, 100.0, 0.1, 2.0, 1000.0).unwrap();
/// // 0.02 * (100/0.1) * 1000 * 4 / 2 = 40_000 Pa.
/// assert!((dp - 40_000.0).abs() < 1e-6);
/// ```
pub fn pressure_drop(
    friction_factor: f64,
    length: f64,
    diameter: f64,
    velocity: f64,
    rho: f64,
) -> Result<f64, PipeFlowError> {
    let f = require_positive("friction_factor", friction_factor)?;
    let length = require_positive("length", length)?;
    let diameter = require_positive("diameter", diameter)?;
    let velocity = require_positive("velocity", velocity)?;
    let rho = require_positive("rho", rho)?;
    Ok(f * (length / diameter) * rho * velocity * velocity / 2.0)
}

/// Mean **wall shear stress** `tau_w = f * rho * v^2 / 8` (Pa).
///
/// The Darcy friction factor `f` is defined so the shear the flow exerts
/// on the pipe wall is `tau_w = (f / 8) rho v^2`. Equivalently it is the
/// force balance over a length `L` of pipe, `tau_w = dP * D / (4 L)`,
/// with `dP` the Darcy-Weisbach [`pressure_drop`]: the pressure force on
/// the cross-section `dP * (pi D^2 / 4)` is carried by the wall shear
/// `tau_w * (pi D L)`. It is therefore a *local* quantity, independent of
/// pipe length and diameter at a fixed friction factor.
///
/// `friction_factor`, `rho` and `velocity` must be finite and strictly
/// positive.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::headloss::wall_shear_stress;
///
/// // f=0.02, rho=1000, v=2: tau = 0.02 * 1000 * 4 / 8 = 10 Pa.
/// let tau = wall_shear_stress(0.02, 1000.0, 2.0).unwrap();
/// assert!((tau - 10.0).abs() < 1e-9);
/// ```
pub fn wall_shear_stress(
    friction_factor: f64,
    rho: f64,
    velocity: f64,
) -> Result<f64, PipeFlowError> {
    let f = require_positive("friction_factor", friction_factor)?;
    let rho = require_positive("rho", rho)?;
    let velocity = require_positive("velocity", velocity)?;
    Ok(f * rho * velocity * velocity / 8.0)
}

/// **Friction velocity** `u* = sqrt(tau_w / rho) = v * sqrt(f / 8)` (m/s).
///
/// The velocity scale of the near-wall turbulence, built from the
/// [`wall_shear_stress`] `tau_w`; it sets the wall units (`y+`, `u+`) of
/// the law of the wall. Equivalently `u* = v sqrt(f / 8)` straight from
/// the Darcy friction factor.
///
/// `friction_factor`, `rho` and `velocity` must be finite and strictly
/// positive.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::headloss::friction_velocity;
///
/// // f=0.02, rho=1000, v=2: u* = 2 * sqrt(0.02/8) = 0.1 m/s.
/// let u_star = friction_velocity(0.02, 1000.0, 2.0).unwrap();
/// assert!((u_star - 0.1).abs() < 1e-9);
/// ```
pub fn friction_velocity(
    friction_factor: f64,
    rho: f64,
    velocity: f64,
) -> Result<f64, PipeFlowError> {
    let tau = wall_shear_stress(friction_factor, rho, velocity)?;
    Ok((tau / rho).sqrt())
}

/// The fully solved state of a straight pipe run: Reynolds number,
/// friction factor + regime, head loss and pressure drop. Returned by
/// [`solve_pipe`].
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PipeFlowResult {
    /// The friction characterisation (Reynolds, factor, regime).
    pub friction: FrictionResult,
    /// Darcy-Weisbach head loss (m of fluid).
    pub head_loss_m: f64,
    /// Frictional pressure drop (Pa).
    pub pressure_drop_pa: f64,
}

/// End-to-end pipe-flow solve from physical inputs.
///
/// Given the fluid (`rho`, dynamic `viscosity`), the geometry
/// (`diameter`, `length`, `relative_roughness = epsilon/D`) and the
/// bulk `velocity`, this:
///
/// 1. computes the Reynolds number,
/// 2. classifies the regime and selects the appropriate friction
///    correlation (`64/Re` laminar, Haaland turbulent),
/// 3. evaluates the Darcy-Weisbach head loss and pressure drop,
///
/// returning everything in a [`PipeFlowResult`]. Standard gravity
/// [`G_STANDARD`] is used for the head/pressure conversion.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::headloss::solve_pipe;
/// use valenx_pipeflow::reynolds::FlowRegime;
///
/// // Water at 20 C through 100 m of 100 mm commercial-steel pipe at 2 m/s.
/// let r = solve_pipe(998.0, 1.002e-3, 0.1, 100.0, 4.6e-4 / 0.1, 2.0).unwrap();
/// assert_eq!(r.friction.regime, FlowRegime::Turbulent);
/// assert!(r.head_loss_m > 0.0);
/// assert!(r.pressure_drop_pa > 0.0);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn solve_pipe(
    rho: f64,
    viscosity: f64,
    diameter: f64,
    length: f64,
    relative_roughness: f64,
    velocity: f64,
) -> Result<PipeFlowResult, PipeFlowError> {
    let re = reynolds_number(rho, velocity, diameter, viscosity)?;
    let friction = friction_factor(re, relative_roughness)?;
    let head_loss_m = head_loss(friction.friction_factor, length, diameter, velocity)?;
    let pressure_drop_pa =
        pressure_drop(friction.friction_factor, length, diameter, velocity, rho)?;
    Ok(PipeFlowResult {
        friction,
        head_loss_m,
        pressure_drop_pa,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Head loss against a hand-computed reference value.
    #[test]
    fn head_loss_matches_hand_calculation() {
        // f=0.02, L=100, D=0.1, v=2, g=9.80665:
        // hf = 0.02 * 1000 * 4 / (2*9.80665) = 80 / 19.6133 = 4.07886 m.
        let hf = head_loss(0.02, 100.0, 0.1, 2.0).unwrap();
        assert!((hf - 4.078_862).abs() < 1e-5, "hf={hf}");
    }

    /// Head loss scales with the square of velocity.
    #[test]
    fn head_loss_scales_with_velocity_squared() {
        let base = head_loss(0.02, 100.0, 0.1, 2.0).unwrap();
        let doubled = head_loss(0.02, 100.0, 0.1, 4.0).unwrap();
        // Doubling v with f fixed quadruples the loss.
        assert!(
            (doubled - 4.0 * base).abs() < 1e-9,
            "base={base}, doubled={doubled}"
        );
    }

    /// Head loss scales linearly with pipe length.
    #[test]
    fn head_loss_scales_linearly_with_length() {
        let base = head_loss(0.02, 100.0, 0.1, 2.0).unwrap();
        let triple = head_loss(0.02, 300.0, 0.1, 2.0).unwrap();
        assert!(
            (triple - 3.0 * base).abs() < 1e-9,
            "base={base}, triple={triple}"
        );
    }

    /// Head loss scales inversely with diameter (at fixed velocity).
    #[test]
    fn head_loss_scales_inversely_with_diameter() {
        let base = head_loss(0.02, 100.0, 0.1, 2.0).unwrap();
        let half_d = head_loss(0.02, 100.0, 0.05, 2.0).unwrap();
        // Halving D doubles L/D, hence doubles the loss.
        assert!(
            (half_d - 2.0 * base).abs() < 1e-9,
            "base={base}, half_d={half_d}"
        );
    }

    /// Custom gravity scales head loss inversely with g.
    #[test]
    fn head_loss_scales_inversely_with_gravity() {
        let earth = head_loss_g(0.02, 100.0, 0.1, 2.0, 9.80665).unwrap();
        let twice_g = head_loss_g(0.02, 100.0, 0.1, 2.0, 2.0 * 9.80665).unwrap();
        assert!((twice_g - earth / 2.0).abs() < 1e-9);
    }

    /// Pressure drop against a hand-computed reference value.
    #[test]
    fn pressure_drop_matches_hand_calculation() {
        // f=0.02, L=100, D=0.1, v=2, rho=1000:
        // dP = 0.02 * 1000 * 1000 * 4 / 2 = 40_000 Pa.
        let dp = pressure_drop(0.02, 100.0, 0.1, 2.0, 1000.0).unwrap();
        assert!((dp - 40_000.0).abs() < 1e-6, "dp={dp}");
    }

    /// `dP = rho * g * h_f` ties the head-loss and pressure-drop forms
    /// together exactly.
    #[test]
    fn pressure_drop_equals_rho_g_head_loss() {
        let rho = 998.0;
        let f = 0.025;
        let hf = head_loss(f, 50.0, 0.08, 1.5).unwrap();
        let dp = pressure_drop(f, 50.0, 0.08, 1.5, rho).unwrap();
        assert!((dp - rho * G_STANDARD * hf).abs() < 1e-6, "dp={dp}");
    }

    /// Pressure drop also scales with v^2 and with length.
    #[test]
    fn pressure_drop_scaling() {
        let base = pressure_drop(0.02, 100.0, 0.1, 2.0, 1000.0).unwrap();
        let v2 = pressure_drop(0.02, 100.0, 0.1, 4.0, 1000.0).unwrap();
        let l2 = pressure_drop(0.02, 200.0, 0.1, 2.0, 1000.0).unwrap();
        assert!((v2 - 4.0 * base).abs() < 1e-6);
        assert!((l2 - 2.0 * base).abs() < 1e-6);
    }

    /// Wall shear stress against a hand-computed reference value.
    #[test]
    fn wall_shear_stress_matches_hand_calculation() {
        // f=0.02, rho=1000, v=2: tau = 0.02 * 1000 * 4 / 8 = 10 Pa.
        let tau = wall_shear_stress(0.02, 1000.0, 2.0).unwrap();
        assert!((tau - 10.0).abs() < 1e-9, "tau={tau}");
    }

    /// The momentum balance `tau_w = dP * D / (4 L)` ties the wall shear
    /// to the Darcy-Weisbach pressure drop exactly, for any `L` and `D`.
    #[test]
    fn wall_shear_equals_pressure_drop_force_balance() {
        let (f, rho, v) = (0.025, 998.0, 1.5);
        for &(l, d) in &[(50.0, 0.08), (100.0, 0.1), (10.0, 0.025)] {
            let tau = wall_shear_stress(f, rho, v).unwrap();
            let dp = pressure_drop(f, l, d, v, rho).unwrap();
            assert!(
                (tau - dp * d / (4.0 * l)).abs() < 1e-9,
                "tau={tau}, dP*D/4L={}",
                dp * d / (4.0 * l)
            );
        }
    }

    /// Wall shear scales with the square of velocity (`f`, `rho` fixed).
    #[test]
    fn wall_shear_scales_with_velocity_squared() {
        let base = wall_shear_stress(0.02, 1000.0, 2.0).unwrap();
        let doubled = wall_shear_stress(0.02, 1000.0, 4.0).unwrap();
        assert!(
            (doubled - 4.0 * base).abs() < 1e-9,
            "base={base}, doubled={doubled}"
        );
    }

    /// Friction velocity against the hand value and both defining
    /// identities `u* = sqrt(tau_w/rho)` and `u* = v sqrt(f/8)`.
    #[test]
    fn friction_velocity_matches_definitions() {
        let (f, rho, v) = (0.02, 1000.0, 2.0);
        let u_star = friction_velocity(f, rho, v).unwrap();
        // u* = v sqrt(f/8) = 2 * sqrt(0.0025) = 0.1 m/s.
        assert!((u_star - 0.1).abs() < 1e-12, "u*={u_star}");
        let tau = wall_shear_stress(f, rho, v).unwrap();
        assert!((u_star - (tau / rho).sqrt()).abs() < 1e-12);
        assert!((u_star - v * (f / 8.0).sqrt()).abs() < 1e-12);
    }

    /// Friction velocity is independent of density (`rho` cancels in
    /// `u* = v sqrt(f/8)`).
    #[test]
    fn friction_velocity_is_independent_of_density() {
        let a = friction_velocity(0.03, 1000.0, 2.5).unwrap();
        let b = friction_velocity(0.03, 13_600.0, 2.5).unwrap(); // mercury
        assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
    }

    /// Bad inputs are rejected by both wall-friction functions.
    #[test]
    fn wall_friction_rejects_bad_inputs() {
        assert!(wall_shear_stress(0.0, 1000.0, 2.0).is_err()); // f
        assert!(wall_shear_stress(0.02, -1.0, 2.0).is_err()); // rho
        assert!(wall_shear_stress(0.02, 1000.0, 0.0).is_err()); // v
        assert!(friction_velocity(0.02, f64::NAN, 2.0).is_err()); // non-finite
        assert!(friction_velocity(0.02, 1000.0, -2.0).is_err()); // v
    }

    /// Bad inputs are rejected.
    #[test]
    fn rejects_bad_inputs() {
        assert!(head_loss(0.0, 100.0, 0.1, 2.0).is_err());
        assert!(head_loss(0.02, -100.0, 0.1, 2.0).is_err());
        assert!(head_loss(0.02, 100.0, 0.0, 2.0).is_err());
        assert!(head_loss_g(0.02, 100.0, 0.1, 2.0, 0.0).is_err());
        assert!(pressure_drop(0.02, 100.0, 0.1, 2.0, f64::NAN).is_err());
    }

    /// The end-to-end solve in the laminar regime reproduces the
    /// analytic `64/Re` friction and a consistent head loss.
    #[test]
    fn solve_pipe_laminar_end_to_end() {
        // Pick conditions that are firmly laminar.
        // rho=900 (oil), mu=0.1 Pa.s, D=0.05, v=1  => Re = 900*1*0.05/0.1 = 450.
        let rho = 900.0;
        let mu = 0.1;
        let d = 0.05;
        let v = 1.0;
        let l = 10.0;
        let r = solve_pipe(rho, mu, d, l, 0.0, v).unwrap();
        let re_expected = rho * v * d / mu;
        assert!((r.friction.reynolds - re_expected).abs() < 1e-6);
        assert!((r.friction.friction_factor - 64.0 / re_expected).abs() < 1e-12);
        // hf = f (L/D) v^2/(2g).
        let f = 64.0 / re_expected;
        let hf_expected = f * (l / d) * v * v / (2.0 * G_STANDARD);
        assert!((r.head_loss_m - hf_expected).abs() < 1e-9);
        // dP = rho g hf.
        assert!((r.pressure_drop_pa - rho * G_STANDARD * r.head_loss_m).abs() < 1e-6);
    }

    /// The end-to-end solve in the turbulent regime returns a positive
    /// head loss and pressure drop and classifies turbulent.
    #[test]
    fn solve_pipe_turbulent_end_to_end() {
        use crate::reynolds::FlowRegime;
        let r = solve_pipe(998.0, 1.002e-3, 0.1, 100.0, 4.6e-4 / 0.1, 2.0).unwrap();
        assert_eq!(r.friction.regime, FlowRegime::Turbulent);
        assert!(r.head_loss_m > 0.0);
        assert!(r.pressure_drop_pa > 0.0);
        // Consistency between head and pressure forms.
        assert!((r.pressure_drop_pa - 998.0 * G_STANDARD * r.head_loss_m).abs() < 1e-3);
    }

    /// The whole result round-trips through JSON.
    #[test]
    fn result_serializes_round_trip() {
        let r = solve_pipe(998.0, 1.002e-3, 0.1, 100.0, 1e-4, 2.0).unwrap();
        let json = serde_json::to_string(&r).unwrap();
        let back: PipeFlowResult = serde_json::from_str(&json).unwrap();
        assert_eq!(r, back);
    }
}

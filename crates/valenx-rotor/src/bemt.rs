//! Blade-element-momentum theory (BEMT) solver.
//!
//! This module couples 2-D blade-element aerodynamics with annular
//! momentum theory and solves, at each radial station, for the local
//! inflow angle `phi` that makes the two consistent. The integrated blade
//! loads give thrust, torque, power, propeller efficiency and the hover
//! figure of merit.
//!
//! ## Per-element model
//!
//! At radius `r` (`hub <= r <= R`) with chord `c(r)`, geometric pitch /
//! twist `beta(r)` (radians), `n` blades, rotor speed `Omega` (rad/s) and
//! axial freestream `V`:
//!
//! - local solidity `sigma = n*c / (2*pi*r)`;
//! - angle of attack `alpha = beta - phi` (`phi` is the unknown inflow
//!   angle measured from the rotor plane);
//! - from the section polar, `Cl(alpha)`, `Cd(alpha)`;
//! - rotor-frame force coefficients in the **propeller** convention,
//!   `Cn = Cl*cos(phi) - Cd*sin(phi)` (thrust: drag subtracts),
//!   `Ct = Cl*sin(phi) + Cd*cos(phi)` (torque: drag adds). (The wind-
//!   turbine decomposition flips both drag signs; using it for a
//!   thrusting rotor pushes the figure of merit above 1.)
//! - Prandtl tip loss `f_tip = (n/2)(R-r)/(r*|sin phi|)`,
//!   `F_tip = (2/pi) acos(exp(-f_tip))`, hub loss analogously, total
//!   `F = F_tip * F_hub` (the acos argument is clamped to `[-1, 1]`,
//!   `sin phi` and the exp argument are guarded);
//! - the velocity triangle is fixed by `phi`: with `U_t = Omega r` the
//!   in-plane speed, `U_a = U_t tan(phi)` the through-disk axial speed,
//!   the induced axial velocity is `v_i = U_a - V` and `W^2 = U_a^2 + U_t^2`;
//! - **inflow solve**: `phi` is the root of the thrust-balance residual
//!   `g(phi) = dT_be/dr - dT_mom/dr`, where the blade element gives
//!   `dT_be/dr = 0.5 rho W^2 (n c) Cn` and annular momentum gives
//!   `dT_mom/dr = 4 pi rho r F U_a v_i`, with a Glauert/Buhl high-thrust
//!   correction in the windmill-brake state (turbine-convention axial
//!   induction `a = 1 - U_a/V` above ~0.4). The root is bracketed by
//!   scanning for the first (smallest-`phi`, physical) sign change and
//!   refined by bisection on `(eps, pi/2 - eps)`. This thrust-balance
//!   form is well-posed in **hover** (`V = 0`), where the bare
//!   `tan(phi) = V(1-a)/(Omega r (1+a'))` rearrangement degenerates to
//!   the trivial `phi = 0`.
//!
//! ## Span integration
//!
//! With the converged `W^2` per station, per-unit-span thrust and torque
//! are `dT/dr = 0.5 rho W^2 (n c) Cn` and `dQ/dr = 0.5 rho W^2 (n c) Ct r`,
//! integrated by the trapezoid rule over the stations to total thrust `T`
//! and torque `Q`. Then `P = Q*Omega`, propeller efficiency
//! `eta = T*V/P` (for `V > 0`) and the hover figure of merit
//! `FM = (T^1.5 / sqrt(2 rho A)) / P` (`A = pi R^2`).
//!
//! ## Honest scope
//!
//! This is 1-D strip theory with the standard engineering corrections
//! (Prandtl loss, Glauert/Buhl high-induction). It is research/educational
//! grade: it omits radial flow and 3-D / vortex effects, unsteady and
//! dynamic-stall aerodynamics, compressibility, Reynolds-number variation
//! of the polar, and rotor-rotor / ground interference. The analytic
//! polar in particular is a thin-airfoil idealisation. Results are a
//! reasonable first estimate, not a substitute for a vortex-method or CFD
//! analysis, and any quantitative use needs measured airfoil data and
//! validation against test data.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::airfoil::Polar;
use crate::error::{require_finite, require_positive, RotorError};

/// Hard cap on bisection iterations for the per-element inflow solve.
///
/// Bisection halves the bracket each step, so ~60 iterations resolves the
/// `(eps, pi/2)` interval to far below `f64` precision; the cap exists
/// purely to guarantee termination, not as a tuning knob.
const MAX_ITERS: u32 = 200;

/// Convergence tolerance on the inflow angle `phi` (radians) for the
/// bisection bracket width.
const PHI_TOL: f64 = 1e-10;

/// Small angular guard keeping the root search strictly inside
/// `(0, pi/2)` so `sin phi`, `cos phi` and the loss terms stay finite.
const ANGLE_EPS: f64 = 1e-6;

/// Induction factor above which the simple momentum relation is replaced
/// by the Glauert/Buhl high-thrust empirical correction.
const A_CRIT: f64 = 0.4;

/// One radial blade station: the local geometry of a strip of the blade.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BladeStation {
    /// Radial position `r` of the station (m), `hub <= r <= R`.
    pub radius_m: f64,
    /// Local chord `c(r)` (m), > 0.
    pub chord_m: f64,
    /// Local geometric twist / pitch `beta(r)` (radians), measured from
    /// the rotor plane. Finite (any sign).
    pub twist_rad: f64,
}

impl BladeStation {
    /// Build a validated blade station.
    ///
    /// # Errors
    ///
    /// Returns [`RotorError::NonPositive`] if radius or chord is not
    /// finite and positive, and [`RotorError::NotFinite`] if the twist is
    /// not finite.
    pub fn new(radius_m: f64, chord_m: f64, twist_rad: f64) -> Result<Self, RotorError> {
        let radius_m = require_positive("station radius", radius_m)?;
        let chord_m = require_positive("chord", chord_m)?;
        let twist_rad = require_finite("twist", twist_rad)?;
        Ok(Self {
            radius_m,
            chord_m,
            twist_rad,
        })
    }
}

/// A rotor / propeller: blade count, tip and hub radius, the radial
/// station table and a section polar. Validated on construction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Rotor {
    /// Number of blades `n` (>= 1).
    pub blade_count: u32,
    /// Tip radius `R` (m).
    pub tip_radius_m: f64,
    /// Hub (root) radius (m), `0 <= hub < R`. Stations inboard of the hub
    /// carry no load.
    pub hub_radius_m: f64,
    /// Radial stations, strictly increasing in radius, all within
    /// `[hub, tip]`. At least two are required for trapezoid integration.
    pub stations: Vec<BladeStation>,
    /// Section airfoil polar shared by all stations.
    pub polar: Polar,
}

impl Rotor {
    /// Build a validated rotor.
    ///
    /// Validates: `blade_count >= 1`; `tip_radius` finite and positive;
    /// `hub_radius` finite, `>= 0` and `< tip_radius`; at least two
    /// stations, strictly increasing in radius and all within
    /// `[hub_radius, tip_radius]`.
    ///
    /// # Errors
    ///
    /// Returns the matching [`RotorError`] variant for any violated
    /// constraint.
    pub fn new(
        blade_count: u32,
        tip_radius_m: f64,
        hub_radius_m: f64,
        stations: Vec<BladeStation>,
        polar: Polar,
    ) -> Result<Self, RotorError> {
        if blade_count == 0 {
            return Err(RotorError::NoBlades);
        }
        let tip_radius_m = require_positive("tip_radius", tip_radius_m)?;
        // Hub may be zero, but must be finite and below the tip.
        let hub_radius_m = require_finite("hub_radius", hub_radius_m)?;
        if hub_radius_m < 0.0 {
            return Err(RotorError::NonPositive {
                name: "hub_radius",
                value: hub_radius_m,
            });
        }
        // Both radii are finite here (validated above), so `>=` is a safe,
        // NaN-free rejection of a hub not strictly below the tip.
        if hub_radius_m >= tip_radius_m {
            return Err(RotorError::HubNotBelowTip {
                hub: hub_radius_m,
                tip: tip_radius_m,
            });
        }
        if stations.len() < 2 {
            return Err(RotorError::TooFewStations {
                count: stations.len(),
            });
        }
        for (i, s) in stations.iter().enumerate() {
            if !(s.radius_m.is_finite() && s.chord_m.is_finite() && s.twist_rad.is_finite()) {
                return Err(RotorError::BadStations {
                    reason: format!("station {i} has a non-finite field"),
                });
            }
            if s.radius_m < hub_radius_m || s.radius_m > tip_radius_m {
                return Err(RotorError::BadStations {
                    reason: format!(
                        "station {i} radius {} outside [hub {hub_radius_m}, tip {tip_radius_m}]",
                        s.radius_m
                    ),
                });
            }
            if i > 0 && s.radius_m <= stations[i - 1].radius_m {
                return Err(RotorError::BadStations {
                    reason: format!(
                        "station {i} radius {} not greater than previous {}",
                        s.radius_m,
                        stations[i - 1].radius_m
                    ),
                });
            }
        }
        Ok(Self {
            blade_count,
            tip_radius_m,
            hub_radius_m,
            stations,
            polar,
        })
    }

    /// Convenience constructor that builds the station table from parallel
    /// slices of radius, chord and twist (all the same length) and uses
    /// the default analytic polar.
    ///
    /// # Errors
    ///
    /// Returns [`RotorError::BadStations`] if the slices differ in length,
    /// plus any error from [`Rotor::new`] / [`BladeStation::new`].
    pub fn from_slices(
        blade_count: u32,
        tip_radius_m: f64,
        hub_radius_m: f64,
        radius_m: &[f64],
        chord_m: &[f64],
        twist_rad: &[f64],
    ) -> Result<Self, RotorError> {
        if radius_m.len() != chord_m.len() || radius_m.len() != twist_rad.len() {
            return Err(RotorError::BadStations {
                reason: format!(
                    "radius/chord/twist lengths differ: {} / {} / {}",
                    radius_m.len(),
                    chord_m.len(),
                    twist_rad.len()
                ),
            });
        }
        let mut stations = Vec::with_capacity(radius_m.len());
        for i in 0..radius_m.len() {
            stations.push(BladeStation::new(radius_m[i], chord_m[i], twist_rad[i])?);
        }
        Rotor::new(
            blade_count,
            tip_radius_m,
            hub_radius_m,
            stations,
            Polar::analytic_default(),
        )
    }

    /// Total rotor disk area `A = pi R^2` (m^2).
    pub fn disk_area(&self) -> f64 {
        PI * self.tip_radius_m * self.tip_radius_m
    }

    /// Solve the rotor at an operating point and integrate the loads.
    ///
    /// `rpm` is the rotor speed (rev/min), `freestream_v` the axial inflow
    /// speed `V` (m/s, `>= 0`; `V = 0` is hover) and `air_density` the air
    /// density `rho` (kg/m^3).
    ///
    /// # Errors
    ///
    /// - [`RotorError::NonPositive`] if `rpm` or `air_density` is not
    ///   finite and positive, or [`RotorError::NotFinite`] / sign error if
    ///   `freestream_v` is not finite or is negative.
    /// - [`RotorError::NoConvergence`] if the inflow-angle root-finder
    ///   fails to bracket or converge at any station — never a silent NaN.
    pub fn solve(
        &self,
        rpm: f64,
        freestream_v: f64,
        air_density: f64,
    ) -> Result<RotorPerformance, RotorError> {
        let rpm = require_positive("rpm", rpm)?;
        let air_density = require_positive("air_density", air_density)?;
        let v = require_finite("freestream_v", freestream_v)?;
        if v < 0.0 {
            return Err(RotorError::NonPositive {
                name: "freestream_v",
                value: v,
            });
        }

        let omega = rpm * 2.0 * PI / 60.0; // rad/s
        let n = f64::from(self.blade_count);
        let r_tip = self.tip_radius_m;
        let r_hub = self.hub_radius_m;

        // Solve each station, collecting per-unit-span dT/dr and dQ/dr.
        let mut elements = Vec::with_capacity(self.stations.len());
        for (i, st) in self.stations.iter().enumerate() {
            let el = solve_element(
                ElementInput {
                    blades: n,
                    r: st.radius_m,
                    chord: st.chord_m,
                    twist: st.twist_rad,
                    r_tip,
                    r_hub,
                    omega,
                    v,
                    rho: air_density,
                    polar: &self.polar,
                },
                i,
            )?;
            elements.push(el);
        }

        // Trapezoid integration of dT/dr and dQ/dr over the stations.
        let mut thrust = 0.0;
        let mut torque = 0.0;
        for w in elements.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            let dr = b.radius_m - a.radius_m; // > 0, stations increasing
            thrust += 0.5 * (a.dt_dr + b.dt_dr) * dr;
            torque += 0.5 * (a.dq_dr + b.dq_dr) * dr;
        }

        let power = torque * omega;
        // Propeller efficiency: useful power T*V over shaft power P. Only
        // meaningful in forward flight (V > 0); defined as 0 in hover.
        let efficiency = if v > 0.0 && power.abs() > f64::MIN_POSITIVE {
            thrust * v / power
        } else {
            0.0
        };
        // Hover figure of merit: ideal hover power over actual shaft power.
        // Only meaningful near hover (V ~ 0); defined as 0 otherwise.
        let area = self.disk_area();
        let figure_of_merit = if v == 0.0 && power.abs() > f64::MIN_POSITIVE && thrust > 0.0 {
            let ideal_power = thrust.powf(1.5) / (2.0 * air_density * area).sqrt();
            ideal_power / power
        } else {
            0.0
        };

        Ok(RotorPerformance {
            rpm,
            omega_rad_s: omega,
            freestream_v_m_s: v,
            air_density,
            disk_area_m2: area,
            thrust_n: thrust,
            torque_nm: torque,
            power_w: power,
            efficiency,
            figure_of_merit,
            elements,
        })
    }
}

/// Inputs to a single blade-element solve (internal).
struct ElementInput<'a> {
    blades: f64,
    r: f64,
    chord: f64,
    twist: f64,
    r_tip: f64,
    r_hub: f64,
    omega: f64,
    v: f64,
    rho: f64,
    polar: &'a Polar,
}

/// The converged aerodynamic state of one radial blade element.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElementResult {
    /// Radial position `r` (m).
    pub radius_m: f64,
    /// Converged inflow angle `phi` (radians).
    pub inflow_angle_rad: f64,
    /// Local angle of attack `alpha = beta - phi` (radians).
    pub alpha_rad: f64,
    /// Axial induction factor `a`.
    pub axial_induction: f64,
    /// Tangential (swirl) induction factor `a'`.
    pub tangential_induction: f64,
    /// Combined Prandtl tip*hub loss factor `F` in `(0, 1]`.
    pub tip_loss: f64,
    /// Section lift coefficient `Cl` at this element.
    pub cl: f64,
    /// Section drag coefficient `Cd` at this element.
    pub cd: f64,
    /// Per-unit-span thrust `dT/dr` at this element (N/m).
    pub dt_dr: f64,
    /// Per-unit-span torque `dQ/dr` at this element (N*m/m).
    pub dq_dr: f64,
}

/// The integrated performance of a rotor at one operating point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RotorPerformance {
    /// Rotor speed (rev/min).
    pub rpm: f64,
    /// Rotor speed `Omega` (rad/s).
    pub omega_rad_s: f64,
    /// Axial freestream `V` (m/s).
    pub freestream_v_m_s: f64,
    /// Air density `rho` (kg/m^3).
    pub air_density: f64,
    /// Disk area `A = pi R^2` (m^2).
    pub disk_area_m2: f64,
    /// Total thrust `T` (N).
    pub thrust_n: f64,
    /// Total torque `Q` (N*m).
    pub torque_nm: f64,
    /// Shaft power `P = Q*Omega` (W).
    pub power_w: f64,
    /// Propeller efficiency `eta = T*V / P` (`0` in hover, `V = 0`).
    pub efficiency: f64,
    /// Hover figure of merit `FM = (T^1.5 / sqrt(2 rho A)) / P` (`0` away
    /// from hover, `V > 0`).
    pub figure_of_merit: f64,
    /// Per-element converged state, one entry per radial station.
    pub elements: Vec<ElementResult>,
}

/// Compute the combined Prandtl tip+hub loss factor `F` at inflow angle
/// `phi`, guarding `sin phi`, the exponential argument and the acos
/// domain. Returns a value in `(0, 1]`.
fn prandtl_loss(blades: f64, r: f64, r_tip: f64, r_hub: f64, phi: f64) -> f64 {
    // Guard sin(phi): the loss factor is only used where the blade is
    // loaded; for a vanishing sin(phi) the geometric factor blows up and
    // F -> 1 (no correction), which is the physically correct limit.
    let s = phi.sin().abs();
    if s < ANGLE_EPS {
        return 1.0;
    }
    let factor = |dist: f64| -> f64 {
        // f = (n/2) * dist / (r * sin phi); F = (2/pi) acos(exp(-f)).
        let f = 0.5 * blades * dist / (r * s);
        // exp(-f) in (0, 1]; clamp the acos argument defensively.
        let arg = (-f).exp().clamp(-1.0, 1.0);
        (2.0 / PI) * arg.acos()
    };
    let f_tip = factor((r_tip - r).max(0.0));
    let f_hub = factor((r - r_hub).max(0.0));
    // Combined loss; never let it reach exactly 0 (it divides later).
    (f_tip * f_hub).max(1e-6)
}

/// The blade-element / momentum thrust-balance residual at trial inflow
/// angle `phi`, returning the residual value and the element's auxiliary
/// quantities at that `phi`.
///
/// We solve for `phi` by balancing the elemental thrust two ways, which
/// — unlike the bare `tan phi = V(1-a)/(Omega r (1+a'))` rearrangement —
/// stays well-posed in **hover** (`V = 0`), where that rearrangement
/// collapses to the trivial `phi = 0`. At trial `phi`:
///
/// - the inflow angle fixes the velocity triangle: with `U_t = Omega r`
///   the in-plane speed, the through-disk axial speed is
///   `U_a = U_t tan phi`, so the induced axial velocity (propeller
///   convention) is `v_i = U_a - V`;
/// - the relative speed is `W = sqrt(U_a^2 + U_t^2)`, the angle of attack
///   `alpha = beta - phi`, and from the polar `Cl, Cd`, giving the
///   rotor-frame `Cn = Cl cos phi + Cd sin phi`;
/// - blade-element thrust per span `dT_be/dr = 0.5 rho W^2 (n c) Cn`;
/// - annular-momentum thrust per span
///   `dT_mom/dr = 4 pi rho r F U_a v_i` (Prandtl loss `F`).
///
/// The residual is `g(phi) = dT_be/dr - dT_mom/dr`. As `phi -> 0+` the
/// blade is at high `alpha` so `dT_be > 0` while `dT_mom -> 0` (positive
/// residual); as `phi -> beta` the lift falls to zero while the momentum
/// sink grows (negative residual) — guaranteeing a sign change to
/// bracket, in hover and in forward flight alike. The everywhere-finite
/// rotor-frame quantities make the residual robust at the tip (`F -> 0`)
/// and for negative-lift elements.
fn residual(inp: &ElementInput, phi: f64) -> (f64, ElementState) {
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let tan_phi = sin_phi / cos_phi.max(ANGLE_EPS);
    let alpha = inp.twist - phi;
    let (cl, cd) = inp.polar.coefficients(alpha);

    // Velocity triangle at this trial inflow angle.
    let u_t = inp.omega * inp.r; // in-plane speed (swirl added later)
    let u_a = u_t * tan_phi; // through-disk axial speed
    let v_i = u_a - inp.v; // induced axial velocity (propeller convention)
    let w2 = u_a * u_a + u_t * u_t; // relative speed squared

    // Rotor-frame force coefficients (propeller convention). With phi the
    // inflow angle from the rotor plane, lift L is perpendicular to the
    // relative wind and drag D parallel to it, so the axial (thrust) and
    // in-plane (torque-absorbing) projections are
    //   Cn = Cl cos phi - Cd sin phi   (thrust: drag subtracts),
    //   Ct = Cl sin phi + Cd cos phi   (torque: drag adds).
    // NOTE: this is the PROPELLER decomposition. The wind-TURBINE form
    // (Cn = Cl cos phi + Cd sin phi, Ct = Cl sin phi - Cd cos phi) flips
    // both drag signs because a turbine extracts power; using the turbine
    // form for a thrusting/hovering rotor over-predicts thrust and
    // under-predicts torque, which drives the figure of merit above 1.
    let cn = cl * cos_phi - cd * sin_phi;
    let ct = cl * sin_phi + cd * cos_phi;

    let f_loss = prandtl_loss(inp.blades, inp.r, inp.r_tip, inp.r_hub, phi);

    // Blade-element vs annular-momentum elemental thrust (per unit span).
    let dt_be = 0.5 * inp.rho * w2 * inp.blades * inp.chord * cn;
    let dt_mom = momentum_thrust_per_span(inp.rho, inp.r, f_loss, inp.v, u_a, v_i);
    let res = dt_be - dt_mom;

    (
        res,
        ElementState {
            phi,
            alpha,
            u_a,
            v_i,
            w2,
            f_loss,
            cl,
            cd,
            cn,
            ct,
        },
    )
}

/// Annular-momentum elemental thrust per unit span (N/m), with the
/// Glauert/Buhl high-induction correction.
///
/// For an axial-flow rotor the elemental momentum thrust is
/// `dT/dr = 4 pi rho r F U_a v_i`, where `U_a` is the through-disk axial
/// speed, `v_i = U_a - V` the induced axial velocity and `F` the Prandtl
/// loss. This simple (Froude / actuator-annulus) relation is exact and
/// well-behaved for a **thrust-producing** disk (`U_a >= V`, `v_i >= 0`)
/// in every axial state from hover (`V = 0`, `U_a = v_i`) through cruise.
///
/// The Glauert/Buhl empirical high-thrust correction addresses the
/// **windmill-brake / turbulent-wake state**, in which the disk extracts
/// energy and the *turbine-convention* axial induction
/// `a = 1 - U_a/V` exceeds the critical `A_CRIT` (~0.4); there the simple
/// momentum parabola `CT = 4 F a (1-a)` over-predicts and is replaced by
/// the Buhl (2005) line `CT = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2`. We
/// therefore gate the correction on that turbine-convention `a` — which
/// for a propeller (`U_a > V`) is negative and for hover (`V = 0`) is
/// `-inf`, so it never misfires on a thrusting rotor — and only when it
/// genuinely applies do we scale the momentum thrust by the ratio of the
/// corrected to the parabolic `CT`. The result is finite for any `phi`
/// and tends smoothly to zero at the tip where `F -> 0`.
fn momentum_thrust_per_span(rho: f64, r: f64, f_loss: f64, v_free: f64, u_a: f64, v_i: f64) -> f64 {
    // Simple annular-momentum thrust (exact for a thrusting disk).
    let simple = 4.0 * PI * rho * r * f_loss * u_a * v_i;

    // Turbine-convention axial induction a = 1 - U_a / V. Only defined for
    // a real freestream; a propeller / hover never enters the corrected
    // regime, so return the simple thrust there.
    if v_free <= 1e-9 {
        return simple;
    }
    // a is finite (v_free > 1e-9, u_a finite), so `<=` is safe here.
    let a = 1.0 - u_a / v_free;
    if a <= A_CRIT {
        return simple;
    }

    // Windmill-brake / turbulent-wake state: blend onto the Buhl
    // empirical CT line (continuous in value and slope with the parabola
    // at a = A_CRIT). Scale the simple thrust by CT_buhl(a)/CT_parab(a).
    let ct_parab = 4.0 * f_loss * a * (1.0 - a);
    let ct_buhl = 8.0 / 9.0 + (4.0 * f_loss - 40.0 / 9.0) * a + (50.0 / 9.0 - 4.0 * f_loss) * a * a;
    if ct_parab.abs() < 1e-12 {
        return simple;
    }
    let scale = (ct_buhl / ct_parab).clamp(0.0, 10.0);
    simple * scale
}

/// Tangential (swirl) induction factor `a'`, in the numerically BOUNDED
/// form. The balance `a'/(1+a') = sigma Ct / (4 F sin phi cos phi)`
/// rearranges to
///
/// `a' = sigma Ct / (4 F sin phi cos phi - sigma Ct)`,
///
/// which stays finite for any `phi` (the naive `kp/(1-kp)` form diverges
/// as the loading group passes 1).
///
/// This swirl factor is REPORT-ONLY: the thrust-balance solver and the
/// load integration use the converged velocity triangle (`W^2`) directly,
/// so `a'` does not feed back into `T`, `Q` or `P`. It is exposed only
/// for inspection. Where the Prandtl loss `F -> 0` (tip / root) the
/// element carries no load and the swirl is physically zero, so we report
/// `0` there rather than the spurious clamp value the bare algebra would
/// give; otherwise the result is bounded to keep `(1 + a')` positive.
fn tangential_induction(sigma: f64, ct: f64, f_loss: f64, sin_phi: f64, cos_phi: f64) -> f64 {
    // Negligible loss -> unloaded element -> no swirl.
    if f_loss < 1e-3 {
        return 0.0;
    }
    let num = sigma * ct;
    let den = 4.0 * f_loss * sin_phi * cos_phi - num;
    if den.abs() < 1e-12 {
        return 0.0;
    }
    let ap = num / den;
    if !ap.is_finite() {
        return 0.0;
    }
    // Keep (1 + a') strictly positive and bounded.
    ap.clamp(-0.99, 5.0)
}

/// Auxiliary converged quantities at a given `phi` (internal).
struct ElementState {
    phi: f64,
    alpha: f64,
    /// Through-disk axial speed `U_a = Omega r tan phi` (m/s).
    u_a: f64,
    /// Induced axial velocity `v_i = U_a - V` (m/s, propeller convention).
    v_i: f64,
    /// Relative speed squared `W^2 = U_a^2 + U_t^2` (m^2/s^2).
    w2: f64,
    f_loss: f64,
    cl: f64,
    cd: f64,
    cn: f64,
    ct: f64,
}

/// Solve one blade element for its inflow angle by bracketing bisection,
/// then assemble its [`ElementResult`].
fn solve_element(inp: ElementInput, station_index: usize) -> Result<ElementResult, RotorError> {
    let lo0 = ANGLE_EPS;
    let hi0 = PI / 2.0 - ANGLE_EPS;

    // The thrust-balance residual is continuous on (eps, pi/2 - eps) but
    // can be NON-monotonic (it may cross zero more than once: the high-
    // loading root at small phi, and spurious crossings at large phi where
    // the section is deeply stalled). Plain bisection on the full interval
    // would converge to an arbitrary one of these. We instead SCAN from
    // small phi upward and bracket the FIRST sign change — the physical
    // operating root for a thrusting rotor is the smallest-phi (highest-
    // loading) one. This makes bracketing deterministic and correct
    // regardless of monotonicity.
    let scan = 400usize;
    let mut lo = lo0;
    let mut hi = hi0;
    let mut r_lo;
    let (mut prev_res, prev_state) = residual(&inp, lo0);
    // An endpoint that is already (numerically) the root.
    if prev_res.abs() < 1e-14 {
        return Ok(assemble(&inp, prev_state));
    }
    let mut prev_phi = lo0;
    let mut found = false;
    r_lo = prev_res;
    for k in 1..=scan {
        let phi = lo0 + (hi0 - lo0) * (k as f64) / (scan as f64);
        let (res, st) = residual(&inp, phi);
        if res == 0.0 {
            // Exact hit on a scan node.
            return Ok(assemble(&inp, st));
        }
        if prev_res * res < 0.0 {
            lo = prev_phi;
            hi = phi;
            r_lo = prev_res;
            found = true;
            break;
        }
        prev_phi = phi;
        prev_res = res;
    }
    if !found {
        return Err(RotorError::NoConvergence {
            station: station_index,
            radius: inp.r,
            reason: "could not bracket the inflow-angle root",
        });
    }

    // Bisection within the bracketed sub-interval: guaranteed convergence,
    // capped iterations.
    let mut state = None;
    for _ in 0..MAX_ITERS {
        let mid = 0.5 * (lo + hi);
        let (r_mid, st) = residual(&inp, mid);
        if (hi - lo) < PHI_TOL || r_mid == 0.0 {
            state = Some(st);
            break;
        }
        if r_lo * r_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            r_lo = r_mid;
        }
        state = Some(st);
    }

    match state {
        Some(st) => Ok(assemble(&inp, st)),
        None => Err(RotorError::NoConvergence {
            station: station_index,
            radius: inp.r,
            reason: "bisection hit the iteration cap without converging",
        }),
    }
}

/// Assemble an [`ElementResult`] from a converged [`ElementState`], using
/// the converged relative speed to form `dT/dr` and `dQ/dr`.
fn assemble(inp: &ElementInput, st: ElementState) -> ElementResult {
    // Per-unit-span loads over all blades from the converged velocity
    // triangle (W^2 already reflects the converged inflow angle).
    let q_dyn = 0.5 * inp.rho * st.w2 * inp.blades * inp.chord;
    let dt_dr = q_dyn * st.cn;
    let dq_dr = q_dyn * st.ct * inp.r;

    // Report induction factors for inspection. Axial fraction
    // a = v_i / U_a (propeller convention); swirl a' from the bounded
    // torque/momentum balance.
    let sigma = inp.blades * inp.chord / (2.0 * PI * inp.r);
    let phi = st.phi;
    let axial_induction = if st.u_a.abs() < 1e-9 {
        0.0
    } else {
        st.v_i / st.u_a
    };
    let tangential_induction = tangential_induction(sigma, st.ct, st.f_loss, phi.sin(), phi.cos());

    ElementResult {
        radius_m: inp.r,
        inflow_angle_rad: st.phi,
        alpha_rad: st.alpha,
        axial_induction,
        tangential_induction,
        tip_loss: st.f_loss,
        cl: st.cl,
        cd: st.cd,
        dt_dr,
        dq_dr,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small 2-blade propeller: R = 0.15 m, hub at 0.15 R, a tapered,
    /// twisted blade with a sensible nose-down twist distribution.
    fn small_prop() -> Rotor {
        let r_tip = 0.15;
        let r_hub = 0.02;
        // Five stations from just outboard of the hub to the tip.
        let radii = [0.03, 0.06, 0.09, 0.12, 0.15];
        // Chord tapers from root to tip (m).
        let chords = [0.025, 0.022, 0.018, 0.014, 0.008];
        // Twist (rad): high pitch at root, washing out toward the tip —
        // a realistic propeller distribution (~25 deg -> ~8 deg).
        let twists = [
            25.0_f64.to_radians(),
            18.0_f64.to_radians(),
            13.0_f64.to_radians(),
            10.0_f64.to_radians(),
            8.0_f64.to_radians(),
        ];
        Rotor::from_slices(2, r_tip, r_hub, &radii, &chords, &twists).unwrap()
    }

    #[test]
    fn solve_returns_finite_converged_loads() {
        // Validation test 4: a normal forward-flight case converges to
        // finite T, Q, P.
        let rotor = small_prop();
        let perf = rotor.solve(5000.0, 5.0, 1.225).unwrap();
        assert!(perf.thrust_n.is_finite());
        assert!(perf.torque_nm.is_finite());
        assert!(perf.power_w.is_finite());
        assert!(perf.thrust_n > 0.0, "expected positive thrust");
        assert!(perf.power_w > 0.0, "expected positive shaft power");
        // Every element converged to finite state.
        for e in &perf.elements {
            assert!(e.inflow_angle_rad.is_finite() && e.inflow_angle_rad > 0.0);
            assert!(e.axial_induction.is_finite());
            assert!(e.tangential_induction.is_finite());
            assert!(e.tip_loss.is_finite() && e.tip_loss > 0.0 && e.tip_loss <= 1.0 + 1e-9);
        }
    }

    #[test]
    fn figure_of_merit_is_physical_in_hover() {
        // Validation test 1: 0 < FM <= 1 for a reasonable hovering rotor.
        // Ideal BEMT cannot beat the actuator-disk (momentum) limit; an
        // FM above 1 would mean the rotor achieves thrust below the
        // theoretical minimum induced power, which is impossible.
        let rotor = small_prop();
        let perf = rotor.solve(6000.0, 0.0, 1.225).unwrap();
        assert!(perf.thrust_n > 0.0, "hover should make positive thrust");
        assert!(perf.power_w > 0.0);
        let fm = perf.figure_of_merit;
        assert!(fm.is_finite(), "FM must be finite");
        assert!(fm > 0.0, "FM must be positive, got {fm}");
        assert!(
            fm <= 1.0,
            "FM must not exceed the actuator-disk limit, got {fm}"
        );
        // Tighter, model-grade sanity band: a real rotor's FM lands in
        // roughly 0.6-0.8; this idealised strip estimate should fall in a
        // plausible band, not be implausibly tiny or right at the limit.
        // (Documents the actual computed value ~0.65; see VALIDATION note
        // in the crate docs.)
        assert!(
            (0.4..0.95).contains(&fm),
            "hover FM out of the plausible 0.4-0.95 band: {fm}"
        );

        // Cross-check the limit directly from the power numbers: the
        // computed shaft power must be at least the ideal induced power
        // for the computed thrust (this is exactly FM <= 1, re-derived
        // independently of the stored field).
        let area = rotor.disk_area();
        let ideal_power = perf.thrust_n.powf(1.5) / (2.0 * 1.225 * area).sqrt();
        assert!(
            perf.power_w >= ideal_power - 1e-9,
            "shaft power {} below the actuator-disk ideal {}",
            perf.power_w,
            ideal_power
        );
    }

    #[test]
    fn thrust_increases_with_rpm() {
        // Validation test 2: thrust grows monotonically with rpm for a
        // fixed rotor at a fixed inflow.
        let rotor = small_prop();
        let v = 4.0;
        let rho = 1.225;
        let mut last = f64::NEG_INFINITY;
        for &rpm in &[2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0] {
            let t = rotor.solve(rpm, v, rho).unwrap().thrust_n;
            assert!(t.is_finite());
            assert!(
                t > last,
                "thrust must increase with rpm: {t} !> {last} at {rpm} rpm"
            );
            last = t;
        }
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        // Validation test 3: zero / negative / non-finite inputs -> Err,
        // never a NaN or a panic.
        let rotor = small_prop();
        assert!(rotor.solve(0.0, 5.0, 1.225).is_err()); // rpm = 0
        assert!(rotor.solve(-100.0, 5.0, 1.225).is_err()); // rpm < 0
        assert!(rotor.solve(5000.0, -1.0, 1.225).is_err()); // V < 0
        assert!(rotor.solve(5000.0, 5.0, 0.0).is_err()); // rho = 0
        assert!(rotor.solve(5000.0, 5.0, -1.0).is_err()); // rho < 0
        assert!(rotor.solve(f64::NAN, 5.0, 1.225).is_err()); // rpm NaN
        assert!(rotor.solve(5000.0, f64::INFINITY, 1.225).is_err()); // V inf
        assert!(rotor.solve(5000.0, 5.0, f64::NAN).is_err()); // rho NaN

        // Bad geometry at construction.
        let st = vec![
            BladeStation::new(0.05, 0.02, 0.2).unwrap(),
            BladeStation::new(0.10, 0.015, 0.15).unwrap(),
        ];
        assert!(matches!(
            Rotor::new(0, 0.15, 0.02, st.clone(), Polar::default()),
            Err(RotorError::NoBlades)
        ));
        assert!(Rotor::new(2, -0.15, 0.02, st.clone(), Polar::default()).is_err());
        assert!(Rotor::new(2, 0.15, 0.20, st.clone(), Polar::default()).is_err()); // hub > tip
        assert!(Rotor::new(2, 0.15, 0.02, vec![st[0]], Polar::default()).is_err());
        // 1 station
    }

    #[test]
    fn degenerate_hub_element_does_not_panic() {
        // Validation test 5: an element right at the hub (where the tip/hub
        // loss geometry and sin(phi) get extreme) must not panic. Put the
        // first station exactly at the hub radius and hover (V = 0, the
        // hardest case for the residual).
        let r_tip = 0.15;
        let r_hub = 0.03;
        let radii = [0.03, 0.075, 0.15]; // first station AT the hub
        let chords = [0.02, 0.016, 0.008];
        let twists = [
            20.0_f64.to_radians(),
            12.0_f64.to_radians(),
            6.0_f64.to_radians(),
        ];
        let rotor = Rotor::from_slices(2, r_tip, r_hub, &radii, &chords, &twists).unwrap();
        // Both hover and forward flight; just assert it returns finite
        // numbers (or a clean Err), never a panic / NaN.
        for &v in &[0.0, 3.0] {
            let perf = rotor.solve(5000.0, v, 1.225).unwrap();
            assert!(perf.thrust_n.is_finite());
            assert!(perf.power_w.is_finite());
            for e in &perf.elements {
                assert!(e.inflow_angle_rad.is_finite());
                assert!(e.tip_loss.is_finite());
            }
        }
    }

    #[test]
    fn higher_rpm_costs_more_power() {
        // Power should also rise with rpm (sanity beyond thrust).
        let rotor = small_prop();
        let p_lo = rotor.solve(3000.0, 4.0, 1.225).unwrap().power_w;
        let p_hi = rotor.solve(6000.0, 4.0, 1.225).unwrap().power_w;
        assert!(p_hi > p_lo, "power must grow with rpm: {p_hi} !> {p_lo}");
    }

    #[test]
    fn efficiency_zero_in_hover_positive_in_cruise() {
        let rotor = small_prop();
        let hover = rotor.solve(6000.0, 0.0, 1.225).unwrap();
        assert_eq!(hover.efficiency, 0.0, "eta is undefined/zero in hover");
        assert!(hover.figure_of_merit > 0.0, "hover defines a positive FM");

        let cruise = rotor.solve(6000.0, 8.0, 1.225).unwrap();
        assert!(cruise.efficiency.is_finite());
        assert!(cruise.efficiency > 0.0, "eta should be positive in cruise");
        // A real prop is well under 1; an ideal strip estimate should at
        // least be physical (0 < eta < 1).
        assert!(
            cruise.efficiency < 1.0,
            "eta should be < 1, got {}",
            cruise.efficiency
        );
        assert_eq!(cruise.figure_of_merit, 0.0, "FM defined only at hover");
    }

    #[test]
    fn denser_air_makes_more_thrust() {
        let rotor = small_prop();
        let sea = rotor.solve(5000.0, 5.0, 1.225).unwrap().thrust_n;
        let alt = rotor.solve(5000.0, 5.0, 0.9).unwrap().thrust_n;
        assert!(sea > alt, "denser air -> more thrust: {sea} !> {alt}");
    }

    #[test]
    fn performance_is_serde_round_trippable() {
        let perf = small_prop().solve(5000.0, 5.0, 1.225).unwrap();
        let json = serde_json::to_string(&perf).unwrap();
        let back: RotorPerformance = serde_json::from_str(&json).unwrap();
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-6 * b.abs().max(1.0);
        assert!(close(back.thrust_n, perf.thrust_n));
        assert!(close(back.torque_nm, perf.torque_nm));
        assert!(close(back.power_w, perf.power_w));
        assert_eq!(back.elements.len(), perf.elements.len());
    }

    #[test]
    fn tabulated_polar_runs() {
        // A coarse symmetric tabulated polar should also produce a finite
        // converged solve, exercising the Table interpolation path.
        let samples: Vec<_> = (-15..=15)
            .step_by(3)
            .map(|deg| {
                let a = (deg as f64).to_radians();
                crate::airfoil::PolarSample {
                    alpha_rad: a,
                    cl: 2.0 * PI * a.sin(),
                    cd: 0.01 + 0.02 * (2.0 * PI * a.sin()).powi(2),
                }
            })
            .collect();
        let polar = Polar::table(samples).unwrap();
        let radii = [0.03, 0.06, 0.09, 0.12, 0.15];
        let chords = [0.025, 0.02, 0.016, 0.012, 0.008];
        let twists = [
            14.0_f64.to_radians(),
            11.0_f64.to_radians(),
            9.0_f64.to_radians(),
            7.0_f64.to_radians(),
            6.0_f64.to_radians(),
        ];
        let mut stations = Vec::new();
        for i in 0..radii.len() {
            stations.push(BladeStation::new(radii[i], chords[i], twists[i]).unwrap());
        }
        let rotor = Rotor::new(2, 0.15, 0.02, stations, polar).unwrap();
        let perf = rotor.solve(5000.0, 5.0, 1.225).unwrap();
        assert!(perf.thrust_n.is_finite() && perf.power_w.is_finite());
    }
}

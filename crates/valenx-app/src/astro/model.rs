//! State, conversions and builders for the Astro / Launch workbench.
//!
//! Everything in this module is **pure, non-UI logic** — the form state
//! structs the panels mutate, the [`valenx_astro`] builders that turn
//! those forms into a [`Vehicle`] + [`AscentConfig`], the unit
//! conversions, and the result formatters. Splitting it out keeps the
//! [`super::panels`] module purely about egui layout and makes the
//! workbench's logic `#[test]`-coverable without standing up an egui
//! context (mirrors the CFD-side [`crate::aero::model`]).

use valenx_astro::config::{AscentConfig, GuidanceMode};
use valenx_astro::guidance::GuidanceProgram;
use valenx_astro::vehicle::{DragModel, Stage, Vehicle};
use valenx_astro::wind::WindModel;
use valenx_astro::AstroError;

// ---------------------------------------------------------------------------
// Enumerations the form exposes
// ---------------------------------------------------------------------------

/// Which workbench sub-view (tab) is selected. The ascent simulator and
/// the closed-form mission planners are the two halves of the panel.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum AstroTab {
    /// The full ascent-to-orbit simulation (vehicle setup + Run + the
    /// trajectory result + flight-profile chart).
    #[default]
    Ascent,
    /// The closed-form mission planners (Hohmann, hoverslam, rendezvous,
    /// launch azimuth) — each an input → output card.
    Planners,
}

impl AstroTab {
    /// All tabs in display order.
    pub const ALL: [AstroTab; 2] = [AstroTab::Ascent, AstroTab::Planners];

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            AstroTab::Ascent => "Ascent to orbit",
            AstroTab::Planners => "Mission planners",
        }
    }
}

/// The steering / cutoff strategy the form offers, mirroring
/// [`GuidanceMode`] without its payload (the target altitude is a
/// separate numeric field so the combo box stays a plain enum).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum GuidanceChoice {
    /// Open-loop gravity turn — burn every stage to depletion.
    #[default]
    OpenLoopGravityTurn,
    /// Closed-loop insertion to the target circular altitude.
    ClosedLoopInsertion,
}

impl GuidanceChoice {
    /// All choices in display order.
    pub const ALL: [GuidanceChoice; 2] = [
        GuidanceChoice::OpenLoopGravityTurn,
        GuidanceChoice::ClosedLoopInsertion,
    ];

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            GuidanceChoice::OpenLoopGravityTurn => "Open-loop gravity turn",
            GuidanceChoice::ClosedLoopInsertion => "Closed-loop insertion",
        }
    }
}

// ---------------------------------------------------------------------------
// Form state
// ---------------------------------------------------------------------------

/// One propulsive stage as the form edits it — the four numbers a user
/// types per stage. Converted to a [`Stage`] by [`StageForm::to_stage`].
#[derive(Clone, Debug, PartialEq)]
pub struct StageForm {
    /// Stage display name.
    pub name: String,
    /// Inert / structural mass jettisoned at staging (kg).
    pub dry_mass: f64,
    /// Usable propellant mass (kg).
    pub propellant_mass: f64,
    /// Vacuum thrust (N).
    pub thrust_vac: f64,
    /// Sea-level thrust (N).
    pub thrust_sl: f64,
    /// Vacuum specific impulse (s).
    pub isp_vac: f64,
    /// Sea-level specific impulse (s).
    pub isp_sl: f64,
}

impl StageForm {
    /// Convert to a backend [`Stage`].
    pub fn to_stage(&self) -> Stage {
        Stage {
            name: self.name.clone(),
            dry_mass: self.dry_mass,
            propellant_mass: self.propellant_mass,
            thrust_vac: self.thrust_vac,
            thrust_sl: self.thrust_sl,
            isp_vac: self.isp_vac,
            isp_sl: self.isp_sl,
        }
    }
}

/// Every Astro / Launch form input.
///
/// The defaults reproduce the `valenx-astro`
/// [`presets::two_stage_medium_lift`](valenx_astro::presets::two_stage_medium_lift)
/// vehicle and a LEO ascent so the panel is runnable the moment it
/// opens — exactly as the CFD workbench's demo box is.
#[derive(Clone, Debug, PartialEq)]
pub struct AscentForm {
    /// The stage stack (index 0 fires first).
    pub stages: Vec<StageForm>,
    /// Payload mass carried to orbit (kg).
    pub payload_mass: f64,
    /// Aerodynamic reference (frontal) area (m²).
    pub reference_area: f64,
    /// Target circular-orbit altitude (km) — used as the closed-loop
    /// insertion target and echoed in the summary.
    pub target_altitude_km: f64,
    /// Launch-site elevation above the equatorial radius (m).
    pub launch_altitude_m: f64,
    /// Steering / cutoff strategy.
    pub guidance: GuidanceChoice,
    /// Vertical-rise hold time (s).
    pub vertical_rise_time: f64,
    /// Pitch-kick angle off vertical (deg).
    pub pitch_kick_deg: f64,
    /// Pitch-kick hold duration (s).
    pub kick_duration: f64,
    /// Steady horizontal wind speed applied to drag (m/s); `0` = calm.
    pub wind_speed_ms: f64,
}

impl Default for AscentForm {
    fn default() -> Self {
        // Mirror presets::two_stage_medium_lift + leo_ascent_config so a
        // fresh panel runs to orbit on the first click.
        Self {
            stages: vec![
                StageForm {
                    name: "first stage".to_string(),
                    dry_mass: 25_000.0,
                    propellant_mass: 410_000.0,
                    thrust_vac: 8_200_000.0,
                    thrust_sl: 7_600_000.0,
                    isp_vac: 311.0,
                    isp_sl: 283.0,
                },
                StageForm {
                    name: "second stage".to_string(),
                    dry_mass: 4_000.0,
                    propellant_mass: 100_000.0,
                    thrust_vac: 980_000.0,
                    thrust_sl: 980_000.0,
                    isp_vac: 348.0,
                    isp_sl: 348.0,
                },
            ],
            payload_mass: 10_000.0,
            reference_area: 10.75,
            target_altitude_km: 300.0,
            launch_altitude_m: 0.0,
            guidance: GuidanceChoice::OpenLoopGravityTurn,
            vertical_rise_time: 20.0,
            pitch_kick_deg: 12.0,
            kick_duration: 5.0,
            wind_speed_ms: 0.0,
        }
    }
}

impl AscentForm {
    /// Build the backend [`Vehicle`] from the form. Pure value
    /// conversion — validation happens inside `valenx-astro` when the
    /// vehicle is simulated, so a non-physical entry surfaces as a
    /// run-time `Err` in the UI rather than a panic here.
    pub fn build_vehicle(&self) -> Vehicle {
        Vehicle {
            stages: self.stages.iter().map(StageForm::to_stage).collect(),
            payload_mass: self.payload_mass,
            reference_area: self.reference_area,
            drag: DragModel::generic_launch_vehicle(),
        }
    }

    /// Build the backend [`AscentConfig`] from the form.
    ///
    /// The integration step / time cap / sample interval mirror the
    /// crate presets (a bounded run); the guidance mode + program come
    /// from the form. The closed-loop target is the `target_altitude_km`
    /// field converted to metres.
    pub fn build_config(&self) -> AscentConfig {
        let mode = match self.guidance {
            GuidanceChoice::OpenLoopGravityTurn => GuidanceMode::OpenLoopGravityTurn,
            GuidanceChoice::ClosedLoopInsertion => GuidanceMode::ClosedLoopInsertion {
                target_altitude_m: self.target_altitude_km * 1_000.0,
            },
        };
        // Closed-loop insertion needs time for the coast + circularise;
        // the open-loop gravity turn finishes sooner. Match the presets.
        let max_time = match self.guidance {
            GuidanceChoice::OpenLoopGravityTurn => 1_500.0,
            GuidanceChoice::ClosedLoopInsertion => 3_000.0,
        };
        let wind = if self.wind_speed_ms.abs() > f64::EPSILON {
            WindModel::Constant(self.wind_speed_ms)
        } else {
            WindModel::None
        };
        AscentConfig {
            launch_altitude_m: self.launch_altitude_m,
            guidance: GuidanceProgram {
                vertical_rise_time: self.vertical_rise_time,
                pitch_kick_deg: self.pitch_kick_deg,
                kick_duration: self.kick_duration,
            },
            time_step: 0.1,
            max_time,
            sample_interval: 2.0,
            mode,
            wind,
        }
    }
}

/// All the inputs the closed-form mission planners take. Each planner is
/// independent; they share this one form so the panel can keep them in
/// adjacent collapsing headers.
#[derive(Clone, Debug, PartialEq)]
pub struct PlannerForm {
    /// Hohmann — departure circular-orbit altitude (km).
    pub hohmann_from_km: f64,
    /// Hohmann — arrival circular-orbit altitude (km).
    pub hohmann_to_km: f64,

    /// Hoverslam — descent speed at the start of the burn (m/s).
    pub hoverslam_descent_speed: f64,
    /// Hoverslam — net deceleration `T/m − g` (m/s²).
    pub hoverslam_net_decel: f64,

    /// Rendezvous — reference circular-orbit altitude (km) of the
    /// target, which sets the mean motion `n`.
    pub rdv_orbit_altitude_km: f64,
    /// Rendezvous — initial chaser radial offset in LVLH (m).
    pub rdv_offset_radial: f64,
    /// Rendezvous — initial chaser along-track offset in LVLH (m).
    pub rdv_offset_along: f64,
    /// Rendezvous — transfer time as a fraction of the orbital period.
    pub rdv_transfer_fraction: f64,

    /// Launch azimuth — launch-site latitude (deg).
    pub azimuth_latitude_deg: f64,
    /// Launch azimuth — target orbit inclination (deg).
    pub azimuth_inclination_deg: f64,

    /// Plane change — circular-orbit altitude (km) where the burn happens.
    pub plane_change_altitude_km: f64,
    /// Plane change — inclination change Δi (deg).
    pub plane_change_delta_inc_deg: f64,

    /// Orbit basics — circular-orbit altitude (km).
    pub basics_altitude_km: f64,

    /// Bi-elliptic — departure circular-orbit altitude (km).
    pub bielliptic_from_km: f64,
    /// Bi-elliptic — arrival circular-orbit altitude (km).
    pub bielliptic_to_km: f64,
    /// Bi-elliptic — intermediate apoapsis altitude (km) the transfer
    /// coasts out to before dropping onto the arrival orbit.
    pub bielliptic_via_km: f64,

    /// Elliptical orbit — perigee altitude (km).
    pub ellipse_perigee_km: f64,
    /// Elliptical orbit — apogee altitude (km).
    pub ellipse_apogee_km: f64,

    /// Synodic period — first circular-orbit altitude (km).
    pub synodic_a_km: f64,
    /// Synodic period — second circular-orbit altitude (km).
    pub synodic_b_km: f64,

    /// Target-period orbit — desired orbital period (hours).
    pub target_period_h: f64,

    /// Injection Δv — parking-orbit altitude (km).
    pub injection_altitude_km: f64,
    /// Injection Δv — target hyperbolic excess speed v∞ (km/s).
    pub injection_v_inf_kms: f64,

    /// Flight-path angle — orbit eccentricity (0–1).
    pub fpa_eccentricity: f64,
    /// Flight-path angle — true anomaly θ (degrees).
    pub fpa_true_anomaly_deg: f64,

    /// Orbital speed at altitude — perigee altitude (km).
    pub speed_perigee_km: f64,
    /// Orbital speed at altitude — apogee altitude (km).
    pub speed_apogee_km: f64,
    /// Orbital speed at altitude — query altitude (km).
    pub speed_query_km: f64,

    /// Time of flight — perigee altitude (km).
    pub tof_perigee_km: f64,
    /// Time of flight — apogee altitude (km).
    pub tof_apogee_km: f64,
    /// Time of flight — true anomaly θ (degrees).
    pub tof_true_anomaly_deg: f64,

    /// Horizon & coverage — satellite altitude (km).
    pub horizon_altitude_km: f64,
}

impl Default for PlannerForm {
    fn default() -> Self {
        Self {
            // LEO -> GEO, the canonical textbook transfer.
            hohmann_from_km: 300.0,
            hohmann_to_km: 35_786.0,
            // A booster decelerating at ~2 g net from 250 m/s.
            hoverslam_descent_speed: 250.0,
            hoverslam_net_decel: 20.0,
            // A ~400 km LEO rendezvous, chaser a few km out, quarter
            // period transfer (sin(nT) = 1, well clear of the singular
            // half-period family).
            rdv_orbit_altitude_km: 400.0,
            rdv_offset_radial: -2_000.0,
            rdv_offset_along: 5_000.0,
            rdv_transfer_fraction: 0.25,
            // KSC latitude -> the ISS inclination.
            azimuth_latitude_deg: 28.5,
            azimuth_inclination_deg: 51.6,
            // A 28.5° plane change at LEO — eye-wateringly expensive.
            plane_change_altitude_km: 400.0,
            plane_change_delta_inc_deg: 28.5,
            basics_altitude_km: 400.0,
            // LEO -> GEO routed via a far apoapsis — the textbook 3-burn case.
            bielliptic_from_km: 300.0,
            bielliptic_to_km: 35_786.0,
            bielliptic_via_km: 250_000.0,
            // A geostationary transfer orbit — the classic eccentric orbit.
            ellipse_perigee_km: 300.0,
            ellipse_apogee_km: 35_786.0,
            // ISS-altitude LEO vs GEO — a wide period spread.
            synodic_a_km: 400.0,
            synodic_b_km: 35_786.0,
            // One sidereal day → the geostationary altitude.
            target_period_h: 23.934,
            // A 300 km LEO with a trans-Mars-class hyperbolic excess.
            injection_altitude_km: 300.0,
            injection_v_inf_kms: 3.0,
            // A moderately eccentric orbit, sampled a third of the way to apogee.
            fpa_eccentricity: 0.5,
            fpa_true_anomaly_deg: 60.0,
            // A GTO queried partway up its climb to apogee.
            speed_perigee_km: 300.0,
            speed_apogee_km: 35_786.0,
            speed_query_km: 2_000.0,
            // A GTO sampled a quarter-turn past perigee.
            tof_perigee_km: 300.0,
            tof_apogee_km: 35_786.0,
            tof_true_anomaly_deg: 90.0,
            // A typical LEO altitude.
            horizon_altitude_km: 400.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Formatting helpers (pure)
// ---------------------------------------------------------------------------

/// Friendly one-line text for an [`AstroError`] surfaced from a run /
/// planner. Pure mapping — never panics, always returns a string.
pub fn friendly_error(err: &AstroError) -> String {
    format!("{err}")
}

/// Format a Δv in m/s with a km/s echo for the large numbers.
pub fn format_delta_v(dv: f64) -> String {
    if dv.abs() >= 1_000.0 {
        format!("{dv:.0} m/s ({:.2} km/s)", dv / 1_000.0)
    } else {
        format!("{dv:.1} m/s")
    }
}

/// Format a duration in seconds as the most readable of s / min / h.
pub fn format_duration(secs: f64) -> String {
    if !secs.is_finite() {
        return "—".to_string();
    }
    if secs >= 3_600.0 {
        format!("{:.2} h", secs / 3_600.0)
    } else if secs >= 120.0 {
        format!("{:.1} min", secs / 60.0)
    } else {
        format!("{secs:.1} s")
    }
}

/// Convert an orbit altitude (km) above the equator to a geocentric
/// radius (m) — the input the maneuver / rendezvous planners take.
pub fn altitude_km_to_radius_m(altitude_km: f64) -> f64 {
    valenx_astro::constants::R_EARTH + altitude_km * 1_000.0
}

/// The orbital speed (m/s) at altitude `query_altitude_km` on the elliptical
/// Earth orbit defined by the given perigee/apogee altitudes (km), from vis-viva
/// `v = √(μ(2/r − 1/a))` with `a = (r_peri + r_apo)/2`. The perigee/apogee inputs
/// may be given in either order. Returns `None` if the query altitude lies
/// outside the orbit's `[perigee, apogee]` span — the orbit never reaches it.
pub fn orbital_speed_at_altitude(
    perigee_km: f64,
    apogee_km: f64,
    query_altitude_km: f64,
) -> Option<f64> {
    let r_peri = altitude_km_to_radius_m(perigee_km.min(apogee_km));
    let r_apo = altitude_km_to_radius_m(perigee_km.max(apogee_km));
    let r = altitude_km_to_radius_m(query_altitude_km);
    // A 1 m tolerance keeps the apsides themselves in range against float error.
    if r < r_peri - 1.0 || r > r_apo + 1.0 {
        return None;
    }
    let mu = valenx_astro::constants::MU_EARTH;
    let a = 0.5 * (r_peri + r_apo);
    Some((mu * (2.0 / r - 1.0 / a)).sqrt())
}

/// The time since perigee (s) at true anomaly `true_anomaly_deg` on the
/// elliptical Earth orbit with the given perigee/apogee altitudes (km), via
/// Kepler's equation: eccentric anomaly
/// `E = 2·atan2(√(1−e)·sin(θ/2), √(1+e)·cos(θ/2))`, mean anomaly
/// `M = E − e·sin E`, and `t = M·T/2π`. The true anomaly wraps to one orbit, so
/// the result is `0` at perigee (θ = 0) and half the period at apogee
/// (θ = 180°) regardless of eccentricity.
pub fn time_since_perigee(perigee_km: f64, apogee_km: f64, true_anomaly_deg: f64) -> f64 {
    let orbit = elliptical_orbit(perigee_km, apogee_km);
    let e = orbit.eccentricity;
    let tau = std::f64::consts::TAU;
    let theta = true_anomaly_deg.to_radians().rem_euclid(tau);
    let (sin_half, cos_half) = (0.5 * theta).sin_cos();
    let ecc_anom = 2.0 * ((1.0 - e).sqrt() * sin_half).atan2((1.0 + e).sqrt() * cos_half);
    let ecc_anom = ecc_anom.rem_euclid(tau);
    let mean_anom = ecc_anom - e * ecc_anom.sin();
    mean_anom * orbit.period_s / tau
}

/// The line-of-sight slant range to the horizon (km) from a satellite at
/// altitude `altitude_km` — `√(h² + 2·R⊕·h)`, the farthest ground point within
/// line of sight (the geometric tangent range). Negative altitudes clamp to 0.
pub fn horizon_slant_range_km(altitude_km: f64) -> f64 {
    let h = (altitude_km * 1000.0).max(0.0); // m
    let r = valenx_astro::constants::R_EARTH;
    (h * h + 2.0 * r * h).sqrt() / 1000.0 // → km
}

/// The Earth-central half-angle (degrees) of a satellite's visibility footprint
/// at altitude `altitude_km` — `acos(R⊕/(R⊕+h))`, the angular radius of the
/// ground area within line of sight (twice this spans the full footprint). `0`
/// at the surface, approaching 90° as the orbit grows.
pub fn coverage_half_angle_deg(altitude_km: f64) -> f64 {
    let h = (altitude_km * 1000.0).max(0.0);
    let r = valenx_astro::constants::R_EARTH;
    (r / (r + h)).acos().to_degrees()
}

/// Δv (m/s) for a pure inclination change of `delta_inc_deg` on a circular
/// orbit at `altitude_km` — a units wrapper (km, deg) over
/// [`valenx_astro::maneuver::circular_plane_change_dv`] (`Δv = 2·v·sin(Δi/2)`).
pub fn plane_change_dv(altitude_km: f64, delta_inc_deg: f64) -> Result<f64, AstroError> {
    let r = altitude_km_to_radius_m(altitude_km);
    valenx_astro::maneuver::circular_plane_change_dv(r, delta_inc_deg.to_radians())
}

/// The Hohmann departure **phase angle** (radians, wrapped to (−π, π]) between
/// circular orbits at altitudes `from_km` and `to_km` — the angular lead the
/// destination must have over the spacecraft at the first burn so the two reach
/// the arrival radius together: `φ = π·(1 − ((r1 + r2)/(2·r2))^{3/2})`. Positive
/// means the target leads (an outward transfer); negative means it trails (an
/// inward transfer); zero for equal orbits.
pub fn hohmann_phase_angle(from_km: f64, to_km: f64) -> f64 {
    let r1 = altitude_km_to_radius_m(from_km);
    let r2 = altitude_km_to_radius_m(to_km);
    let a_t = 0.5 * (r1 + r2);
    let raw = std::f64::consts::PI * (1.0 - (a_t / r2).powf(1.5));
    // Wrap into (−π, π] so an inward transfer's large lead reads as a lead/lag.
    (raw + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
}

/// Circular-orbit speed, escape speed (m/s) and orbital period (s) at altitude
/// `altitude_km`, from `r = R_⊕ + h` and Earth's μ:
/// `v_circ = √(μ/r)`, `v_esc = √(2μ/r)`, `T = 2π·√(r³/μ)`.
pub fn circular_orbit_basics(altitude_km: f64) -> (f64, f64, f64) {
    let r = altitude_km_to_radius_m(altitude_km);
    let mu = valenx_astro::constants::MU_EARTH;
    let v_circ = (mu / r).sqrt();
    let v_esc = (2.0 * mu / r).sqrt();
    let period = std::f64::consts::TAU * (r * r * r / mu).sqrt();
    (v_circ, v_esc, period)
}

/// The specific orbital energy (J/kg) and specific angular momentum (m²/s) of a circular
/// orbit at altitude `altitude_km` — `ε = −μ/(2a)`
/// ([`valenx_astro::orbit::specific_orbital_energy`]) and `h = √(μ·r)`
/// ([`valenx_astro::orbit::circular_angular_momentum`]), with `a = r` for a circle. Returns
/// `None` for a non-physical altitude (radius non-positive).
pub fn circular_orbit_energy_momentum(altitude_km: f64) -> Option<(f64, f64)> {
    let r = altitude_km_to_radius_m(altitude_km);
    let energy = valenx_astro::orbit::specific_orbital_energy(r).ok()?;
    let ang_mom = valenx_astro::orbit::circular_angular_momentum(r).ok()?;
    Some((energy, ang_mom))
}

/// The prograde Δv (m/s) to escape from a circular orbit at altitude
/// `altitude_km` — the burn from circular speed up to escape speed,
/// `Δv = v_esc − v_circ = (√2 − 1)·v_circ ≈ 0.414·v_circ`.
pub fn escape_delta_v_from_circular(altitude_km: f64) -> f64 {
    let (v_circ, v_esc, _) = circular_orbit_basics(altitude_km);
    v_esc - v_circ
}

/// The circular-orbit altitude (km) whose period equals `period_s` — the
/// inverse of [`circular_orbit_basics`]'s period via Kepler's third law,
/// `r = (μ·T²/4π²)^(1/3)`, then altitude `= r − R_⊕`. One sidereal day
/// (86 164 s) gives the geostationary altitude (≈ 35 786 km); half a sidereal
/// day (≈ 43 082 s) the GPS orbit (≈ 20 200 km). Negative for a period too
/// short to clear the surface.
pub fn orbit_altitude_for_period_km(period_s: f64) -> f64 {
    let mu = valenx_astro::constants::MU_EARTH;
    let four_pi_sq = 4.0 * std::f64::consts::PI * std::f64::consts::PI;
    let r = (mu * period_s * period_s / four_pi_sq).cbrt();
    (r - valenx_astro::constants::R_EARTH) / 1000.0
}

/// The departure Δv (m/s) from a circular parking orbit at altitude
/// `altitude_km` onto a hyperbolic escape trajectory with hyperbolic excess
/// speed `v_inf_ms` (m/s): `Δv = √(v_esc² + v_inf²) − v_circ`. With `v_inf = 0`
/// it reduces to the bare escape Δv; a trans-Mars `v_inf ≈ 3 km/s` from LEO
/// gives the classic ~3.6 km/s trans-Mars-injection burn.
pub fn injection_delta_v(altitude_km: f64, v_inf_ms: f64) -> f64 {
    let (v_circ, v_esc, _) = circular_orbit_basics(altitude_km);
    (v_esc * v_esc + v_inf_ms * v_inf_ms).sqrt() - v_circ
}

/// The flight-path angle `γ = atan( e·sinθ / (1 + e·cosθ) )` (radians) at true
/// anomaly `theta_rad` on an orbit of eccentricity `eccentricity` — the angle of
/// the velocity above the local horizontal (perpendicular to the radius). It is
/// zero at both apsides (`θ = 0, π`, where the motion is purely tangential) and
/// everywhere on a circular orbit (`e = 0`); positive while climbing from
/// perigee toward apogee and negative while descending.
pub fn flight_path_angle(eccentricity: f64, theta_rad: f64) -> f64 {
    let (sin_t, cos_t) = theta_rad.sin_cos();
    (eccentricity * sin_t).atan2(1.0 + eccentricity * cos_t)
}

/// The three-burn Δv budget and timing of a **bi-elliptic** transfer between
/// circular orbits at altitudes `from_km` and `to_km`, routed via an
/// intermediate apoapsis at altitude `via_km` — a units wrapper (km) over
/// [`valenx_astro::bielliptic_transfer`]. For large radius ratios a high
/// enough `via_km` beats the Hohmann total Δv, at the cost of a far longer
/// flight time.
pub fn bielliptic_transfer_altitudes(
    from_km: f64,
    to_km: f64,
    via_km: f64,
) -> Result<valenx_astro::Transfer, AstroError> {
    let r1 = altitude_km_to_radius_m(from_km);
    let r2 = altitude_km_to_radius_m(to_km);
    let r_b = altitude_km_to_radius_m(via_km);
    valenx_astro::bielliptic_transfer(r1, r2, r_b)
}

/// A summary of an elliptical Earth orbit specified by its perigee and apogee.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticalOrbit {
    /// Semi-major axis `a` (m).
    pub semi_major_axis_m: f64,
    /// Eccentricity `e` (0 = circular).
    pub eccentricity: f64,
    /// Orbital period `T` (s).
    pub period_s: f64,
    /// Speed at perigee (m/s) — the orbit's fastest point.
    pub perigee_speed_ms: f64,
    /// Speed at apogee (m/s) — the orbit's slowest point.
    pub apogee_speed_ms: f64,
    /// Specific orbital energy `ε = −μ/(2a)` (J/kg) — the conserved total energy
    /// per unit mass (kinetic + gravitational potential). Negative for a bound
    /// ellipse; it rises toward zero as the orbit nears escape (`a → ∞`).
    pub specific_energy_j_per_kg: f64,
    /// Specific angular momentum `h = √(μ·a·(1−e²))` (m²/s) — conserved along the
    /// orbit (Kepler's 2nd law: the areal sweep rate is `h/2`). Equals the
    /// product `r·v` at both apsides.
    pub specific_angular_momentum_m2_s: f64,
    /// Semi-latus rectum `p = a·(1−e²)` (m) — the conic parameter in the orbit
    /// equation `r = p/(1 + e·cosθ)`; the orbital radius at ±90° true anomaly, and
    /// the harmonic mean of the apsis radii.
    pub semi_latus_rectum_m: f64,
}

/// Summarize the elliptical Earth orbit with the two given altitudes (km),
/// taking the lower as perigee (so the inputs may be in either order). Uses
/// vis-viva `v = √(μ(2/r − 1/a))` and Kepler's third law `T = 2π√(a³/μ)` with
/// Earth's μ, where `a = (r_peri + r_apo)/2` and
/// `e = (r_apo − r_peri)/(r_apo + r_peri)`. Equal altitudes give a circle
/// (`e = 0`, equal peri/apogee speeds).
pub fn elliptical_orbit(altitude_a_km: f64, altitude_b_km: f64) -> EllipticalOrbit {
    let ra = altitude_km_to_radius_m(altitude_a_km);
    let rb = altitude_km_to_radius_m(altitude_b_km);
    let r_peri = ra.min(rb);
    let r_apo = ra.max(rb);
    let mu = valenx_astro::constants::MU_EARTH;
    let a = 0.5 * (r_peri + r_apo);
    let eccentricity = (r_apo - r_peri) / (r_apo + r_peri);
    let period_s = std::f64::consts::TAU * (a * a * a / mu).sqrt();
    let perigee_speed_ms = (mu * (2.0 / r_peri - 1.0 / a)).sqrt();
    let apogee_speed_ms = (mu * (2.0 / r_apo - 1.0 / a)).sqrt();
    let specific_energy_j_per_kg = -mu / (2.0 * a);
    let specific_angular_momentum_m2_s = (mu * a * (1.0 - eccentricity * eccentricity)).sqrt();
    let semi_latus_rectum_m =
        valenx_astro::orbit::semi_latus_rectum_from_elements(a, eccentricity).unwrap_or(0.0);
    EllipticalOrbit {
        semi_major_axis_m: a,
        eccentricity,
        period_s,
        perigee_speed_ms,
        apogee_speed_ms,
        specific_energy_j_per_kg,
        specific_angular_momentum_m2_s,
        semi_latus_rectum_m,
    }
}

/// Synodic period (s) between two circular Earth orbits at altitudes
/// `alt_a_km` and `alt_b_km` — the time between successive
/// same-relative-geometry alignments, `T_syn = 1 / |1/T_a − 1/T_b|`, which
/// sets how often a transfer or rendezvous window recurs. Returns
/// `f64::INFINITY` when the two periods are equal (the orbits never drift
/// apart).
pub fn synodic_period(alt_a_km: f64, alt_b_km: f64) -> f64 {
    let (_, _, t_a) = circular_orbit_basics(alt_a_km);
    let (_, _, t_b) = circular_orbit_basics(alt_b_km);
    let rate = (1.0 / t_a - 1.0 / t_b).abs();
    if rate > 0.0 {
        1.0 / rate
    } else {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_astro::simulate_ascent;

    #[test]
    fn default_form_reproduces_the_medium_lift_preset() {
        // The default vehicle + config must be a runnable case that
        // reaches orbit, so the panel is useful the instant it opens.
        let form = AscentForm::default();
        let vehicle = form.build_vehicle();
        let config = form.build_config();
        vehicle.validate().expect("default vehicle is valid");
        config.validate().expect("default config is valid");
        let r = simulate_ascent(&vehicle, &config).expect("default case runs");
        assert!(r.reached_space, "apoapsis only {:.1} km", r.apoapsis_km());
        assert!(r.ideal_delta_v > 9_000.0, "Δv budget {}", r.ideal_delta_v);
    }

    #[test]
    fn closed_loop_choice_sets_the_target_altitude() {
        let form = AscentForm {
            guidance: GuidanceChoice::ClosedLoopInsertion,
            target_altitude_km: 250.0,
            ..Default::default()
        };
        let config = form.build_config();
        match config.mode {
            GuidanceMode::ClosedLoopInsertion { target_altitude_m } => {
                assert!((target_altitude_m - 250_000.0).abs() < 1e-6);
            }
            other => panic!("expected closed-loop, got {other:?}"),
        }
        // Closed-loop gets the longer time budget for the coast + burn.
        assert!(config.max_time >= 3_000.0);
    }

    #[test]
    fn wind_speed_maps_to_a_constant_wind_model() {
        let mut form = AscentForm::default();
        assert_eq!(form.build_config().wind, WindModel::None);
        form.wind_speed_ms = 120.0;
        assert_eq!(form.build_config().wind, WindModel::Constant(120.0));
    }

    #[test]
    fn altitude_to_radius_adds_earth_radius() {
        let r = altitude_km_to_radius_m(300.0);
        assert!((r - (valenx_astro::constants::R_EARTH + 300_000.0)).abs() < 1e-6);
    }

    #[test]
    fn plane_change_dv_wraps_the_backend() {
        // Zero inclination change is free; Δv grows with the angle.
        assert_eq!(plane_change_dv(400.0, 0.0).unwrap(), 0.0);
        let small = plane_change_dv(400.0, 15.0).unwrap();
        let large = plane_change_dv(400.0, 45.0).unwrap();
        assert!(large > small && small > 0.0, "Δv grows with Δi: {small} → {large}");
        // A non-physical altitude (below the Earth's centre) errors, not panics.
        assert!(plane_change_dv(-10_000.0, 30.0).is_err());
    }

    #[test]
    fn circular_orbit_basics_match_textbook_values() {
        // At Earth's radius: ~7.9 km/s circular, escape = √2× that, ~84 min.
        let (v_circ, v_esc, period) = circular_orbit_basics(0.0);
        assert!((v_circ - 7910.0).abs() < 120.0, "v_circ {v_circ}");
        assert!((v_esc - 2.0_f64.sqrt() * v_circ).abs() < 1e-6, "v_esc = √2·v_circ");
        assert!((period / 60.0 - 84.4).abs() < 2.0, "T {} min", period / 60.0);
        // 300 km LEO: ~7.73 km/s, ~90 min period.
        let (v300, _, t300) = circular_orbit_basics(300.0);
        assert!((v300 - 7730.0).abs() < 120.0, "v300 {v300}");
        assert!((t300 / 60.0 - 90.5).abs() < 3.0, "T300 {} min", t300 / 60.0);
    }

    #[test]
    fn bielliptic_transfer_altitudes_is_a_three_burn_budget() {
        // LEO 300 km -> GEO 35 786 km, routed via a far 250 000 km apoapsis.
        let t = bielliptic_transfer_altitudes(300.0, 35_786.0, 250_000.0)
            .expect("valid altitudes");
        // Three finite, non-negative burns that sum to the reported total.
        assert!(t.delta_v1 >= 0.0 && t.delta_v2 >= 0.0 && t.delta_v3 >= 0.0);
        assert!(t.total_delta_v.is_finite());
        assert!((t.total_delta_v - (t.delta_v1 + t.delta_v2 + t.delta_v3)).abs() < 1e-6);
        // Unlike a two-burn Hohmann, the bi-elliptic uses a real third burn.
        assert!(t.delta_v3 > 0.0, "the third (circularising) burn is non-zero");
        // Coasting out past the target and back takes longer than the direct
        // Hohmann between the same two orbits.
        let h = valenx_astro::hohmann_transfer(
            altitude_km_to_radius_m(300.0),
            altitude_km_to_radius_m(35_786.0),
        )
        .expect("valid");
        assert!(
            t.transfer_time > h.transfer_time,
            "bi-elliptic is the slow way round"
        );
        // A non-finite altitude surfaces as an error, not a panic.
        assert!(bielliptic_transfer_altitudes(f64::NAN, 35_786.0, 250_000.0).is_err());
    }

    #[test]
    fn elliptical_orbit_matches_a_gto_and_reduces_to_circular() {
        // A geostationary transfer orbit: 300 km perigee, 35 786 km apogee.
        let gto = elliptical_orbit(300.0, 35_786.0);
        assert!((0.70..=0.75).contains(&gto.eccentricity), "GTO e {}", gto.eccentricity);
        let hours = gto.period_s / 3600.0;
        assert!((10.0..=11.1).contains(&hours), "GTO period {hours} h");
        // Vis-viva: ~10.2 km/s at perigee, ~1.6 km/s at apogee; perigee is
        // always the faster end of an eccentric orbit.
        assert!(
            (9_800.0..=10_500.0).contains(&gto.perigee_speed_ms),
            "v_peri {}",
            gto.perigee_speed_ms
        );
        assert!(
            (1_450.0..=1_750.0).contains(&gto.apogee_speed_ms),
            "v_apo {}",
            gto.apogee_speed_ms
        );
        assert!(gto.perigee_speed_ms > gto.apogee_speed_ms);
        // The two altitudes may be given in either order.
        assert_eq!(gto, elliptical_orbit(35_786.0, 300.0));
        // Equal altitudes ⇒ a circle: e = 0, equal speeds, and both period and
        // speed agree with circular_orbit_basics at that altitude.
        let circ = elliptical_orbit(500.0, 500.0);
        assert!(circ.eccentricity < 1e-12, "circular e {}", circ.eccentricity);
        assert!((circ.perigee_speed_ms - circ.apogee_speed_ms).abs() < 1e-9);
        let (v_circ, _, period) = circular_orbit_basics(500.0);
        assert!((circ.perigee_speed_ms - v_circ).abs() < 1e-6, "circular speed");
        assert!((circ.period_s - period).abs() < 1e-3, "circular period");
        // Specific orbital energy ε = −μ/(2a) equals the vis-viva energy
        // ½v² − μ/r evaluated anywhere on the orbit — cross-check at perigee.
        let mu = valenx_astro::constants::MU_EARTH;
        let r_peri = altitude_km_to_radius_m(300.0);
        let eps_visviva = 0.5 * gto.perigee_speed_ms.powi(2) - mu / r_peri;
        assert!(
            (gto.specific_energy_j_per_kg - eps_visviva).abs() / eps_visviva.abs() < 1e-9,
            "ε {} vs vis-viva {}",
            gto.specific_energy_j_per_kg,
            eps_visviva
        );
        // A bound ellipse has negative specific energy.
        assert!(gto.specific_energy_j_per_kg < 0.0, "bound orbit ε < 0");
        // Specific angular momentum h = √(μ·a·(1−e²)) is conserved, so it equals
        // the kinematic product r·v at both apsides.
        let r_apo = altitude_km_to_radius_m(35_786.0);
        let h_peri = r_peri * gto.perigee_speed_ms;
        let h_apo = r_apo * gto.apogee_speed_ms;
        assert!(
            (gto.specific_angular_momentum_m2_s - h_peri).abs() / h_peri < 1e-9,
            "h {} vs r_peri·v_peri {}",
            gto.specific_angular_momentum_m2_s,
            h_peri
        );
        assert!((h_peri - h_apo).abs() / h_peri < 1e-9, "h conserved peri↔apo");
    }

    #[test]
    fn synodic_period_combines_two_orbital_periods() {
        let alt_a = 400.0;
        let alt_b = 35_786.0;
        let (_, _, t_a) = circular_orbit_basics(alt_a);
        let (_, _, t_b) = circular_orbit_basics(alt_b);
        let syn = synodic_period(alt_a, alt_b);
        // T_syn = 1 / |1/T_a − 1/T_b|.
        let expected = 1.0 / (1.0 / t_a - 1.0 / t_b).abs();
        assert!((syn - expected).abs() / expected < 1e-9, "syn {syn} vs {expected}");
        // Always longer than the shorter of the two periods: the faster craft
        // needs more than one of its own orbits to lap the slower one.
        assert!(syn > t_a.min(t_b));
        // Identical orbits never realign → an infinite synodic period.
        assert!(synodic_period(400.0, 400.0).is_infinite());
        // Order-independent.
        assert!((synodic_period(alt_a, alt_b) - synodic_period(alt_b, alt_a)).abs() < 1e-6);
    }

    #[test]
    fn hohmann_phase_angle_matches_known_transfers() {
        // A 1 : 1.524 radius ratio (the Earth→Mars orbital ratio) gives the
        // famous ~44° Hohmann phase angle, independent of absolute scale —
        // 300 km → 3795 km altitude is r2/r1 ≈ 1.524 about Earth.
        let mars_like = hohmann_phase_angle(300.0, 3795.0).to_degrees();
        assert!((mars_like - 44.0).abs() < 2.0, "Mars-ratio phase {mars_like}°");
        // LEO 300 km → GEO 35 786 km: the target leads by ~101° at departure.
        let leo_geo = hohmann_phase_angle(300.0, 35_786.0).to_degrees();
        assert!((leo_geo - 100.7).abs() < 2.0, "LEO→GEO phase {leo_geo}°");
        assert!(hohmann_phase_angle(300.0, 35_786.0) > 0.0, "outward transfer leads");
        // Equal orbits need no phasing.
        assert!(hohmann_phase_angle(500.0, 500.0).abs() < 1e-9);
        // The result is always a valid wrapped angle in (−π, π].
        let inward = hohmann_phase_angle(20_000.0, 300.0);
        assert!(inward > -std::f64::consts::PI - 1e-9 && inward <= std::f64::consts::PI + 1e-9);
    }

    #[test]
    fn escape_delta_v_is_root_two_minus_one_times_circular_speed() {
        let (v_circ, _, _) = circular_orbit_basics(300.0);
        let dv = escape_delta_v_from_circular(300.0);
        // Δv = (√2 − 1)·v_circ exactly.
        assert!((dv - (2.0_f64.sqrt() - 1.0) * v_circ).abs() < 1e-6, "Δv {dv}");
        // From a ~300 km LEO it is ~3.2 km/s.
        assert!((dv / 1000.0 - 3.2).abs() < 0.2, "LEO escape Δv {} km/s", dv / 1000.0);
        // Positive, and less than the circular speed itself.
        assert!(dv > 0.0 && dv < v_circ);
    }

    #[test]
    fn orbit_altitude_for_period_matches_geostationary_and_gps() {
        // One sidereal day → geostationary, ≈ 35 786 km.
        let geo = orbit_altitude_for_period_km(86_164.0);
        assert!((geo - 35_786.0).abs() < 100.0, "geostationary {geo} km");
        // Half a sidereal day → the GPS orbit, ≈ 20 200 km.
        let gps = orbit_altitude_for_period_km(43_082.0);
        assert!((gps - 20_200.0).abs() < 200.0, "GPS {gps} km");
        // Inverse of circular_orbit_basics: the period at the computed altitude
        // round-trips back to the input period.
        let (_, _, t) = circular_orbit_basics(geo);
        assert!((t - 86_164.0).abs() / 86_164.0 < 1e-6, "round-trip period");
    }

    #[test]
    fn injection_delta_v_reduces_to_escape_and_matches_tmi() {
        // v∞ = 0 is the bare escape Δv (no energy left at infinity).
        let esc = escape_delta_v_from_circular(300.0);
        assert!(
            (injection_delta_v(300.0, 0.0) - esc).abs() < 1e-6,
            "v∞=0 → escape Δv"
        );
        // A trans-Mars v∞ ≈ 3 km/s from a 300 km LEO → the classic ~3.6 km/s TMI.
        let tmi = injection_delta_v(300.0, 3000.0);
        assert!((tmi / 1000.0 - 3.6).abs() < 0.2, "TMI {} km/s", tmi / 1000.0);
        // A higher required hyperbolic excess always costs more Δv.
        assert!(injection_delta_v(300.0, 4000.0) > tmi);
    }

    #[test]
    fn flight_path_angle_is_zero_at_apsides_and_peaks_between() {
        use std::f64::consts::PI;
        // A circular orbit (e = 0) is purely tangential everywhere → γ = 0.
        assert!(flight_path_angle(0.0, 1.234).abs() < 1e-12);
        // Both apsides (θ = 0, π) → γ = 0.
        assert!(flight_path_angle(0.6, 0.0).abs() < 1e-12);
        assert!(flight_path_angle(0.6, PI).abs() < 1e-12);
        // e = 0.5 at θ = 90° → atan(0.5) ≈ 26.565°.
        let g = flight_path_angle(0.5, PI / 2.0);
        assert!((g - 0.5_f64.atan()).abs() < 1e-12, "γ {g}");
        assert!((g.to_degrees() - 26.565).abs() < 1e-2);
        // Positive while climbing (0 < θ < π), negative while descending.
        assert!(flight_path_angle(0.5, PI / 4.0) > 0.0);
        assert!(flight_path_angle(0.5, -PI / 4.0) < 0.0);
    }

    #[test]
    fn orbital_speed_at_altitude_matches_apsis_speeds_and_vis_viva() {
        // GTO: 300 km perigee, 35 786 km apogee.
        let gto = elliptical_orbit(300.0, 35_786.0);
        // At the perigee altitude the speed equals the perigee speed.
        let v_peri = orbital_speed_at_altitude(300.0, 35_786.0, 300.0).unwrap();
        assert!((v_peri - gto.perigee_speed_ms).abs() < 1e-6, "v_peri {v_peri}");
        // At the apogee altitude → apogee speed.
        let v_apo = orbital_speed_at_altitude(300.0, 35_786.0, 35_786.0).unwrap();
        assert!((v_apo - gto.apogee_speed_ms).abs() < 1e-6, "v_apo {v_apo}");
        // Speed decreases monotonically from perigee to apogee.
        let v_mid = orbital_speed_at_altitude(300.0, 35_786.0, 18_000.0).unwrap();
        assert!(v_apo < v_mid && v_mid < v_peri, "mid {v_mid}");
        // The perigee/apogee inputs are order-independent.
        assert!(
            (orbital_speed_at_altitude(35_786.0, 300.0, 300.0).unwrap() - v_peri).abs() < 1e-6
        );
        // An altitude the orbit never reaches → None.
        assert!(orbital_speed_at_altitude(300.0, 35_786.0, 50_000.0).is_none());
    }

    #[test]
    fn time_since_perigee_follows_keplers_equation() {
        let gto = elliptical_orbit(300.0, 35_786.0);
        // Perigee (θ=0) → t = 0.
        assert!(time_since_perigee(300.0, 35_786.0, 0.0).abs() < 1e-6);
        // Apogee (θ=180°) → exactly half the orbital period (any eccentricity).
        let t_apo = time_since_perigee(300.0, 35_786.0, 180.0);
        assert!((t_apo - gto.period_s / 2.0).abs() / gto.period_s < 1e-9, "apogee at T/2");
        // A circular orbit sweeps uniformly: θ=90° → a quarter period.
        let circ = elliptical_orbit(500.0, 500.0);
        let t_q = time_since_perigee(500.0, 500.0, 90.0);
        assert!((t_q - circ.period_s / 4.0).abs() / circ.period_s < 1e-9, "circular θ=90 → T/4");
        // On the eccentric GTO the satellite lingers near apogee: the 90°→180°
        // leg takes longer than 0°→90° (Kepler's 2nd law / areal velocity).
        let t90 = time_since_perigee(300.0, 35_786.0, 90.0);
        assert!(t90 < t_apo - t90, "more time spent approaching apogee: {t90} vs {}", t_apo - t90);
    }

    #[test]
    fn horizon_and_coverage_match_a_leo_satellite() {
        // A 400 km LEO satellite: horizon ≈ 2293 km, coverage half-angle ≈ 19.85°.
        let slant = horizon_slant_range_km(400.0);
        assert!((slant - 2293.0).abs() < 5.0, "horizon slant {slant} km");
        let half = coverage_half_angle_deg(400.0);
        assert!((half - 19.85).abs() < 0.1, "coverage half-angle {half}°");
        // At the surface both vanish.
        assert!(horizon_slant_range_km(0.0).abs() < 1e-6);
        assert!(coverage_half_angle_deg(0.0).abs() < 1e-9);
        // A higher orbit sees a farther horizon and covers a wider footprint.
        assert!(horizon_slant_range_km(35_786.0) > slant);
        assert!(coverage_half_angle_deg(35_786.0) > half);
        // Geometric identity slant = √(h²+2Rh).
        let r = valenx_astro::constants::R_EARTH;
        let h = 1_000.0e3;
        let expected = (h * h + 2.0 * r * h).sqrt() / 1000.0;
        assert!((horizon_slant_range_km(1000.0) - expected).abs() < 1e-6);
    }

    #[test]
    fn formatters_never_panic_on_edge_values() {
        // Pure formatters must tolerate non-finite / extreme inputs.
        let _ = format_delta_v(f64::INFINITY);
        let _ = format_delta_v(0.0);
        let _ = format_duration(f64::NAN);
        let _ = format_duration(f64::INFINITY);
        let _ = format_duration(45.0);
        assert_eq!(format_duration(f64::NAN), "—");
    }

    #[test]
    fn empty_stage_stack_builds_a_vehicle_that_fails_validation_not_a_panic() {
        // A user who deletes every stage must get a clean validation
        // error from the backend, never a panic in the builder.
        let mut form = AscentForm::default();
        form.stages.clear();
        let vehicle = form.build_vehicle();
        assert!(vehicle.validate().is_err());
    }
}

//! # valenx-marine-hydro
//!
//! Calm-water ship/boat hull **resistance** and **powering** estimation:
//! how much drag a displacement hull feels at a given speed, and the
//! effective (towed) power needed to push it through the water.
//!
//! This is the in-house companion to [`valenx-marine`] (which does the
//! box-form *hydrostatics* — displacement, `KB`, `BM`, `GM`). Where that
//! crate answers "does it float and is it stable?", this crate answers
//! "how much does it cost, in drag and power, to drive it forward?".
//!
//! ## What it computes
//!
//! For a hull at speed `V` in calm water it splits the total resistance into
//! physically-meaningful components and reports the effective power:
//!
//! - **Frictional resistance** `R_f = 1/2 rho V^2 S C_f`, where the skin-friction
//!   coefficient `C_f` comes from the **ITTC-1957 model-ship correlation line**
//!   `C_f = 0.075 / (log10(Re) - 2)^2`, with `Re = V L / nu`
//!   ([`ittc57_friction_coefficient`], [`frictional_resistance`]).
//! - **Form (viscous-pressure) factor** `(1 + k)` that scales the flat-plate
//!   friction up to the real three-dimensional hull, via the
//!   **Holtrop-Mennen** regression ([`Hull::form_factor`]).
//! - **Wave-making (residuary) resistance** `R_w`, the energy radiated into the
//!   wave system the hull makes, via the **Holtrop-Mennen** statistical method
//!   ([`Hull::wave_resistance`]).
//! - **Model-ship correlation allowance** `R_a = 1/2 rho V^2 S C_a` (a small
//!   roughness / scale-effect term).
//! - **Total resistance** `R_t = R_f (1 + k) + R_w + R_a`
//!   ([`Hull::total_resistance`]).
//! - **Effective power** `P_e = R_t V` ([`Hull::effective_power`]).
//! - The dimensionless **Froude number** `Fn = V / sqrt(g L)` and **Reynolds
//!   number** `Re = V L / nu` at each speed.
//!
//! [`Hull::resistance_at`] returns all of the above at one speed as a
//! [`ResistancePoint`]; [`resistance_curve`] sweeps a speed range and returns a
//! [`ResistanceCurve`] suitable for a future plot or design-speed readout.
//!
//! ## Wetted surface
//!
//! Resistance needs the **wetted surface** `S`. You can pass a measured `S`
//! directly, or let the crate estimate it from the principal dimensions with
//! the Holtrop-Mennen wetted-surface regression ([`Hull::wetted_surface`]).
//!
//! ## Honest scope and accuracy
//!
//! Research / **preliminary-design** grade. This is a *statistical / empirical*
//! estimator, **not** a free-surface Navier-Stokes (CFD) solver: it computes no
//! flow field, resolves no boundary layer, and models no actual free surface.
//! It is the kind of first-cut powering estimate used at the concept stage,
//! before a model test or a full RANS/free-surface computation.
//!
//! What is **rigorously validated** (**digit-exact**, to 4+ significant
//! figures, against the published Holtrop & Mennen (1982) worked example — see
//! the `validation` tests):
//!
//! - the ITTC-1957 friction line `C_f` (also a pinnable closed form);
//! - the Holtrop wetted-surface estimate `S` (paper: `7381.45 m^2`);
//! - the frictional resistance `R_f` (paper: `869.63 kN`);
//! - the geometry coefficients (`C_b`, `C_p`) and the Reynolds/Froude numbers.
//!
//! What is validated **to a stated few-per-cent tolerance** (the Holtrop
//! regression is statistical, so a digit-exact claim would be dishonest, but
//! the implementation *does* reproduce the paper's published intermediate and
//! component values):
//!
//! - the half-angle of entrance `i_E` (computed 11.84 deg vs paper 12.08 deg,
//!   ~2 %);
//! - the **form factor `(1 + k)`** (computed 1.162 vs the value implied by the
//!   paper's own resistance balance, 1.156, ~0.6 %);
//! - the **wave-making resistance `R_w`** (computed ~555 kN vs paper 557.11 kN,
//!   ~0.4 %);
//! - the reconstructed bare-hull total and effective power (within ~1-2 % of
//!   the paper once its separate appendage term is accounted for).
//!
//! Even so, treat the wave term as **preliminary**: the Holtrop method carries
//! the method's own scatter (a few per cent on a hull inside its calibration
//! envelope, more outside it), so `R_w` is an estimate, not a guarantee, and we
//! deliberately do **not** claim digit-exactness for it. See
//! [`Hull::wave_resistance`] for the per-method caveats.
//!
//! Excluded entirely (a real powering job needs them, this does not provide
//! them): appendage resistance beyond a single lumped factor, air/wind
//! resistance, added resistance in waves, shallow-water effects, the propeller
//! and the propulsive coefficients (wake, thrust deduction, open-water and
//! hull efficiencies) that turn effective power into delivered/shaft power, and
//! any trial-correction or class margin. **Not** a substitute for a licensed
//! naval architect, model-test data, or class-approved software, and not for
//! design, construction, or any decision affecting safety.
//!
//! Units are SI and must be self-consistent: metres for lengths, m/s for
//! speed, kg/m^3 for density, m^2/s for kinematic viscosity. Resistances are
//! newtons (N), powers watts (W). Helpers are provided for knots and kW.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::{Deserialize, Serialize};

// Re-export the hydrostatics hull so callers can build geometry once and feed
// it to both crates without a second `use`.
pub use valenx_marine::{
    Hull as HydrostaticHull, FRESHWATER_DENSITY, GRAVITY, SEAWATER_DENSITY,
};

/// One knot in metres per second (1 nmi/h = 1852 m / 3600 s).
pub const KNOT_MS: f64 = 1852.0 / 3600.0;

/// Standard ITTC fresh-water kinematic viscosity at 15 degrees C (m^2/s).
pub const FRESHWATER_NU_15C: f64 = 1.139_02e-6;

/// Standard ITTC sea-water kinematic viscosity at 15 degrees C (m^2/s).
///
/// This is the value implied by the Holtrop & Mennen (1982) worked example
/// (it reproduces that example's `Re` and `C_f` to four significant figures).
pub const SEAWATER_NU_15C: f64 = 1.188_31e-6;

/// A default model-ship correlation allowance `C_a` for a medium hull.
///
/// Holtrop's own `C_a` regression depends on length and form; this constant is
/// a representative mid-range value (about `4e-4`) for when you have nothing
/// better. Pass your own via [`WaterProperties::with_correlation_allowance`].
pub const DEFAULT_CORRELATION_ALLOWANCE: f64 = 0.0004;

/// Convert a speed in knots to metres per second.
#[must_use]
pub fn knots_to_ms(knots: f64) -> f64 {
    knots * KNOT_MS
}

/// Convert a speed in metres per second to knots.
#[must_use]
pub fn ms_to_knots(ms: f64) -> f64 {
    ms / KNOT_MS
}

/// An out-of-domain resistance/powering input.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum HydroError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// A quantity that must be finite (but may be zero/negative) was not.
    #[error("{quantity} must be finite, got {value}")]
    NonFinite {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// A coefficient was outside its physical range.
    #[error("{quantity} must be in {lo}..={hi}, got {value}")]
    OutOfRange {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The lower bound (inclusive).
        lo: f64,
        /// The upper bound (inclusive).
        hi: f64,
        /// The offending value.
        value: f64,
    },
    /// A speed range was empty or had a non-positive step.
    #[error("invalid speed range: start={start}, end={end}, step={step}")]
    BadRange {
        /// Range start (m/s).
        start: f64,
        /// Range end (m/s).
        end: f64,
        /// Range step (m/s).
        step: f64,
    },
}

fn require_positive(quantity: &'static str, value: f64) -> Result<f64, HydroError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(HydroError::NonPositive { quantity, value })
    }
}

fn require_finite(quantity: &'static str, value: f64) -> Result<f64, HydroError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HydroError::NonFinite { quantity, value })
    }
}

fn require_range(
    quantity: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
) -> Result<f64, HydroError> {
    if value.is_finite() && value >= lo && value <= hi {
        Ok(value)
    } else {
        Err(HydroError::OutOfRange {
            quantity,
            lo,
            hi,
            value,
        })
    }
}

// ===========================================================================
// Water + the ITTC-1957 frictional line (the rigorously-validated core)
// ===========================================================================

/// The fluid the hull moves through: density, kinematic viscosity and the
/// model-ship correlation allowance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WaterProperties {
    /// Density `rho` (kg/m^3).
    pub density: f64,
    /// Kinematic viscosity `nu` (m^2/s).
    pub kinematic_viscosity: f64,
    /// Model-ship correlation allowance `C_a` (dimensionless), the small extra
    /// flat-plate-like term `R_a = 1/2 rho V^2 S C_a` for roughness/scale.
    pub correlation_allowance: f64,
}

impl WaterProperties {
    /// Standard ITTC sea water at 15 degrees C
    /// (`rho = 1025`, `nu = 1.18831e-6`) with the default correlation
    /// allowance. These are the values behind the Holtrop & Mennen (1982)
    /// worked example used in this crate's validation tests.
    #[must_use]
    pub fn seawater() -> Self {
        Self {
            density: SEAWATER_DENSITY,
            kinematic_viscosity: SEAWATER_NU_15C,
            correlation_allowance: DEFAULT_CORRELATION_ALLOWANCE,
        }
    }

    /// Standard ITTC fresh water at 15 degrees C
    /// (`rho = 1000`, `nu = 1.13902e-6`) with the default correlation
    /// allowance.
    #[must_use]
    pub fn freshwater() -> Self {
        Self {
            density: FRESHWATER_DENSITY,
            kinematic_viscosity: FRESHWATER_NU_15C,
            correlation_allowance: DEFAULT_CORRELATION_ALLOWANCE,
        }
    }

    /// Return a copy with a different correlation allowance `C_a`.
    #[must_use]
    pub fn with_correlation_allowance(mut self, c_a: f64) -> Self {
        self.correlation_allowance = c_a;
        self
    }

    /// Validate the water properties (density and viscosity positive;
    /// correlation allowance finite).
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] when any field is out of its physical domain.
    pub fn validate(&self) -> Result<(), HydroError> {
        require_positive("water density", self.density)?;
        require_positive("kinematic viscosity", self.kinematic_viscosity)?;
        require_finite("correlation allowance", self.correlation_allowance)?;
        Ok(())
    }
}

impl Default for WaterProperties {
    fn default() -> Self {
        Self::seawater()
    }
}

/// Reynolds number `Re = V L / nu` for a hull of length `L` at speed `V`.
///
/// # Errors
///
/// Returns [`HydroError`] if any argument is not finite and positive.
pub fn reynolds_number(
    speed_ms: f64,
    length_m: f64,
    kinematic_viscosity: f64,
) -> Result<f64, HydroError> {
    let v = require_positive("speed", speed_ms)?;
    let l = require_positive("length", length_m)?;
    let nu = require_positive("kinematic viscosity", kinematic_viscosity)?;
    Ok(v * l / nu)
}

/// Froude number `Fn = V / sqrt(g L)` for a hull of length `L` at speed `V`.
///
/// # Errors
///
/// Returns [`HydroError`] if any argument is not finite and positive.
pub fn froude_number(speed_ms: f64, length_m: f64) -> Result<f64, HydroError> {
    let v = require_positive("speed", speed_ms)?;
    let l = require_positive("length", length_m)?;
    Ok(v / (GRAVITY * l).sqrt())
}

/// The **ITTC-1957 model-ship correlation line** skin-friction coefficient
/// `C_f = 0.075 / (log10(Re) - 2)^2`.
///
/// This is the international-standard flat-plate friction line adopted by the
/// 8th ITTC (1957) and used essentially universally for extrapolating model
/// resistance to full scale. It is an exact closed form, so it is pinnable:
/// e.g. `Re = 1e9` gives `C_f = 0.0015306`, and `Re = 1e7` gives exactly
/// `C_f = 0.003`.
///
/// # Errors
///
/// Returns [`HydroError`] if `Re` is not finite and positive, or if
/// `log10(Re) <= 2` (i.e. `Re <= 100`, where the line is undefined / singular).
pub fn ittc57_friction_coefficient(reynolds: f64) -> Result<f64, HydroError> {
    let re = require_positive("Reynolds number", reynolds)?;
    let denom = re.log10() - 2.0;
    if denom <= 0.0 {
        return Err(HydroError::OutOfRange {
            quantity: "Reynolds number (log10(Re) must exceed 2)",
            lo: 100.0,
            hi: f64::INFINITY,
            value: re,
        });
    }
    Ok(0.075 / (denom * denom))
}

/// Frictional resistance `R_f = 1/2 rho V^2 S C_f` (N) for a flat-plate-like
/// skin-friction coefficient `C_f` and wetted surface `S`.
///
/// # Errors
///
/// Returns [`HydroError`] if density, speed or wetted surface are not finite
/// and positive, or if `C_f` is not finite and non-negative.
pub fn frictional_resistance(
    density: f64,
    speed_ms: f64,
    wetted_surface_m2: f64,
    friction_coefficient: f64,
) -> Result<f64, HydroError> {
    let rho = require_positive("water density", density)?;
    let v = require_positive("speed", speed_ms)?;
    let s = require_positive("wetted surface", wetted_surface_m2)?;
    let cf = require_finite("friction coefficient", friction_coefficient)?;
    if cf < 0.0 {
        return Err(HydroError::NonPositive {
            quantity: "friction coefficient",
            value: cf,
        });
    }
    Ok(0.5 * rho * v * v * s * cf)
}

// ===========================================================================
// The resistance hull (Holtrop-Mennen geometry)
// ===========================================================================

/// A hull described well enough for a Holtrop-Mennen resistance estimate.
///
/// The principal dimensions are the same ones [`valenx_marine::Hull`] uses;
/// the extra form coefficients (`C_m`, `C_wp`, `lcb`, bulb/transom areas, the
/// stern-shape parameter) are what the Holtrop wave-making and form-factor
/// regressions need. Build one directly, or with [`Hull::from_hydrostatic`]
/// from a [`valenx_marine::Hull`] plus the extra coefficients.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Hull {
    /// Waterline length `L_wl` (m).
    pub length_m: f64,
    /// Waterline beam / breadth `B` (m).
    pub beam_m: f64,
    /// Draft `T` (m). For a hull with different fore/aft drafts use the mean.
    pub draft_m: f64,
    /// Displaced volume `nabla` (m^3).
    pub volume_m3: f64,
    /// Midship-section coefficient `C_m` in `(0, 1]` (immersed midship area /
    /// `B*T`).
    pub midship_coefficient: f64,
    /// Waterplane-area coefficient `C_wp` in `(0, 1]` (waterplane area /
    /// `L*B`).
    pub waterplane_coefficient: f64,
    /// Longitudinal centre of buoyancy as a percentage of `L` forward of
    /// amidships (`+` forward, `-` aft). Holtrop's `lcb`.
    pub lcb_percent: f64,
    /// Transverse area of the bulbous bow at the forward perpendicular
    /// `A_bt` (m^2); `0.0` for no bulb.
    pub bulb_area_m2: f64,
    /// Height of the centre of the bulb area above the keel `h_b` (m); ignored
    /// when there is no bulb.
    pub bulb_centre_m: f64,
    /// Immersed transom area at rest `A_t` (m^2); `0.0` for no transom.
    pub transom_area_m2: f64,
    /// Stern-shape parameter `C_stern` (Holtrop): about `-25` for a pram /
    /// barge stern, `-10` for a V-shaped stern, `0` for a normal stern, `+10`
    /// for a U-shaped stern with a Hogner form.
    pub stern_parameter: f64,
    /// Optional measured wetted surface `S` (m^2). When `None`, the Holtrop
    /// wetted-surface estimate ([`Hull::wetted_surface`]) is used.
    pub wetted_surface_m2: Option<f64>,
}

impl Hull {
    /// Build a validated resistance hull from principal dimensions and form
    /// coefficients.
    ///
    /// `length_m`, `beam_m`, `draft_m`, `volume_m3` must be finite and
    /// positive; `midship_coefficient` and `waterplane_coefficient` must lie in
    /// `(0, 1]`; `lcb_percent` finite; bulb/transom areas finite and
    /// non-negative; `stern_parameter` finite. A supplied `wetted_surface_m2`
    /// (if any) must be finite and positive.
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] when any input is out of its physical domain.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        length_m: f64,
        beam_m: f64,
        draft_m: f64,
        volume_m3: f64,
        midship_coefficient: f64,
        waterplane_coefficient: f64,
        lcb_percent: f64,
        bulb_area_m2: f64,
        bulb_centre_m: f64,
        transom_area_m2: f64,
        stern_parameter: f64,
        wetted_surface_m2: Option<f64>,
    ) -> Result<Self, HydroError> {
        let length_m = require_positive("length", length_m)?;
        let beam_m = require_positive("beam", beam_m)?;
        let draft_m = require_positive("draft", draft_m)?;
        let volume_m3 = require_positive("displaced volume", volume_m3)?;
        let midship_coefficient =
            require_range("midship coefficient", midship_coefficient, f64::MIN_POSITIVE, 1.0)?;
        let waterplane_coefficient = require_range(
            "waterplane coefficient",
            waterplane_coefficient,
            f64::MIN_POSITIVE,
            1.0,
        )?;
        let lcb_percent = require_finite("lcb percent", lcb_percent)?;
        if !(bulb_area_m2.is_finite() && bulb_area_m2 >= 0.0) {
            return Err(HydroError::NonFinite {
                quantity: "bulb area",
                value: bulb_area_m2,
            });
        }
        let bulb_centre_m = require_finite("bulb centre", bulb_centre_m)?;
        if !(transom_area_m2.is_finite() && transom_area_m2 >= 0.0) {
            return Err(HydroError::NonFinite {
                quantity: "transom area",
                value: transom_area_m2,
            });
        }
        let stern_parameter = require_finite("stern parameter", stern_parameter)?;
        if let Some(s) = wetted_surface_m2 {
            require_positive("wetted surface", s)?;
        }
        Ok(Self {
            length_m,
            beam_m,
            draft_m,
            volume_m3,
            midship_coefficient,
            waterplane_coefficient,
            lcb_percent,
            bulb_area_m2,
            bulb_centre_m,
            transom_area_m2,
            stern_parameter,
            wetted_surface_m2,
        })
    }

    /// Build a resistance hull from a [`valenx_marine::Hull`] (which already
    /// carries `L`, `B`, `T` and the block coefficient `C_b`, from which the
    /// displaced volume is taken) plus the extra Holtrop coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] when any input is out of its physical domain.
    #[allow(clippy::too_many_arguments)]
    pub fn from_hydrostatic(
        hull: &valenx_marine::Hull,
        midship_coefficient: f64,
        waterplane_coefficient: f64,
        lcb_percent: f64,
        bulb_area_m2: f64,
        bulb_centre_m: f64,
        transom_area_m2: f64,
        stern_parameter: f64,
        wetted_surface_m2: Option<f64>,
    ) -> Result<Self, HydroError> {
        Self::new(
            hull.length_m,
            hull.beam_m,
            hull.draft_m,
            hull.displaced_volume(),
            midship_coefficient,
            waterplane_coefficient,
            lcb_percent,
            bulb_area_m2,
            bulb_centre_m,
            transom_area_m2,
            stern_parameter,
            wetted_surface_m2,
        )
    }

    /// Block coefficient `C_b = nabla / (L B T)`.
    #[must_use]
    pub fn block_coefficient(&self) -> f64 {
        self.volume_m3 / (self.length_m * self.beam_m * self.draft_m)
    }

    /// Prismatic coefficient `C_p = C_b / C_m = nabla / (A_m L)`.
    #[must_use]
    pub fn prismatic_coefficient(&self) -> f64 {
        self.block_coefficient() / self.midship_coefficient
    }

    /// Holtrop **length of run** `L_R = L (1 - C_p + 0.06 C_p lcb / (4 C_p - 1))`
    /// (m), the effective length over which the afterbody runs out.
    #[must_use]
    pub fn length_of_run(&self) -> f64 {
        let cp = self.prismatic_coefficient();
        self.length_m * (1.0 - cp + 0.06 * cp * self.lcb_percent / (4.0 * cp - 1.0))
    }

    /// **Holtrop-Mennen wetted-surface estimate** `S` (m^2) from the principal
    /// dimensions and form coefficients:
    ///
    /// ```text
    /// S = L (2T + B) sqrt(C_m)
    ///       (0.453 + 0.4425 C_b - 0.2862 C_m - 0.003467 B/T + 0.3696 C_wp)
    ///     + 2.38 A_bt / C_b
    /// ```
    ///
    /// If a measured [`Hull::wetted_surface_m2`] was supplied it is returned
    /// instead. This estimate is validated digit-exact against the Holtrop &
    /// Mennen (1982) worked example (`S = 7381.45 m^2`).
    #[must_use]
    pub fn wetted_surface(&self) -> f64 {
        if let Some(s) = self.wetted_surface_m2 {
            return s;
        }
        let l = self.length_m;
        let b = self.beam_m;
        let t = self.draft_m;
        let cb = self.block_coefficient();
        let cm = self.midship_coefficient;
        let cwp = self.waterplane_coefficient;
        l * (2.0 * t + b)
            * cm.sqrt()
            * (0.453 + 0.4425 * cb - 0.2862 * cm - 0.003_467 * b / t + 0.3696 * cwp)
            + 2.38 * self.bulb_area_m2 / cb
    }

    /// **Holtrop-Mennen form factor** `(1 + k)` — the viscous-pressure
    /// multiplier on the flat-plate friction (1984 re-analysis form):
    ///
    /// ```text
    /// 1 + k = c13 { 0.93 + c12 (B/L_R)^0.92497 (0.95 - C_p)^-0.521448
    ///                       (1 - C_p + 0.0225 lcb)^0.6906 }
    /// ```
    ///
    /// with `c13 = 1 + 0.003 C_stern` and `c12` a function of `T/L`.
    ///
    /// This is the genuine Holtrop regression and gives physically-sensible
    /// factors (typically ~1.1-1.3 for normal forms). For the Holtrop & Mennen
    /// (1982) worked example it returns `1.162`, within ~0.6 % of the `1.156`
    /// implied by that paper's own resistance balance. It is validated to that
    /// stated tolerance (not digit-exact, as befits a regression).
    #[must_use]
    pub fn form_factor(&self) -> f64 {
        let l = self.length_m;
        let b = self.beam_m;
        let t = self.draft_m;
        let cp = self.prismatic_coefficient();
        let lr = self.length_of_run();
        let lcb = self.lcb_percent;

        // c12: depends on T/L.
        let t_over_l = t / l;
        let c12 = if t_over_l > 0.05 {
            t_over_l.powf(0.228_844_6)
        } else if t_over_l > 0.02 {
            48.20 * (t_over_l - 0.02).powf(2.078) + 0.479_948
        } else {
            0.479_948
        };
        // c13: stern-shape correction.
        let c13 = 1.0 + 0.003 * self.stern_parameter;

        c13 * (0.93
            + c12
                * (b / lr).powf(0.924_97)
                * (0.95 - cp).powf(-0.521_448)
                * (1.0 - cp + 0.0225 * lcb).powf(0.6906))
    }

    /// **Holtrop-Mennen wave-making (residuary) resistance** `R_w` (N) at speed
    /// `V`, for the low/moderate-speed branch (`Fn <= ~0.4`, the form most ship
    /// hulls operate in):
    ///
    /// ```text
    /// R_w = c1 c2 c5 nabla rho g exp( m1 Fn^d + m2 cos(lambda Fn^-2) )
    /// ```
    ///
    /// with the standard Holtrop `c1, c2, c3, c5, c7, c15, c16, m1, m2, lambda`
    /// sub-coefficients derived from the hull form, the half-angle of entrance,
    /// and the bulb/transom geometry.
    ///
    /// **Preliminary** (see the crate-level honest-scope note). For the
    /// Holtrop & Mennen (1982) worked example this returns ~555 kN, within
    /// ~0.4 % of the paper's published `R_w = 557.11 kN`, so the implementation
    /// is faithful to the regression. But the wave term is the dominant source
    /// of uncertainty in any statistical powering estimate: the Holtrop method
    /// carries several per cent scatter inside its calibration envelope and
    /// more outside it. We deliberately do **not** claim a *digit-exact* wave
    /// match — treat `R_w` as an estimate, not a guarantee.
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] if `speed_ms` is not finite and positive or if
    /// the water properties are invalid.
    pub fn wave_resistance(
        &self,
        speed_ms: f64,
        water: &WaterProperties,
    ) -> Result<f64, HydroError> {
        let v = require_positive("speed", speed_ms)?;
        water.validate()?;
        let l = self.length_m;
        let b = self.beam_m;
        let t = self.draft_m;
        let cp = self.prismatic_coefficient();
        let cb = self.block_coefficient();
        let cm = self.midship_coefficient;
        let cwp = self.waterplane_coefficient;
        let nabla = self.volume_m3;
        let fn_ = froude_number(v, l)?;

        // c7 depends on B/L.
        let b_over_l = b / l;
        let c7 = if b_over_l < 0.11 {
            0.229_577 * b_over_l.powf(0.333_333)
        } else if b_over_l <= 0.25 {
            b_over_l
        } else {
            0.5 - 0.0625 * l / b
        };

        // Half-angle of entrance i_E (degrees) via Holtrop's regression. A
        // tiny floor keeps c1 finite for the degenerate i_E -> 90 deg edge.
        let lr = self.length_of_run();
        let i_e =
            canonical_half_entrance_angle(l, b, t, cp, cwp, cb, lr, self.lcb_percent).min(89.0);

        // c1 (forebody / entrance) coefficient.
        let c1 = 2_223_105.0
            * c7.powf(3.786_13)
            * (t / b).powf(1.078_961)
            * (90.0 - i_e).powf(-1.377_565);

        // c3, c2 (bulbous-bow influence on wave resistance), c5 (transom).
        let c3 = 0.566_15 * self.bulb_area_m2.powf(1.5)
            / (b * t * (0.31 * self.bulb_area_m2.sqrt() + t - self.bulb_centre_m));
        let c2 = (-1.89 * c3.sqrt()).exp();
        let c5 = 1.0 - 0.8 * self.transom_area_m2 / (b * t * cm);

        // lambda, d, m1, m2.
        let l_over_b = l / b;
        let lambda = if l_over_b < 12.0 {
            1.446 * cp - 0.03 * l_over_b
        } else {
            1.446 * cp - 0.36
        };
        let d = -0.9;
        // c16 and m1.
        let c16 = if cp < 0.8 {
            8.078_1 * cp - 13.873_7 * cp * cp + 6.984_388 * cp * cp * cp
        } else {
            1.732_5 - 0.7067 * cp
        };
        let m1 = 0.014_04 * l / t - 1.752_54 * nabla.powf(1.0 / 3.0) / l - 4.793_23 * b / l
            - c16;
        // c15 and m2.
        let l3_over_v = l * l * l / nabla;
        let c15 = if l3_over_v < 512.0 {
            -1.693_85
        } else if l3_over_v <= 1726.91 {
            -1.693_85 + (l / nabla.powf(1.0 / 3.0) - 8.0) / 2.36
        } else {
            0.0
        };
        let m2 = c15 * cp * cp * (-0.1 * fn_.powi(-2)).exp();

        let rho = water.density;
        let g = GRAVITY;
        Ok(c1 * c2 * c5
            * nabla
            * rho
            * g
            * (m1 * fn_.powf(d) + m2 * (lambda * fn_.powi(-2)).cos()).exp())
    }

    /// Resistance breakdown and effective power at a single speed `V`.
    ///
    /// `R_t = R_f (1 + k) + R_w + R_a`, `P_e = R_t V`. The wetted surface is
    /// the measured value if supplied, otherwise the Holtrop estimate.
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] if the speed is not finite and positive or the
    /// water properties are invalid.
    pub fn resistance_at(
        &self,
        speed_ms: f64,
        water: &WaterProperties,
    ) -> Result<ResistancePoint, HydroError> {
        let v = require_positive("speed", speed_ms)?;
        water.validate()?;

        let s = self.wetted_surface();
        let re = reynolds_number(v, self.length_m, water.kinematic_viscosity)?;
        let fn_ = froude_number(v, self.length_m)?;
        let cf = ittc57_friction_coefficient(re)?;
        let one_plus_k = self.form_factor();

        let r_f = frictional_resistance(water.density, v, s, cf)?;
        let r_viscous = r_f * one_plus_k;
        let r_w = self.wave_resistance(v, water)?;
        let r_a = 0.5 * water.density * v * v * s * water.correlation_allowance;
        let r_t = r_viscous + r_w + r_a;
        let p_e = r_t * v;

        Ok(ResistancePoint {
            speed_ms: v,
            speed_knots: ms_to_knots(v),
            froude_number: fn_,
            reynolds_number: re,
            wetted_surface_m2: s,
            friction_coefficient: cf,
            form_factor: one_plus_k,
            frictional_resistance_n: r_f,
            viscous_resistance_n: r_viscous,
            wave_resistance_n: r_w,
            correlation_resistance_n: r_a,
            total_resistance_n: r_t,
            effective_power_w: p_e,
        })
    }

    /// Convenience: [`Hull::total_resistance`] in newtons at speed `V`.
    ///
    /// # Errors
    ///
    /// See [`Hull::resistance_at`].
    pub fn total_resistance(
        &self,
        speed_ms: f64,
        water: &WaterProperties,
    ) -> Result<f64, HydroError> {
        Ok(self.resistance_at(speed_ms, water)?.total_resistance_n)
    }

    /// Convenience: effective power `P_e = R_t V` (W) at speed `V`.
    ///
    /// # Errors
    ///
    /// See [`Hull::resistance_at`].
    pub fn effective_power(
        &self,
        speed_ms: f64,
        water: &WaterProperties,
    ) -> Result<f64, HydroError> {
        Ok(self.resistance_at(speed_ms, water)?.effective_power_w)
    }
}

/// Canonical Holtrop half-angle of entrance `i_E` (degrees).
///
/// `i_E = 1 + 89 exp{ -(L/B)^0.80856 (1-C_wp)^0.30484 (1 - C_p - 0.0225 lcb)^0.6367
/// (L_R/B)^0.34574 (100 nabla / L^3)^0.16302 }` (Holtrop & Mennen 1982/1984).
///
/// The argument list mirrors the published formula's variables one-to-one;
/// bundling them into a struct would only obscure that correspondence.
#[allow(clippy::too_many_arguments)]
fn canonical_half_entrance_angle(
    l: f64,
    b: f64,
    t: f64,
    cp: f64,
    cwp: f64,
    cb: f64,
    lr: f64,
    lcb: f64,
) -> f64 {
    let nabla = cb * l * b * t; // = C_b L B T
    1.0 + 89.0
        * (-(l / b).powf(0.808_56)
            * (1.0 - cwp).powf(0.304_84)
            * (1.0 - cp - 0.0225 * lcb).powf(0.6367)
            * (lr / b).powf(0.345_74)
            * (100.0 * nabla / (l * l * l)).powf(0.163_02))
        .exp()
}

// ===========================================================================
// Result types + the curve API
// ===========================================================================

/// The full resistance breakdown and effective power at one speed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ResistancePoint {
    /// Speed `V` (m/s).
    pub speed_ms: f64,
    /// Speed `V` (knots).
    pub speed_knots: f64,
    /// Froude number `Fn = V / sqrt(g L)`.
    pub froude_number: f64,
    /// Reynolds number `Re = V L / nu`.
    pub reynolds_number: f64,
    /// Wetted surface `S` used (m^2) — measured or Holtrop-estimated.
    pub wetted_surface_m2: f64,
    /// ITTC-57 skin-friction coefficient `C_f`.
    pub friction_coefficient: f64,
    /// Form factor `(1 + k)`.
    pub form_factor: f64,
    /// Bare flat-plate frictional resistance `R_f` (N).
    pub frictional_resistance_n: f64,
    /// Viscous resistance `R_f (1 + k)` (N).
    pub viscous_resistance_n: f64,
    /// Wave-making (residuary) resistance `R_w` (N) — **preliminary**.
    pub wave_resistance_n: f64,
    /// Model-ship correlation-allowance resistance `R_a` (N).
    pub correlation_resistance_n: f64,
    /// Total resistance `R_t = R_f (1+k) + R_w + R_a` (N).
    pub total_resistance_n: f64,
    /// Effective (towed) power `P_e = R_t V` (W).
    pub effective_power_w: f64,
}

impl ResistancePoint {
    /// Total resistance in kilonewtons.
    #[must_use]
    pub fn total_resistance_kn(&self) -> f64 {
        self.total_resistance_n / 1000.0
    }

    /// Effective power in kilowatts.
    #[must_use]
    pub fn effective_power_kw(&self) -> f64 {
        self.effective_power_w / 1000.0
    }
}

/// A resistance/powering curve over a swept speed range, plus the hull and
/// water it was computed for. Ready to feed a plot or a design-speed readout.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResistanceCurve {
    /// The hull the curve was computed for.
    pub hull: Hull,
    /// The water the curve was computed in.
    pub water: WaterProperties,
    /// One [`ResistancePoint`] per speed, in ascending speed order.
    pub points: Vec<ResistancePoint>,
}

impl ResistanceCurve {
    /// The point in the curve closest to a target speed (m/s), if any.
    #[must_use]
    pub fn nearest(&self, speed_ms: f64) -> Option<&ResistancePoint> {
        self.points.iter().min_by(|a, b| {
            (a.speed_ms - speed_ms)
                .abs()
                .total_cmp(&(b.speed_ms - speed_ms).abs())
        })
    }

    /// The point of maximum total resistance in the curve, if any.
    #[must_use]
    pub fn peak_resistance(&self) -> Option<&ResistancePoint> {
        self.points
            .iter()
            .max_by(|a, b| a.total_resistance_n.total_cmp(&b.total_resistance_n))
    }
}

/// An inclusive speed range to sweep, in m/s.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SpeedRange {
    /// First speed (m/s).
    pub start_ms: f64,
    /// Last speed (m/s), inclusive (within a small tolerance of the step).
    pub end_ms: f64,
    /// Step (m/s).
    pub step_ms: f64,
}

impl SpeedRange {
    /// A range from `start` to `end` (inclusive) in `n` equal steps (`n >= 1`).
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] if `start`/`end` are not finite, `end < start`,
    /// or `n == 0`.
    pub fn linspace(start_ms: f64, end_ms: f64, n: usize) -> Result<Self, HydroError> {
        let start = require_positive("range start", start_ms)?;
        let end = require_finite("range end", end_ms)?;
        if end < start || n == 0 {
            return Err(HydroError::BadRange {
                start,
                end,
                step: f64::NAN,
            });
        }
        let step = if n == 1 {
            (end - start).max(1.0)
        } else {
            (end - start) / (n as f64 - 1.0)
        };
        Ok(Self {
            start_ms: start,
            end_ms: end,
            step_ms: step,
        })
    }

    /// A range from `start` to `end` (inclusive) in steps of `step`.
    ///
    /// # Errors
    ///
    /// Returns [`HydroError`] if any value is not finite, `step <= 0`, or
    /// `end < start`.
    pub fn stepped(start_ms: f64, end_ms: f64, step_ms: f64) -> Result<Self, HydroError> {
        let start = require_positive("range start", start_ms)?;
        let end = require_finite("range end", end_ms)?;
        let step = require_positive("range step", step_ms)?;
        if end < start {
            return Err(HydroError::BadRange { start, end, step });
        }
        Ok(Self {
            start_ms: start,
            end_ms: end,
            step_ms: step,
        })
    }

    /// The speeds (m/s) this range yields, ascending and inclusive of `end`
    /// (within half a step).
    #[must_use]
    pub fn speeds(&self) -> Vec<f64> {
        let mut out = Vec::new();
        let mut v = self.start_ms;
        // Include end within a half-step tolerance to avoid float drift.
        while v <= self.end_ms + self.step_ms * 0.5 {
            out.push(v.min(self.end_ms));
            v += self.step_ms;
        }
        // Guarantee the exact end point is present.
        if let Some(last) = out.last() {
            if (*last - self.end_ms).abs() > self.step_ms * 1e-9 {
                out.push(self.end_ms);
            }
        }
        out
    }
}

/// Sweep a [`SpeedRange`] and return the full [`ResistanceCurve`] for a hull in
/// the given water. This is the high-level entry point for a plot/readout.
///
/// # Errors
///
/// Returns [`HydroError`] if the water properties are invalid or any speed in
/// the range is out of domain.
pub fn resistance_curve(
    hull: &Hull,
    range: &SpeedRange,
    water: &WaterProperties,
) -> Result<ResistanceCurve, HydroError> {
    water.validate()?;
    let speeds = range.speeds();
    if speeds.is_empty() {
        return Err(HydroError::BadRange {
            start: range.start_ms,
            end: range.end_ms,
            step: range.step_ms,
        });
    }
    let mut points = Vec::with_capacity(speeds.len());
    for v in speeds {
        points.push(hull.resistance_at(v, water)?);
    }
    Ok(ResistanceCurve {
        hull: *hull,
        water: *water,
        points,
    })
}

#[cfg(test)]
mod tests;

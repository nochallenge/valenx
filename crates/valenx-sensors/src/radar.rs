//! Simulated **radar / RF** ranging sensor (monostatic, point-target).
//!
//! A [`Radar`] models the classic **monostatic radar range equation** for a
//! single point target in free space: it takes the target's position and radial
//! velocity plus a radar cross-section (RCS), and reports the received power as a
//! [`RadarReturn`] — slant range, signal-to-noise ratio (SNR), Doppler shift,
//! and a threshold detection decision. The pieces are the textbook closed forms:
//!
//! * **Range equation** — received power
//!   `Pr = Pt·G²·λ²·σ / ((4π)³·R⁴)`, so received power falls as `1/R⁴`
//!   (doubling the range drops `Pr` by exactly `12 dB`). Here `Pt` is the
//!   transmit power (W), `G` the antenna gain *toward the target* (linear), `λ`
//!   the wavelength (m), `σ` the target RCS (m²), and `R` the slant range (m).
//! * **RCS** ([`Rcs`]) — closed-form optical-region cross-sections for a few
//!   canonical shapes: a sphere (`σ = π·r²`), a flat plate at normal incidence
//!   (`σ = 4π·A²/λ²`), a dihedral / right-angle corner reflector
//!   (`σ = 8π·a⁴/λ²`), and a long cylinder broadside
//!   (`σ = 2π·r·L²/λ`).
//! * **Antenna gain** ([`AntennaPattern`]) — a simple main-lobe model whose
//!   gain is maximal at boresight and falls off symmetrically with the
//!   off-boresight angle (a Gaussian beam matched to a half-power beamwidth, or
//!   the textbook `sinc²` aperture pattern).
//! * **Doppler** — the radial-velocity shift `fd = −2·vr/λ` (a closing target,
//!   `vr < 0` by the sign convention below, gives a positive `fd`).
//! * **Detection** — the SNR is formed against a thermal-noise floor
//!   (`kTB` scaled by the receiver noise figure) and compared to a configured
//!   threshold; `SNR ≥ threshold` is a detection. Because `Pr ∝ R⁻⁴`, SNR is
//!   strictly decreasing in range, so detection holds out to a single maximum
//!   range and not beyond.
//!
//! ## Sign / frame conventions
//!
//! * **Range** `R = ‖target_position‖` is the slant range from the radar at the
//!   origin to the target (m), always `> 0` for a real target.
//! * **Boresight** is the radar's `+x` axis; the off-boresight angle of a target
//!   is the angle between `+x` and the target bearing.
//! * **Radial velocity** `vr` is the component of the target velocity *along the
//!   line of sight from radar to target* (`vr = v · r̂`). A **receding** target
//!   has `vr > 0`; a **closing** target has `vr < 0`. With `fd = −2·vr/λ` a
//!   closing target therefore produces a **positive** Doppler shift, matching the
//!   usual "approaching ⇒ frequency up" convention.
//!
//! ## Determinism — seeded thermal noise
//!
//! Like the rest of the crate, the optional receiver thermal noise is drawn from
//! the in-crate seeded [`crate::SplitMix64`] (no `rand` dependency), so a given
//! seed reproduces the same returns bit-for-bit. With `noise_std_db == 0` the SNR
//! is the exact deterministic value the tests pin.
//!
//! ## Honesty / scope caveats
//!
//! This is an **analytic radar-equation, point-target, free-space model** — the
//! same fidelity as the rest of `valenx-sensors`: it reproduces the *geometry and
//! first-order link budget* of a radar, exactly what an autonomy / tracking
//! pipeline needs to be built and V&V'd in the loop, but it is deliberately not a
//! high-fidelity electromagnetics model. In particular it does **not** model:
//!
//! * **clutter** (ground / sea / rain returns) or any background other than
//!   thermal noise;
//! * **multipath**, ducting, or terrain masking;
//! * **atmospheric / propagation loss** (the `1/R⁴` free-space term only);
//! * **extended-target effects** — scintillation / glint, the Swerling
//!   fluctuation models, or range/Doppler spread (the target is a single point
//!   scatterer with a fixed RCS);
//! * **electronic countermeasures** (jamming, chaff, deception) or any EW;
//! * **waveform / processing detail** — no pulse compression, CFAR, MTI/STAP,
//!   ambiguity function, range-gate or Doppler-bin structure (the "detection" is
//!   a single SNR-vs-threshold test);
//! * antenna **sidelobes / back-lobes** beyond the single main-lobe model, and no
//!   polarisation.
//!
//! The RCS formulae are the standard **optical-region** (high-frequency)
//! approximations and are only valid when the body is large compared to the
//! wavelength; they do not cover the Rayleigh or resonance regions. **DEFENSIVE
//! sensing only** — this models detection / ranging / tracking *input* (defense
//! thrust **M8 sensor-RF**); it is not weapons cueing or targeting-for-lethality.
//! Calibrating against a specific real radar is out of scope; this is the
//! research/educational-grade, reproducible testbed such work would build on.

use std::f64::consts::PI;

use nalgebra::Vector3;

use crate::error::SensorError;
use crate::rng::SplitMix64;

/// Boltzmann constant `k` (J/K), for the thermal-noise floor `kTB`.
const BOLTZMANN: f64 = 1.380_649e-23;
/// Standard reference temperature `T₀` (K) for receiver noise (290 K, IEEE).
const REFERENCE_TEMP_K: f64 = 290.0;
/// Speed of light in vacuum (m/s), for the frequency ↔ wavelength conversion.
const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// A small positive tolerance for "near-zero" guards.
const EPS: f64 = 1e-12;

/// Convert a linear power ratio to decibels (`10·log₁₀`).
///
/// A non-positive ratio has no finite dB value; this returns `f64::NEG_INFINITY`
/// for `ratio == 0` (and for any `ratio ≤ 0`) rather than a `NaN`, so a
/// zero-power return reads as "infinitely far below" the floor instead of
/// poisoning downstream arithmetic.
#[must_use]
pub fn linear_to_db(ratio: f64) -> f64 {
    if ratio > 0.0 {
        10.0 * ratio.log10()
    } else {
        f64::NEG_INFINITY
    }
}

/// Convert decibels back to a linear power ratio (`10^(db/10)`).
#[must_use]
pub fn db_to_linear(db: f64) -> f64 {
    10.0_f64.powf(db / 10.0)
}

/// Wavelength `λ = c / f` (m) for a frequency `f` (Hz).
///
/// # Errors
/// Returns [`SensorError::InvalidConfig`] if `frequency_hz` is not finite and
/// strictly positive.
pub fn wavelength_from_frequency(frequency_hz: f64) -> Result<f64, SensorError> {
    if !(frequency_hz.is_finite() && frequency_hz > 0.0) {
        return Err(SensorError::InvalidConfig(format!(
            "frequency must be finite and > 0 Hz, got {frequency_hz}"
        )));
    }
    Ok(SPEED_OF_LIGHT / frequency_hz)
}

/// A canonical shape whose **optical-region** radar cross-section has a textbook
/// closed form. Build one and read its [`Rcs::sigma`] (m²).
///
/// These are the high-frequency (body ≫ wavelength) approximations; the sphere
/// is wavelength-independent, the others scale with `λ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Rcs {
    /// A conducting **sphere** of `radius` (m) in the optical region:
    /// `σ = π·r²` (independent of wavelength).
    Sphere {
        /// Sphere radius (m).
        radius: f64,
    },
    /// A **flat plate** of area `area` (m²) at normal incidence:
    /// `σ = 4π·A²/λ²`.
    FlatPlate {
        /// Plate area (m²).
        area: f64,
        /// Wavelength (m).
        wavelength: f64,
    },
    /// A **right-angle dihedral corner reflector** with square faces of side
    /// `side` (m), at its peak (boresight) response: `σ = 8π·a⁴/λ²`.
    CornerReflector {
        /// Face side length `a` (m).
        side: f64,
        /// Wavelength (m).
        wavelength: f64,
    },
    /// A long **circular cylinder** of `radius` (m) and length `length` (m)
    /// illuminated broadside: `σ = 2π·r·L²/λ`.
    Cylinder {
        /// Cylinder radius (m).
        radius: f64,
        /// Cylinder length (m).
        length: f64,
        /// Wavelength (m).
        wavelength: f64,
    },
}

impl Rcs {
    /// The radar cross-section `σ` (m²) of this shape.
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if any dimension is not finite and
    /// strictly positive, or a wavelength is not finite and `> 0`.
    pub fn sigma(&self) -> Result<f64, SensorError> {
        match *self {
            Rcs::Sphere { radius } => {
                require_positive("sphere radius", radius)?;
                Ok(PI * radius * radius)
            }
            Rcs::FlatPlate { area, wavelength } => {
                require_positive("plate area", area)?;
                require_positive("wavelength", wavelength)?;
                Ok(4.0 * PI * area * area / (wavelength * wavelength))
            }
            Rcs::CornerReflector { side, wavelength } => {
                require_positive("corner side", side)?;
                require_positive("wavelength", wavelength)?;
                let a4 = side * side * side * side;
                Ok(8.0 * PI * a4 / (wavelength * wavelength))
            }
            Rcs::Cylinder {
                radius,
                length,
                wavelength,
            } => {
                require_positive("cylinder radius", radius)?;
                require_positive("cylinder length", length)?;
                require_positive("wavelength", wavelength)?;
                Ok(2.0 * PI * radius * length * length / wavelength)
            }
        }
    }
}

/// Validate that `value` is finite and strictly positive, naming the field.
fn require_positive(name: &str, value: f64) -> Result<(), SensorError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SensorError::InvalidConfig(format!(
            "{name} must be finite and > 0, got {value}"
        )))
    }
}

/// A simple antenna **main-lobe** gain pattern: maximal at boresight, falling off
/// symmetrically with the off-boresight angle.
///
/// Two shapes are offered, both parameterised by a peak (boresight) gain and a
/// half-power (−3 dB) beamwidth:
///
/// * [`AntennaPattern::Gaussian`] — a Gaussian main lobe
///   `G(θ) = G₀·exp(−4 ln2 · (θ/θ₃ ᵈᴮ)²)`, which is exactly `−3 dB` at
///   `θ = ±θ₃ ᵈᴮ/2`. Smooth, no nulls or sidelobes.
/// * [`AntennaPattern::SincSquared`] — the textbook uniform-aperture pattern
///   `G(θ) = G₀·sinc²(c·θ/θ₃ ᵈᴮ)` (with `sinc x = sin(πx)/(πx)`), which has the
///   characteristic main-lobe nulls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AntennaPattern {
    /// Gaussian main lobe.
    Gaussian {
        /// Peak (boresight) gain, **linear** (≥ 1, dimensionless).
        peak_gain: f64,
        /// Half-power (−3 dB) full beamwidth `θ₃ ᵈᴮ` (rad, > 0).
        beamwidth: f64,
    },
    /// Uniform-aperture `sinc²` main lobe.
    SincSquared {
        /// Peak (boresight) gain, **linear** (≥ 1, dimensionless).
        peak_gain: f64,
        /// Half-power (−3 dB) full beamwidth `θ₃ ᵈᴮ` (rad, > 0).
        beamwidth: f64,
    },
}

impl AntennaPattern {
    /// The peak (boresight) gain (linear).
    #[must_use]
    pub fn peak_gain(&self) -> f64 {
        match *self {
            AntennaPattern::Gaussian { peak_gain, .. }
            | AntennaPattern::SincSquared { peak_gain, .. } => peak_gain,
        }
    }

    /// The half-power beamwidth `θ₃ ᵈᴮ` (rad).
    #[must_use]
    pub fn beamwidth(&self) -> f64 {
        match *self {
            AntennaPattern::Gaussian { beamwidth, .. }
            | AntennaPattern::SincSquared { beamwidth, .. } => beamwidth,
        }
    }

    /// Validate the pattern's parameters.
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if the peak gain is not finite and
    /// `≥ 1`, or the beamwidth is not finite and `> 0`.
    pub fn validate(&self) -> Result<(), SensorError> {
        let g = self.peak_gain();
        if !(g.is_finite() && g >= 1.0) {
            return Err(SensorError::InvalidConfig(format!(
                "antenna peak gain must be finite and ≥ 1 (linear), got {g}"
            )));
        }
        let bw = self.beamwidth();
        if !(bw.is_finite() && bw > 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "antenna beamwidth must be finite and > 0 rad, got {bw}"
            )));
        }
        Ok(())
    }

    /// Linear gain at off-boresight angle `theta` (rad).
    ///
    /// The pattern is symmetric in `θ` (only `|θ|` matters), maximal at
    /// `θ = 0` (returns the peak gain), and monotonically decreasing across the
    /// main lobe. A non-finite `theta` yields `0` (no gain) rather than a `NaN`.
    #[must_use]
    pub fn gain(&self, theta: f64) -> f64 {
        if !theta.is_finite() {
            return 0.0;
        }
        match *self {
            AntennaPattern::Gaussian {
                peak_gain,
                beamwidth,
            } => {
                if beamwidth <= 0.0 {
                    return peak_gain;
                }
                let x = theta / beamwidth;
                // −3 dB (factor 1/2) at θ = ±beamwidth/2 ⇒ exponent 4·ln2·x².
                peak_gain * (-4.0 * std::f64::consts::LN_2 * x * x).exp()
            }
            AntennaPattern::SincSquared {
                peak_gain,
                beamwidth,
            } => {
                if beamwidth <= 0.0 {
                    return peak_gain;
                }
                // c chosen so sinc²(c/2) = 1/2 at θ = ±beamwidth/2: c ≈ 0.8859.
                const C: f64 = 0.885_894_4;
                let arg = C * theta / beamwidth;
                peak_gain * sinc(arg) * sinc(arg)
            }
        }
    }
}

/// The normalised cardinal sine `sinc x = sin(πx)/(πx)`, with `sinc 0 = 1`.
fn sinc(x: f64) -> f64 {
    if x.abs() < EPS {
        1.0
    } else {
        let pix = PI * x;
        pix.sin() / pix
    }
}

/// Doppler frequency shift `fd` (Hz) for a target with radial velocity `vr`
/// (m/s, **receding positive**) at wavelength `lambda` (m).
///
/// Uses `fd = −2·vr/λ`, so a **closing** target (`vr < 0`) yields a **positive**
/// `fd` ("approaching ⇒ frequency up"). The factor of two is the two-way
/// (transmit + receive) path of a monostatic radar.
///
/// # Errors
/// Returns [`SensorError::InvalidConfig`] if `lambda` is not finite and `> 0`, or
/// [`SensorError::NonFinite`] if `vr` is not finite.
pub fn doppler_shift(vr: f64, lambda: f64) -> Result<f64, SensorError> {
    require_positive("wavelength", lambda)?;
    if !vr.is_finite() {
        return Err(SensorError::NonFinite(format!(
            "radial velocity must be finite, got {vr}"
        )));
    }
    Ok(-2.0 * vr / lambda)
}

/// Configuration for a [`Radar`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadarConfig {
    /// Transmit power `Pt` (W, > 0).
    pub tx_power_w: f64,
    /// The antenna main-lobe gain pattern (its peak gain is `G` at boresight).
    pub antenna: AntennaPattern,
    /// Carrier wavelength `λ` (m, > 0). Build from a frequency with
    /// [`wavelength_from_frequency`].
    pub wavelength: f64,
    /// Receiver **noise bandwidth** `B` (Hz, > 0), used in the `kTB` noise floor.
    pub bandwidth_hz: f64,
    /// Receiver **noise figure** `F` (dB, ≥ 0): the excess noise above the ideal
    /// `kTB₀` floor.
    pub noise_figure_db: f64,
    /// Detection threshold on SNR (dB). `SNR ≥ threshold` is a detection.
    pub detection_threshold_db: f64,
    /// Standard deviation of additive zero-mean Gaussian noise on the SNR (dB,
    /// ≥ 0), drawn from the seeded PRNG. `0` makes the SNR deterministic.
    pub noise_std_db: f64,
}

impl Default for RadarConfig {
    /// A modest X-band-ish surveillance radar: 1 kW transmit, a 30 dB (linear
    /// 1000) Gaussian pencil beam 3° wide, 10 GHz carrier (`λ ≈ 3 cm`), 1 MHz
    /// noise bandwidth, 3 dB noise figure, a 13 dB detection threshold, no SNR
    /// noise.
    fn default() -> Self {
        Self {
            tx_power_w: 1_000.0,
            antenna: AntennaPattern::Gaussian {
                peak_gain: 1_000.0,
                beamwidth: 3.0_f64.to_radians(),
            },
            wavelength: SPEED_OF_LIGHT / 10.0e9,
            bandwidth_hz: 1.0e6,
            noise_figure_db: 3.0,
            detection_threshold_db: 13.0,
            noise_std_db: 0.0,
        }
    }
}

/// A monostatic radar measurement of a single target.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadarReturn {
    /// Slant range to the target (m).
    pub range: f64,
    /// Off-boresight angle of the target (rad).
    pub bearing: f64,
    /// Received power at the antenna (W), before the noise floor.
    pub received_power_w: f64,
    /// Signal-to-noise ratio (dB), including any seeded SNR noise.
    pub snr_db: f64,
    /// Doppler frequency shift (Hz); positive for a closing target.
    pub doppler_hz: f64,
    /// Whether the (possibly noisy) SNR met the detection threshold.
    pub detected: bool,
}

/// A simulated monostatic radar at the origin looking down its `+x` boresight.
#[derive(Debug, Clone)]
pub struct Radar {
    config: RadarConfig,
    /// Precomputed thermal-noise floor `k·T₀·B·F` (W).
    noise_floor_w: f64,
    rng: SplitMix64,
}

impl Radar {
    /// Build a radar from a config and a noise seed, validating the config.
    ///
    /// # Errors
    /// - [`SensorError::InvalidConfig`] if the transmit power, wavelength, or
    ///   noise bandwidth is not finite and `> 0`, the noise figure is not finite
    ///   or negative, the detection threshold is non-finite, or the antenna
    ///   pattern is invalid (see [`AntennaPattern::validate`]).
    /// - [`SensorError::InvalidNoise`] if `noise_std_db` is negative or
    ///   non-finite.
    pub fn new(config: RadarConfig, seed: u64) -> Result<Self, SensorError> {
        require_positive("tx_power_w", config.tx_power_w)?;
        require_positive("wavelength", config.wavelength)?;
        require_positive("bandwidth_hz", config.bandwidth_hz)?;
        config.antenna.validate()?;
        if !(config.noise_figure_db.is_finite() && config.noise_figure_db >= 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "noise_figure_db must be finite and ≥ 0, got {}",
                config.noise_figure_db
            )));
        }
        if !config.detection_threshold_db.is_finite() {
            return Err(SensorError::InvalidConfig(format!(
                "detection_threshold_db must be finite, got {}",
                config.detection_threshold_db
            )));
        }
        if !(config.noise_std_db.is_finite() && config.noise_std_db >= 0.0) {
            return Err(SensorError::InvalidNoise(format!(
                "noise_std_db must be finite and ≥ 0, got {}",
                config.noise_std_db
            )));
        }

        // Thermal-noise floor kT₀B scaled by the (linear) noise figure.
        let noise_floor_w = BOLTZMANN
            * REFERENCE_TEMP_K
            * config.bandwidth_hz
            * db_to_linear(config.noise_figure_db);

        Ok(Self {
            config,
            noise_floor_w,
            rng: SplitMix64::new(seed),
        })
    }

    /// The configuration this radar was built with.
    #[must_use]
    pub fn config(&self) -> &RadarConfig {
        &self.config
    }

    /// The thermal-noise floor `k·T₀·B·F` (W).
    #[must_use]
    pub fn noise_floor_w(&self) -> f64 {
        self.noise_floor_w
    }

    /// Received power `Pr` (W) from the monostatic range equation, for a target
    /// at slant range `range` (m) seen through antenna gain `gain` (linear) with
    /// cross-section `rcs` (m²):
    /// `Pr = Pt·G²·λ²·σ / ((4π)³·R⁴)`.
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if `range`, `gain`, or `rcs` is not
    /// finite and `> 0` (a real target is in front of the radar with positive
    /// gain and a positive cross-section). Guarding `range > 0` makes the `R⁴`
    /// divide safe.
    pub fn received_power(&self, range: f64, gain: f64, rcs: f64) -> Result<f64, SensorError> {
        require_positive("range", range)?;
        require_positive("gain", gain)?;
        require_positive("rcs", rcs)?;
        let four_pi_cubed = (4.0 * PI).powi(3);
        let num = self.config.tx_power_w
            * gain
            * gain
            * self.config.wavelength
            * self.config.wavelength
            * rcs;
        let den = four_pi_cubed * range.powi(4);
        Ok(num / den)
    }

    /// The clean (noise-free) SNR in dB for a received power `received_power_w`
    /// (W) against this radar's thermal-noise floor.
    #[must_use]
    pub fn snr_db(&self, received_power_w: f64) -> f64 {
        linear_to_db(received_power_w / self.noise_floor_w)
    }

    /// Measure a single target at `target_position` (m, radar-frame, boresight
    /// `+x`) moving at `target_velocity` (m/s) with radar cross-section `rcs`
    /// (m²).
    ///
    /// Returns the full [`RadarReturn`] (range, bearing, received power, SNR,
    /// Doppler, detection) or `None` if the target is at (or essentially at) the
    /// radar origin so the range / line-of-sight is undefined.
    ///
    /// The antenna gain is evaluated at the target's off-boresight angle, the
    /// received power from the range equation, the SNR against the noise floor
    /// (plus optional seeded dB noise), the Doppler from the line-of-sight
    /// velocity component, and `detected` is `SNR ≥ detection_threshold_db`.
    ///
    /// # Errors
    /// Returns [`SensorError::NonFinite`] if any component of `target_position`
    /// or `target_velocity`, or `rcs`, is not finite. (A non-positive or
    /// non-finite `rcs` is rejected; a degenerate position returns `None`.)
    pub fn measure(
        &mut self,
        target_position: Vector3<f64>,
        target_velocity: Vector3<f64>,
        rcs: f64,
    ) -> Result<Option<RadarReturn>, SensorError> {
        if !target_position.iter().all(|c| c.is_finite()) {
            return Err(SensorError::NonFinite(
                "target_position must be finite".into(),
            ));
        }
        if !target_velocity.iter().all(|c| c.is_finite()) {
            return Err(SensorError::NonFinite(
                "target_velocity must be finite".into(),
            ));
        }
        require_positive("rcs", rcs)?;

        let range = target_position.norm();
        if range < EPS {
            // Target at the origin: range and bearing undefined.
            return Ok(None);
        }
        let los = target_position / range; // unit line of sight, radar → target

        // Off-boresight angle: angle between +x and the line of sight.
        let cos_theta = los.x.clamp(-1.0, 1.0);
        let bearing = cos_theta.acos();

        // Antenna gain toward the target, then the range equation.
        let gain = self.config.antenna.gain(bearing);
        // gain can be ~0 far out on the pattern; guard the SNR rather than the
        // power equation (received_power requires gain > 0).
        let received_power_w = if gain > 0.0 {
            self.received_power(range, gain, rcs)?
        } else {
            0.0
        };

        // SNR against the thermal-noise floor, plus optional seeded dB noise.
        let clean_snr_db = self.snr_db(received_power_w);
        let snr_db = if self.config.noise_std_db > 0.0 && clean_snr_db.is_finite() {
            clean_snr_db + self.rng.next_normal(0.0, self.config.noise_std_db)
        } else {
            clean_snr_db
        };

        // Doppler from the radial (line-of-sight) velocity component.
        let vr = target_velocity.dot(&los);
        let doppler_hz = doppler_shift(vr, self.config.wavelength)?;

        let detected = snr_db >= self.config.detection_threshold_db;

        Ok(Some(RadarReturn {
            range,
            bearing,
            received_power_w,
            snr_db,
            doppler_hz,
            detected,
        }))
    }

    /// The **maximum detection range** (m) for a target of cross-section `rcs`
    /// (m²) on boresight, i.e. the range at which the boresight SNR equals the
    /// detection threshold. Beyond this range the target is not detected; within
    /// it, it is (the SNR is strictly decreasing in `R`).
    ///
    /// Inverts the range equation at the threshold:
    /// `R_max = [ Pt·G²·λ²·σ / ((4π)³ · N · SNR_thresh) ]^(1/4)`.
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if `rcs` is not finite and `> 0`.
    pub fn max_detection_range(&self, rcs: f64) -> Result<f64, SensorError> {
        require_positive("rcs", rcs)?;
        let gain = self.config.antenna.peak_gain();
        let four_pi_cubed = (4.0 * PI).powi(3);
        let snr_thresh_linear = db_to_linear(self.config.detection_threshold_db);
        let num = self.config.tx_power_w
            * gain
            * gain
            * self.config.wavelength
            * self.config.wavelength
            * rcs;
        let den = four_pi_cubed * self.noise_floor_w * snr_thresh_linear;
        Ok((num / den).powf(0.25))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    // ---- BENCHMARK PIN (1): R⁴ law — doubling range drops Pr by 12 dB. ----

    #[test]
    fn range_equation_is_inverse_fourth_power() {
        let radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let g = radar.config().antenna.peak_gain();
        let p_near = radar.received_power(1_000.0, g, 1.0).unwrap();
        let p_far = radar.received_power(2_000.0, g, 1.0).unwrap();
        // Pr ∝ 1/R⁴: doubling R ⇒ ratio = 1/16 ⇒ −12.0412 dB exactly.
        let drop_db = linear_to_db(p_far / p_near);
        assert!(
            (drop_db - (-10.0 * 16.0_f64.log10())).abs() < 1e-9,
            "drop = {drop_db} dB (expected ≈ −12.04)"
        );
        assert!(
            (p_near / p_far - 16.0).abs() < 1e-6,
            "ratio {}",
            p_near / p_far
        );
    }

    #[test]
    fn range_equation_matches_closed_form() {
        // Hand-evaluate Pr = Pt·G²·λ²·σ / ((4π)³·R⁴) for simple round numbers.
        let cfg = RadarConfig {
            tx_power_w: 1.0,
            antenna: AntennaPattern::Gaussian {
                peak_gain: 1.0,
                beamwidth: 1.0,
            },
            wavelength: 1.0,
            bandwidth_hz: 1.0,
            noise_figure_db: 0.0,
            detection_threshold_db: 0.0,
            noise_std_db: 0.0,
        };
        let radar = Radar::new(cfg, 0).unwrap();
        let pr = radar.received_power(1.0, 1.0, 1.0).unwrap();
        let expected = 1.0 / (4.0 * PI).powi(3);
        assert!(
            (pr - expected).abs() < 1e-18,
            "Pr = {pr}, expected {expected}"
        );
    }

    // ---- BENCHMARK PIN (2): RCS closed forms at textbook values. ----

    #[test]
    fn sphere_rcs_is_pi_r_squared() {
        let r = 2.0;
        let sigma = Rcs::Sphere { radius: r }.sigma().unwrap();
        assert!((sigma - PI * r * r).abs() < 1e-12, "σ = {sigma}");
        // And it is wavelength-independent: σ at r=1 is exactly π.
        let unit = Rcs::Sphere { radius: 1.0 }.sigma().unwrap();
        assert!((unit - PI).abs() < 1e-12);
    }

    #[test]
    fn flat_plate_rcs_matches_textbook() {
        // σ = 4π·A²/λ². Pick A = 1 m², λ = 0.03 m ⇒ σ = 4π/0.0009.
        let area = 1.0;
        let lambda = 0.03;
        let sigma = Rcs::FlatPlate {
            area,
            wavelength: lambda,
        }
        .sigma()
        .unwrap();
        let expected = 4.0 * PI * area * area / (lambda * lambda);
        assert!(
            (sigma - expected).abs() < 1e-6,
            "σ = {sigma}, expected {expected}"
        );
        // A 1 m square plate at X-band has a huge RCS (~10^4–10^5 m²).
        assert!(sigma > 1.0e4, "plate σ should be very large, got {sigma}");
    }

    #[test]
    fn corner_reflector_and_cylinder_match_closed_form() {
        let lambda = 0.03;
        let side = 0.2;
        let corner = Rcs::CornerReflector {
            side,
            wavelength: lambda,
        }
        .sigma()
        .unwrap();
        let corner_expected = 8.0 * PI * side.powi(4) / (lambda * lambda);
        assert!(
            (corner - corner_expected).abs() < 1e-9,
            "corner σ = {corner}"
        );

        let cyl = Rcs::Cylinder {
            radius: 0.1,
            length: 2.0,
            wavelength: lambda,
        }
        .sigma()
        .unwrap();
        let cyl_expected = 2.0 * PI * 0.1 * 2.0 * 2.0 / lambda;
        assert!((cyl - cyl_expected).abs() < 1e-9, "cylinder σ = {cyl}");
    }

    // ---- BENCHMARK PIN (3): Doppler fd = 2·vr/λ for a known closing speed. ----

    #[test]
    fn doppler_matches_two_vr_over_lambda() {
        let lambda = 0.03; // 10 GHz
                           // Closing at 300 m/s (vr negative by convention) ⇒ positive fd.
        let fd = doppler_shift(-300.0, lambda).unwrap();
        let expected = 2.0 * 300.0 / lambda; // = 20 kHz
        assert!(
            (fd - expected).abs() < 1e-6,
            "fd = {fd}, expected {expected}"
        );
        assert!(fd > 0.0, "closing target ⇒ positive Doppler");
        // Receding at the same speed ⇒ equal magnitude, opposite sign.
        let fd_recede = doppler_shift(300.0, lambda).unwrap();
        assert!((fd_recede + expected).abs() < 1e-6);
    }

    #[test]
    fn measure_doppler_uses_line_of_sight_component() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let lambda = radar.config().wavelength;
        // Target on +x boresight, closing straight in at 100 m/s (velocity −x).
        let ret = radar
            .measure(v(1_000.0, 0.0, 0.0), v(-100.0, 0.0, 0.0), 1.0)
            .unwrap()
            .unwrap();
        let expected = 2.0 * 100.0 / lambda;
        assert!(
            (ret.doppler_hz - expected).abs() < 1e-3,
            "fd = {}",
            ret.doppler_hz
        );

        // Pure cross-range velocity (along +y for an +x target) ⇒ zero Doppler.
        let cross = radar
            .measure(v(1_000.0, 0.0, 0.0), v(0.0, 100.0, 0.0), 1.0)
            .unwrap()
            .unwrap();
        assert!(
            cross.doppler_hz.abs() < 1e-6,
            "cross fd = {}",
            cross.doppler_hz
        );
    }

    // ---- BENCHMARK PIN (4): detection crosses threshold at the expected max
    //      range and is monotone in R. ----

    #[test]
    fn detection_holds_within_max_range_and_fails_beyond() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let rcs = 1.0;
        let r_max = radar.max_detection_range(rcs).unwrap();
        assert!(r_max.is_finite() && r_max > 0.0, "r_max = {r_max}");

        // Just inside: detected. Just outside: not.
        let inside = radar
            .measure(v(r_max * 0.99, 0.0, 0.0), v(0.0, 0.0, 0.0), rcs)
            .unwrap()
            .unwrap();
        let outside = radar
            .measure(v(r_max * 1.01, 0.0, 0.0), v(0.0, 0.0, 0.0), rcs)
            .unwrap()
            .unwrap();
        assert!(
            inside.detected,
            "should detect inside r_max (SNR {})",
            inside.snr_db
        );
        assert!(
            !outside.detected,
            "should not detect beyond r_max (SNR {})",
            outside.snr_db
        );

        // At exactly r_max the SNR equals the threshold (to floating point).
        let at = radar
            .measure(v(r_max, 0.0, 0.0), v(0.0, 0.0, 0.0), rcs)
            .unwrap()
            .unwrap();
        assert!(
            (at.snr_db - radar.config().detection_threshold_db).abs() < 1e-6,
            "SNR at r_max = {} vs threshold {}",
            at.snr_db,
            radar.config().detection_threshold_db
        );
    }

    #[test]
    fn snr_is_monotonically_decreasing_in_range() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let mut last = f64::INFINITY;
        for i in 1..=50 {
            let r = i as f64 * 1_000.0;
            let ret = radar
                .measure(v(r, 0.0, 0.0), v(0.0, 0.0, 0.0), 1.0)
                .unwrap()
                .unwrap();
            assert!(ret.snr_db < last, "SNR not decreasing at R = {r}");
            last = ret.snr_db;
        }
    }

    // ---- BENCHMARK PIN (5): antenna gain maximal at boresight and symmetric. ----

    #[test]
    fn antenna_gain_is_peak_at_boresight_and_symmetric() {
        for pat in [
            AntennaPattern::Gaussian {
                peak_gain: 1_000.0,
                beamwidth: 0.05,
            },
            AntennaPattern::SincSquared {
                peak_gain: 1_000.0,
                beamwidth: 0.05,
            },
        ] {
            let g0 = pat.gain(0.0);
            assert!((g0 - pat.peak_gain()).abs() < 1e-9, "boresight gain");
            // Symmetric: g(+θ) == g(−θ).
            for theta in [0.005, 0.01, 0.02, 0.03] {
                let gp = pat.gain(theta);
                let gm = pat.gain(-theta);
                assert!((gp - gm).abs() < 1e-9, "asymmetric at {theta}");
                // Off boresight is strictly less than the peak.
                assert!(gp < g0, "gain not falling off at {theta}");
            }
        }
    }

    #[test]
    fn gaussian_beamwidth_is_minus_three_db_at_half_beamwidth() {
        let pat = AntennaPattern::Gaussian {
            peak_gain: 100.0,
            beamwidth: 0.1,
        };
        // At ±beamwidth/2 the gain should be exactly half (−3.0103 dB).
        let g_half = pat.gain(0.05);
        assert!(
            (g_half / pat.peak_gain() - 0.5).abs() < 1e-9,
            "ratio = {}",
            g_half / pat.peak_gain()
        );
    }

    #[test]
    fn sinc_squared_beamwidth_is_minus_three_db_at_half_beamwidth() {
        let pat = AntennaPattern::SincSquared {
            peak_gain: 100.0,
            beamwidth: 0.1,
        };
        let g_half = pat.gain(0.05);
        // The C constant is tuned for −3 dB at the half-beamwidth; allow a small
        // band for the approximation.
        assert!(
            (g_half / pat.peak_gain() - 0.5).abs() < 1e-3,
            "ratio = {}",
            g_half / pat.peak_gain()
        );
    }

    // ---- fail-loud on bad config / NaN, guarded divides. ----

    #[test]
    fn invalid_configs_are_rejected() {
        // Each case perturbs one field of an otherwise-valid default config.
        let bad = [
            // Pt ≤ 0.
            RadarConfig {
                tx_power_w: 0.0,
                ..Default::default()
            },
            // λ ≤ 0.
            RadarConfig {
                wavelength: -1.0,
                ..Default::default()
            },
            // Bandwidth ≤ 0.
            RadarConfig {
                bandwidth_hz: 0.0,
                ..Default::default()
            },
            // Negative noise figure.
            RadarConfig {
                noise_figure_db: -1.0,
                ..Default::default()
            },
            // Negative SNR noise.
            RadarConfig {
                noise_std_db: -0.1,
                ..Default::default()
            },
            // NaN threshold.
            RadarConfig {
                detection_threshold_db: f64::NAN,
                ..Default::default()
            },
            // Bad antenna: gain < 1.
            RadarConfig {
                antenna: AntennaPattern::Gaussian {
                    peak_gain: 0.5,
                    beamwidth: 0.1,
                },
                ..Default::default()
            },
            // Bad antenna: zero beamwidth.
            RadarConfig {
                antenna: AntennaPattern::Gaussian {
                    peak_gain: 10.0,
                    beamwidth: 0.0,
                },
                ..Default::default()
            },
        ];
        for (i, cfg) in bad.into_iter().enumerate() {
            assert!(
                Radar::new(cfg, 0).is_err(),
                "config #{i} should be rejected"
            );
        }
    }

    #[test]
    fn negative_or_nonfinite_rcs_is_rejected() {
        assert!(Rcs::Sphere { radius: -1.0 }.sigma().is_err());
        assert!(Rcs::Sphere { radius: f64::NAN }.sigma().is_err());
        assert!(Rcs::FlatPlate {
            area: 1.0,
            wavelength: 0.0
        }
        .sigma()
        .is_err());

        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        // Negative RCS into received_power and measure.
        assert!(radar.received_power(100.0, 10.0, -1.0).is_err());
        assert!(radar
            .measure(v(100.0, 0.0, 0.0), v(0.0, 0.0, 0.0), -1.0)
            .is_err());
    }

    #[test]
    fn nonfinite_inputs_are_rejected() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        assert!(radar
            .measure(v(f64::NAN, 0.0, 0.0), v(0.0, 0.0, 0.0), 1.0)
            .is_err());
        assert!(radar
            .measure(v(100.0, 0.0, 0.0), v(f64::INFINITY, 0.0, 0.0), 1.0)
            .is_err());
        assert!(doppler_shift(f64::NAN, 0.03).is_err());
        assert!(doppler_shift(100.0, 0.0).is_err());
        assert!(wavelength_from_frequency(0.0).is_err());
        assert!(wavelength_from_frequency(-1.0).is_err());
    }

    #[test]
    fn target_at_origin_returns_none() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let ret = radar
            .measure(v(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0), 1.0)
            .unwrap();
        assert!(ret.is_none(), "degenerate range should be a no-measure");
    }

    // ---- determinism, integration, conversions. ----

    #[test]
    fn snr_noise_is_deterministic_for_a_seed() {
        let cfg = RadarConfig {
            noise_std_db: 1.5,
            ..Default::default()
        };
        let mut a = Radar::new(cfg, 99).unwrap();
        let mut b = Radar::new(cfg, 99).unwrap();
        let ra = a
            .measure(v(5_000.0, 0.0, 0.0), v(-50.0, 0.0, 0.0), 1.0)
            .unwrap();
        let rb = b
            .measure(v(5_000.0, 0.0, 0.0), v(-50.0, 0.0, 0.0), 1.0)
            .unwrap();
        assert_eq!(ra, rb, "same seed ⇒ identical return");
    }

    #[test]
    fn off_boresight_target_has_lower_received_power_than_on_boresight() {
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        // Same range, one on boresight, one well off it.
        let on = radar
            .measure(v(2_000.0, 0.0, 0.0), v(0.0, 0.0, 0.0), 1.0)
            .unwrap()
            .unwrap();
        // ~5° off boresight (well outside the 3° beam).
        let off_angle = 5.0_f64.to_radians();
        let pos = v(2_000.0 * off_angle.cos(), 2_000.0 * off_angle.sin(), 0.0);
        let off = radar.measure(pos, v(0.0, 0.0, 0.0), 1.0).unwrap().unwrap();
        assert!((on.range - off.range).abs() < 1e-6, "ranges should match");
        assert!(
            off.received_power_w < on.received_power_w,
            "off-boresight should receive less power"
        );
        assert!(off.bearing > on.bearing);
    }

    #[test]
    fn wavelength_frequency_roundtrip() {
        let lambda = wavelength_from_frequency(10.0e9).unwrap();
        assert!((lambda - SPEED_OF_LIGHT / 10.0e9).abs() < 1e-15);
        // ~3 cm at X-band.
        assert!((lambda - 0.0299792458).abs() < 1e-9, "λ = {lambda}");
    }

    #[test]
    fn db_conversions_roundtrip() {
        for &lin in &[1.0, 2.0, 10.0, 1_000.0, 0.5] {
            let db = linear_to_db(lin);
            let back = db_to_linear(db);
            assert!((back - lin).abs() < 1e-9, "roundtrip {lin} → {db} → {back}");
        }
        // 10× is +10 dB, 2× is ~+3.0103 dB.
        assert!((linear_to_db(10.0) - 10.0).abs() < 1e-12);
        assert!((linear_to_db(2.0) - 3.010_299_956_639_812).abs() < 1e-9);
        // Non-positive ratio → −∞, not NaN.
        assert_eq!(linear_to_db(0.0), f64::NEG_INFINITY);
        assert_eq!(linear_to_db(-1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn end_to_end_realistic_detection() {
        // Default X-band radar vs a 1 m² target: should comfortably detect a
        // close, on-boresight target and report a sane SNR and Doppler.
        let mut radar = Radar::new(RadarConfig::default(), 0).unwrap();
        let sigma = Rcs::Sphere { radius: 0.5641896 }.sigma().unwrap(); // ≈ 1 m²
        assert!((sigma - 1.0).abs() < 1e-3, "σ ≈ 1 m², got {sigma}");
        let ret = radar
            .measure(v(10_000.0, 0.0, 0.0), v(-200.0, 0.0, 0.0), sigma)
            .unwrap()
            .unwrap();
        assert!((ret.range - 10_000.0).abs() < 1e-6);
        assert!(ret.bearing.abs() < 1e-9, "on boresight");
        assert!(ret.snr_db.is_finite());
        assert!(ret.doppler_hz > 0.0, "closing ⇒ positive Doppler");
        assert!(ret.received_power_w > 0.0);
    }
}

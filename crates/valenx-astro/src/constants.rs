//! Physical constants used across the ascent simulator.
//!
//! Values are the standard geodetic / atmospheric references (WGS-84
//! Earth, US Standard Atmosphere 1976). They are `pub const` so callers
//! and tests can reference the exact same numbers the solver uses.

/// Earth gravitational parameter `μ = G·M_⊕` (m³/s²), the WGS-84 value.
pub const MU_EARTH: f64 = 3.986_004_418e14;

/// Earth second dynamic form factor (J2) — the dominant oblateness term
/// of the geopotential (dimensionless), WGS-84.
pub const J2_EARTH: f64 = 1.082_626_68e-3;

/// Earth equatorial radius (m), WGS-84. This is the semi-major axis `a` of
/// the WGS-84 reference ellipsoid (and the equatorial radius used by the
/// spherical-Earth models elsewhere in the crate).
pub const R_EARTH: f64 = 6_378_137.0;

/// WGS-84 ellipsoid flattening `f = (a − b)/a` (dimensionless) — the defining
/// shape parameter of the reference geoid. The polar semi-axis is
/// `b = a(1 − f)` ([`R_EARTH_POLAR`]). Used by the geodetic ↔ Earth-fixed
/// (ECEF) coordinate transforms in [`crate::frames`].
pub const WGS84_FLATTENING: f64 = 1.0 / 298.257_223_563;

/// WGS-84 polar (semi-minor) Earth radius `b = a(1 − f)` (m).
pub const R_EARTH_POLAR: f64 = R_EARTH * (1.0 - WGS84_FLATTENING);

/// First eccentricity *squared* of the WGS-84 ellipsoid,
/// `e² = f(2 − f) = 1 − (b/a)²` (dimensionless). The core quantity in the
/// geodetic-latitude ↔ ECEF conversion (the prime-vertical radius of
/// curvature `N = a/√(1 − e²·sin²φ)`).
pub const WGS84_ECC_SQ: f64 = WGS84_FLATTENING * (2.0 - WGS84_FLATTENING);

/// Earth sidereal rotation rate (rad/s).
pub const OMEGA_EARTH: f64 = 7.292_115_9e-5;

/// Earth's mean orbital angular rate about the Sun (rad/s) — `2π` per tropical
/// year (365.242190 days), ≈ 0.9856°/day. A **sun-synchronous** orbit's nodal
/// precession must match this so the orbit plane keeps a fixed Sun angle.
pub const EARTH_ORBITAL_RATE: f64 = 1.991_063e-7;

/// Standard gravitational acceleration used to convert specific impulse
/// (seconds) into effective exhaust velocity (m/s): `v_e = Isp · g₀`.
pub const G0: f64 = 9.806_65;

/// Specific gas constant for dry air (J/(kg·K)).
pub const R_AIR: f64 = 287.052_8;

/// Ratio of specific heats for air (dimensionless), used for the local
/// speed of sound `a = √(γ·R·T)`.
pub const GAMMA_AIR: f64 = 1.4;

/// Effective Earth radius used in the geopotential-altitude conversion
/// of the US Standard Atmosphere 1976 (m).
pub const ATMOS_EARTH_RADIUS: f64 = 6_356_766.0;

/// Top of the modelled atmosphere (geometric altitude, m). Above this
/// the density is treated as zero — drag is negligible there.
pub const ATMOS_TOP_M: f64 = 86_000.0;

//! Gerstner (trochoidal) directional wave field.
//!
//! The sea surface is modelled as a **sum of `N` directional Gerstner
//! (trochoidal) waves**. A single Gerstner wave does not merely raise and lower
//! the surface like a sine; it moves each surface particle in a *circle* (deep
//! water), so crests are sharpened and troughs are flattened — the classic
//! trochoidal profile. For a horizontal rest position `x = (x, z)` and time `t`,
//! one wave of amplitude `A`, wave vector `k_vec = k * d` (unit direction `d`,
//! wavenumber `k = 2π / L` for wavelength `L`), angular frequency `ω`, steepness
//! `Q ∈ [0, 1]`, and phase offset `φ` displaces the surface point to
//!
//! ```text
//! θ      = k_vec · x − ω t + φ
//! Δx     = (Q A) d.x cos θ
//! Δz     = (Q A) d.y cos θ      (d.y is the z-component of the horizontal dir)
//! Δheight=        A      sin θ
//! ```
//!
//! Summing `N` such waves gives a horizontally displaced point whose height is
//! the sum of the `A sin θ` terms. The **deep-water dispersion relation**
//! `ω = sqrt(g k)` ties frequency to wavenumber, so the phase speed of each
//! component is `c = ω / k = sqrt(g / k)` — longer waves travel faster.
//!
//! ## Steepness and the gimbal limit
//!
//! The horizontal displacement amplitude is `Q A` per component. If the summed
//! steepness is too high the surface self-intersects (the trochoid forms a
//! loop). This model **does not** clamp `Q` for you beyond requiring
//! `0 ≤ Q ≤ 1` per wave; choosing physically sane amplitudes/steepness is the
//! caller's responsibility, exactly as in a real-time ocean renderer.
//!
//! ## Surface normal
//!
//! The analytic normal of the Gerstner surface is obtained from the partial
//! derivatives of the displaced position with respect to the two horizontal rest
//! coordinates (Tessendorf 2001, "Simulating Ocean Water", eq. for the Gerstner
//! Jacobian). [`OceanWaveField::normal_at`] evaluates that closed form.

use crate::error::OceanError;
use nalgebra::Vector3;

/// Standard gravity, m/s². Used as the default for the dispersion relation.
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// Density of sea water at ~15 °C, kg/m³. Provided for convenience; the buoyancy
/// model takes the density explicitly so fresh water (≈ 1000) is equally valid.
pub const SEAWATER_DENSITY: f64 = 1025.0;

/// A single directional Gerstner (trochoidal) wave component.
///
/// Construct with [`GerstnerWave::new`], which validates the parameters and
/// derives the wavenumber and (deep-water) angular frequency. All fields are
/// read-only after construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GerstnerWave {
    /// Vertical amplitude `A` (m). Crest-to-trough height of this component
    /// alone is `2 A`.
    amplitude: f64,
    /// Wavelength `L` (m), the crest-to-crest spacing along the direction of
    /// travel.
    wavelength: f64,
    /// Steepness / phase-bunching factor `Q ∈ [0, 1]`. `0` is a pure (sinusoidal
    /// height, no horizontal motion) wave; `1` is the steepest non-self-
    /// intersecting single trochoid.
    steepness: f64,
    /// Phase offset `φ` (rad).
    phase: f64,
    /// Unit horizontal direction of travel `d = (dx, dz)`.
    direction: [f64; 2],
    /// Wavenumber `k = 2π / L` (rad/m). Cached; always `> 0`.
    k: f64,
    /// Angular frequency `ω = sqrt(g k)` (rad/s). Cached; always `> 0`.
    omega: f64,
}

impl GerstnerWave {
    /// Build a single Gerstner wave.
    ///
    /// * `amplitude` — vertical amplitude `A` (m), must be `> 0` and finite.
    /// * `wavelength` — `L` (m), must be `> 0` and finite.
    /// * `steepness` — `Q`, must lie in `[0, 1]`.
    /// * `phase` — phase offset `φ` (rad), must be finite.
    /// * `direction` — horizontal direction `(dx, dz)`; must be finite and have
    ///   non-zero length (it is normalised internally).
    /// * `gravity` — `g` (m/s²) for the dispersion relation, must be `> 0`.
    ///
    /// # Errors
    ///
    /// Returns [`OceanError::InvalidConfig`] if `amplitude`, `wavelength`, or
    /// `gravity` is non-positive, if `steepness ∉ [0, 1]`, or if `direction` is
    /// the zero vector; [`OceanError::NonFinite`] if any input is `NaN`/`±∞`.
    pub fn new(
        amplitude: f64,
        wavelength: f64,
        steepness: f64,
        phase: f64,
        direction: (f64, f64),
        gravity: f64,
    ) -> Result<Self, OceanError> {
        let (dx, dz) = direction;
        for (name, v) in [
            ("amplitude", amplitude),
            ("wavelength", wavelength),
            ("steepness", steepness),
            ("phase", phase),
            ("direction.x", dx),
            ("direction.z", dz),
            ("gravity", gravity),
        ] {
            if !v.is_finite() {
                return Err(OceanError::NonFinite(format!("{name} = {v}")));
            }
        }
        if amplitude <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "amplitude must be > 0, got {amplitude}"
            )));
        }
        if wavelength <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "wavelength must be > 0, got {wavelength}"
            )));
        }
        if gravity <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "gravity must be > 0, got {gravity}"
            )));
        }
        if !(0.0..=1.0).contains(&steepness) {
            return Err(OceanError::InvalidConfig(format!(
                "steepness must be in [0, 1], got {steepness}"
            )));
        }
        let len = (dx * dx + dz * dz).sqrt();
        if len <= 0.0 {
            return Err(OceanError::InvalidConfig(
                "direction must be a non-zero vector".to_string(),
            ));
        }
        // k > 0 because wavelength > 0; omega > 0 because g, k > 0.
        let k = 2.0 * std::f64::consts::PI / wavelength;
        let omega = (gravity * k).sqrt();
        Ok(Self {
            amplitude,
            wavelength,
            steepness,
            phase,
            direction: [dx / len, dz / len],
            k,
            omega,
        })
    }

    /// Vertical amplitude `A` (m).
    #[must_use]
    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }

    /// Wavelength `L` (m).
    #[must_use]
    pub fn wavelength(&self) -> f64 {
        self.wavelength
    }

    /// Steepness `Q`.
    #[must_use]
    pub fn steepness(&self) -> f64 {
        self.steepness
    }

    /// Wavenumber `k = 2π / L` (rad/m).
    #[must_use]
    pub fn wavenumber(&self) -> f64 {
        self.k
    }

    /// Angular frequency `ω = sqrt(g k)` (rad/s).
    #[must_use]
    pub fn angular_frequency(&self) -> f64 {
        self.omega
    }

    /// Temporal period `T = 2π / ω` (s).
    #[must_use]
    pub fn period(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.omega
    }

    /// Deep-water phase speed `c = ω / k = sqrt(g / k)` (m/s).
    ///
    /// `k > 0` by construction, so the divide is safe.
    #[must_use]
    pub fn phase_speed(&self) -> f64 {
        self.omega / self.k
    }

    /// The scalar phase `θ = k_vec · x − ω t + φ` at horizontal rest position
    /// `(x, z)` and time `t`.
    #[inline]
    fn theta(&self, x: f64, z: f64, t: f64) -> f64 {
        self.k * (self.direction[0] * x + self.direction[1] * z) - self.omega * t + self.phase
    }
}

/// A sum-of-`N`-Gerstner-waves ocean surface.
///
/// Holds the component waves and the gravity used by them. Evaluate the surface
/// with [`OceanWaveField::displacement_at`] (full 3-D offset of a rest point),
/// [`OceanWaveField::height_at`] (water height directly, ignoring the horizontal
/// trochoidal shift — the quantity a buoyancy probe at a fixed `(x, z)` needs),
/// and [`OceanWaveField::normal_at`] (analytic unit surface normal).
#[derive(Debug, Clone)]
pub struct OceanWaveField {
    waves: Vec<GerstnerWave>,
    /// Mean (still-water) sea level (m). Heights are reported relative to the
    /// world origin, i.e. `mean_level + Σ A sin θ`.
    mean_level: f64,
}

impl OceanWaveField {
    /// Build a wave field from a set of [`GerstnerWave`] components and a mean
    /// (still-water) level.
    ///
    /// An empty wave set is allowed and yields a perfectly flat sea at
    /// `mean_level` (useful as a degenerate baseline and for the
    /// zero-net-force-at-equilibrium benchmark).
    ///
    /// # Errors
    ///
    /// Returns [`OceanError::NonFinite`] if `mean_level` is not finite.
    pub fn new(waves: Vec<GerstnerWave>, mean_level: f64) -> Result<Self, OceanError> {
        if !mean_level.is_finite() {
            return Err(OceanError::NonFinite(format!("mean_level = {mean_level}")));
        }
        Ok(Self { waves, mean_level })
    }

    /// A flat sea at the given mean level (no waves).
    ///
    /// # Errors
    ///
    /// Returns [`OceanError::NonFinite`] if `mean_level` is not finite.
    pub fn flat(mean_level: f64) -> Result<Self, OceanError> {
        Self::new(Vec::new(), mean_level)
    }

    /// The component waves.
    #[must_use]
    pub fn waves(&self) -> &[GerstnerWave] {
        &self.waves
    }

    /// The mean (still-water) sea level (m).
    #[must_use]
    pub fn mean_level(&self) -> f64 {
        self.mean_level
    }

    /// The full Gerstner displacement of the rest point `(x, z)` at time `t`:
    /// the world-space offset `(Δx, Δheight, Δz)` to add to `(x, mean_level, z)`
    /// to get the displaced surface point.
    ///
    /// The `y` component is the height offset; `x`/`z` are the horizontal
    /// trochoidal shift. Returns a zero vector for a flat sea.
    #[must_use]
    pub fn displacement_at(&self, x: f64, z: f64, t: f64) -> Vector3<f64> {
        let mut d = Vector3::zeros();
        for w in &self.waves {
            let theta = w.theta(x, z, t);
            let (s, c) = theta.sin_cos();
            let qa = w.steepness * w.amplitude;
            d.x += qa * w.direction[0] * c;
            d.z += qa * w.direction[1] * c;
            d.y += w.amplitude * s;
        }
        d
    }

    /// The water **height** at the horizontal position `(x, z)` and time `t`,
    /// relative to the world origin — i.e. `mean_level + Σ A sin θ`.
    ///
    /// This evaluates the height contribution directly at the queried `(x, z)`
    /// (it does **not** invert the horizontal trochoidal shift). This is exactly
    /// the quantity a vertical buoyancy probe at a fixed horizontal location
    /// samples, and the one whose crest-to-trough amplitude and periodicity the
    /// benchmarks pin.
    #[must_use]
    pub fn height_at(&self, x: f64, z: f64, t: f64) -> f64 {
        let mut h = self.mean_level;
        for w in &self.waves {
            h += w.amplitude * w.theta(x, z, t).sin();
        }
        h
    }

    /// The analytic **unit surface normal** at the rest point `(x, z)` and time
    /// `t`.
    ///
    /// Derived from the tangents of the displaced Gerstner surface
    /// (`∂P/∂x × ∂P/∂z`, normalised). For a flat sea this is `+y`. The result is
    /// guaranteed unit-length and points upward (`y > 0`).
    #[must_use]
    pub fn normal_at(&self, x: f64, z: f64, t: f64) -> Vector3<f64> {
        // Tangent vectors of P(x,z) = (x + Σ QA dx cosθ, Σ A sinθ, z + Σ QA dz cosθ).
        // ∂θ/∂x = k dx ; ∂θ/∂z = k dz.
        // Accumulate the Jacobian terms, then cross.
        let mut dpdx = Vector3::new(1.0, 0.0, 0.0);
        let mut dpdz = Vector3::new(0.0, 0.0, 1.0);
        for w in &self.waves {
            let theta = w.theta(x, z, t);
            let (s, c) = theta.sin_cos();
            let qa = w.steepness * w.amplitude;
            let kdx = w.k * w.direction[0];
            let kdz = w.k * w.direction[1];
            // d/dx of the components:
            dpdx.x += -qa * w.direction[0] * s * kdx;
            dpdx.y += w.amplitude * c * kdx;
            dpdx.z += -qa * w.direction[1] * s * kdx;
            // d/dz of the components:
            dpdz.x += -qa * w.direction[0] * s * kdz;
            dpdz.y += w.amplitude * c * kdz;
            dpdz.z += -qa * w.direction[1] * s * kdz;
        }
        let mut n = dpdx.cross(&dpdz);
        // Cross of the two horizontal tangents points up for sane steepness;
        // guard the (degenerate) zero case and the sign.
        let len = n.norm();
        if len <= 0.0 {
            return Vector3::new(0.0, 1.0, 0.0);
        }
        n /= len;
        if n.y < 0.0 {
            n = -n;
        }
        n
    }

    /// A short, deterministic preset: `n` waves spread over a range of
    /// wavelengths and headings from a fixed in-crate seed (no `rand`).
    ///
    /// Wavelengths geometrically span `[base_wavelength, base_wavelength * 2^(n-1)]`
    /// roughly; amplitudes scale `∝ wavelength` at the given `steepness`;
    /// directions fan out around `+x`; phases come from a SplitMix64-style
    /// integer hash of the index. This is a convenience for demos/tests, not a
    /// calibrated sea spectrum (a Phillips/JONSWAP spectral field is a documented
    /// follow-up).
    ///
    /// # Errors
    ///
    /// Propagates [`GerstnerWave::new`] errors (so `n == 0`, a non-positive
    /// `base_wavelength`/`amplitude`/`gravity`, or `steepness ∉ [0, 1]` fail
    /// loud).
    pub fn deterministic_sea(
        n: usize,
        base_wavelength: f64,
        base_amplitude: f64,
        steepness: f64,
        gravity: f64,
        mean_level: f64,
    ) -> Result<Self, OceanError> {
        if n == 0 {
            return Err(OceanError::InvalidConfig(
                "deterministic_sea needs at least one wave".to_string(),
            ));
        }
        let mut waves = Vec::with_capacity(n);
        for i in 0..n {
            // Deterministic, seed-free spread.
            let frac = i as f64 / n as f64;
            let wavelength = base_wavelength * (1.0 + frac); // 1×..~2×
            let amplitude = base_amplitude * (wavelength / base_wavelength);
            // Fan the heading ±0.6 rad around +x.
            let angle = (frac - 0.5) * 1.2;
            let direction = (angle.cos(), angle.sin());
            // Integer-hash phase in [0, 2π) — SplitMix64 finaliser on the index.
            let phase = splitmix_phase(i as u64);
            waves.push(GerstnerWave::new(
                amplitude, wavelength, steepness, phase, direction, gravity,
            )?);
        }
        Self::new(waves, mean_level)
    }
}

/// Map an integer index to a deterministic phase in `[0, 2π)` using the
/// SplitMix64 finaliser (no `rand`, no per-crate state).
fn splitmix_phase(i: u64) -> f64 {
    let mut z = i.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    // Top 53 bits → [0,1) → scale to [0, 2π).
    let u = (z >> 11) as f64 / (1u64 << 53) as f64;
    u * std::f64::consts::TAU
}

#[cfg(test)]
mod tests {
    use super::*;

    const G: f64 = STANDARD_GRAVITY;

    fn single(amplitude: f64, wavelength: f64) -> OceanWaveField {
        let w = GerstnerWave::new(amplitude, wavelength, 0.0, 0.0, (1.0, 0.0), G).unwrap();
        OceanWaveField::new(vec![w], 0.0).unwrap()
    }

    // ---- BENCHMARK-PIN (1): single Gerstner wave ----

    #[test]
    fn crest_to_trough_is_twice_amplitude() {
        let a = 1.7;
        let l = 30.0;
        let field = single(a, l);
        // Sample the height over a full wavelength; max - min must equal 2A.
        let (mut hi, mut lo) = (f64::NEG_INFINITY, f64::INFINITY);
        let n = 100_000;
        for i in 0..n {
            let x = l * i as f64 / n as f64;
            let h = field.height_at(x, 0.0, 0.0);
            hi = hi.max(h);
            lo = lo.min(h);
        }
        // Dense sampling reaches the analytic ±A extrema to well under 1e-6.
        assert!((hi - a).abs() < 1e-6, "crest {hi} != {a}");
        assert!((lo + a).abs() < 1e-6, "trough {lo} != {}", -a);
        assert!(
            ((hi - lo) - 2.0 * a).abs() < 1e-6,
            "height {} != 2A",
            hi - lo
        );
    }

    #[test]
    fn deep_water_phase_speed_matches_sqrt_g_over_k() {
        let l = 42.0;
        let w = GerstnerWave::new(1.0, l, 0.0, 0.0, (1.0, 0.0), G).unwrap();
        let k = 2.0 * std::f64::consts::PI / l;
        let c_expected = (G / k).sqrt();
        assert!((w.phase_speed() - c_expected).abs() < 1e-9);
        // Equivalent: omega/k and sqrt(g*k)/k agree.
        assert!((w.angular_frequency() - (G * k).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn periodic_in_space_to_1e_9() {
        let l = 25.0;
        let field = single(0.9, l);
        for &t in &[0.0, 1.3, 7.7] {
            for &x in &[0.0, 3.0, 11.5, 19.2] {
                let a = field.height_at(x, 0.0, t);
                let b = field.height_at(x + l, 0.0, t);
                assert!((a - b).abs() < 1e-9, "space period: {a} vs {b}");
            }
        }
    }

    #[test]
    fn periodic_in_time_to_1e_9() {
        let l = 18.0;
        let w = GerstnerWave::new(0.6, l, 0.0, 0.0, (1.0, 0.0), G).unwrap();
        let period = w.period();
        let field = OceanWaveField::new(vec![w], 0.0).unwrap();
        for &x in &[0.0, 4.0, 9.5] {
            for &t in &[0.0, 0.7, 2.2] {
                let a = field.height_at(x, 0.0, t);
                let b = field.height_at(x, 0.0, t + period);
                assert!((a - b).abs() < 1e-9, "time period: {a} vs {b}");
            }
        }
    }

    // ---- normal / displacement sanity ----

    #[test]
    fn flat_sea_normal_is_up_and_height_is_mean() {
        let field = OceanWaveField::flat(3.0).unwrap();
        assert_eq!(field.height_at(123.0, -45.0, 9.9), 3.0);
        let n = field.normal_at(1.0, 2.0, 3.0);
        assert!((n - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
        assert_eq!(field.displacement_at(1.0, 2.0, 3.0), Vector3::zeros());
    }

    #[test]
    fn normal_is_unit_and_upward() {
        let field = OceanWaveField::deterministic_sea(5, 20.0, 0.6, 0.5, G, 0.0).unwrap();
        for &(x, z, t) in &[(0.0, 0.0, 0.0), (3.3, -7.1, 2.0), (50.0, 12.0, 9.0)] {
            let n = field.normal_at(x, z, t);
            assert!(
                (n.norm() - 1.0).abs() < 1e-9,
                "normal not unit: {}",
                n.norm()
            );
            assert!(n.y > 0.0, "normal points down: {n:?}");
        }
    }

    #[test]
    fn normal_matches_finite_difference() {
        // Compare the analytic Gerstner normal to a central finite-difference of
        // the displaced surface tangents.
        let field = single(0.4, 30.0); // gentle, well-behaved
        let (x, z, t) = (5.0, 0.0, 1.0);
        let h = 1e-4;
        let p = |xx: f64, zz: f64| {
            let d = field.displacement_at(xx, zz, t);
            Vector3::new(xx + d.x, field.mean_level() + d.y, zz + d.z)
        };
        let tx = (p(x + h, z) - p(x - h, z)) / (2.0 * h);
        let tz = (p(x, z + h) - p(x, z - h)) / (2.0 * h);
        let mut fd = tx.cross(&tz).normalize();
        // `normal_at` documents an upward-oriented (`y > 0`) normal; orient the
        // finite-difference normal the same way so we compare direction, not the
        // arbitrary handedness sign of this particular tangent cross product.
        if fd.y < 0.0 {
            fd = -fd;
        }
        let an = field.normal_at(x, z, t);
        assert!((fd - an).norm() < 1e-5, "fd {fd:?} vs analytic {an:?}");
    }

    // ---- fail-loud config ----

    #[test]
    fn rejects_bad_config() {
        assert!(GerstnerWave::new(0.0, 10.0, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(-1.0, 10.0, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 0.0, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, -5.0, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 10.0, 1.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 10.0, -0.1, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 10.0, 0.5, 0.0, (0.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 10.0, 0.5, 0.0, (1.0, 0.0), 0.0).is_err());
        assert!(GerstnerWave::new(f64::NAN, 10.0, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, f64::INFINITY, 0.5, 0.0, (1.0, 0.0), G).is_err());
        assert!(GerstnerWave::new(1.0, 10.0, 0.5, 0.0, (f64::NAN, 0.0), G).is_err());
        assert!(OceanWaveField::flat(f64::NAN).is_err());
        assert!(OceanWaveField::deterministic_sea(0, 10.0, 0.5, 0.5, G, 0.0).is_err());
    }

    #[test]
    fn deterministic_sea_is_reproducible() {
        let a = OceanWaveField::deterministic_sea(6, 15.0, 0.4, 0.5, G, 0.0).unwrap();
        let b = OceanWaveField::deterministic_sea(6, 15.0, 0.4, 0.5, G, 0.0).unwrap();
        for i in 0..50 {
            let x = i as f64;
            assert_eq!(a.height_at(x, 0.0, 0.0), b.height_at(x, 0.0, 0.0));
        }
    }
}

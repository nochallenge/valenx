//! Aerodynamic drag model: a Mach-dependent drag-coefficient curve plus
//! the dynamic-pressure and drag-force relations that turn an airspeed
//! and atmosphere state into a force.
//!
//! The drag coefficient of a slender launch vehicle is strongly
//! Mach-dependent: a low, roughly constant **subsonic plateau**, a sharp
//! **transonic drag rise** peaking near Mach 1 (the "sound barrier"), and
//! a gradual **supersonic decline** as the bow shock lays back. This
//! module captures that with a documented piecewise-linear
//! `Cd(Mach)` table ([`DragCurve`]) and linear interpolation between
//! nodes — exact at the table nodes, monotone within each modelled
//! region.
//!
//! The dynamic pressure `q = ½ρv²` and drag force `D = q·C_d·A` are exact
//! definitions. Mach number is obtained from the local speed of sound via
//! the shared [`crate::atmosphere`] model, and **max-Q** — the peak
//! dynamic pressure along a trajectory, the structural design driver — is
//! computed over a supplied `(ρ, v)` series.
//!
//! The closed-form pieces (`q`, `D`, the node values) are pinned directly
//! by the unit tests against the textbook formulas.

use serde::{Deserialize, Serialize};

use crate::atmosphere;
use crate::error::{AstroError, Result};

/// A Mach-dependent drag-coefficient curve for a slender launch vehicle.
///
/// The default curve ([`DragCurve::launch_vehicle`]) is the canonical
/// shape: a ~0.20 subsonic plateau, a transonic peak of ~0.50 at Mach 1,
/// and a supersonic decline back toward ~0.20. Construct your own with
/// [`DragCurve::new`] from a strictly-increasing Mach table.
///
/// `Serialize`/`Deserialize` are derived for mission persistence. Note that
/// deserialization restores the stored nodes directly and does **not**
/// re-run [`DragCurve::new`]'s table validation; build untrusted curves via
/// [`DragCurve::new`] rather than deserializing them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DragCurve {
    /// `(Mach, Cd)` nodes, sorted by strictly-increasing Mach. Below the
    /// first node and above the last the endpoint `Cd` is held flat.
    nodes: Vec<(f64, f64)>,
}

impl DragCurve {
    /// The canonical slender-launch-vehicle `Cd(Mach)` curve:
    /// subsonic plateau → transonic peak at Mach 1 → supersonic decline.
    pub fn launch_vehicle() -> Self {
        // Strictly increasing in Mach; physically representative shape.
        Self {
            nodes: vec![
                (0.0, 0.20),
                (0.6, 0.20),
                (0.8, 0.28),
                (1.0, 0.50), // transonic peak
                (1.2, 0.48),
                (2.0, 0.34),
                (3.0, 0.26),
                (5.0, 0.22),
                (8.0, 0.20),
            ],
        }
    }

    /// Build a curve from an explicit `(Mach, Cd)` table.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidAero`] if the table has fewer than two
    /// nodes, contains a non-finite or negative value, or its Mach column
    /// is not strictly increasing (a non-monotone table would make the
    /// interpolation pick an arbitrary bracketing segment).
    pub fn new(nodes: Vec<(f64, f64)>) -> Result<Self> {
        if nodes.len() < 2 {
            return Err(AstroError::InvalidAero("drag curve needs >= 2 nodes"));
        }
        for &(m, cd) in &nodes {
            if !m.is_finite() || m < 0.0 {
                return Err(AstroError::InvalidAero("Mach node must be finite and >= 0"));
            }
            if !cd.is_finite() || cd < 0.0 {
                return Err(AstroError::InvalidAero("Cd node must be finite and >= 0"));
            }
        }
        for w in nodes.windows(2) {
            if w[1].0 <= w[0].0 {
                return Err(AstroError::InvalidAero(
                    "Mach column must be strictly increasing",
                ));
            }
        }
        Ok(Self { nodes })
    }

    /// Drag coefficient at the given Mach number by linear interpolation.
    ///
    /// Exact at the table nodes; below the first / above the last node the
    /// endpoint value is held flat. A non-finite Mach clamps to the
    /// subsonic endpoint so the result is always finite (drag is then
    /// resolved by the dynamic pressure, which a non-finite airspeed would
    /// itself reject upstream).
    pub fn cd(&self, mach: f64) -> f64 {
        if !mach.is_finite() {
            return self.nodes[0].1;
        }
        if mach <= self.nodes[0].0 {
            return self.nodes[0].1;
        }
        let last = self.nodes[self.nodes.len() - 1];
        if mach >= last.0 {
            return last.1;
        }
        for w in self.nodes.windows(2) {
            let (m0, c0) = w[0];
            let (m1, c1) = w[1];
            if mach >= m0 && mach <= m1 {
                let t = (mach - m0) / (m1 - m0);
                return c0 + t * (c1 - c0);
            }
        }
        last.1
    }
}

impl Default for DragCurve {
    fn default() -> Self {
        Self::launch_vehicle()
    }
}

/// Dynamic pressure `q = ½·ρ·v²` (Pa).
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `density` or `speed` is
/// non-finite, or `density` is negative — so a bad input cannot produce a
/// silent `NaN` `q` that then poisons a max-Q reduction.
pub fn dynamic_pressure(density: f64, speed: f64) -> Result<f64> {
    if !density.is_finite() || density < 0.0 {
        return Err(AstroError::InvalidParameter(
            "density must be finite and >= 0",
        ));
    }
    if !speed.is_finite() {
        return Err(AstroError::InvalidParameter("speed must be finite"));
    }
    Ok(0.5 * density * speed * speed)
}

/// Aerodynamic drag force `D = q·C_d·A` (N) from a dynamic pressure (Pa),
/// drag coefficient and reference area (m²).
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if any argument is non-finite
/// or negative.
pub fn drag_force(dynamic_pressure: f64, cd: f64, area: f64) -> Result<f64> {
    for value in [dynamic_pressure, cd, area] {
        if !value.is_finite() || value < 0.0 {
            return Err(AstroError::InvalidParameter(
                "drag_force inputs must be finite and >= 0",
            ));
        }
    }
    Ok(dynamic_pressure * cd * area)
}

/// Mach number `M = v / a(altitude)` from an airspeed (m/s) and a
/// geometric altitude (m), using the shared [`crate::atmosphere`] local
/// speed of sound.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `speed` or `altitude_m` is
/// non-finite. (The standard atmosphere always returns a strictly
/// positive speed of sound, so the division is safe.)
pub fn mach(speed: f64, altitude_m: f64) -> Result<f64> {
    if !speed.is_finite() {
        return Err(AstroError::InvalidParameter("speed must be finite"));
    }
    if !altitude_m.is_finite() {
        return Err(AstroError::InvalidParameter("altitude must be finite"));
    }
    let a = atmosphere::sample(altitude_m).speed_of_sound;
    Ok(speed / a)
}

/// Peak dynamic pressure ("max-Q", Pa) over a `(density, speed)` series,
/// and the index of the sample at which it occurs.
///
/// # Errors
///
/// - [`AstroError::InvalidAero`] (`"empty series"`) for an empty series.
/// - [`AstroError::InvalidParameter`] if any sample is non-physical (see
///   [`dynamic_pressure`]).
pub fn max_q(series: &[(f64, f64)]) -> Result<(f64, usize)> {
    if series.is_empty() {
        return Err(AstroError::InvalidAero("empty (density, speed) series"));
    }
    let mut best_q = f64::NEG_INFINITY;
    let mut best_i = 0usize;
    for (i, &(rho, v)) in series.iter().enumerate() {
        let q = dynamic_pressure(rho, v)?;
        if q > best_q {
            best_q = q;
            best_i = i;
        }
    }
    Ok((best_q, best_i))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dynamic_pressure_is_exact() {
        // q = ½·ρ·v². ρ=1.0, v=100 -> 5000 Pa.
        let q = dynamic_pressure(1.0, 100.0).expect("ok");
        assert!((q - 5_000.0).abs() < 1e-9, "q = {q}");
        // ρ=1.225 (sea level), v=250 -> 0.5·1.225·62500 = 38281.25.
        let q2 = dynamic_pressure(1.225, 250.0).expect("ok");
        assert!((q2 - 38_281.25).abs() < 1e-6, "q2 = {q2}");
    }

    #[test]
    fn drag_force_is_exact() {
        // D = q·Cd·A. q=5000, Cd=0.5, A=10 -> 25000 N.
        let d = drag_force(5_000.0, 0.5, 10.0).expect("ok");
        assert!((d - 25_000.0).abs() < 1e-9, "D = {d}");
    }

    #[test]
    fn cd_curve_exact_at_nodes() {
        let c = DragCurve::launch_vehicle();
        assert!((c.cd(0.0) - 0.20).abs() < 1e-12);
        assert!((c.cd(1.0) - 0.50).abs() < 1e-12); // transonic peak node
        assert!((c.cd(2.0) - 0.34).abs() < 1e-12);
        assert!((c.cd(8.0) - 0.20).abs() < 1e-12);
    }

    #[test]
    fn cd_curve_linear_interpolation_pinned() {
        let c = DragCurve::launch_vehicle();
        // Midpoint of (0.8,0.28)->(1.0,0.50): t=0.5 -> 0.39.
        assert!((c.cd(0.9) - 0.39).abs() < 1e-12, "cd(0.9) = {}", c.cd(0.9));
        // (1.2,0.48)->(2.0,0.34): M=1.5 -> t=0.375 -> 0.4275.
        assert!(
            (c.cd(1.5) - 0.4275).abs() < 1e-12,
            "cd(1.5) = {}",
            c.cd(1.5)
        );
    }

    #[test]
    fn cd_curve_shape_subsonic_transonic_supersonic() {
        let c = DragCurve::launch_vehicle();
        let sub = c.cd(0.3); // subsonic plateau
        let peak = c.cd(1.0); // transonic peak
        let suprsnc = c.cd(4.0); // supersonic decline
                                 // Transonic peak is the maximum; subsonic and high-supersonic are
                                 // both well below it.
        assert!(peak > sub, "peak {peak} should exceed subsonic {sub}");
        assert!(
            peak > suprsnc,
            "peak {peak} should exceed supersonic {suprsnc}"
        );
        // Subsonic is a low plateau ~0.2.
        assert!((sub - 0.20).abs() < 1e-9);
        // Supersonic declines monotonically past the peak across the
        // table's supersonic nodes.
        let mut prev = c.cd(1.2);
        for &m in &[1.5, 2.0, 3.0, 5.0, 8.0] {
            let v = c.cd(m);
            assert!(
                v <= prev + 1e-12,
                "Cd rose in supersonic at M={m}: {v} > {prev}"
            );
            prev = v;
        }
    }

    #[test]
    fn cd_held_flat_outside_table() {
        let c = DragCurve::launch_vehicle();
        assert!((c.cd(-1.0) - 0.20).abs() < 1e-12); // below first node
        assert!((c.cd(20.0) - 0.20).abs() < 1e-12); // above last node
    }

    #[test]
    fn mach_uses_atmosphere_speed_of_sound() {
        // Sea-level speed of sound ≈ 340.3 m/s; at that speed M ≈ 1.
        let a = atmosphere::sample(0.0).speed_of_sound;
        let m = mach(a, 0.0).expect("ok");
        assert!((m - 1.0).abs() < 1e-9, "M = {m}");
        // Half the speed of sound -> M = 0.5.
        let m2 = mach(0.5 * a, 0.0).expect("ok");
        assert!((m2 - 0.5).abs() < 1e-9);
    }

    #[test]
    fn max_q_finds_the_peak() {
        // A synthetic ascent: density falls, speed rises; the product
        // peaks in the middle. Hand-built so the max is unambiguous.
        let series = [
            (1.225, 0.0),  // q = 0
            (0.9, 150.0),  // q = 10125
            (0.6, 300.0),  // q = 27000  <- peak
            (0.3, 400.0),  // q = 24000
            (0.05, 600.0), // q = 9000
        ];
        let (q, i) = max_q(&series).expect("ok");
        assert_eq!(i, 2, "peak index");
        assert!((q - 27_000.0).abs() < 1e-6, "max-Q = {q}");
    }

    #[test]
    fn new_rejects_bad_tables() {
        assert!(DragCurve::new(vec![(0.0, 0.2)]).is_err()); // < 2 nodes
                                                            // Non-monotone Mach column.
        assert!(DragCurve::new(vec![(0.0, 0.2), (0.5, 0.3), (0.5, 0.4)]).is_err());
        assert!(DragCurve::new(vec![(0.0, 0.2), (-0.1, 0.3)]).is_err());
        assert!(DragCurve::new(vec![(0.0, 0.2), (1.0, -0.3)]).is_err()); // negative Cd
                                                                         // A valid custom table works.
        assert!(DragCurve::new(vec![(0.0, 0.25), (1.0, 0.6), (3.0, 0.3)]).is_ok());
    }

    #[test]
    fn rejects_non_physical_inputs() {
        assert!(dynamic_pressure(f64::NAN, 100.0).is_err());
        assert!(dynamic_pressure(-1.0, 100.0).is_err());
        assert!(dynamic_pressure(1.0, f64::INFINITY).is_err());
        assert!(drag_force(-1.0, 0.5, 10.0).is_err());
        assert!(mach(f64::NAN, 0.0).is_err());
        assert!(max_q(&[]).is_err());
    }
}

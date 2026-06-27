//! Two-parameter cubic equations of state: Peng–Robinson and Soave–Redlich–Kwong.
//!
//! Both can be written in the unified pressure-explicit form
//!
//! ```text
//! P = R T / (V - b)  -  a(T) / ((V + δ1 b)(V + δ2 b))
//! ```
//!
//! where `(δ1, δ2) = (1 + √2, 1 - √2)` for Peng–Robinson and `(1, 0)` for SRK.
//! In dimensionless form this becomes a cubic in the compressibility factor
//! `Z = P V / (R T)`:
//!
//! ```text
//! Z³ + c2 Z² + c1 Z + c0 = 0
//! ```
//!
//! with coefficients built from the reduced parameters `A = a P / (R T)²` and
//! `B = b P / (R T)`. This module solves that cubic in closed form (Cardano),
//! returns the liquid and vapor compressibility roots, and evaluates the
//! analytic fugacity coefficient used for phase-equilibrium calculations.
//!
//! References: Peng & Robinson, *Ind. Eng. Chem. Fundam.* **15** (1976) 59;
//! Soave, *Chem. Eng. Sci.* **27** (1972) 1197; Poling, Prausnitz & O'Connell,
//! *The Properties of Gases and Liquids*, 5th ed., chs. 4–6.

use crate::error::{require_positive, Result, ThermoError, GAS_CONSTANT};
use crate::fluid::Fluid;

/// Which cubic equation of state to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EosModel {
    /// Peng–Robinson (1976). Good for vapor pressures and liquid densities of
    /// non-polar and slightly polar fluids; the default for most engineering work.
    PengRobinson,
    /// Soave–Redlich–Kwong (1972). Good for vapor-phase fugacities of light gases.
    Srk,
}

impl EosModel {
    /// The constant `Ωa` in `a_c = Ωa R² Tc² / Pc`.
    fn omega_a(self) -> f64 {
        match self {
            EosModel::PengRobinson => 0.457_235_528_921_382_3,
            EosModel::Srk => 0.427_480_233_540_245_6,
        }
    }

    /// The constant `Ωb` in `b = Ωb R Tc / Pc`.
    fn omega_b(self) -> f64 {
        match self {
            EosModel::PengRobinson => 0.077_796_073_903_888_4,
            EosModel::Srk => 0.086_640_349_964_980_3,
        }
    }

    /// The cubic-EoS volume-translation constants `(δ1, δ2)`.
    fn deltas(self) -> (f64, f64) {
        match self {
            EosModel::PengRobinson => (
                1.0 + std::f64::consts::SQRT_2,
                1.0 - std::f64::consts::SQRT_2,
            ),
            EosModel::Srk => (1.0, 0.0),
        }
    }

    /// The `κ` coefficient of the Soave/PR `α(T)` correlation,
    /// `α = [1 + κ(1 - √Tr)]²`.
    fn kappa(self, omega: f64) -> f64 {
        match self {
            EosModel::PengRobinson => 0.374_64 + 1.542_26 * omega - 0.269_92 * omega * omega,
            EosModel::Srk => 0.480 + 1.574 * omega - 0.176 * omega * omega,
        }
    }
}

/// A cubic equation of state bound to a particular fluid.
///
/// Construct with [`Eos::new`], then evaluate pressure, compressibility roots,
/// and fugacity coefficients at a chosen temperature/pressure.
#[derive(Debug, Clone, Copy)]
pub struct Eos {
    /// The fluid whose critical constants parametrize this EoS.
    pub fluid: Fluid,
    /// Which cubic model is in use.
    pub model: EosModel,
}

/// The two physically meaningful compressibility roots of the cubic at a state.
///
/// When the cubic has three real roots the smallest corresponds to the liquid
/// and the largest to the vapor (the middle root is non-physical). When it has
/// a single real root (e.g. a supercritical state) both fields hold that root.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZRoots {
    /// Liquid-like (smallest) compressibility factor.
    pub liquid: f64,
    /// Vapor-like (largest) compressibility factor.
    pub vapor: f64,
    /// `true` if the cubic had three distinct real roots (two-phase region
    /// possible); `false` if only one real root exists.
    pub three_real: bool,
}

impl Eos {
    /// Bind a cubic equation of state to a fluid.
    #[must_use]
    pub fn new(fluid: Fluid, model: EosModel) -> Self {
        Eos { fluid, model }
    }

    /// The temperature-dependent attraction parameter `a(T)` in SI units
    /// (Pa·m⁶/mol²).
    #[must_use]
    pub fn a(&self, t: f64) -> f64 {
        let f = &self.fluid;
        let ac = self.model.omega_a() * GAS_CONSTANT * GAS_CONSTANT * f.tc * f.tc / f.pc;
        let tr = t / f.tc;
        let kappa = self.model.kappa(f.omega);
        let sqrt_alpha = 1.0 + kappa * (1.0 - tr.sqrt());
        ac * sqrt_alpha * sqrt_alpha
    }

    /// The (temperature-independent) co-volume parameter `b` in SI units
    /// (m³/mol).
    #[must_use]
    pub fn b(&self) -> f64 {
        let f = &self.fluid;
        self.model.omega_b() * GAS_CONSTANT * f.tc / f.pc
    }

    /// Pressure `P` (Pa) predicted by the EoS at temperature `t` (K) and molar
    /// volume `v` (m³/mol).
    ///
    /// # Errors
    ///
    /// Returns [`ThermoError::NonPositive`] if `t` or `v` is not strictly
    /// positive, or [`ThermoError::OutOfRange`] if `v <= b` (inside the
    /// hard-core volume, where the EoS is undefined).
    pub fn pressure(&self, t: f64, v: f64) -> Result<f64> {
        require_positive("temperature", t)?;
        require_positive("molar volume", v)?;
        let b = self.b();
        if v <= b {
            return Err(ThermoError::OutOfRange {
                name: "molar volume",
                value: v,
                expected: "must exceed the co-volume b",
            });
        }
        let (d1, d2) = self.model.deltas();
        let a = self.a(t);
        Ok(GAS_CONSTANT * t / (v - b) - a / ((v + d1 * b) * (v + d2 * b)))
    }

    /// Reduced cubic parameters `(A, B)` at temperature `t` (K) and pressure
    /// `p` (Pa).
    fn reduced_ab(&self, t: f64, p: f64) -> (f64, f64) {
        let rt = GAS_CONSTANT * t;
        let a = self.a(t) * p / (rt * rt);
        let b = self.b() * p / rt;
        (a, b)
    }

    /// Solve the cubic in `Z` at the given temperature and pressure, returning
    /// the liquid and vapor compressibility roots.
    ///
    /// # Errors
    ///
    /// Returns [`ThermoError::NonPositive`] if `t` or `p` is not strictly
    /// positive.
    pub fn z_roots(&self, t: f64, p: f64) -> Result<ZRoots> {
        require_positive("temperature", t)?;
        require_positive("pressure", p)?;
        let (a, b) = self.reduced_ab(t, p);
        let (d1, d2) = self.model.deltas();

        // General cubic Z³ + c2 Z² + c1 Z + c0 = 0 for
        // P = RT/(V-b) - a/((V+δ1 b)(V+δ2 b)). With u = δ1+δ2, w = δ1·δ2:
        //   c2 = (u - 1) B - 1
        //   c1 = A + (w - u) B² - u B
        //   c0 = -(A B + w B² + w B³)
        let u = d1 + d2;
        let w = d1 * d2;
        let c2 = (u - 1.0) * b - 1.0;
        let c1 = a + (w - u) * b * b - u * b;
        let c0 = -(a * b + w * b * b + w * b * b * b);

        let roots = real_cubic_roots(c2, c1, c0);
        // Keep only physically admissible roots (Z must exceed B, i.e. V > b).
        let mut phys: Vec<f64> = roots.iter().copied().filter(|&z| z > b).collect();
        if phys.is_empty() {
            // Fall back to the full set if the B-filter removed everything
            // (can happen at extreme states from round-off); pick positives.
            phys = roots.iter().copied().filter(|&z| z > 0.0).collect();
        }
        if phys.is_empty() {
            return Err(ThermoError::NoSuchRoot {
                what: "compressibility factor Z",
            });
        }
        phys.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let liquid = *phys.first().unwrap();
        let vapor = *phys.last().unwrap();
        Ok(ZRoots {
            liquid,
            vapor,
            three_real: phys.len() >= 3,
        })
    }

    /// Natural log of the fugacity coefficient `ln φ` for the phase with
    /// compressibility `z` at `(t, p)`.
    ///
    /// The closed-form expression for a generic two-parameter cubic is
    ///
    /// ```text
    /// ln φ = (Z - 1) - ln(Z - B) - A/(B(δ1-δ2)) · ln((Z+δ1 B)/(Z+δ2 B))
    /// ```
    ///
    /// For SRK, where `δ2 = 0`, the `1/(δ1-δ2)` prefactor and the log reduce to
    /// the familiar `-(A/B) ln((Z+B)/Z)` form.
    ///
    /// # Errors
    ///
    /// Returns [`ThermoError::OutOfRange`] if `z` does not exceed the reduced
    /// co-volume `B` (a non-physical root).
    pub fn ln_phi(&self, t: f64, p: f64, z: f64) -> Result<f64> {
        let (a, b) = self.reduced_ab(t, p);
        let (d1, d2) = self.model.deltas();
        if z <= b {
            return Err(ThermoError::OutOfRange {
                name: "compressibility factor",
                value: z,
                expected: "must exceed reduced co-volume B",
            });
        }
        let term = a / (b * (d1 - d2)) * ((z + d1 * b) / (z + d2 * b)).ln();
        Ok((z - 1.0) - (z - b).ln() - term)
    }

    /// Fugacity coefficient `φ = f/P` of the requested phase at `(t, p)`.
    ///
    /// `vapor = true` selects the vapor (largest-Z) root; `false` selects the
    /// liquid (smallest-Z) root.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`Eos::z_roots`] and returns
    /// [`ThermoError::OutOfRange`] for a non-physical root.
    pub fn fugacity_coefficient(&self, t: f64, p: f64, vapor: bool) -> Result<f64> {
        let roots = self.z_roots(t, p)?;
        let z = if vapor { roots.vapor } else { roots.liquid };
        Ok(self.ln_phi(t, p, z)?.exp())
    }
}

/// Return the real roots of the depressed-or-general monic cubic
/// `x³ + a2 x² + a1 x + a0 = 0` using Cardano's formula.
///
/// Returns a vector of 1 or 3 real roots (a double root is reported with
/// multiplicity, yielding 3 entries).
fn real_cubic_roots(a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    // Depress: x = y - a2/3  ⇒  y³ + p y + q = 0.
    let shift = a2 / 3.0;
    let p = a1 - a2 * a2 / 3.0;
    let q = 2.0 * a2 * a2 * a2 / 27.0 - a2 * a1 / 3.0 + a0;

    let disc = (q * q) / 4.0 + (p * p * p) / 27.0;
    let eps = 1e-14;

    if disc > eps {
        // One real root.
        let sqrt_disc = disc.sqrt();
        let u = cbrt(-q / 2.0 + sqrt_disc);
        let v = cbrt(-q / 2.0 - sqrt_disc);
        vec![u + v - shift]
    } else if disc < -eps {
        // Three distinct real roots — trigonometric form.
        let m = 2.0 * (-p / 3.0).sqrt();
        let theta = (3.0 * q / (p * m)).clamp(-1.0, 1.0).acos();
        let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;
        vec![
            m * (theta / 3.0).cos() - shift,
            m * (theta / 3.0 - two_pi_3).cos() - shift,
            m * (theta / 3.0 + two_pi_3).cos() - shift,
        ]
    } else {
        // Discriminant ≈ 0: a multiple real root.
        if p.abs() < eps && q.abs() < eps {
            let r = -shift;
            vec![r, r, r]
        } else {
            let single = 3.0 * q / p;
            let double = -3.0 * q / (2.0 * p);
            vec![single - shift, double - shift, double - shift]
        }
    }
}

/// Real cube root that handles negative arguments (unlike `f64::powf(1/3)`).
fn cbrt(x: f64) -> f64 {
    x.cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_known_roots() {
        // (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6.
        let mut r = real_cubic_roots(-6.0, 11.0, -6.0);
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(r.len(), 3);
        for (got, want) in r.iter().zip([1.0, 2.0, 3.0]) {
            assert!((got - want).abs() < 1e-9, "got {got}, want {want}");
        }
    }

    #[test]
    fn cubic_single_real_root() {
        // x³ + x + 1 has one real root near -0.6823.
        let r = real_cubic_roots(0.0, 1.0, 1.0);
        assert_eq!(r.len(), 1);
        assert!((r[0] - (-0.682_327_8)).abs() < 1e-6, "got {}", r[0]);
    }

    #[test]
    fn ideal_gas_limit_z_to_one() {
        // As P → 0 at fixed T, Z(vapor) → 1 for any real fluid.
        let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
        let t = 350.0;
        let z_hi = eos.z_roots(t, 1.0e5).unwrap().vapor;
        let z_lo = eos.z_roots(t, 1.0e2).unwrap().vapor;
        let z_vlo = eos.z_roots(t, 1.0e0).unwrap().vapor;
        assert!(z_vlo > z_lo && z_lo > z_hi, "Z should rise toward 1 as P→0");
        assert!((z_vlo - 1.0).abs() < 1e-4, "Z(P→0) = {z_vlo}, want ≈ 1");
    }

    #[test]
    fn fugacity_coefficient_to_one_at_low_pressure() {
        // φ → 1 as P → 0 (ideal-gas limit).
        let eos = Eos::new(Fluid::nitrogen(), EosModel::Srk);
        let phi = eos.fugacity_coefficient(300.0, 1.0, true).unwrap();
        assert!((phi - 1.0).abs() < 1e-4, "φ(P→0) = {phi}");
    }

    #[test]
    fn pressure_rejects_bad_input() {
        let eos = Eos::new(Fluid::methane(), EosModel::PengRobinson);
        assert!(matches!(
            eos.pressure(-300.0, 1e-3),
            Err(ThermoError::NonPositive { .. })
        ));
        assert!(matches!(
            eos.z_roots(300.0, -1e5),
            Err(ThermoError::NonPositive { .. })
        ));
    }
}

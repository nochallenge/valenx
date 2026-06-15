//! Free (unforced) time response `x(t)` of an SDOF system.
//!
//! Given initial displacement `x0` and initial velocity `v0`, the
//! solution of `m x'' + c x' + k x = 0` is one of three closed forms
//! selected by the damping regime (Rao, *Mechanical Vibrations*, ch. 2):
//!
//! - **Underdamped** (`0 <= zeta < 1`):
//!   ```text
//!   x(t) = e^(-zeta*wn*t) * [ x0*cos(wd*t)
//!                             + (v0 + zeta*wn*x0)/wd * sin(wd*t) ]
//!   ```
//!   a decaying oscillation at the damped frequency `wd`, with the
//!   undamped case `zeta = 0` recovering pure `x0*cos(wn*t) +
//!   (v0/wn)*sin(wn*t)`.
//!
//! - **Critically damped** (`zeta = 1`, double root `s = -wn`):
//!   ```text
//!   x(t) = [ x0 + (v0 + wn*x0) * t ] * e^(-wn*t)
//!   ```
//!
//! - **Overdamped** (`zeta > 1`, two real roots `s1, s2`):
//!   ```text
//!   x(t) = A * e^(s1*t) + B * e^(s2*t)
//!   ```
//!   with `s1,2 = wn*(-zeta +/- sqrt(zeta^2 - 1))` and `A, B` fixed by
//!   the initial conditions.
//!
//! Each closed form is constructed once into a [`FreeResponse`] and then
//! evaluated cheaply at any time `t`.

use crate::error::VibrationError;
use crate::model::{DampingRegime, SdofSystem};
use serde::{Deserialize, Serialize};

/// A pre-solved free-vibration response, ready to evaluate at any time.
///
/// Build one with [`FreeResponse::new`] from a system and its initial
/// conditions, then call [`displacement`](FreeResponse::displacement) to
/// sample `x(t)`. The variant captures the closed-form coefficients for
/// the system's damping regime so repeated evaluation is just a handful
/// of `exp`/`sin`/`cos` calls.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum FreeResponse {
    /// Decaying (or pure, when `zeta = 0`) oscillation.
    Underdamped {
        /// Decay rate `zeta*wn` (1/s) in the `e^(-decay*t)` envelope.
        decay: f64,
        /// Damped angular frequency `wd` (rad/s).
        damped_freq: f64,
        /// Coefficient of `cos(wd*t)` — equals the initial displacement.
        cos_coeff: f64,
        /// Coefficient of `sin(wd*t)`.
        sin_coeff: f64,
    },
    /// Non-oscillating fastest decay: `(a + b*t) e^(-wn*t)`.
    CriticallyDamped {
        /// Undamped natural frequency `wn` (rad/s); the decay rate.
        natural_freq: f64,
        /// Constant term `a` (= initial displacement).
        a: f64,
        /// Linear-in-`t` term `b`.
        b: f64,
    },
    /// Non-oscillating slow decay: `A e^(s1*t) + B e^(s2*t)`.
    Overdamped {
        /// First (less negative) real root `s1` (1/s).
        s1: f64,
        /// Second (more negative) real root `s2` (1/s).
        s2: f64,
        /// Weight `A` on the `s1` mode.
        a: f64,
        /// Weight `B` on the `s2` mode.
        b: f64,
    },
}

impl FreeResponse {
    /// Solve the free response of `system` for initial displacement
    /// `x0` (m) and initial velocity `v0` (m/s).
    ///
    /// # Errors
    ///
    /// Returns [`VibrationError::BadParameter`] if `x0` or `v0` is not a
    /// finite number.
    pub fn new(system: &SdofSystem, x0: f64, v0: f64) -> Result<Self, VibrationError> {
        if !x0.is_finite() {
            return Err(VibrationError::BadParameter {
                name: "x0",
                reason: format!("initial displacement must be finite, got {x0}"),
            });
        }
        if !v0.is_finite() {
            return Err(VibrationError::BadParameter {
                name: "v0",
                reason: format!("initial velocity must be finite, got {v0}"),
            });
        }

        let wn = system.natural_freq_rad_s();
        let zeta = system.damping_ratio();

        match system.regime() {
            DampingRegime::Undamped | DampingRegime::Underdamped => {
                // wd is real and > 0 here.
                let wd = system.damped_freq_rad_s()?;
                Ok(FreeResponse::Underdamped {
                    decay: zeta * wn,
                    damped_freq: wd,
                    cos_coeff: x0,
                    sin_coeff: (v0 + zeta * wn * x0) / wd,
                })
            }
            DampingRegime::CriticallyDamped => Ok(FreeResponse::CriticallyDamped {
                natural_freq: wn,
                a: x0,
                b: v0 + wn * x0,
            }),
            DampingRegime::Overdamped => {
                let root = wn * (zeta * zeta - 1.0).sqrt();
                let s1 = -zeta * wn + root; // less negative
                let s2 = -zeta * wn - root; // more negative
                                            // Solve  A + B = x0 ; s1*A + s2*B = v0.
                let a = (v0 - s2 * x0) / (s1 - s2);
                let b = x0 - a;
                Ok(FreeResponse::Overdamped { s1, s2, a, b })
            }
        }
    }

    /// Evaluate the displacement `x(t)` (m) at time `t` (s).
    ///
    /// Time may be any finite value; for `t = 0` this returns the
    /// initial displacement `x0` exactly (up to floating-point round-off
    /// in the overdamped weight split).
    pub fn displacement(&self, t: f64) -> f64 {
        match *self {
            FreeResponse::Underdamped {
                decay,
                damped_freq,
                cos_coeff,
                sin_coeff,
            } => {
                let env = (-decay * t).exp();
                env * (cos_coeff * (damped_freq * t).cos() + sin_coeff * (damped_freq * t).sin())
            }
            FreeResponse::CriticallyDamped { natural_freq, a, b } => {
                (a + b * t) * (-natural_freq * t).exp()
            }
            FreeResponse::Overdamped { s1, s2, a, b } => a * (s1 * t).exp() + b * (s2 * t).exp(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Loose tolerance for analytic-vs-evaluated comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn initial_displacement_recovered_at_t0() {
        // Across all three regimes, x(0) == x0.
        for zeta in [0.0_f64, 0.3, 1.0, 2.5] {
            let sys = SdofSystem::from_modal(12.0, zeta).expect("valid");
            let resp = FreeResponse::new(&sys, 0.37, -1.1).expect("valid");
            assert!((resp.displacement(0.0) - 0.37).abs() < EPS, "zeta = {zeta}");
        }
    }

    #[test]
    fn undamped_is_pure_cosine_for_zero_velocity() {
        // zeta = 0, v0 = 0  =>  x(t) = x0 cos(wn t).
        let wn = 5.0;
        let sys = SdofSystem::from_modal(wn, 0.0).expect("valid");
        let x0 = 2.0;
        let resp = FreeResponse::new(&sys, x0, 0.0).expect("valid");

        // Quarter period: cos(pi/2) = 0.
        let t_quarter = std::f64::consts::FRAC_PI_2 / wn;
        assert!(resp.displacement(t_quarter).abs() < EPS);

        // Half period: cos(pi) = -1  =>  x = -x0.
        let t_half = std::f64::consts::PI / wn;
        assert!((resp.displacement(t_half) - (-x0)).abs() < EPS);

        // Full period: cos(2pi) = 1  =>  x = x0.
        let t_full = std::f64::consts::TAU / wn;
        assert!((resp.displacement(t_full) - x0).abs() < EPS);
    }

    #[test]
    fn undamped_amplitude_is_conserved() {
        // With v0 != 0 the amplitude is sqrt(x0^2 + (v0/wn)^2); the
        // response should never exceed it and should reach it.
        let wn = 8.0;
        let sys = SdofSystem::from_modal(wn, 0.0).expect("valid");
        let (x0, v0) = (1.0, 4.0);
        let amp = (x0 * x0 + (v0 / wn).powi(2)).sqrt();
        let resp = FreeResponse::new(&sys, x0, v0).expect("valid");

        let mut peak = 0.0_f64;
        for i in 0..2000 {
            let t = i as f64 * 0.001;
            peak = peak.max(resp.displacement(t).abs());
            assert!(resp.displacement(t).abs() <= amp + 1e-6);
        }
        // Over more than one full period the peak essentially reaches amp.
        assert!((peak - amp).abs() < 1e-3);
    }

    #[test]
    fn underdamped_oscillates_and_decays() {
        // A lightly-damped system should change sign (oscillate) and
        // have a strictly smaller envelope at successive same-phase
        // times (decay).
        let sys = SdofSystem::from_modal(20.0, 0.05).expect("valid");
        let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");
        let td = sys.damped_period_s().expect("ud");

        // Sign change within the first half-period => it oscillates.
        let half = resp.displacement(td / 2.0);
        assert!(half < 0.0, "expected sign change, got {half}");

        // Peaks one damped period apart decay by the envelope ratio
        // e^(-zeta*wn*Td) < 1.
        let p0 = resp.displacement(0.0);
        let p1 = resp.displacement(td);
        assert!(p1 < p0 && p1 > 0.0);
        assert!(p1 < p0); // strict decay
    }

    #[test]
    fn underdamped_envelope_matches_exponential() {
        // At t = n*Td the response equals x0 * e^(-zeta*wn*n*Td).
        let wn = 15.0;
        let zeta = 0.1;
        let sys = SdofSystem::from_modal(wn, zeta).expect("valid");
        let x0 = 1.0;
        let resp = FreeResponse::new(&sys, x0, 0.0).expect("valid");
        let td = sys.damped_period_s().expect("ud");

        for n in 1..=4 {
            let t = n as f64 * td;
            let expected = x0 * (-zeta * wn * t).exp();
            // At whole damped periods sin term vanishes, cos = 1.
            assert!((resp.displacement(t) - expected).abs() < 1e-6, "n = {n}");
        }
    }

    #[test]
    fn critically_damped_does_not_oscillate() {
        // Released from rest, a critically-damped system decays
        // monotonically toward zero without crossing it.
        let sys = SdofSystem::from_modal(10.0, 1.0).expect("valid");
        let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");

        let mut prev = resp.displacement(0.0);
        for i in 1..=500 {
            let t = i as f64 * 0.005;
            let x = resp.displacement(t);
            assert!(x >= -EPS, "overshot below zero: {x}");
            assert!(x <= prev + EPS, "not monotone decreasing");
            prev = x;
        }
        // And it has essentially decayed after several time constants.
        assert!(resp.displacement(3.0) < 1e-3);
    }

    #[test]
    fn overdamped_decays_without_oscillating() {
        let sys = SdofSystem::from_modal(6.0, 2.0).expect("valid");
        let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");

        let mut prev = resp.displacement(0.0);
        for i in 1..=1000 {
            let t = i as f64 * 0.005;
            let x = resp.displacement(t);
            assert!(x >= -EPS, "overdamped should not undershoot: {x}");
            assert!(x <= prev + EPS, "overdamped should be monotone");
            prev = x;
        }
    }

    #[test]
    fn overdamped_matches_two_exponentials_known_roots() {
        // Pick wn and zeta giving clean integer roots.
        // s = wn(-zeta +/- sqrt(zeta^2 - 1)). With wn = 5, zeta = 1.25:
        // sqrt(zeta^2 - 1) = sqrt(0.5625) = 0.75
        // s1 = 5*(-1.25 + 0.75) = -2.5 ; s2 = 5*(-1.25 - 0.75) = -10.
        let sys = SdofSystem::from_modal(5.0, 1.25).expect("valid");
        let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");

        // Released from rest: A + B = 1, s1 A + s2 B = 0.
        // B = A*(-s1/s2) ... solve: A = -s2/(s1 - s2) = 10/7.5 = 1.3333..
        let s1 = -2.5;
        let s2 = -10.0;
        let a = -s2 / (s1 - s2);
        let b = 1.0 - a;
        for i in 0..50 {
            let t = i as f64 * 0.02;
            let expected = a * (s1 * t).exp() + b * (s2 * t).exp();
            assert!((resp.displacement(t) - expected).abs() < 1e-9, "i = {i}");
        }
    }

    #[test]
    fn rejects_non_finite_initial_conditions() {
        let sys = SdofSystem::from_modal(10.0, 0.1).expect("valid");
        assert!(FreeResponse::new(&sys, f64::NAN, 0.0).is_err());
        assert!(FreeResponse::new(&sys, 0.0, f64::INFINITY).is_err());
    }
}

//! The Ziegler-Nichols ultimate-gain tuning table.
//!
//! Given a validated [`UltimateMeasurement`] `(Ku, Tu)`, the classic
//! 1942 Ziegler-Nichols closed-loop rules produce controller settings
//! for three structures:
//!
//! | Controller | `Kp`       | `Ti`       | `Td`       |
//! | ---------- | ---------- | ---------- | ---------- |
//! | P          | `0.5 Ku`   | infinite   | `0`        |
//! | PI         | `0.45 Ku`  | `Tu / 1.2` | `0`        |
//! | PID        | `0.6 Ku`   | `Tu / 2`   | `Tu / 8`   |
//!
//! Here `Kp` is the proportional gain, `Ti` the integral (reset) time
//! and `Td` the derivative time of the *standard* (a.k.a. ideal /
//! interacting-free) controller
//!
//! ```text
//! u(t) = Kp * ( e(t) + (1/Ti) integral e dt + Td de/dt ).
//! ```
//!
//! A P controller has no integral action, so its reset time is infinite
//! (an infinite `Ti` makes the `1/Ti` integral term vanish); P and PI
//! controllers have no derivative action, so their `Td` is zero. Both of
//! those degenerate cases are represented faithfully:
//! [`Gains::integral_time`] is `f64::INFINITY` for the P controller, and
//! its parallel-form integral gain [`Gains::ki`] is `0`.
//!
//! ## Honest scope
//!
//! These are the textbook closed-form constants and nothing more. The
//! rules target roughly quarter-amplitude decay and are known to give
//! aggressive, lightly damped responses; they presuppose a clean
//! ultimate-gain experiment and offer no robustness or noise guarantees.
//! Treat the output as an educational starting point, not a deployable
//! tune.

use serde::{Deserialize, Serialize};

use crate::ultimate::UltimateMeasurement;

/// Which controller structure a tuning rule targets.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum ControllerKind {
    /// Proportional only.
    P,
    /// Proportional + integral.
    Pi,
    /// Proportional + integral + derivative (classic Z-N).
    Pid,
}

impl ControllerKind {
    /// Whether this structure includes integral (reset) action.
    pub fn has_integral(self) -> bool {
        matches!(self, ControllerKind::Pi | ControllerKind::Pid)
    }

    /// Whether this structure includes derivative action.
    pub fn has_derivative(self) -> bool {
        matches!(self, ControllerKind::Pid)
    }
}

/// A complete set of controller settings in standard time-constant form.
///
/// The triple `(Kp, Ti, Td)` describes the standard controller
///
/// ```text
/// u(t) = Kp * ( e(t) + (1/Ti) integral e dt + Td de/dt ).
/// ```
///
/// For controllers without integral action `Ti` is [`f64::INFINITY`];
/// for controllers without derivative action `Td` is `0`. The
/// equivalent *parallel* (independent-gain) form
///
/// ```text
/// u(t) = Kp e + Ki integral e dt + Kd de/dt
/// ```
///
/// is available through [`Gains::ki`] and [`Gains::kd`], with
/// `Ki = Kp / Ti` and `Kd = Kp * Td`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Gains {
    /// Which controller structure these settings target.
    pub kind: ControllerKind,
    /// Proportional gain `Kp` (dimensionless).
    pub kp: f64,
    /// Integral (reset) time `Ti`, in seconds; [`f64::INFINITY`] when
    /// the controller has no integral action.
    pub ti: f64,
    /// Derivative time `Td`, in seconds; `0` when the controller has no
    /// derivative action.
    pub td: f64,
}

impl Gains {
    /// The proportional gain `Kp`.
    pub fn kp(&self) -> f64 {
        self.kp
    }

    /// The integral (reset) time `Ti`, in seconds.
    ///
    /// Returns [`f64::INFINITY`] for a pure-P controller, faithfully
    /// encoding "no integral action".
    pub fn integral_time(&self) -> f64 {
        self.ti
    }

    /// The derivative time `Td`, in seconds.
    ///
    /// Returns `0` for P and PI controllers.
    pub fn derivative_time(&self) -> f64 {
        self.td
    }

    /// The parallel-form integral gain `Ki = Kp / Ti`.
    ///
    /// Because a pure-P controller has `Ti = INFINITY`, this evaluates
    /// to exactly `0` there, matching the absence of integral action.
    pub fn ki(&self) -> f64 {
        self.kp / self.ti
    }

    /// The parallel-form derivative gain `Kd = Kp * Td`.
    ///
    /// Evaluates to `0` for P and PI controllers.
    pub fn kd(&self) -> f64 {
        self.kp * self.td
    }
}

/// Ziegler-Nichols ultimate-gain rules over a validated measurement.
///
/// Construct from an [`UltimateMeasurement`] and read off the settings
/// for whichever controller structure you need. Each accessor applies
/// the corresponding row of the closed-loop table; the constants are
/// the classic 1942 values and are not configurable.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ZieglerNichols {
    measurement: UltimateMeasurement,
}

impl ZieglerNichols {
    /// Wrap a validated measurement in the Ziegler-Nichols tuner.
    pub fn new(measurement: UltimateMeasurement) -> Self {
        Self { measurement }
    }

    /// The measurement these rules are derived from.
    pub fn measurement(&self) -> UltimateMeasurement {
        self.measurement
    }

    /// Proportional-only settings: `Kp = 0.5 Ku`.
    ///
    /// `Ti` is [`f64::INFINITY`] and `Td` is `0`.
    pub fn p(&self) -> Gains {
        let ku = self.measurement.ultimate_gain();
        Gains {
            kind: ControllerKind::P,
            kp: 0.5 * ku,
            ti: f64::INFINITY,
            td: 0.0,
        }
    }

    /// PI settings: `Kp = 0.45 Ku`, `Ti = Tu / 1.2`.
    ///
    /// `Td` is `0`.
    pub fn pi(&self) -> Gains {
        let ku = self.measurement.ultimate_gain();
        let tu = self.measurement.ultimate_period();
        Gains {
            kind: ControllerKind::Pi,
            kp: 0.45 * ku,
            ti: tu / 1.2,
            td: 0.0,
        }
    }

    /// Classic PID settings: `Kp = 0.6 Ku`, `Ti = Tu / 2`, `Td = Tu / 8`.
    pub fn pid(&self) -> Gains {
        let ku = self.measurement.ultimate_gain();
        let tu = self.measurement.ultimate_period();
        Gains {
            kind: ControllerKind::Pid,
            kp: 0.6 * ku,
            ti: tu / 2.0,
            td: tu / 8.0,
        }
    }

    /// Settings for an arbitrary [`ControllerKind`].
    ///
    /// Dispatches to [`ZieglerNichols::p`], [`ZieglerNichols::pi`] or
    /// [`ZieglerNichols::pid`] according to `kind`.
    pub fn gains(&self, kind: ControllerKind) -> Gains {
        match kind {
            ControllerKind::P => self.p(),
            ControllerKind::Pi => self.pi(),
            ControllerKind::Pid => self.pid(),
        }
    }
}

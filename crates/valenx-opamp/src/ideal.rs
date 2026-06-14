//! Ideal closed-loop op-amp configurations.
//!
//! All formulae here assume an **ideal op-amp**: infinite open-loop
//! gain, infinite input impedance, zero output impedance, and the two
//! resulting "golden rules" — no current flows into the inputs, and
//! negative feedback drives the differential input voltage to zero
//! (`V+ = V-`). Under those assumptions each topology reduces to an
//! algebraic resistor ratio.

use crate::error::{ensure_positive, Result};
use serde::{Deserialize, Serialize};

/// An inverting amplifier: input through `r_in` into the virtual-ground
/// summing node, feedback `r_f`, non-inverting input grounded.
///
/// Closed-loop voltage gain is `G = -r_f / r_in`, always negative — the
/// output is inverted relative to the input.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Inverting {
    /// Input resistance `Rin` (ohms). Must be `> 0`.
    pub r_in: f64,
    /// Feedback resistance `Rf` (ohms). Must be `> 0`.
    pub r_f: f64,
}

impl Inverting {
    /// Construct an inverting stage, validating both resistors are
    /// finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if either resistance is
    /// non-finite or `<= 0`.
    pub fn new(r_in: f64, r_f: f64) -> Result<Self> {
        let r_in = ensure_positive("r_in", r_in)?;
        let r_f = ensure_positive("r_f", r_f)?;
        Ok(Self { r_in, r_f })
    }

    /// Closed-loop voltage gain `G = -Rf / Rin` (dimensionless).
    ///
    /// Always negative for the ideal inverting topology.
    pub fn gain(&self) -> f64 {
        -self.r_f / self.r_in
    }

    /// Magnitude of the closed-loop gain, `|Rf / Rin|`.
    ///
    /// This is the value that enters the gain-bandwidth relations in
    /// [`crate::bandwidth`].
    pub fn gain_magnitude(&self) -> f64 {
        self.gain().abs()
    }

    /// Output voltage for a given input voltage, `Vout = G · Vin`.
    pub fn output(&self, v_in: f64) -> f64 {
        self.gain() * v_in
    }
}

/// A non-inverting amplifier: input drives `V+` directly, the feedback
/// divider `r_f` (top) / `r_in` (bottom, to ground) sets `V-`.
///
/// Closed-loop voltage gain is `G = 1 + r_f / r_in`, always `>= 1` — the
/// output is in phase with, and never smaller than, the input.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NonInverting {
    /// Lower divider resistance `Rin` to ground (ohms). Must be `> 0`.
    pub r_in: f64,
    /// Feedback resistance `Rf` from output to the `V-` node (ohms).
    /// Must be `> 0`.
    pub r_f: f64,
}

impl NonInverting {
    /// Construct a non-inverting stage, validating both resistors are
    /// finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if either resistance is
    /// non-finite or `<= 0`.
    pub fn new(r_in: f64, r_f: f64) -> Result<Self> {
        let r_in = ensure_positive("r_in", r_in)?;
        let r_f = ensure_positive("r_f", r_f)?;
        Ok(Self { r_in, r_f })
    }

    /// Closed-loop voltage gain `G = 1 + Rf / Rin` (dimensionless).
    ///
    /// Strictly greater than `1` for any positive resistor pair; it
    /// approaches `1` only in the limit `Rf -> 0`.
    pub fn gain(&self) -> f64 {
        1.0 + self.r_f / self.r_in
    }

    /// Magnitude of the closed-loop gain.
    ///
    /// Identical to [`gain`](Self::gain) here since the non-inverting
    /// gain is always positive; provided for symmetry with
    /// [`Inverting::gain_magnitude`].
    pub fn gain_magnitude(&self) -> f64 {
        self.gain().abs()
    }

    /// Output voltage for a given input voltage, `Vout = G · Vin`.
    pub fn output(&self, v_in: f64) -> f64 {
        self.gain() * v_in
    }
}

/// A unity-gain buffer (voltage follower): the special case of a
/// non-inverting amplifier with `Rf = 0`, giving exactly `G = 1`.
///
/// Modelled as its own zero-parameter type because the follower needs
/// no resistors and its gain is fixed.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct VoltageFollower;

impl VoltageFollower {
    /// Construct a voltage follower. Infallible — there is nothing to
    /// validate.
    pub fn new() -> Self {
        Self
    }

    /// Closed-loop gain — exactly `1.0` by definition.
    pub fn gain(&self) -> f64 {
        1.0
    }

    /// Output voltage equals the input voltage.
    pub fn output(&self, v_in: f64) -> f64 {
        v_in
    }
}

/// One input branch of a [`SummingAmplifier`]: a source voltage `v`
/// driving the virtual-ground node through resistance `r`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SummingInput {
    /// Branch input voltage `Vᵢ` (volts).
    pub v: f64,
    /// Branch input resistance `Rᵢ` (ohms). Must be `> 0`.
    pub r: f64,
}

/// An inverting summing amplifier: several inputs share one
/// virtual-ground summing node and a common feedback resistor `r_f`.
///
/// The output is the negated, resistor-weighted sum of the inputs:
/// `Vout = -Rf · Σ(Vᵢ / Rᵢ)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SummingAmplifier {
    /// Common feedback resistance `Rf` (ohms). Must be `> 0`.
    pub r_f: f64,
    /// The input branches; at least one is required.
    pub inputs: Vec<SummingInput>,
}

impl SummingAmplifier {
    /// Construct a summing amplifier from a feedback resistance and a
    /// non-empty set of `(voltage, resistance)` input branches.
    ///
    /// Every resistance (the feedback resistor and each branch
    /// resistor) is validated finite and `> 0`; each branch voltage is
    /// validated finite.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError::NoInputs`](crate::OpAmpError::NoInputs) for
    /// an empty input set, or a numeric variant if any value fails the
    /// finite / positive checks.
    pub fn new(r_f: f64, inputs: impl IntoIterator<Item = (f64, f64)>) -> Result<Self> {
        let r_f = ensure_positive("r_f", r_f)?;
        let inputs: Vec<SummingInput> = inputs
            .into_iter()
            .enumerate()
            .map(|(i, (v, r))| {
                let _ = i; // index reserved for future per-branch error context
                let v = crate::error::ensure_finite("v", v)?;
                let r = ensure_positive("r", r)?;
                Ok(SummingInput { v, r })
            })
            .collect::<Result<_>>()?;
        if inputs.is_empty() {
            return Err(crate::OpAmpError::NoInputs);
        }
        Ok(Self { r_f, inputs })
    }

    /// Output voltage `Vout = -Rf · Σ(Vᵢ / Rᵢ)`.
    pub fn output(&self) -> f64 {
        let weighted: f64 = self.inputs.iter().map(|i| i.v / i.r).sum();
        -self.r_f * weighted
    }

    /// Per-branch closed-loop gain `-Rf / Rᵢ` for branch `index`.
    ///
    /// Returns `None` if `index` is out of range.
    pub fn branch_gain(&self, index: usize) -> Option<f64> {
        self.inputs.get(index).map(|i| -self.r_f / i.r)
    }
}

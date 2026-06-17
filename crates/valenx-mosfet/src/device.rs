//! N-channel enhancement MOSFET square-law (Shockley level-1) model.
//!
//! This module implements the idealized long-channel drain-current
//! equations for an n-channel enhancement-mode MOSFET, expressed with a
//! single lumped transconductance parameter `k` (units: A/V², equal to
//! `μ_n · C_ox · W / L` in the device-physics decomposition). All
//! voltages are referenced to the source terminal and given in volts.
//!
//! # Conventions
//!
//! - `vgs` — gate-to-source voltage (V).
//! - `vds` — drain-to-source voltage (V), assumed `≥ 0` for an NMOS.
//! - `vth` — threshold voltage (V), positive for an enhancement NMOS.
//! - `vov = vgs − vth` — gate overdrive (a.k.a. effective voltage).
//!
//! # Regions
//!
//! With `vov = vgs − vth`:
//!
//! - **Cutoff** (`vov ≤ 0`): the channel is off, `Id = 0`.
//! - **Triode / linear** (`vov > 0` and `vds < vov`): the channel
//!   conducts as a voltage-controlled resistor,
//!   `Id = k · (vov · vds − ½ · vds²)`.
//! - **Saturation** (`vov > 0` and `vds ≥ vov`): the channel pinches
//!   off, the current is (ideally) independent of `vds`,
//!   `Id = ½ · k · vov²`.
//!
//! The triode and saturation expressions agree at the boundary
//! `vds = vov`: substituting `vds = vov` into the triode formula gives
//! `k · (vov² − ½ vov²) = ½ k vov²`, the saturation value. This
//! continuity is asserted in the tests.
//!
//! # Transconductance
//!
//! The small-signal transconductance is `gm = ∂Id/∂Vgs`. In
//! saturation, differentiating `Id = ½ k (vgs − vth)²` with respect to
//! `vgs` yields `gm = k · (vgs − vth) = k · vov`. This crate reports the
//! saturation-region `gm` from [`Mosfet::gm`]; it is zero in cutoff
//! because `Id` is identically zero there.

use serde::{Deserialize, Serialize};

use crate::error::{MosfetError, Result};

/// Operating region of the MOSFET for a given bias point.
///
/// Returned by [`Mosfet::region`] and reported alongside the current by
/// [`Mosfet::operating_point`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Region {
    /// Channel off (`vgs − vth ≤ 0`): `Id = 0`.
    Cutoff,
    /// Linear / ohmic region (`vds < vgs − vth`): channel acts as a
    /// voltage-controlled resistor.
    Triode,
    /// Saturation / active region (`vds ≥ vgs − vth`): channel pinched
    /// off, current ideally flat in `vds`.
    Saturation,
}

impl Region {
    /// Short human-readable label for UI / logs.
    pub fn label(self) -> &'static str {
        match self {
            Region::Cutoff => "cutoff",
            Region::Triode => "triode",
            Region::Saturation => "saturation",
        }
    }
}

/// A bias point: the resolved [`Region`], drain current `id` (A), and
/// small-signal transconductance `gm` (S = A/V).
///
/// Produced by [`Mosfet::operating_point`] so a single call yields all
/// three quantities consistently for one `(vgs, vds)` pair.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Resolved operating region.
    pub region: Region,
    /// Drain current `Id` in amperes.
    pub id: f64,
    /// Transconductance `gm = dId/dVgs` in siemens (saturation value).
    pub gm: f64,
}

/// N-channel enhancement-mode MOSFET described by the square-law model.
///
/// Construct with [`Mosfet::new`] (validates parameters) or
/// [`Mosfet::nmos`] for a convenient default. The struct holds only the
/// two parameters the level-1 model needs: the lumped transconductance
/// parameter `k` and the threshold voltage `vth`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Mosfet {
    /// Transconductance parameter `k = μ_n · C_ox · W / L` (A/V²).
    /// Strictly positive.
    k: f64,
    /// Threshold voltage `vth` (V). Positive for an enhancement NMOS;
    /// any finite value is accepted so depletion-style thresholds can be
    /// modeled too.
    vth: f64,
}

impl Mosfet {
    /// Build a MOSFET from its transconductance parameter `k` (A/V²) and
    /// threshold voltage `vth` (V).
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Invalid`] if `k` is not strictly positive,
    /// and [`MosfetError::Domain`] if either argument is non-finite
    /// (`NaN` / `±∞`).
    pub fn new(k: f64, vth: f64) -> Result<Self> {
        if !k.is_finite() {
            return Err(MosfetError::domain("k", "must be finite"));
        }
        if !vth.is_finite() {
            return Err(MosfetError::domain("vth", "must be finite"));
        }
        if k <= 0.0 {
            return Err(MosfetError::invalid(
                "k",
                format!("transconductance parameter must be > 0, got {k}"),
            ));
        }
        Ok(Self { k, vth })
    }

    /// Convenience constructor for a textbook NMOS: `k = 0.5 mA/V²`
    /// (`0.5e-3` A/V²) and `vth = 1.0 V`.
    ///
    /// Equivalent to `Mosfet::new(0.5e-3, 1.0).unwrap()`; provided so the
    /// common case needs no error handling at the call site.
    pub fn nmos() -> Self {
        // Safe: both literals are finite and k > 0.
        Self {
            k: 0.5e-3,
            vth: 1.0,
        }
    }

    /// Transconductance parameter `k` (A/V²).
    pub fn k(&self) -> f64 {
        self.k
    }

    /// Threshold voltage `vth` (V).
    pub fn vth(&self) -> f64 {
        self.vth
    }

    /// Gate overdrive `vov = vgs − vth` (V) for the supplied gate bias.
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `vgs` is non-finite.
    pub fn overdrive(&self, vgs: f64) -> Result<f64> {
        check_finite("vgs", vgs)?;
        Ok(vgs - self.vth)
    }

    /// Classify the operating [`Region`] for a bias point.
    ///
    /// Uses the standard inequalities with `vov = vgs − vth`:
    /// cutoff when `vov ≤ 0`; triode when `vds < vov`; saturation when
    /// `vds ≥ vov`. The saturation/triode boundary `vds == vov` is
    /// assigned to saturation (the conventional choice; both formulas
    /// give the same current there).
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `vgs` or `vds` is non-finite.
    pub fn region(&self, vgs: f64, vds: f64) -> Result<Region> {
        check_finite("vgs", vgs)?;
        check_finite("vds", vds)?;
        let vov = vgs - self.vth;
        if vov <= 0.0 {
            return Ok(Region::Cutoff);
        }
        if vds < vov {
            Ok(Region::Triode)
        } else {
            Ok(Region::Saturation)
        }
    }

    /// Drain current `Id` (A) for a bias point under the square-law
    /// model.
    ///
    /// Selects the region with [`Mosfet::region`] and applies:
    ///
    /// - cutoff: `Id = 0`;
    /// - triode: `Id = k · (vov · vds − ½ vds²)`;
    /// - saturation: `Id = ½ · k · vov²`.
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `vgs` or `vds` is non-finite.
    pub fn drain_current(&self, vgs: f64, vds: f64) -> Result<f64> {
        let region = self.region(vgs, vds)?;
        let vov = vgs - self.vth;
        let id = match region {
            Region::Cutoff => 0.0,
            Region::Triode => self.k * (vov * vds - 0.5 * vds * vds),
            Region::Saturation => 0.5 * self.k * vov * vov,
        };
        Ok(id)
    }

    /// Small-signal transconductance `gm = dId/dVgs` (S) at the gate
    /// bias `vgs`, evaluated for the saturation region.
    ///
    /// Returns `k · (vgs − vth) = k · vov` when `vov > 0`, and `0` when
    /// `vov ≤ 0` (the device is in cutoff, where `Id ≡ 0`).
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `vgs` is non-finite.
    pub fn gm(&self, vgs: f64) -> Result<f64> {
        let vov = self.overdrive(vgs)?;
        if vov <= 0.0 {
            Ok(0.0)
        } else {
            Ok(self.k * vov)
        }
    }

    /// Resolve the full [`OperatingPoint`] — region, drain current, and
    /// saturation transconductance — for one bias point in a single
    /// call.
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `vgs` or `vds` is non-finite.
    pub fn operating_point(&self, vgs: f64, vds: f64) -> Result<OperatingPoint> {
        Ok(OperatingPoint {
            region: self.region(vgs, vds)?,
            id: self.drain_current(vgs, vds)?,
            gm: self.gm(vgs)?,
        })
    }

    /// Gate **overdrive** `vov = vgs − vth` (V) needed to carry a target
    /// saturation drain current — the design inverse of the saturation
    /// branch of [`drain_current`](Mosfet::drain_current).
    ///
    /// Inverting `Id = ½ · k · vov²` gives `vov = sqrt(2 · Id / k)`. A
    /// zero target maps to zero overdrive (the cutoff edge `vgs = vth`).
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `target_id` is non-finite, and
    /// [`MosfetError::Invalid`] if it is negative (an NMOS drain current
    /// is non-negative).
    pub fn overdrive_for_saturation_current(&self, target_id: f64) -> Result<f64> {
        if !target_id.is_finite() {
            return Err(MosfetError::domain("target_id", "must be finite"));
        }
        if target_id < 0.0 {
            return Err(MosfetError::invalid(
                "target_id",
                format!("drain current must be >= 0, got {target_id}"),
            ));
        }
        Ok((2.0 * target_id / self.k).sqrt())
    }

    /// Gate-to-source voltage `vgs` (V) that biases the device to a target
    /// saturation drain current — the analog-design inverse that picks the
    /// bias for a desired current.
    ///
    /// `vgs = vth + sqrt(2 · Id / k)`, i.e. the threshold plus the
    /// [overdrive](Mosfet::overdrive_for_saturation_current) the current
    /// demands. Holding `vds ≥ vgs − vth` keeps the device in saturation,
    /// where feeding this `vgs` back into
    /// [`drain_current`](Mosfet::drain_current) reproduces `target_id`.
    ///
    /// # Errors
    ///
    /// Returns [`MosfetError::Domain`] if `target_id` is non-finite, and
    /// [`MosfetError::Invalid`] if it is negative.
    pub fn vgs_for_saturation_current(&self, target_id: f64) -> Result<f64> {
        Ok(self.vth + self.overdrive_for_saturation_current(target_id)?)
    }
}

/// Reject non-finite bias voltages with a [`MosfetError::Domain`].
fn check_finite(what: &'static str, v: f64) -> Result<()> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(MosfetError::domain(what, "must be finite"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for current comparisons (A). The
    /// example currents are O(1e-3); 1e-12 leaves ~9 orders of margin.
    const EPS: f64 = 1e-12;

    fn assert_close(a: f64, b: f64, ctx: &str) {
        assert!((a - b).abs() < EPS, "{ctx}: |{a} - {b}| >= {EPS}");
    }

    // ---- construction / validation ----------------------------------

    #[test]
    fn new_rejects_nonpositive_k() {
        assert!(matches!(
            Mosfet::new(0.0, 1.0),
            Err(MosfetError::Invalid { what: "k", .. })
        ));
        assert!(matches!(
            Mosfet::new(-1.0e-3, 1.0),
            Err(MosfetError::Invalid { what: "k", .. })
        ));
    }

    #[test]
    fn new_rejects_nonfinite_params() {
        assert!(matches!(
            Mosfet::new(f64::NAN, 1.0),
            Err(MosfetError::Domain { what: "k", .. })
        ));
        assert!(matches!(
            Mosfet::new(1.0e-3, f64::INFINITY),
            Err(MosfetError::Domain { what: "vth", .. })
        ));
    }

    #[test]
    fn new_accepts_valid_params_and_exposes_them() {
        let m = Mosfet::new(2.0e-3, 0.7).expect("valid");
        assert_close(m.k(), 2.0e-3, "k accessor");
        assert_close(m.vth(), 0.7, "vth accessor");
    }

    #[test]
    fn nmos_default_matches_new() {
        let a = Mosfet::nmos();
        let b = Mosfet::new(0.5e-3, 1.0).expect("valid");
        assert_eq!(a, b);
    }

    // ---- VALIDATE: Id = 0 for Vgs < Vth (cutoff) --------------------

    #[test]
    fn cutoff_current_is_zero_below_threshold() {
        let m = Mosfet::nmos(); // vth = 1.0
                                // vgs strictly below vth across a sweep of vds.
        for &vds in &[0.0, 0.5, 1.0, 3.0, 10.0] {
            let id = m.drain_current(0.5, vds).expect("finite");
            assert_eq!(m.region(0.5, vds).expect("finite"), Region::Cutoff);
            assert_close(id, 0.0, "cutoff current");
        }
    }

    #[test]
    fn cutoff_at_exact_threshold() {
        // vgs == vth => vov == 0 => cutoff (boundary uses vov <= 0).
        let m = Mosfet::nmos();
        assert_eq!(m.region(1.0, 2.0).expect("finite"), Region::Cutoff);
        assert_close(m.drain_current(1.0, 2.0).expect("finite"), 0.0, "vov=0");
    }

    // ---- VALIDATE: saturation is quadratic in overdrive -------------

    #[test]
    fn saturation_current_is_half_k_vov_squared() {
        // k = 2 A/V^2 keeps the arithmetic exact and easy to read.
        let m = Mosfet::new(2.0, 1.0).expect("valid");
        // vgs = 4 => vov = 3; vds = 10 >> vov so saturation.
        let id = m.drain_current(4.0, 10.0).expect("finite");
        assert_eq!(m.region(4.0, 10.0).expect("finite"), Region::Saturation);
        // 0.5 * 2 * 3^2 = 9.
        assert_close(id, 9.0, "saturation Id");
    }

    #[test]
    fn saturation_scales_as_overdrive_squared() {
        // Doubling the overdrive should quadruple the saturation current.
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let id1 = m.drain_current(2.0, 5.0).expect("finite"); // vov = 1
        let id2 = m.drain_current(3.0, 5.0).expect("finite"); // vov = 2
        assert_eq!(m.region(2.0, 5.0).expect("finite"), Region::Saturation);
        assert_eq!(m.region(3.0, 5.0).expect("finite"), Region::Saturation);
        // id1 = 0.5e-3, id2 = 2.0e-3 => ratio 4.
        assert_close(id1, 0.5e-3, "vov=1 sat");
        assert_close(id2, 2.0e-3, "vov=2 sat");
        assert_close(id2 / id1, 4.0, "quadratic scaling");
    }

    #[test]
    fn saturation_is_flat_in_vds() {
        // Ideal level-1 saturation current is independent of vds.
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let a = m.drain_current(3.0, 2.0).expect("finite"); // vds == vov boundary
        let b = m.drain_current(3.0, 5.0).expect("finite");
        let c = m.drain_current(3.0, 50.0).expect("finite");
        assert_close(a, b, "sat flat a==b");
        assert_close(b, c, "sat flat b==c");
    }

    // ---- VALIDATE: triode at low Vds --------------------------------

    #[test]
    fn triode_formula_at_low_vds() {
        let m = Mosfet::new(2.0, 1.0).expect("valid");
        // vgs = 4 => vov = 3; vds = 1 < vov so triode.
        let id = m.drain_current(4.0, 1.0).expect("finite");
        assert_eq!(m.region(4.0, 1.0).expect("finite"), Region::Triode);
        // k*(vov*vds - 0.5*vds^2) = 2*(3*1 - 0.5*1) = 2*2.5 = 5.
        assert_close(id, 5.0, "triode Id");
    }

    #[test]
    fn triode_small_vds_is_nearly_ohmic() {
        // For vds << vov, Id ~= k*vov*vds (the 0.5*vds^2 term is tiny):
        // a linear (resistive) channel. Check the on-resistance slope.
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let vgs = 3.0; // vov = 2
        let vds = 1.0e-3; // tiny
        let id = m.drain_current(vgs, vds).expect("finite");
        // k*vov*vds = 1e-3 * 2 * 1e-3 = 2e-6; the -0.5*k*vds^2 term is
        // 5e-10, negligible at EPS=1e-12? No — 5e-10 > 1e-12, so compare
        // against the exact triode value instead.
        let exact = 1.0e-3 * (2.0 * vds - 0.5 * vds * vds);
        assert_close(id, exact, "triode exact");
        // And confirm it is close to the ohmic approximation.
        let ohmic = 1.0e-3 * 2.0 * vds;
        assert!((id - ohmic).abs() < 1.0e-9, "near-ohmic: {id} vs {ohmic}");
    }

    // ---- VALIDATE: region boundary (vds vs vov) ---------------------

    #[test]
    fn region_boundary_saturation_when_vds_ge_overdrive() {
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let vgs = 3.0; // vov = 2
                       // vds just below vov => triode.
        assert_eq!(
            m.region(vgs, 1.999_999).expect("finite"),
            Region::Triode,
            "just below boundary"
        );
        // vds exactly vov => saturation (conventional boundary choice).
        assert_eq!(
            m.region(vgs, 2.0).expect("finite"),
            Region::Saturation,
            "at boundary"
        );
        // vds above vov => saturation.
        assert_eq!(
            m.region(vgs, 2.000_001).expect("finite"),
            Region::Saturation,
            "above boundary"
        );
    }

    #[test]
    fn triode_and_saturation_agree_at_pinchoff() {
        // Continuity: at vds = vov the triode and saturation closed
        // forms must give the same current.
        let m = Mosfet::new(3.3e-3, 0.8).expect("valid");
        for &vgs in &[1.5, 2.0, 3.0, 5.0] {
            let vov = vgs - m.vth();
            // Evaluate the triode expression at exactly the pinch-off
            // point and the saturation expression, independent of the
            // region selector, then compare.
            let triode = m.k() * (vov * vov - 0.5 * vov * vov);
            let sat = 0.5 * m.k() * vov * vov;
            assert_close(triode, sat, "pinchoff continuity");
            // The public API at the boundary returns the saturation value.
            let id = m.drain_current(vgs, vov).expect("finite");
            assert_close(id, sat, "api at boundary");
        }
    }

    // ---- VALIDATE: gm = dId/dVgs (finite-difference) ----------------

    #[test]
    fn gm_equals_k_times_overdrive_in_saturation() {
        let m = Mosfet::new(2.0e-3, 1.0).expect("valid");
        // vgs = 3 => vov = 2 => gm = k*vov = 4e-3.
        let gm = m.gm(3.0).expect("finite");
        assert_close(gm, 4.0e-3, "gm closed form");
    }

    #[test]
    fn gm_is_zero_in_cutoff() {
        let m = Mosfet::nmos(); // vth = 1
        assert_close(m.gm(0.5).expect("finite"), 0.0, "gm below vth");
        assert_close(m.gm(1.0).expect("finite"), 0.0, "gm at vth");
    }

    #[test]
    fn gm_matches_numeric_derivative_of_saturation_current() {
        // gm = dId/dVgs verified by a central finite difference of the
        // saturation drain current. Keep vds large so both bias points
        // stay in saturation.
        let m = Mosfet::new(1.2e-3, 0.9).expect("valid");
        let vgs = 2.5;
        let vds = 20.0;
        let h = 1.0e-6;
        let id_plus = m.drain_current(vgs + h, vds).expect("finite");
        let id_minus = m.drain_current(vgs - h, vds).expect("finite");
        let numeric = (id_plus - id_minus) / (2.0 * h);
        let analytic = m.gm(vgs).expect("finite");
        // Central difference error on a quadratic is ~0 to round-off;
        // use a slightly looser tolerance to absorb f64 cancellation.
        assert!(
            (numeric - analytic).abs() < 1.0e-9,
            "gm fd mismatch: numeric {numeric}, analytic {analytic}"
        );
    }

    // ---- operating_point aggregates the three quantities ------------

    #[test]
    fn operating_point_is_consistent_with_individual_calls() {
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let vgs = 3.0;
        let vds = 5.0;
        let op = m.operating_point(vgs, vds).expect("finite");
        assert_eq!(op.region, m.region(vgs, vds).expect("finite"));
        assert_close(op.id, m.drain_current(vgs, vds).expect("finite"), "op id");
        assert_close(op.gm, m.gm(vgs).expect("finite"), "op gm");
    }

    // ---- domain rejection at the IV boundary ------------------------

    #[test]
    fn nonfinite_bias_is_rejected() {
        let m = Mosfet::nmos();
        assert!(matches!(
            m.drain_current(f64::NAN, 1.0),
            Err(MosfetError::Domain { what: "vgs", .. })
        ));
        assert!(matches!(
            m.drain_current(1.0, f64::INFINITY),
            Err(MosfetError::Domain { what: "vds", .. })
        ));
        assert!(matches!(
            m.gm(f64::NAN),
            Err(MosfetError::Domain { what: "vgs", .. })
        ));
        assert!(matches!(
            m.region(1.0, f64::NAN),
            Err(MosfetError::Domain { what: "vds", .. })
        ));
    }

    // ---- region labels ----------------------------------------------

    #[test]
    fn region_labels_are_stable() {
        assert_eq!(Region::Cutoff.label(), "cutoff");
        assert_eq!(Region::Triode.label(), "triode");
        assert_eq!(Region::Saturation.label(), "saturation");
    }

    // ---- VALIDATE: saturation bias-design inverse -------------------

    #[test]
    fn vgs_for_saturation_current_inverts_drain_current() {
        let m = Mosfet::new(1.0e-3, 1.0).expect("valid");
        let target = 2.0e-3;
        let vgs = m.vgs_for_saturation_current(target).expect("finite");
        // Bias well into saturation (large vds) and recover the target.
        assert_eq!(m.region(vgs, 50.0).expect("finite"), Region::Saturation);
        assert_close(
            m.drain_current(vgs, 50.0).expect("finite"),
            target,
            "round-trip Id",
        );
    }

    #[test]
    fn overdrive_for_saturation_current_matches_hand_value() {
        // k = 2, Id = 9 => vov = sqrt(2*9/2) = 3; vgs = vth + 3 = 4.
        let m = Mosfet::new(2.0, 1.0).expect("valid");
        let vov = m.overdrive_for_saturation_current(9.0).expect("finite");
        assert_close(vov, 3.0, "vov hand value");
        let vgs = m.vgs_for_saturation_current(9.0).expect("finite");
        assert_close(vgs, 4.0, "vgs hand value");
        // The forward overdrive() recovers the same vov from that vgs.
        assert_close(
            m.overdrive(vgs).expect("finite"),
            vov,
            "overdrive round-trip",
        );
    }

    #[test]
    fn saturation_bias_inverse_is_consistent_with_gm() {
        // gm at the computed bias equals sqrt(2*k*Id), the gm-Id identity.
        let m = Mosfet::new(1.2e-3, 0.7).expect("valid");
        let target = 0.6e-3;
        let vgs = m.vgs_for_saturation_current(target).expect("finite");
        let gm = m.gm(vgs).expect("finite");
        let expected = (2.0 * m.k() * target).sqrt();
        assert_close(gm, expected, "gm = sqrt(2 k Id)");
    }

    #[test]
    fn zero_target_current_biases_at_threshold() {
        let m = Mosfet::nmos(); // vth = 1
        assert_close(
            m.overdrive_for_saturation_current(0.0).expect("finite"),
            0.0,
            "vov at Id=0",
        );
        assert_close(
            m.vgs_for_saturation_current(0.0).expect("finite"),
            m.vth(),
            "vgs at Id=0",
        );
    }

    #[test]
    fn saturation_bias_inverse_rejects_bad_target() {
        let m = Mosfet::nmos();
        assert!(matches!(
            m.vgs_for_saturation_current(-1.0e-3),
            Err(MosfetError::Invalid {
                what: "target_id",
                ..
            })
        ));
        assert!(matches!(
            m.overdrive_for_saturation_current(f64::NAN),
            Err(MosfetError::Domain {
                what: "target_id",
                ..
            })
        ));
    }
}

//! Preload, stiffness sharing and separation of a bolted joint.
//!
//! ## Preload from torque
//!
//! Tightening a nut to a torque `T` stretches the bolt to a preload
//! tension `F`. The two are linked by the **nut-factor** (or
//! "torque-coefficient") relation
//!
//! ```text
//! T = K * F * d
//! ```
//!
//! where `d` is the nominal bolt diameter and `K` a dimensionless
//! friction coefficient — the famous `K ≈ 0.2` for as-received steel.
//! Inverting gives the achieved preload `F = T / (K d)`. A higher `K`
//! (more friction) means more of the wrench torque is lost to friction,
//! so a higher torque is needed to reach the same preload.
//!
//! ## Stiffness sharing under external load
//!
//! Once preloaded, the bolt (stiffness `kb`) and the clamped members
//! (stiffness `km`) act as two springs in parallel. An external tensile
//! load `P` applied to the joint does **not** all go into the bolt: it
//! is shared in proportion to stiffness through the **joint stiffness
//! constant**
//!
//! ```text
//! C = kb / (kb + km)
//! ```
//!
//! The bolt tension rises by `C * P` (so the total bolt force is
//! `F + C P`) while the clamping force between the members drops by
//! `(1 - C) P`. Because members are usually much stiffer than the
//! slender bolt, `C` is typically small (~0.2–0.3): most of an external
//! load is carried by *unloading the joint*, not by stretching the bolt.
//!
//! ## Separation
//!
//! The members stay clamped only while the clamping force is positive.
//! Setting the residual clamp `F - (1 - C) P` to zero gives the
//! **separation load**
//!
//! ```text
//! P_sep = F / (1 - C)
//! ```
//!
//! Beyond `P_sep` the joint gaps open and the bolt alone carries the
//! full external load — the failure mode this analysis exists to avoid.

use serde::{Deserialize, Serialize};

use crate::error::BoltError;

/// A nut factor `K` validated to lie in `(0, 1)`.
///
/// `K` lumps thread friction, under-head friction and the thread-helix
/// geometry into a single torque coefficient. Common values:
/// ~0.10 (well-lubricated), ~0.20 (plain steel, "as received"),
/// ~0.30 (dry, zinc-plated).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NutFactor(f64);

impl NutFactor {
    /// The textbook default for plain steel fasteners, `K = 0.2`.
    pub const STEEL_AS_RECEIVED: f64 = 0.2;

    /// Build a validated nut factor.
    ///
    /// # Errors
    ///
    /// Returns [`BoltError::NutFactorRange`] if `k` is not in the open
    /// interval `(0, 1)`, or [`BoltError::NotFinite`] if it is NaN /
    /// infinite.
    pub fn new(k: f64) -> Result<Self, BoltError> {
        if !k.is_finite() {
            return Err(BoltError::NotFinite {
                name: "K",
                value: k,
            });
        }
        if k <= 0.0 || k >= 1.0 {
            return Err(BoltError::NutFactorRange { value: k });
        }
        Ok(Self(k))
    }

    /// The underlying coefficient.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// The dimensionless joint stiffness constant `C = kb / (kb + km)`,
/// validated to lie in `(0, 1)`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StiffnessRatio(f64);

impl StiffnessRatio {
    /// Build the stiffness constant from the bolt and member
    /// stiffnesses, `C = kb / (kb + km)`.
    ///
    /// Both stiffnesses share a unit (e.g. N/m); the result is
    /// dimensionless. Because both are strictly positive, the result is
    /// guaranteed to land in `(0, 1)`.
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if either stiffness is not strictly
    /// positive and finite.
    pub fn from_stiffnesses(bolt: f64, member: f64) -> Result<Self, BoltError> {
        let kb = BoltError::require_positive("bolt_stiffness", bolt)?;
        let km = BoltError::require_positive("member_stiffness", member)?;
        Ok(Self(kb / (kb + km)))
    }

    /// Wrap an already-computed ratio, validating it lies in `(0, 1)`.
    ///
    /// # Errors
    ///
    /// Returns [`BoltError::StiffnessRatioRange`] if `c` is not in
    /// `(0, 1)`, or [`BoltError::NotFinite`] if it is NaN / infinite.
    pub fn new(c: f64) -> Result<Self, BoltError> {
        if !c.is_finite() {
            return Err(BoltError::NotFinite {
                name: "C",
                value: c,
            });
        }
        if c <= 0.0 || c >= 1.0 {
            return Err(BoltError::StiffnessRatioRange { value: c });
        }
        Ok(Self(c))
    }

    /// The constant `C` itself (the fraction of external load the bolt
    /// picks up).
    pub fn value(self) -> f64 {
        self.0
    }

    /// The member share `1 - C` (the fraction by which the external load
    /// *relieves* the member clamping force).
    pub fn member_fraction(self) -> f64 {
        1.0 - self.0
    }
}

/// A fully specified bolted joint: a known preload, a known stiffness
/// share, and the nominal diameter used to relate torque and preload.
///
/// Build with [`BoltedJoint::from_torque`] (specify the tightening
/// torque and let the joint compute the preload) or
/// [`BoltedJoint::with_preload`] (specify the preload directly).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BoltedJoint {
    /// Achieved bolt preload tension `F` (N).
    preload_n: f64,
    /// Joint stiffness constant `C`.
    stiffness: StiffnessRatio,
    /// Nominal bolt diameter `d` (m).
    nominal_diameter_m: f64,
}

impl BoltedJoint {
    /// Build a joint by specifying the tightening torque, computing the
    /// preload `F = T / (K d)`.
    ///
    /// `torque_nm` is in newton-metres, `nominal_diameter_m` in metres;
    /// the resulting preload is in newtons.
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if the torque is negative / non-finite or
    /// the diameter is not strictly positive.
    pub fn from_torque(
        torque_nm: f64,
        nut_factor: NutFactor,
        nominal_diameter_m: f64,
        stiffness: StiffnessRatio,
    ) -> Result<Self, BoltError> {
        let t = BoltError::require_non_negative("torque_nm", torque_nm)?;
        let d = BoltError::require_positive("nominal_diameter_m", nominal_diameter_m)?;
        let preload_n = t / (nut_factor.value() * d);
        Ok(Self {
            preload_n,
            stiffness,
            nominal_diameter_m: d,
        })
    }

    /// Build a joint from a directly specified preload tension `F` (N).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if the preload is negative / non-finite or
    /// the diameter is not strictly positive.
    pub fn with_preload(
        preload_n: f64,
        nominal_diameter_m: f64,
        stiffness: StiffnessRatio,
    ) -> Result<Self, BoltError> {
        let f = BoltError::require_non_negative("preload_n", preload_n)?;
        let d = BoltError::require_positive("nominal_diameter_m", nominal_diameter_m)?;
        Ok(Self {
            preload_n: f,
            stiffness,
            nominal_diameter_m: d,
        })
    }

    /// Achieved preload tension `F` (N).
    pub fn preload_n(&self) -> f64 {
        self.preload_n
    }

    /// Joint stiffness constant `C`.
    pub fn stiffness(&self) -> StiffnessRatio {
        self.stiffness
    }

    /// Nominal bolt diameter `d` (m).
    pub fn nominal_diameter_m(&self) -> f64 {
        self.nominal_diameter_m
    }

    /// The tightening torque required to *reach* this joint's preload at
    /// a given nut factor, `T = K F d` (N·m). Inverse of
    /// [`BoltedJoint::from_torque`].
    pub fn required_torque_nm(&self, nut_factor: NutFactor) -> f64 {
        nut_factor.value() * self.preload_n * self.nominal_diameter_m
    }

    /// The *increment* in bolt tension caused by an external tensile
    /// load `P`: `ΔF_bolt = C P` (N).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if `external_load_n` is negative /
    /// non-finite.
    pub fn bolt_load_increment_n(&self, external_load_n: f64) -> Result<f64, BoltError> {
        let p = BoltError::require_non_negative("external_load_n", external_load_n)?;
        Ok(self.stiffness.value() * p)
    }

    /// The total tension carried by the bolt under an external load `P`,
    /// `F_bolt = F + C P` (N).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if `external_load_n` is negative /
    /// non-finite.
    pub fn bolt_load_n(&self, external_load_n: f64) -> Result<f64, BoltError> {
        Ok(self.preload_n + self.bolt_load_increment_n(external_load_n)?)
    }

    /// The residual clamping force holding the members together under an
    /// external load `P`, `F_clamp = F - (1 - C) P` (N).
    ///
    /// A positive value means the joint is still clamped; zero or
    /// negative means it has separated (the value is *not* clamped to
    /// zero here so callers can see how far past separation they are).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if `external_load_n` is negative /
    /// non-finite.
    pub fn clamping_force_n(&self, external_load_n: f64) -> Result<f64, BoltError> {
        let p = BoltError::require_non_negative("external_load_n", external_load_n)?;
        Ok(self.preload_n - self.stiffness.member_fraction() * p)
    }

    /// The external tensile load at which the members just separate,
    /// `P_sep = F / (1 - C)` (N).
    ///
    /// Loads above this gap the joint open. Because `C < 1`, the
    /// denominator `1 - C` is strictly positive, so this is always
    /// finite.
    pub fn separation_load_n(&self) -> f64 {
        self.preload_n / self.stiffness.member_fraction()
    }

    /// Whether an external load `P` keeps the joint clamped (`true`) or
    /// opens it (`false`). Separation is the boundary `P = P_sep`, which
    /// is reported as *not* clamped (residual clamp is zero).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if `external_load_n` is negative /
    /// non-finite.
    pub fn stays_clamped(&self, external_load_n: f64) -> Result<bool, BoltError> {
        Ok(self.clamping_force_n(external_load_n)? > 0.0)
    }

    /// The factor of safety against joint separation under a service
    /// load `P`, `n_sep = P_sep / P`.
    ///
    /// # Errors
    ///
    /// Returns [`BoltError`] if `external_load_n` is not strictly
    /// positive (a zero load has infinite margin and no meaningful
    /// ratio) or non-finite.
    pub fn separation_safety_factor(&self, external_load_n: f64) -> Result<f64, BoltError> {
        let p = BoltError::require_positive("external_load_n", external_load_n)?;
        Ok(self.separation_load_n() / p)
    }
}

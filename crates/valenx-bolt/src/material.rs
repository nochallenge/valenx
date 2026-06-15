//! Bolt material strength data.
//!
//! A bolt grade is captured by two stresses: the **proof strength**
//! `S_p` (the stress the bolt can carry with no measurable permanent
//! set — the design ceiling for preload) and the **ultimate tensile
//! strength** `S_u` (the stress at fracture). Multiplying a strength by
//! the tensile-stress area (see [`crate::stress`]) gives the
//! corresponding axial load.
//!
//! The built-in [`BoltGrade`] values are the nominal table strengths
//! from the metric property-class system (ISO 898-1). They are nominal
//! design figures, not measured lot certificates.

use serde::{Deserialize, Serialize};

use crate::error::BoltError;

/// A bolt material strength, in pascals.
///
/// Construct from a named [`BoltGrade`], or build a custom grade with
/// [`BoltMaterial::new`] (which validates that both strengths are
/// strictly positive and finite and that proof does not exceed
/// ultimate).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BoltMaterial {
    /// Proof strength `S_p` (Pa) — the preload design ceiling.
    pub proof_strength_pa: f64,
    /// Ultimate tensile strength `S_u` (Pa) — the fracture stress.
    pub tensile_strength_pa: f64,
}

impl BoltMaterial {
    /// Build a custom material from a proof and an ultimate strength
    /// (both in pascals).
    ///
    /// # Errors
    ///
    /// Returns [`BoltError::NonPositive`] / [`BoltError::NotFinite`] if
    /// either strength is not strictly positive and finite, and
    /// [`BoltError::NonPositive`] (named `proof_minus_tensile`) if the
    /// proof strength exceeds the ultimate strength — a physically
    /// impossible grade.
    pub fn new(proof_strength_pa: f64, tensile_strength_pa: f64) -> Result<Self, BoltError> {
        let proof = BoltError::require_positive("proof_strength_pa", proof_strength_pa)?;
        let tensile = BoltError::require_positive("tensile_strength_pa", tensile_strength_pa)?;
        if proof > tensile {
            // Encoded as a non-positive margin so the message reads
            // naturally (the ultimate-minus-proof margin must be >= 0).
            return Err(BoltError::NonPositive {
                name: "tensile_minus_proof",
                value: tensile - proof,
            });
        }
        Ok(Self {
            proof_strength_pa: proof,
            tensile_strength_pa: tensile,
        })
    }

    /// The recommended preload stress for a *reused* connection,
    /// `0.75 S_p` (Shigley's rule of thumb), in pascals.
    pub fn recommended_preload_stress_pa(&self) -> f64 {
        0.75 * self.proof_strength_pa
    }

    /// The recommended preload stress for a *permanent* connection,
    /// `0.90 S_p`, in pascals.
    pub fn recommended_permanent_preload_stress_pa(&self) -> f64 {
        0.90 * self.proof_strength_pa
    }
}

/// Common metric bolt property classes (ISO 898-1), plus a custom
/// escape hatch.
///
/// The strengths are the nominal class values; `8.8` etc. is read as
/// "ultimate ≈ first-number × 100 MPa, yield ≈ product of the two
/// numbers × 10 MPa", and proof strength is the tabulated value close
/// to yield.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum BoltGrade {
    /// Class 4.6 — low-carbon steel. `S_p` = 225 MPa, `S_u` = 400 MPa.
    Class4_6,
    /// Class 5.8 — `S_p` = 380 MPa, `S_u` = 520 MPa.
    Class5_8,
    /// Class 8.8 — the workhorse medium-carbon grade.
    /// `S_p` = 600 MPa, `S_u` = 830 MPa.
    Class8_8,
    /// Class 10.9 — alloy steel, quenched & tempered.
    /// `S_p` = 830 MPa, `S_u` = 1040 MPa.
    Class10_9,
    /// Class 12.9 — high-strength alloy steel.
    /// `S_p` = 970 MPa, `S_u` = 1220 MPa.
    Class12_9,
    /// A user-supplied grade with explicit proof / ultimate stresses
    /// (pascals).
    Custom {
        /// Proof strength `S_p` (Pa).
        proof_strength_pa: f64,
        /// Ultimate tensile strength `S_u` (Pa).
        tensile_strength_pa: f64,
    },
}

impl BoltGrade {
    /// Resolve this grade to a validated [`BoltMaterial`].
    ///
    /// The named classes always succeed (their tabulated values are
    /// valid by construction). [`BoltGrade::Custom`] is validated through
    /// [`BoltMaterial::new`] and may return a [`BoltError`].
    pub fn material(self) -> Result<BoltMaterial, BoltError> {
        // MPa → Pa.
        const M: f64 = 1.0e6;
        match self {
            BoltGrade::Class4_6 => BoltMaterial::new(225.0 * M, 400.0 * M),
            BoltGrade::Class5_8 => BoltMaterial::new(380.0 * M, 520.0 * M),
            BoltGrade::Class8_8 => BoltMaterial::new(600.0 * M, 830.0 * M),
            BoltGrade::Class10_9 => BoltMaterial::new(830.0 * M, 1040.0 * M),
            BoltGrade::Class12_9 => BoltMaterial::new(970.0 * M, 1220.0 * M),
            BoltGrade::Custom {
                proof_strength_pa,
                tensile_strength_pa,
            } => BoltMaterial::new(proof_strength_pa, tensile_strength_pa),
        }
    }

    /// Short label for UI / logs.
    pub fn label(self) -> &'static str {
        match self {
            BoltGrade::Class4_6 => "4.6",
            BoltGrade::Class5_8 => "5.8",
            BoltGrade::Class8_8 => "8.8",
            BoltGrade::Class10_9 => "10.9",
            BoltGrade::Class12_9 => "12.9",
            BoltGrade::Custom { .. } => "custom",
        }
    }
}

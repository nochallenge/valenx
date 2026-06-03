//! Provenance — who produced a `Results`, when, against what inputs.
//!
//! Every `Results` carries one `Provenance` value. It's what makes
//! reproducibility a mechanical check instead of folklore.

use serde::{Deserialize, Serialize};

/// SHA-256 digest as a hex string.
///
/// The actual hashing lives in `valenx-core`; this crate just stores
/// the result. Keeping the type opaque avoids a `sha2` dep leaking
/// into every downstream that touches results.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Sha256Hex(pub String);

impl Sha256Hex {
    /// Wrap a 64-character hex string as a [`Sha256Hex`].
    pub fn new(hex: impl Into<String>) -> Self {
        Self(hex.into())
    }
}

/// A reference to another `Results`' provenance — used for derived
/// fields whose source is another run.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProvenanceRef {
    /// UUID of the ancestor run.
    pub run_id: String,
    /// Human label ("raw pressure from run 2024-11-02").
    pub label: Option<String>,
}

/// Everything we captured about how a `Results` was produced.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Provenance {
    /// Identifier of the Valenx adapter that produced the result.
    pub adapter: String,
    /// Adapter crate version at the time of the run.
    pub adapter_version: String,
    /// Underlying external tool name (e.g. `"OpenFOAM"`, `"CalculiX"`).
    pub tool: String,
    /// Detected tool version string.
    pub tool_version: String,
    /// SHA-256 of the canonical case configuration.
    pub case_hash: Sha256Hex,
    /// SHA-256 of the input mesh.
    pub mesh_hash: Sha256Hex,
    /// SHA-256 of the generated solver input deck.
    pub input_hash: Sha256Hex,
    /// SHA-256 of the `tools.lock` snapshot.
    pub tools_lock_hash: Sha256Hex,
    /// UUID assigned to this run.
    pub run_id: String,
    /// Wall time the solve took, in seconds. Kept as `f64` rather
    /// than `Duration` so serialization is straightforward.
    pub wall_time_seconds: f64,
    /// ISO-8601 completion timestamp.
    pub completed_at: String,
    /// Ancestors for derived results.
    pub ancestors: Vec<ProvenanceRef>,
}

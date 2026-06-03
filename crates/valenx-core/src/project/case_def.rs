//! In-memory shape of per-case `case.toml` files (RFC 0001).
//!
//! Cases live in `cases/<name>/case.toml`. The schema below covers
//! the structure every adapter agrees on; physics-specific sections
//! (`[flow]`, `[structural]`, `[em]`, …) are stored as opaque TOML
//! tables so adapters can validate them against their own shapes.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// The root of `case.toml`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CaseDef {
    pub case: CaseHeader,

    /// Physics-specific sections stored verbatim. An adapter reads
    /// the sub-tables it cares about and ignores the rest. Keeping
    /// this as a `toml::Value`-like map means we don't have to
    /// exhaustively enumerate every physics schema here.
    #[serde(flatten)]
    pub sections: BTreeMap<String, toml::Value>,
}

/// Required `[case]` header of a `case.toml` file.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CaseHeader {
    /// Case-file SemVer.
    pub format: String,
    /// Display name shown in the UI's case picker.
    pub name: String,
    /// Physics domain: "cfd", "fea", "em", "chemistry", …
    pub physics: String,
    /// Adapter + native solver identifier, e.g. `"openfoam.simpleFoam"`.
    pub solver: String,
    /// Referenced mesh key (matches an entry under `[mesh.*]` in
    /// `project.toml`). `"default"` is the convention.
    pub mesh: String,
    /// Optional free-text description.
    #[serde(default)]
    pub description: Option<String>,
}

impl CaseDef {
    /// Access a physics section by name (e.g. "flow", "structural").
    pub fn section(&self, name: &str) -> Option<&toml::Value> {
        self.sections.get(name)
    }

    /// True if this case has physics section `name`.
    pub fn has_section(&self, name: &str) -> bool {
        self.sections.contains_key(name)
    }
}

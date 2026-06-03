//! In-memory shape of the per-project `tools.lock` file (RFC 0001).
//!
//! `tools.lock` pins the exact version + checksum + integration mode
//! of every external tool the project depends on. Opening a project
//! compares these entries against what's installed on the user's
//! machine; a mismatch becomes a user-facing warning, not a silent
//! substitution.

use serde::{Deserialize, Serialize};

use crate::license::LicenseMode;

/// Deserialised `tools.lock` snapshot — pins external tool versions
/// for reproducibility.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ToolsLock {
    /// Lock-file format SemVer.
    pub format: String,
    /// Identifier of the harness that wrote the lock (free text).
    #[serde(default)]
    pub generated_by: Option<String>,
    /// ISO-8601 timestamp when the lock was written.
    #[serde(default)]
    pub generated_at: Option<String>,

    /// Pinned tool entries (one per detected tool).
    #[serde(default, rename = "tool")]
    pub tools: Vec<ToolEntry>,
}

/// One pinned entry inside [`ToolsLock`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ToolEntry {
    /// Tool name (matches [`crate::adapter::AdapterInfo::id`]).
    pub name: String,
    /// Pinned version string.
    pub version: String,
    /// Hex SHA-256 prefixed `"sha256:"` per RFC 0001 convention.
    #[serde(default)]
    pub checksum: Option<String>,
    /// Distribution channel ("stable", "nightly", site-specific tag).
    #[serde(default)]
    pub channel: Option<String>,
    /// How Valenx talks to the tool (in-process, subprocess, remote).
    #[serde(default)]
    pub integration_mode: Option<LockedIntegrationMode>,
    /// Free-text license note recorded at lock time.
    #[serde(default)]
    pub license: Option<String>,
}

/// String-valued integration mode, mirrored as enum for type-safety.
/// Kept separate from [`crate::license::LicenseMode`] because this
/// one is loaded from files and uses kebab-case strings.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum LockedIntegrationMode {
    Bundled,
    DynamicLinked,
    Subprocess,
}

impl From<LockedIntegrationMode> for LicenseMode {
    fn from(m: LockedIntegrationMode) -> Self {
        match m {
            LockedIntegrationMode::Bundled => LicenseMode::Bundled,
            LockedIntegrationMode::DynamicLinked => LicenseMode::DynamicLinked,
            LockedIntegrationMode::Subprocess => LicenseMode::Subprocess,
        }
    }
}

impl ToolsLock {
    /// Look up a tool by name.
    pub fn get(&self, name: &str) -> Option<&ToolEntry> {
        self.tools.iter().find(|t| t.name == name)
    }
}

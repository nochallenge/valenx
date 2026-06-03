//! In-memory shape of `project.toml` (per RFC 0001).
//!
//! Every path in this module is **relative to the project root**
//! (the directory containing `project.toml`). The loader in
//! `super::loader` rejects absolute paths and paths escaping the
//! project.

use std::collections::BTreeMap;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// The canonical `project.toml` root. Additive fields get defaults
/// so old files written by newer apps still parse — per RFC 0001's
/// "unknown keys are preserved" commitment this is a best-effort
/// for known keys; unknown keys round-trip through the raw `extra`
/// map.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Project {
    pub project: ProjectHeader,
    #[serde(default)]
    pub units: UnitsConfig,
    #[serde(default)]
    pub geometry: GeometrySection,
    #[serde(default)]
    pub mesh: BTreeMap<String, MeshEntry>,
    #[serde(default)]
    pub cases: CasesSection,
    #[serde(default)]
    pub ui: UiSection,
}

/// `[project]` header block of `project.toml`.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectHeader {
    /// File-format SemVer (e.g. "1.0", "1.2").
    pub format: String,
    /// Display name of the project.
    pub name: String,
    /// Minimum Valenx version that can open this project.
    #[serde(default)]
    pub valenx_min: Option<String>,
    /// ISO-8601 creation timestamp.
    #[serde(default)]
    pub created: Option<String>,
    /// ISO-8601 last-modified timestamp.
    #[serde(default)]
    pub modified: Option<String>,
    /// Free-text author / owner field.
    #[serde(default)]
    pub author: Option<String>,
    /// Free-text description shown on the project landing page.
    #[serde(default)]
    pub description: Option<String>,
}

/// `[units]` block — the canonical unit symbols the project uses.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnitsConfig {
    /// Length unit symbol (default `"m"`).
    #[serde(default = "default_length")]
    pub length: String,
    /// Mass unit symbol (default `"kg"`).
    #[serde(default = "default_mass")]
    pub mass: String,
    /// Time unit symbol (default `"s"`).
    #[serde(default = "default_time")]
    pub time: String,
    /// Temperature unit symbol (default `"K"`).
    #[serde(default = "default_temperature")]
    pub temperature: String,
}

fn default_length() -> String {
    "m".to_string()
}
fn default_mass() -> String {
    "kg".to_string()
}
fn default_time() -> String {
    "s".to_string()
}
fn default_temperature() -> String {
    "K".to_string()
}

/// `[geometry]` block — the project's registered geometry sources.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeometrySection {
    /// Geometry entries (one per source file).
    #[serde(default)]
    pub entries: Vec<GeometryEntry>,
}

/// One geometry registration.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeometryEntry {
    /// Stable identifier the case files reference.
    pub id: String,
    /// Relative path to the source geometry file.
    pub source: PathBuf,
    /// Format identifier (`"step"`, `"iges"`, `"stl"`, …).
    pub format: String,
}

/// `[mesh.<key>]` block — one named mesh asset.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeshEntry {
    /// Relative path to the mesh file under the project root.
    pub source: PathBuf,
    /// SHA-256 of the meshing config that produced this mesh, if
    /// known.
    #[serde(default)]
    pub config_hash: Option<String>,
}

/// `[cases]` block — declared case order + bookkeeping.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CasesSection {
    /// Stable display order for cases in the UI.
    #[serde(default)]
    pub order: Vec<String>,
}

/// User-intent UI state that survives save/load (per RFC 0001 § 6).
/// Transient UI (panel layout, scroll position) does not live here.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UiSection {
    #[serde(default)]
    pub last_camera_pivot: Option<[f64; 3]>,
}

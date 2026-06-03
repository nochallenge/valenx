//! The top-level `Results` type returned by every adapter's
//! `collect()`. Pulls together metadata, field catalog, scalars,
//! artifacts, and provenance.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::artifact::Artifact;
use crate::catalog::{FieldCatalog, ScalarCatalog};
use crate::provenance::Provenance;

/// Top-line descriptive metadata for a results bundle.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ResultMeta {
    /// Case this result was produced from (stable project-local ID).
    pub case_id: String,
    /// Short human label ("CFD steady, airfoil, Re=6e6").
    pub description: Option<String>,
}

/// Canonical results produced by one adapter run.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Results {
    pub meta: ResultMeta,
    pub fields: FieldCatalog,
    pub scalars: ScalarCatalog,
    pub artifacts: Vec<Artifact>,
    pub provenance: Provenance,
    /// Path to the on-disk manifest (`results/manifest.toml`), used
    /// by the UI to resolve relative artifact paths.
    pub manifest_path: Option<PathBuf>,
}

impl Results {
    /// A minimal constructor useful for tests and adapter scaffolds.
    pub fn empty(case_id: &str, provenance: Provenance) -> Self {
        Self {
            meta: ResultMeta {
                case_id: case_id.to_string(),
                description: None,
            },
            fields: FieldCatalog::new(),
            scalars: ScalarCatalog::new(),
            artifacts: Vec::new(),
            provenance,
            manifest_path: None,
        }
    }
}

//! Streaming docking events for LLM-driven orchestration.
//!
//! Every long-running step in [`crate::runner::dock_with_events`] emits
//! one of these. Serialized as NDJSON (one JSON object per line) when
//! piped to stdout; the desktop UI subscribes to the in-process callback
//! variant and never touches the disk.

use serde::Serialize;

/// One step of the docking run. `tag` is the discriminator the LLM keys on.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum DockEvent {
    /// Inputs validated, search about to begin.
    Started {
        /// Receptor PDBQT path (or `"<inline>"`).
        receptor: String,
        /// Ligand PDBQT path (or `"<inline>"`).
        ligand: String,
        /// Ligand heavy-atom count (excludes H/HD).
        ligand_heavy_atoms: usize,
        /// Number of rotatable bonds.
        n_torsions: usize,
    },
    /// Receptor grid bundle finished precomputing.
    GridBuilt {
        /// Number of distinct atom types gridded.
        n_grids: usize,
        /// Total grid voxel count across all atom-type maps.
        total_voxels: usize,
    },
    /// Periodic progress during the search.
    SearchProgress {
        /// Fraction in [0, 1].
        fraction: f64,
        /// Best score so far (kcal/mol).
        best_score: f64,
    },
    /// One pose has been added to the result set.
    PoseFound {
        /// 1-based rank (1 = best).
        rank: usize,
        /// Score in kcal/mol.
        score: f64,
        /// Heavy-atom RMSD to the top pose (0.0 if rank=1).
        rmsd_to_top: f64,
    },
    /// Search complete; ready for output.
    Complete {
        /// Total poses kept after clustering + cutoff.
        n_poses: usize,
        /// Wall-clock seconds.
        wall_seconds: f64,
    },
    /// Fatal error — the run will not produce results.
    Error {
        /// Stable code from [`crate::error::DockError::code`].
        code: &'static str,
        /// Category from [`crate::error::DockError::category`].
        category: String,
        /// Human-readable detail.
        message: String,
    },
}

impl DockEvent {
    /// Serialize to a single line of NDJSON (no trailing newline).
    ///
    /// The fallback path on serialize failure still produces valid
    /// JSON: we wrap the underlying error message through
    /// `serde_json::Value::String` so any quotes or backslashes are
    /// escaped properly. Hand-formatting the fallback with naive
    /// string interpolation (the previous behaviour) would emit
    /// malformed JSON whenever the error contained a `"`.
    pub fn to_ndjson(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|e| {
            let safe = serde_json::Value::String(e.to_string()).to_string();
            format!(
                r#"{{"type":"error","code":"dock.serialize","category":"Internal","message":{safe}}}"#
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn started_event_serializes_to_expected_tag() {
        let e = DockEvent::Started {
            receptor: "rec.pdbqt".into(),
            ligand: "lig.pdbqt".into(),
            ligand_heavy_atoms: 10,
            n_torsions: 3,
        };
        let line = e.to_ndjson();
        assert!(line.starts_with("{"));
        assert!(line.contains("\"type\":\"started\""));
        assert!(line.contains("\"ligand_heavy_atoms\":10"));
        assert!(line.contains("\"n_torsions\":3"));
    }

    #[test]
    fn pose_found_serializes_rank_and_score() {
        let e = DockEvent::PoseFound {
            rank: 1,
            score: -9.345,
            rmsd_to_top: 0.0,
        };
        let line = e.to_ndjson();
        assert!(line.contains("\"type\":\"pose_found\""));
        assert!(line.contains("\"rank\":1"));
        assert!(line.contains("-9.345"));
    }
}

//! # valenx-adapter-prody
//!
//! Adapter for [ProDy](http://prody.csb.pitt.edu/) — the canonical
//! Python library for protein dynamics. ProDy ships elastic-network
//! models (ENM / GNM / ANM), normal-mode analysis, ensemble PCA, the
//! NMD trajectory format consumed by VMD's NMWiz plug-in, and
//! integrations with the BLAST / DALI / PDB databases.
//!
//! **Phase 5.5 — subprocess wrapper for user-provided ProDy
//! scripts.** The user authors an `analyse.py` referenced from
//! `[bio.prody].script` in `case.toml`, plus the input `.pdb` to
//! operate on, an output basename, the number of modes, and the
//! contact cutoff. `prepare()` stages the script + input PDB into
//! the workdir and writes a flat `valenx_params.json` so the
//! script can read the parsed knobs without re-parsing case.toml.
//! `run()` invokes `python <script>` via the shared subprocess
//! runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "input_pdb":       "1ake.pdb",
//!   "output_basename": "modes",
//!   "num_modes":       20,
//!   "cutoff":          15.0
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to ProDy themselves. This keeps the
//! adapter free of upstream API churn.
//!
//! On `collect()` we walk the workdir for `<basename>*.npz`
//! (ProDy's NumPy mode storage), `<basename>*.nmd` (NMD trajectory
//! consumed by NMWiz / VMD), and `<basename>*.csv` (analysis
//! tables).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ProdyInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ProdyAdapter::new())
}

pub struct ProdyAdapter;

impl ProdyAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ProdyAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "prody";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ProdyAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ProDy",
            // ProDy 2.4+ (2023) ships the modern ANM / GNM API
            // surface and the NMD writer we lean on. Upper bound
            // 3.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 4, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "http://prody.csb.pitt.edu/",
            homepage_url: "http://prody.csb.pitt.edu/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import prody` to confirm the library is
                // installed in the chosen interpreter (vs. just
                // having Python on PATH). Failure is non-fatal —
                // surfaced via a warning so the user can still
                // configure a different interpreter via `python`.
                let importable = prody_importable(&binary_path);
                let mut warnings = Vec::new();
                if !importable {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `prody` — install with `pip install prody` (or \
                         `conda install -c conda-forge prody`) for runs to \
                         succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with ProDy installed; `pip install prody` \
                       after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ProdyInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.prody].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script into the workdir.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.prody].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }
        let script_filename =
            input
                .script
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.prody].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the input `.pdb` so the script can resolve it via
        // a bare relative filename inside the workdir.
        let source_pdb = confined_join(&case.path, &input.input_pdb)?;
        if !source_pdb.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.prody].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }
        let pdb_filename =
            input
                .input_pdb
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.prody].input_pdb path `{}` has no filename",
                        input.input_pdb.display()
                    ),
                })?;
        let dest_pdb = workdir.join(pdb_filename);
        if source_pdb != dest_pdb {
            fs::copy(&source_pdb, &dest_pdb)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's script can read the parsed `[bio.prody]` knobs
        // without having to reparse case.toml itself. Built by
        // hand to avoid pulling in serde_json for a 4-key flat
        // object.
        let params_json = format!(
            "{{\n  \"input_pdb\": {},\n  \"output_basename\": {},\n  \"num_modes\": {},\n  \"cutoff\": {}\n}}\n",
            json_string(&pdb_filename.to_string_lossy()),
            json_string(&input.output_basename),
            input.num_modes,
            input.cutoff,
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Python-script adapter: bare `python` / `python3` walks
        // PATH; absolute paths or pinned interpreters are honored
        // verbatim.
        // Round-4 security: validate python interpreter spec
        // against the allow-list AND resolve to a real binary
        // in one step. Closes the arbitrary-binary-exec class
        // that round-3 only patched in 8 of the 48 affected
        // adapters.
        let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
            &input.python,
            PYTHON_BINARIES,
        )
        // Round-5: do NOT rewrap as ToolNotInstalled — the resolver
        // returns InvalidCase for allow-list rejections (which a hint
        // string would have hidden) and a clear Other for PATH lookup
        // failures. Pass the error through unchanged.
        ?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // ENM diagonalisation on a small protein finishes in
            // seconds; whole-ribosome ANM / large-ensemble PCA can
            // run for an hour or more. 2 hours is a generous
            // default that covers the long tail.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ProDy", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] prody done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] prody done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("Error") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash the staged Python script as the
        // canonical input descriptor. Falls back to case.toml
        // when no script is staged yet (partial / failed runs).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ProDy",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict NPZ / NMD / CSV outputs to those whose stem
        // starts with the configured `output_basename`. Failure
        // to read params is non-fatal — accept everything.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-prody", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_ok = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            if !stem_ok {
                continue;
            }
            let (kind, label) = match ext.as_deref() {
                Some("npz") => (ArtifactKind::Native, "ProDy ENM modes".to_string()),
                Some("nmd") => (ArtifactKind::Native, "ProDy NMD trajectory".to_string()),
                Some("csv") => (ArtifactKind::Tabular, "ProDy analysis".to_string()),
                _ => continue,
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.prody.analyze"],
        }
    }
}

/// Probe whether `import prody` succeeds in the given Python
/// interpreter. Returns true only when Python exits 0.
fn prody_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import prody")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output();
    matches!(output, Ok(o) if o.status.success())
}

/// Escape a string for embedding inside a JSON string literal.
fn json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000C}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("py"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Pull `output_basename` out of our own hand-emitted
/// `valenx_params.json` for collect()-time output filtering.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. We wrote the file ourselves so we know
/// its shape; a full JSON parser would be overkill.
fn extract_json_string(text: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\":");
    let idx = text.find(&needle)?;
    let rest = &text[idx + needle.len()..];
    let start = rest.find('"')? + 1;
    let body = &rest[start..];
    let mut out = String::new();
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                '"' => out.push('"'),
                '\\' => out.push('\\'),
                'n' => out.push('\n'),
                'r' => out.push('\r'),
                't' => out.push('\t'),
                'b' => out.push('\u{0008}'),
                'f' => out.push('\u{000C}'),
                other => out.push(other),
            },
            c => out.push(c),
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = ProdyAdapter::new().info();
        assert_eq!(info.id, "prody");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "ProDy");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ProdyAdapter::new().info();
        // ProDy 2.4+ is the modern stable line; upper bound 3.0
        // reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ProdyAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.prody.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ProdyAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}

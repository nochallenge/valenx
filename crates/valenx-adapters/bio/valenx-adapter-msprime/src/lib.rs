//! # valenx-adapter-msprime
//!
//! Adapter for [msprime](https://tskit.dev/msprime/) — Jerome
//! Kelleher's coalescent backwards-in-time population-genetics
//! simulator. msprime simulates the ancestry of a sample under a
//! configurable demography and recombination map, then layers
//! mutations onto the resulting tree sequence. It is the speed-of-
//! light coalescent simulator (millions of samples per minute on a
//! workstation) and the canonical companion to SLiM (forward-time)
//! and tskit (tree-sequence analysis).
//!
//! **Phase 29 — subprocess wrapper for user-provided msprime
//! scripts.** The user authors a `simulate.py` referenced from
//! `[bio.msprime].script` in `case.toml`, plus the demographic
//! knobs (`population_size`, `num_samples`, `recombination_rate`,
//! `mutation_rate`, `output_basename`). `prepare()` stages the
//! script into the workdir and writes a flat `valenx_params.json`
//! so the script can read the parsed knobs without re-parsing
//! `case.toml`. `run()` invokes `python <script>` via the shared
//! subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "population_size":    10000,
//!   "num_samples":        100,
//!   "recombination_rate": 1e-08,
//!   "mutation_rate":      1.5e-08,
//!   "output_basename":    "sim"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to msprime themselves. This keeps
//! the adapter free of upstream API churn — msprime's Python
//! signature has shifted between minor releases.
//!
//! On `collect()` we walk the workdir for `<basename>.trees`
//! (tskit tree-sequence native), `<basename>.vcf` (genotype calls),
//! and `<basename>.csv` (per-sample summaries the user's script
//! emitted via tskit's tabular APIs).

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

use crate::case_input::MsprimeInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MsprimeAdapter::new())
}

pub struct MsprimeAdapter;

impl MsprimeAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MsprimeAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "msprime";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for MsprimeAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "msprime",
            // msprime 1.3 (2024) is the floor we test against — it
            // ships the modern `sim_ancestry()` / `sim_mutations()`
            // split and the tskit 0.5+ tree-sequence format we
            // surface in collect(). Upper bound 2.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 3, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://tskit.dev/msprime/",
            homepage_url: "https://tskit.dev/msprime/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import msprime` to confirm the library is
                // installed in the chosen interpreter (vs. just
                // having Python on PATH). Failure is non-fatal —
                // surfaced via a warning so the user can still
                // configure a different interpreter via `python`.
                let importable = msprime_importable(&binary_path);
                let mut warnings = Vec::new();
                if !importable {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `msprime` — install with `pip install msprime` (or \
                         `conda install -c conda-forge msprime`) for runs \
                         to succeed"
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
                hint: "Python 3.9+ with msprime installed; `pip install \
                       msprime` after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MsprimeInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.msprime].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against
        // the case directory; `confined_join` rejects absolute paths
        // and `..` traversal so the staged copy stays confined to the
        // case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.msprime].script `{}` not found (resolved {})",
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
                        "[bio.msprime].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's script can read the parsed `[bio.msprime]` knobs
        // without having to reparse case.toml itself. Built by
        // hand to avoid pulling in serde_json for a 5-key flat
        // object. Floats use `{:e}` so Python's `json.load` parses
        // them back as floats rather than strings.
        let params_json = format!(
            "{{\n  \"population_size\": {},\n  \"num_samples\": {},\n  \"recombination_rate\": {:e},\n  \"mutation_rate\": {:e},\n  \"output_basename\": {}\n}}\n",
            input.population_size,
            input.num_samples,
            input.recombination_rate,
            input.mutation_rate,
            json_string(&input.output_basename),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Phase 17 Python-script adapter: bare `python` / `python3`
        // walks PATH; absolute paths or pinned interpreters are
        // honored verbatim.
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
            // Coalescent simulations are fast — typical runs
            // finish in seconds-to-minutes on a workstation.
            // Whole-genome runs with millions of samples can run
            // for a few hours. 4 hours is a generous default that
            // covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting msprime", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] msprime done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] msprime done") {
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
            "msprime",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged params back out so we can restrict the
        // collected outputs to those whose stem matches the
        // configured `output_basename`. Failure to read the params
        // is non-fatal — accept everything in that case.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-msprime", ?e, "workdir read failed");
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
            // Restrict to outputs whose stem matches the
            // configured `output_basename`. If params couldn't be
            // read, accept everything (best-effort).
            let stem_ok = match basename.as_deref() {
                Some(b) => stem == b,
                None => true,
            };
            if !stem_ok {
                continue;
            }
            let (kind, label) = match ext.as_deref() {
                Some("trees") => (ArtifactKind::Native, "msprime tree sequence".to_string()),
                Some("vcf") => (ArtifactKind::Tabular, "msprime VCF genotypes".to_string()),
                Some("csv") => (ArtifactKind::Tabular, "msprime sample summary".to_string()),
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
            ribbon_contributions: vec!["bio.msprime.simulate"],
        }
    }
}

/// Probe whether `import msprime` succeeds in the given Python
/// interpreter. Returns true only when Python exits 0 — any
/// failure (interpreter missing, library not on sys.path, import-
/// time exception) is treated as "not importable" and surfaced via
/// a probe warning rather than blocking adapter use.
fn msprime_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import msprime")
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
/// `valenx_params.json`. Trivially small — we wrote the file
/// ourselves so we know its shape; a full JSON parser would be
/// overkill.
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
        let info = MsprimeAdapter::new().info();
        assert_eq!(info.id, "msprime");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "msprime");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MsprimeAdapter::new().info();
        // msprime 1.3+ (2024) ships the modern `sim_ancestry()` /
        // `sim_mutations()` split; upper bound 2.0 reserves room
        // for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 3, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MsprimeAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.msprime.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MsprimeAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}

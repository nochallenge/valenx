//! # valenx-adapter-rfantibody
//!
//! Adapter for [RFantibody](https://github.com/RosettaCommons/RFantibody)
//! — RosettaCommons' antibody-specific diffusion design model. Where
//! RFdiffusion handles general protein backbones, RFantibody specializes
//! in the antibody case: take an antibody framework + an antigen target,
//! redesign one or more CDR loops (H1/H2/H3/L1/L2/L3) so the antibody
//! binds the target.
//!
//! **Phase 27.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `design_rfantibody.py` (or whatever filename)
//! referenced from `[bio.rfantibody].script` in `case.toml` plus a
//! framework PDB and a target PDB. `prepare()` stages the script + both
//! PDBs into the workdir and `run()` invokes `python <script>` via the
//! shared subprocess runner. The script is responsible for invoking
//! RFantibody's inference entry point with the parsed knobs and writing
//! PDB outputs under `<output_basename>*.pdb`.
//!
//! ## `valenx_params.json`
//!
//! RFantibody's CLI is a Hydra config tree built on top of
//! RFdiffusion's — there's no stable flag shape we can pin without
//! breaking on the next refactor. Instead, `prepare()` writes a flat
//! JSON file `valenx_params.json` into the workdir containing the
//! parsed `[bio.rfantibody]` knobs:
//!
//! ```json
//! {
//!   "framework_pdb":   "framework.pdb",
//!   "target_pdb":      "antigen.pdb",
//!   "design_loops":    ["H3"],
//!   "num_designs":     8,
//!   "diffusion_steps": 50,
//!   "output_basename": "design"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to RFantibody themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.pdb`
//! files and parse each via [`valenx_bio::format::pdb::read`]. Each is
//! surfaced as a typed [`ArtifactKind::Native`] artifact with an
//! "RFantibody design" label.

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

use crate::case_input::RfAntibodyInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RfAntibodyAdapter::new())
}

pub struct RfAntibodyAdapter;

impl RfAntibodyAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RfAntibodyAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rfantibody";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for RfAntibodyAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RFantibody",
            // RFantibody's first tagged release is the 1.0 weights /
            // inference code drop. Upper bound 2.0 reserves room for
            // an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://github.com/RosettaCommons/RFantibody",
            homepage_url: "https://github.com/RosettaCommons/RFantibody",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import rfantibody; print(rfantibody.__version__)`
                // — confirms the package is importable from the chosen
                // interpreter (vs. just having Python on PATH). Some
                // RFantibody checkouts don't expose `__version__` —
                // fall back to a "couldn't import" warning so the probe
                // still surfaces a useful state.
                let found_version = detect_rfantibody_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `rfantibody` — install RFantibody from \
                         https://github.com/RosettaCommons/RFantibody and \
                         ensure it's importable from the chosen interpreter \
                         for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with RFantibody installed; clone \
                       https://github.com/RosettaCommons/RFantibody and \
                       follow the install steps after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RfAntibodyInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.rfantibody].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. `confined_join`
        // rejects absolute paths and `..` traversal so the staged copy
        // stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rfantibody].script `{}` not found (resolved {})",
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
                        "[bio.rfantibody].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the framework PDB.
        let framework_filename = stage_pdb(case, workdir, &input.framework_pdb, "framework_pdb")?;
        // Stage the target PDB.
        let target_filename = stage_pdb(case, workdir, &input.target_pdb, "target_pdb")?;

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's design script can read the parsed `[bio.rfantibody]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in a serde_json dep for a 6-key flat
        // object.
        let loops_array = format!(
            "[{}]",
            input
                .design_loops
                .iter()
                .map(|s| json_string(s))
                .collect::<Vec<_>>()
                .join(", ")
        );
        let params_json = format!(
            "{{\n  \"framework_pdb\": {},\n  \"target_pdb\": {},\n  \"design_loops\": {},\n  \"num_designs\": {},\n  \"diffusion_steps\": {},\n  \"output_basename\": {}\n}}\n",
            json_string(&framework_filename.to_string_lossy()),
            json_string(&target_filename.to_string_lossy()),
            loops_array,
            input.num_designs,
            input.diffusion_steps,
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
            // RFantibody sampling can take hours per design on a
            // consumer GPU; with `num_designs` typically 8 – 64 the
            // tail is long. 4 hours mirrors the RFdiffusion default
            // — generous enough that long runs aren't pre-empted.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RFantibody", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] rfantibody done` to signal completion
            // before exit; lift to a 95% progress tick.
            if line.contains("[valenx] rfantibody done") {
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
        // Provenance: hash the staged framework PDB. Falls back to
        // the script, then case.toml, when the PDB isn't present yet.
        let framework_path = read_params_pdb(&job.workdir, "framework_pdb")
            .map(|name| job.workdir.join(name))
            .filter(|p| p.is_file());
        let target_path = read_params_pdb(&job.workdir, "target_pdb")
            .map(|name| job.workdir.join(name))
            .filter(|p| p.is_file());
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = framework_path
            .clone()
            .or_else(|| target_path.clone())
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RFantibody",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged params back out so we can restrict the
        // collected output PDBs to those matching `output_basename`
        // and label the input PDBs distinctly. Failure to read the
        // params is non-fatal — collect still surfaces every PDB.
        let basename = read_params(&job.workdir);

        if let Some(p) = framework_path.clone() {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RFantibody framework PDB".to_string(),
            });
        }
        if let Some(p) = target_path.clone() {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RFantibody target PDB".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RFantibody script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-rfantibody", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut design_paths: Vec<PathBuf> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            match ext.as_deref() {
                Some("pdb") => {
                    // Skip the staged input PDBs — they're already
                    // surfaced above.
                    if Some(&path) == framework_path.as_ref() || Some(&path) == target_path.as_ref()
                    {
                        continue;
                    }
                    // Restrict to designs whose stem starts with the
                    // configured `output_basename`. If params couldn't
                    // be read, accept everything (best-effort).
                    let stem_ok = match basename.as_ref() {
                        Some(b) => {
                            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                            stem.starts_with(b.as_str())
                        }
                        None => true,
                    };
                    if stem_ok {
                        design_paths.push(path);
                    }
                }
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "RFantibody log".to_string(),
                }),
                _ => continue,
            }
        }
        design_paths.sort();
        for path in design_paths {
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("design")
                .to_string();
            // Round-22 M2: cap the per-PDB read at MAX_PDB_FILE_BYTES
            // (256 MiB) so a poisoned workdir with a multi-GB `.pdb`
            // can't OOM `collect()` before the parser runs.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
            ) {
                Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                    Ok(structure) => format!(
                        "RFantibody design `{}` ({} atoms, {} residues)",
                        stem,
                        structure.atom_count(),
                        structure.residue_count()
                    ),
                    Err(_) => "RFantibody design".to_string(),
                },
                Err(_) => "RFantibody design".to_string(),
            };
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
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
            ribbon_contributions: vec!["bio.rfantibody.design"],
        }
    }
}

/// Stage one of the input PDBs into the workdir. Returns the
/// destination filename (just the leaf, for embedding in
/// `valenx_params.json`).
fn stage_pdb(
    case: &Case,
    workdir: &Path,
    rel_path: &Path,
    field_name: &str,
) -> Result<PathBuf, AdapterError> {
    let source = confined_join(&case.path, rel_path)?;
    if !source.is_file() {
        return Err(AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!(
                "[bio.rfantibody].{field_name} `{}` not found (resolved {})",
                rel_path.display(),
                source.display()
            ),
        });
    }
    let filename = rel_path
        .file_name()
        .ok_or_else(|| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!(
                "[bio.rfantibody].{field_name} path `{}` has no filename",
                rel_path.display()
            ),
        })?;
    let dest = workdir.join(filename);
    if source != dest {
        fs::copy(&source, &dest)?;
    }
    Ok(PathBuf::from(filename))
}

/// Escape a string for embedding inside a JSON string literal.
/// Avoids pulling in a serde_json dependency for the trivially small
/// `valenx_params.json` we emit.
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

/// Read just the `output_basename` field from `valenx_params.json` so
/// `collect()` can restrict design-PDB pickup to the configured stem.
fn read_params(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Read one of the staged PDB filenames (`framework_pdb` or
/// `target_pdb`) from `valenx_params.json` so `collect()` can label
/// them distinctly.
fn read_params_pdb(workdir: &Path, key: &str) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, key)
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. Trivially small — we wrote the file
/// ourselves so we know its shape.
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

/// Run `python -c "import rfantibody; print(rfantibody.__version__)"`
/// and parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, rfantibody not importable, version
/// string malformed); `probe()` falls back to a "rfantibody not
/// importable" warning in that case.
fn detect_rfantibody_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import rfantibody; print(rfantibody.__version__)")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    let raw = valenx_core::adapter_helpers::extract_semver(&stdout)?;
    let dots = raw.chars().filter(|c| *c == '.').count();
    let normalised: String = match dots {
        0 => format!("{raw}.0.0"),
        1 => format!("{raw}.0"),
        _ => raw,
    };
    Version::parse(&normalised).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal valid PDB record covering one residue.
    const SAMPLE_PDB: &str = "\
ATOM      1  N   ALA A   1      11.104  13.207   2.063  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.804  13.793   3.215  1.00  0.00           C
ATOM      3  C   ALA A   1      11.072  15.058   3.668  1.00  0.00           C
ATOM      4  O   ALA A   1       9.835  15.117   3.586  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.916  12.789   4.357  1.00  0.00           C
END
";

    #[test]
    fn info_is_bio_domain() {
        let info = RfAntibodyAdapter::new().info();
        assert_eq!(info.id, "rfantibody");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "RFantibody");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RfAntibodyAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RfAntibodyAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rfantibody.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RfAntibodyAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` should walk the workdir for `<output_basename>*.pdb`
    /// design files, surface them as Native artifacts with parsed
    /// labels, and tag the framework + target PDBs as auxiliary
    /// `Other` artifacts with distinct labels.
    #[test]
    fn collect_walks_workdir_and_filters_designs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rfantibody-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Stage the params + input PDBs + script + designs.
        fs::write(
            tmp.join("valenx_params.json"),
            "{\n  \"framework_pdb\": \"framework.pdb\",\n  \"target_pdb\": \"antigen.pdb\",\n  \"design_loops\": [\"H3\"],\n  \"num_designs\": 2,\n  \"diffusion_steps\": 50,\n  \"output_basename\": \"design\"\n}\n",
        )
        .unwrap();
        fs::write(tmp.join("framework.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("antigen.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("design_0.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design_1.pdb"), SAMPLE_PDB).unwrap();
        // A stray pdb that shouldn't be picked up as a design.
        fs::write(tmp.join("unrelated.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("run.log"), b"RFantibody run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = RfAntibodyAdapter::new().collect(&job).unwrap();

        // Two design PDBs picked up.
        let designs: Vec<_> = results
            .artifacts
            .iter()
            .filter(|a| a.kind == ArtifactKind::Native)
            .collect();
        assert_eq!(
            designs.len(),
            2,
            "expected 2 design PDBs, got {}: {:?}",
            designs.len(),
            results.artifacts
        );
        for d in &designs {
            assert!(
                d.label.contains("RFantibody design"),
                "label was: {}",
                d.label
            );
        }

        // Framework + target tagged as Other with distinct labels.
        let framework_art = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| s == "framework.pdb")
                    .unwrap_or(false)
            })
            .expect("framework PDB artifact present");
        assert_eq!(framework_art.kind, ArtifactKind::Other);
        assert_eq!(framework_art.label, "RFantibody framework PDB");

        let target_art = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| s == "antigen.pdb")
                    .unwrap_or(false)
            })
            .expect("target PDB artifact present");
        assert_eq!(target_art.kind, ArtifactKind::Other);
        assert_eq!(target_art.label, "RFantibody target PDB");

        let py_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "py"))
            .expect("script artifact present");
        assert_eq!(py_art.kind, ArtifactKind::Other);
        assert_eq!(py_art.label, "RFantibody script");

        let log_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "log"))
            .expect("log artifact present");
        assert_eq!(log_art.kind, ArtifactKind::Log);

        let _ = fs::remove_dir_all(&tmp);
    }
}

//! # valenx-adapter-tcoffee
//!
//! Adapter for [T-Coffee](https://tcoffee.org/) — Notredame &
//! Higgins's consensus-based multiple-sequence aligner. Where MAFFT and
//! Clustal Omega are pure progressive aligners, T-Coffee builds a
//! library of pairwise alignments from multiple sources (local +
//! global, structure-aware, profile-driven) and uses that library as a
//! consistency objective for the progressive step. This makes it
//! noticeably more accurate on hard cases — distantly related
//! sequences, low coverage — at the cost of being substantially slower.
//!
//! **Phase 18.7 — subprocess wrapper around `t_coffee`.** The user
//! supplies a multi-FASTA via `[bio.tcoffee]` in `case.toml`.
//! `prepare()` resolves it against the case directory and composes
//! `t_coffee <input> -output=<outfmt> -outfile=<basename>.aln
//! [-mode=<mode>] [extras...]`. T-Coffee uses `=`-style flags
//! exclusively (not space-separated), and the output filename always
//! ends in `.aln` regardless of the requested format — `outfmt`
//! controls the file *contents*, not the extension.
//!
//! Alongside the alignment, T-Coffee writes a guide tree
//! (`*.dnd`, NEWICK) and a log; `collect()` surfaces all three
//! categories.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::TCoffeeInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(TCoffeeAdapter::new())
}

pub struct TCoffeeAdapter;

impl TCoffeeAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for TCoffeeAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "tcoffee";
/// T-Coffee's binary name. The canonical name is `t_coffee`
/// (underscore, not hyphen) everywhere — Bioconda, Homebrew, Debian,
/// and the upstream tarball all install under that exact name.
const BINARIES: &[&str] = &["t_coffee"];

impl Adapter for TCoffeeAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "T-Coffee",
            // T-Coffee's modern series is the 13.x line, which has been
            // the upstream stable since ~2019 and is what every distro
            // and Bioconda ship today. Floor at 13.0.0; upper bound
            // 14.0.0 reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(13, 0, 0),
                max_exclusive: Version::new(14, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "https://tcoffee.readthedocs.io/",
            homepage_url: "https://tcoffee.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `t_coffee -version` prints something like
                // "PROGRAM: T-COFFEE Version_13.46.0.919e8c6b" — the
                // semver detector pulls the leading 13.46.0 cleanly via
                // its built-in version-string regex. If the parse fails
                // the probe still succeeds (binary on PATH is enough);
                // the upstream banner format has been stable for years.
                let found_version = detect_tool_version_semver(&binary_path, &["-version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "T-Coffee 13.0.0+ required; install via `apt install t-coffee`, \
                       `brew install t-coffee`, or `conda install -c bioconda t-coffee`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = TCoffeeInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.tcoffee].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the input FASTA against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "seqs.fa"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.tcoffee].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "T-Coffee 13.0.0+ required; install via `apt install t-coffee`, \
                       `brew install t-coffee`, or `conda install -c bioconda t-coffee`"
                .into(),
        })?;

        // Compose the `t_coffee` invocation. T-Coffee uses `=`-style
        // flags exclusively (not space-separated), and the output file
        // is always `<basename>.aln` regardless of `outfmt` — the
        // format affects the file *content*, not the extension. That
        // makes `collect()` straightforward: walk for `<basename>*` and
        // pick up the alignment plus anything else T-Coffee dropped
        // alongside it.
        let output_filename = format!("{}.aln", input.output_basename);

        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            source_input.into_os_string(),
            OsString::from(format!("-output={}", input.outfmt)),
            OsString::from(format!("-outfile={output_filename}")),
        ];
        if let Some(mode) = &input.mode {
            native_command.push(OsString::from(format!("-mode={mode}")));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // T-Coffee is meaningfully slower than Clustal Omega or
            // MAFFT — its consistency-library construction can dominate
            // wall time on inputs in the hundreds of sequences, and
            // structure-aware modes (`expresso`, `psicoffee`) issue
            // network calls to PDB / NCBI. 4 hours covers the long tail
            // without being absurd, matching the other MSA adapters.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting T-Coffee", |line| {
            let mut hint = subprocess::Hint::default();
            // T-Coffee is chatty by default; explicit `ERROR` /
            // `FATAL` markers are the ones worth lifting as warnings.
            // Anything else is progress reporting we don't need to
            // surface to the run report.
            if line.contains("ERROR") || line.contains("FATAL") {
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
        // Surface three categories: anything starting with the
        // configured output basename (the alignment, plus any
        // `<basename>.html` / `<basename>.score_ascii` etc. T-Coffee
        // emits adjacent), `*.dnd` guide trees, and `*.log` files.
        // Provenance hashes the case.toml since the basename lives in
        // [bio.tcoffee].output_basename and we can recover it from the
        // command without re-parsing.
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "T-Coffee",
            "unknown",
            &job.workdir.join("case.toml"),
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Recover the basename from the prepared command. `prepare()`
        // pushed `-outfile=<basename>.aln` at index 3; strip the prefix
        // and the trailing `.aln` to get the user's basename.
        let basename = job
            .native_command
            .get(3)
            .and_then(|s| s.to_str())
            .and_then(|s| s.strip_prefix("-outfile="))
            .and_then(|s| Path::new(s).file_stem().and_then(|stem| stem.to_str()))
            .map(str::to_string)
            .unwrap_or_default();

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-tcoffee", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());

            let starts_with_basename = matches!(
                (&file_name, basename.as_str()),
                (Some(name), bn) if !bn.is_empty() && name.starts_with(bn)
            );

            let (kind, label) = if starts_with_basename {
                (ArtifactKind::Tabular, "T-Coffee alignment".to_string())
            } else if ext.as_deref() == Some("dnd") {
                (ArtifactKind::Native, "T-Coffee guide tree".to_string())
            } else if ext.as_deref() == Some("log") {
                (ArtifactKind::Log, "T-Coffee log".to_string())
            } else {
                continue;
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
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry to
        // surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.tcoffee.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = TCoffeeAdapter::new().info();
        assert_eq!(info.id, "tcoffee");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "T-Coffee");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = TCoffeeAdapter::new().info();
        // 13.0.0 is the floor we test against; 14.0 reserves room for
        // an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(13, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(14, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = TCoffeeAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.tcoffee.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = TCoffeeAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}

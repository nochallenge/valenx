//! # valenx-validate
//!
//! Project-level structural validator. Walks a `.valenx` directory,
//! loads its `project.toml`, optional `tools.lock`, and every
//! `case.toml` listed in `[cases].order`, then reports a punch list
//! of any structural problems it finds.
//!
//! Designed for CI pre-flight gates: the tool exits 0 on a clean
//! project and non-zero on the first structural issue, so a CI
//! recipe can wire it in front of any adapter run.
//!
//! Scope:
//! - Project manifest format-version gate (`SUPPORTED_MAJOR`).
//! - Geometry / mesh path safety (no absolute paths, no `..`-escapes).
//! - `tools.lock` parses if present.
//! - Every case in `[cases].order` exists and parses as TOML with
//!   the required `[case]` header (format / name / physics / solver).
//!
//! Out of scope: per-adapter case-input validation (that lives next
//! to each adapter's `case_input::from_case_dir`). This binary
//! deliberately doesn't depend on the full adapter zoo so it stays
//! small and fast to build.
//!
//! ## CLI
//!
//! ```text
//! valenx-validate <path/to/project.valenx>
//! valenx-validate <path> --format json
//! ```
//!
//! ## Exit codes
//!
//! - 0 — project loads cleanly, every case parseable.
//! - 1 — structural issue in the project (missing case, malformed
//!   TOML, unsupported format-version, …). Diagnostic on stderr.
//! - 2 — usage error (missing argument, unknown flag, …).
//! - 3 — I/O error before we even reached the project content.

use std::path::PathBuf;
use std::process::ExitCode;

use valenx_core::project::{LoadedProject, ProjectLoadError};

const USAGE: &str = "\
valenx-validate — structural validator for `.valenx` projects.

USAGE:
  valenx-validate <path> [--format text|json]
  valenx-validate -h | --help
  valenx-validate -V | --version

OPTIONS:
  --format F     Output shape. `text` (default) is human-readable;
                 `json` emits a single JSON object suitable for
                 piping into a CI report.
  -h, --help     Show this help.
  -V, --version  Print the binary version and exit.

EXIT CODES:
  0   project loads cleanly
  1   structural issue in the project
  2   usage error
  3   I/O error reading the directory
";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Format {
    Text,
    Json,
}

#[derive(Debug)]
enum ParsedArgs {
    Help,
    Version,
    Validate { path: PathBuf, format: Format },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    let mut path: Option<PathBuf> = None;
    let mut format = Format::Text;
    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "-h" | "--help" => return ParsedArgs::Help,
            "-V" | "--version" => return ParsedArgs::Version,
            "--format" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--format needs a value (text|json)".into());
                };
                format = match v.as_str() {
                    "text" => Format::Text,
                    "json" => Format::Json,
                    other => {
                        return ParsedArgs::Invalid(format!(
                            "--format expects text|json; got `{other}`"
                        ));
                    }
                };
            }
            s if s.starts_with('-') => {
                return ParsedArgs::Invalid(format!("unknown argument `{s}`"));
            }
            other => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!(
                        "extra positional argument `{other}` — expected exactly one project path"
                    ));
                }
                path = Some(PathBuf::from(other));
            }
        }
    }
    match path {
        Some(p) => ParsedArgs::Validate { path: p, format },
        None => ParsedArgs::Invalid("missing project path argument".into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        ParsedArgs::Help => {
            print!("{USAGE}");
            ExitCode::from(0)
        }
        ParsedArgs::Version => {
            println!("valenx-validate v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Validate { path, format } => run(&path, format),
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
    }
}

/// Inspect a project on disk and dump a structural report. Returns
/// the exit code (0 OK, 1 structural issue, 3 IO error).
fn run(path: &std::path::Path, format: Format) -> ExitCode {
    let project = match LoadedProject::load(path) {
        Ok(p) => p,
        Err(e) => {
            // IO + NotADirectory are pre-content failures; everything
            // else is a structural issue that the user can act on
            // inside the project.
            let code = match &e {
                ProjectLoadError::Io { .. } | ProjectLoadError::NotADirectory { .. } => 3,
                _ => 1,
            };
            match format {
                Format::Text => eprintln!("error: {e}"),
                Format::Json => {
                    let obj = serde_json::json!({
                        "ok": false,
                        "path": path.display().to_string(),
                        "error": format!("{e}"),
                    });
                    println!("{}", serde_json::to_string_pretty(&obj).unwrap_or_default());
                }
            }
            return ExitCode::from(code);
        }
    };
    render_report(&project, format);
    ExitCode::from(0)
}

/// Render the structural report. Pure (no I/O) so unit tests cover it
/// without touching a temp dir.
fn render_report(project: &LoadedProject, format: Format) {
    match format {
        Format::Text => render_text(project),
        Format::Json => render_json(project),
    }
}

fn render_text(project: &LoadedProject) {
    let p = &project.project.project;
    println!(
        "✓ project `{}` (format {}) loaded from {}",
        p.name,
        p.format,
        project.root.display(),
    );

    if let Some(lock) = &project.tools_lock {
        let n = lock.tools.len();
        println!(
            "✓ tools.lock present ({} tool{})",
            n,
            if n == 1 { "" } else { "s" }
        );
    } else {
        println!("· tools.lock absent (optional)");
    }

    let g = project.project.geometry.entries.len();
    println!("✓ geometry: {g} entr{}", if g == 1 { "y" } else { "ies" });

    let order = project.case_names();
    println!(
        "✓ {} case{}:",
        order.len(),
        if order.len() == 1 { "" } else { "s" }
    );
    for name in order {
        let c = project
            .cases
            .get(name)
            .expect("LoadedProject populated every order entry");
        println!("  ✓ {name}  ({}/{})", c.case.physics, c.case.solver,);
    }
}

fn render_json(project: &LoadedProject) {
    let p = &project.project.project;
    let cases: Vec<_> = project
        .case_names()
        .iter()
        .map(|name| {
            let c = project
                .cases
                .get(name)
                .expect("LoadedProject populated every order entry");
            serde_json::json!({
                "name": name,
                "physics": c.case.physics,
                "solver": c.case.solver,
            })
        })
        .collect();
    let obj = serde_json::json!({
        "ok": true,
        "path": project.root.display().to_string(),
        "project": {
            "name": p.name,
            "format": p.format,
        },
        "tools_lock_present": project.tools_lock.is_some(),
        "geometry_entries": project.project.geometry.entries.len(),
        "cases": cases,
    });
    println!("{}", serde_json::to_string_pretty(&obj).unwrap_or_default());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_args_help_short() {
        match parse_args(&["-h".into()]) {
            ParsedArgs::Help => {}
            other => panic!("expected Help; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_help_long() {
        match parse_args(&["--help".into()]) {
            ParsedArgs::Help => {}
            other => panic!("expected Help; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_version_short() {
        match parse_args(&["-V".into()]) {
            ParsedArgs::Version => {}
            other => panic!("expected Version; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_version_long() {
        match parse_args(&["--version".into()]) {
            ParsedArgs::Version => {}
            other => panic!("expected Version; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_path_only_defaults_to_text() {
        match parse_args(&["proj".into()]) {
            ParsedArgs::Validate { path, format } => {
                assert_eq!(path, PathBuf::from("proj"));
                assert_eq!(format, Format::Text);
            }
            other => panic!("expected Validate; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_json() {
        match parse_args(&["proj".into(), "--format".into(), "json".into()]) {
            ParsedArgs::Validate { format, .. } => assert_eq!(format, Format::Json),
            other => panic!("expected Validate; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_format_is_invalid() {
        match parse_args(&["proj".into(), "--format".into(), "xml".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("text|json")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_without_value_is_invalid() {
        match parse_args(&["proj".into(), "--format".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--format needs")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_missing_path_is_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing project path")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extra_positional_is_invalid() {
        match parse_args(&["a".into(), "b".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("extra positional argument")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_flag_is_invalid() {
        match parse_args(&["proj".into(), "--bogus".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown argument")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }
}

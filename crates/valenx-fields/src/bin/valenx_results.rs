//! # valenx-results
//!
//! Headless inspector for the `results.json` sidecar that
//! `ValenxApp` writes next to every finished run's workdir.
//!
//! Reads the JSON, deserialises into [`valenx_fields::Results`], and
//! prints a punch list of fields, scalars, artifacts, and
//! provenance. Useful for post-run debugging without firing up the
//! GUI — `cargo run --bin valenx-results -- /tmp/valenx-run-…/results.json`.
//!
//! Two output modes:
//! - `text` (default) — human-readable summary.
//! - `json` — pretty-printed JSON envelope of the same fields,
//!   suitable for piping into a CI report or downstream tool.
//!
//! ## Exit codes
//!
//! - 0 — file loaded and printed.
//! - 1 — file present but not a valid `Results` JSON.
//! - 2 — usage error (missing argument, unknown flag, bad
//!   `--format`, …).
//! - 3 — I/O error reading the file.

use std::io::Read;
use std::path::PathBuf;
use std::process::ExitCode;

use valenx_fields::Results;

/// Round-23 named finding: per-file / per-stdin cap on the
/// `valenx-results` CLI inspector. Pre-fix both branches did bare
/// `std::io::stdin().read_to_string(...)` / `std::fs::read_to_string(...)`
/// against caller-supplied input. A `cat /dev/zero | valenx-results -`
/// or stale path pointing at a multi-GB file would OOM the process
/// before serde_json saw the first byte. 256 MiB matches the bio
/// CLI cap in `valenx_core::io_caps::MAX_BIO_CLI_BYTES` (which this
/// crate cannot depend on — valenx-core depends on valenx-fields,
/// so we duplicate the cap inline here).
const MAX_RESULTS_INPUT_BYTES: u64 = 256 * 1024 * 1024;

const USAGE: &str = "\
valenx-results — headless inspector for adapter `results.json` files.

USAGE:
  valenx-results <path/to/results.json> [--format text|json]
  valenx-results -                        # read from stdin
  valenx-results -h | --help
  valenx-results -V | --version

OPTIONS:
  --format F     Output shape. `text` (default) is human-readable;
                 `json` re-emits the file pretty-printed.
  -h, --help     Show this help.
  -V, --version  Print the binary version and exit.

EXIT CODES:
  0   results loaded and printed
  1   file present but not a valid Results JSON
  2   usage error
  3   I/O error reading the file
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
    Inspect { path: PathBuf, format: Format },
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
            s if s.starts_with('-') && s != "-" => {
                return ParsedArgs::Invalid(format!("unknown argument `{s}`"));
            }
            other => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!(
                        "extra positional argument `{other}` — expected exactly one results.json path"
                    ));
                }
                path = Some(PathBuf::from(other));
            }
        }
    }
    match path {
        Some(p) => ParsedArgs::Inspect { path: p, format },
        None => ParsedArgs::Invalid("missing path argument".into()),
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
            println!("valenx-results v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Inspect { path, format } => run(&path, format),
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
    }
}

fn run(path: &std::path::Path, format: Format) -> ExitCode {
    // Treat `-` as a directive to read from stdin, matching standard
    // Unix CLI convention. Lets users pipe a results.json from
    // anywhere (`kubectl get …`, `aws s3 cp s3://… -`, in-flight
    // output of a sweep harness) without staging to a temp file.
    let (text, display_path): (String, String) = if path == std::path::Path::new("-") {
        // Round-23: bounded stdin read — `cat /dev/zero | valenx-results -`
        // pre-fix would slurp unbounded. Take cap+1 and reject on overflow.
        let mut buf = Vec::new();
        if let Err(e) = std::io::stdin()
            .lock()
            .take(MAX_RESULTS_INPUT_BYTES.saturating_add(1))
            .read_to_end(&mut buf)
        {
            eprintln!("error: read stdin: {e}");
            return ExitCode::from(3);
        }
        if buf.len() as u64 > MAX_RESULTS_INPUT_BYTES {
            eprintln!("error: stdin exceeded {MAX_RESULTS_INPUT_BYTES}-byte cap");
            return ExitCode::from(3);
        }
        match String::from_utf8(buf) {
            Ok(s) => (s, "<stdin>".to_string()),
            Err(e) => {
                eprintln!("error: stdin not valid UTF-8: {e}");
                return ExitCode::from(3);
            }
        }
    } else {
        // Round-23: bounded file read — stat then bounded-take so a
        // file growing past the cap mid-read still produces a clean
        // error rather than slurping unbounded.
        match std::fs::metadata(path) {
            Ok(md) if md.len() > MAX_RESULTS_INPUT_BYTES => {
                eprintln!(
                    "error: {} exceeds {MAX_RESULTS_INPUT_BYTES}-byte cap (actual: {})",
                    path.display(),
                    md.len()
                );
                return ExitCode::from(3);
            }
            Ok(_) => {}
            Err(e) => {
                eprintln!("error: stat {}: {e}", path.display());
                return ExitCode::from(3);
            }
        }
        let mut buf = Vec::new();
        match std::fs::File::open(path).and_then(|f| {
            f.take(MAX_RESULTS_INPUT_BYTES.saturating_add(1))
                .read_to_end(&mut buf)
                .map(|_| ())
        }) {
            Ok(()) => {}
            Err(e) => {
                eprintln!("error: read {}: {e}", path.display());
                return ExitCode::from(3);
            }
        }
        if buf.len() as u64 > MAX_RESULTS_INPUT_BYTES {
            eprintln!(
                "error: {} grew past {MAX_RESULTS_INPUT_BYTES}-byte cap mid-read",
                path.display()
            );
            return ExitCode::from(3);
        }
        match String::from_utf8(buf) {
            Ok(s) => (s, path.display().to_string()),
            Err(e) => {
                eprintln!("error: {} not valid UTF-8: {e}", path.display());
                return ExitCode::from(3);
            }
        }
    };
    let results: Results = match serde_json::from_str(&text) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: parse {display_path}: {e}");
            return ExitCode::from(1);
        }
    };
    match format {
        Format::Text => render_text(&results, std::path::Path::new(&display_path)),
        Format::Json => render_json(&results),
    }
    ExitCode::from(0)
}

fn render_text(results: &Results, path: &std::path::Path) {
    println!("results from {}", path.display());
    println!(
        "  case: {} (run_id={})",
        results.meta.case_id, results.provenance.run_id
    );
    println!(
        "  adapter: {} v{}  (tool {} v{})",
        results.provenance.adapter,
        results.provenance.adapter_version,
        results.provenance.tool,
        results.provenance.tool_version,
    );

    let n_fields = results.fields.len();
    println!("  {n_fields} field{}", if n_fields == 1 { "" } else { "s" });
    let mut field_names: Vec<&str> = results.fields.names().collect();
    field_names.sort();
    field_names.dedup();
    for name in field_names {
        let timesteps = results.fields.time_series(name);
        let n_times = timesteps.len();
        println!(
            "    · {name}{}",
            if n_times > 1 {
                format!(" ({n_times} timesteps)")
            } else {
                String::new()
            }
        );
    }

    let n_scalars = results.scalars.len();
    println!(
        "  {n_scalars} scalar{}",
        if n_scalars == 1 { "" } else { "s" }
    );
    let mut scalar_names: Vec<&str> = results.scalars.names().collect();
    scalar_names.sort();
    scalar_names.dedup();
    // Cap the per-name list so a 10k-row PyBaMM time series doesn't
    // blast the terminal. Anyone wanting the full firehose drops to
    // `--format json`.
    const MAX_LISTED_SCALARS: usize = 12;
    for name in scalar_names.iter().take(MAX_LISTED_SCALARS) {
        let n = results.scalars.all(name).len();
        if let Some(rec) = results.scalars.get(name) {
            println!(
                "    · {name} = {} ({n} sample{})",
                rec.value,
                if n == 1 { "" } else { "s" }
            );
        }
    }
    if scalar_names.len() > MAX_LISTED_SCALARS {
        println!(
            "    … {} more scalar names hidden — `--format json` for the full list",
            scalar_names.len() - MAX_LISTED_SCALARS
        );
    }

    let n_artifacts = results.artifacts.len();
    println!(
        "  {n_artifacts} artifact{}",
        if n_artifacts == 1 { "" } else { "s" }
    );
    for a in &results.artifacts {
        println!(
            "    · {} ({}, {})",
            a.path.display(),
            a.label,
            artifact_kind_label(&a.kind),
        );
    }
}

fn artifact_kind_label(kind: &valenx_fields::artifact::ArtifactKind) -> &'static str {
    use valenx_fields::artifact::ArtifactKind::*;
    match kind {
        Native => "native",
        Tabular => "tabular",
        VizData => "viz",
        Log => "log",
        Image => "image",
        Other => "other",
    }
}

fn render_json(results: &Results) {
    // Re-emit the file pretty-printed; the Results type is
    // serde-roundtripping per its existing tests, so this is just a
    // nice CLI surface on top of `serde_json::to_string_pretty`.
    let s = serde_json::to_string_pretty(results)
        .unwrap_or_else(|e| format!("{{\"error\":\"serialise: {e}\"}}"));
    println!("{s}");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_help_short_long() {
        match parse_args(&["-h".into()]) {
            ParsedArgs::Help => {}
            other => panic!("expected Help; got {other:?}"),
        }
        match parse_args(&["--help".into()]) {
            ParsedArgs::Help => {}
            other => panic!("expected Help; got {other:?}"),
        }
    }

    #[test]
    fn parse_version_short_long() {
        match parse_args(&["-V".into()]) {
            ParsedArgs::Version => {}
            other => panic!("expected Version; got {other:?}"),
        }
        match parse_args(&["--version".into()]) {
            ParsedArgs::Version => {}
            other => panic!("expected Version; got {other:?}"),
        }
    }

    #[test]
    fn parse_path_only_defaults_text() {
        match parse_args(&["x.json".into()]) {
            ParsedArgs::Inspect { path, format } => {
                assert_eq!(path, PathBuf::from("x.json"));
                assert_eq!(format, Format::Text);
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_dash_is_stdin_directive() {
        // `-` as the path means read from stdin per Unix convention.
        // The arg parser must accept it and not flag it as
        // "unknown argument" via the dash-prefix branch.
        match parse_args(&["-".into()]) {
            ParsedArgs::Inspect { path, .. } => {
                assert_eq!(path, PathBuf::from("-"));
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_format_json() {
        match parse_args(&["x.json".into(), "--format".into(), "json".into()]) {
            ParsedArgs::Inspect { format, .. } => assert_eq!(format, Format::Json),
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_unknown_format_invalid() {
        match parse_args(&["x.json".into(), "--format".into(), "xml".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("text|json")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_format_without_value_invalid() {
        match parse_args(&["x.json".into(), "--format".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--format needs")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_missing_path_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing path")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_extra_positional_invalid() {
        match parse_args(&["a.json".into(), "b.json".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("extra positional")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_unknown_flag_invalid() {
        match parse_args(&["x.json".into(), "--bogus".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown argument")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn artifact_kind_label_round_trip() {
        use valenx_fields::artifact::ArtifactKind::*;
        assert_eq!(artifact_kind_label(&Native), "native");
        assert_eq!(artifact_kind_label(&Tabular), "tabular");
        assert_eq!(artifact_kind_label(&VizData), "viz");
        assert_eq!(artifact_kind_label(&Log), "log");
        assert_eq!(artifact_kind_label(&Image), "image");
        assert_eq!(artifact_kind_label(&Other), "other");
    }
}

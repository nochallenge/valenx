//! `valenx-mesh-info` — headless CLI for mesh quality inspection.
//!
//! Reads a canonical `valenx_mesh::Mesh` from a JSON file and
//! prints the `QualityReport` (via its `Display` impl) plus the
//! aspect-ratio + skewness histograms. Equivalent of the GUI's
//! Quality panel for users who don't have the GUI to hand —
//! CI mesh-quality gates, SSH'd-into-a-cluster debugging, etc.
//!
//! Usage:
//!
//! ```text
//! valenx-mesh-info <mesh.json>
//! valenx-mesh-info help
//! ```
//!
//! Exit codes:
//! - 0: report printed
//! - 1: I/O failure (file missing / unreadable)
//! - 2: invalid CLI usage (missing arg, unknown flag)
//! - 3: parse failure (file isn't a Mesh in canonical JSON form)

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use valenx_mesh::{
    aspect_ratio_histogram, skewness_histogram, AspectRatioHistogram, Mesh, SkewnessHistogram,
    DEFAULT_AR_BUCKETS, DEFAULT_SKEW_BUCKETS,
};

const USAGE: &str = "\
valenx-mesh-info — print quality report + histograms for a mesh.

USAGE:
  valenx-mesh-info <mesh.json> [--format text|json]
                                [--check METRIC=THRESHOLD ...]
  valenx-mesh-info -                       # read mesh JSON from stdin
  valenx-mesh-info help
  valenx-mesh-info -V | --version

OPTIONS:
  --format F   Output format: `text` (default, human-readable
               with Unicode-bar histograms) or `json` (single
               serde object — quality + both histograms with
               their bucket edges + counts).

  --check K=V  Quality gate. Repeatable. Exit 4 if any check
               fails. Supported metrics:
                 max-skew=<float>           upper bound on max_skewness
                 max-aspect=<float>         upper bound on max_aspect_ratio
                 inverted=<int>             upper bound on inverted_count
                 min-orthogonality=<float>  lower bound on min_orthogonality
               A check that wants a value the report doesn't have
               (e.g. `min-orthogonality` on a single-element mesh)
               counts as a failure — fail-loud rather than skip.

EXIT CODES:
  0   report printed (and all --check passed, if any)
  1   I/O failure
  2   invalid CLI usage
  3   parse failure (not a valid Mesh JSON)
  4   one or more --check thresholds failed
";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OutputFormat {
    Text,
    Json,
}

impl OutputFormat {
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "text" | "txt" => Some(Self::Text),
            "json" => Some(Self::Json),
            _ => None,
        }
    }
}

/// One CI-gate predicate parsed from a `--check K=V` flag.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Check {
    /// `max_skewness <= v`
    MaxSkewLeq(f64),
    /// `max_aspect_ratio <= v`
    MaxAspectLeq(f64),
    /// `inverted_count <= v`
    InvertedLeq(u64),
    /// `min_orthogonality >= v`
    MinOrthogonalityGeq(f64),
}

impl Check {
    /// Parse a `metric=value` string. Returns the parsed Check or
    /// an error string suitable for ParsedArgs::Invalid.
    fn parse(spec: &str) -> Result<Self, String> {
        let Some((metric, value)) = spec.split_once('=') else {
            return Err(format!("--check `{spec}` needs a `metric=value` form"));
        };
        match metric {
            "max-skew" => value
                .parse::<f64>()
                .map(Check::MaxSkewLeq)
                .map_err(|_| format!("--check max-skew=`{value}` (not a number)")),
            "max-aspect" => value
                .parse::<f64>()
                .map(Check::MaxAspectLeq)
                .map_err(|_| format!("--check max-aspect=`{value}` (not a number)")),
            "inverted" => value
                .parse::<u64>()
                .map(Check::InvertedLeq)
                .map_err(|_| format!("--check inverted=`{value}` (not an integer)")),
            "min-orthogonality" => value
                .parse::<f64>()
                .map(Check::MinOrthogonalityGeq)
                .map_err(|_| format!("--check min-orthogonality=`{value}` (not a number)")),
            other => Err(format!(
                "--check unknown metric `{other}` (try max-skew / max-aspect / inverted / min-orthogonality)"
            )),
        }
    }
}

#[derive(Debug, PartialEq)]
enum ParsedArgs {
    Help,
    Version,
    Inspect {
        path: PathBuf,
        format: OutputFormat,
        checks: Vec<Check>,
    },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    if args.is_empty() {
        return ParsedArgs::Invalid("missing mesh file path".into());
    }
    if matches!(args[0].as_str(), "help" | "-h" | "--help") {
        return ParsedArgs::Help;
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        return ParsedArgs::Version;
    }
    let mut path: Option<PathBuf> = None;
    let mut format = OutputFormat::Text;
    let mut checks: Vec<Check> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--format" | "-f" => {
                i += 1;
                let Some(v) = args.get(i) else {
                    return ParsedArgs::Invalid("--format needs a value (text|json)".into());
                };
                let Some(parsed) = OutputFormat::from_str(v) else {
                    return ParsedArgs::Invalid(format!("unknown format `{v}` (try text or json)"));
                };
                format = parsed;
            }
            "--check" | "-c" => {
                i += 1;
                let Some(spec) = args.get(i) else {
                    return ParsedArgs::Invalid("--check needs a `metric=value` argument".into());
                };
                match Check::parse(spec) {
                    Ok(c) => checks.push(c),
                    Err(msg) => return ParsedArgs::Invalid(msg),
                }
            }
            other if other.starts_with('-') && other != "-" => {
                return ParsedArgs::Invalid(format!("unknown flag `{other}`"));
            }
            _ => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!("unexpected extra argument `{}`", args[i]));
                }
                path = Some(PathBuf::from(&args[i]));
            }
        }
        i += 1;
    }
    let Some(path) = path else {
        return ParsedArgs::Invalid("missing mesh file path".into());
    };
    ParsedArgs::Inspect {
        path,
        format,
        checks,
    }
}

/// Run every Check against the report and return a list of failure
/// messages (empty if everything passed). The format of each message
/// is stable enough for grep/awk pipelines: `<metric>: actual=<X>
/// threshold=<Y>`.
fn evaluate_checks(report: &valenx_mesh::QualityReport, checks: &[Check]) -> Vec<String> {
    let mut failures = Vec::new();
    for check in checks {
        match *check {
            Check::MaxSkewLeq(threshold) => {
                let Some(actual) = report.max_skewness else {
                    failures.push(format!(
                        "max-skew: not measurable (no elements with computable skewness); threshold={threshold}"
                    ));
                    continue;
                };
                if actual > threshold {
                    failures.push(format!(
                        "max-skew: actual={actual:.6} threshold={threshold}"
                    ));
                }
            }
            Check::MaxAspectLeq(threshold) => {
                let Some(actual) = report.max_aspect_ratio else {
                    failures.push(format!("max-aspect: not measurable; threshold={threshold}"));
                    continue;
                };
                if actual > threshold {
                    failures.push(format!(
                        "max-aspect: actual={actual:.6} threshold={threshold}"
                    ));
                }
            }
            Check::InvertedLeq(threshold) => {
                if report.inverted_count > threshold {
                    failures.push(format!(
                        "inverted: actual={} threshold={threshold}",
                        report.inverted_count
                    ));
                }
            }
            Check::MinOrthogonalityGeq(threshold) => {
                let Some(actual) = report.min_orthogonality else {
                    failures.push(format!(
                        "min-orthogonality: not measurable (no interior faces); threshold={threshold}"
                    ));
                    continue;
                };
                if actual < threshold {
                    failures.push(format!(
                        "min-orthogonality: actual={actual:.6} threshold={threshold}"
                    ));
                }
            }
        }
    }
    failures
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        ParsedArgs::Help => {
            print!("{USAGE}");
            ExitCode::from(0)
        }
        ParsedArgs::Version => {
            println!("valenx-mesh-info v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
        ParsedArgs::Inspect {
            path,
            format,
            checks,
        } => match inspect(&path, format, &checks) {
            Ok(()) => ExitCode::from(0),
            Err(e) => {
                eprintln!("error: {e}");
                e.exit_code()
            }
        },
    }
}

#[derive(Debug)]
enum InspectError {
    Io(String),
    Parse(String),
    /// One or more --check thresholds failed. Each string is one
    /// failure line (already formatted by `evaluate_checks`).
    CheckFailures(Vec<String>),
}

impl InspectError {
    fn exit_code(&self) -> ExitCode {
        match self {
            InspectError::Io(_) => ExitCode::from(1),
            InspectError::Parse(_) => ExitCode::from(3),
            InspectError::CheckFailures(_) => ExitCode::from(4),
        }
    }
}

impl std::fmt::Display for InspectError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InspectError::Io(s) | InspectError::Parse(s) => f.write_str(s),
            InspectError::CheckFailures(lines) => {
                writeln!(f, "{} quality check(s) failed:", lines.len())?;
                for l in lines {
                    writeln!(f, "  - {l}")?;
                }
                Ok(())
            }
        }
    }
}

/// Cap on the bytes the inspector CLI will read from stdin OR from
/// the file path argument before refusing to allocate further.
///
/// Round-20 L2: pre-fix both the stdin and file-path branches were
/// unbounded — `read_to_end` would grow the `Vec` to whatever the
/// input was, so a user piping `/dev/zero` (or pointing the CLI at a
/// multi-GB stale mesh file from a runaway adapter) could OOM the
/// inspector before `serde_json::from_slice` saw the first token.
/// 256 MiB is generous for legitimate canonical mesh JSON (even a
/// 10M-element tet mesh fits inside that envelope with serde_json's
/// representation overhead) while refusing the slurp-the-pipe DoS.
const MAX_CLI_MESH_BYTES: u64 = 256 * 1024 * 1024;

fn inspect(path: &Path, format: OutputFormat, checks: &[Check]) -> Result<(), InspectError> {
    use std::io::Read;
    // `-` reads the mesh JSON from stdin per Unix convention. Lets
    // users pipe a canonical mesh straight from gmsh / netgen
    // collect() output without staging to a temp file.
    let (bytes, display_path): (Vec<u8>, String) = if path == Path::new("-") {
        // Round-20 L2: bounded stdin read. `take(cap + 1)` so a file
        // that exceeds the cap by even one byte trips the check
        // below rather than silently truncating to `cap` bytes (which
        // would produce a malformed-JSON parse error masking the real
        // problem).
        let mut buf = Vec::new();
        std::io::stdin()
            .take(MAX_CLI_MESH_BYTES + 1)
            .read_to_end(&mut buf)
            .map_err(|e| InspectError::Io(format!("read stdin: {e}")))?;
        if buf.len() as u64 > MAX_CLI_MESH_BYTES {
            return Err(InspectError::Io(format!(
                "read stdin: input exceeds {MAX_CLI_MESH_BYTES}-byte cap",
            )));
        }
        (buf, "<stdin>".to_string())
    } else {
        // Round-20 L2: stat-then-bounded-read on the file path.
        let meta = std::fs::metadata(path)
            .map_err(|e| InspectError::Io(format!("stat {}: {e}", path.display())))?;
        if meta.len() > MAX_CLI_MESH_BYTES {
            return Err(InspectError::Io(format!(
                "read {}: file is {} bytes, exceeds {MAX_CLI_MESH_BYTES}-byte cap",
                path.display(),
                meta.len(),
            )));
        }
        let mut buf = Vec::new();
        std::fs::File::open(path)
            .map_err(|e| InspectError::Io(format!("open {}: {e}", path.display())))?
            .take(MAX_CLI_MESH_BYTES + 1)
            .read_to_end(&mut buf)
            .map_err(|e| InspectError::Io(format!("read {}: {e}", path.display())))?;
        if buf.len() as u64 > MAX_CLI_MESH_BYTES {
            return Err(InspectError::Io(format!(
                "read {}: file grew past {MAX_CLI_MESH_BYTES}-byte cap mid-read",
                path.display(),
            )));
        }
        (buf, path.display().to_string())
    };
    let mut mesh: Mesh = serde_json::from_slice(&bytes)
        .map_err(|e| InspectError::Parse(format!("parse {display_path} as Mesh JSON: {e}",)))?;
    // recompute_quality_stats populates mesh.stats AND returns the
    // report — single walk over elements, no double-computation.
    let report = mesh.recompute_quality_stats();
    let ar_hist = aspect_ratio_histogram(&mesh, DEFAULT_AR_BUCKETS);
    let sk_hist = skewness_histogram(&mesh, DEFAULT_SKEW_BUCKETS);

    match format {
        OutputFormat::Text => {
            print_report(&display_path, &report, &ar_hist, &sk_hist);
        }
        OutputFormat::Json => {
            println!(
                "{}",
                render_json(&display_path, &report, &ar_hist, &sk_hist)
            );
        }
    }
    // Quality gates run after the report so users always see the
    // numbers, even when a check fails. The error path is exit
    // code 4 with the failure list on stderr.
    let failures = evaluate_checks(&report, checks);
    if !failures.is_empty() {
        return Err(InspectError::CheckFailures(failures));
    }
    Ok(())
}

/// Build the JSON output as a string. Uses serde_json's pretty
/// printer so the output reads in a terminal but is also a single
/// parseable object for CI gates / `jq` pipelines.
fn render_json(
    path: &str,
    report: &valenx_mesh::QualityReport,
    ar: &AspectRatioHistogram,
    sk: &SkewnessHistogram,
) -> String {
    let v = serde_json::json!({
        "path": path,
        "quality": report,
        "aspect_histogram": ar,
        "skewness_histogram": sk,
    });
    serde_json::to_string_pretty(&v).expect("static structure serialises")
}

/// Print the multi-section text report. Public-via-`pub(crate)` so
/// integration tests in `tests/` can drive it without spawning the
/// binary; not exported beyond this binary because the format is
/// stable for grep/awk but not part of a real public API.
fn print_report(
    path: &str,
    report: &valenx_mesh::QualityReport,
    ar: &AspectRatioHistogram,
    sk: &SkewnessHistogram,
) {
    println!("# Mesh quality report — {path}");
    print!("{report}");
    println!();
    println!("# Aspect-ratio histogram");
    print_ar_histogram(ar);
    println!();
    println!("# Skewness histogram");
    print_sk_histogram(sk);
}

fn print_ar_histogram(hist: &AspectRatioHistogram) {
    let max_count = hist
        .counts
        .iter()
        .copied()
        .chain(std::iter::once(hist.overflow))
        .max()
        .unwrap_or(0);
    for (i, &edge) in hist.buckets.iter().enumerate() {
        let count = hist.counts[i];
        println!("  ≤ {:>6.2}: {} {count}", edge, bar(count, max_count, 24));
    }
    if hist.overflow > 0 {
        println!(
            "  >  max: {} {n}",
            bar(hist.overflow, max_count, 24),
            n = hist.overflow
        );
    }
    if hist.uncategorised > 0 {
        println!("  uncategorised: {}", hist.uncategorised);
    }
}

fn print_sk_histogram(hist: &SkewnessHistogram) {
    let labels = ["excellent", "good", "acceptable", "poor", "very poor"];
    let max_count = hist.counts.iter().copied().max().unwrap_or(0);
    for (i, &edge) in hist.buckets.iter().enumerate() {
        let count = hist.counts[i];
        let label = labels.get(i).copied().unwrap_or("");
        println!(
            "  ≤ {edge:.2} ({label:<11}): {} {count}",
            bar(count, max_count, 24)
        );
    }
    if hist.uncategorised > 0 {
        println!("  uncategorised: {}", hist.uncategorised);
    }
}

/// Unicode-block bar of width up to `width` proportional to `count/max`.
/// Always emits at least one block when count > 0 so non-empty buckets
/// stay visible.
fn bar(count: u64, max: u64, width: usize) -> String {
    if max == 0 || count == 0 {
        return " ".repeat(width);
    }
    let n = ((count as f64 / max as f64) * width as f64).ceil() as usize;
    let n = n.max(1).min(width);
    let mut s = String::with_capacity(width * 3);
    for _ in 0..n {
        s.push('█');
    }
    for _ in n..width {
        s.push(' ');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_args_version_variants() {
        assert_eq!(parse_args(&["-V".into()]), ParsedArgs::Version);
        assert_eq!(parse_args(&["--version".into()]), ParsedArgs::Version);
    }

    #[test]
    fn parse_args_help_variants() {
        let cases = [
            vec!["help".into()],
            vec!["-h".into()],
            vec!["--help".into()],
        ];
        for args in &cases {
            assert_eq!(parse_args(args), ParsedArgs::Help, "for args {args:?}");
        }
    }

    #[test]
    fn parse_args_missing_path_is_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("missing mesh file")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_picks_up_path_with_default_format() {
        match parse_args(&["./mesh.json".into()]) {
            ParsedArgs::Inspect { path, format, .. } => {
                assert_eq!(path, PathBuf::from("./mesh.json"));
                assert_eq!(format, OutputFormat::Text);
            }
            other => panic!("expected Inspect, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extra_arg_is_invalid() {
        match parse_args(&["a.json".into(), "b.json".into()]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("unexpected extra")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_json_flag() {
        match parse_args(&["./mesh.json".into(), "--format".into(), "json".into()]) {
            ParsedArgs::Inspect { path, format, .. } => {
                assert_eq!(path, PathBuf::from("./mesh.json"));
                assert_eq!(format, OutputFormat::Json);
            }
            other => panic!("expected Inspect, got {other:?}"),
        }
        // Order doesn't matter — flag before path also valid.
        match parse_args(&["--format".into(), "json".into(), "./mesh.json".into()]) {
            ParsedArgs::Inspect { path, format, .. } => {
                assert_eq!(path, PathBuf::from("./mesh.json"));
                assert_eq!(format, OutputFormat::Json);
            }
            other => panic!("expected Inspect, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_unknown_value_is_invalid() {
        match parse_args(&["./mesh.json".into(), "--format".into(), "xml".into()]) {
            ParsedArgs::Invalid(s) => {
                assert!(s.contains("xml"), "expected error to mention `xml`: {s}");
            }
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_missing_value_is_invalid() {
        match parse_args(&["--format".into()]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("--format")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_check_flag_collects_multiple_metrics() {
        match parse_args(&[
            "./mesh.json".into(),
            "--check".into(),
            "max-skew=0.9".into(),
            "--check".into(),
            "inverted=0".into(),
            "--check".into(),
            "max-aspect=100".into(),
        ]) {
            ParsedArgs::Inspect { checks, .. } => {
                assert_eq!(checks.len(), 3);
                assert!(matches!(checks[0], Check::MaxSkewLeq(v) if (v - 0.9).abs() < 1e-12));
                assert!(matches!(checks[1], Check::InvertedLeq(0)));
                assert!(matches!(checks[2], Check::MaxAspectLeq(v) if (v - 100.0).abs() < 1e-12));
            }
            other => panic!("expected Inspect, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_check_min_orthogonality_uses_lower_bound() {
        match parse_args(&[
            "./mesh.json".into(),
            "--check".into(),
            "min-orthogonality=0.7".into(),
        ]) {
            ParsedArgs::Inspect { checks, .. } => {
                assert_eq!(checks.len(), 1);
                assert!(
                    matches!(checks[0], Check::MinOrthogonalityGeq(v) if (v - 0.7).abs() < 1e-12)
                );
            }
            other => panic!("expected Inspect, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_check_unknown_metric_is_invalid() {
        match parse_args(&[
            "./mesh.json".into(),
            "--check".into(),
            "fancy-metric=42".into(),
        ]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("fancy-metric")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_check_malformed_value_is_invalid() {
        match parse_args(&[
            "./mesh.json".into(),
            "--check".into(),
            "max-skew=not-a-number".into(),
        ]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("not-a-number")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn parse_args_check_missing_equals_is_invalid() {
        match parse_args(&["./mesh.json".into(), "--check".into(), "max-skew".into()]) {
            ParsedArgs::Invalid(s) => assert!(s.contains("max-skew")),
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn evaluate_checks_passes_within_thresholds() {
        let r = valenx_mesh::QualityReport {
            element_count: 1,
            min_size: Some(0.5),
            max_size: Some(0.5),
            mean_size: Some(0.5),
            max_aspect_ratio: Some(1.4),
            max_skewness: Some(0.25),
            min_orthogonality: Some(0.95),
            inverted_count: 0,
        };
        let checks = vec![
            Check::MaxSkewLeq(0.5),
            Check::MaxAspectLeq(2.0),
            Check::InvertedLeq(0),
            Check::MinOrthogonalityGeq(0.9),
        ];
        let failures = evaluate_checks(&r, &checks);
        assert!(failures.is_empty(), "expected no failures: {failures:?}");
    }

    #[test]
    fn evaluate_checks_reports_each_failure() {
        let r = valenx_mesh::QualityReport {
            element_count: 10,
            min_size: Some(0.0),
            max_size: Some(1.0),
            mean_size: Some(0.5),
            max_aspect_ratio: Some(150.0), // exceeds threshold 100
            max_skewness: Some(0.95),      // exceeds 0.9
            min_orthogonality: Some(0.4),  // below 0.5
            inverted_count: 3,             // exceeds 0
        };
        let checks = vec![
            Check::MaxSkewLeq(0.9),
            Check::MaxAspectLeq(100.0),
            Check::InvertedLeq(0),
            Check::MinOrthogonalityGeq(0.5),
        ];
        let failures = evaluate_checks(&r, &checks);
        assert_eq!(failures.len(), 4, "expected all to fail: {failures:?}");
        // Spot-check the message format includes the metric + actual + threshold.
        assert!(failures[0].contains("max-skew"));
        assert!(failures[0].contains("0.95"));
        assert!(failures[0].contains("0.9"));
    }

    #[test]
    fn evaluate_checks_missing_metric_counts_as_failure() {
        // Mesh has no interior faces -> min_orthogonality is None.
        // A --check min-orthogonality=0.5 should fail because we
        // can't verify the threshold was met.
        let r = valenx_mesh::QualityReport {
            element_count: 1,
            min_size: Some(0.5),
            max_size: Some(0.5),
            mean_size: Some(0.5),
            max_aspect_ratio: None,
            max_skewness: None,
            min_orthogonality: None,
            inverted_count: 0,
        };
        let checks = vec![Check::MinOrthogonalityGeq(0.5)];
        let failures = evaluate_checks(&r, &checks);
        assert_eq!(failures.len(), 1);
        assert!(failures[0].contains("not measurable"));
    }

    #[test]
    fn render_json_round_trips_through_serde() {
        // The JSON output must be a parseable object with the
        // expected top-level keys. Build a tiny report + histograms
        // and verify the rendered string deserialises back.
        let report = valenx_mesh::QualityReport {
            element_count: 3,
            min_size: Some(0.5),
            max_size: Some(1.5),
            mean_size: Some(1.0),
            max_aspect_ratio: Some(1.4),
            max_skewness: Some(0.25),
            min_orthogonality: Some(0.95),
            inverted_count: 0,
        };
        let ar = valenx_mesh::AspectRatioHistogram {
            buckets: vec![1.5, 2.0],
            counts: vec![3, 0],
            overflow: 0,
            uncategorised: 0,
        };
        let sk = valenx_mesh::SkewnessHistogram {
            buckets: vec![0.25, 0.5],
            counts: vec![3, 0],
            uncategorised: 0,
        };
        let out = render_json("test.json", &report, &ar, &sk);
        let v: serde_json::Value = serde_json::from_str(&out).expect("parseable JSON");
        assert_eq!(v["path"], "test.json");
        assert_eq!(v["quality"]["element_count"], 3);
        assert_eq!(v["quality"]["max_skewness"], 0.25);
        assert_eq!(v["aspect_histogram"]["counts"][0], 3);
        assert_eq!(v["skewness_histogram"]["counts"][0], 3);
    }

    #[test]
    fn bar_proportional() {
        // count = max -> full width; count = max/2 -> ~half.
        let s = bar(10, 10, 8);
        assert_eq!(s.chars().filter(|c| *c == '█').count(), 8);
        let h = bar(5, 10, 8);
        let blocks = h.chars().filter(|c| *c == '█').count();
        assert_eq!(blocks, 4);
    }

    #[test]
    fn bar_zero_count_is_all_spaces() {
        let s = bar(0, 100, 4);
        assert_eq!(s, "    ");
    }

    #[test]
    fn bar_small_count_shows_at_least_one_block() {
        let s = bar(1, 1000, 8);
        assert!(s.starts_with('█'));
    }

    #[test]
    fn inspect_reports_on_valid_mesh_json() {
        // Round-trip: serialise a Mesh, write it, run inspect.
        // Inspect prints to stdout; we just verify the function
        // doesn't error on a well-formed mesh.
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};

        let mut m = Mesh::new("right-tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);

        let tmp = std::env::temp_dir().join(format!(
            "valenx-mesh-info-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, serde_json::to_string(&m).unwrap()).unwrap();
        let r = inspect(&tmp, OutputFormat::Text, &[]);
        assert!(r.is_ok(), "expected Ok, got {r:?}");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn inspect_io_error_on_missing_file() {
        let r = inspect(
            Path::new("/tmp/this-does-not-exist-valenx-test.json"),
            OutputFormat::Text,
            &[],
        );
        match r {
            Err(InspectError::Io(_)) => {}
            other => panic!("expected Io error, got {other:?}"),
        }
    }

    #[test]
    fn inspect_parse_error_on_garbage() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mesh-info-bad-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, b"not a mesh").unwrap();
        let r = inspect(&tmp, OutputFormat::Text, &[]);
        match r {
            Err(InspectError::Parse(_)) => {}
            other => panic!("expected Parse error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-20 L2 RED→GREEN: an over-cap mesh JSON file must be
    /// rejected as an IO error with a cap-mention message BEFORE
    /// `serde_json::from_slice` allocates the full file into memory.
    /// Pre-fix the file-path branch did a bare `std::fs::read(path)`
    /// — a 1 GiB stale mesh file from a runaway adapter would have
    /// allocated 1 GiB before any reasonability check ran.
    #[test]
    fn inspect_rejects_oversize_file() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mesh-info-r20l2-oversize-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // 257 MiB — past the 256 MiB MAX_CLI_MESH_BYTES cap. Sparse
        // file via set_len so the test isn't writing 257 MiB.
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.set_len(MAX_CLI_MESH_BYTES + 1).unwrap();
        f.write_all(b"{").unwrap();
        drop(f);
        let r = inspect(&tmp, OutputFormat::Text, &[]);
        match r {
            Err(InspectError::Io(msg)) => {
                assert!(
                    msg.contains("cap") || msg.contains("exceed"),
                    "expected cap-mention in IO error, got: {msg}"
                );
            }
            other => panic!("expected Io error for oversize file, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}

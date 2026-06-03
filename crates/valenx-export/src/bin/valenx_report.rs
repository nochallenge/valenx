//! # valenx-report
//!
//! Headless CSV / HTML exporter. Reads a `results.json` produced by
//! any adapter's `collect()` pipeline and writes whichever artifacts
//! the caller asks for:
//!
//! - `--html <path>` → self-contained HTML report (provenance, fields,
//!   scalars, artifacts) via [`valenx_export::write_html_report`].
//! - `--markdown <path>` → GitHub-flavoured Markdown report via
//!   [`valenx_export::write_markdown_report`]. Same five sections as
//!   the HTML, ideal for PR comments and release notes.
//! - `--csv <path>` → flat scalar history CSV via
//!   [`valenx_export::write_scalars_csv`].
//!
//! At least one of `--html` / `--markdown` / `--csv` must be supplied;
//! the binary refuses to no-op silently. Useful in CI jobs that want
//! to attach a human-readable report to a pull request without
//! booting the GUI.
//!
//! ## CLI
//!
//! ```text
//! valenx-report <path/to/results.json> --html out.html
//! valenx-report <path/to/results.json> --csv scalars.csv
//! valenx-report <path/to/results.json> --html out.html --csv scalars.csv
//! ```
//!
//! ## Exit codes
//!
//! - 0 — every requested artifact was written.
//! - 1 — the input JSON is malformed.
//! - 2 — usage error (missing argument, unknown flag, no destination
//!   given, …).
//! - 3 — I/O error reading the input or writing an output.

use std::path::PathBuf;
use std::process::ExitCode;

use valenx_export::{write_html_report, write_markdown_report, write_scalars_csv};
use valenx_fields::Results;

const USAGE: &str = "\
valenx-report — write CSV / HTML / Markdown reports from a `results.json`.

USAGE:
  valenx-report <path/to/results.json> --html <out.html>
  valenx-report <path/to/results.json> --markdown <out.md>
  valenx-report <path/to/results.json> --csv  <out.csv>
  valenx-report <path/to/results.json> --html out.html --markdown out.md --csv scalars.csv
  valenx-report - --html out.html        # read results.json from stdin
  valenx-report -h | --help
  valenx-report -V | --version

OPTIONS:
  --html      PATH   Write a self-contained HTML report.
  --markdown  PATH   Write a GitHub-flavoured Markdown report.
  --csv       PATH   Write a flat scalar history CSV.
  -h, --help         Show this help.
  -V, --version      Print the binary version and exit.

At least one of --html / --markdown / --csv must be supplied.

EXIT CODES:
  0   every requested artifact written
  1   input JSON is malformed
  2   usage error
  3   I/O error
";

#[derive(Debug)]
enum ParsedArgs {
    Help,
    Version,
    Export {
        input: PathBuf,
        html: Option<PathBuf>,
        markdown: Option<PathBuf>,
        csv: Option<PathBuf>,
    },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    let mut input: Option<PathBuf> = None;
    let mut html: Option<PathBuf> = None;
    let mut markdown: Option<PathBuf> = None;
    let mut csv: Option<PathBuf> = None;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "-h" | "--help" => return ParsedArgs::Help,
            "-V" | "--version" => return ParsedArgs::Version,
            "--html" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--html needs a path".into());
                };
                if html.is_some() {
                    return ParsedArgs::Invalid("--html given twice".into());
                }
                html = Some(PathBuf::from(v));
            }
            "--markdown" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--markdown needs a path".into());
                };
                if markdown.is_some() {
                    return ParsedArgs::Invalid("--markdown given twice".into());
                }
                markdown = Some(PathBuf::from(v));
            }
            "--csv" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--csv needs a path".into());
                };
                if csv.is_some() {
                    return ParsedArgs::Invalid("--csv given twice".into());
                }
                csv = Some(PathBuf::from(v));
            }
            s if s.starts_with('-') && s != "-" => {
                return ParsedArgs::Invalid(format!("unknown argument `{s}`"));
            }
            other => {
                if input.is_some() {
                    return ParsedArgs::Invalid(format!(
                        "extra positional argument `{other}` — expected exactly one input path"
                    ));
                }
                input = Some(PathBuf::from(other));
            }
        }
    }
    let Some(input) = input else {
        return ParsedArgs::Invalid("missing input path argument".into());
    };
    if html.is_none() && markdown.is_none() && csv.is_none() {
        return ParsedArgs::Invalid(
            "at least one of --html / --markdown / --csv must be supplied".into(),
        );
    }
    ParsedArgs::Export {
        input,
        html,
        markdown,
        csv,
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
            println!("valenx-report v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Export {
            input,
            html,
            markdown,
            csv,
        } => run(&input, html.as_deref(), markdown.as_deref(), csv.as_deref()),
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
    }
}

fn run(
    input: &std::path::Path,
    html: Option<&std::path::Path>,
    markdown: Option<&std::path::Path>,
    csv: Option<&std::path::Path>,
) -> ExitCode {
    // `-` reads from stdin per Unix convention. Lets users feed
    // results from any source (kubectl / aws s3 cp / a sweep
    // harness) without staging to a temp file.
    //
    // Round-23 sweep: both stdin and path reads are bounded at
    // MAX_BIO_CLI_BYTES (256 MiB) — sister to the round-21 bio CLI
    // cap. A `cat /dev/zero | valenx-report -` would have OOMed
    // pre-fix; stale paths against a multi-GB file have the same
    // shape.
    let (text, display_path): (String, String) = if input == std::path::Path::new("-") {
        match valenx_core::io_caps::read_capped_stdin_to_string(
            valenx_core::io_caps::MAX_BIO_CLI_BYTES,
        ) {
            Ok(s) => (s, "<stdin>".to_string()),
            Err(e) => {
                eprintln!("error: read stdin: {e}");
                return ExitCode::from(3);
            }
        }
    } else {
        match valenx_core::io_caps::read_capped_to_string(
            input,
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        ) {
            Ok(t) => (t, input.display().to_string()),
            Err(e) => {
                eprintln!("error: read {}: {e}", input.display());
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

    if let Some(html_path) = html {
        if let Err(e) = write_html_report(&results, html_path) {
            eprintln!("error: write html {}: {e}", html_path.display());
            return ExitCode::from(3);
        }
        println!("wrote {}", html_path.display());
    }
    if let Some(md_path) = markdown {
        if let Err(e) = write_markdown_report(&results, md_path) {
            eprintln!("error: write markdown {}: {e}", md_path.display());
            return ExitCode::from(3);
        }
        println!("wrote {}", md_path.display());
    }
    if let Some(csv_path) = csv {
        if let Err(e) = write_scalars_csv(&results, csv_path) {
            eprintln!("error: write csv {}: {e}", csv_path.display());
            return ExitCode::from(3);
        }
        println!("wrote {}", csv_path.display());
    }
    ExitCode::from(0)
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
    fn parse_dash_is_stdin_directive() {
        // `-` as the input path means read from stdin per Unix
        // convention. The arg parser must accept it without flagging
        // it as an unknown flag via the dash-prefix branch.
        match parse_args(&["-".into(), "--html".into(), "o.html".into()]) {
            ParsedArgs::Export { input, html, .. } => {
                assert_eq!(input, PathBuf::from("-"));
                assert!(html.is_some());
            }
            other => panic!("expected Export; got {other:?}"),
        }
    }

    #[test]
    fn parse_html_only() {
        match parse_args(&["r.json".into(), "--html".into(), "o.html".into()]) {
            ParsedArgs::Export {
                input,
                html,
                markdown,
                csv,
            } => {
                assert_eq!(input, PathBuf::from("r.json"));
                assert_eq!(html, Some(PathBuf::from("o.html")));
                assert_eq!(markdown, None);
                assert_eq!(csv, None);
            }
            other => panic!("expected Export; got {other:?}"),
        }
    }

    #[test]
    fn parse_markdown_only() {
        match parse_args(&["r.json".into(), "--markdown".into(), "o.md".into()]) {
            ParsedArgs::Export {
                html,
                markdown,
                csv,
                ..
            } => {
                assert_eq!(html, None);
                assert_eq!(markdown, Some(PathBuf::from("o.md")));
                assert_eq!(csv, None);
            }
            other => panic!("expected Export; got {other:?}"),
        }
    }

    #[test]
    fn parse_csv_only() {
        match parse_args(&["r.json".into(), "--csv".into(), "o.csv".into()]) {
            ParsedArgs::Export {
                input,
                html,
                markdown,
                csv,
            } => {
                assert_eq!(input, PathBuf::from("r.json"));
                assert_eq!(html, None);
                assert_eq!(markdown, None);
                assert_eq!(csv, Some(PathBuf::from("o.csv")));
            }
            other => panic!("expected Export; got {other:?}"),
        }
    }

    #[test]
    fn parse_all_three() {
        match parse_args(&[
            "r.json".into(),
            "--html".into(),
            "o.html".into(),
            "--markdown".into(),
            "o.md".into(),
            "--csv".into(),
            "o.csv".into(),
        ]) {
            ParsedArgs::Export {
                html,
                markdown,
                csv,
                ..
            } => {
                assert!(html.is_some());
                assert!(markdown.is_some());
                assert!(csv.is_some());
            }
            other => panic!("expected Export; got {other:?}"),
        }
    }

    #[test]
    fn parse_neither_destination_invalid() {
        match parse_args(&["r.json".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("at least one")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_markdown_without_value_invalid() {
        match parse_args(&["r.json".into(), "--markdown".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--markdown needs")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_markdown_twice_invalid() {
        match parse_args(&[
            "r.json".into(),
            "--markdown".into(),
            "a.md".into(),
            "--markdown".into(),
            "b.md".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--markdown given twice")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_html_without_value_invalid() {
        match parse_args(&["r.json".into(), "--html".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--html needs")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_csv_without_value_invalid() {
        match parse_args(&["r.json".into(), "--csv".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--csv needs")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_html_twice_invalid() {
        match parse_args(&[
            "r.json".into(),
            "--html".into(),
            "a.html".into(),
            "--html".into(),
            "b.html".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--html given twice")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_missing_input_invalid() {
        match parse_args(&["--html".into(), "a.html".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing input path")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_extra_positional_invalid() {
        match parse_args(&[
            "a.json".into(),
            "b.json".into(),
            "--html".into(),
            "o.html".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("extra positional")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_unknown_flag_invalid() {
        match parse_args(&["r.json".into(), "--bogus".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown argument")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }
}

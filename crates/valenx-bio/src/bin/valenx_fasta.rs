//! `valenx-fasta` — headless CLI for FASTA inspection / validation /
//! extraction.
//!
//! Three subcommands wrap [`valenx_bio::format::fasta::read`]:
//!
//! - `inspect` — list every record's name, alphabet, and residue
//!   count. Text output is a fixed-width punch list; `--format json`
//!   re-emits the same data as a single parseable object.
//! - `validate` — parse the file under a chosen alphabet and exit 0
//!   if every byte is valid, or 1 with the offending byte / position
//!   if any record fails validation.
//! - `extract` — print one record (selected by `--name`) as a single
//!   `>name\n<body>\n` chunk, suitable for piping into another tool.
//!
//! The path argument accepts `-` to read from stdin per Unix
//! convention, matching the rest of the Valenx inspector CLIs
//! (`valenx-mesh-info`, `valenx-results`, `valenx-report`).
//!
//! ## Exit codes
//!
//! - 0 — operation succeeded.
//! - 1 — content error (parse failure / invalid alphabet byte / no
//!   record matched `--name`).
//! - 2 — usage error (missing argument, unknown flag, unknown
//!   subcommand, …).
//! - 3 — I/O error reading the input.

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use valenx_bio::format::fasta;
use valenx_bio::{Alphabet, Sequence};

const USAGE: &str = "\
valenx-fasta — inspect / validate / extract FASTA records.

USAGE:
  valenx-fasta inspect <path|->  [--format text|json]
  valenx-fasta validate <path|-> [--alphabet protein|dna|rna]
  valenx-fasta extract <path|->  --name <record-name>
  valenx-fasta -h | --help
  valenx-fasta -V | --version

SUBCOMMANDS:
  inspect    Print one row per record: name | alphabet | residue count.
             `--format json` emits a parseable JSON object.
  validate   Parse the file under `--alphabet` (default `dna`) and exit
             1 with an actionable error on the first invalid byte.
  extract    Print just the record matching `--name` as a one-line-body
             `>name\\n<body>\\n` chunk for piping into another tool.

EXIT CODES:
  0   operation succeeded
  1   content error (parse failure / no record matched / invalid byte)
  2   usage error (missing argument, unknown flag, …)
  3   I/O error reading the input
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

/// Parse the alphabet flag value. `protein` / `dna` / `rna` only —
/// these are the three alphabets [`valenx_bio::Alphabet`] currently
/// supports.
fn parse_alphabet(s: &str) -> Option<Alphabet> {
    match s {
        "protein" => Some(Alphabet::Protein),
        "dna" => Some(Alphabet::Dna),
        "rna" => Some(Alphabet::Rna),
        _ => None,
    }
}

#[derive(Debug, PartialEq)]
enum ParsedArgs {
    Help,
    Version,
    Inspect { path: PathBuf, format: OutputFormat },
    Validate { path: PathBuf, alphabet: Alphabet },
    Extract { path: PathBuf, name: String },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    if args.is_empty() {
        return ParsedArgs::Invalid("missing subcommand (inspect|validate|extract)".into());
    }
    match args[0].as_str() {
        "-h" | "--help" | "help" => return ParsedArgs::Help,
        "-V" | "--version" => return ParsedArgs::Version,
        _ => {}
    }
    let sub = &args[0];
    let rest = &args[1..];
    match sub.as_str() {
        "inspect" => parse_inspect(rest),
        "validate" => parse_validate(rest),
        "extract" => parse_extract(rest),
        other => ParsedArgs::Invalid(format!(
            "unknown subcommand `{other}` (try inspect|validate|extract)"
        )),
    }
}

fn parse_inspect(rest: &[String]) -> ParsedArgs {
    let mut path: Option<PathBuf> = None;
    let mut format = OutputFormat::Text;
    let mut iter = rest.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--format" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--format needs a value (text|json)".into());
                };
                let Some(parsed) = OutputFormat::from_str(v) else {
                    return ParsedArgs::Invalid(format!("unknown format `{v}` (try text or json)"));
                };
                format = parsed;
            }
            s if s.starts_with('-') && s != "-" => {
                return ParsedArgs::Invalid(format!("unknown flag `{s}`"));
            }
            other => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!("unexpected extra argument `{other}`"));
                }
                path = Some(PathBuf::from(other));
            }
        }
    }
    let Some(path) = path else {
        return ParsedArgs::Invalid("inspect: missing FASTA path".into());
    };
    ParsedArgs::Inspect { path, format }
}

fn parse_validate(rest: &[String]) -> ParsedArgs {
    let mut path: Option<PathBuf> = None;
    let mut alphabet = Alphabet::Dna;
    let mut iter = rest.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--alphabet" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid(
                        "--alphabet needs a value (protein|dna|rna)".into(),
                    );
                };
                let Some(parsed) = parse_alphabet(v) else {
                    return ParsedArgs::Invalid(format!(
                        "unknown alphabet `{v}` (try protein, dna, or rna)"
                    ));
                };
                alphabet = parsed;
            }
            s if s.starts_with('-') && s != "-" => {
                return ParsedArgs::Invalid(format!("unknown flag `{s}`"));
            }
            other => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!("unexpected extra argument `{other}`"));
                }
                path = Some(PathBuf::from(other));
            }
        }
    }
    let Some(path) = path else {
        return ParsedArgs::Invalid("validate: missing FASTA path".into());
    };
    ParsedArgs::Validate { path, alphabet }
}

fn parse_extract(rest: &[String]) -> ParsedArgs {
    let mut path: Option<PathBuf> = None;
    let mut name: Option<String> = None;
    let mut iter = rest.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--name" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--name needs a value".into());
                };
                if name.is_some() {
                    return ParsedArgs::Invalid("--name given twice".into());
                }
                name = Some(v.clone());
            }
            s if s.starts_with('-') && s != "-" => {
                return ParsedArgs::Invalid(format!("unknown flag `{s}`"));
            }
            other => {
                if path.is_some() {
                    return ParsedArgs::Invalid(format!("unexpected extra argument `{other}`"));
                }
                path = Some(PathBuf::from(other));
            }
        }
    }
    let Some(path) = path else {
        return ParsedArgs::Invalid("extract: missing FASTA path".into());
    };
    let Some(name) = name else {
        return ParsedArgs::Invalid("extract: missing --name argument".into());
    };
    ParsedArgs::Extract { path, name }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        ParsedArgs::Help => {
            print!("{USAGE}");
            ExitCode::from(0)
        }
        ParsedArgs::Version => {
            println!("valenx-fasta v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
        ParsedArgs::Inspect { path, format } => match run_inspect(&path, format) {
            Ok(()) => ExitCode::from(0),
            Err(e) => {
                eprintln!("error: {e}");
                e.exit_code()
            }
        },
        ParsedArgs::Validate { path, alphabet } => match run_validate(&path, alphabet) {
            Ok(()) => ExitCode::from(0),
            Err(e) => {
                eprintln!("error: {e}");
                e.exit_code()
            }
        },
        ParsedArgs::Extract { path, name } => match run_extract(&path, &name) {
            Ok(()) => ExitCode::from(0),
            Err(e) => {
                eprintln!("error: {e}");
                e.exit_code()
            }
        },
    }
}

#[derive(Debug)]
enum CliError {
    Io(String),
    Content(String),
}

impl CliError {
    fn exit_code(&self) -> ExitCode {
        match self {
            CliError::Io(_) => ExitCode::from(3),
            CliError::Content(_) => ExitCode::from(1),
        }
    }
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CliError::Io(s) | CliError::Content(s) => f.write_str(s),
        }
    }
}

/// Read the FASTA bytes from `path`, treating `-` as stdin per the
/// Unix convention. Returns the text plus a display-friendly path.
fn read_input(path: &Path) -> Result<(String, String), CliError> {
    // Round-21 M4: see valenx_blast for the cap rationale.
    if path == Path::new("-") {
        let buf = valenx_core::io_caps::read_capped_stdin_to_string(
            valenx_core::io_caps::MAX_BIO_CLI_BYTES,
        )
        .map_err(|e| CliError::Io(format!("read stdin: {e}")))?;
        Ok((buf, "<stdin>".to_string()))
    } else {
        let text = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        )
        .map_err(|e| CliError::Io(format!("read {}: {e}", path.display())))?;
        Ok((text, path.display().to_string()))
    }
}

/// Auto-detect alphabet by parsing the text against each candidate in
/// turn. Falls back to protein when nothing else fits — protein is a
/// strict superset of DNA / RNA in practice (every nucleotide letter
/// is also a valid amino-acid one-letter code), so a parse-as-protein
/// only fails on bytes outside the IUPAC alphabet entirely.
///
/// Returns the detected alphabet AND the parsed records; the caller
/// avoids re-parsing.
fn detect_and_parse(text: &str) -> Result<(Alphabet, Vec<Sequence>), CliError> {
    // Try DNA first — narrowest alphabet, so a clean DNA parse means
    // the data really is DNA.
    if let Ok(seqs) = fasta::read(text, Alphabet::Dna) {
        return Ok((Alphabet::Dna, seqs));
    }
    // Then RNA.
    if let Ok(seqs) = fasta::read(text, Alphabet::Rna) {
        return Ok((Alphabet::Rna, seqs));
    }
    // Then protein, the broadest.
    match fasta::read(text, Alphabet::Protein) {
        Ok(seqs) => Ok((Alphabet::Protein, seqs)),
        Err(e) => Err(CliError::Content(format!("parse FASTA: {e}"))),
    }
}

fn run_inspect(path: &Path, format: OutputFormat) -> Result<(), CliError> {
    let (text, display_path) = read_input(path)?;
    // Inspect runs in auto-detect mode: show whatever alphabet each
    // file matched. Validation is the user's call via `validate`.
    let (_alphabet, seqs) = detect_and_parse(&text)?;
    match format {
        OutputFormat::Text => {
            print_inspect_text(&display_path, &seqs);
        }
        OutputFormat::Json => {
            println!("{}", render_inspect_json(&display_path, &seqs));
        }
    }
    Ok(())
}

fn print_inspect_text(path: &str, seqs: &[Sequence]) {
    let basename = Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    println!("{basename}: {} records", seqs.len());
    let max_name = seqs.iter().map(|s| s.name.len()).max().unwrap_or(0);
    let max_alpha = seqs
        .iter()
        .map(|s| s.alphabet.id().len())
        .max()
        .unwrap_or(0);
    let max_len_digits = seqs
        .iter()
        .map(|s| s.len().to_string().len())
        .max()
        .unwrap_or(1);
    for s in seqs {
        let unit = if s.alphabet == Alphabet::Protein {
            "residues"
        } else {
            "bases"
        };
        println!(
            "  {name:nw$} | {alpha:aw$} | {len:>lw$} {unit}",
            name = s.name,
            alpha = s.alphabet.id(),
            len = s.len(),
            unit = unit,
            nw = max_name,
            aw = max_alpha,
            lw = max_len_digits,
        );
    }
}

fn render_inspect_json(path: &str, seqs: &[Sequence]) -> String {
    let records: Vec<serde_json::Value> = seqs
        .iter()
        .map(|s| {
            serde_json::json!({
                "name": s.name,
                "alphabet": s.alphabet.id(),
                "length": s.len(),
            })
        })
        .collect();
    let v = serde_json::json!({
        "path": path,
        "records": records,
    });
    serde_json::to_string_pretty(&v).expect("static structure serialises")
}

fn run_validate(path: &Path, alphabet: Alphabet) -> Result<(), CliError> {
    let (text, _display_path) = read_input(path)?;
    // Bubble fasta::read's error up verbatim — it already includes
    // byte / position / alphabet / record-name info as the spec asks.
    fasta::read(&text, alphabet).map_err(|e| CliError::Content(format!("{e}")))?;
    Ok(())
}

fn run_extract(path: &Path, name: &str) -> Result<(), CliError> {
    let (text, _display_path) = read_input(path)?;
    let (_alphabet, seqs) = detect_and_parse(&text)?;
    let Some(found) = seqs.iter().find(|s| s.name == name) else {
        return Err(CliError::Content(format!(
            "no record named `{name}` in input ({} records found)",
            seqs.len()
        )));
    };
    // Single-line body output — easier to pipe into other tools that
    // expect FASTA-shaped input. `fasta::write` wraps at 60 chars; we
    // emit a one-line body here on purpose.
    println!(">{}", found.name);
    println!("{}", found.as_str());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_args_help_variants() {
        for arg in ["-h", "--help", "help"] {
            assert_eq!(parse_args(&[arg.into()]), ParsedArgs::Help);
        }
    }

    #[test]
    fn parse_args_version_variants() {
        for arg in ["-V", "--version"] {
            assert_eq!(parse_args(&[arg.into()]), ParsedArgs::Version);
        }
    }

    #[test]
    fn parse_args_no_subcommand_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(msg) => {
                assert!(msg.contains("missing subcommand"), "got: {msg}");
            }
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_subcommand_invalid() {
        match parse_args(&["frobnicate".into(), "x".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("frobnicate")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_default_format() {
        match parse_args(&["inspect".into(), "x.fasta".into()]) {
            ParsedArgs::Inspect { path, format } => {
                assert_eq!(path, PathBuf::from("x.fasta"));
                assert_eq!(format, OutputFormat::Text);
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_format_json() {
        match parse_args(&[
            "inspect".into(),
            "x.fasta".into(),
            "--format".into(),
            "json".into(),
        ]) {
            ParsedArgs::Inspect { format, .. } => assert_eq!(format, OutputFormat::Json),
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_format_unknown_invalid() {
        match parse_args(&[
            "inspect".into(),
            "x.fasta".into(),
            "--format".into(),
            "yaml".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("yaml"), "got: {msg}"),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_missing_path_invalid() {
        match parse_args(&["inspect".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing FASTA path")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_inspect_dash_is_stdin() {
        // `-` for stdin must NOT be treated as an unknown flag.
        match parse_args(&["inspect".into(), "-".into()]) {
            ParsedArgs::Inspect { path, .. } => {
                assert_eq!(path, PathBuf::from("-"));
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_validate_default_alphabet_is_dna() {
        match parse_args(&["validate".into(), "x.fasta".into()]) {
            ParsedArgs::Validate { alphabet, .. } => assert_eq!(alphabet, Alphabet::Dna),
            other => panic!("expected Validate; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_validate_alphabet_protein() {
        match parse_args(&[
            "validate".into(),
            "x.fasta".into(),
            "--alphabet".into(),
            "protein".into(),
        ]) {
            ParsedArgs::Validate { alphabet, .. } => assert_eq!(alphabet, Alphabet::Protein),
            other => panic!("expected Validate; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_validate_unknown_alphabet_invalid() {
        match parse_args(&[
            "validate".into(),
            "x.fasta".into(),
            "--alphabet".into(),
            "klingon".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("klingon")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extract_with_name() {
        match parse_args(&[
            "extract".into(),
            "x.fasta".into(),
            "--name".into(),
            "p53".into(),
        ]) {
            ParsedArgs::Extract { path, name } => {
                assert_eq!(path, PathBuf::from("x.fasta"));
                assert_eq!(name, "p53");
            }
            other => panic!("expected Extract; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extract_missing_name_invalid() {
        match parse_args(&["extract".into(), "x.fasta".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--name")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extract_missing_path_invalid() {
        match parse_args(&["extract".into(), "--name".into(), "x".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing FASTA path")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_alphabet_recognises_each() {
        assert_eq!(parse_alphabet("dna"), Some(Alphabet::Dna));
        assert_eq!(parse_alphabet("rna"), Some(Alphabet::Rna));
        assert_eq!(parse_alphabet("protein"), Some(Alphabet::Protein));
        assert_eq!(parse_alphabet("klingon"), None);
    }

    #[test]
    fn detect_and_parse_picks_dna_for_dna_input() {
        let (alpha, seqs) = detect_and_parse(">x\nACGTACGT\n").unwrap();
        assert_eq!(alpha, Alphabet::Dna);
        assert_eq!(seqs.len(), 1);
    }

    #[test]
    fn detect_and_parse_picks_protein_for_protein_only_residues() {
        // `E` (glutamate) is not a valid IUPAC nucleotide code, so DNA
        // and RNA both reject — protein wins.
        let (alpha, _seqs) = detect_and_parse(">x\nMEEPQSDPSV\n").unwrap();
        assert_eq!(alpha, Alphabet::Protein);
    }

    #[test]
    fn render_inspect_json_round_trips() {
        let s = Sequence::new("p53", Alphabet::Protein, "MEEP").unwrap();
        let out = render_inspect_json("x.fasta", &[s]);
        let v: serde_json::Value = serde_json::from_str(&out).expect("parseable JSON");
        assert_eq!(v["path"], "x.fasta");
        assert_eq!(v["records"][0]["name"], "p53");
        assert_eq!(v["records"][0]["alphabet"], "protein");
        assert_eq!(v["records"][0]["length"], 4);
    }

    #[test]
    fn run_validate_passes_for_clean_dna() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fasta-validate-ok-{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, ">x\nACGTACGT\n").unwrap();
        let r = run_validate(&tmp, Alphabet::Dna);
        assert!(r.is_ok(), "expected Ok, got {r:?}");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn run_validate_rejects_protein_under_dna_alphabet() {
        // Glutamate `E` is invalid in DNA — error has alphabet+position.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fasta-validate-bad-{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, ">p53\nMEEPQSDPSV\n").unwrap();
        let r = run_validate(&tmp, Alphabet::Dna);
        match r {
            Err(CliError::Content(msg)) => {
                assert!(msg.contains("position"), "expected position in: {msg}");
            }
            other => panic!("expected Content error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn run_extract_returns_content_error_on_missing_name() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fasta-extract-{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, ">p53\nMEEP\n").unwrap();
        let r = run_extract(&tmp, "ubq");
        match r {
            Err(CliError::Content(msg)) => assert!(msg.contains("ubq")),
            other => panic!("expected Content error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn cli_error_exit_codes_map_correctly() {
        assert_eq!(CliError::Io("x".into()).exit_code(), ExitCode::from(3));
        assert_eq!(CliError::Content("x".into()).exit_code(), ExitCode::from(1));
    }
}

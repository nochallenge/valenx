//! `valenx-blast` — headless CLI wrapping the BLAST+ binaries.
//!
//! Wraps `blastp` (protein) and `blastn` (nucleotide) from the NCBI
//! BLAST+ distribution. We do NOT bind libblast — the BLAST+ tools
//! are well-tested as standalone CLIs and shipping a static link
//! would balloon binary size and complicate the build matrix.
//! Instead, the wrapper:
//!
//! 1. Reads the query (FASTA, or `-` for stdin) and auto-detects the
//!    alphabet by trying DNA → RNA → protein. RNA queries route to
//!    `blastn` (BLAST+ accepts U-as-T uracil/thymine for RNA
//!    queries against DNA databases).
//! 2. Locates `blastp` / `blastn` on PATH via
//!    [`valenx_core::adapter_helpers::find_on_path`]. If neither is
//!    installed, exits 3 with an actionable hint.
//! 3. Spawns the chosen binary with `-outfmt 6` (tab-separated, the
//!    BLAST+ canonical machine-readable form) and captures stdout.
//! 4. In text mode prints the captured stdout verbatim (the user
//!    gets exactly what `blastp` / `blastn` would have printed); in
//!    `--format json` mode parses each line into a structured
//!    record and emits a single parseable JSON document.
//!
//! Tests can't actually run BLAST+ in CI (the runner won't have the
//! binary or a database). Integration coverage stops at help /
//! version / missing-binary path / empty-query rejection — the
//! happy-path BLAST+ run lives in adapter-side smoke tests once the
//! BLAST+ adapter lands.
//!
//! ## Exit codes
//!
//! - 0 — BLAST+ ran and stdout was captured (or printed).
//! - 1 — content error (empty query / unparseable query / BLAST+
//!   exited non-zero with a parse failure on its stderr).
//! - 2 — usage error (missing argument, unknown flag, …).
//! - 3 — I/O error (PATH probe failed / spawn failed / can't read
//!   the input).

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use valenx_bio::format::fasta;
use valenx_bio::{Alphabet, Sequence};
use valenx_core::adapter_helpers::find_on_path;

const USAGE: &str = "\
valenx-blast — wrap NCBI BLAST+ binaries (blastp / blastn).

USAGE:
  valenx-blast <query.fasta|-> --db <db_prefix>
  valenx-blast <query.fasta|-> --db <db> --evalue 1e-5 --max-target-seqs 10
  valenx-blast <query.fasta|-> --db <db> --format json
  valenx-blast -h | --help
  valenx-blast -V | --version

OPTIONS:
  --db PATH               BLAST+ database prefix. Passed verbatim to
                          `blastp`/`blastn -db`. We don't try to
                          detect or stage databases — the user is
                          responsible for the format/index.
  --evalue FLOAT          E-value threshold. Surfaced as `-evalue` in
                          the BLAST+ command. Optional; BLAST+'s own
                          default applies when omitted.
  --max-target-seqs N     Maximum number of target sequences to
                          report. Surfaced as `-max_target_seqs`.
                          Optional.
  --format F              Output shape. `text` (default) prints
                          BLAST+'s `-outfmt 6` tab-separated stdout
                          verbatim; `json` parses each line into a
                          structured `{query, subject, percent_identity,
                          length, mismatch, gapopen, qstart, qend,
                          sstart, send, evalue, bitscore}` record and
                          emits a single JSON document.
  -h, --help              Show this help.
  -V, --version           Print the binary version and exit.

The query alphabet auto-detects: DNA / RNA queries route to `blastn`,
protein queries to `blastp`. Override is not supported — the BLAST+
binary that matches your sequence type is always the right call.

EXIT CODES:
  0   BLAST+ ran and stdout was captured (or printed)
  1   content error (empty query / unparseable query / BLAST+ stderr
      surfaced an input-data problem)
  2   usage error
  3   I/O error (BLAST+ not installed / spawn failed / can't read input)
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

#[derive(Debug, PartialEq)]
enum ParsedArgs {
    Help,
    Version,
    Run {
        query: PathBuf,
        db: String,
        evalue: Option<f64>,
        max_target_seqs: Option<u64>,
        format: OutputFormat,
    },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    if args.is_empty() {
        return ParsedArgs::Invalid("missing query path".into());
    }
    if matches!(args[0].as_str(), "-h" | "--help" | "help") {
        return ParsedArgs::Help;
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        return ParsedArgs::Version;
    }
    let mut query: Option<PathBuf> = None;
    let mut db: Option<String> = None;
    let mut evalue: Option<f64> = None;
    let mut max_target_seqs: Option<u64> = None;
    let mut format = OutputFormat::Text;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--db" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--db needs a value".into());
                };
                if db.is_some() {
                    return ParsedArgs::Invalid("--db given twice".into());
                }
                db = Some(v.clone());
            }
            "--evalue" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--evalue needs a value".into());
                };
                let Ok(parsed) = v.parse::<f64>() else {
                    return ParsedArgs::Invalid(format!("--evalue `{v}` is not a number"));
                };
                evalue = Some(parsed);
            }
            "--max-target-seqs" => {
                let Some(v) = iter.next() else {
                    return ParsedArgs::Invalid("--max-target-seqs needs a value".into());
                };
                let Ok(parsed) = v.parse::<u64>() else {
                    return ParsedArgs::Invalid(format!(
                        "--max-target-seqs `{v}` is not a positive integer"
                    ));
                };
                max_target_seqs = Some(parsed);
            }
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
                if query.is_some() {
                    return ParsedArgs::Invalid(format!("unexpected extra argument `{other}`"));
                }
                query = Some(PathBuf::from(other));
            }
        }
    }
    let Some(query) = query else {
        return ParsedArgs::Invalid("missing query path".into());
    };
    let Some(db) = db else {
        return ParsedArgs::Invalid("missing --db argument".into());
    };
    ParsedArgs::Run {
        query,
        db,
        evalue,
        max_target_seqs,
        format,
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
            println!("valenx-blast v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
        ParsedArgs::Run {
            query,
            db,
            evalue,
            max_target_seqs,
            format,
        } => match run(&query, &db, evalue, max_target_seqs, format) {
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

/// Read the query bytes (file or stdin), then parse to detect the
/// alphabet. The `Alphabet` choice picks the BLAST+ binary; the
/// returned `String` is the raw FASTA text we hand to BLAST+ on
/// stdin (preserving the original record structure).
fn load_query(path: &Path) -> Result<(String, Alphabet, Vec<Sequence>), CliError> {
    // Round-21 M4: bound both stdin and file reads at
    // MAX_BIO_CLI_BYTES (256 MiB). Pre-fix the bare
    // `std::io::stdin().read_to_string` would slurp a
    // `cat /dev/zero | valenx-blast` pipe unbounded.
    let text: String = if path == Path::new("-") {
        valenx_core::io_caps::read_capped_stdin_to_string(valenx_core::io_caps::MAX_BIO_CLI_BYTES)
            .map_err(|e| CliError::Io(format!("read stdin: {e}")))?
    } else {
        valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        )
        .map_err(|e| CliError::Io(format!("read {}: {e}", path.display())))?
    };
    if text.trim().is_empty() {
        return Err(CliError::Content("query is empty".into()));
    }
    let (alphabet, seqs) = detect_alphabet(&text)?;
    if seqs.is_empty() {
        return Err(CliError::Content("query has zero records".into()));
    }
    Ok((text, alphabet, seqs))
}

/// Try DNA → RNA → protein in turn. The first alphabet that parses
/// cleanly wins. Protein is the broadest catch-all; DNA is checked
/// first so a clean nucleotide query routes to `blastn` rather than
/// `blastp`.
fn detect_alphabet(text: &str) -> Result<(Alphabet, Vec<Sequence>), CliError> {
    if let Ok(seqs) = fasta::read(text, Alphabet::Dna) {
        return Ok((Alphabet::Dna, seqs));
    }
    if let Ok(seqs) = fasta::read(text, Alphabet::Rna) {
        return Ok((Alphabet::Rna, seqs));
    }
    match fasta::read(text, Alphabet::Protein) {
        Ok(seqs) => Ok((Alphabet::Protein, seqs)),
        Err(e) => Err(CliError::Content(format!("parse query: {e}"))),
    }
}

/// Pick `blastp` (Protein) or `blastn` (Dna / Rna). RNA routes to
/// blastn — BLAST+ accepts U-as-T translation transparently.
fn binary_for(alphabet: Alphabet) -> &'static [&'static str] {
    match alphabet {
        Alphabet::Protein => &["blastp"],
        Alphabet::Dna | Alphabet::Rna => &["blastn"],
    }
}

fn run(
    query_path: &Path,
    db: &str,
    evalue: Option<f64>,
    max_target_seqs: Option<u64>,
    format: OutputFormat,
) -> Result<(), CliError> {
    let (query_text, alphabet, _seqs) = load_query(query_path)?;
    let candidates = binary_for(alphabet);
    let binary = find_on_path(candidates).ok_or_else(|| {
        CliError::Io(format!(
            "BLAST+ binary not found on PATH (tried {:?}); \
             install NCBI BLAST+ and ensure {} is on PATH",
            candidates, candidates[0]
        ))
    })?;

    let mut cmd = std::process::Command::new(&binary);
    cmd.arg("-db")
        .arg(db)
        .arg("-outfmt")
        .arg("6")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    if let Some(e) = evalue {
        cmd.arg("-evalue").arg(e.to_string());
    }
    if let Some(n) = max_target_seqs {
        cmd.arg("-max_target_seqs").arg(n.to_string());
    }

    let mut child = cmd
        .spawn()
        .map_err(|e| CliError::Io(format!("spawn {}: {e}", binary.display())))?;
    {
        use std::io::Write;
        let stdin = child
            .stdin
            .as_mut()
            .ok_or_else(|| CliError::Io("BLAST+ stdin pipe missing".into()))?;
        stdin
            .write_all(query_text.as_bytes())
            .map_err(|e| CliError::Io(format!("write query stdin: {e}")))?;
    }
    let out = child
        .wait_with_output()
        .map_err(|e| CliError::Io(format!("wait BLAST+: {e}")))?;

    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(CliError::Content(format!(
            "{} exited with {:?}: {}",
            binary.display(),
            out.status.code(),
            stderr.trim()
        )));
    }

    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    match format {
        OutputFormat::Text => {
            print!("{stdout}");
        }
        OutputFormat::Json => {
            println!("{}", render_json(&stdout));
        }
    }
    Ok(())
}

/// Parse BLAST+'s `-outfmt 6` tab-separated stdout into a JSON
/// document. The 12 columns are the BLAST+ default: query, subject,
/// percent_identity, length, mismatch, gapopen, qstart, qend, sstart,
/// send, evalue, bitscore. Lines with fewer columns are skipped (the
/// BLAST+ "no hits" output is one informational line that doesn't
/// match — quietly ignored).
fn render_json(stdout: &str) -> String {
    let hits: Vec<serde_json::Value> = stdout
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            if cols.len() < 12 {
                return None;
            }
            Some(serde_json::json!({
                "query": cols[0],
                "subject": cols[1],
                "percent_identity": cols[2].parse::<f64>().unwrap_or(0.0),
                "length": cols[3].parse::<u64>().unwrap_or(0),
                "mismatch": cols[4].parse::<u64>().unwrap_or(0),
                "gapopen": cols[5].parse::<u64>().unwrap_or(0),
                "qstart": cols[6].parse::<u64>().unwrap_or(0),
                "qend": cols[7].parse::<u64>().unwrap_or(0),
                "sstart": cols[8].parse::<u64>().unwrap_or(0),
                "send": cols[9].parse::<u64>().unwrap_or(0),
                "evalue": cols[10].parse::<f64>().unwrap_or(0.0),
                "bitscore": cols[11].parse::<f64>().unwrap_or(0.0),
            }))
        })
        .collect();
    let v = serde_json::json!({ "hits": hits });
    serde_json::to_string_pretty(&v).expect("static structure serialises")
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
    fn parse_args_missing_query_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing query")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_missing_db_invalid() {
        match parse_args(&["q.fasta".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--db")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_minimal_run() {
        match parse_args(&["q.fasta".into(), "--db".into(), "nr".into()]) {
            ParsedArgs::Run {
                query,
                db,
                evalue,
                max_target_seqs,
                format,
            } => {
                assert_eq!(query, PathBuf::from("q.fasta"));
                assert_eq!(db, "nr");
                assert!(evalue.is_none());
                assert!(max_target_seqs.is_none());
                assert_eq!(format, OutputFormat::Text);
            }
            other => panic!("expected Run; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_full_set_of_flags() {
        match parse_args(&[
            "q.fasta".into(),
            "--db".into(),
            "nr".into(),
            "--evalue".into(),
            "1e-5".into(),
            "--max-target-seqs".into(),
            "10".into(),
            "--format".into(),
            "json".into(),
        ]) {
            ParsedArgs::Run {
                evalue,
                max_target_seqs,
                format,
                ..
            } => {
                assert!(evalue.is_some());
                assert_eq!(max_target_seqs, Some(10));
                assert_eq!(format, OutputFormat::Json);
            }
            other => panic!("expected Run; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_evalue_not_a_number_invalid() {
        match parse_args(&[
            "q.fasta".into(),
            "--db".into(),
            "nr".into(),
            "--evalue".into(),
            "not-a-number".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("not a number")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_flag_invalid() {
        match parse_args(&[
            "q.fasta".into(),
            "--db".into(),
            "nr".into(),
            "--bogus".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown flag")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_dash_is_stdin() {
        match parse_args(&["-".into(), "--db".into(), "nr".into()]) {
            ParsedArgs::Run { query, .. } => assert_eq!(query, PathBuf::from("-")),
            other => panic!("expected Run; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_unknown_invalid() {
        match parse_args(&[
            "q.fasta".into(),
            "--db".into(),
            "nr".into(),
            "--format".into(),
            "yaml".into(),
        ]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("yaml")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn detect_alphabet_protein_routes_to_blastp() {
        let (alpha, _seqs) = detect_alphabet(">x\nMEEPQSDPSV\n").unwrap();
        assert_eq!(alpha, Alphabet::Protein);
        assert_eq!(binary_for(alpha), &["blastp"]);
    }

    #[test]
    fn detect_alphabet_dna_routes_to_blastn() {
        let (alpha, _seqs) = detect_alphabet(">x\nACGTACGT\n").unwrap();
        assert_eq!(alpha, Alphabet::Dna);
        assert_eq!(binary_for(alpha), &["blastn"]);
    }

    #[test]
    fn render_json_parses_outfmt6_tab_lines() {
        // BLAST+ outfmt 6 sample: 12 tab-separated columns.
        let stdout = "Query_1\tsubj_a\t99.5\t100\t1\t0\t1\t100\t1\t100\t1e-50\t180.0\n\
                      Query_1\tsubj_b\t75.2\t98\t10\t1\t1\t98\t5\t102\t1e-20\t90.0\n";
        let s = render_json(stdout);
        let v: serde_json::Value = serde_json::from_str(&s).expect("parseable JSON");
        let hits = v["hits"].as_array().expect("hits is array");
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0]["query"], "Query_1");
        assert_eq!(hits[0]["subject"], "subj_a");
        assert_eq!(hits[0]["length"], 100);
        assert_eq!(hits[1]["subject"], "subj_b");
    }

    #[test]
    fn render_json_skips_short_lines() {
        // BLAST+'s "BLAST 0 hits found" message is a one-line
        // informational header that doesn't have 12 tab-separated
        // columns. render_json should drop it without panicking.
        let stdout = "# BLAST 0 hits found\n";
        let s = render_json(stdout);
        let v: serde_json::Value = serde_json::from_str(&s).expect("parseable JSON");
        let hits = v["hits"].as_array().expect("hits is array");
        assert!(hits.is_empty());
    }

    #[test]
    fn binary_for_rna_is_blastn() {
        assert_eq!(binary_for(Alphabet::Rna), &["blastn"]);
    }

    #[test]
    fn cli_error_exit_codes_map_correctly() {
        assert_eq!(CliError::Io("x".into()).exit_code(), ExitCode::from(3));
        assert_eq!(CliError::Content("x".into()).exit_code(), ExitCode::from(1));
    }

    #[test]
    fn load_query_rejects_empty_input() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-blast-empty-{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, "").unwrap();
        let r = load_query(&tmp);
        match r {
            Err(CliError::Content(msg)) => assert!(msg.contains("empty")),
            other => panic!("expected Content error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-21 M4 RED→GREEN: an oversize query file is rejected
    /// with `CliError::Io` from the capped reader rather than being
    /// slurped into memory. Pre-fix the bare `fs::read_to_string`
    /// would have allocated the whole file before the alphabet
    /// detector ran. Use a small cap-shaped scratch file rather
    /// than allocating 257 MiB on disk in CI.
    #[test]
    fn load_query_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-blast-r21-oversize-{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // Sparse-allocate just past the 256 MiB cap.
        let f = std::fs::File::create(&tmp).unwrap();
        f.set_len(valenx_core::io_caps::MAX_BIO_CLI_BYTES + 1)
            .unwrap();
        drop(f);
        let r = load_query(&tmp);
        match r {
            Err(CliError::Io(msg)) => {
                assert!(
                    msg.contains("read") || msg.contains("exceed") || msg.contains("cap"),
                    "expected cap-related io error, got: {msg}"
                );
            }
            other => panic!("expected Io error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}

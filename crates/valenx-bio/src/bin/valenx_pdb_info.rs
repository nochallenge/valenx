//! `valenx-pdb-info` — headless CLI for PDB structure summaries.
//!
//! Reads a PDB file (ATOM / HETATM records) via
//! [`valenx_bio::format::pdb::read`] and prints a quick punch list of
//! per-chain residue / atom counts plus a per-chain element-frequency
//! tally. Useful for confirming a PDB landed correctly in CI, for
//! `gh codespace` debugging where the GUI isn't reachable, or as a
//! quick `head` on a multi-MB structure dump.
//!
//! ## Output
//!
//! Text (default):
//!
//! ```text
//! PDB id: 1ubq-tiny (1 chain, 2 residues, 9 atoms)
//!   chain A: 2 residues, 9 atoms
//!     residue range: 1..=2
//!     elements: 4 C, 2 N, 2 O
//! ```
//!
//! JSON (`--format json`):
//!
//! ```json
//! {
//!   "id": "1ubq-tiny",
//!   "chains": [{"id":"A","residues":2,"atoms":9}],
//!   "atom_count": 9,
//!   "residue_count": 2
//! }
//! ```
//!
//! ## Exit codes
//!
//! - 0 — file loaded and printed.
//! - 1 — content error (PDB parse failure).
//! - 2 — usage error (missing argument, unknown flag, …).
//! - 3 — I/O error reading the input.

#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use valenx_bio::format::pdb;
use valenx_bio::Structure;

const USAGE: &str = "\
valenx-pdb-info — print a structure summary for a PDB file.

USAGE:
  valenx-pdb-info <path|->  [--format text|json]
  valenx-pdb-info -                       # read from stdin
  valenx-pdb-info -h | --help
  valenx-pdb-info -V | --version

OPTIONS:
  --format F     Output format: `text` (default, human-readable) or
                 `json` (single parseable object — id + per-chain
                 counts + total atom / residue counts).
  -h, --help     Show this help.
  -V, --version  Print the binary version and exit.

EXIT CODES:
  0   summary printed
  1   PDB parse failure
  2   usage error
  3   I/O error
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
    Inspect { path: PathBuf, format: OutputFormat },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    if args.is_empty() {
        return ParsedArgs::Invalid("missing PDB file path".into());
    }
    if matches!(args[0].as_str(), "-h" | "--help" | "help") {
        return ParsedArgs::Help;
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        return ParsedArgs::Version;
    }
    let mut path: Option<PathBuf> = None;
    let mut format = OutputFormat::Text;
    let mut iter = args.iter();
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
        return ParsedArgs::Invalid("missing PDB file path".into());
    };
    ParsedArgs::Inspect { path, format }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        ParsedArgs::Help => {
            print!("{USAGE}");
            ExitCode::from(0)
        }
        ParsedArgs::Version => {
            println!("valenx-pdb-info v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
        ParsedArgs::Inspect { path, format } => match inspect(&path, format) {
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
    Content(String),
}

impl InspectError {
    fn exit_code(&self) -> ExitCode {
        match self {
            InspectError::Io(_) => ExitCode::from(3),
            InspectError::Content(_) => ExitCode::from(1),
        }
    }
}

impl std::fmt::Display for InspectError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InspectError::Io(s) | InspectError::Content(s) => f.write_str(s),
        }
    }
}

/// Derive a default structure id from the input path's stem, falling
/// back to "stdin" when reading from `-`. Mirrors how the test fixture
/// is named (`1ubq-tiny.pdb` -> id "1ubq-tiny").
fn id_from_path(path: &Path) -> String {
    if path == Path::new("-") {
        return "stdin".to_string();
    }
    path.file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| "structure".to_string())
}

fn inspect(path: &Path, format: OutputFormat) -> Result<(), InspectError> {
    // `-` reads from stdin per Unix convention, matching every other
    // Valenx inspector CLI.
    //
    // Round-21 M4: bound stdin and file reads at MAX_BIO_CLI_BYTES
    // so a piped `/dev/zero` or stale path can't OOM the inspector.
    let text: String = if path == Path::new("-") {
        valenx_core::io_caps::read_capped_stdin_to_string(valenx_core::io_caps::MAX_BIO_CLI_BYTES)
            .map_err(|e| InspectError::Io(format!("read stdin: {e}")))?
    } else {
        valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        )
        .map_err(|e| InspectError::Io(format!("read {}: {e}", path.display())))?
    };
    let id = id_from_path(path);
    let structure = pdb::read(&id, &text).map_err(|e| InspectError::Content(format!("{e}")))?;
    match format {
        OutputFormat::Text => print_text(&structure),
        OutputFormat::Json => println!("{}", render_json(&structure)),
    }
    Ok(())
}

/// Tally element symbols across a chain into a deterministic
/// `(symbol, count)` list (alphabetically sorted on symbol so the
/// output is stable for grep/awk pipelines).
fn element_tally(chain: &valenx_bio::Chain) -> Vec<(String, usize)> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for residue in &chain.residues {
        for atom in &residue.atoms {
            *counts.entry(atom.element.clone()).or_insert(0) += 1;
        }
    }
    counts.into_iter().collect()
}

fn print_text(structure: &Structure) {
    let chain_count = structure.chains.len();
    println!(
        "PDB id: {} ({} chain{}, {} residues, {} atoms)",
        structure.id,
        chain_count,
        if chain_count == 1 { "" } else { "s" },
        structure.residue_count(),
        structure.atom_count(),
    );
    for chain in &structure.chains {
        let n_atoms: usize = chain.residues.iter().map(|r| r.atoms.len()).sum();
        println!(
            "  chain {}: {} residues, {n_atoms} atoms",
            chain.id,
            chain.residues.len(),
        );
        if let (Some(first), Some(last)) = (chain.residues.first(), chain.residues.last()) {
            println!("    residue range: {}..={}", first.seq_id, last.seq_id);
        }
        let tally = element_tally(chain);
        if !tally.is_empty() {
            let parts: Vec<String> = tally.iter().map(|(sym, n)| format!("{n} {sym}")).collect();
            println!("    elements: {}", parts.join(", "));
        }
    }
}

fn render_json(structure: &Structure) -> String {
    let chains: Vec<serde_json::Value> = structure
        .chains
        .iter()
        .map(|c| {
            let n_atoms: usize = c.residues.iter().map(|r| r.atoms.len()).sum();
            serde_json::json!({
                "id": c.id.to_string(),
                "residues": c.residues.len(),
                "atoms": n_atoms,
            })
        })
        .collect();
    let v = serde_json::json!({
        "id": structure.id,
        "chains": chains,
        "atom_count": structure.atom_count(),
        "residue_count": structure.residue_count(),
    });
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
    fn parse_args_missing_path_invalid() {
        match parse_args(&[]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("missing PDB file")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_picks_up_path_with_default_format() {
        match parse_args(&["x.pdb".into()]) {
            ParsedArgs::Inspect { path, format } => {
                assert_eq!(path, PathBuf::from("x.pdb"));
                assert_eq!(format, OutputFormat::Text);
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_json_flag() {
        match parse_args(&["x.pdb".into(), "--format".into(), "json".into()]) {
            ParsedArgs::Inspect { format, .. } => assert_eq!(format, OutputFormat::Json),
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_format_invalid() {
        match parse_args(&["x.pdb".into(), "--format".into(), "yaml".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("yaml")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_extra_arg_invalid() {
        match parse_args(&["a.pdb".into(), "b.pdb".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unexpected extra")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_dash_is_stdin() {
        // `-` for stdin must NOT be treated as an unknown flag.
        match parse_args(&["-".into()]) {
            ParsedArgs::Inspect { path, .. } => {
                assert_eq!(path, PathBuf::from("-"));
            }
            other => panic!("expected Inspect; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_format_missing_value_invalid() {
        match parse_args(&["x.pdb".into(), "--format".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--format")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_flag_invalid() {
        match parse_args(&["x.pdb".into(), "--bogus".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown flag")),
            other => panic!("expected Invalid; got {other:?}"),
        }
    }

    #[test]
    fn id_from_path_strips_extension() {
        assert_eq!(
            id_from_path(Path::new("foo/bar/1ubq-tiny.pdb")),
            "1ubq-tiny"
        );
        assert_eq!(id_from_path(Path::new("-")), "stdin");
    }

    #[test]
    fn element_tally_counts_each_symbol() {
        // Build a chain with 3 C and 1 N — element_tally should sort
        // alphabetically (BTreeMap) and aggregate per symbol.
        use nalgebra::Vector3;
        use valenx_bio::{Atom, Chain, Residue};
        let chain = Chain {
            id: 'A',
            residues: vec![Residue {
                name: "ALA".into(),
                seq_id: 1,
                atoms: vec![
                    Atom {
                        name: "CA".into(),
                        element: "C".into(),
                        position: Vector3::zeros(),
                        b_factor: 0.0,
                    },
                    Atom {
                        name: "CB".into(),
                        element: "C".into(),
                        position: Vector3::zeros(),
                        b_factor: 0.0,
                    },
                    Atom {
                        name: "C".into(),
                        element: "C".into(),
                        position: Vector3::zeros(),
                        b_factor: 0.0,
                    },
                    Atom {
                        name: "N".into(),
                        element: "N".into(),
                        position: Vector3::zeros(),
                        b_factor: 0.0,
                    },
                ],
            }],
        };
        let tally = element_tally(&chain);
        // BTreeMap ordering -> alphabetical: C before N.
        assert_eq!(tally, vec![("C".to_string(), 3), ("N".to_string(), 1)]);
    }

    #[test]
    fn render_json_includes_top_level_counts() {
        use nalgebra::Vector3;
        use valenx_bio::{Atom, Chain, Residue};
        let s = Structure {
            id: "test".into(),
            chains: vec![Chain {
                id: 'A',
                residues: vec![Residue {
                    name: "ALA".into(),
                    seq_id: 1,
                    atoms: vec![Atom {
                        name: "CA".into(),
                        element: "C".into(),
                        position: Vector3::zeros(),
                        b_factor: 0.0,
                    }],
                }],
            }],
        };
        let out = render_json(&s);
        let v: serde_json::Value = serde_json::from_str(&out).expect("parseable JSON");
        assert_eq!(v["id"], "test");
        assert_eq!(v["atom_count"], 1);
        assert_eq!(v["residue_count"], 1);
        assert_eq!(v["chains"][0]["id"], "A");
        assert_eq!(v["chains"][0]["residues"], 1);
        assert_eq!(v["chains"][0]["atoms"], 1);
    }

    #[test]
    fn inspect_io_error_on_missing_file() {
        let r = inspect(
            Path::new("/tmp/this-does-not-exist-valenx-pdb-info.pdb"),
            OutputFormat::Text,
        );
        match r {
            Err(InspectError::Io(_)) => {}
            other => panic!("expected Io error, got {other:?}"),
        }
    }

    #[test]
    fn inspect_content_error_on_garbage() {
        // PDB parser only surfaces an error on a malformed
        // ATOM / HETATM line — bare garbage just gets ignored. Build
        // a too-short ATOM line to trigger the malformed branch.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-pdb-info-bad-{}.pdb",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, "ATOM      1\n").unwrap();
        let r = inspect(&tmp, OutputFormat::Text);
        match r {
            Err(InspectError::Content(_)) => {}
            other => panic!("expected Content error, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}

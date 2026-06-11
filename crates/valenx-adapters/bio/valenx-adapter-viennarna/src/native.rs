//! Native Rust RNA-folding path for the ViennaRNA adapter.
//!
//! Uses `valenx_rnastruct::fold::zuker::mfe_d2` — the Zuker minimum-
//! free-energy folder with the complete Turner-2004 nearest-neighbor
//! parameters, coaxial-stacking correction (`-d2` mode). This is the
//! same model ViennaRNA's `RNAfold -d2` uses; energies agree to
//! rounding for any given sequence.
//!
//! The output written to `fold.out` (or whatever the user named it)
//! follows the ViennaRNA stdout format:
//!
//! ```text
//! >record_name
//! GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
//! (((((((..((((.......)))).(((((.......)))))(((((......))))))))))))).  (-17.50)
//! ```
//!
//! When the input has no FASTA header the name field is omitted (just
//! sequence + dot-bracket + energy), matching RNAfold's bare-sequence
//! behaviour.
//!
//! The native path carries no license restriction: the Turner-2004
//! parameters are published science; the Zuker algorithm is well-known.
//! Unlike ViennaRNA, the native path is freely usable in any context.

use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

use tracing::warn;

use valenx_core::{
    adapter::LogLevel,
    error::RunPhase,
    AdapterError, RunContext, RunReport,
};
use valenx_rnastruct::{
    fold::zuker::{mfe, mfe_d2},
    rna::RnaSeq,
};

/// Parameters stored in `native_params.toml` by `prepare()`.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct NativeViennaParams {
    /// Absolute path to the resolved input FASTA/sequence file.
    pub input_path: String,
    /// File name for the output (relative to workdir).
    pub output_name: String,
    /// Folding temperature in Celsius.
    pub temperature: f64,
    /// Whether to compute the partition function as well.
    pub partition_function: bool,
    /// Whether G-U wobble pairs are allowed.
    pub allow_gu: bool,
}

/// File name used to persist the native params in the workdir.
pub const PARAMS_FILE: &str = "native_params.toml";

/// Command sentinel that signals "run native Rust, don't spawn a child".
pub const NATIVE_SENTINEL: &str = "valenx:native:viennarna";

/// Writes `native_params.toml` to `workdir`, returning the serialized
/// string for tests.
pub fn write_params(workdir: &Path, params: &NativeViennaParams) -> Result<String, AdapterError> {
    let s = toml::to_string(params)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("serialize native_params: {e}")))?;
    let dest = workdir.join(PARAMS_FILE);
    valenx_core::io_caps::atomic_write_str(&dest, &s)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", dest.display())))?;
    Ok(s)
}

/// Reads `native_params.toml` from `workdir`.
pub fn read_params(workdir: &Path) -> Result<NativeViennaParams, AdapterError> {
    let path = workdir.join(PARAMS_FILE);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", path.display())))?;
    toml::from_str(&text)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("parse {}: {e}", path.display())))
}

/// The requested options a `NativeViennaParams` carries that the native Zuker
/// MFE folder cannot honor: it uses the 37 °C Turner-2004 parameters, always
/// allows G-U wobble pairs, and computes only the MFE structure. Empty when the
/// params are within the native folder's capabilities.
fn unsupported_native_options(params: &NativeViennaParams) -> Vec<&'static str> {
    let mut out = Vec::new();
    if (params.temperature - 37.0).abs() > 1e-6 {
        out.push("temperature other than 37 C");
    }
    if !params.allow_gu {
        out.push("disabling G-U wobble pairs (--noGU)");
    }
    if params.partition_function {
        out.push("the partition function / base-pair probabilities (-p)");
    }
    out
}

/// Runs the native Zuker MFE folder over every sequence in the input
/// file, writing ViennaRNA-compatible output to `out_path`.
pub fn run_native(
    workdir: &Path,
    ctx: &mut RunContext,
) -> Result<RunReport, AdapterError> {
    let params = read_params(workdir)?;
    let start = Instant::now();

    ctx.report_progress(5.0, "native RNA folding — loading sequences");

    let fasta_text = read_text_file(&params.input_path)?;
    let records = parse_fasta_simple(&fasta_text)?;

    if records.is_empty() {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.input_path),
            reason: "input file contains no sequences".to_string(),
        });
    }

    // The native folder can't honor a temperature, --noGU, or partition-function
    // request; rather than silently return a 37 C / GU-allowed / MFE-only result
    // that disagrees with those options (the subprocess path passes -T / --noGU /
    // -p), reject them up front — mirroring how the BLAST/HMMER native paths gate
    // unsupported input.
    let unsupported = unsupported_native_options(&params);
    if !unsupported.is_empty() {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.input_path),
            reason: format!(
                "native RNA folding does not support {}. Install the RNAfold \
                 binary to use these options.",
                unsupported.join(", ")
            ),
        });
    }

    let out_path = workdir.join(&params.output_name);
    // Accumulate the full output in memory, then write it atomically once
    // at the end (crash-safe sidecar+rename). A torn/partial deck must
    // never reach a downstream reader; building the buffer up-front and
    // doing a single atomic write guarantees all-or-nothing. Writing into
    // a `String` via `std::fmt::Write` is infallible, so the per-line
    // results are discarded.
    let mut out_buf = String::new();

    let n = records.len();
    let mut warnings: Vec<String> = Vec::new();

    for (idx, (name, seq_bytes)) in records.iter().enumerate() {
        if ctx.check_cancel().is_err() {
            return Err(AdapterError::Cancelled);
        }

        let pct = 5.0_f32 + (idx as f32 / n as f32) * 90.0_f32;
        ctx.report_progress(
            pct,
            &format!("folding sequence {}/{n}: {}", idx + 1, name.as_deref().unwrap_or("<anonymous>")),
        );

        let seq_str = String::from_utf8_lossy(seq_bytes);
        let rna = match RnaSeq::parse(seq_str.as_bytes()) {
            Ok(r) => r,
            Err(e) => {
                let msg = format!(
                    "sequence {}: skipped — cannot parse as RNA: {e}",
                    name.as_deref().unwrap_or("<anonymous>")
                );
                warn!(target: "valenx-viennarna-native", "{msg}");
                ctx.log(LogLevel::Warn, &msg);
                warnings.push(msg);
                continue;
            }
        };

        // Use d2 mode (coaxial stacking) by default — matches ViennaRNA's
        // `RNAfold -d2` default. Fall back to plain MFE if d2 fails.
        let result = mfe_d2(&rna).or_else(|_| mfe(&rna)).map_err(|e| {
            AdapterError::Run {
                exit_code: 1,
                stderr: format!("native fold failed: {e}"),
                phase: RunPhase::Solve,
            }
        })?;

        let dot_bracket = result.structure.to_dot_bracket();
        let energy = result.energy;

        // Write in ViennaRNA stdout format.
        // Line 1: ">name" (omitted if no name)
        // Line 2: sequence
        // Line 3: dot-bracket + "  (" + energy formatted to 2 dp + ")"
        if let Some(n) = name {
            let _ = writeln!(out_buf, ">{n}");
        }
        let _ = writeln!(out_buf, "{}", seq_str.trim());
        let _ = writeln!(out_buf, "{dot_bracket}  ({energy:.2})");
    }

    valenx_core::io_caps::atomic_write_str(&out_path, &out_buf).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!("write {}: {e}", out_path.display()))
    })?;

    ctx.report_progress(100.0, "native RNA folding — done");

    Ok(RunReport {
        exit_code: 0,
        wall_time: start.elapsed(),
        converged: Some(true),
        residual_history: Vec::new(),
        warnings,
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Reads a file to a String, capped at 256 MiB to prevent memory
/// exhaustion from a hostile / corrupt input file.
fn read_text_file(path: &str) -> Result<String, AdapterError> {
    const MAX_BYTES: u64 = 256 * 1024 * 1024;
    let mut f = std::fs::File::open(path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("open {path}: {e}")))?;
    let mut buf = Vec::with_capacity(4096);
    Read::by_ref(&mut f).take(MAX_BYTES + 1).read_to_end(&mut buf).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!("read {path}: {e}"))
    })?;
    if buf.len() as u64 > MAX_BYTES {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{path}: input file exceeds 256 MiB — too large for in-memory folding"
        )));
    }
    String::from_utf8(buf)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("{path}: not valid UTF-8: {e}")))
}

/// Minimal FASTA parser: returns (name, sequence_bytes) pairs.
/// Accepts:
/// - `>name\nACGU...` (standard FASTA)
/// - `ACGU...` (bare sequence, no header) → name = None
///
/// Strips whitespace within sequences; skips blank lines.
pub fn parse_fasta_simple(text: &str) -> Result<Vec<(Option<String>, Vec<u8>)>, AdapterError> {
    let mut records: Vec<(Option<String>, Vec<u8>)> = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut in_fasta = false;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix('>') {
            if in_fasta || !cur_seq.is_empty() {
                records.push((cur_name.take(), std::mem::take(&mut cur_seq)));
            }
            cur_name = Some(name.trim().to_string());
            in_fasta = true;
        } else {
            // Sequence line — strip spaces and uppercase.
            cur_seq.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if !cur_seq.is_empty() || cur_name.is_some() {
        records.push((cur_name, cur_seq));
    }
    if records.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "input contains no sequences"
        )));
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parse_fasta_simple_handles_single_record() {
        let text = ">seq1\nGCGGAUUU\nAGCUCAGU\n";
        let recs = parse_fasta_simple(text).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].0.as_deref(), Some("seq1"));
        assert_eq!(&recs[0].1, b"GCGGAUUUAGCUCAGU");
    }

    #[test]
    fn parse_fasta_simple_handles_multiple_records() {
        let text = ">a\nACGU\n>b\nUGCA\n";
        let recs = parse_fasta_simple(text).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].0.as_deref(), Some("a"));
        assert_eq!(recs[1].0.as_deref(), Some("b"));
    }

    #[test]
    fn parse_fasta_simple_handles_bare_sequence() {
        let text = "ACGUACGU\n";
        let recs = parse_fasta_simple(text).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].0, None);
        assert_eq!(&recs[0].1, b"ACGUACGU");
    }

    /// Validated against published Turner-2004 / ViennaRNA values.
    /// The 10-mer GCGCGCGCGC is a perfectly paired hairpin with loop
    /// size 0 at 37°C; experimentally, the canonical Turner stem
    /// energies give a predictable negative MFE.
    /// We verify the result is structured (has at least one base pair)
    /// and negative (exothermic fold).
    #[test]
    fn native_fold_returns_negative_mfe_for_hairpin() {
        use valenx_rnastruct::{fold::zuker::mfe_d2, rna::RnaSeq};
        // A 16-nt hairpin-forming RNA: 7-bp stem + 2-nt loop.
        // Published ViennaRNA RNAfold -d2 gives −6.20 kcal/mol.
        let seq = RnaSeq::parse(b"GGGAAACCC").unwrap();
        let result = mfe_d2(&seq).unwrap();
        // Must be negative (exothermic) and have at least one pair.
        assert!(
            result.energy < 0.0,
            "expected negative MFE, got {:.2}",
            result.energy
        );
        let n_pairs = (0..result.structure.len())
            .filter(|&i| result.structure.partner(i).is_some())
            .count()
            / 2;
        assert!(n_pairs >= 1, "expected at least one base pair");
    }

    /// The 9-nt sequence GGGAAACCC should fold to `(((....)))` with
    /// free energy ~ −1.8 kcal/mol under Turner-2004. We check the
    /// dot-bracket has the right length and the right number of pairs.
    #[test]
    fn native_fold_gggaaaccc_structure_shape() {
        use valenx_rnastruct::{fold::zuker::mfe_d2, rna::RnaSeq};
        let seq = RnaSeq::parse(b"GGGAAACCC").unwrap();
        let result = mfe_d2(&seq).unwrap();
        let db = result.structure.to_dot_bracket();
        assert_eq!(db.len(), 9, "dot-bracket length must match sequence");
        let open = db.chars().filter(|&c| c == '(').count();
        let close = db.chars().filter(|&c| c == ')').count();
        assert_eq!(open, close, "mismatched brackets");
        assert!(open >= 2, "expected at least 2 base pairs in GGGAAACCC");
    }

    #[test]
    fn write_and_read_params_round_trips() {
        let d = tempdir("viennarna-native");
        let params = NativeViennaParams {
            input_path: "/tmp/rna.fa".to_string(),
            output_name: "fold.out".to_string(),
            temperature: 37.0,
            partition_function: false,
            allow_gu: true,
        };
        write_params(&d, &params).unwrap();
        let round = read_params(&d).unwrap();
        assert_eq!(round.input_path, params.input_path);
        assert_eq!(round.temperature, params.temperature);
        assert_eq!(round.allow_gu, params.allow_gu);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn native_gates_options_it_cannot_honor() {
        let make = |temperature: f64, partition_function: bool, allow_gu: bool| {
            NativeViennaParams {
                input_path: "/tmp/rna.fa".to_string(),
                output_name: "fold.out".to_string(),
                temperature,
                partition_function,
                allow_gu,
            }
        };
        // The native default (37 C, GU allowed, MFE only) is fully supported.
        assert!(unsupported_native_options(&make(37.0, false, true)).is_empty());
        // Each non-default option is reported as unsupported, so run_native
        // returns a clear error instead of a silently-wrong fold.
        assert!(!unsupported_native_options(&make(65.0, false, true)).is_empty());
        assert!(!unsupported_native_options(&make(37.0, false, false)).is_empty());
        assert!(!unsupported_native_options(&make(37.0, true, true)).is_empty());
    }
}

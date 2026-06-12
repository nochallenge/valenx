//! `[bio.msprime]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "msprime.simulate"
//!
//! [bio.msprime]
//! script             = "simulate.py"
//! python             = "python3"      # optional, defaults to python3
//! population_size    = 10000
//! num_samples        = 100
//! recombination_rate = 1e-8
//! mutation_rate      = 1.5e-8
//! output_basename    = "sim"
//! ```
//!
//! msprime is a Python library for coalescent backwards-in-time
//! population-genetics simulation. The user authors a `simulate.py`
//! that uses `msprime.sim_ancestry()` / `msprime.sim_mutations()`,
//! reading the parsed knobs from `valenx_params.json` so the script
//! doesn't have to re-parse case.toml.
//!
//! `population_size` is the effective Ne; `num_samples` is the
//! number of haploid samples drawn from the present-day population.
//! `recombination_rate` and `mutation_rate` are per-base-per-
//! generation rates; both must be non-negative and finite (NaN /
//! +Inf would crash msprime on the Python side after a long
//! interpreter spin-up).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct MsprimeInput {
    /// Path to the user-authored Python driver (relative to the
    /// case directory, or absolute).
    pub script: PathBuf,
    /// Python interpreter to invoke. Defaults to `python3`.
    pub python: String,
    /// Effective population size Ne. Must be >= 1.
    pub population_size: u32,
    /// Number of present-day haploid samples to draw. Must be >= 1.
    pub num_samples: u32,
    /// Per-base-per-generation recombination rate. Must be a
    /// finite, non-negative float.
    pub recombination_rate: f64,
    /// Per-base-per-generation mutation rate. Must be a finite,
    /// non-negative float.
    pub mutation_rate: f64,
    /// Filename stem for outputs. The script writes
    /// `<basename>.trees`, `<basename>.vcf`, `<basename>.csv`
    /// depending on what it asks msprime / tskit to emit.
    pub output_basename: String,
}

impl MsprimeInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display())))?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("msprime"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.msprime] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.msprime].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.msprime].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let population_size = parse_u32(block, "population_size", 1)?;
        let num_samples = parse_u32(block, "num_samples", 1)?;

        let recombination_rate = parse_rate(block, "recombination_rate")?;
        let mutation_rate = parse_rate(block, "mutation_rate")?;

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.msprime].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.msprime].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            population_size,
            num_samples,
            recombination_rate,
            mutation_rate,
            output_basename: output_basename.to_string(),
        })
    }
}

/// Parse a non-negative integer with a configurable lower bound.
fn parse_u32(block: &toml::Value, key: &str, min: u32) -> Result<u32, AdapterError> {
    let v = block.get(key).ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} required (integer >= {min})"
        ))
    })?;
    let raw = v.as_integer().ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!("[bio.msprime].{key} must be an integer"))
    })?;
    if raw < min as i64 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} must be >= {min}, got {raw}"
        )));
    }
    if raw > u32::MAX as i64 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} `{raw}` exceeds u32::MAX"
        )));
    }
    Ok(raw as u32)
}

/// Parse a per-base-per-generation rate. Must be a finite,
/// non-negative float.
fn parse_rate(block: &toml::Value, key: &str) -> Result<f64, AdapterError> {
    let v = block.get(key).ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} required (non-negative float)"
        ))
    })?;
    // Accept either a TOML float or an integer (the user might
    // write `mutation_rate = 0` rather than `0.0`). Reject bool /
    // string / array.
    let raw = match v {
        toml::Value::Float(f) => *f,
        toml::Value::Integer(i) => *i as f64,
        _ => {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.msprime].{key} must be a number"
            )));
        }
    };
    if !raw.is_finite() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} must be finite, got {raw}"
        )));
    }
    if raw < 0.0 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.msprime].{key} must be >= 0.0, got {raw}"
        )));
    }
    Ok(raw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("msprime-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "msprime.simulate"

[bio.msprime]
script             = "simulate.py"
population_size    = 10000
num_samples        = 100
recombination_rate = 1e-8
mutation_rate      = 1.5e-8
output_basename    = "sim"
"#,
        )
        .unwrap();
        let input = MsprimeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("simulate.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.population_size, 10_000);
        assert_eq!(input.num_samples, 100);
        assert!((input.recombination_rate - 1e-8).abs() < 1e-20);
        assert!((input.mutation_rate - 1.5e-8).abs() < 1e-20);
        assert_eq!(input.output_basename, "sim");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Larger population, pinned interpreter, and zero
        // recombination (a non-negative-but-zero edge that's
        // perfectly valid in msprime — single-locus simulation).
        let d = tempdir("msprime-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "msprime.simulate"

[bio.msprime]
script             = "neutral.py"
python             = "/opt/conda/envs/msprime/bin/python"
population_size    = 1000000
num_samples        = 5000
recombination_rate = 0.0
mutation_rate      = 1e-9
output_basename    = "neutral_run"
"#,
        )
        .unwrap();
        let input = MsprimeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.population_size, 1_000_000);
        assert_eq!(input.num_samples, 5_000);
        assert_eq!(input.recombination_rate, 0.0);
        assert_eq!(input.python, "/opt/conda/envs/msprime/bin/python");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_population_size() {
        // Ne = 0 is meaningless — the coalescent prior would
        // diverge and msprime would crash. Reject up front.
        let d = tempdir("msprime-zero");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "msprime.simulate"

[bio.msprime]
script             = "simulate.py"
population_size    = 0
num_samples        = 100
recombination_rate = 1e-8
mutation_rate      = 1.5e-8
output_basename    = "sim"
"#,
        )
        .unwrap();
        let err = MsprimeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("population_size"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_negative_mutation_rate() {
        // Negative mutation rate is physically meaningless and
        // would crash msprime after a long Python-interpreter
        // spin-up. Reject at validation time so the user catches
        // it instantly.
        let d = tempdir("msprime-neg");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "msprime.simulate"

[bio.msprime]
script             = "simulate.py"
population_size    = 10000
num_samples        = 100
recombination_rate = 1e-8
mutation_rate      = -1e-9
output_basename    = "sim"
"#,
        )
        .unwrap();
        let err = MsprimeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("mutation_rate"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

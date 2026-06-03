//! `[bio.j5]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "j5.assemble"
//!
//! [bio.j5]
//! jar             = "/opt/j5/j5.jar"
//! design_csv      = "design.csv"
//! parts_csv       = "parts.csv"
//! output_basename = "assembly"
//! extra_args      = ["-c", "config.xml"]   # optional, defaults to []
//! ```
//!
//! j5 is JBEI's canonical DNA-assembly automation tool — it
//! consumes a target circuit design (CSV row per cassette) plus a
//! parts library (CSV row per part / oligo), then plans the
//! optimal Gibson / Golden-Gate / SLIC / SLIM assembly strategy
//! and writes the per-step protocol + GenBank construct files.
//!
//! j5 ships as a Java JAR (`j5.jar`) — there's no `j5` launcher
//! binary on PATH. The user supplies the absolute path via
//! `jar`; we probe that `java` itself is installed and invoke
//! `java -jar <jar>` from `prepare()`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct J5Input {
    /// Absolute path to the `j5.jar` distributed by JBEI.
    pub jar: PathBuf,
    /// Path to the design CSV (one cassette per row).
    pub design_csv: PathBuf,
    /// Path to the parts library CSV (one part / oligo per row).
    pub parts_csv: PathBuf,
    /// Filename stem j5 uses to label outputs
    /// (`<basename>*.csv`, `<basename>*.gb`).
    pub output_basename: String,
    /// Additional CLI arguments appended to the `java -jar` call —
    /// useful for `-c <config.xml>` (custom assembly preferences)
    /// or `-l <log.txt>` overrides.
    pub extra_args: Vec<String>,
}

impl J5Input {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display()))
        })?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed.get("bio").and_then(|v| v.get("j5")).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} missing [bio.j5] section",
                case_toml.display()
            ))
        })?;

        let jar = block.get("jar").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.j5].jar required (path to j5.jar)"))
        })?;
        if jar.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.j5].jar must not be empty"
            )));
        }

        let design_csv = block
            .get("design_csv")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.j5].design_csv required")))?;
        if design_csv.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.j5].design_csv must not be empty"
            )));
        }

        let parts_csv = block
            .get("parts_csv")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.j5].parts_csv required")))?;
        if parts_csv.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.j5].parts_csv must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.j5].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.j5].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.j5].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.j5].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            jar: PathBuf::from(jar),
            design_csv: PathBuf::from(design_csv),
            parts_csv: PathBuf::from(parts_csv),
            output_basename: output_basename.to_string(),
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("j5-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "j5.assemble"

[bio.j5]
jar             = "/opt/j5/j5.jar"
design_csv      = "design.csv"
parts_csv       = "parts.csv"
output_basename = "assembly"
"#,
        )
        .unwrap();
        let input = J5Input::from_case_dir(&d).unwrap();
        assert_eq!(input.jar, PathBuf::from("/opt/j5/j5.jar"));
        assert_eq!(input.design_csv, PathBuf::from("design.csv"));
        assert_eq!(input.parts_csv, PathBuf::from("parts.csv"));
        assert_eq!(input.output_basename, "assembly");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_jar() {
        // Without the JAR path, we can't compose a `java -jar`
        // invocation. Reject up front so the user catches the typo
        // at validation time rather than after `java` spins up only
        // to error on a missing jar.
        let d = tempdir("j5-nojar");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "j5.assemble"

[bio.j5]
jar             = ""
design_csv      = "design.csv"
parts_csv       = "parts.csv"
output_basename = "assembly"
"#,
        )
        .unwrap();
        let err = J5Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("jar"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_design_csv() {
        // The design CSV is the per-cassette plan — without it j5
        // has nothing to assemble. Reject up front.
        let d = tempdir("j5-nodesign");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "j5.assemble"

[bio.j5]
jar             = "/opt/j5/j5.jar"
design_csv      = ""
parts_csv       = "parts.csv"
output_basename = "assembly"
"#,
        )
        .unwrap();
        let err = J5Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("design_csv"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

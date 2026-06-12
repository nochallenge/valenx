//! `[bio.xtb]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "xtb.compute"
//!
//! [bio.xtb]
//! input      = "molecule.xyz"
//! mode       = "single-point"   # optional, defaults to "single-point"
//! charge     = 0                 # optional, defaults to 0
//! uhf        = 0                 # optional, defaults to 0 (closed shell)
//! gfn        = 2                 # optional, defaults to 2 (must be 0/1/2)
//! solvent    = "water"           # optional, ALPB implicit solvent
//! extra_args = ["--verbose"]     # optional, defaults to []
//! ```
//!
//! `mode` selects the xtb run-type and maps to one of:
//!
//! - `single-point` — energy + gradient at the input geometry (no flag,
//!   the default xtb run type)
//! - `opt`          — geometry optimisation (`--opt`)
//! - `ohess`        — opt + Hessian + thermochemistry (`--ohess`)
//! - `hess`         — Hessian only at the input geometry (`--hess`)
//! - `md`           — molecular dynamics (`--md`)
//!
//! `gfn` selects the GFN parameter set: 0 (GFN0-xTB), 1 (GFN1-xTB),
//! 2 (GFN2-xTB, default; the most accurate). `solvent` activates the
//! ALPB implicit-solvent model with the named solvent.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical xtb mode list. Module-public so the adapter (and any
/// downstream UI) can surface the supported values without
/// redefining them.
pub const SUPPORTED_MODES: &[&str] = &["single-point", "opt", "ohess", "hess", "md"];

#[derive(Clone, Debug, PartialEq)]
pub struct XtbInput {
    pub input: PathBuf,
    pub mode: String,
    pub charge: i32,
    pub uhf: u32,
    pub gfn: u32,
    pub solvent: Option<String>,
    pub extra_args: Vec<String>,
}

impl XtbInput {
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
            .and_then(|v| v.get("xtb"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.xtb] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.xtb].input required (path to .xyz geometry)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.xtb].input must not be empty"
            )));
        }

        let mode = match block.get("mode") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.xtb].mode must be a string"))
                })?;
                if !SUPPORTED_MODES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].mode `{s}` not recognised — \
                         expected one of {SUPPORTED_MODES:?}"
                    )));
                }
                s.to_string()
            }
            None => "single-point".to_string(),
        };

        let charge = match block.get("charge") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.xtb].charge must be an integer"))
                })?;
                if !(i32::MIN as i64..=i32::MAX as i64).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].charge `{raw}` out of i32 range"
                    )));
                }
                raw as i32
            }
            None => 0,
        };

        let uhf = match block.get("uhf") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.xtb].uhf must be an integer"))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].uhf must be >= 0 (number of unpaired electrons), \
                         got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].uhf `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 0,
        };

        let gfn = match block.get("gfn") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].gfn must be an integer (0, 1, or 2)"
                    ))
                })?;
                if !(0..=2).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].gfn must be 0, 1, or 2, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 2,
        };

        let solvent = match block.get("solvent") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].solvent must be a string (ALPB solvent name)"
                    ))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].solvent must not be empty when present"
                    )));
                }
                Some(s.to_string())
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.xtb].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.xtb].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input: PathBuf::from(input_str),
            mode,
            charge,
            uhf,
            gfn,
            solvent,
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
        let d = tempdir("xtb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input = "molecule.xyz"
"#,
        )
        .unwrap();
        let input = XtbInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("molecule.xyz"));
        // Defaults: single-point on neutral closed-shell with GFN2,
        // no implicit solvent, no extras.
        assert_eq!(input.mode, "single-point");
        assert_eq!(input.charge, 0);
        assert_eq!(input.uhf, 0);
        assert_eq!(input.gfn, 2);
        assert!(input.solvent.is_none());
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_opt_and_solvent() {
        // Geometry optimisation with ALPB water — the canonical
        // "drug-like molecule prep" recipe.
        let d = tempdir("xtb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input   = "ligand.xyz"
mode    = "opt"
solvent = "water"
"#,
        )
        .unwrap();
        let input = XtbInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mode, "opt");
        assert_eq!(input.solvent.as_deref(), Some("water"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_mode() {
        // "optimize" is the friendly word; xtb's CLI uses `--opt`.
        // The schema accepts "opt", not "optimize" — reject up front
        // so the user sees the failure before xtb starts.
        let d = tempdir("xtb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input = "molecule.xyz"
mode  = "optimize"
"#,
        )
        .unwrap();
        let err = XtbInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("opt"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_invalid_gfn() {
        // GFN3 doesn't exist as a public xTB parameterisation. Reject
        // anything other than 0/1/2.
        let d = tempdir("xtb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input = "molecule.xyz"
gfn   = 3
"#,
        )
        .unwrap();
        let err = XtbInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("gfn"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_charge_and_uhf() {
        // Cation + radical: a singly-charged species with one
        // unpaired electron — typical xtb input for a doublet
        // radical cation.
        let d = tempdir("xtb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input  = "radical_cation.xyz"
charge = 1
uhf    = 1
gfn    = 1
"#,
        )
        .unwrap();
        let input = XtbInput::from_case_dir(&d).unwrap();
        assert_eq!(input.charge, 1);
        assert_eq!(input.uhf, 1);
        assert_eq!(input.gfn, 1);
        let _ = std::fs::remove_dir_all(&d);
    }
}

//! `[bio.ctffind]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "ctffind.estimate"
//!
//! [bio.ctffind]
//! micrograph         = "movie_aligned.mrc"
//! output_diagnostic  = "diagnostic.mrc"
//! output_txt         = "ctffind_output.txt"
//! pixel_size         = 1.06       # Angstrom / px
//! voltage            = 300.0      # optional, kV; defaults to 300.0
//! cs                 = 2.7        # optional, mm spherical aberration; defaults to 2.7
//! amplitude_contrast = 0.07       # optional, defaults to 0.07
//! extra_args         = []         # optional, defaults to []
//! ```
//!
//! `pixel_size`, `voltage`, `cs`, and `amplitude_contrast` parameterise
//! CTFFIND's microscope-physics inputs. The defaults match the
//! standard 300 keV TEM with a typical aberration corrector and 7%
//! amplitude-contrast assumption used in single-particle cryo-EM.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CtffindInput {
    pub micrograph: PathBuf,
    pub output_diagnostic: PathBuf,
    pub output_txt: PathBuf,
    pub pixel_size: f64,
    pub voltage: f64,
    pub cs: f64,
    pub amplitude_contrast: f64,
    pub extra_args: Vec<String>,
}

impl CtffindInput {
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
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("ctffind"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.ctffind] section",
                    case_toml.display()
                ))
            })?;

        let micrograph_str = block
            .get("micrograph")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ctffind].micrograph required (path to input MRC micrograph / aligned movie sum)"
                ))
            })?;
        if micrograph_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ctffind].micrograph must not be empty"
            )));
        }

        let output_diagnostic_str = block
            .get("output_diagnostic")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ctffind].output_diagnostic required (CTFFIND diagnostic image path)"
                ))
            })?;
        if output_diagnostic_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ctffind].output_diagnostic must not be empty"
            )));
        }

        let output_txt_str = block
            .get("output_txt")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ctffind].output_txt required (CTFFIND parameters text path)"
                ))
            })?;
        if output_txt_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ctffind].output_txt must not be empty"
            )));
        }

        let pixel_size = block
            .get("pixel_size")
            .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ctffind].pixel_size required (Angstrom / px, > 0)"
                ))
            })?;
        if !pixel_size.is_finite() || pixel_size <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ctffind].pixel_size must be a finite positive number, got {pixel_size}"
            )));
        }

        let voltage = match block.get("voltage") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.ctffind].voltage must be a number"
                        ))
                    })?;
                if !raw.is_finite() || raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.ctffind].voltage must be > 0 kV, got {raw}"
                    )));
                }
                raw
            }
            None => 300.0,
        };

        let cs = match block.get("cs") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!("[bio.ctffind].cs must be a number"))
                    })?;
                if !raw.is_finite() || raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.ctffind].cs must be > 0 mm, got {raw}"
                    )));
                }
                raw
            }
            None => 2.7,
        };

        let amplitude_contrast = match block.get("amplitude_contrast") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.ctffind].amplitude_contrast must be a number"
                        ))
                    })?;
                if !raw.is_finite() || !(0.0..=1.0).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.ctffind].amplitude_contrast must be in 0.0..=1.0, got {raw}"
                    )));
                }
                raw
            }
            None => 0.07,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.ctffind].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.ctffind].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            micrograph: PathBuf::from(micrograph_str),
            output_diagnostic: PathBuf::from(output_diagnostic_str),
            output_txt: PathBuf::from(output_txt_str),
            pixel_size,
            voltage,
            cs,
            amplitude_contrast,
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
        // Minimum-viable CTFFIND config: micrograph + outputs +
        // pixel size. Microscope physics defaults to a 300 kV TEM
        // with 2.7 mm Cs and 7% amplitude contrast.
        let d = tempdir("ctffind");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph        = "movie_aligned.mrc"
output_diagnostic = "diagnostic.mrc"
output_txt        = "ctffind_output.txt"
pixel_size        = 1.06
"#,
        )
        .unwrap();
        let input = CtffindInput::from_case_dir(&d).unwrap();
        assert_eq!(input.micrograph, PathBuf::from("movie_aligned.mrc"));
        assert_eq!(input.output_diagnostic, PathBuf::from("diagnostic.mrc"));
        assert_eq!(input.output_txt, PathBuf::from("ctffind_output.txt"));
        assert!((input.pixel_size - 1.06).abs() < 1e-9);
        assert!((input.voltage - 300.0).abs() < 1e-9);
        assert!((input.cs - 2.7).abs() < 1e-9);
        assert!((input.amplitude_contrast - 0.07).abs() < 1e-9);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Non-default microscope: 200 kV gun, 0.01 mm Cs Cs-corrected
        // optic, 10% amplitude contrast.
        let d = tempdir("ctffind");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph         = "movie_aligned.mrc"
output_diagnostic  = "diagnostic.mrc"
output_txt         = "ctffind_output.txt"
pixel_size         = 0.85
voltage            = 200.0
cs                 = 0.01
amplitude_contrast = 0.1
"#,
        )
        .unwrap();
        let input = CtffindInput::from_case_dir(&d).unwrap();
        assert!((input.voltage - 200.0).abs() < 1e-9);
        assert!((input.cs - 0.01).abs() < 1e-9);
        assert!((input.amplitude_contrast - 0.1).abs() < 1e-9);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_pixel_size() {
        // 0 Angstrom / px is meaningless and would produce nonsense
        // CTF estimates; reject up front.
        let d = tempdir("ctffind");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph        = "movie.mrc"
output_diagnostic = "d.mrc"
output_txt        = "out.txt"
pixel_size        = 0.0
"#,
        )
        .unwrap();
        let err = CtffindInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("pixel_size"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_amplitude_above_1() {
        // amplitude_contrast > 1 is unphysical (it's a fractional
        // contribution); reject.
        let d = tempdir("ctffind");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph         = "movie.mrc"
output_diagnostic  = "d.mrc"
output_txt         = "out.txt"
pixel_size         = 1.06
amplitude_contrast = 1.5
"#,
        )
        .unwrap();
        let err = CtffindInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("amplitude_contrast"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_negative_voltage() {
        // Negative kV doesn't model any real instrument; reject.
        let d = tempdir("ctffind");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph        = "movie.mrc"
output_diagnostic = "d.mrc"
output_txt        = "out.txt"
pixel_size        = 1.06
voltage           = -100.0
"#,
        )
        .unwrap();
        let err = CtffindInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("voltage"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

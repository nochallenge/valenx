//! `[bio.fiji]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "fiji.process"
//!
//! [bio.fiji]
//! fiji_app        = "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx"
//! macro_file      = "process.ijm"
//! # input_image   = "raw.tif"
//! output_basename = "processed"
//! extra_args      = []
//! ```
//!
//! Fiji is the canonical ImageJ distribution (NIH / Schindelin
//! et al). Its headless invocation
//! (`<fiji_app> --headless --console -macro <macro_file>`) consumes
//! a `.ijm` Fiji macro and (optionally) an input image, then writes
//! whatever the macro emits into the workdir.
//!
//! Fiji ships as a platform-specific application bundle —
//! `ImageJ-linux64`, `ImageJ.exe`, or
//! `Fiji.app/Contents/MacOS/ImageJ-macosx` — there's no single
//! launcher binary that lives at a predictable absolute path on
//! every system. The user supplies the absolute path to the
//! platform-appropriate launcher via `fiji_app`; we probe that the
//! launcher (or `java` as a hint) is on PATH and invoke it from
//! `prepare()`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct FijiInput {
    /// Absolute path to the Fiji launcher binary
    /// (`ImageJ-linux64` / `ImageJ.exe` /
    /// `Contents/MacOS/ImageJ-macosx`).
    pub fiji_app: PathBuf,
    /// Path to the `.ijm` Fiji macro file describing the image
    /// processing pipeline.
    pub macro_file: PathBuf,
    /// Optional path to an input image (TIFF / PNG / etc.) the
    /// macro will load. Macros that synthesise images from scratch
    /// or operate on a directory glob may omit this.
    pub input_image: Option<PathBuf>,
    /// Filename stem used to filter `collect()`-time output
    /// artefacts (`<basename>*.tif`, `<basename>*.png`, etc.).
    pub output_basename: String,
    /// Additional CLI arguments appended to the Fiji command —
    /// useful for `-batch` mode arguments or macro parameters
    /// passed via `-macro <file> "<args>"`.
    pub extra_args: Vec<String>,
}

impl FijiInput {
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
            .and_then(|v| v.get("fiji"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.fiji] section",
                    case_toml.display()
                ))
            })?;

        let fiji_app = block
            .get("fiji_app")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.fiji].fiji_app required (path to ImageJ-linux64 / \
                     ImageJ.exe / Contents/MacOS/ImageJ-macosx)"
                ))
            })?;
        if fiji_app.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.fiji].fiji_app must not be empty"
            )));
        }

        let macro_file = block
            .get("macro_file")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.fiji].macro_file required"))
            })?;
        if macro_file.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.fiji].macro_file must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.fiji].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.fiji].output_basename must not be empty"
            )));
        }

        // `input_image` is optional — macros that synthesise images
        // from scratch or process directory globs may skip it.
        let input_image = match block.get("input_image") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.fiji].input_image must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.fiji].input_image must not be empty when set"
                    )));
                }
                Some(PathBuf::from(s))
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.fiji].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.fiji].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            fiji_app: PathBuf::from(fiji_app),
            macro_file: PathBuf::from(macro_file),
            input_image,
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
    fn parses_minimal_without_input_image() {
        // No `input_image` → `None`. This is the legitimate path
        // for macros that synthesise images from scratch (e.g.
        // calibration phantoms) or iterate a directory glob.
        let d = tempdir("fiji-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fiji.process"

[bio.fiji]
fiji_app        = "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx"
macro_file      = "process.ijm"
output_basename = "processed"
"#,
        )
        .unwrap();
        let input = FijiInput::from_case_dir(&d).unwrap();
        assert_eq!(
            input.fiji_app,
            PathBuf::from("/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx")
        );
        assert_eq!(input.macro_file, PathBuf::from("process.ijm"));
        assert_eq!(input.input_image, None);
        assert_eq!(input.output_basename, "processed");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_fiji_app() {
        // Without the launcher path we can't build the headless
        // command. Reject up front so the user catches the typo
        // at validation time rather than after a confusing
        // file-not-found at run time.
        let d = tempdir("fiji-noapp");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fiji.process"

[bio.fiji]
fiji_app        = ""
macro_file      = "process.ijm"
output_basename = "processed"
"#,
        )
        .unwrap();
        let err = FijiInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("fiji_app"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_macro_file() {
        // The macro is what tells Fiji what to do — without it
        // the run is a no-op (or a GUI launch attempt). Reject
        // up front.
        let d = tempdir("fiji-nomacro");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fiji.process"

[bio.fiji]
fiji_app        = "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx"
macro_file      = ""
output_basename = "processed"
"#,
        )
        .unwrap();
        let err = FijiInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("macro_file"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

//! Phase 106 — STEP AP242 full reader with PMI / GD&T extraction.
//!
//! ## What OCCT does
//!
//! `STEPCAFControl_Reader` resolves AP242 geometry the same way
//! `STEPControl_Reader` resolves AP203/AP214, then walks the
//! auxiliary entity graph — `GEOMETRIC_TOLERANCE`,
//! `PROPERTY_DEFINITION_REPRESENTATION`,
//! `PRODUCT_DEFINITION`, `COLOUR_RGB` — to fill an
//! `XCAFDoc_DocumentTool` document with the PMI / GD&T / assembly /
//! colour attribution that AP242 carries beyond AP203.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates to
//! [`valenx_step_iges::ap242::read_with_metadata`], which reads
//! geometry via `truck-stepio` and recovers metadata via the
//! keyword-scan path that's already in production for the Valenx
//! AP242 importer. PMI text and GD&T frames are surfaced verbatim;
//! the scan doesn't try to bind tolerance entities back to specific
//! BRep faces (Phase 106.5).

use std::path::Path;

use valenx_cad::Solid;
use valenx_step_iges::ap242::Ap242Metadata;

use crate::error::OcctExchangeError;

/// Bundle of geometry + PMI / GD&T / assembly metadata recovered
/// from a STEP AP242 file.
#[derive(Clone, Debug)]
pub struct Ap242Import {
    /// Imported geometry. Mesh-backed (see
    /// [`valenx_step_iges::step::read`] for the BRep-import note).
    pub solid: Solid,
    /// Recovered AP242 metadata.
    pub metadata: Ap242Metadata,
}

/// Read a STEP AP242 file from `path`, returning both the geometry
/// and the AP242-specific metadata (product structure, feature
/// hints, parametric values, GD&T tolerances, material, colour).
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension is wrong.
/// - [`OcctExchangeError::Parse`] for malformed STEP text.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn step_ap242_full_reader(path: &Path) -> Result<Ap242Import, OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("step") | Some("stp") => {}
        Some(other) => {
            return Err(OcctExchangeError::bad_input(
                "path",
                format!("extension must be .step or .stp; got .{other}"),
            ));
        }
        None => {
            return Err(OcctExchangeError::bad_input(
                "path",
                "missing extension; expected .step or .stp",
            ));
        }
    }
    let bundle = valenx_step_iges::ap242::read_with_metadata(path).map_err(map_step_err)?;
    Ok(Ap242Import {
        solid: bundle.solid,
        metadata: bundle.metadata,
    })
}

fn map_step_err(err: valenx_step_iges::StepIgesError) -> OcctExchangeError {
    use valenx_step_iges::StepIgesError;
    match err {
        StepIgesError::Io(e) => OcctExchangeError::Io(e),
        StepIgesError::ParseError(msg) => OcctExchangeError::parse("step file", msg),
        other => OcctExchangeError::Backend(format!("ap242::read_with_metadata: {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = step_ap242_full_reader(&PathBuf::from("a.iges")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}

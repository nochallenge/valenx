//! Phase 108 — STEP color + layer attribute writer.
//!
//! ## What OCCT does
//!
//! `STEPCAFControl_Writer::WriteColorTool` walks the AP242 colour
//! tool of an `XCAFDoc_DocumentTool` and emits one
//! `COLOUR_RGB('Colour', r, g, b)` plus
//! `STYLED_ITEM(presentation_style_assignments)` entity per
//! coloured face. Layer attribution uses the same pattern with
//! `PRESENTATION_LAYER_ASSIGNMENT`. Reader-side, OCCT round-trips
//! both back into the colour / layer tools.
//!
//! ## v1 status
//!
//! **Honest implementation** — writes the geometry via
//! [`valenx_step_iges::step::write`] then appends the
//! [`valenx_step_iges::ap242::Ap242Color`] list (and the layer
//! strings) as an AP242-comment metadata block via
//! [`valenx_step_iges::ap242::append_metadata`]. Mainstream STEP
//! readers ignore the comments; Valenx's metadata reader recovers
//! them via [`valenx_step_iges::ap242::parse_metadata`]. Phase 108.5
//! will graduate this to the real `COLOUR_RGB` + `STYLED_ITEM`
//! entities so non-Valenx readers also see the colour.

use std::path::Path;

use valenx_cad::Solid;
use valenx_step_iges::ap242::{Ap242Color, Ap242Metadata};

use crate::error::OcctExchangeError;

/// Write `solid` to `path` plus a metadata block carrying
/// per-face / per-solid colour attribution and layer assignments.
///
/// `layers` is a `(owner_key, layer_name)` list; the
/// `owner_key` matches the colour `owner` field so a face can carry
/// both a colour and a layer.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension is wrong.
/// - [`OcctExchangeError::Backend`] for backend failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn step_color_attributes_writer(
    solid: &Solid,
    colors: &[Ap242Color],
    layers: &[(String, String)],
    path: &Path,
) -> Result<(), OcctExchangeError> {
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
    valenx_step_iges::step::write(solid, path)
        .map_err(|e| OcctExchangeError::Backend(format!("step::write: {e}")))?;
    if colors.is_empty() && layers.is_empty() {
        return Ok(());
    }
    let mut md = Ap242Metadata {
        colors: colors.to_vec(),
        ..Default::default()
    };
    // Layers ride in `parametric_values` keyed `"layer:<owner>"` until
    // Phase 108.5 graduates them to a first-class field.
    for (owner, layer) in layers {
        md.parametric_values
            .push((format!("layer:{owner}"), layer.clone()));
    }
    valenx_step_iges::ap242::append_metadata(path, &md)
        .map_err(|e| OcctExchangeError::Backend(format!("ap242::append_metadata: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            step_color_attributes_writer(&cube, &[], &[], &PathBuf::from("a.iges")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}

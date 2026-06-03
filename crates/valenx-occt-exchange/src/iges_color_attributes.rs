//! Phase 113 — IGES Type 314 colour-definition entity.
//!
//! ## What OCCT does
//!
//! `IGESGraph_Color` carries a `(red, green, blue, color_name)`
//! definition that other entities reference by directory-entry index
//! via the "Color" field of their directory record. The IGES 5.3
//! spec encodes the colour as percentages (`0.0..=100.0`) plus an
//! optional case-sensitive colour-name HString.
//!
//! ## v1 status
//!
//! **Honest implementation** for the record-string construction.
//! This phase is the helper that callers feed into a future
//! IGES-with-colour writer; the full integration (inserting Type 314
//! records into the directory + parameter sections of an IGES file
//! produced by [`valenx_step_iges::iges::write`]) waits on Phase
//! 113.5, since the writer doesn't yet expose post-hoc entity
//! injection. The pure record-format function is testable and
//! correct on its own.

use crate::error::OcctExchangeError;

/// An RGB colour bound to a name, ready for IGES Type 314 encoding.
///
/// Components are stored in the 0..=1 range — the
/// [`iges_color_attributes()`] emitter rescales them into the
/// 0..=100 percent range the IGES spec wants on the wire.
#[derive(Clone, Debug, PartialEq)]
pub struct Color314 {
    /// Red channel (0..=1).
    pub r: f32,
    /// Green channel (0..=1).
    pub g: f32,
    /// Blue channel (0..=1).
    pub b: f32,
    /// Optional colour name — passed through to the HString record
    /// field. Empty string omits the name (parameter count drops by
    /// one).
    pub name: String,
}

/// Render a single Type 314 parameter-data record.
///
/// Format (IGES 5.3, Section 4.43):
/// `314, R%, G%, B%, [Hname];`
///
/// where `R%` / `G%` / `B%` are percentages and `Hname` follows the
/// IGES HString convention `<len>H<text>`. The record is terminated
/// by `;` and is meant to be inserted into the Parameter Data
/// section verbatim (the host writer slots the directory-entry
/// pointer in the trailing columns 66-72 itself).
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] when any component is non-finite
///   or outside the `0..=1` range.
pub fn iges_color_attributes(color: &Color314) -> Result<String, OcctExchangeError> {
    for (component, label) in [(color.r, "r"), (color.g, "g"), (color.b, "b")] {
        if !component.is_finite() {
            return Err(OcctExchangeError::bad_input(
                "color",
                format!("{label} component is non-finite"),
            ));
        }
        if !(0.0..=1.0).contains(&component) {
            return Err(OcctExchangeError::bad_input(
                "color",
                format!("{label} component {component} outside 0..=1"),
            ));
        }
    }
    let r_pct = color.r * 100.0;
    let g_pct = color.g * 100.0;
    let b_pct = color.b * 100.0;
    let body = if color.name.is_empty() {
        format!("314,{r_pct},{g_pct},{b_pct};")
    } else {
        format!(
            "314,{r_pct},{g_pct},{b_pct},{}H{};",
            color.name.len(),
            color.name,
        )
    };
    Ok(body)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_without_name_omits_hstring() {
        let c = Color314 {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            name: String::new(),
        };
        let s = iges_color_attributes(&c).unwrap();
        assert!(s.starts_with("314,"));
        assert!(s.contains("100"));
        assert!(s.ends_with(';'));
        assert!(!s.contains('H'));
    }

    #[test]
    fn record_with_name_emits_hstring() {
        let c = Color314 {
            r: 0.25,
            g: 0.5,
            b: 0.75,
            name: "Red".into(),
        };
        let s = iges_color_attributes(&c).unwrap();
        assert!(s.contains("3HRed"));
    }

    #[test]
    fn rejects_out_of_range() {
        let c = Color314 {
            r: 1.5,
            g: 0.0,
            b: 0.0,
            name: String::new(),
        };
        let err = iges_color_attributes(&c).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn rejects_nan() {
        let c = Color314 {
            r: f32::NAN,
            g: 0.0,
            b: 0.0,
            name: String::new(),
        };
        let err = iges_color_attributes(&c).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}

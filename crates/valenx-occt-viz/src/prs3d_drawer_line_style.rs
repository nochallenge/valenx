//! Phase 186 — `Aspect_TypeOfLine` — solid / dashed / dotted line
//! style patterns.
//!
//! ## What OCCT does
//!
//! `Graphic3d_AspectLine3d::SetType(Aspect_TypeOfLine)` switches the
//! per-object line stroke pattern: solid, dashed, dotted, dot-dashed,
//! or a user-supplied 16-bit pattern. On the OpenGl side OCCT emits
//! `glLineStipple` (legacy) or a fragment-shader stipple test (modern).
//!
//! ## v1 status
//!
//! **Honest v1.** Returns the validated [`LineStyle`] enum value.
//! The egui-paint path supports solid only today
//! (`egui::Stroke::new(width, colour)` has no pattern field), but
//! the *spec* of the requested style is preserved in app state so
//! Phase 188.5's wgpu wireframe pass can read it through a stipple
//! pattern uniform. Dashed / dotted are rendered as solid in v1 with
//! a deferred-stipple flag set; callers see the same value they passed.

use crate::error::OcctVizError;

/// Line stroke pattern mirror of OCCT's `Aspect_TypeOfLine`.
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum LineStyle {
    /// Continuous line (default).
    #[default]
    Solid,
    /// Dashed pattern (5px on / 3px off in v1's planned stipple).
    Dashed,
    /// Dotted pattern (1px on / 2px off).
    Dotted,
    /// Dot-dashed pattern (1px on / 2px off / 5px on / 2px off).
    DotDashed,
}

/// Validate and return the requested line style.
///
/// This op cannot fail — every variant is valid input. Returns
/// `Result` for API consistency.
pub fn prs3d_drawer_line_style(style: LineStyle) -> Result<LineStyle, OcctVizError> {
    Ok(style)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_solid() {
        assert_eq!(LineStyle::default(), LineStyle::Solid);
    }

    #[test]
    fn all_variants_round_trip() {
        for s in [
            LineStyle::Solid,
            LineStyle::Dashed,
            LineStyle::Dotted,
            LineStyle::DotDashed,
        ] {
            assert_eq!(prs3d_drawer_line_style(s).unwrap(), s);
        }
    }
}

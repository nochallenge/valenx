//! Phase 196 — `Prs3d_Drawer::SetSectionAspect` — draw section-plane
//! outlines + section curves.
//!
//! ## What OCCT does
//!
//! When a clipping plane is active ([`crate::v3d_viewer_clipping_plane()`]),
//! OCCT optionally renders the *section curves* — the 1D curves where
//! the plane intersects each face — as a coloured outline on top of
//! the clipped geometry. Plus an optional translucent fill of the
//! plane itself (the "cap") for visual reference. Configured via
//! `Prs3d_Drawer::SetSectionAspect(Prs3d_LineAspect)` for the curve
//! and a separate fill aspect for the cap.
//!
//! ## v1 status
//!
//! **Honest v1.** Returns a validated [`SectionDisplay`] config struct
//! that the caller (valenx_app's Mesh Toolbox section-plane panel)
//! stores in app state. The section curves themselves are already
//! computed by Valenx's `cut_overlay` path in
//! `valenx_app::viewport`; this op standardises the *display
//! aspect* (curve colour, cap fill colour, cap visible) so the
//! section view looks like OCCT's. Wires through to the clipping
//! plane infrastructure that Phase 169.5 will land — without that,
//! `cap_visible=true` is recorded but cap rendering remains pending
//! until clipping planes go live.

use crate::error::OcctVizError;
use crate::prs3d_drawer_face_color::{prs3d_drawer_face_color, FaceColorRgba};

/// Display configuration for a section plane.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SectionDisplay {
    /// Colour of the section curves (the plane-face intersection
    /// polylines).
    pub curve_color: FaceColorRgba,
    /// Colour of the optional translucent cap (the plane disk that
    /// fills the section). Alpha < 1 recommended.
    pub cap_color: FaceColorRgba,
    /// Whether to draw the cap. False = curves only (matches
    /// AutoCAD's section-line convention).
    pub cap_visible: bool,
}

/// Build a [`SectionDisplay`] from its three components.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] propagated from
///   [`crate::prs3d_drawer_face_color()`] if either colour is invalid.
pub fn prs3d_drawer_section_plane_display(
    curve_color: FaceColorRgba,
    cap_color: FaceColorRgba,
    cap_visible: bool,
) -> Result<SectionDisplay, OcctVizError> {
    // Re-validate via the constructor in case the caller hand-rolled a
    // `FaceColorRgba` literal with out-of-range channels.
    let curve = prs3d_drawer_face_color(
        curve_color.r,
        curve_color.g,
        curve_color.b,
        curve_color.a,
    )?;
    let cap = prs3d_drawer_face_color(cap_color.r, cap_color.g, cap_color.b, cap_color.a)?;
    Ok(SectionDisplay {
        curve_color: curve,
        cap_color: cap,
        cap_visible,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn yellow() -> FaceColorRgba {
        prs3d_drawer_face_color(1.0, 1.0, 0.0, 1.0).unwrap()
    }
    fn translucent_grey() -> FaceColorRgba {
        prs3d_drawer_face_color(0.5, 0.5, 0.5, 0.3).unwrap()
    }

    #[test]
    fn builds_full_display() {
        let s = prs3d_drawer_section_plane_display(yellow(), translucent_grey(), true).unwrap();
        assert!(s.cap_visible);
        assert_eq!(s.curve_color.r, 1.0);
        assert!(s.cap_color.a < 1.0);
    }

    #[test]
    fn rejects_invalid_curve_color() {
        let bad = FaceColorRgba {
            r: 2.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let err = prs3d_drawer_section_plane_display(bad, yellow(), false).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }
}

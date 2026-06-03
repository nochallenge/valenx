//! Phase 160 — `BRepAlgoAPI_Splitter` driven by a plane — split a
//! solid by a plane into two solids.
//!
//! ## What OCCT does
//!
//! `BRepAlgoAPI_Splitter(arguments, tools)` is the general
//! shape-splitter: it takes a list of "arguments" (the shapes to
//! split) and a list of "tools" (the shapes to split them with) and
//! produces a compound of all the resulting fragments. The plane-as-
//! tool variant uses a `BRepBuilderAPI_MakeFace(infinite_plane)` as
//! the cutter, producing two half-solids per input solid.
//!
//! Used for: making an exploded view (cut the assembly in half),
//! generating cross-sections, building "construction history" decals
//! from a model.
//!
//! Block 1's `valenx_occt_surface::algo_section()` does the
//! "section curve" variant (returns a wire, not a solid pair); this
//! phase returns the two halves as full solids.
//!
//! ## v1 status
//!
//! Stub — implementing the plane-split requires truck-shapeops to
//! accept a plane (or infinite face) as one of the boolean operands
//! and to split rather than union/intersect/subtract. truck's API
//! doesn't expose that primitive. Phase 160.5 ships when
//! valenx-cad lands a `split_solid_with_plane` op (likely via
//! a "construct large box, subtract, reconstruct" trick until
//! truck gains a native splitter).

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Two halves produced by [`local_op_split_solid_with_plane`].
#[derive(Debug)]
pub struct SplitHalves {
    /// Half on the positive side of the plane normal.
    pub positive_half: Solid,
    /// Half on the negative side of the plane normal.
    pub negative_half: Solid,
}

/// Split `solid` by the plane through `plane_origin` with normal
/// `plane_normal`.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a zero `plane_normal`.
/// - [`OcctAdvancedError::NotYetImplemented`] otherwise in v1.
pub fn local_op_split_solid_with_plane(
    _solid: &Solid,
    _plane_origin: [f64; 3],
    plane_normal: [f64; 3],
) -> Result<SplitHalves, OcctAdvancedError> {
    let n =
        (plane_normal[0].powi(2) + plane_normal[1].powi(2) + plane_normal[2].powi(2)).sqrt();
    if n < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "plane_normal",
            "must be non-zero",
        ));
    }
    Err(OcctAdvancedError::not_yet("local_op_split_solid_with_plane"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn rejects_zero_normal() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_split_solid_with_plane(&cube, [0.5; 3], [0.0; 3]).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn stub_with_valid_inputs() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_split_solid_with_plane(&cube, [0.5; 3], [0.0, 0.0, 1.0]).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}

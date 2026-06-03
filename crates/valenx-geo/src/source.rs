//! Where a piece of geometry came from.

use serde::{Deserialize, Serialize};

/// Geometry file formats we import from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SourceFormat {
    /// ISO 10303 STEP (default protocol: AP242 when possible,
    /// AP214 fallback).
    StepAp242,
    StepAp214,
    Iges,
    /// Binary or ASCII STL — the lowest common denominator.
    Stl,
    /// OCCT's native BRep serialization.
    BRep,
    /// FreeCAD's .FCStd archive (internally uses BRep).
    FcStd,
    /// Wavefront OBJ, for quick imports.
    Obj,
    /// Native `.valenx` parametric — the project's own format once
    /// the native CAD kernel matures.
    ValenxNative,
}

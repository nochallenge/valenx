//! Parse the `[geometry]` section of a FreeCAD case into a typed
//! spec the Python-script generator consumes.
//!
//! Scope (Phase 2 MVP):
//!
//! - Import STEP (AP242/AP214), IGES, STL, BREP, or native FCStd.
//! - Export baked STL + BREP + a JSON summary (parts, bbox, volume).
//! - Optional feature-tree summary for the Valenx browser.
//!
//! Parameter editing (edit linear-pattern count, see model rebuild)
//! lands once the scripting layer gets a real parameter surface.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Everything the FreeCAD Python script needs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeometryImportInput {
    pub source: PathBuf,
    pub source_format: SourceFormat,
    pub exports: Vec<ExportFormat>,
    /// Refine tessellation when exporting STL (`DeviationMM`). Lower
    /// means finer mesh.
    pub stl_deviation_mm: f64,
    /// Include a feature-tree summary per object (names, types,
    /// visibility).
    pub emit_feature_tree: bool,
    /// When `Some`, the FreeCAD script generates the primitive
    /// procedurally instead of opening `source`. Lets users build
    /// quick test geometry from `case.toml` without hand-modelling
    /// in FreeCAD's GUI first. The `source` field is ignored when
    /// this is set.
    #[serde(default)]
    pub primitive: Option<Primitive>,
}

/// Procedural primitive specifications. Each variant corresponds to a
/// `Part.makeXxx()` call in FreeCAD's Python API; the generated
/// script picks the right call when `GeometryImportInput.primitive`
/// is `Some`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Primitive {
    /// Axis-aligned box. Dimensions in millimetres.
    Cube { x: f64, y: f64, z: f64 },
    /// Cylinder along the +Z axis. Radius and height in millimetres.
    Cylinder { radius: f64, height: f64 },
    /// Sphere centred on the origin. Radius in millimetres.
    Sphere { radius: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SourceFormat {
    Step,
    Iges,
    Stl,
    Brep,
    Fcstd,
}

impl SourceFormat {
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_ascii_lowercase().as_str() {
            "stp" | "step" => Some(Self::Step),
            "igs" | "iges" => Some(Self::Iges),
            "stl" => Some(Self::Stl),
            "brep" | "brp" => Some(Self::Brep),
            "fcstd" => Some(Self::Fcstd),
            _ => None,
        }
    }

    /// FreeCAD importer module name — `Part.open` handles every
    /// geometry format we support; FCStd goes through the
    /// `FreeCAD.open` door.
    pub fn freecad_opener(self) -> &'static str {
        match self {
            Self::Fcstd => "FreeCAD.open",
            _ => "Part.open",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ExportFormat {
    Stl,
    Brep,
    Step,
    /// IGES — Initial Graphics Exchange Spec. Older
    /// interchange format than STEP, but still required by some
    /// legacy CAD/CAM toolchains (Fanuc CAM, older versions of
    /// SolidWorks). Lossier than STEP for assemblies + colour
    /// metadata; use STEP first when the downstream tool
    /// supports it.
    Iges,
}

impl ExportFormat {
    pub fn from_str_lenient(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "stl" => Some(Self::Stl),
            "brep" | "brp" => Some(Self::Brep),
            "step" | "stp" => Some(Self::Step),
            "iges" | "igs" => Some(Self::Iges),
            _ => None,
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Self::Stl => "stl",
            Self::Brep => "brep",
            Self::Step => "step",
            Self::Iges => "iges",
        }
    }
}

impl GeometryImportInput {
    pub fn from_case_dir(case_dir: &Path) -> Result<(CaseHeader, Self), AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )?;
        let case_def: CaseDef = toml::from_str(&text).map_err(|e| AdapterError::InvalidCase {
            case_path: case_toml.clone(),
            reason: format!("parse: {e}"),
        })?;
        let input = Self::from_case_def(&case_def).map_err(|e| with_case_path(e, &case_toml))?;
        Ok((case_def.case, input))
    }

    pub fn from_case_def(case_def: &CaseDef) -> Result<Self, AdapterError> {
        if case_def.case.physics != "geometry" {
            return Err(invalid(format!(
                "freecad adapter only handles physics=\"geometry\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let geom = case_def
            .section("geometry")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [geometry] section"))?;

        // Optional [geometry.primitive] block. When set, the script
        // generates the primitive procedurally and `source` becomes
        // optional / ignored — useful for quick-and-dirty test
        // geometry without having to model anything in the FreeCAD
        // GUI first.
        let primitive = parse_primitive(geom.get("primitive"))?;

        let source = geom
            .get("source")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        let source_format_explicit = geom
            .get("source_format")
            .and_then(|v| v.as_str())
            .and_then(SourceFormat::from_extension);

        let (source, source_format) = match (primitive, source) {
            (Some(_), maybe_src) => {
                // Primitive case: source/source_format are placeholders
                // since the script doesn't read them. Default to a
                // dummy STEP value so existing serialisers keep
                // working.
                (
                    maybe_src.unwrap_or_else(|| PathBuf::from("primitive")),
                    source_format_explicit.unwrap_or(SourceFormat::Step),
                )
            }
            (None, Some(src)) => {
                let fmt = source_format_explicit
                    .or_else(|| {
                        src.extension()
                            .and_then(|e| e.to_str())
                            .and_then(SourceFormat::from_extension)
                    })
                    .ok_or_else(|| {
                        invalid(format!(
                            "can't infer source format for {} — set [geometry] \
                             source_format = \"step|iges|stl|brep|fcstd\"",
                            src.display()
                        ))
                    })?;
                (src, fmt)
            }
            (None, None) => {
                return Err(invalid(
                    "[geometry] needs either `source = \"file.step\"` or a \
                     [geometry.primitive] block",
                ));
            }
        };

        let exports: Vec<ExportFormat> = geom
            .get("exports")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_str())
                    .filter_map(ExportFormat::from_str_lenient)
                    .collect()
            })
            .unwrap_or_else(|| vec![ExportFormat::Stl, ExportFormat::Brep]);

        let stl_deviation_mm = geom
            .get("stl_deviation_mm")
            .and_then(|v| v.as_float())
            .or_else(|| {
                geom.get("stl_deviation_mm")
                    .and_then(|v| v.as_integer())
                    .map(|i| i as f64)
            })
            .unwrap_or(0.1);

        let emit_feature_tree = geom
            .get("emit_feature_tree")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        Ok(Self {
            source,
            source_format,
            exports,
            stl_deviation_mm,
            emit_feature_tree,
            primitive,
        })
    }
}

fn as_f64(v: &toml::Value) -> Option<f64> {
    v.as_float().or_else(|| v.as_integer().map(|i| i as f64))
}

/// Parse the optional `[geometry.primitive]` sub-table.
fn parse_primitive(v: Option<&toml::Value>) -> Result<Option<Primitive>, AdapterError> {
    let Some(tbl) = v.and_then(|v| v.as_table()) else {
        return Ok(None);
    };
    let kind = tbl
        .get("kind")
        .and_then(|v| v.as_str())
        .ok_or_else(|| invalid("[geometry.primitive] needs `kind = \"cube|cylinder|sphere\"`"))?;
    let prim = match kind {
        "cube" | "box" => {
            let x = tbl.get("x").and_then(as_f64).unwrap_or(10.0);
            let y = tbl.get("y").and_then(as_f64).unwrap_or(10.0);
            let z = tbl.get("z").and_then(as_f64).unwrap_or(10.0);
            Primitive::Cube { x, y, z }
        }
        "cylinder" => {
            let radius = tbl.get("radius").and_then(as_f64).unwrap_or(5.0);
            let height = tbl.get("height").and_then(as_f64).unwrap_or(10.0);
            Primitive::Cylinder { radius, height }
        }
        "sphere" => {
            let radius = tbl.get("radius").and_then(as_f64).unwrap_or(5.0);
            Primitive::Sphere { radius }
        }
        other => {
            return Err(invalid(format!(
                "[geometry.primitive] unknown kind \"{other}\" — \
                 supported: cube, cylinder, sphere"
            )));
        }
    };
    Ok(Some(prim))
}

fn invalid(reason: impl Into<String>) -> AdapterError {
    AdapterError::InvalidCase {
        case_path: PathBuf::new(),
        reason: reason.into(),
    }
}

fn with_case_path(err: AdapterError, path: &Path) -> AdapterError {
    if let AdapterError::InvalidCase { case_path, reason } = err {
        if case_path.as_os_str().is_empty() {
            AdapterError::InvalidCase {
                case_path: path.to_path_buf(),
                reason,
            }
        } else {
            AdapterError::InvalidCase { case_path, reason }
        }
    } else {
        err
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_step_case() {
        let text = r#"
[case]
format  = "1.0"
name    = "bracket"
physics = "geometry"
solver  = "freecad.import"
mesh    = "(none)"

[geometry]
source  = "bracket.step"
exports = ["stl", "brep"]
stl_deviation_mm = 0.05
emit_feature_tree = true
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = GeometryImportInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.source_format, SourceFormat::Step);
        assert_eq!(input.exports.len(), 2);
        assert!((input.stl_deviation_mm - 0.05).abs() < 1e-9);
        assert!(input.emit_feature_tree);
    }

    #[test]
    fn infers_format_from_extension() {
        let text = r#"
[case]
format  = "1.0"
name    = "x"
physics = "geometry"
solver  = "freecad.import"
mesh    = "(none)"

[geometry]
source = "part.igs"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = GeometryImportInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.source_format, SourceFormat::Iges);
        // Also exercise the new Iges export variant via from_str_lenient
        assert_eq!(
            ExportFormat::from_str_lenient("iges"),
            Some(ExportFormat::Iges)
        );
        assert_eq!(
            ExportFormat::from_str_lenient("igs"),
            Some(ExportFormat::Iges)
        );
    }

    #[test]
    fn rejects_non_geometry_physics() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "cfd"
solver = "freecad.import"
mesh = "(none)"

[geometry]
source = "x.step"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        assert!(matches!(
            GeometryImportInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn unknown_extension_errors() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "geometry"
solver = "freecad.import"
mesh = "(none)"

[geometry]
source = "bracket.xyz"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = GeometryImportInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("source format"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn parses_primitive_cube_case() {
        let text = r#"
[case]
format = "1.0"
name = "test-cube"
physics = "geometry"
solver = "freecad.primitive"
mesh = "(none)"

[geometry]
exports = ["stl"]

[geometry.primitive]
kind = "cube"
x = 10.0
y = 20.0
z = 30.0
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = GeometryImportInput::from_case_def(&cd).expect("parse");
        assert_eq!(
            input.primitive,
            Some(Primitive::Cube {
                x: 10.0,
                y: 20.0,
                z: 30.0
            })
        );
    }

    #[test]
    fn parses_primitive_cylinder_with_defaults() {
        let text = r#"
[case]
format = "1.0"
name = "cyl"
physics = "geometry"
solver = "freecad.primitive"
mesh = "(none)"

[geometry]
[geometry.primitive]
kind = "cylinder"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = GeometryImportInput::from_case_def(&cd).expect("parse");
        // Defaults: radius=5, height=10.
        assert_eq!(
            input.primitive,
            Some(Primitive::Cylinder {
                radius: 5.0,
                height: 10.0
            })
        );
    }

    #[test]
    fn primitive_unknown_kind_errors() {
        let text = r#"
[case]
format = "1.0"
name = "bad"
physics = "geometry"
solver = "freecad.primitive"
mesh = "(none)"

[geometry]
[geometry.primitive]
kind = "torus"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = GeometryImportInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("torus"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn neither_source_nor_primitive_errors() {
        let text = r#"
[case]
format = "1.0"
name = "empty"
physics = "geometry"
solver = "freecad.import"
mesh = "(none)"

[geometry]
exports = ["stl"]
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = GeometryImportInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("source") || reason.contains("primitive"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }
}

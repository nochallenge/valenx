//! Parse the `[mesh]` section of a meshing case's `case.toml` into
//! a typed spec the `.geo` writer can consume.
//!
//! Phase 2 starts narrow on purpose — we support procedural primitive
//! domains (box, sphere) and external merges from STL/BREP files.
//! Feature-tree-driven CAD meshes land when the FreeCAD / OCCT
//! adapters graduate beyond scaffolding.

use std::path::{Path, PathBuf};

use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Everything the `.geo` writer needs to emit a meshing script.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshSpec {
    pub domain: Domain,
    pub sizes: MeshSizes,
    pub algorithm_2d: Algorithm2D,
    pub algorithm_3d: Algorithm3D,
    pub dim: MeshDim,
    pub physical_volume_name: String,
    pub physical_surface_name: String,
    /// Optional boundary-layer (prism-layer) spec. `None` disables
    /// inflation and produces a purely tetrahedral mesh. `Some(_)`
    /// stacks prism cells on the named surface for y+-resolving
    /// RANS work.
    pub boundary_layer: Option<BoundaryLayer>,
}

/// Gmsh `BoundaryLayer` field configuration. Maps onto gmsh's native
/// `Field[N] = BoundaryLayer` with Ratio / Size / NbLayers. Users
/// set the target y+ externally — converting y+ into a first-layer
/// height needs the flow's friction velocity, which isn't known at
/// meshing time — so the spec is over physical distance rather than
/// y+ directly.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryLayer {
    /// First-cell normal thickness (metres). 1 µm–1 mm is typical
    /// for external aero; set this to match your target y+ via a
    /// separate wall-units calculation.
    pub first_layer_thickness: f64,
    /// Growth ratio between successive layers (>1.0). 1.2 is a
    /// conservative default.
    pub growth_rate: f64,
    /// Total number of prism layers to stack.
    pub num_layers: u32,
    /// Names of physical surfaces inflation is applied to. Empty
    /// means the whole `Physical Surface` of the domain.
    pub surfaces: Vec<String>,
}

/// Source geometry for the mesh. `Box` and `Sphere` are procedural;
/// `MergeFile` pulls in an external surface mesh (STL, BRep, STEP)
/// that gmsh reclassifies and volume-meshes.
#[derive(Clone, Debug, PartialEq)]
pub enum Domain {
    Box {
        origin: [f64; 3],
        size: [f64; 3],
    },
    Sphere {
        center: [f64; 3],
        radius: f64,
    },
    MergeFile {
        /// Path relative to the case workdir.
        path: PathBuf,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshSizes {
    pub char_length_min: f64,
    pub char_length_max: f64,
}

impl Default for MeshSizes {
    fn default() -> Self {
        Self {
            char_length_min: 0.05,
            char_length_max: 0.2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Algorithm2D {
    /// MeshAdapt.
    MeshAdapt = 1,
    /// Automatic.
    Automatic = 2,
    /// Delaunay.
    Delaunay = 5,
    /// Frontal-Delaunay (default).
    FrontalDelaunay = 6,
    /// Frontal (quad-dominant).
    FrontalQuad = 8,
}

impl Algorithm2D {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "mesh-adapt" | "meshadapt" => Self::MeshAdapt,
            "automatic" | "auto" => Self::Automatic,
            "delaunay" => Self::Delaunay,
            "frontal-delaunay" | "frontal" => Self::FrontalDelaunay,
            "frontal-quad" | "quad" => Self::FrontalQuad,
            _ => Self::FrontalDelaunay,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Algorithm3D {
    /// Delaunay (default).
    Delaunay = 1,
    /// Frontal.
    Frontal = 4,
    /// HXT parallel.
    Hxt = 10,
}

impl Algorithm3D {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "delaunay" => Self::Delaunay,
            "frontal" => Self::Frontal,
            "hxt" | "parallel" => Self::Hxt,
            _ => Self::Delaunay,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MeshDim {
    Two,
    #[default]
    Three,
}

impl MeshDim {
    pub fn from_int(n: i64) -> Self {
        match n {
            2 => Self::Two,
            _ => Self::Three,
        }
    }

    pub fn as_int(self) -> i32 {
        match self {
            Self::Two => 2,
            Self::Three => 3,
        }
    }
}

impl MeshSpec {
    /// Load + parse a case directory's `case.toml` into a meshing spec.
    pub fn from_case_dir(case_dir: &Path) -> Result<(CaseHeader, Self), AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(&case_toml, valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize)?;
        let case_def: CaseDef = toml::from_str(&text).map_err(|e| AdapterError::InvalidCase {
            case_path: case_toml.clone(),
            reason: format!("parse: {e}"),
        })?;
        let spec = Self::from_case_def(&case_def).map_err(|e| with_case_path(e, &case_toml))?;
        Ok((case_def.case, spec))
    }

    pub fn from_case_def(case_def: &CaseDef) -> Result<Self, AdapterError> {
        if case_def.case.physics != "meshing" {
            return Err(invalid(format!(
                "gmsh adapter only handles physics=\"meshing\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let mesh = case_def
            .section("mesh")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [mesh] section"))?;

        let domain = parse_domain(mesh)?;

        let sizes = MeshSizes {
            char_length_min: mesh
                .get("characteristic_length_min")
                .and_then(|v| v.as_float())
                .or_else(|| mesh.get("char_length_min").and_then(|v| v.as_float()))
                .unwrap_or(0.05),
            char_length_max: mesh
                .get("characteristic_length_max")
                .and_then(|v| v.as_float())
                .or_else(|| mesh.get("characteristic_length").and_then(|v| v.as_float()))
                .or_else(|| mesh.get("char_length_max").and_then(|v| v.as_float()))
                .unwrap_or(0.2),
        };

        let algorithm_2d = mesh
            .get("algorithm_2d")
            .and_then(|v| v.as_str())
            .map(Algorithm2D::from_str_lenient)
            .unwrap_or(Algorithm2D::FrontalDelaunay);

        let algorithm_3d = mesh
            .get("algorithm_3d")
            .and_then(|v| v.as_str())
            .map(Algorithm3D::from_str_lenient)
            .unwrap_or(Algorithm3D::Delaunay);

        let dim = mesh
            .get("dim")
            .and_then(|v| v.as_integer())
            .map(MeshDim::from_int)
            .unwrap_or_default();

        let physical_volume_name = mesh
            .get("volume_name")
            .and_then(|v| v.as_str())
            .unwrap_or("domain")
            .to_string();
        let physical_surface_name = mesh
            .get("surface_name")
            .and_then(|v| v.as_str())
            .unwrap_or("boundary")
            .to_string();

        let boundary_layer = mesh
            .get("boundary_layer")
            .and_then(|v| v.as_table())
            .map(|tbl| BoundaryLayer {
                first_layer_thickness: tbl
                    .get("first_layer_thickness")
                    .and_then(|v| v.as_float())
                    .or_else(|| tbl.get("thickness").and_then(|v| v.as_float()))
                    .unwrap_or(1e-4),
                growth_rate: tbl
                    .get("growth_rate")
                    .and_then(|v| v.as_float())
                    .or_else(|| tbl.get("ratio").and_then(|v| v.as_float()))
                    .unwrap_or(1.2),
                num_layers: tbl
                    .get("num_layers")
                    .and_then(|v| v.as_integer())
                    .map(|i| i.max(1) as u32)
                    .unwrap_or(10),
                surfaces: tbl
                    .get("surfaces")
                    .and_then(|v| v.as_array())
                    .map(|arr| {
                        arr.iter()
                            .filter_map(|v| v.as_str())
                            .map(|s| s.to_string())
                            .collect()
                    })
                    .unwrap_or_default(),
            });

        Ok(Self {
            domain,
            sizes,
            algorithm_2d,
            algorithm_3d,
            dim,
            physical_volume_name,
            physical_surface_name,
            boundary_layer,
        })
    }
}

fn parse_domain(mesh: &toml::value::Table) -> Result<Domain, AdapterError> {
    let kind = mesh.get("type").and_then(|v| v.as_str()).ok_or_else(|| {
        invalid("[mesh] missing `type` (expected \"box\", \"sphere\", or \"merge\")")
    })?;
    match kind {
        "box" => {
            let origin = parse_vec3(mesh.get("origin"))
                .ok_or_else(|| invalid("[mesh] box domain requires `origin = [x, y, z]`"))?;
            let size = parse_vec3(mesh.get("size"))
                .ok_or_else(|| invalid("[mesh] box domain requires `size = [sx, sy, sz]`"))?;
            Ok(Domain::Box { origin, size })
        }
        "sphere" => {
            let center = parse_vec3(mesh.get("center"))
                .ok_or_else(|| invalid("[mesh] sphere domain requires `center = [x, y, z]`"))?;
            let radius = mesh
                .get("radius")
                .and_then(|v| v.as_float())
                .ok_or_else(|| invalid("[mesh] sphere domain requires numeric `radius`"))?;
            Ok(Domain::Sphere { center, radius })
        }
        "merge" | "merge-file" => {
            let path = mesh
                .get("source")
                .and_then(|v| v.as_str())
                .ok_or_else(|| invalid("[mesh] merge domain requires `source = \"file.stl\"`"))?;
            Ok(Domain::MergeFile {
                path: PathBuf::from(path),
            })
        }
        other => Err(invalid(format!(
            "[mesh] unsupported type \"{other}\" \
             (supported: box, sphere, merge)"
        ))),
    }
}

fn parse_vec3(v: Option<&toml::Value>) -> Option<[f64; 3]> {
    let arr = v?.as_array()?;
    if arr.len() != 3 {
        return None;
    }
    let mut out = [0.0f64; 3];
    for (i, el) in arr.iter().enumerate() {
        out[i] = el
            .as_float()
            .or_else(|| el.as_integer().map(|x| x as f64))?;
    }
    Some(out)
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

    fn sample_box_case() -> CaseDef {
        let text = r#"
[case]
format  = "1.0"
name    = "box-mesh"
physics = "meshing"
solver  = "gmsh.delaunay"
mesh    = "default"

[mesh]
type   = "box"
origin = [0.0, 0.0, 0.0]
size   = [1.0, 1.0, 1.0]
characteristic_length = 0.1
algorithm_3d = "delaunay"
"#;
        toml::from_str(text).unwrap()
    }

    #[test]
    fn parses_box_spec() {
        let cd = sample_box_case();
        let spec = MeshSpec::from_case_def(&cd).expect("parse");
        assert_eq!(
            spec.domain,
            Domain::Box {
                origin: [0.0, 0.0, 0.0],
                size: [1.0, 1.0, 1.0],
            }
        );
        assert!((spec.sizes.char_length_max - 0.1).abs() < 1e-9);
        assert_eq!(spec.algorithm_3d, Algorithm3D::Delaunay);
        assert_eq!(spec.dim, MeshDim::Three);
    }

    #[test]
    fn parses_sphere_spec() {
        let text = r#"
[case]
format = "1.0"
name = "sphere"
physics = "meshing"
solver = "gmsh.delaunay"
mesh = "default"

[mesh]
type = "sphere"
center = [0, 0, 0]
radius = 0.5
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let spec = MeshSpec::from_case_def(&cd).unwrap();
        match spec.domain {
            Domain::Sphere { center, radius } => {
                assert_eq!(center, [0.0, 0.0, 0.0]);
                assert!((radius - 0.5).abs() < 1e-9);
            }
            other => panic!("wrong domain: {other:?}"),
        }
    }

    #[test]
    fn parses_merge_spec() {
        let text = r#"
[case]
format = "1.0"
name = "merge"
physics = "meshing"
solver = "gmsh.delaunay"
mesh = "default"

[mesh]
type = "merge"
source = "inlet.stl"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let spec = MeshSpec::from_case_def(&cd).unwrap();
        match spec.domain {
            Domain::MergeFile { path } => assert_eq!(path, PathBuf::from("inlet.stl")),
            other => panic!("wrong domain: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_meshing_physics() {
        let mut cd = sample_box_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            MeshSpec::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn algorithm_parsing_is_lenient() {
        assert_eq!(
            Algorithm2D::from_str_lenient("Frontal-Delaunay"),
            Algorithm2D::FrontalDelaunay
        );
        assert_eq!(
            Algorithm2D::from_str_lenient("junk"),
            Algorithm2D::FrontalDelaunay
        );
        assert_eq!(Algorithm3D::from_str_lenient("HXT"), Algorithm3D::Hxt);
    }
}

//! Post-load validation for [`super::LoadedProject`].
//!
//! The `loader` module enforces the structural / safety rules that
//! must hold for the project to be usable at all (manifest exists,
//! parses, paths stay inside the root). This module adds *advisory*
//! checks that surface common operator mistakes without blocking
//! the load:
//!
//! - case `solver` strings should follow the `<adapter>.<flavour>`
//!   convention (e.g. `openfoam.simpleFoam`)
//! - case `physics` should be one of the known categories
//! - geometry entries should have a non-empty `id`
//! - mesh entries should reference a geometry id that exists
//!
//! Returns a `Vec<ProjectWarning>` — empty when everything checks
//! out. The UI surfaces the warnings in the project pane so users
//! see them on load.

use serde::{Deserialize, Serialize};

use super::LoadedProject;

/// Known physics categories. Cases with a `physics` value outside
/// this set still load (adapters may add their own categories) but
/// the validator emits a warning so typos are visible.
const KNOWN_PHYSICS: &[&str] = &[
    "cfd",
    "fea",
    "structural", // legacy alias for "fea"
    "thermal",
    "chemistry",
    "molecular-dynamics",
    "md", // shorthand
    "electromagnetics",
    "em", // shorthand
    "battery",
    "dynamics",
    "multibody",
    "multi-physics",
    "coupling", // alias for "multi-physics"
    "geometry",
    "mesh",
];

/// One advisory warning surfaced by [`validate`]. Each carries a
/// short `code` (machine-stable) + a human-readable `message`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProjectWarning {
    pub code: &'static str,
    pub message: String,
}

/// Run every advisory check against a loaded project. Returns the
/// collected warnings in declaration order — empty when the
/// project is clean.
pub fn validate(project: &LoadedProject) -> Vec<ProjectWarning> {
    let mut warnings: Vec<ProjectWarning> = Vec::new();

    // Geometry entries with empty ids.
    for (i, g) in project.project.geometry.entries.iter().enumerate() {
        if g.id.trim().is_empty() {
            warnings.push(ProjectWarning {
                code: "geometry.empty_id",
                message: format!(
                    "geometry entry #{i} has an empty `id` — every geometry needs a stable name to be referenced from cases"
                ),
            });
        }
    }

    // Cases that reference an unknown mesh id (mesh declared in
    // [mesh.<name>] but referenced from `[case].mesh = "..."` —
    // typos here surface as silent "missing mesh" runtime errors
    // unless flagged on load).
    let mesh_ids: std::collections::HashSet<String> =
        project.project.mesh.keys().cloned().collect();
    // The literal `(none)` is a valid sentinel meaning "no mesh
    // applicable" (chemistry / molecular dynamics / coupling cases
    // use it). Only warn for non-sentinel mesh refs not in the map.
    for (case_name, case_def) in project.cases.iter() {
        let m = &case_def.case.mesh;
        if !m.is_empty() && m != "(none)" && !mesh_ids.contains(m) {
            warnings.push(ProjectWarning {
                code: "case.unknown_mesh",
                message: format!("case `{case_name}` mesh `{m}` is not declared in [mesh.*]"),
            });
        }
    }

    // Cases with non-conventional solver / physics strings.
    for (case_name, case_def) in project.cases.iter() {
        let solver = &case_def.case.solver;
        if !solver.is_empty() && !solver.contains('.') {
            warnings.push(ProjectWarning {
                code: "case.solver_no_dot",
                message: format!(
                    "case `{case_name}` solver `{solver}` doesn't follow the `<adapter>.<flavour>` convention — \
                     adapter id resolution falls back to the bare string"
                ),
            });
        }
        let physics = &case_def.case.physics;
        if !physics.is_empty() && !KNOWN_PHYSICS.contains(&physics.as_str()) {
            warnings.push(ProjectWarning {
                code: "case.physics_unknown",
                message: format!(
                    "case `{case_name}` physics `{physics}` is not in the known set ({}); \
                     this may be a typo or a new physics category",
                    KNOWN_PHYSICS.join(", ")
                ),
            });
        }
        // Empty case name in `cases.order` would have been caught by
        // the loader; but a case file with an empty `name` field is
        // an advisory issue.
        if case_def.case.name.trim().is_empty() {
            warnings.push(ProjectWarning {
                code: "case.empty_name",
                message: format!(
                    "case `{case_name}` has an empty `[case].name` field — UI labels will fall back to the directory name"
                ),
            });
        }
    }

    warnings
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project::case_def::CaseHeader;
    use crate::project::manifest::*;
    use crate::CaseDef;
    use std::collections::BTreeMap;
    use std::path::PathBuf;

    fn empty_project() -> LoadedProject {
        LoadedProject {
            root: PathBuf::from("/tmp/test"),
            project: Project {
                project: ProjectHeader {
                    format: "1.0".into(),
                    name: "test".into(),
                    valenx_min: None,
                    created: None,
                    modified: None,
                    author: None,
                    description: None,
                },
                geometry: GeometrySection {
                    entries: Vec::new(),
                },
                mesh: BTreeMap::new(),
                cases: CasesSection { order: Vec::new() },
                ui: UiSection::default(),
                units: UnitsConfig::default(),
            },
            tools_lock: None,
            cases: BTreeMap::new(),
        }
    }

    fn case_with(solver: &str, physics: &str, name: &str) -> CaseDef {
        case_with_mesh(solver, physics, name, "(none)")
    }

    fn case_with_mesh(solver: &str, physics: &str, name: &str, mesh: &str) -> CaseDef {
        CaseDef {
            case: CaseHeader {
                format: "1.0".into(),
                name: name.into(),
                physics: physics.into(),
                solver: solver.into(),
                mesh: mesh.into(),
                description: None,
            },
            sections: BTreeMap::new(),
        }
    }

    #[test]
    fn empty_project_has_no_warnings() {
        let p = empty_project();
        assert!(validate(&p).is_empty());
    }

    #[test]
    fn warns_when_case_solver_lacks_dot() {
        let mut p = empty_project();
        p.cases
            .insert("smoke".into(), case_with("openfoam", "cfd", "smoke"));
        let w = validate(&p);
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, "case.solver_no_dot");
    }

    #[test]
    fn warns_when_case_physics_unknown() {
        let mut p = empty_project();
        p.cases.insert(
            "weird".into(),
            case_with("openfoam.simpleFoam", "magic", "weird"),
        );
        let w = validate(&p);
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, "case.physics_unknown");
        assert!(w[0].message.contains("magic"));
    }

    #[test]
    fn warns_when_case_name_is_empty() {
        let mut p = empty_project();
        p.cases.insert(
            "anon".into(),
            case_with("openfoam.simpleFoam", "cfd", "   "),
        );
        let w = validate(&p);
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, "case.empty_name");
    }

    #[test]
    fn warns_when_geometry_id_is_empty() {
        let mut p = empty_project();
        p.project.geometry.entries.push(GeometryEntry {
            id: "".into(),
            source: PathBuf::from("airfoil.step"),
            format: "step".into(),
        });
        let w = validate(&p);
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, "geometry.empty_id");
    }

    #[test]
    fn warns_when_case_refs_unknown_mesh() {
        let mut p = empty_project();
        // Mesh declared as "fluid"...
        p.project.mesh.insert(
            "fluid".into(),
            MeshEntry {
                source: PathBuf::from("fluid.msh"),
                config_hash: None,
            },
        );
        // ...case references the typo'd "flud".
        p.cases.insert(
            "smoke".into(),
            case_with_mesh("openfoam.simpleFoam", "cfd", "smoke", "flud"),
        );
        let w = validate(&p);
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, "case.unknown_mesh");
        assert!(w[0].message.contains("flud"));
    }

    #[test]
    fn case_mesh_sentinel_none_does_not_warn() {
        // (none) is the convention for cases without a mesh
        // (chemistry / MD / coupling). Must NOT warn.
        let mut p = empty_project();
        p.cases.insert(
            "ch".into(),
            case_with("cantera.equilibrium", "chemistry", "ch"),
        );
        assert!(validate(&p).is_empty());
    }

    #[test]
    fn known_physics_aliases_dont_warn() {
        // `structural` (legacy alias for `fea`), `md`, `em`,
        // `coupling` (alias for `multi-physics`) — all in
        // KNOWN_PHYSICS — should NOT trigger warnings.
        let mut p = empty_project();
        for (name, physics) in [
            ("a", "structural"),
            ("b", "md"),
            ("c", "em"),
            ("d", "coupling"),
        ] {
            p.cases
                .insert(name.into(), case_with("openfoam.simpleFoam", physics, name));
        }
        let w = validate(&p);
        assert!(w.is_empty(), "got: {w:?}");
    }

    #[test]
    fn multiple_issues_are_reported_in_order() {
        let mut p = empty_project();
        p.project.geometry.entries.push(GeometryEntry {
            id: "".into(),
            source: PathBuf::from("a.step"),
            format: "step".into(),
        });
        p.cases
            .insert("smoke".into(), case_with("openfoam", "magic", "smoke"));
        let w = validate(&p);
        assert!(w.len() >= 3);
        // Geometry warnings come first, then mesh / cases.
        assert_eq!(w[0].code, "geometry.empty_id");
        // Followed by the case warnings (solver_no_dot,
        // physics_unknown — order within case checks is stable).
        assert!(w.iter().any(|x| x.code == "case.solver_no_dot"));
        assert!(w.iter().any(|x| x.code == "case.physics_unknown"));
    }
}

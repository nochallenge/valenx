//! **Per-file registry of agent-bridge 3-D mesh producers.**
//!
//! The Workbench+Agent bridge (`show_3d{kind}`) used to build each model with a
//! per-kind `else if kind == "<x>"` arm inlined into
//! [`crate::agent_commands`]'s shared reducer — so wiring a new model meant
//! editing that one shared `match`, which serialises parallel work and is a
//! merge-conflict magnet. This module replaces those arms with a single
//! generic lookup ([`lookup`](crate::products_registry::lookup)) keyed by the
//! wire `kind` string, so the reducer
//! never grows a new branch again.
//!
//! ## What lives where
//!
//! The *substance* of each product — how its [`crate::WorkspaceProduct`] is
//! assembled from that tool's canonical mesh / camera / readout-row producers —
//! lives in **that tool's own module** as a pure `pub(crate) fn …() ->
//! WorkspaceProduct` builder (e.g. `crate::rocket_workbench::rocket_product`).
//! This file only holds the tiny `kind → builder` table in
//! [`lookup`](crate::products_registry::lookup). Adding a
//! tool is therefore a *one-line* edit to the table here plus a self-contained
//! builder in the new tool's file — the per-tool code (the part that actually
//! varies and that two contributors might touch at once) is fully isolated, so
//! parallel wiring conflicts shrink to one trivial line in this table.
//!
//! ## Why a shared-`match` table and not the `inventory` crate
//!
//! The plan offered `inventory` (link-time distributed-slice registration, no
//! central list at all) as the first choice. We deliberately use the fallback
//! shared-`match` table instead, because this workspace links Windows builds
//! with **`rust-lld` / `lld-link`** (see `.cargo/config.toml`, committed
//! 2026-06-20). `inventory` populates its slice via `ctor`-style
//! life-before-main static registration, which is exactly the construct a
//! non-default linker can dead-strip when the registering module isn't
//! otherwise referenced — and the gate here is a `cargo test -p valenx-app
//! --lib` build, the case most prone to that stripping. A `match` is
//! linker-agnostic and resolved at compile time, so the registry can never
//! silently lose a kind. (It also avoids adding a third-party dependency +
//! `deny.toml` / `Cargo.lock` churn for a five-entry table.)
//!
//! ## Builder contract
//!
//! Every builder is a `fn() -> WorkspaceProduct`
//! ([`MeshProducerEntry::build`](crate::products_registry::MeshProducerEntry::build))
//! — pure, app-state-free, built only from that tool's canonical inputs — so
//! the reducer can call it with nothing but the channel it already knows.
//! Behaviour is byte-for-byte what the old inline arms produced.
//!
//! ## Adding a new 3-D tool — the one-liner pattern
//!
//! Copy this into the new tool's own module (PHASE-C wiring subagents: this is
//! the whole change on the producer side — no edit to `agent_commands`):
//!
//! ```ignore
//! // 1. In `crate::foo_workbench` (the tool's own file): a pure builder that
//! //    assembles the WorkspaceProduct from the tool's canonical producers.
//! pub(crate) fn foo_product() -> crate::WorkspaceProduct {
//!     let (mesh, lines) = foo_loaded_mesh();              // your existing producer
//!     let camera = foo_camera(&mesh.mesh);               // your existing camera
//!     crate::WorkspaceProduct {
//!         title: "Foo".into(), lines, mesh: Some(mesh),
//!         vertex_colors: None, camera, kind2d: None, last_export: None,
//!     }
//! }
//! ```
//!
//! …then add the single table line in
//! [`lookup`](crate::products_registry::lookup) below:
//! `"foo" => Some(crate::foo_workbench::foo_product),`. That is the only shared
//! edit, and it is a one-liner — everything else lives in the tool's file.

use crate::types::LoadedMesh;
use crate::WorkspaceProduct;

/// Wrap a freshly-built triangle [`valenx_mesh::Mesh`] into a
/// fully-populated [`LoadedMesh`] (mesh + quality report + aspect-ratio /
/// skewness histograms on the default buckets), tagged `path`.
///
/// The shared helper the per-tool `*_product` builders use so the
/// quality-metric plumbing (identical in every workbench's own
/// `load_*_3d`) lives in exactly one place instead of being copy-pasted 18×.
/// `path` is the synthetic `<kind>/valenx-…` tag the tile shows.
pub(crate) fn loaded_mesh_from(mesh: valenx_mesh::Mesh, path: &str) -> LoadedMesh {
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    LoadedMesh {
        path: std::path::PathBuf::from(path),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    }
}

/// A fixed 3/4-view [`valenx_viz::OrbitCamera`] framing `mesh`'s AABB at the
/// same hero angle (azimuth 35°, elevation 22°) as
/// [`crate::rocket_workbench::lv1_camera`] / the gear-train camera, for a
/// Workbench+Agent product's per-tile 3-D view. The shared camera builder the
/// per-tool `*_product` builders use (so framing stays identical across tiles).
pub(crate) fn camera_for(mesh: &valenx_mesh::Mesh) -> valenx_viz::OrbitCamera {
    let mut camera = valenx_viz::OrbitCamera::default();
    if let Some((min, max)) = crate::mesh_loader::mesh_bounding_box(mesh) {
        camera.frame_bounds(min, max);
    }
    camera.azimuth_deg = 35.0;
    camera.elevation_deg = 22.0;
    camera
}

/// Split a workbench's monospace readout (the `compute(&state)` string each
/// tool already formats for its panel) into the per-row result lines a
/// [`WorkspaceProduct`] carries: one trimmed line per source line, with the
/// blank spacer rows dropped. Lets a `*_product` builder reuse the tool's own
/// computed headline numbers as the tile's result rows rather than
/// re-deriving them.
pub(crate) fn lines_from_readout(readout: &str) -> Vec<String> {
    readout
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(str::to_owned)
        .collect()
}

/// One registry entry: a wire `kind` string mapped to a pure builder that
/// returns the [`WorkspaceProduct`] for that 3-D model.
///
/// The builder is app-state-free (the existing producers build from canonical
/// inputs only), so the bridge can invoke it knowing nothing but the channel it
/// publishes into. Returned by [`lookup`] so callers can read [`Self::kind`]
/// (e.g. for diagnostics) as well as invoke [`Self::build`].
#[derive(Clone, Copy)]
pub struct MeshProducerEntry {
    /// The `show_3d` wire `kind` this entry answers to (e.g. `"rocket"`).
    pub kind: &'static str,
    /// Pure builder for this kind's product — same output as the old inline
    /// reducer arm.
    pub build: fn() -> WorkspaceProduct,
}

/// Resolve a `show_3d` `kind` to its registry entry, or `None` for an unknown
/// kind (the reducer then skips it safely — no panic, no placeholder churn,
/// matching the rest of its bad-input handling).
///
/// **This `match` is the single shared edit point for 3-D mesh tools**: each
/// arm is one line pairing a wire `kind` with the per-tool builder in that
/// tool's own module. The substantive per-tool code lives in those builders,
/// not here — so adding a kind is a one-line addition (see the module docs for
/// the copy-paste pattern). Note `dna` is intentionally absent: `show_3d:dna`
/// is a *text* card handled directly in the reducer (no mesh), and the 2-D
/// `show_2d` drawings (`rcbeam` / `dna`) have their own separate path.
pub fn lookup(kind: &str) -> Option<MeshProducerEntry> {
    let build: fn() -> WorkspaceProduct = match kind {
        "rocket" => crate::rocket_workbench::rocket_product,
        "gear" => crate::gears_workbench::gear_product,
        "bracket" => crate::bracket_product::bracket_workspace_product,
        "rcbeam" => crate::rcbeam_workbench::rcbeam_product,
        "fem" => crate::fem_workbench::fem_product,
        // Machine-design family (each builder lives in its own workbench
        // module; see that module's `*_product`).
        "geartooth" => crate::geartooth_workbench::geartooth_product,
        "gearbox" => crate::gearbox_workbench::gearbox_product,
        "bearing" => crate::bearing_workbench::bearing_product,
        "clutch" => crate::clutch_workbench::clutch_product,
        "brake" => crate::brake_workbench::brake_product,
        "pulley" => crate::pulley_workbench::pulley_product,
        "flywheel" => crate::flywheel_workbench::flywheel_product,
        "leadscrew" => crate::leadscrew_workbench::leadscrew_product,
        "screwthread" => crate::screwthread_workbench::screwthread_product,
        "shaftdesign" => crate::shaftdesign_workbench::shaftdesign_product,
        "camdynamics" => crate::camdynamics_workbench::camdynamics_product,
        "conveyor" => crate::conveyor_workbench::conveyor_product,
        "bolt" => crate::bolt_workbench::bolt_product,
        "rivet" => crate::rivet_workbench::rivet_product,
        // Structural family.
        "beam" => crate::beam_workbench::beam_product,
        "truss" => crate::truss_workbench::truss_product,
        "plate" => crate::plate_workbench::plate_product,
        "buckling" => crate::buckling_workbench::buckling_product,
        // Civil / strength-of-materials / mechanisms / vibration family (each
        // builder lives in its own workbench module; see that module's
        // `*_product`).
        "columnsteel" => crate::columnsteel_workbench::columnsteel_product,
        "retainingwall" => crate::retainingwall_workbench::retainingwall_product,
        "soilbearing" => crate::soilbearing_workbench::soilbearing_product,
        "statics" => crate::statics_workbench::statics_product,
        "mohr" => crate::mohr_workbench::mohr_product,
        "torsion" => crate::torsion_workbench::torsion_product,
        "fatigue" => crate::fatigue_workbench::fatigue_product,
        "fracture" => crate::fracture_workbench::fracture_product,
        "creep" => crate::creep_workbench::creep_product,
        "pressurevessel" => crate::pressurevessel_workbench::pressurevessel_product,
        "straingauge" => crate::straingauge_workbench::straingauge_product,
        "strainrosette" => crate::strainrosette_workbench::strainrosette_product,
        "rail" => crate::rail_workbench::rail_product,
        "springdesign" => crate::springdesign_workbench::springdesign_product,
        "springs" => crate::springs_workbench::springs_product,
        "springcombination" => crate::springcombination_workbench::springcombination_product,
        "leverage" => crate::leverage_workbench::leverage_product,
        "inclinedplane" => crate::inclinedplane_workbench::inclinedplane_product,
        "vibration" => crate::vibration_workbench::vibration_product,
        // Electric-machines / power-transmission / mechanisms family (each
        // builder lives in its own workbench module; see that module's
        // `*_product`).
        "dcmotor" => crate::dcmotor_workbench::dcmotor_product,
        "inductionmotor" => crate::inductionmotor_workbench::inductionmotor_product,
        "beltdrive" => crate::beltdrive_workbench::beltdrive_product,
        "chaindrive" => crate::chaindrive_workbench::chaindrive_product,
        "fourbar" => crate::fourbar_workbench::fourbar_product,
        // Thermal / HVAC / thermodynamics / measurement family.
        "thermalexpansion" => crate::thermalexpansion_workbench::thermalexpansion_product,
        "dimensional" => crate::dimensional_workbench::dimensional_product,
        "projectile" => crate::projectile_workbench::projectile_product,
        "heattransfer" => crate::heattransfer_workbench::heattransfer_product,
        "insulation" => crate::insulation_workbench::insulation_product,
        "heatexchanger" => crate::heatexchanger_workbench::heatexchanger_product,
        "heatpump" => crate::heatpump_workbench::heatpump_product,
        "refrigeration" => crate::refrigeration_workbench::refrigeration_product,
        "psychrometrics" => crate::psychrometrics_workbench::psychrometrics_product,
        "thermocouple" => crate::thermocouple_workbench::thermocouple_product,
        "thermistor" => crate::thermistor_workbench::thermistor_product,
        "thermocycle" => crate::thermocycle_workbench::thermocycle_product,
        "fanlaws" => crate::fanlaws_workbench::fanlaws_product,
        // Fluid-mechanics / hydraulics / thermo-fluids / electronics family
        // (each builder lives in its own workbench module; see that module's
        // `*_product`).
        "pump" => crate::pump_workbench::pump_product,
        "pipeflow" => crate::pipeflow_workbench::pipeflow_product,
        "pipenetwork" => crate::pipenetwork_workbench::pipenetwork_product,
        "hydraulics" => crate::hydraulics_workbench::hydraulics_product,
        "pneumatics" => crate::pneumatics_workbench::pneumatics_product,
        "fluidstatics" => crate::fluidstatics_workbench::fluidstatics_product,
        "openchannel" => crate::openchannel_workbench::openchannel_product,
        "weir" => crate::weir_workbench::weir_product,
        "orifice" => crate::orifice_workbench::orifice_product,
        "combustion" => crate::combustion_workbench::combustion_product,
        "diffusion" => crate::diffusion_workbench::diffusion_product,
        "marine" => crate::marine_workbench::marine_product,
        "resistornetwork" => crate::resistornetwork_workbench::resistornetwork_product,
        "capacitor" => crate::capacitor_workbench::capacitor_product,
        "opamp" => crate::opamp_workbench::opamp_product,
        "bjt" => crate::bjt_workbench::bjt_product,
        "mosfet" => crate::mosfet_workbench::mosfet_product,
        "rectifier" => crate::rectifier_workbench::rectifier_product,
        _ => return None,
    };
    Some(MeshProducerEntry {
        kind: kind_static(kind)?,
        build,
    })
}

/// Map a looked-up `kind` to its `'static` spelling for [`MeshProducerEntry::kind`].
/// Kept in lockstep with [`lookup`]'s arms so the returned `kind` is always the
/// canonical literal (not a borrow of the caller's input).
fn kind_static(kind: &str) -> Option<&'static str> {
    Some(match kind {
        "rocket" => "rocket",
        "gear" => "gear",
        "bracket" => "bracket",
        "rcbeam" => "rcbeam",
        "fem" => "fem",
        "geartooth" => "geartooth",
        "gearbox" => "gearbox",
        "bearing" => "bearing",
        "clutch" => "clutch",
        "brake" => "brake",
        "pulley" => "pulley",
        "flywheel" => "flywheel",
        "leadscrew" => "leadscrew",
        "screwthread" => "screwthread",
        "shaftdesign" => "shaftdesign",
        "camdynamics" => "camdynamics",
        "conveyor" => "conveyor",
        "bolt" => "bolt",
        "rivet" => "rivet",
        "beam" => "beam",
        "truss" => "truss",
        "plate" => "plate",
        "buckling" => "buckling",
        "columnsteel" => "columnsteel",
        "retainingwall" => "retainingwall",
        "soilbearing" => "soilbearing",
        "statics" => "statics",
        "mohr" => "mohr",
        "torsion" => "torsion",
        "fatigue" => "fatigue",
        "fracture" => "fracture",
        "creep" => "creep",
        "pressurevessel" => "pressurevessel",
        "straingauge" => "straingauge",
        "strainrosette" => "strainrosette",
        "rail" => "rail",
        "springdesign" => "springdesign",
        "springs" => "springs",
        "springcombination" => "springcombination",
        "leverage" => "leverage",
        "inclinedplane" => "inclinedplane",
        "vibration" => "vibration",
        "dcmotor" => "dcmotor",
        "inductionmotor" => "inductionmotor",
        "beltdrive" => "beltdrive",
        "chaindrive" => "chaindrive",
        "fourbar" => "fourbar",
        "thermalexpansion" => "thermalexpansion",
        "dimensional" => "dimensional",
        "projectile" => "projectile",
        "heattransfer" => "heattransfer",
        "insulation" => "insulation",
        "heatexchanger" => "heatexchanger",
        "heatpump" => "heatpump",
        "refrigeration" => "refrigeration",
        "psychrometrics" => "psychrometrics",
        "thermocouple" => "thermocouple",
        "thermistor" => "thermistor",
        "thermocycle" => "thermocycle",
        "fanlaws" => "fanlaws",
        "pump" => "pump",
        "pipeflow" => "pipeflow",
        "pipenetwork" => "pipenetwork",
        "hydraulics" => "hydraulics",
        "pneumatics" => "pneumatics",
        "fluidstatics" => "fluidstatics",
        "openchannel" => "openchannel",
        "weir" => "weir",
        "orifice" => "orifice",
        "combustion" => "combustion",
        "diffusion" => "diffusion",
        "marine" => "marine",
        "resistornetwork" => "resistornetwork",
        "capacitor" => "capacitor",
        "opamp" => "opamp",
        "bjt" => "bjt",
        "mosfet" => "mosfet",
        "rectifier" => "rectifier",
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact set of 3-D mesh kinds the registry is expected to resolve, so a
    /// future addition that forgets the table line trips the count assertion.
    const KNOWN_3D_KINDS: &[&str] = &[
        "rocket",
        "gear",
        "bracket",
        "rcbeam",
        "fem",
        // Machine-design + structural workbenches wired into the bridge.
        "geartooth",
        "gearbox",
        "bearing",
        "clutch",
        "brake",
        "pulley",
        "flywheel",
        "leadscrew",
        "screwthread",
        "shaftdesign",
        "camdynamics",
        "conveyor",
        "bolt",
        "rivet",
        "beam",
        "truss",
        "plate",
        "buckling",
        // Civil / strength-of-materials / mechanisms / vibration workbenches
        // wired into the bridge.
        "columnsteel",
        "retainingwall",
        "soilbearing",
        "statics",
        "mohr",
        "torsion",
        "fatigue",
        "fracture",
        "creep",
        "pressurevessel",
        "straingauge",
        "strainrosette",
        "rail",
        "springdesign",
        "springs",
        "springcombination",
        "leverage",
        "inclinedplane",
        "vibration",
        // Electric-machines / power-transmission / mechanisms / thermal / HVAC /
        // thermodynamics / measurement workbenches wired into the bridge.
        "dcmotor",
        "inductionmotor",
        "beltdrive",
        "chaindrive",
        "fourbar",
        "thermalexpansion",
        "dimensional",
        "projectile",
        "heattransfer",
        "insulation",
        "heatexchanger",
        "heatpump",
        "refrigeration",
        "psychrometrics",
        "thermocouple",
        "thermistor",
        "thermocycle",
        "fanlaws",
        // Fluid-mechanics / hydraulics / thermo-fluids / electronics
        // workbenches wired into the bridge.
        "pump",
        "pipeflow",
        "pipenetwork",
        "hydraulics",
        "pneumatics",
        "fluidstatics",
        "openchannel",
        "weir",
        "orifice",
        "combustion",
        "diffusion",
        "marine",
        "resistornetwork",
        "capacitor",
        "opamp",
        "bjt",
        "mosfet",
        "rectifier",
    ];

    /// The machine-design / structural / civil / strength-of-materials /
    /// mechanisms / vibration workbench kinds wired into the bridge — every one
    /// must resolve via [`lookup`] and build a non-empty `Tri3` mesh product
    /// carrying a title and at least one readout row (the agent-bridge
    /// `show_3d{kind}` payload).
    const WIRED_WORKBENCH_KINDS: &[&str] = &[
        "geartooth",
        "gearbox",
        "bearing",
        "clutch",
        "brake",
        "pulley",
        "flywheel",
        "leadscrew",
        "screwthread",
        "shaftdesign",
        "camdynamics",
        "conveyor",
        "bolt",
        "rivet",
        "beam",
        "truss",
        "plate",
        "buckling",
        // The 19 civil / strength-of-materials / mechanisms / vibration
        // workbenches wired in this change.
        "columnsteel",
        "retainingwall",
        "soilbearing",
        "statics",
        "mohr",
        "torsion",
        "fatigue",
        "fracture",
        "creep",
        "pressurevessel",
        "straingauge",
        "strainrosette",
        "rail",
        "springdesign",
        "springs",
        "springcombination",
        "leverage",
        "inclinedplane",
        "vibration",
        // The 18 electric-machines / power-transmission / mechanisms / thermal /
        // HVAC / thermodynamics / measurement workbenches wired in this change.
        "dcmotor",
        "inductionmotor",
        "beltdrive",
        "chaindrive",
        "fourbar",
        "thermalexpansion",
        "dimensional",
        "projectile",
        "heattransfer",
        "insulation",
        "heatexchanger",
        "heatpump",
        "refrigeration",
        "psychrometrics",
        "thermocouple",
        "thermistor",
        "thermocycle",
        "fanlaws",
        // The 18 fluid-mechanics / hydraulics / thermo-fluids / electronics
        // workbenches wired in this change.
        "pump",
        "pipeflow",
        "pipenetwork",
        "hydraulics",
        "pneumatics",
        "fluidstatics",
        "openchannel",
        "weir",
        "orifice",
        "combustion",
        "diffusion",
        "marine",
        "resistornetwork",
        "capacitor",
        "opamp",
        "bjt",
        "mosfet",
        "rectifier",
    ];

    #[test]
    fn every_wired_workbench_kind_resolves_and_builds_a_tri3_mesh() {
        // Each wired engineering-workbench kind resolves to a registry entry
        // whose pure builder yields a non-empty triangle-mesh product — so a
        // missing table line (or a builder that stops producing geometry) is
        // caught here without the file-poll bridge plumbing.
        for &k in WIRED_WORKBENCH_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(
                entry.kind, k,
                "entry.kind echoes the looked-up kind for {k}"
            );
            let product = (entry.build)();
            let mesh = product
                .mesh
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: a 3-D product carries a mesh"));
            assert!(!mesh.mesh.nodes.is_empty(), "{k}: mesh has vertices");
            assert!(mesh.mesh.total_elements() > 0, "{k}: mesh has triangles");
            // A title and at least one headline result row from the tool's own
            // compute() — the tile is more than bare geometry.
            assert!(!product.title.is_empty(), "{k}: product has a title");
            assert!(
                !product.lines.is_empty(),
                "{k}: product carries result rows"
            );
        }
    }

    #[test]
    fn registry_resolves_all_known_3d_kinds() {
        // Every migrated kind resolves to an entry whose `kind` echoes the query
        // (so the table and the `kind_static` spelling stay in lockstep).
        for &k in KNOWN_3D_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(entry.kind, k, "entry.kind echoes the looked-up kind");
        }
    }

    #[test]
    fn registry_builds_a_live_mesh_for_each_3d_kind() {
        // Each builder is pure and yields a non-empty 3-D mesh product (the FEM
        // one additionally carries the von-Mises vertex colours). This exercises
        // the builders the reducer dispatches to, so a regression in any of them
        // is caught here without the file-poll plumbing.
        for &k in KNOWN_3D_KINDS {
            let entry = lookup(k).unwrap();
            let product = (entry.build)();
            let mesh = product
                .mesh
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: a 3-D product carries a mesh"));
            assert!(!mesh.mesh.nodes.is_empty(), "{k}: mesh has vertices");
            assert!(mesh.mesh.total_elements() > 0, "{k}: mesh has triangles");
        }
        // The FEM cantilever is the only kind that ships per-vertex colours.
        assert!(
            (lookup("fem").unwrap().build)().vertex_colors.is_some(),
            "fem product carries von-Mises vertex colours"
        );
    }

    #[test]
    fn registry_returns_none_for_an_unknown_kind() {
        // An unknown kind resolves to None → the reducer skips it safely.
        assert!(lookup("not-a-model").is_none());
        assert!(lookup("").is_none());
        // `dna` is a text card / 2-D drawing, NOT a registry 3-D mesh kind.
        assert!(lookup("dna").is_none());
    }
}

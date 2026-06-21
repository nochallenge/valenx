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
//!         image: None, image_texture: None, animation: None,
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

/// Promote a [`valenx_viz::TriangleMesh`] triangle soup into a
/// [`valenx_mesh::Mesh`] of one `Tri3` element **per source triangle, in the
/// source order**, tagged `id`.
///
/// Used by the molecular-geometry product builders (`molecule` / `reactdyn`),
/// whose canonical geometry comes out of
/// [`crate::genetics::molecule_view::ball_and_stick`] as a `TriangleMesh`
/// rather than a `valenx_mesh::Mesh`. Keeping one element per source triangle
/// **in order** is what lets a paired per-triangle colour list expand cleanly to
/// the triangle-major per-vertex `vertex_colors` the tile renderer expects:
/// [`crate::viewport::mesh_to_triangle_surface`] re-emits this mesh's `Tri3`
/// elements in the same order, so triangle *k* of the returned mesh is triangle
/// *k* of the input (see [`per_triangle_to_vertex_colors`]). Vertices are not
/// welded — each triangle owns three fresh nodes — which matches how the source
/// soup is laid out and keeps the index-alignment trivial.
pub(crate) fn mesh_from_triangle_soup(
    soup: &valenx_viz::TriangleMesh,
    id: &str,
) -> valenx_mesh::Mesh {
    let mut nodes: Vec<nalgebra::Vector3<f64>> = Vec::with_capacity(soup.triangles.len() * 3);
    let mut conn: Vec<u32> = Vec::with_capacity(soup.triangles.len() * 3);
    for tri in &soup.triangles {
        for v in &tri.vertices {
            conn.push(nodes.len() as u32);
            nodes.push(nalgebra::Vector3::new(
                v[0] as f64,
                v[1] as f64,
                v[2] as f64,
            ));
        }
    }
    let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
    block.connectivity = conn;
    let mut mesh = valenx_mesh::Mesh::new(id);
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// Expand a **per-triangle** colour list (one `[r, g, b]` per triangle, in the
/// triangle order of a [`mesh_from_triangle_soup`] mesh) into the triangle-major
/// **per-vertex** `vertex_colors` a [`WorkspaceProduct`] carries: each
/// triangle's colour is repeated three times (once per corner), in the same
/// order [`crate::wgpu_renderer::triangles_to_vertices`] walks the surface.
///
/// The molecular builders pair this with [`mesh_from_triangle_soup`] (one
/// element per source triangle, in order) so the returned vec lines up 1:1 with
/// the renderer's emitted surface vertices.
pub(crate) fn per_triangle_to_vertex_colors(per_tri: &[[f32; 3]]) -> Vec<[f32; 3]> {
    let mut out = Vec::with_capacity(per_tri.len() * 3);
    for &c in per_tri {
        out.push(c);
        out.push(c);
        out.push(c);
    }
    out
}

/// Build the triangle-major **per-vertex** `vertex_colors` for a `Tri3`
/// `valenx_mesh::Mesh` whose nodes carry a scalar field (one value per node),
/// mapped through the shared `valenx_fields` cool-to-warm divergent ramp over
/// `[min, max]`.
///
/// Walks the mesh's `Tri3` elements in the exact order
/// [`crate::viewport::mesh_to_triangle_surface`] re-emits them (block-major,
/// then the three corners of each triangle) so the result is index-aligned with
/// the surface vertices [`crate::wgpu_renderer::triangles_to_vertices_colored`]
/// draws — i.e. its length is `3 × (number of Tri3 elements)`. A node index past
/// the field length, or a non-finite value, degrades that corner to mid-ramp
/// rather than panicking. Used by the aero product (per-face `Cp` painted on the
/// voxelized body shell).
pub(crate) fn node_field_to_vertex_colors(
    mesh: &valenx_mesh::Mesh,
    field: &[f64],
    min: f64,
    max: f64,
) -> Vec<[f32; 3]> {
    let to_rgb = |node: u32| -> [f32; 3] {
        let v = field.get(node as usize).copied().unwrap_or(f64::NAN);
        let v = if v.is_finite() { v } else { 0.5 * (min + max) };
        let [r, g, b] = valenx_fields::colormap::cool_to_warm_in_range(v, min, max);
        [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0]
    };
    let mut out = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type == valenx_mesh::ElementType::Tri3 {
            for f in block.connectivity.chunks_exact(3) {
                out.push(to_rgb(f[0]));
                out.push(to_rgb(f[1]));
                out.push(to_rgb(f[2]));
            }
        }
    }
    out
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
/// **This `match` is the single shared edit point for bridge product tools**:
/// each arm is one line pairing a wire `kind` with the per-tool builder in that
/// tool's own module. The substantive per-tool code lives in those builders,
/// not here — so adding a kind is a one-line addition (see the module docs for
/// the copy-paste pattern). Most arms are 3-D mesh tools; the trailing block is
/// the **DATA-ONLY** tools whose builders return a mesh-less *text-card*
/// `WorkspaceProduct` (`mesh: None`, populated `lines`) — the reducer dispatches
/// those through this same table and the tile renders the card (the mesh-less
/// path the `dna` card uses) instead of a 3-D view. Note `dna` itself is
/// intentionally absent: `show_3d:dna` is a text card handled directly in the
/// reducer, and the 2-D `show_2d` drawings (`rcbeam` / `dna`) have their own
/// separate path.
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
        // Electrical / EM / power-systems / electrochemistry / photonics /
        // acoustics / nuclear / signals / propulsion family (each builder lives
        // in its own workbench module; see that module's `*_product`).
        "filter" => crate::filter_workbench::filter_product,
        "antenna" => crate::antenna_workbench::antenna_product,
        "transmissionline" => crate::transmissionline_workbench::transmissionline_product,
        "coil" => crate::coil_workbench::coil_product,
        "led" => crate::led_workbench::led_product,
        "transformer" => crate::transformer_workbench::transformer_product,
        "threephase" => crate::threephase_workbench::threephase_product,
        "powerfactor" => crate::powerfactor_workbench::powerfactor_product,
        "electrochem" => crate::electrochem_workbench::electrochem_product,
        "batterypack" => crate::batterypack_workbench::batterypack_product,
        "batteryecm" => crate::batteryecm_workbench::batteryecm_product,
        "solarpv" => crate::solarpv_workbench::solarpv_product,
        "optics" => crate::optics_workbench::optics_product,
        "acoustics" => crate::acoustics_workbench::acoustics_product,
        "radioactivity" => crate::radioactivity_workbench::radioactivity_product,
        "queueing" => crate::queueing_workbench::queueing_product,
        "fft" => crate::fft_workbench::fft_product,
        "engine" => crate::engine_workbench::engine_product,
        // Aerospace + bio family — the last of the 3-D mesh workbenches
        // (each builder lives in its own workbench module; see that module's
        // `*_product`).
        "fixedwing" => crate::fixedwing_workbench::fixedwing_product,
        "drone" => crate::drone_workbench::drone_product,
        "windturbine" => crate::windturbine_workbench::windturbine_product,
        "pharmacokinetics" => crate::pharmacokinetics_workbench::pharmacokinetics_product,
        "enzymekinetics" => crate::enzymekinetics_workbench::enzymekinetics_product,
        "hemodynamics" => crate::hemodynamics_workbench::hemodynamics_product,
        "bonemech" => crate::bonemech_workbench::bonemech_product,
        "bmr" => crate::bmr_workbench::bmr_product,
        "thermoreg" => crate::thermoreg_workbench::thermoreg_product,
        "osmosis" => crate::osmosis_workbench::osmosis_product,
        "acidbase" => crate::acidbase_workbench::acidbase_product,
        "popdynamics" => crate::popdynamics_workbench::popdynamics_product,
        // CAD / reconstruction / reinforcement / molecular / sheet-metal / aero
        // family — workbenches that produce real geometry from their own kernels
        // (not the mechanical `*_solid_mesh` pattern); each builder lives in its
        // own module and extracts a `valenx_mesh::Mesh` from that kernel. The
        // molecular (`molecule` / `reactdyn`) and `aero` products additionally
        // carry per-vertex colours (CPK element palette / Cp field).
        "cad" => crate::cad_workbench::cad_product,
        "reverse" => crate::reverse_workbench::reverse_product,
        "reinforcement" => crate::reinforcement_workbench::reinforcement_product,
        "molecule" => crate::genetics::molecule_view::molecule_product,
        "reactdyn" => crate::reactdyn_workbench::reactdyn_product,
        "sheetmetal" => crate::sheetmetal_workbench::sheetmetal_product,
        "aero" => crate::aero_workbench::aero_product,
        // DATA-ONLY workbenches — their output is numbers / analysis, not 3-D
        // geometry, so each builder returns a TEXT-CARD `WorkspaceProduct`
        // (`mesh: None`, populated `lines`). The reducer dispatches these
        // through the same `lookup` as the mesh kinds, and the tile renders the
        // card (the mesh-less path the `dna` card uses) rather than a 3-D view.
        // Each builder lives in its own workbench module and formats the genuine
        // computed result rows from that tool's canonical default state.
        "fields" => crate::fields_workbench::fields_product,
        "fasteners" => crate::fasteners_workbench::fasteners_product,
        "frames" => crate::frames_workbench::frames_product,
        "collision" => crate::collision_workbench::collision_product,
        "geomatics" => crate::geomatics_workbench::geomatics_product,
        "hvac" => crate::hvac_workbench::hvac_product,
        "gasdynamics" => crate::gasdynamics_workbench::gasdynamics_product,
        "piping" => crate::piping_workbench::piping_product,
        "cfd" => crate::cfd_workbench::cfd_product,
        "astro" => crate::astro_workbench::astro_product,
        "car" => crate::car_workbench::car_product,
        "neuro" => crate::neuro_workbench::neuro_product,
        "variant_effect" => crate::variant_effect_workbench::variant_effect_product,
        // IMAGE / 2-D-DRAWING / extra-card workbenches — these resolve through
        // the same `lookup`, but their builders return a product that is
        // neither a `Tri3` mesh nor a bare text card:
        //   • `render`  → an IMAGE product (`image: Some(ColorImage)`): a small
        //     synchronous path-trace of a Cornell box (the pane's image branch
        //     uploads + draws it).
        //   • `animate` → a DATA-ONLY text CARD (`mesh: None`, `lines`): the
        //     keyframe-timeline summary (no posed body to rasterise — see
        //     `animate_product`).
        //   • `draft2d` / `interior` → 2-D DRAWING products
        //     (`kind2d: Some(..)`): the 2-D CAD drawing / floor plan, painted by
        //     the same egui 2-D branch as `rcbeam` / `dna`.
        // They are asserted by their own registry tests (image.is_some() /
        // kind2d.is_some()), not the Tri3-mesh assertion.
        "render" => crate::render_workbench::render_product,
        "animate" => crate::animate_workbench::animate_product,
        "draft2d" => crate::draft2d_workbench::draft2d_product,
        "interior" => crate::interior_workbench::interior_product,
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
        "filter" => "filter",
        "antenna" => "antenna",
        "transmissionline" => "transmissionline",
        "coil" => "coil",
        "led" => "led",
        "transformer" => "transformer",
        "threephase" => "threephase",
        "powerfactor" => "powerfactor",
        "electrochem" => "electrochem",
        "batterypack" => "batterypack",
        "batteryecm" => "batteryecm",
        "solarpv" => "solarpv",
        "optics" => "optics",
        "acoustics" => "acoustics",
        "radioactivity" => "radioactivity",
        "queueing" => "queueing",
        "fft" => "fft",
        "engine" => "engine",
        // Aerospace + bio family — the last of the 3-D mesh workbenches.
        "fixedwing" => "fixedwing",
        "drone" => "drone",
        "windturbine" => "windturbine",
        "pharmacokinetics" => "pharmacokinetics",
        "enzymekinetics" => "enzymekinetics",
        "hemodynamics" => "hemodynamics",
        "bonemech" => "bonemech",
        "bmr" => "bmr",
        "thermoreg" => "thermoreg",
        "osmosis" => "osmosis",
        "acidbase" => "acidbase",
        "popdynamics" => "popdynamics",
        // CAD / reconstruction / reinforcement / molecular / sheet-metal / aero
        // family.
        "cad" => "cad",
        "reverse" => "reverse",
        "reinforcement" => "reinforcement",
        "molecule" => "molecule",
        "reactdyn" => "reactdyn",
        "sheetmetal" => "sheetmetal",
        "aero" => "aero",
        // DATA-ONLY card workbenches (mesh-less text-card products).
        "fields" => "fields",
        "fasteners" => "fasteners",
        "frames" => "frames",
        "collision" => "collision",
        "geomatics" => "geomatics",
        "hvac" => "hvac",
        "gasdynamics" => "gasdynamics",
        "piping" => "piping",
        "cfd" => "cfd",
        "astro" => "astro",
        "car" => "car",
        "neuro" => "neuro",
        "variant_effect" => "variant_effect",
        // IMAGE / 2-D-drawing / extra-card workbenches (see `lookup`).
        "render" => "render",
        "animate" => "animate",
        "draft2d" => "draft2d",
        "interior" => "interior",
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
        // The 18 electrical / EM / power-systems / electrochemistry /
        // photonics / acoustics / nuclear / signals / propulsion workbenches
        // wired in this change.
        "filter",
        "antenna",
        "transmissionline",
        "coil",
        "led",
        "transformer",
        "threephase",
        "powerfactor",
        "electrochem",
        "batterypack",
        "batteryecm",
        "solarpv",
        "optics",
        "acoustics",
        "radioactivity",
        "queueing",
        "fft",
        "engine",
        // The 12 aerospace + bio workbenches — the last of the 3-D mesh
        // tools — wired in this change.
        "fixedwing",
        "drone",
        "windturbine",
        "pharmacokinetics",
        "enzymekinetics",
        "hemodynamics",
        "bonemech",
        "bmr",
        "thermoreg",
        "osmosis",
        "acidbase",
        "popdynamics",
        // The 7 real-geometry-extraction workbenches wired in this change — CAD
        // CSG, point-cloud reconstruction, rebar cage, coloured molecule,
        // reaction-dynamics frame, folded sheet metal, aero Cp on a demo body.
        "cad",
        "reverse",
        "reinforcement",
        "molecule",
        "reactdyn",
        "sheetmetal",
        "aero",
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
        // The 18 electrical / EM / power-systems / electrochemistry /
        // photonics / acoustics / nuclear / signals / propulsion workbenches
        // wired in this change.
        "filter",
        "antenna",
        "transmissionline",
        "coil",
        "led",
        "transformer",
        "threephase",
        "powerfactor",
        "electrochem",
        "batterypack",
        "batteryecm",
        "solarpv",
        "optics",
        "acoustics",
        "radioactivity",
        "queueing",
        "fft",
        "engine",
        // The 12 aerospace + bio workbenches — the last of the 3-D mesh
        // tools — wired in this change.
        "fixedwing",
        "drone",
        "windturbine",
        "pharmacokinetics",
        "enzymekinetics",
        "hemodynamics",
        "bonemech",
        "bmr",
        "thermoreg",
        "osmosis",
        "acidbase",
        "popdynamics",
        // The 7 real-geometry-extraction workbenches wired in this change. Each
        // must resolve and build a non-empty `Tri3` product with a title and at
        // least one readout row — `molecule` / `reactdyn` / `aero` additionally
        // carry per-vertex colours (asserted separately below).
        "cad",
        "reverse",
        "reinforcement",
        "molecule",
        "reactdyn",
        "sheetmetal",
        "aero",
    ];

    /// The DATA-ONLY card workbench kinds wired into the bridge. These differ
    /// from [`KNOWN_3D_KINDS`] / [`WIRED_WORKBENCH_KINDS`] in one essential way:
    /// their output is numbers / analysis, **not** 3-D geometry, so each builder
    /// returns a *text-card* [`WorkspaceProduct`] (`mesh: None`, populated
    /// `lines`) rather than a `Tri3` mesh. They are therefore asserted by the
    /// card-specific [`every_card_kind_resolves_and_builds_a_text_card`] test
    /// below — which checks `mesh.is_none()` + non-empty `lines` — and are kept
    /// out of the mesh arrays so the Tri3-mesh assertions never run against them.
    const CARD_KINDS: &[&str] = &[
        "fields",
        "fasteners",
        "frames",
        "collision",
        "geomatics",
        "hvac",
        "gasdynamics",
        "piping",
        "cfd",
        "astro",
        "car",
        "neuro",
        "variant_effect",
        // The animation timeline is a DATA-ONLY card too — it summarises the
        // keyframe timeline (count / duration / sampled joint values); there is
        // no posed body to rasterise, so it carries no mesh and no 2-D drawing.
        "animate",
    ];

    /// The IMAGE workbench kinds — their builder returns a product carrying a
    /// raster [`crate::WorkspaceProduct::image`] (a CPU `egui::ColorImage`)
    /// instead of a mesh / card / 2-D drawing. Asserted by
    /// [`image_kinds_resolve_and_build_an_image_product`] (`image.is_some()` +
    /// `mesh.is_none()`), so they are kept out of the mesh / card / 2-D arrays.
    const IMAGE_KINDS: &[&str] = &["render"];

    /// The 2-D-DRAWING workbench kinds whose builder returns a product with a
    /// [`crate::WorkspaceProduct::kind2d`] (`Some(..)`) painted by the egui 2-D
    /// branch — the same path `rcbeam` / `dna` use, but resolved through the
    /// registry. Asserted by
    /// [`drawing_2d_kinds_resolve_and_build_a_2d_drawing`] (`kind2d.is_some()` +
    /// `mesh.is_none()` + `image.is_none()`).
    const DRAWING_2D_KINDS: &[&str] = &["draft2d", "interior"];

    #[test]
    fn every_card_kind_resolves_and_builds_a_text_card() {
        // Each DATA-ONLY workbench kind resolves to a registry entry whose pure
        // builder yields a *text-card* product: NO mesh (so the tile takes the
        // mesh-less card path the `dna` card uses, not the 3-D viewport), a title,
        // and at least one genuine computed result row. This is the card-kind
        // counterpart to the Tri3-mesh assertion above — the distinguishing check
        // is `mesh.is_none()`, which is exactly what separates a card kind from a
        // mesh kind, so a card builder that accidentally grew a mesh (or a mesh
        // builder mis-listed here) is caught.
        for &k in CARD_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(
                entry.kind, k,
                "entry.kind echoes the looked-up kind for {k}"
            );
            let product = (entry.build)();
            assert!(
                product.mesh.is_none(),
                "{k}: a DATA-ONLY card carries no mesh"
            );
            // No 2-D drawing either — these render as the plain text card, not the
            // egui 2-D-drawing branch (which is reserved for `rcbeam` / `dna` /
            // `draft2d` / `interior`).
            assert!(
                product.kind2d.is_none(),
                "{k}: a DATA-ONLY card is a text card, not a 2-D drawing"
            );
            // And no raster image — a text card is text, not the IMAGE branch
            // (which is `render`). Guards a card builder that mis-listed an
            // image kind here.
            assert!(
                product.image.is_none(),
                "{k}: a DATA-ONLY card is a text card, not an image"
            );
            assert!(!product.title.is_empty(), "{k}: product has a title");
            assert!(
                !product.lines.is_empty(),
                "{k}: card carries genuine result rows"
            );
        }
    }

    #[test]
    fn image_kinds_resolve_and_build_an_image_product() {
        // Each IMAGE workbench kind resolves to a registry entry whose builder
        // yields a product carrying a non-empty raster `image` (a CPU
        // `egui::ColorImage`) and NO mesh — so the tile takes the image branch,
        // not the 3-D viewport / card / 2-D drawing. The image's pixel buffer
        // is sized `w·h` Color32 entries (the `from_rgb` invariant), which we
        // assert is non-empty and matches its declared size.
        for &k in IMAGE_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(
                entry.kind, k,
                "entry.kind echoes the looked-up kind for {k}"
            );
            let product = (entry.build)();
            assert!(
                product.mesh.is_none(),
                "{k}: an image product carries no mesh"
            );
            assert!(
                product.kind2d.is_none(),
                "{k}: an image product is not a 2-D drawing"
            );
            let image = product
                .image
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: an image product carries an image"));
            let [w, h] = image.size;
            assert!(w > 0 && h > 0, "{k}: image has non-zero dimensions");
            assert_eq!(
                image.pixels.len(),
                w * h,
                "{k}: ColorImage pixel count equals w·h"
            );
            assert!(!product.title.is_empty(), "{k}: product has a title");
        }
    }

    #[test]
    fn drawing_2d_kinds_resolve_and_build_a_2d_drawing() {
        // Each 2-D-DRAWING workbench kind resolves to a registry entry whose
        // builder yields a product with a `kind2d` (the egui 2-D branch paints
        // it) and NO mesh / NO image — the distinguishing check is
        // `kind2d.is_some()`, exactly what routes it to the 2-D painter rather
        // than the 3-D viewport / image branch / text card. The drawing's view
        // must carry geometry (a non-empty entity / room list) so the painter
        // has something to draw.
        for &k in DRAWING_2D_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(
                entry.kind, k,
                "entry.kind echoes the looked-up kind for {k}"
            );
            let product = (entry.build)();
            assert!(product.mesh.is_none(), "{k}: a 2-D drawing carries no mesh");
            assert!(
                product.image.is_none(),
                "{k}: a 2-D drawing carries no image"
            );
            match product
                .kind2d
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: a 2-D drawing carries a kind2d"))
            {
                crate::Workspace2dKind::Draft2d(view) => {
                    assert!(!view.entities.is_empty(), "{k}: drawing has entities");
                }
                crate::Workspace2dKind::FloorPlan(plan) => {
                    assert!(!plan.rooms.is_empty(), "{k}: floor plan has a room");
                }
                other => panic!("{k}: unexpected 2-D kind {other:?}"),
            }
            assert!(!product.title.is_empty(), "{k}: product has a title");
            assert!(
                !product.lines.is_empty(),
                "{k}: drawing carries readout rows"
            );
        }
    }

    #[test]
    fn card_kinds_are_disjoint_from_mesh_kinds() {
        // A kind is exactly one of: a 3-D mesh kind, a DATA-ONLY card kind, an
        // IMAGE kind, or a 2-D-DRAWING kind — never two at once, so each
        // test-suite's invariant (mesh-Some / mesh-None+card / image / kind2d)
        // can't collide. (Guards against a future edit that lists a kind in two
        // of these arrays.)
        let mesh = |c: &str| KNOWN_3D_KINDS.contains(&c) || WIRED_WORKBENCH_KINDS.contains(&c);
        for &c in CARD_KINDS {
            assert!(!mesh(c), "{c} is a card kind and must not be a mesh kind");
            assert!(
                !IMAGE_KINDS.contains(&c),
                "{c} is a card kind and must not be an image kind"
            );
            assert!(
                !DRAWING_2D_KINDS.contains(&c),
                "{c} is a card kind and must not be a 2-D-drawing kind"
            );
        }
        for &c in IMAGE_KINDS {
            assert!(!mesh(c), "{c} is an image kind and must not be a mesh kind");
            assert!(
                !DRAWING_2D_KINDS.contains(&c),
                "{c} is an image kind and must not be a 2-D-drawing kind"
            );
        }
        for &c in DRAWING_2D_KINDS {
            assert!(
                !mesh(c),
                "{c} is a 2-D-drawing kind and must not be a mesh kind"
            );
        }
    }

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
        // The FEM cantilever ships von-Mises per-vertex colours.
        assert!(
            (lookup("fem").unwrap().build)().vertex_colors.is_some(),
            "fem product carries von-Mises vertex colours"
        );
    }

    #[test]
    fn coloured_products_carry_a_vertex_aligned_colour_vec() {
        // The molecular (CPK by element) and aero (Cp field) products ship
        // per-vertex colours. The tile renderer only takes the coloured path
        // when the colour vec length EXACTLY equals the emitted surface-vertex
        // stream — three per `Tri3` element (triangle-major, vertex-within-
        // triangle). Assert that invariant here so a future mesh/colour drift
        // silently falls back to grey rather than mis-indexing (and this catches
        // the alignment without the GPU).
        for k in ["molecule", "reactdyn", "aero"] {
            let product = (lookup(k).unwrap().build)();
            let mesh = product
                .mesh
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: a 3-D product carries a mesh"));
            let colors = product
                .vertex_colors
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: coloured product carries vertex_colors"));
            let expected = mesh.mesh.total_elements() * 3;
            assert_eq!(
                colors.len(),
                expected,
                "{k}: vertex_colors length ({}) must equal 3 \u{00D7} triangle count ({expected}) \
                 so the renderer takes the coloured path",
                colors.len(),
            );
            // Every channel is a sane [0, 1] colour.
            for c in colors {
                for ch in c {
                    assert!(
                        ch.is_finite() && (0.0..=1.0).contains(ch),
                        "{k}: colour channel {ch} out of [0, 1]"
                    );
                }
            }
        }
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

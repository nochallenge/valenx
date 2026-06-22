//! **Per-product base material colour.**
//!
//! Most catalogue products build a plain `Tri3` surface mesh and leave
//! [`WorkspaceProduct::vertex_colors`](crate::WorkspaceProduct) `None`, so the
//! tile renderer paints them in the single neutral brushed-metal constant
//! ([`crate::wgpu_renderer::METAL_BASE`]) — every rocket, beam, motor and pump
//! comes out the same flat grey. A few rebuilt products (the machinery /
//! aero / marine / fasteners families, plus the FEM stress map and the
//! molecular CPK / aero-Cp overlays) already author *per-part* colours and so
//! read with material variety; the rest do not.
//!
//! This module closes that gap centrally: [`base_color_for`] maps a product
//! `kind` to one sensible, mid-saturation **base** material tone by category
//! (steel for the machine-design family, concrete for structural / civil,
//! copper for electronics, a warm tone for thermal, aluminium for aerospace,
//! a teal organic tone for bio / chem, a fluid blue-grey for hydraulics /
//! marine, and a pleasant neutral default for everything else, including
//! unknown kinds).
//!
//! It is wired into [`materialize_pending`](crate::agent_commands) — the single
//! place every product is built and read out — which, *after* the build, fills
//! a uniform base-colour `vertex_colors` for any product that has a `mesh` but
//! authored **no** colours of its own. Products that already set
//! `vertex_colors` (the rebuilt per-part-coloured ones, the FEM von-Mises map,
//! the molecule / aero overlays) are left untouched, so they keep their richer
//! shading; card / 2-D / image products (`mesh: None`) are never touched.
//!
//! ## Why a base colour and not a real material model
//!
//! The renderer's coloured path ([`crate::wgpu_renderer::triangles_to_vertices_colored`])
//! takes a flat `[r, g, b]` albedo per surface vertex. A single base colour
//! repeated across the whole product is therefore the smallest change that
//! turns the uniform-grey catalogue into a varied, professional-looking one
//! without disturbing the products that already colour their own parts. Real
//! per-part materials remain each producer's job (see the gearbox / aero
//! builders); this is the sensible *floor*.

/// Broad material category a product belongs to, used to pick its
/// [`base_color_for`] tone. Each variant maps to one professional, mid-
/// saturation base colour (see [`Category::base_color`]); the mapping from a
/// product `kind` to a category lives in [`category_for`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    /// Machine-design / mechanical hardware (gears, bearings, shafts, springs,
    /// fasteners, motors, mechanisms): brushed structural steel.
    Mechanical,
    /// Structural / civil members (beams, trusses, columns, walls, reinforced
    /// concrete, brackets, pressure vessels): a light concrete grey.
    Structural,
    /// Electrical / electronics (passives, semiconductors, antennas, coils,
    /// transformers, batteries, PV, power systems): copper.
    Electronics,
    /// Thermal / thermodynamic / HVAC / combustion equipment: a warm
    /// terracotta tone.
    Thermal,
    /// Aerospace airframes (rockets, fixed-wing, drones, turbines, aero/astro
    /// bodies): light aluminium / off-white.
    Aerospace,
    /// Bio / chemistry / physiology (molecules, kinetics, haemodynamics, bone,
    /// metabolism, neuro): a soft organic teal.
    BioChem,
    /// Fluid mechanics / hydraulics / marine / piping: a wet blue-grey.
    Fluid,
    /// Anything not otherwise categorised (CAD parts, reverse-engineered
    /// scans, sheet metal, and any unknown kind): a pleasant neutral grey.
    Neutral,
}

impl Category {
    /// The base material colour for this category, as a linear `[r, g, b]` in
    /// `[0, 1]` — the same space [`crate::wgpu_renderer::triangles_to_vertices_colored`]
    /// consumes. Deliberately mid-saturation and a touch desaturated so the lit
    /// render reads as a real engineering material, never neon.
    pub const fn base_color(self) -> [f32; 3] {
        match self {
            // Brushed structural steel — cool, neutral, slightly blue.
            Category::Mechanical => [0.62, 0.64, 0.68],
            // Light cured concrete — warm-neutral grey.
            Category::Structural => [0.70, 0.70, 0.67],
            // Copper — warm metallic orange-brown (PCB / windings / busbar).
            Category::Electronics => [0.72, 0.45, 0.20],
            // Warm terracotta — heat / thermal equipment.
            Category::Thermal => [0.78, 0.42, 0.30],
            // Aluminium / off-white — aerospace skin.
            Category::Aerospace => [0.82, 0.84, 0.86],
            // Soft organic teal — bio / chemistry.
            Category::BioChem => [0.30, 0.62, 0.58],
            // Wet blue-grey — fluids / marine / piping.
            Category::Fluid => [0.40, 0.54, 0.62],
            // Pleasant neutral grey — default / CAD / unknown.
            Category::Neutral => [0.66, 0.66, 0.68],
        }
    }
}

/// Map a product `kind` (the agent-bridge `show_3d` wire string, identical to a
/// [`crate::products_registry`] table key) to its material [`Category`].
///
/// Grouped by the same families as the registry table. An unknown kind — or any
/// kind whose honest material is "just a part" (CAD / scans / sheet metal) —
/// falls through to [`Category::Neutral`]. Products that author their own
/// per-part colours (`fem`, `molecule`, `reactdyn`, `aero`, `gearbox`,
/// `bearing`, …) still get a category here, but it is never applied to them:
/// [`materialize_pending`](crate::agent_commands) only fills a base colour when
/// the product left `vertex_colors` `None`.
pub fn category_for(kind: &str) -> Category {
    match kind {
        // ---- Machine-design / mechanical hardware → steel ----
        "gear" | "geartooth" | "gearbox" | "bearing" | "clutch" | "brake" | "pulley"
        | "flywheel" | "leadscrew" | "screwthread" | "shaftdesign" | "camdynamics"
        | "conveyor" | "bolt" | "rivet" | "fasteners" | "springs" | "springdesign"
        | "springcombination" | "beltdrive" | "chaindrive" | "fourbar" | "leverage"
        | "inclinedplane" | "vibration" | "dcmotor" | "inductionmotor" | "engine" => {
            Category::Mechanical
        }

        // ---- Structural / civil members → concrete grey ----
        "beam" | "truss" | "plate" | "buckling" | "columnsteel" | "retainingwall"
        | "soilbearing" | "statics" | "mohr" | "torsion" | "fatigue" | "fracture"
        | "creep" | "pressurevessel" | "straingauge" | "strainrosette" | "rail"
        | "rcbeam" | "reinforcement" | "bracket" | "frames" => Category::Structural,

        // ---- Electrical / electronics → copper ----
        "resistornetwork" | "capacitor" | "opamp" | "bjt" | "mosfet" | "rectifier"
        | "filter" | "antenna" | "transmissionline" | "coil" | "led" | "transformer"
        | "threephase" | "powerfactor" | "electrochem" | "batterypack" | "batteryecm"
        | "solarpv" | "fft" => Category::Electronics,

        // ---- Thermal / thermodynamics / HVAC / combustion → warm tone ----
        "heattransfer" | "insulation" | "heatexchanger" | "heatpump" | "refrigeration"
        | "psychrometrics" | "thermocouple" | "thermistor" | "thermocycle" | "fanlaws"
        | "thermalexpansion" | "combustion" | "hvac" => Category::Thermal,

        // ---- Aerospace airframes → aluminium / white ----
        "rocket" | "fixedwing" | "drone" | "windturbine" | "aero" | "astro"
        | "projectile" => Category::Aerospace,

        // ---- Bio / chemistry / physiology → organic teal ----
        "molecule" | "reactdyn" | "pharmacokinetics" | "enzymekinetics" | "hemodynamics"
        | "bonemech" | "bmr" | "thermoreg" | "osmosis" | "acidbase" | "popdynamics"
        | "neuro" | "variant_effect" | "diffusion" => Category::BioChem,

        // ---- Fluid mechanics / hydraulics / marine / piping → blue-grey ----
        "pump" | "pipeflow" | "pipenetwork" | "hydraulics" | "pneumatics" | "fluidstatics"
        | "openchannel" | "weir" | "orifice" | "marine" | "piping" | "cfd" | "gasdynamics" => {
            Category::Fluid
        }

        // ---- Everything else (CAD parts, scans, sheet metal, car, fem,
        //      fields, collision, geomatics, …) + unknown → neutral ----
        _ => Category::Neutral,
    }
}

/// The base material colour for a product `kind`, as a linear `[r, g, b]` in
/// `[0, 1]`. A thin convenience over [`category_for`] + [`Category::base_color`];
/// the one entry point [`materialize_pending`](crate::agent_commands) calls.
///
/// Unknown kinds resolve to the neutral default rather than an error, so a new
/// catalogue kind always renders in a sensible tone even before it gets its own
/// category line.
pub fn base_color_for(kind: &str) -> [f32; 3] {
    category_for(kind).base_color()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_categories_get_distinct_base_colours() {
        // A representative kind from each mandated family resolves to its own
        // category and a distinct base colour, so the catalogue actually reads
        // with material variety instead of one uniform grey.
        let mech = base_color_for("gear");
        let elec = base_color_for("capacitor");
        let struc = base_color_for("beam");
        let unknown = base_color_for("not-a-kind");

        assert_eq!(category_for("gear"), Category::Mechanical);
        assert_eq!(category_for("capacitor"), Category::Electronics);
        assert_eq!(category_for("beam"), Category::Structural);
        assert_eq!(category_for("not-a-kind"), Category::Neutral);

        // The three named families are mutually distinct …
        assert_ne!(mech, elec, "mechanical vs electrical differ");
        assert_ne!(mech, struc, "mechanical vs structural differ");
        assert_ne!(elec, struc, "electrical vs structural differ");
        // … and none of them collides with the neutral default.
        assert_ne!(mech, unknown, "mechanical is not the neutral default");
        assert_ne!(elec, unknown, "electrical is not the neutral default");
        assert_ne!(struc, unknown, "structural is not the neutral default");
    }

    #[test]
    fn unknown_kind_falls_back_to_the_neutral_default() {
        // Any unrecognised kind (and the empty string) is the neutral grey, never
        // a panic — a brand-new catalogue kind still renders in a sensible tone.
        assert_eq!(base_color_for("not-a-kind"), Category::Neutral.base_color());
        assert_eq!(base_color_for(""), Category::Neutral.base_color());
    }

    #[test]
    fn every_base_colour_is_a_sane_professional_albedo() {
        // Every category's base colour is a finite [0, 1] albedo (the renderer's
        // colour space) and not pure white/black, so nothing reads as a neon or
        // clipped material.
        for cat in [
            Category::Mechanical,
            Category::Structural,
            Category::Electronics,
            Category::Thermal,
            Category::Aerospace,
            Category::BioChem,
            Category::Fluid,
            Category::Neutral,
        ] {
            let c = cat.base_color();
            for ch in c {
                assert!(
                    ch.is_finite() && (0.0..=1.0).contains(&ch),
                    "{cat:?} channel {ch} out of [0, 1]"
                );
            }
            assert!(
                c.iter().any(|&ch| ch > 0.05),
                "{cat:?} is not pure black"
            );
            assert!(
                c.iter().any(|&ch| ch < 0.95),
                "{cat:?} is not pure white"
            );
        }
    }
}

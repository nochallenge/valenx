//! Generic AI-drivability guard: **every interactive widget on every workbench
//! panel carries an accessible name.**
//!
//! Closing AI-drivability gap 2. The per-panel headless tests (and the older
//! `uq_workbench` test) only blanket-verified that numeric `DragValue`s
//! (`Role::SpinButton`) were named. This module generalises that to the full
//! interactive widget set — buttons, combo boxes, checkboxes, sliders, text
//! inputs, toggle/selectable labels and radio buttons — across **all** of the
//! `show_*_workbench` panels, so an AI driver (or a screen reader) can address
//! any control by name.
//!
//! Mechanism (mirrors `uq_workbench::headless_ui_tests`): each panel is mounted
//! once in a headless [`egui::Context`] with `enable_accesskit()`, its show
//! flag set, and the resulting accesskit node tree is swept. Every node whose
//! role is in the INTERACTIVE set must have a non-empty `name()` **or** a
//! `labelled_by()` edge that resolves to a node which itself has a name (the
//! `let lbl = ui.label(..); widget.labelled_by(lbl.id)` idiom egui recommends —
//! a bare `DragValue` / `from_id_source` `ComboBox` has an empty own-name, so it
//! must borrow one).
//!
//! egui 0.28 widget -> accesskit role map (see `egui/src/response.rs`):
//!   DragValue -> SpinButton · Button/ImageButton/CollapsingHeader -> Button ·
//!   ComboBox -> ComboBox · Checkbox -> CheckBox · Slider -> Slider ·
//!   TextEdit -> TextInput · SelectableLabel -> ToggleButton ·
//!   RadioButton -> RadioButton.

#![allow(clippy::field_reassign_with_default)]

use eframe::egui;
use egui::accesskit::{Node, NodeId, Role};

use crate::ValenxApp;

/// Roles that represent a user-operable control which an AI driver must be able
/// to address by name. Static-text, links, color wells and progress indicators
/// are intentionally excluded (not operated by name in the drive-by-name model;
/// a `ColorWell` carries no caption by construction).
const INTERACTIVE_ROLES: &[Role] = &[
    Role::SpinButton,   // DragValue
    Role::Button,       // Button / ImageButton / CollapsingHeader
    Role::ComboBox,     // ComboBox
    Role::CheckBox,     // Checkbox
    Role::Slider,       // Slider
    Role::TextInput,    // TextEdit
    Role::ToggleButton, // SelectableLabel (egui 0.28 maps these here)
    Role::RadioButton,  // RadioButton
];

/// A node has an accessible name if its own `name()` is non-empty, or it is
/// `labelled_by` at least one node that resolves (in `by_id`) to a non-empty
/// name. (egui clears a `DragValue`'s own name and gives a `from_id_source`
/// `ComboBox` an empty name, so those must borrow a caption's name.)
fn has_accessible_name(node: &Node, by_id: &std::collections::HashMap<NodeId, &Node>) -> bool {
    if node.name().is_some_and(|s| !s.trim().is_empty()) {
        return true;
    }
    node.labelled_by().iter().any(|lid| {
        by_id
            .get(lid)
            .and_then(|n| n.name())
            .is_some_and(|s| !s.trim().is_empty())
    })
}

/// Assert that EVERY interactive widget in `nodes` (an accesskit tree captured
/// from one headless frame of `panel`) carries an accessible name. On failure
/// the panic names the `panel`, the offending role, and any nearby text /
/// value the node does expose, so the unlabelled control can be located fast.
///
/// A widget that is `set_disabled()` (e.g. an `add_enabled_ui(false, ..)`
/// coefficient row) is still required to be named — egui keeps its caption /
/// `labelled_by` edge when greyed, and the codebase deliberately renders such
/// rows "always (greyed when unused) so the controls keep stable accessible
/// names" — so we do NOT exempt disabled nodes.
pub fn assert_all_interactive_widgets_named(panel: &str, nodes: &[(NodeId, Node)]) {
    let by_id: std::collections::HashMap<NodeId, &Node> =
        nodes.iter().map(|(id, n)| (*id, n)).collect();

    let mut offenders: Vec<String> = Vec::new();
    for (_, node) in nodes {
        let role = node.role();
        if !INTERACTIVE_ROLES.contains(&role) {
            continue;
        }
        if has_accessible_name(node, &by_id) {
            continue;
        }
        // Build a targeted hint from whatever the node DOES expose.
        let value_hint = node
            .value()
            .map(|v| format!(" value={v:?}"))
            .unwrap_or_default();
        let labelled = if node.labelled_by().is_empty() {
            " (no labelled_by edge)".to_string()
        } else {
            " (labelled_by edge present but its target has no name)".to_string()
        };
        offenders.push(format!("{role:?}{value_hint}{labelled}"));
    }

    assert!(
        offenders.is_empty(),
        "panel `{panel}` has {} unnamed interactive widget(s) — every control must \
         carry an accessible name (.labelled_by(caption.id) for DragValues / \
         from_id_source ComboBoxes / unlabelled Sliders, or label text on \
         buttons / checkboxes): [{}]",
        offenders.len(),
        offenders.join("; "),
    );
}

/// Draw one panel headlessly with accesskit enabled and return its node tree.
fn draw_and_collect(
    app: &mut ValenxApp,
    draw: fn(&mut ValenxApp, &egui::Context),
) -> Vec<(NodeId, Node)> {
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    let out = ctx.run(egui::RawInput::default(), |ctx| {
        draw(app, ctx);
    });
    out.platform_output
        .accesskit_update
        .expect("accesskit tree is produced when enabled")
        .nodes
}

/// Every `show_*_workbench` panel, paired with a closure that sets its show
/// flag and draws it. Mirrors the authoritative dispatch in `update.rs`
/// (one entry per `crate::*::draw_*_workbench(self, ctx)` call there, minus the
/// tab strip). Kept in lock-step with that block: if a panel is added there,
/// add it here too and this guard will start covering it.
#[allow(clippy::type_complexity)]
const PANELS: &[(&str, fn(&mut ValenxApp, &egui::Context))] = &[
    ("acidbase_workbench", |app, ctx| {
        app.show_acidbase_workbench = true;
        crate::acidbase_workbench::draw_acidbase_workbench(app, ctx);
    }),
    ("acoustics_workbench", |app, ctx| {
        app.show_acoustics_workbench = true;
        crate::acoustics_workbench::draw_acoustics_workbench(app, ctx);
    }),
    ("aero_workbench", |app, ctx| {
        app.show_aero_workbench = true;
        crate::aero_workbench::draw_aero_workbench(app, ctx);
    }),
    ("animate_workbench", |app, ctx| {
        app.show_animate_workbench = true;
        crate::animate_workbench::draw_animate_workbench(app, ctx);
    }),
    ("antenna_workbench", |app, ctx| {
        app.show_antenna_workbench = true;
        crate::antenna_workbench::draw_antenna_workbench(app, ctx);
    }),
    ("assistant_workbench", |app, ctx| {
        app.show_assistant_panel = true;
        crate::assistant_workbench::draw_assistant_workbench(app, ctx);
    }),
    ("astro_workbench", |app, ctx| {
        app.show_astro_workbench = true;
        crate::astro_workbench::draw_astro_workbench(app, ctx);
    }),
    ("autonomy_workbench", |app, ctx| {
        app.show_autonomy_workbench = true;
        crate::autonomy_workbench::draw_autonomy_workbench(app, ctx);
    }),
    ("batteryecm_workbench", |app, ctx| {
        app.show_batteryecm_workbench = true;
        crate::batteryecm_workbench::draw_batteryecm_workbench(app, ctx);
    }),
    ("batterypack_workbench", |app, ctx| {
        app.show_batterypack_workbench = true;
        crate::batterypack_workbench::draw_batterypack_workbench(app, ctx);
    }),
    ("beam_workbench", |app, ctx| {
        app.show_beam_workbench = true;
        crate::beam_workbench::draw_beam_workbench(app, ctx);
    }),
    ("bearing_workbench", |app, ctx| {
        app.show_bearing_workbench = true;
        crate::bearing_workbench::draw_bearing_workbench(app, ctx);
    }),
    ("beltdrive_workbench", |app, ctx| {
        app.show_beltdrive_workbench = true;
        crate::beltdrive_workbench::draw_beltdrive_workbench(app, ctx);
    }),
    ("bjt_workbench", |app, ctx| {
        app.show_bjt_workbench = true;
        crate::bjt_workbench::draw_bjt_workbench(app, ctx);
    }),
    ("blackhole_workbench", |app, ctx| {
        app.show_blackhole_workbench = true;
        crate::blackhole_workbench::draw_blackhole_workbench(app, ctx);
    }),
    ("bmr_workbench", |app, ctx| {
        app.show_bmr_workbench = true;
        crate::bmr_workbench::draw_bmr_workbench(app, ctx);
    }),
    ("bolt_workbench", |app, ctx| {
        app.show_bolt_workbench = true;
        crate::bolt_workbench::draw_bolt_workbench(app, ctx);
    }),
    ("bonemech_workbench", |app, ctx| {
        app.show_bonemech_workbench = true;
        crate::bonemech_workbench::draw_bonemech_workbench(app, ctx);
    }),
    ("brake_workbench", |app, ctx| {
        app.show_brake_workbench = true;
        crate::brake_workbench::draw_brake_workbench(app, ctx);
    }),
    ("buckling_workbench", |app, ctx| {
        app.show_buckling_workbench = true;
        crate::buckling_workbench::draw_buckling_workbench(app, ctx);
    }),
    ("cad_workbench", |app, ctx| {
        app.show_cad_workbench = true;
        crate::cad_workbench::draw_cad_workbench(app, ctx);
    }),
    ("camdynamics_workbench", |app, ctx| {
        app.show_camdynamics_workbench = true;
        crate::camdynamics_workbench::draw_camdynamics_workbench(app, ctx);
    }),
    ("capacitor_workbench", |app, ctx| {
        app.show_capacitor_workbench = true;
        crate::capacitor_workbench::draw_capacitor_workbench(app, ctx);
    }),
    ("car_workbench", |app, ctx| {
        app.show_car_workbench = true;
        crate::car_workbench::draw_car_workbench(app, ctx);
    }),
    ("cfd_workbench", |app, ctx| {
        app.show_cfd_workbench = true;
        crate::cfd_workbench::draw_cfd_workbench(app, ctx);
    }),
    ("chaindrive_workbench", |app, ctx| {
        app.show_chaindrive_workbench = true;
        crate::chaindrive_workbench::draw_chaindrive_workbench(app, ctx);
    }),
    ("clutch_workbench", |app, ctx| {
        app.show_clutch_workbench = true;
        crate::clutch_workbench::draw_clutch_workbench(app, ctx);
    }),
    ("coil_workbench", |app, ctx| {
        app.show_coil_workbench = true;
        crate::coil_workbench::draw_coil_workbench(app, ctx);
    }),
    ("collision_workbench", |app, ctx| {
        app.show_collision_workbench = true;
        crate::collision_workbench::draw_collision_workbench(app, ctx);
    }),
    ("columnsteel_workbench", |app, ctx| {
        app.show_columnsteel_workbench = true;
        crate::columnsteel_workbench::draw_columnsteel_workbench(app, ctx);
    }),
    ("combustion_workbench", |app, ctx| {
        app.show_combustion_workbench = true;
        crate::combustion_workbench::draw_combustion_workbench(app, ctx);
    }),
    ("conveyor_workbench", |app, ctx| {
        app.show_conveyor_workbench = true;
        crate::conveyor_workbench::draw_conveyor_workbench(app, ctx);
    }),
    ("cosim_workbench", |app, ctx| {
        app.show_cosim_workbench = true;
        crate::cosim_workbench::draw_cosim_workbench(app, ctx);
    }),
    ("creep_workbench", |app, ctx| {
        app.show_creep_workbench = true;
        crate::creep_workbench::draw_creep_workbench(app, ctx);
    }),
    ("dcmotor_workbench", |app, ctx| {
        app.show_dcmotor_workbench = true;
        crate::dcmotor_workbench::draw_dcmotor_workbench(app, ctx);
    }),
    ("diffusion_workbench", |app, ctx| {
        app.show_diffusion_workbench = true;
        crate::diffusion_workbench::draw_diffusion_workbench(app, ctx);
    }),
    ("dimensional_workbench", |app, ctx| {
        app.show_dimensional_workbench = true;
        crate::dimensional_workbench::draw_dimensional_workbench(app, ctx);
    }),
    ("draft2d_workbench", |app, ctx| {
        app.show_draft2d_workbench = true;
        crate::draft2d_workbench::draw_draft2d_workbench(app, ctx);
    }),
    ("drone_workbench", |app, ctx| {
        app.show_drone_workbench = true;
        crate::drone_workbench::draw_drone_workbench(app, ctx);
    }),
    ("electrochem_workbench", |app, ctx| {
        app.show_electrochem_workbench = true;
        crate::electrochem_workbench::draw_electrochem_workbench(app, ctx);
    }),
    ("engine_workbench", |app, ctx| {
        app.show_engine_workbench = true;
        crate::engine_workbench::draw_engine_workbench(app, ctx);
    }),
    ("enzymekinetics_workbench", |app, ctx| {
        app.show_enzymekinetics_workbench = true;
        crate::enzymekinetics_workbench::draw_enzymekinetics_workbench(app, ctx);
    }),
    ("fanlaws_workbench", |app, ctx| {
        app.show_fanlaws_workbench = true;
        crate::fanlaws_workbench::draw_fanlaws_workbench(app, ctx);
    }),
    ("fasteners_workbench", |app, ctx| {
        app.show_fasteners_workbench = true;
        crate::fasteners_workbench::draw_fasteners_workbench(app, ctx);
    }),
    ("fatigue_workbench", |app, ctx| {
        app.show_fatigue_workbench = true;
        crate::fatigue_workbench::draw_fatigue_workbench(app, ctx);
    }),
    ("fem_workbench", |app, ctx| {
        app.show_fem_workbench = true;
        crate::fem_workbench::draw_fem_workbench(app, ctx);
    }),
    ("fft_workbench", |app, ctx| {
        app.show_fft_workbench = true;
        crate::fft_workbench::draw_fft_workbench(app, ctx);
    }),
    ("fields_workbench", |app, ctx| {
        app.show_fields_workbench = true;
        crate::fields_workbench::draw_fields_workbench(app, ctx);
    }),
    ("filter_workbench", |app, ctx| {
        app.show_filter_workbench = true;
        crate::filter_workbench::draw_filter_workbench(app, ctx);
    }),
    ("fixedwing_workbench", |app, ctx| {
        app.show_fixedwing_workbench = true;
        crate::fixedwing_workbench::draw_fixedwing_workbench(app, ctx);
    }),
    ("fluids_workbench", |app, ctx| {
        app.show_fluids_workbench = true;
        crate::fluids_workbench::draw_fluids_workbench(app, ctx);
    }),
    ("fluidstatics_workbench", |app, ctx| {
        app.show_fluidstatics_workbench = true;
        crate::fluidstatics_workbench::draw_fluidstatics_workbench(app, ctx);
    }),
    ("flywheel_workbench", |app, ctx| {
        app.show_flywheel_workbench = true;
        crate::flywheel_workbench::draw_flywheel_workbench(app, ctx);
    }),
    ("fourbar_workbench", |app, ctx| {
        app.show_fourbar_workbench = true;
        crate::fourbar_workbench::draw_fourbar_workbench(app, ctx);
    }),
    ("fracture_workbench", |app, ctx| {
        app.show_fracture_workbench = true;
        crate::fracture_workbench::draw_fracture_workbench(app, ctx);
    }),
    ("frames_workbench", |app, ctx| {
        app.show_frames_workbench = true;
        crate::frames_workbench::draw_frames_workbench(app, ctx);
    }),
    ("gasdynamics_workbench", |app, ctx| {
        app.show_gasdynamics_workbench = true;
        crate::gasdynamics_workbench::draw_gasdynamics_workbench(app, ctx);
    }),
    ("gearbox_workbench", |app, ctx| {
        app.show_gearbox_workbench = true;
        crate::gearbox_workbench::draw_gearbox_workbench(app, ctx);
    }),
    ("gears_workbench", |app, ctx| {
        app.show_gears_workbench = true;
        crate::gears_workbench::draw_gears_workbench(app, ctx);
    }),
    ("geartooth_workbench", |app, ctx| {
        app.show_geartooth_workbench = true;
        crate::geartooth_workbench::draw_geartooth_workbench(app, ctx);
    }),
    ("genetics_workbench", |app, ctx| {
        app.show_genetics_workbench = true;
        crate::genetics_workbench::draw_genetics_workbench(app, ctx);
    }),
    ("geomatics_workbench", |app, ctx| {
        app.show_geomatics_workbench = true;
        crate::geomatics_workbench::draw_geomatics_workbench(app, ctx);
    }),
    ("heatexchanger_workbench", |app, ctx| {
        app.show_heatexchanger_workbench = true;
        crate::heatexchanger_workbench::draw_heatexchanger_workbench(app, ctx);
    }),
    ("heatpump_workbench", |app, ctx| {
        app.show_heatpump_workbench = true;
        crate::heatpump_workbench::draw_heatpump_workbench(app, ctx);
    }),
    ("heattransfer_workbench", |app, ctx| {
        app.show_heattransfer_workbench = true;
        crate::heattransfer_workbench::draw_heattransfer_workbench(app, ctx);
    }),
    ("hemodynamics_workbench", |app, ctx| {
        app.show_hemodynamics_workbench = true;
        crate::hemodynamics_workbench::draw_hemodynamics_workbench(app, ctx);
    }),
    ("hvac_workbench", |app, ctx| {
        app.show_hvac_workbench = true;
        crate::hvac_workbench::draw_hvac_workbench(app, ctx);
    }),
    ("hydraulics_workbench", |app, ctx| {
        app.show_hydraulics_workbench = true;
        crate::hydraulics_workbench::draw_hydraulics_workbench(app, ctx);
    }),
    ("inclinedplane_workbench", |app, ctx| {
        app.show_inclinedplane_workbench = true;
        crate::inclinedplane_workbench::draw_inclinedplane_workbench(app, ctx);
    }),
    ("inductionmotor_workbench", |app, ctx| {
        app.show_inductionmotor_workbench = true;
        crate::inductionmotor_workbench::draw_inductionmotor_workbench(app, ctx);
    }),
    ("insulation_workbench", |app, ctx| {
        app.show_insulation_workbench = true;
        crate::insulation_workbench::draw_insulation_workbench(app, ctx);
    }),
    ("interior_workbench", |app, ctx| {
        app.show_interior_workbench = true;
        crate::interior_workbench::draw_interior_workbench(app, ctx);
    }),
    ("leadscrew_workbench", |app, ctx| {
        app.show_leadscrew_workbench = true;
        crate::leadscrew_workbench::draw_leadscrew_workbench(app, ctx);
    }),
    ("led_workbench", |app, ctx| {
        app.show_led_workbench = true;
        crate::led_workbench::draw_led_workbench(app, ctx);
    }),
    ("leverage_workbench", |app, ctx| {
        app.show_leverage_workbench = true;
        crate::leverage_workbench::draw_leverage_workbench(app, ctx);
    }),
    ("marine_workbench", |app, ctx| {
        app.show_marine_workbench = true;
        crate::marine_workbench::draw_marine_workbench(app, ctx);
    }),
    ("mbd_workbench", |app, ctx| {
        app.show_mbd_workbench = true;
        crate::mbd_workbench::draw_mbd_workbench(app, ctx);
    }),
    ("mesh_toolbox", |app, ctx| {
        app.show_mesh_toolbox = true;
        crate::mesh_toolbox::draw_mesh_toolbox(app, ctx);
    }),
    ("missionsim_workbench", |app, ctx| {
        app.show_missionsim_workbench = true;
        crate::missionsim_workbench::draw_missionsim_workbench(app, ctx);
    }),
    ("mohr_workbench", |app, ctx| {
        app.show_mohr_workbench = true;
        crate::mohr_workbench::draw_mohr_workbench(app, ctx);
    }),
    ("mosfet_workbench", |app, ctx| {
        app.show_mosfet_workbench = true;
        crate::mosfet_workbench::draw_mosfet_workbench(app, ctx);
    }),
    ("neuro_workbench", |app, ctx| {
        app.show_neuro_workbench = true;
        crate::neuro_workbench::draw_neuro_workbench(app, ctx);
    }),
    ("ocean_workbench", |app, ctx| {
        app.show_ocean_workbench = true;
        crate::ocean_workbench::draw_ocean_workbench(app, ctx);
    }),
    ("opamp_workbench", |app, ctx| {
        app.show_opamp_workbench = true;
        crate::opamp_workbench::draw_opamp_workbench(app, ctx);
    }),
    ("openchannel_workbench", |app, ctx| {
        app.show_openchannel_workbench = true;
        crate::openchannel_workbench::draw_openchannel_workbench(app, ctx);
    }),
    ("optics_workbench", |app, ctx| {
        app.show_optics_workbench = true;
        crate::optics_workbench::draw_optics_workbench(app, ctx);
    }),
    ("orifice_workbench", |app, ctx| {
        app.show_orifice_workbench = true;
        crate::orifice_workbench::draw_orifice_workbench(app, ctx);
    }),
    ("osmosis_workbench", |app, ctx| {
        app.show_osmosis_workbench = true;
        crate::osmosis_workbench::draw_osmosis_workbench(app, ctx);
    }),
    ("param_sketch_panel", |app, ctx| {
        app.show_param_sketch = true;
        crate::param_sketch_panel::draw_param_sketch_workbench(app, ctx);
    }),
    ("pharmacokinetics_workbench", |app, ctx| {
        app.show_pharmacokinetics_workbench = true;
        crate::pharmacokinetics_workbench::draw_pharmacokinetics_workbench(app, ctx);
    }),
    ("photogrammetry_workbench", |app, ctx| {
        app.show_photogrammetry_workbench = true;
        crate::photogrammetry_workbench::draw_photogrammetry_workbench(app, ctx);
    }),
    ("pipeflow_workbench", |app, ctx| {
        app.show_pipeflow_workbench = true;
        crate::pipeflow_workbench::draw_pipeflow_workbench(app, ctx);
    }),
    ("pipenetwork_workbench", |app, ctx| {
        app.show_pipenetwork_workbench = true;
        crate::pipenetwork_workbench::draw_pipenetwork_workbench(app, ctx);
    }),
    ("piping_workbench", |app, ctx| {
        app.show_piping_workbench = true;
        crate::piping_workbench::draw_piping_workbench(app, ctx);
    }),
    ("plate_workbench", |app, ctx| {
        app.show_plate_workbench = true;
        crate::plate_workbench::draw_plate_workbench(app, ctx);
    }),
    ("pneumatics_workbench", |app, ctx| {
        app.show_pneumatics_workbench = true;
        crate::pneumatics_workbench::draw_pneumatics_workbench(app, ctx);
    }),
    ("popdynamics_workbench", |app, ctx| {
        app.show_popdynamics_workbench = true;
        crate::popdynamics_workbench::draw_popdynamics_workbench(app, ctx);
    }),
    ("powerfactor_workbench", |app, ctx| {
        app.show_powerfactor_workbench = true;
        crate::powerfactor_workbench::draw_powerfactor_workbench(app, ctx);
    }),
    ("ppi_workbench", |app, ctx| {
        app.show_ppi_workbench = true;
        crate::ppi_workbench::draw_ppi_workbench(app, ctx);
    }),
    ("pressurevessel_workbench", |app, ctx| {
        app.show_pressurevessel_workbench = true;
        crate::pressurevessel_workbench::draw_pressurevessel_workbench(app, ctx);
    }),
    ("projectile_workbench", |app, ctx| {
        app.show_projectile_workbench = true;
        crate::projectile_workbench::draw_projectile_workbench(app, ctx);
    }),
    ("psychrometrics_workbench", |app, ctx| {
        app.show_psychrometrics_workbench = true;
        crate::psychrometrics_workbench::draw_psychrometrics_workbench(app, ctx);
    }),
    ("pulley_workbench", |app, ctx| {
        app.show_pulley_workbench = true;
        crate::pulley_workbench::draw_pulley_workbench(app, ctx);
    }),
    ("pump_workbench", |app, ctx| {
        app.show_pump_workbench = true;
        crate::pump_workbench::draw_pump_workbench(app, ctx);
    }),
    ("queueing_workbench", |app, ctx| {
        app.show_queueing_workbench = true;
        crate::queueing_workbench::draw_queueing_workbench(app, ctx);
    }),
    ("radioactivity_workbench", |app, ctx| {
        app.show_radioactivity_workbench = true;
        crate::radioactivity_workbench::draw_radioactivity_workbench(app, ctx);
    }),
    ("rail_workbench", |app, ctx| {
        app.show_rail_workbench = true;
        crate::rail_workbench::draw_rail_workbench(app, ctx);
    }),
    ("rcbeam_workbench", |app, ctx| {
        app.show_rcbeam_workbench = true;
        crate::rcbeam_workbench::draw_rcbeam_workbench(app, ctx);
    }),
    ("reactdyn_workbench", |app, ctx| {
        app.show_reactdyn_workbench = true;
        crate::reactdyn_workbench::draw_reactdyn_workbench(app, ctx);
    }),
    ("rectifier_workbench", |app, ctx| {
        app.show_rectifier_workbench = true;
        crate::rectifier_workbench::draw_rectifier_workbench(app, ctx);
    }),
    ("refrigeration_workbench", |app, ctx| {
        app.show_refrigeration_workbench = true;
        crate::refrigeration_workbench::draw_refrigeration_workbench(app, ctx);
    }),
    ("reinforcement_workbench", |app, ctx| {
        app.show_reinforcement_workbench = true;
        crate::reinforcement_workbench::draw_reinforcement_workbench(app, ctx);
    }),
    ("render_workbench", |app, ctx| {
        app.show_render_workbench = true;
        crate::render_workbench::draw_render_workbench(app, ctx);
    }),
    ("resistornetwork_workbench", |app, ctx| {
        app.show_resistornetwork_workbench = true;
        crate::resistornetwork_workbench::draw_resistornetwork_workbench(app, ctx);
    }),
    ("retainingwall_workbench", |app, ctx| {
        app.show_retainingwall_workbench = true;
        crate::retainingwall_workbench::draw_retainingwall_workbench(app, ctx);
    }),
    ("reverse_workbench", |app, ctx| {
        app.show_reverse_workbench = true;
        crate::reverse_workbench::draw_reverse_workbench(app, ctx);
    }),
    ("rivet_workbench", |app, ctx| {
        app.show_rivet_workbench = true;
        crate::rivet_workbench::draw_rivet_workbench(app, ctx);
    }),
    ("rocket_workbench", |app, ctx| {
        app.show_rocket_workbench = true;
        crate::rocket_workbench::draw_rocket_workbench(app, ctx);
    }),
    ("rom_workbench", |app, ctx| {
        app.show_rom_workbench = true;
        crate::rom_workbench::draw_rom_workbench(app, ctx);
    }),
    ("rotor_workbench", |app, ctx| {
        app.show_rotor_workbench = true;
        crate::rotor_workbench::draw_rotor_workbench(app, ctx);
    }),
    ("screwthread_workbench", |app, ctx| {
        app.show_screwthread_workbench = true;
        crate::screwthread_workbench::draw_screwthread_workbench(app, ctx);
    }),
    ("shaftdesign_workbench", |app, ctx| {
        app.show_shaftdesign_workbench = true;
        crate::shaftdesign_workbench::draw_shaftdesign_workbench(app, ctx);
    }),
    ("sheetmetal_workbench", |app, ctx| {
        app.show_sheetmetal_workbench = true;
        crate::sheetmetal_workbench::draw_sheetmetal_workbench(app, ctx);
    }),
    ("soilbearing_workbench", |app, ctx| {
        app.show_soilbearing_workbench = true;
        crate::soilbearing_workbench::draw_soilbearing_workbench(app, ctx);
    }),
    ("solarpv_workbench", |app, ctx| {
        app.show_solarpv_workbench = true;
        crate::solarpv_workbench::draw_solarpv_workbench(app, ctx);
    }),
    ("springcombination_workbench", |app, ctx| {
        app.show_springcombination_workbench = true;
        crate::springcombination_workbench::draw_springcombination_workbench(app, ctx);
    }),
    ("springdesign_workbench", |app, ctx| {
        app.show_springdesign_workbench = true;
        crate::springdesign_workbench::draw_springdesign_workbench(app, ctx);
    }),
    ("springs_workbench", |app, ctx| {
        app.show_springs_workbench = true;
        crate::springs_workbench::draw_springs_workbench(app, ctx);
    }),
    ("statics_workbench", |app, ctx| {
        app.show_statics_workbench = true;
        crate::statics_workbench::draw_statics_workbench(app, ctx);
    }),
    ("straingauge_workbench", |app, ctx| {
        app.show_straingauge_workbench = true;
        crate::straingauge_workbench::draw_straingauge_workbench(app, ctx);
    }),
    ("strainrosette_workbench", |app, ctx| {
        app.show_strainrosette_workbench = true;
        crate::strainrosette_workbench::draw_strainrosette_workbench(app, ctx);
    }),
    ("survivability_workbench", |app, ctx| {
        app.show_survivability_workbench = true;
        crate::survivability_workbench::draw_survivability_workbench(app, ctx);
    }),
    ("thermalexpansion_workbench", |app, ctx| {
        app.show_thermalexpansion_workbench = true;
        crate::thermalexpansion_workbench::draw_thermalexpansion_workbench(app, ctx);
    }),
    ("thermistor_workbench", |app, ctx| {
        app.show_thermistor_workbench = true;
        crate::thermistor_workbench::draw_thermistor_workbench(app, ctx);
    }),
    ("thermocouple_workbench", |app, ctx| {
        app.show_thermocouple_workbench = true;
        crate::thermocouple_workbench::draw_thermocouple_workbench(app, ctx);
    }),
    ("thermocycle_workbench", |app, ctx| {
        app.show_thermocycle_workbench = true;
        crate::thermocycle_workbench::draw_thermocycle_workbench(app, ctx);
    }),
    ("thermoreg_workbench", |app, ctx| {
        app.show_thermoreg_workbench = true;
        crate::thermoreg_workbench::draw_thermoreg_workbench(app, ctx);
    }),
    ("threephase_workbench", |app, ctx| {
        app.show_threephase_workbench = true;
        crate::threephase_workbench::draw_threephase_workbench(app, ctx);
    }),
    ("torsion_workbench", |app, ctx| {
        app.show_torsion_workbench = true;
        crate::torsion_workbench::draw_torsion_workbench(app, ctx);
    }),
    ("transformer_workbench", |app, ctx| {
        app.show_transformer_workbench = true;
        crate::transformer_workbench::draw_transformer_workbench(app, ctx);
    }),
    ("transmissionline_workbench", |app, ctx| {
        app.show_transmissionline_workbench = true;
        crate::transmissionline_workbench::draw_transmissionline_workbench(app, ctx);
    }),
    ("truss_workbench", |app, ctx| {
        app.show_truss_workbench = true;
        crate::truss_workbench::draw_truss_workbench(app, ctx);
    }),
    ("uas_workbench", |app, ctx| {
        app.show_uas_workbench = true;
        crate::uas_workbench::draw_uas_workbench(app, ctx);
    }),
    ("uq_workbench", |app, ctx| {
        app.show_uq_workbench = true;
        crate::uq_workbench::draw_uq_workbench(app, ctx);
    }),
    ("variant_effect_workbench", |app, ctx| {
        app.show_variant_effect_workbench = true;
        crate::variant_effect_workbench::draw_variant_effect_workbench(app, ctx);
    }),
    ("vibration_workbench", |app, ctx| {
        app.show_vibration_workbench = true;
        crate::vibration_workbench::draw_vibration_workbench(app, ctx);
    }),
    ("weir_workbench", |app, ctx| {
        app.show_weir_workbench = true;
        crate::weir_workbench::draw_weir_workbench(app, ctx);
    }),
    ("windturbine_workbench", |app, ctx| {
        app.show_windturbine_workbench = true;
        crate::windturbine_workbench::draw_windturbine_workbench(app, ctx);
    }),
];

#[test]
fn every_workbench_panel_names_all_its_interactive_widgets() {
    // Table-driven: draw each panel once headlessly and assert that no
    // interactive widget is left without an accessible name. A regression
    // (a new bare DragValue / from_id_source ComboBox / unlabelled Slider on
    // any panel) fails here with the panel + role called out.
    for (panel, draw) in PANELS {
        let mut app = ValenxApp::default();
        let nodes = draw_and_collect(&mut app, *draw);
        assert_all_interactive_widgets_named(panel, &nodes);
    }
}

#[test]
fn panel_table_covers_every_workbench_in_the_dispatch() {
    // Guard the guard: the table must stay as large as the real dispatch so a
    // newly-added panel can't silently escape the naming check. The count is
    // the number of `crate::*::draw_*(self, ctx)` workbench calls in
    // `update.rs` (the tab strip is not a workbench and is excluded).
    assert_eq!(
        PANELS.len(),
        148,
        "PANELS must list every show_*_workbench panel drawn in update.rs"
    );
    // No duplicate panel ids.
    let mut seen = std::collections::HashSet::new();
    for (p, _) in PANELS {
        assert!(seen.insert(*p), "duplicate panel id in PANELS: {p}");
    }
}

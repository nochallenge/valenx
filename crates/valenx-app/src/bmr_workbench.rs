//! The right-side **BMR / TDEE Workbench** panel — native
//! basal-metabolic-rate and total-daily-energy-expenditure estimation
//! over `valenx-bmr`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_bmr_workbench`,
//! toggled from the View menu. The form takes a person's biological sex,
//! age, height and mass and a physical-activity level; "Analyze" evaluates
//! the chosen BMR regression ([`valenx_bmr::BmrEquation`] —
//! Mifflin-St Jeor or Harris-Benedict) and scales it to a daily
//! expenditure ([`valenx_bmr::tdee_for_level`]), reporting BMR, the
//! activity multiplier and TDEE in kcal/day, and "Show 3-D figure" loads a
//! representative stylized human figure into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_bmr::{ActivityLevel, BmrEquation, Sex};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the BMR / TDEE Workbench.
pub struct BmrWorkbenchState {
    /// Biological sex used by the regression.
    sex: Sex,
    /// Age (years).
    age_years: f64,
    /// Standing height (cm).
    height_cm: f64,
    /// Body mass (kg).
    mass_kg: f64,
    /// Which predictive BMR equation to evaluate.
    equation: BmrEquation,
    /// Physical-activity level whose factor scales BMR into TDEE.
    activity: ActivityLevel,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D figure (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for BmrWorkbenchState {
    fn default() -> Self {
        // A 30-year-old, 180 cm, 80 kg male at moderate activity:
        // Mifflin-St Jeor BMR = 1780 kcal/day, TDEE = 1780 * 1.55 = 2759.
        Self {
            sex: Sex::Male,
            age_years: 30.0,
            height_cm: 180.0,
            mass_kg: 80.0,
            equation: BmrEquation::MifflinStJeor,
            activity: ActivityLevel::ModeratelyActive,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the BMR / TDEE Workbench right-side panel. A no-op when the
/// `show_bmr_workbench` toggle is off.
pub fn draw_bmr_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bmr_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_bmr_workbench",
        "BMR / TDEE",
        |app, ui| {
            ui.label(
                egui::RichText::new("native resting + daily energy expenditure · valenx-bmr")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.bmr;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Sex").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.sex, Sex::Male, "male");
                        ui.radio_value(&mut s.sex, Sex::Female, "female");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Anthropometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("age (yr)");
                        ui.add(egui::DragValue::new(&mut s.age_years).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("height (cm)");
                        ui.add(egui::DragValue::new(&mut s.height_cm).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Equation").strong());
                    ui.radio_value(
                        &mut s.equation,
                        BmrEquation::MifflinStJeor,
                        "Mifflin-St Jeor (1990)",
                    );
                    ui.radio_value(
                        &mut s.equation,
                        BmrEquation::HarrisBenedict,
                        "Harris-Benedict (Roza & Shizgal 1984)",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Activity level").strong());
                    ui.radio_value(&mut s.activity, ActivityLevel::Sedentary, "sedentary ×1.2");
                    ui.radio_value(
                        &mut s.activity,
                        ActivityLevel::LightlyActive,
                        "lightly active ×1.375",
                    );
                    ui.radio_value(
                        &mut s.activity,
                        ActivityLevel::ModeratelyActive,
                        "moderately active ×1.55",
                    );
                    ui.radio_value(
                        &mut s.activity,
                        ActivityLevel::VeryActive,
                        "very active ×1.725",
                    );
                    ui.radio_value(
                        &mut s.activity,
                        ActivityLevel::ExtraActive,
                        "extra active ×1.9",
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_bmr(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D figure").strong())
                        .on_hover_text(
                            "Build a representative stylized human figure (head, torso and four limbs) as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Energy expenditure").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_bmr_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.bmr` borrow is released
    // here): build the figure's 3-D solid and load it.
    if app.bmr.show_3d_request {
        app.bmr.show_3d_request = false;
        load_figure_3d(app);
    }
}

/// Validate the form, evaluate the BMR/TDEE and format the readout.
fn run_bmr(s: &mut BmrWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the chosen BMR equation and scale it to a daily expenditure,
/// formatting the full readout and mapping any domain error to a display
/// string. Extracted so it is unit-testable.
fn compute(s: &BmrWorkbenchState) -> Result<String, String> {
    let bmr = s
        .equation
        .evaluate(s.sex, s.mass_kg, s.height_cm, s.age_years)
        .map_err(|e| e.to_string())?;
    let factor = s.activity.factor();
    let tdee = valenx_bmr::tdee_for_level(bmr, s.activity).map_err(|e| e.to_string())?;

    let sex = match s.sex {
        Sex::Male => "male",
        Sex::Female => "female",
    };
    let equation = match s.equation {
        BmrEquation::MifflinStJeor => "Mifflin-St Jeor",
        BmrEquation::HarrisBenedict => "Harris-Benedict",
    };

    Ok(format!(
        "sex             : {sex}\n\
         age             : {:.0} yr\n\
         height / mass   : {:.1} cm / {:.1} kg\n\
         equation        : {equation}\n\n\
         BMR             : {bmr:.1} kcal/day\n\
         activity factor : ×{factor:.3}\n\
         TDEE            : {tdee:.1} kcal/day",
        s.age_years, s.height_cm, s.mass_kg,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Build a representative stylized human figure as a triangle [`Mesh`] —
/// a torso box, a head box atop it, two arm boxes at the shoulders and two
/// leg boxes below. Representative geometry only (not to scale, not
/// person-specific); the BMR / TDEE numbers are the `valenx-bmr` result.
/// `None` for an invalid configuration (the form fails to evaluate a BMR).
fn figure_solid_mesh(s: &BmrWorkbenchState) -> Option<Mesh> {
    // Gate on the crate actually accepting the anthropometry.
    s.equation
        .evaluate(s.sex, s.mass_kg, s.height_cm, s.age_years)
        .ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Torso.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.1),
        Vector3::new(0.18, 0.11, 0.32),
    );
    // Head.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.56),
        Vector3::new(0.11, 0.1, 0.13),
    );
    // Left arm.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.26, 0.0, 1.05),
        Vector3::new(0.06, 0.06, 0.3),
    );
    // Right arm.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.26, 0.0, 1.05),
        Vector3::new(0.06, 0.06, 0.3),
    );
    // Left leg.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.09, 0.0, 0.4),
        Vector3::new(0.07, 0.07, 0.4),
    );
    // Right leg.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.09, 0.0, 0.4),
        Vector3::new(0.07, 0.07, 0.4),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-bmr");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D figure solid and load it into the central viewport.
fn load_figure_3d(app: &mut ValenxApp) {
    let Some(mesh) = figure_solid_mesh(&app.bmr) else {
        app.bmr.error = Some("inputs are invalid — cannot build the 3-D figure".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<figure>/valenx-bmr"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"bmr"}`** product: a DATA-ONLY text card
/// of the workbench's own BMR / TDEE headline numbers. The metabolic result
/// has no characteristic shape — the panel's stylised human figure is a
/// representative decorative solid, not a real object — so the bridge product
/// is right-sized to a card (`mesh: None`) carrying just the readout (the
/// confidence badge is appended centrally). The panel's "Show 3-D figure"
/// button still builds that representative solid into the central viewport.
/// Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`BmrWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` readout.
pub(crate) fn bmr_product() -> crate::WorkspaceProduct {
    let s = BmrWorkbenchState::default();
    let readout = compute(&s).expect("default body inputs ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    crate::WorkspaceProduct {
        title: "BMR / TDEE".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = BmrWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_bmr_and_tdee() {
        let mut s = BmrWorkbenchState::default();
        run_bmr(&mut s);
        assert!(s.error.is_none(), "default should analyze: {:?}", s.error);
        assert!(s.result.contains("BMR"));
        assert!(s.result.contains("TDEE"));
        // 30 y, 180 cm, 80 kg male, Mifflin-St Jeor = 1780 kcal/day,
        // moderate activity TDEE = 1780 * 1.55 = 2759 kcal/day.
        assert!(s.result.contains("1780.0"));
        assert!(s.result.contains("2759.0"));
    }

    #[test]
    fn analyze_rejects_zero_mass() {
        let mut s = BmrWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        run_bmr(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn mifflin_male_ground_truth_and_tdee_scaling() {
        // Ground truth: Mifflin-St Jeor male BMR is exactly
        // 10*mass + 6.25*height - 5*age + 5, and TDEE is that BMR times
        // the activity factor. Hand-computed for the default body.
        let mass = 80.0;
        let height = 180.0;
        let age = 30.0;
        let bmr: f64 = 10.0 * mass + 6.25 * height - 5.0 * age + 5.0;
        assert!((bmr - 1780.0).abs() < 1e-9, "hand BMR {bmr}");
        let factor = ActivityLevel::ModeratelyActive.factor();
        assert!((factor - 1.55).abs() < 1e-12);
        let tdee = bmr * factor;
        assert!((tdee - 2759.0).abs() < 1e-9, "hand TDEE {tdee}");
        // The crate path must agree with the hand computation.
        let via_crate = valenx_bmr::mifflin_st_jeor(Sex::Male, mass, height, age).unwrap();
        assert!((via_crate - bmr).abs() < 1e-9);
    }

    #[test]
    fn female_bmr_below_male_for_same_body() {
        // The Mifflin-St Jeor female intercept is 166 kcal lower, so a
        // female TDEE must come out below the male one for the same body.
        let mut male = BmrWorkbenchState::default();
        run_bmr(&mut male);
        let mut female = BmrWorkbenchState {
            sex: Sex::Female,
            ..Default::default()
        };
        run_bmr(&mut female);
        assert!(male.error.is_none() && female.error.is_none());
        // 1780 male vs 1614 female BMR -> the female readout must differ.
        assert!(female.result.contains("1614.0"));
        assert!(!female.result.contains("1780.0"));
    }

    #[test]
    fn figure_mesh_for_default_is_nonempty_and_in_range() {
        let s = BmrWorkbenchState::default();
        let mesh = figure_solid_mesh(&s).expect("default figure yields a solid");
        assert!(mesh.nodes.len() > 8, "expected torso + head + four limbs");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn figure_mesh_none_for_invalid() {
        let s = BmrWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        assert!(figure_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bmr_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_bmr_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bmr_workbench = true;
        run_bmr(&mut app.bmr);
        draw_workbench(&mut app);
    }
}

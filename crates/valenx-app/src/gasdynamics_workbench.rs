//! The right-side **Gas Dynamics** workbench panel — a 1-D
//! compressible-flow calculator over `valenx-gasdynamics`.
//!
//! Mirrors the frames / fasteners / … workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_gasdynamics_workbench`,
//! toggled from the View menu and off by default. The form takes a Mach
//! number `M` and specific-heat ratio `gamma`; the panel reports — live, in
//! process — the full classic compressible-flow set for a calorically-perfect
//! gas:
//!
//! - **Isentropic** stagnation ratios `T0/T`, `p0/p`, `rho0/rho` and the
//!   area-Mach ratio `A/A*`.
//! - **Normal shock** (supersonic only): downstream `M2` and the `p2/p1`,
//!   `T2/T1`, `rho2/rho1`, `p02/p01` jumps.
//! - **Prandtl-Meyer** (supersonic only): the expansion angle `nu(M)` and the
//!   Mach angle `mu = asin(1/M)`.
//! - **Fanno** flow (adiabatic + friction): the sonic-referenced `T/T*`,
//!   `p/p*`, `rho/rho*`, `V/V*`, `p0/p0*` and the friction parameter `4fL*/D`.
//! - **Rayleigh** flow (heat addition): the sonic-referenced `T/T*`, `p/p*`,
//!   `rho/rho*`, `V/V*`, `T0/T0*` and `p0/p0*`.
//!
//! Research / educational grade — the same calorically-perfect ideal-gas scope
//! and caveats as `valenx-gasdynamics` itself.

use eframe::egui;

use valenx_gasdynamics::{
    area_mach_ratio, fanno_state, normal_shock, prandtl_meyer_angle, rayleigh_state,
    stagnation_ratios, FannoState, GasError, NormalShock, RayleighState, StagnationRatios,
};

use crate::ValenxApp;

/// Persistent form + result state for the Gas Dynamics workbench.
pub struct GasDynamicsWorkbenchState {
    /// Free-stream / local Mach number `M` (`> 0`).
    mach: f64,
    /// Specific-heat ratio `gamma = cp / cv` (`> 1`; air `1.4`).
    gamma: f64,
    /// Formatted multi-section readout (empty until the first compute).
    result: String,
    /// Validation / domain error, if any.
    error: Option<String>,
}

impl Default for GasDynamicsWorkbenchState {
    fn default() -> Self {
        // A representative supersonic air state (the NACA-1135 worked example).
        Self {
            mach: 2.0,
            gamma: 1.4,
            result: String::new(),
            error: None,
        }
    }
}

/// The full compressible-flow report at one `(M, gamma)`. The shock /
/// expansion sections are `None` for subsonic flow, where they are undefined.
struct FlowReport {
    iso: StagnationRatios,
    area_mach: f64,
    fanno: FannoState,
    rayleigh: RayleighState,
    /// Normal-shock jump — `Some` only for `M >= 1`.
    shock: Option<NormalShock>,
    /// Prandtl-Meyer angle `nu(M)` in degrees — `Some` only for `M >= 1`.
    nu_deg: Option<f64>,
    /// Mach angle `mu = asin(1/M)` in degrees — `Some` only for `M >= 1`.
    mu_deg: Option<f64>,
}

/// Evaluate every defined relation at `(mach, gamma)`. The subsonic-defined
/// quantities (isentropic, area-Mach, Fanno, Rayleigh) are always computed; the
/// normal shock and Prandtl-Meyer expansion are only added for `M >= 1`.
fn compute_report(mach: f64, gamma: f64) -> Result<FlowReport, GasError> {
    let iso = stagnation_ratios(mach, gamma)?;
    let area_mach = area_mach_ratio(mach, gamma)?;
    let fanno = fanno_state(mach, gamma)?;
    let rayleigh = rayleigh_state(mach, gamma)?;
    let (shock, nu_deg, mu_deg) = if mach >= 1.0 {
        let s = normal_shock(mach, gamma)?;
        let nu = prandtl_meyer_angle(mach, gamma)?.to_degrees();
        let mu = (1.0 / mach).asin().to_degrees();
        (Some(s), Some(nu), Some(mu))
    } else {
        (None, None, None)
    };
    Ok(FlowReport {
        iso,
        area_mach,
        fanno,
        rayleigh,
        shock,
        nu_deg,
        mu_deg,
    })
}

/// Validate the form, evaluate the relations and format the monospace readout.
/// Extracted from the draw closure so it is unit-testable. On any domain error
/// (`M <= 0`, `gamma <= 1`, non-finite) it sets `error` and clears `result`.
fn run_gasdynamics(s: &mut GasDynamicsWorkbenchState) {
    s.error = None;
    match compute_report(s.mach, s.gamma) {
        Err(e) => {
            s.error = Some(e.to_string());
            s.result.clear();
        }
        Ok(r) => {
            let regime = if s.mach > 1.0 {
                "supersonic"
            } else if s.mach < 1.0 {
                "subsonic"
            } else {
                "sonic (M = 1)"
            };

            let mut out = format!(
                "Mach M             : {:.4}\n\
                 gamma              : {:.4}\n\
                 flow regime        : {regime}\n\n\
                 ISENTROPIC (stagnation / static)\n\
                 \x20 T0/T             : {:.6}\n\
                 \x20 p0/p             : {:.6}\n\
                 \x20 rho0/rho         : {:.6}\n\
                 \x20 A/A*             : {:.6}\n",
                s.mach, s.gamma, r.iso.t0_over_t, r.iso.p0_over_p, r.iso.rho0_over_rho, r.area_mach,
            );

            if let (Some(sh), Some(nu), Some(mu)) = (r.shock, r.nu_deg, r.mu_deg) {
                out.push_str(&format!(
                    "\nNORMAL SHOCK (M1 = {:.4})\n\
                     \x20 M2               : {:.6}\n\
                     \x20 p2/p1            : {:.6}\n\
                     \x20 T2/T1            : {:.6}\n\
                     \x20 rho2/rho1        : {:.6}\n\
                     \x20 p02/p01          : {:.6}\n\
                     \nPRANDTL-MEYER\n\
                     \x20 nu(M)            : {:.4} deg\n\
                     \x20 Mach angle mu    : {:.4} deg\n",
                    s.mach,
                    sh.downstream_mach,
                    sh.pressure_ratio,
                    sh.temperature_ratio,
                    sh.density_ratio,
                    sh.stagnation_pressure_ratio,
                    nu,
                    mu,
                ));
            } else {
                out.push_str("\n(subsonic — no normal shock or Prandtl-Meyer expansion fan)\n");
            }

            out.push_str(&format!(
                "\nFANNO (adiabatic + friction, sonic-ref)\n\
                 \x20 T/T*             : {:.6}\n\
                 \x20 p/p*             : {:.6}\n\
                 \x20 rho/rho*         : {:.6}\n\
                 \x20 V/V*             : {:.6}\n\
                 \x20 p0/p0*           : {:.6}\n\
                 \x20 4fL*/D           : {:.6}\n\
                 \nRAYLEIGH (heat addition, sonic-ref)\n\
                 \x20 T/T*             : {:.6}\n\
                 \x20 p/p*             : {:.6}\n\
                 \x20 rho/rho*         : {:.6}\n\
                 \x20 V/V*             : {:.6}\n\
                 \x20 T0/T0*           : {:.6}\n\
                 \x20 p0/p0*           : {:.6}\n",
                r.fanno.temperature_ratio,
                r.fanno.pressure_ratio,
                r.fanno.density_ratio,
                r.fanno.velocity_ratio,
                r.fanno.stagnation_pressure_ratio,
                r.fanno.friction_length,
                r.rayleigh.temperature_ratio,
                r.rayleigh.pressure_ratio,
                r.rayleigh.density_ratio,
                r.rayleigh.velocity_ratio,
                r.rayleigh.stagnation_temperature_ratio,
                r.rayleigh.stagnation_pressure_ratio,
            ));

            s.result = out;
        }
    }
}

/// Draw the Gas Dynamics right-side panel. A no-op when the
/// `show_gasdynamics_workbench` toggle is off.
pub fn draw_gasdynamics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_gasdynamics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_gasdynamics_workbench",
        "Gas Dynamics",
        |app, ui| {
            ui.label(
                egui::RichText::new("1-D compressible-flow relations · valenx-gasdynamics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.gasdynamics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Inputs").strong());
                    let mut changed = false;
                    ui.horizontal(|ui| {
                        ui.label("Mach M");
                        changed |= ui
                            .add(egui::DragValue::new(&mut s.mach).speed(0.01))
                            .changed();
                    });
                    ui.horizontal(|ui| {
                        ui.label("gamma");
                        changed |= ui
                            .add(egui::DragValue::new(&mut s.gamma).speed(0.005))
                            .changed();
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        changed = true;
                    }

                    // Live recompute on any input change, and once on first
                    // open so the panel is never blank.
                    if changed || (s.result.is_empty() && s.error.is_none()) {
                        run_gasdynamics(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_gasdynamics_workbench = false;
    }
}

/// Build the **Gas Dynamics** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine compressible-flow ratios ([`run_gasdynamics`]) for the canonical
/// default condition (M = 2.0, gamma = 1.4 — the supersonic NACA-1135 worked
/// example: isentropic, normal-shock, Prandtl-Meyer, Fanno and Rayleigh ratios).
/// Registered as the `"gasdynamics"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view.
pub(crate) fn gasdynamics_product() -> crate::WorkspaceProduct {
    let mut s = GasDynamicsWorkbenchState::default();
    run_gasdynamics(&mut s);
    crate::WorkspaceProduct {
        title: "Gas Dynamics".into(),
        lines: crate::products_registry::lines_from_readout(&s.result),
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn default_state_is_idle() {
        let s = GasDynamicsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
        assert!(close(s.mach, 2.0, 1e-12));
        assert!(close(s.gamma, 1.4, 1e-12));
    }

    #[test]
    fn supersonic_report_matches_naca_1135() {
        // gamma = 1.4, M = 2 — the canonical NACA-1135 worked state.
        let r = compute_report(2.0, 1.4).unwrap();
        // Isentropic: T0/T = 1.8 exact, A/A* = 1.6875 exact.
        assert!(close(r.iso.t0_over_t, 1.8, 1e-9));
        assert!(close(r.area_mach, 1.6875, 1e-9));
        // Normal shock at M1 = 2: M2 = 0.57735, p2/p1 = 4.5, p02/p01 = 0.72087.
        let sh = r.shock.expect("M >= 1 has a normal shock");
        assert!(
            close(sh.downstream_mach, 0.577_35, 1e-4),
            "M2 {}",
            sh.downstream_mach
        );
        assert!(
            close(sh.pressure_ratio, 4.5, 1e-9),
            "p2/p1 {}",
            sh.pressure_ratio
        );
        assert!(
            close(sh.stagnation_pressure_ratio, 0.720_87, 1e-4),
            "p02/p01 {}",
            sh.stagnation_pressure_ratio
        );
        // Prandtl-Meyer: nu(2) = 26.380 deg, Mach angle = 30 deg.
        assert!(
            close(r.nu_deg.unwrap(), 26.3798, 1e-3),
            "nu {}",
            r.nu_deg.unwrap()
        );
        assert!(
            close(r.mu_deg.unwrap(), 30.0, 1e-6),
            "mu {}",
            r.mu_deg.unwrap()
        );
        // Fanno + Rayleigh sonic-ref values (cross-check a couple).
        assert!(
            close(r.fanno.pressure_ratio, 0.4082, 1e-3),
            "fanno p/p* {}",
            r.fanno.pressure_ratio
        );
        assert!(
            close(r.rayleigh.stagnation_temperature_ratio, 0.7934, 1e-3),
            "rayleigh T0/T0* {}",
            r.rayleigh.stagnation_temperature_ratio
        );
    }

    #[test]
    fn subsonic_has_no_shock_or_expansion() {
        let r = compute_report(0.5, 1.4).unwrap();
        assert!(r.shock.is_none());
        assert!(r.nu_deg.is_none());
        assert!(r.mu_deg.is_none());
        // The subsonic-defined relations are still present and sane.
        assert!(r.iso.t0_over_t > 1.0);
        assert!(r.fanno.pressure_ratio > 1.0); // p/p* > 1 below sonic
        assert!(r.rayleigh.pressure_ratio > 1.0);
    }

    #[test]
    fn run_formats_a_supersonic_readout() {
        let mut s = GasDynamicsWorkbenchState::default();
        run_gasdynamics(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("supersonic"));
        assert!(s.result.contains("ISENTROPIC"));
        assert!(s.result.contains("NORMAL SHOCK"));
        assert!(s.result.contains("PRANDTL-MEYER"));
        assert!(s.result.contains("FANNO"));
        assert!(s.result.contains("RAYLEIGH"));
    }

    #[test]
    fn run_subsonic_omits_shock_section() {
        let mut s = GasDynamicsWorkbenchState {
            mach: 0.5,
            ..Default::default()
        };
        run_gasdynamics(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("subsonic"));
        assert!(!s.result.contains("NORMAL SHOCK"));
        assert!(s.result.contains("FANNO"));
        assert!(s.result.contains("RAYLEIGH"));
    }

    #[test]
    fn run_rejects_bad_inputs() {
        // gamma must be > 1.
        let mut s = GasDynamicsWorkbenchState {
            gamma: 1.0,
            ..Default::default()
        };
        run_gasdynamics(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());

        // Mach must be > 0.
        let mut s = GasDynamicsWorkbenchState {
            mach: 0.0,
            ..Default::default()
        };
        run_gasdynamics(&mut s);
        assert!(s.error.is_some());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_gasdynamics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_gasdynamics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gasdynamics_workbench = true;
        // First draw auto-computes the default supersonic report.
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_an_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gasdynamics_workbench = true;
        app.gasdynamics.gamma = 1.0; // invalid → error path
        draw_workbench(&mut app);
    }
}

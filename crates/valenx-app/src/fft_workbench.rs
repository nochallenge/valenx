//! The right-side **FFT / Spectrum Workbench** panel — native discrete
//! Fourier transform of a synthesized test tone over `valenx-fft`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fft_workbench`,
//! toggled from the View menu. The form sets a sample count `N`, a sample
//! rate `fs` and a tone frequency `f0` (with amplitude); "Analyze"
//! synthesizes `x[n] = A sin(2 pi f0 n / fs)` in the workbench (a plain
//! `Vec<f64>`), runs the [`dft_real`] direct-summation transform, and
//! reports the bin resolution, Nyquist frequency, peak magnitude bin
//! (which lands at `f0`) and a few leading magnitude-spectrum values.
//! "Show 3-D" loads a representative bar-spectrum solid (one bar per
//! magnitude bin) into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fft::{bin_frequency, dft_real, magnitude_spectrum, nyquist_frequency};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the FFT / Spectrum Workbench.
pub struct FftWorkbenchState {
    /// Number of samples `N` in the synthesized signal.
    n_samples: usize,
    /// Sample rate `fs` (Hz).
    fs_hz: f64,
    /// Tone frequency `f0` (Hz) of the synthesized sinusoid.
    tone_hz: f64,
    /// Tone amplitude `A`.
    amplitude: f64,
    /// Formatted spectrum readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bar-spectrum solid (serviced
    /// after the panel draws).
    show_3d_request: bool,
}

impl Default for FftWorkbenchState {
    fn default() -> Self {
        // A 64-sample window at fs = 1000 Hz carrying a single 125 Hz
        // tone: the bin resolution is fs/N = 15.625 Hz, so the tone
        // lands exactly on bin 8 (8 * 15.625 = 125 Hz) with magnitude
        // A*N/2 = 32.
        Self {
            n_samples: 64,
            fs_hz: 1000.0,
            tone_hz: 125.0,
            amplitude: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the FFT / Spectrum Workbench right-side panel. A no-op when the
/// `show_fft_workbench` toggle is off.
pub fn draw_fft_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fft_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fft_workbench",
        "FFT / Spectrum",
        |app, ui| {
            ui.label(
                egui::RichText::new("native DFT of a synthesized test tone · valenx-fft")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fft;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Signal window").strong());
                    ui.horizontal(|ui| {
                        ui.label("samples N");
                        ui.add(
                            egui::DragValue::new(&mut s.n_samples)
                                .speed(1.0)
                                .range(1..=4096),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("sample rate fs (Hz)");
                        ui.add(egui::DragValue::new(&mut s.fs_hz).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Test tone").strong());
                    ui.horizontal(|ui| {
                        ui.label("frequency f0 (Hz)");
                        ui.add(egui::DragValue::new(&mut s.tone_hz).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("amplitude A");
                        ui.add(egui::DragValue::new(&mut s.amplitude).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fft(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative bar-spectrum (one bar per magnitude-spectrum bin, heights proportional to |X[k]|) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Spectrum").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fft_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fft` borrow is
    // released here): build the bar-spectrum 3-D solid and load it.
    if app.fft.show_3d_request {
        app.fft.show_3d_request = false;
        load_spectrum_3d(app);
    }
}

/// Validate the form, run the transform and format the readout.
fn run_fft(s: &mut FftWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Synthesize the test signal `x[n] = A sin(2 pi f0 n / fs)` for the
/// current form. The crate does the DFT; this signal synthesis is done
/// here in the workbench as a plain `Vec<f64>`. Extracted so it is
/// unit-testable and shared with the 3-D gate.
fn synth_signal(s: &FftWorkbenchState) -> Vec<f64> {
    (0..s.n_samples)
        .map(|n| s.amplitude * (TAU * s.tone_hz * (n as f64) / s.fs_hz).sin())
        .collect()
}

/// The magnitude spectrum of the synthesized signal together with the
/// derived bin resolution and Nyquist frequency — the quantities both
/// the readout and the 3-D gate need. Maps any domain error to a display
/// string. Extracted so it is unit-testable and shared.
fn spectrum(s: &FftWorkbenchState) -> Result<(Vec<f64>, f64, f64), String> {
    let x = synth_signal(s);
    let transform = dft_real(&x).map_err(|e| e.to_string())?;
    let mag = magnitude_spectrum(&transform);
    let nyquist = nyquist_frequency(s.fs_hz).map_err(|e| e.to_string())?;
    let delta_f = s.fs_hz / (s.n_samples as f64);
    Ok((mag, delta_f, nyquist))
}

/// Synthesize, transform and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &FftWorkbenchState) -> Result<String, String> {
    let (mag, delta_f, nyquist) = spectrum(s)?;

    // Peak magnitude bin over the physically meaningful first half
    // [0, N/2] (bins above N/2 are the conjugate mirror of a real
    // signal). Fall back to the whole spectrum for very short windows.
    let half = (mag.len() / 2).max(1);
    let (peak_bin, &peak_mag) = mag[..half]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .ok_or_else(|| "empty spectrum".to_string())?;
    let peak_freq = bin_frequency(peak_bin, mag.len(), s.fs_hz).map_err(|e| e.to_string())?;

    // A few leading magnitude-spectrum values (DC onward).
    let lead = mag.len().min(6);
    let leading: String = mag[..lead]
        .iter()
        .enumerate()
        .map(|(k, m)| format!("  |X[{k}]| = {m:.3}"))
        .collect::<Vec<_>>()
        .join("\n");

    Ok(format!(
        "samples N       : {n}\n\
         sample rate fs  : {fs:.1} Hz\n\
         tone f0 / A     : {f0:.1} Hz / {a:.2}\n\n\
         bin resolution  : {df:.3} Hz\n\
         Nyquist fs/2    : {nyq:.1} Hz\n\
         peak bin        : {pk}\n\
         peak frequency  : {pf:.1} Hz\n\
         peak magnitude  : {pm:.3}\n\n\
         leading |X[k]|  :\n{leading}",
        n = s.n_samples,
        fs = s.fs_hz,
        f0 = s.tone_hz,
        a = s.amplitude,
        df = delta_f,
        nyq = nyquist,
        pk = peak_bin,
        pf = peak_freq,
        pm = peak_mag,
    ))
}

/// Append a single outward-facing bar (an axis-aligned box of half-extents
/// `h` centred at `c`) to the node / triangle buffers.
fn push_bar(
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

/// Build the magnitude spectrum as a triangle [`Mesh`] — a representative
/// bar chart, one bar per magnitude-spectrum bin (over the meaningful
/// first half), bar heights proportional to `|X[k]|`, laid out along the
/// `x` axis. Representative geometry (not to scale; the spectrum values
/// are the `valenx-fft` result). `None` for an invalid configuration.
fn spectrum_bars_mesh(s: &FftWorkbenchState) -> Option<Mesh> {
    let (mag, _delta_f, _nyquist) = spectrum(s).ok()?;

    // Plot the physically meaningful first half [0, N/2].
    let half = (mag.len() / 2).max(1);
    let bars = &mag[..half];
    let peak = bars.iter().copied().fold(0.0_f64, f64::max);
    // Normalise heights to a fixed plot height; guard a flat (all-zero)
    // spectrum so every bar is at least visible.
    let scale = if peak > 0.0 { 1.0 / peak } else { 0.0 };

    let width = 0.6_f64; // per-bar footprint along x
    let depth = 0.25_f64; // half-depth along y
    let min_h = 0.01_f64; // floor so empty bins still render

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    for (k, &m) in bars.iter().enumerate() {
        let height = (m * scale).clamp(min_h, 1.0);
        let cx = (k as f64) * width;
        // Box sits on the z = 0 plane, growing upward in z.
        push_bar(
            &mut nodes,
            &mut tris,
            Vector3::new(cx, 0.0, height * 0.5),
            Vector3::new(width * 0.4, depth, height * 0.5),
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fft");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bar-spectrum solid and load it into the central viewport.
fn load_spectrum_3d(app: &mut ValenxApp) {
    let Some(mesh) = spectrum_bars_mesh(&app.fft) else {
        app.fft.error =
            Some("signal parameters are invalid — cannot build the 3-D spectrum".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<spectrum>/valenx-fft"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = FftWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_resolution_nyquist_and_peak() {
        let mut s = FftWorkbenchState::default();
        run_fft(&mut s);
        assert!(
            s.error.is_none(),
            "default tone should analyze: {:?}",
            s.error
        );
        // N=64, fs=1000 -> df = 15.625 Hz, Nyquist = 500.0 Hz.
        assert!(s.result.contains("15.625"));
        assert!(s.result.contains("500.0"));
        // f0 = 125 Hz lands on bin 8 -> peak frequency reads 125.0 Hz.
        assert!(s.result.contains("peak bin        : 8"));
        assert!(s.result.contains("125.0"));
    }

    #[test]
    fn analyze_rejects_zero_sample_rate() {
        let mut s = FftWorkbenchState {
            fs_hz: 0.0,
            ..Default::default()
        };
        run_fft(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_pure_tone_peaks_at_f0_bin() {
        // Hand-computed ground truth: a pure tone at f0 = 125 Hz sampled
        // at fs = 1000 Hz over N = 64 samples has bin resolution
        // df = fs/N = 15.625 Hz, so the tone sits exactly on
        // bin = round(f0/df) = round(125/15.625) = 8, whose physical
        // frequency bin_frequency(8) = 8 * 1000 / 64 = 125.0 Hz. A real
        // sine of amplitude A on a bin has magnitude A*N/2 = 32.
        let s = FftWorkbenchState::default();
        let (mag, delta_f, nyquist) = spectrum(&s).expect("default tone analyzes");
        assert!((delta_f - 15.625).abs() < 1e-12);
        assert!((nyquist - 500.0).abs() < 1e-12);

        let half = mag.len() / 2;
        let (peak_bin, &peak_mag) = mag[..half]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        assert_eq!(peak_bin, 8);
        assert!((bin_frequency(peak_bin, mag.len(), s.fs_hz).unwrap() - 125.0).abs() < 1e-12);
        assert!((peak_mag - 32.0).abs() < 1e-7);
    }

    #[test]
    fn synth_signal_has_expected_length_and_zero_first_sample() {
        // x[0] = A sin(0) = 0, and the window holds exactly N samples.
        let s = FftWorkbenchState::default();
        let x = synth_signal(&s);
        assert_eq!(x.len(), 64);
        assert!((x[0] - 0.0_f64).abs() < 1e-12);
    }

    #[test]
    fn spectrum_bars_mesh_for_default_is_nonempty_and_in_range() {
        let s = FftWorkbenchState::default();
        let mesh = spectrum_bars_mesh(&s).expect("default tone yields a bar mesh");
        assert!(mesh.nodes.len() > 8, "expected one box per plotted bin");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fft_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fft_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fft_workbench = true;
        run_fft(&mut app.fft);
        draw_workbench(&mut app);
    }
}

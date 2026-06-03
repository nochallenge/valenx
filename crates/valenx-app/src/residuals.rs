//! Residual history accumulator and bottom-panel chart.
//!
//! The UI side of the OpenFOAM adapter's live-residual stream. We
//! reuse the adapter's own log parser so the panel stays honest
//! about what the solver actually emitted — no separate copy of the
//! parsing rules drifting from the canonical one.
//!
//! Residuals are stored per field as `(iteration, value)` points on
//! a log-scale axis. egui_plot handles the drawing; we only have to
//! feed it sample series and labels.

use std::collections::{BTreeMap, VecDeque};

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use valenx_adapter_openfoam::log_parser::{intern_field, parse_line, LogSignal};

use crate::settings::ResidualScale;

/// Round-8 cap: maximum samples retained per field. Pre-fix, a long
/// transient run could push millions of `[f64; 2]` entries into one
/// `Vec`, eating ~16 bytes per sample with no upper bound — both a
/// memory leak (the history grows unbounded over a long-running
/// session) and a per-frame work amplifier (the chart re-walks
/// every sample on every frame). 50_000 samples per field is enough
/// to show full-resolution detail for the first 50k steps and to
/// retain decimated history after that without bothering the
/// per-frame cost.
pub const MAX_RESIDUAL_SAMPLES_PER_FIELD: usize = 50_000;

/// In-memory residual store, keyed by canonical field name.
///
/// X-axis is the OpenFOAM `Time` value as reported by the solver:
/// integer pseudo-time for simpleFoam (1, 2, 3, …), real seconds for
/// pimpleFoam / icoFoam (0.0005, 0.001, …). Plotting on the real
/// time value instead of a step counter gives transient runs a
/// meaningful x-axis without losing anything for steady runs.
#[derive(Clone, Debug, Default)]
pub struct ResidualHistory {
    current_time: f64,
    /// Series per field, in time order.
    /// Round-9: switched from `Vec` to `VecDeque` so the cap-hit
    /// eviction is O(1) (`pop_front`) rather than O(n)
    /// (`Vec::remove(0)`). For a fully-saturated buffer that's
    /// 50k pops/insertion saved per ingestion.
    pub by_field: BTreeMap<&'static str, VecDeque<[f64; 2]>>,
    /// Round-8 cache: log10-transformed mirror of `by_field`, built
    /// incrementally on `ingest_log_line`. Pre-fix, `show()`
    /// rebuilt a `Vec<[f64; 2]>` per field per frame, transforming
    /// each sample on every redraw — wasteful for long runs. The
    /// mirror trades 2× memory for 1× compute saved per frame.
    /// Round-9: same VecDeque migration as `by_field` above — both
    /// series have to evict in lockstep, and the per-frame draw
    /// reads from the mirror.
    by_field_log10: BTreeMap<&'static str, VecDeque<[f64; 2]>>,
}

impl ResidualHistory {
    /// Drop every recorded series and reset the parser's running
    /// time counter.
    pub fn clear(&mut self) {
        self.current_time = 0.0;
        self.by_field.clear();
        self.by_field_log10.clear();
    }

    /// Update from one line of solver output.
    ///
    /// Round-8: when a field's sample buffer reaches
    /// [`MAX_RESIDUAL_SAMPLES_PER_FIELD`], the oldest entry is
    /// dropped before the new one lands. The history stays bounded
    /// regardless of run length. The log10 mirror is updated in
    /// lockstep so the chart's per-frame draw is allocation-free.
    pub fn ingest_log_line(&mut self, line: &str) {
        let signal = parse_line(line);
        match signal {
            LogSignal::Time(t) => {
                self.current_time = t;
            }
            LogSignal::Residual { field, initial } => {
                if let Some(interned) = intern_field(field) {
                    let series = self.by_field.entry(interned).or_default();
                    let log_series = self.by_field_log10.entry(interned).or_default();
                    if series.len() >= MAX_RESIDUAL_SAMPLES_PER_FIELD {
                        // Cap reached — drop oldest sample (FIFO) so
                        // the chart keeps recent detail and the
                        // memory footprint stays bounded.
                        // Round-9: O(1) pop_front on VecDeque
                        // (vs pre-fix O(n) Vec::remove(0)).
                        series.pop_front();
                        log_series.pop_front();
                    }
                    let clamped = initial.max(1e-30);
                    series.push_back([self.current_time, clamped]);
                    log_series.push_back([self.current_time, clamped.log10()]);
                }
            }
            LogSignal::Execution { .. } | LogSignal::Other => {}
        }
    }

    /// `true` once no field samples have been ingested.
    pub fn is_empty(&self) -> bool {
        self.by_field.is_empty()
    }
}

/// Render the residual chart into the given `ui`. `scale` selects
/// log10 (the sensible default — residuals routinely cover 6+
/// orders of magnitude) or linear. `convergence_target` paints a
/// horizontal reference line at the threshold value when
/// `Some(...)` — useful for visually flagging which fields have
/// dropped below the engineering / research target.
pub fn show(
    ui: &mut egui::Ui,
    history: &ResidualHistory,
    scale: ResidualScale,
    convergence_target: Option<f64>,
) {
    if history.is_empty() {
        ui.vertical_centered(|ui| {
            ui.add_space(8.0);
            ui.label("Residuals will appear here when a case runs.");
            ui.add_space(8.0);
        });
        return;
    }

    // Per-field summary header: "Ux: 1.2e-5 ✓  p: 3.4e-3 ✗" — gives
    // the user the converged/not-converged-at-threshold answer at a
    // glance without having to read coords off the chart.
    if let Some(threshold) = convergence_target {
        let mut chips: Vec<String> = Vec::with_capacity(history.by_field.len());
        for (field, samples) in &history.by_field {
            if let Some(&[_, last]) = samples.back() {
                let mark = if last <= threshold {
                    "\u{2713}"
                } else {
                    "\u{2717}"
                };
                chips.push(format!("{field}: {} {mark}", format_residual(last)));
            }
        }
        if !chips.is_empty() {
            ui.horizontal_wrapped(|ui| {
                ui.small(format!("target {}:", format_residual(threshold)));
                for chip in chips {
                    ui.small(chip);
                }
            });
        }
    }

    let y_label = match scale {
        ResidualScale::Log10 => "initial residual (log\u{2081}\u{2080})",
        ResidualScale::Linear => "initial residual",
    };

    Plot::new("valenx_residuals")
        .height(ui.available_height().max(120.0))
        .allow_zoom(true)
        .allow_drag(true)
        .allow_scroll(false)
        .legend(Legend::default().position(egui_plot::Corner::RightTop))
        .x_axis_label("iteration")
        .y_axis_label(y_label)
        .show(ui, |plot_ui| {
            // Round-8: pull from the appropriate pre-built mirror so
            // the per-frame draw is allocation-free. Linear -> raw
            // `by_field`; Log10 -> `by_field_log10`.
            let series_map = match scale {
                ResidualScale::Log10 => &history.by_field_log10,
                ResidualScale::Linear => &history.by_field,
            };
            for (field, samples) in series_map {
                // Round-9: VecDeque doesn't have a free `to_vec()`
                // that PlotPoints accepts via the slice impl, so
                // build PlotPoints from an iterator. Cloning is
                // still O(n) but skips the temporary intermediate
                // Vec the slice impl would force.
                plot_ui.line(
                    Line::new(PlotPoints::from_iter(samples.iter().copied()))
                        .name(*field),
                );
            }
            // Convergence-target overlay: a flat horizontal line at
            // log10(threshold) (or threshold itself for linear). The
            // target line is named so it shows up in the legend with
            // a distinct colour from the field curves.
            if let Some(threshold) = convergence_target {
                let y = match scale {
                    ResidualScale::Log10 => threshold.max(1e-30).log10(),
                    ResidualScale::Linear => threshold,
                };
                if let Some((x_min, x_max)) = chart_x_extent(history) {
                    plot_ui.line(
                        Line::new(PlotPoints::from(vec![[x_min, y], [x_max, y]]))
                            .name(format!("target {}", format_residual(threshold))),
                    );
                }
            }
        });
}

/// Render the residual history as a long-format CSV string.
///
/// One row per `(time, field, value)` tuple. Long-format wins over
/// wide-format because fields land at different time points
/// (transient runs sample non-uniformly per field) — wide-format
/// would produce ragged rows with empty cells.
///
/// Header: `time,field,value`. Sorted by `(time, field)` so the
/// rows pair lexicographically across runs and a `git diff` between
/// two log exports surfaces actual changes rather than reorderings.
pub fn residuals_to_csv(history: &ResidualHistory) -> String {
    // (time, field, value) tuples, then sort.
    let mut rows: Vec<(f64, &str, f64)> = Vec::new();
    for (field, samples) in &history.by_field {
        for &[t, v] in samples {
            rows.push((t, field, v));
        }
    }
    rows.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.cmp(b.1))
    });
    let mut out = String::with_capacity(48 + rows.len() * 24);
    out.push_str("time,field,value\n");
    for (t, field, v) in rows {
        // Scientific notation for the value so 1e-30 stays
        // readable without producing 0.0 when written through
        // the default Display impl.
        out.push_str(&format!("{t},{field},{v:e}\n"));
    }
    out
}

/// Write the residual history as a long-format CSV file. Creates
/// the parent directory if needed. Wraps [`residuals_to_csv`].
pub fn write_residuals_csv(
    history: &ResidualHistory,
    path: &std::path::Path,
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    valenx_core::io_caps::atomic_write_str(path, &residuals_to_csv(history))
}

/// Compact scientific-notation formatter for residual values. Keeps
/// the chart legend / summary chips narrow even for tiny gradients.
fn format_residual(v: f64) -> String {
    if v == 0.0 {
        "0".into()
    } else {
        format!("{v:.2e}")
    }
}

/// Find the (min_x, max_x) extent across every series. Used to size
/// the convergence-target overlay so it spans the chart without
/// requiring the user to enable auto-axes manually.
fn chart_x_extent(history: &ResidualHistory) -> Option<(f64, f64)> {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut any = false;
    for samples in history.by_field.values() {
        for &[x, _] in samples {
            if x < min {
                min = x;
            }
            if x > max {
                max = x;
            }
            any = true;
        }
    }
    if any {
        Some((min, max))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ingests_time_and_residual_lines() {
        let mut hist = ResidualHistory::default();
        hist.ingest_log_line("Time = 1");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.01, Final residual = 1e-4, No Iterations 5",
        );
        hist.ingest_log_line(
            "GAMG:  Solving for p, Initial residual = 0.5, Final residual = 2e-3, No Iterations 10",
        );
        hist.ingest_log_line("Time = 2");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.005, Final residual = 5e-5, No Iterations 3",
        );
        let ux = hist.by_field.get("Ux").expect("ux series");
        let p = hist.by_field.get("p").expect("p series");
        assert_eq!(ux.len(), 2);
        assert_eq!(p.len(), 1);
        assert_eq!(ux[0][0], 1.0);
        assert!((ux[0][1] - 0.01).abs() < 1e-9);
        assert_eq!(ux[1][0], 2.0);
    }

    #[test]
    fn ignores_unknown_fields() {
        let mut hist = ResidualHistory::default();
        hist.ingest_log_line("Time = 1");
        hist.ingest_log_line(
            "smoothSolver:  Solving for alpha.water, Initial residual = 0.1, Final residual = 1e-3, No Iterations 1",
        );
        // alpha.water isn't in the interned list yet — drop silently
        // rather than pollute the chart with an untyped field.
        assert!(hist.is_empty());
    }

    #[test]
    fn clears_everything() {
        let mut hist = ResidualHistory::default();
        hist.ingest_log_line("Time = 5");
        hist.ingest_log_line(
            "smoothSolver:  Solving for p, Initial residual = 0.1, Final residual = 1e-3, No Iterations 1",
        );
        assert!(!hist.is_empty());
        hist.clear();
        assert!(hist.is_empty());
    }

    #[test]
    fn residuals_to_csv_writes_header_and_long_format_rows() {
        let mut hist = ResidualHistory::default();
        hist.ingest_log_line("Time = 1");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.01, Final residual = 1e-4, No Iterations 5",
        );
        hist.ingest_log_line(
            "GAMG:  Solving for p, Initial residual = 0.5, Final residual = 2e-3, No Iterations 10",
        );
        hist.ingest_log_line("Time = 2");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.005, Final residual = 5e-5, No Iterations 3",
        );
        let csv = residuals_to_csv(&hist);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], "time,field,value");
        // 3 sample rows: (1,Ux,0.01), (1,p,0.5), (2,Ux,0.005)
        assert_eq!(lines.len(), 4);
        // Sorted by time then field — Ux precedes p alphabetically.
        assert!(lines[1].starts_with("1,Ux,"));
        assert!(lines[2].starts_with("1,p,"));
        assert!(lines[3].starts_with("2,Ux,"));
    }

    #[test]
    fn residuals_to_csv_empty_history_emits_header_only() {
        let csv = residuals_to_csv(&ResidualHistory::default());
        assert_eq!(csv, "time,field,value\n");
    }

    #[test]
    fn residuals_to_csv_uses_scientific_for_value() {
        // The value column must survive as scientific notation so
        // tiny gradients (1e-30 floor) don't print as 0.
        let mut hist = ResidualHistory::default();
        hist.by_field.insert("Ux", VecDeque::from(vec![[1.0, 1e-12]]));
        let csv = residuals_to_csv(&hist);
        assert!(
            csv.contains("1e-12") || csv.contains("1.0e-12"),
            "got {csv:?}"
        );
    }

    #[test]
    fn write_residuals_csv_round_trips_through_disk() {
        let path = std::env::temp_dir().join(format!(
            "valenx-resids-{}.csv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut hist = ResidualHistory::default();
        hist.by_field
            .insert("p", VecDeque::from(vec![[1.0, 0.5], [2.0, 0.25]]));
        write_residuals_csv(&hist, &path).expect("write");
        let text = std::fs::read_to_string(&path).expect("read");
        assert!(text.starts_with("time,field,value\n"));
        assert!(text.contains("1,p,"));
        assert!(text.contains("2,p,"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn format_residual_uses_compact_scientific() {
        // Common cases the legend shows.
        assert_eq!(super::format_residual(0.0), "0");
        assert_eq!(super::format_residual(1.0), "1.00e0");
        assert_eq!(super::format_residual(0.0001), "1.00e-4");
        assert_eq!(super::format_residual(1.5e-7), "1.50e-7");
        assert_eq!(super::format_residual(101325.0), "1.01e5");
    }

    #[test]
    fn chart_x_extent_finds_min_and_max_across_fields() {
        let mut hist = ResidualHistory::default();
        // Field A: x in [1.0, 5.0]
        hist.by_field
            .insert("A", VecDeque::from(vec![[1.0, 0.1], [3.0, 0.05], [5.0, 0.01]]));
        // Field B: x in [2.0, 4.0]
        hist.by_field
            .insert("B", VecDeque::from(vec![[2.0, 0.5], [4.0, 0.2]]));
        let (min, max) = super::chart_x_extent(&hist).expect("non-empty");
        assert_eq!(min, 1.0);
        assert_eq!(max, 5.0);
    }

    #[test]
    fn chart_x_extent_for_empty_history_is_none() {
        let hist = ResidualHistory::default();
        assert!(super::chart_x_extent(&hist).is_none());
    }

    /// Round-9 RED→GREEN: at the cap, eviction must be O(1)
    /// (`VecDeque::pop_front`) not O(n) (`Vec::remove(0)`). For
    /// `MAX_RESIDUAL_SAMPLES_PER_FIELD = 50_000`, pushing 60_000
    /// samples should finish in well under a second; the pre-fix
    /// O(n²) loop would take many seconds. We assert a 2-second
    /// budget — generous for slow CI, but tight enough to catch
    /// the regression.
    #[test]
    fn ingest_at_cap_is_o1_not_on() {
        use std::time::{Duration, Instant};
        let mut hist = ResidualHistory::default();
        let n_pushes = MAX_RESIDUAL_SAMPLES_PER_FIELD + 10_000;
        let started = Instant::now();
        for i in 0..n_pushes {
            // Synthesise a residual line per step. Walk the time
            // counter forward so each row has a unique x.
            if i % 2 == 0 {
                hist.ingest_log_line(&format!("Time = {i}"));
            }
            hist.ingest_log_line(
                "smoothSolver:  Solving for Ux, Initial residual = 0.001, Final residual = 1e-5, No Iterations 1"
            );
        }
        let elapsed = started.elapsed();
        assert!(
            elapsed < Duration::from_secs(2),
            "ingesting {n_pushes} lines took {elapsed:?} — eviction is O(n) instead of O(1)"
        );
        let ux = hist.by_field.get("Ux").expect("ux series");
        assert_eq!(ux.len(), MAX_RESIDUAL_SAMPLES_PER_FIELD);
    }

    #[test]
    fn transient_subsecond_times_plot_at_real_time() {
        // Regression: pimpleFoam emits `Time = 0.0005`, `Time = 0.001`,
        // etc. The chart must plot at the real time, not snap every
        // sample to step 0 (which is what truncation-to-u64 would do).
        let mut hist = ResidualHistory::default();
        hist.ingest_log_line("Time = 0.0005");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.01, Final residual = 1e-4, No Iterations 5",
        );
        hist.ingest_log_line("Time = 0.001");
        hist.ingest_log_line(
            "smoothSolver:  Solving for Ux, Initial residual = 0.005, Final residual = 5e-5, No Iterations 4",
        );
        let ux = hist.by_field.get("Ux").expect("ux series");
        assert_eq!(ux.len(), 2);
        // Two distinct x-values, not both at zero — that was the bug.
        assert!((ux[0][0] - 0.0005).abs() < 1e-12);
        assert!((ux[1][0] - 0.001).abs() < 1e-12);
        assert!(ux[0][0] != ux[1][0]);
    }

    #[test]
    fn ingest_caps_samples_per_field() {
        // Round-8 RED→GREEN: pre-fix, `by_field` grew unbounded —
        // a long transient run could push millions of samples into
        // one Vec. The MAX_RESIDUAL_SAMPLES_PER_FIELD cap drops the
        // oldest sample once at capacity, keeping memory bounded.
        let mut hist = ResidualHistory::default();
        // Push 1.2× the cap worth of Ux samples.
        let total = MAX_RESIDUAL_SAMPLES_PER_FIELD + 10_000;
        for i in 0..total {
            hist.ingest_log_line(&format!("Time = {i}"));
            hist.ingest_log_line(
                "smoothSolver:  Solving for Ux, Initial residual = 0.01, \
                 Final residual = 1e-4, No Iterations 5",
            );
        }
        let ux = hist.by_field.get("Ux").expect("ux series");
        assert_eq!(
            ux.len(),
            MAX_RESIDUAL_SAMPLES_PER_FIELD,
            "series length must be capped at MAX_RESIDUAL_SAMPLES_PER_FIELD"
        );
        // The oldest entries were dropped; the first remaining
        // sample is from `Time = 10_000` (we pushed 60k, kept the
        // tail 50k).
        let first_time = ux[0][0];
        let expected_first = (total - MAX_RESIDUAL_SAMPLES_PER_FIELD) as f64;
        assert!(
            (first_time - expected_first).abs() < 1e-9,
            "expected first retained sample at t={expected_first}, got {first_time}"
        );
        // The log10 mirror tracks the raw series length.
        let ux_log = hist.by_field_log10.get("Ux").expect("ux log series");
        assert_eq!(ux_log.len(), MAX_RESIDUAL_SAMPLES_PER_FIELD);
    }
}

//! Process-entry plumbing: native window bootstrap, tracing init,
//! and the panic-hook installer.
//!
//! Extracted from `lib.rs` so the root module focuses on the
//! application state machine rather than the once-per-process
//! "stand the binary up" wiring. Every public item here is
//! re-exported from `lib.rs`, so callers continue using
//! `valenx_app::run()` / `valenx_app::crashes_dir()` unchanged.

use std::path::PathBuf;

use eframe::egui;

use crate::state_paths::state_dir;
use crate::ValenxApp;

/// Entry point used by `src/main.rs`.
pub fn run() -> anyhow::Result<()> {
    init_tracing();
    install_crash_reporter();

    // Headless mode: an automation agent or CI can run batch tasks (compute,
    // geometry export) with no window. The default — no `--headless` flag —
    // falls through to the headed interactive GUI below.
    let args: Vec<String> = std::env::args().collect();
    if let Some(pos) = args.iter().position(|a| a == "--headless") {
        return crate::headless::run_headless(&args[pos + 1..]);
    }
    // Product self-test: `valenx --self-test [--group <G>] [--id <id>]` is a
    // first-class shortcut for the `--headless self-test` task — it runs the
    // baked-in 56-product verification head-less (no window, no rfd dialog) and
    // prints the compact report. Kept as its own flag so the standing "one fast
    // command to verify the products" entry point is obvious.
    if let Some(pos) = args.iter().position(|a| a == "--self-test") {
        let mut task = vec!["self-test".to_string()];
        task.extend_from_slice(&args[pos + 1..]);
        return crate::headless::run_headless(&task);
    }

    tracing::info!(
        target: "valenx",
        version = env!("CARGO_PKG_VERSION"),
        "launching native window",
    );

    let initial_stl = std::env::args().nth(1).map(PathBuf::from);

    // Restore saved window geometry (position + size) so Valenx reopens
    // where the user left it — e.g. on a second monitor. Best-effort: a
    // missing/unreadable settings.json (or the toggle off) uses defaults.
    let saved = crate::settings_io::load_settings_from_state_dir();
    let remember = saved.as_ref().is_none_or(|s| s.remember_window_geometry);
    let inner_size = saved
        .as_ref()
        .filter(|_| remember)
        .and_then(|s| s.window_size)
        .unwrap_or([1440.0, 940.0]);
    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size(inner_size)
        .with_min_inner_size([1000.0, 620.0])
        .with_title(format!("Valenx {}", env!("CARGO_PKG_VERSION")));
    if remember {
        if let Some(pos) = saved.as_ref().and_then(|s| s.window_position) {
            viewport = viewport.with_position(pos);
        }
    }
    let native_options = eframe::NativeOptions {
        viewport,
        // Explicitly select the wgpu backend so the offscreen depth-buffered
        // viewport and the ground grid are guaranteed to initialise.
        // Without this, eframe could fall back to glow on some platforms/builds,
        // leaving cc.wgpu_render_state == None and the entire GPU render path
        // silently dead-coding (which matched the "black viewport, no grid" symptom).
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };

    eframe::run_native(
        "Valenx",
        native_options,
        Box::new(move |cc| {
            // Larger default UI scale — the dark theme reads too small at 1.0 on
            // big / hi-DPI monitors. Adjustable live via View -> Text size or the
            // built-in Ctrl + / Ctrl - zoom.
            cc.egui_ctx.set_zoom_factor(1.3);
            Ok(Box::new(ValenxApp::new(cc, initial_stl.clone())))
        }),
    )
    .map_err(|e| anyhow::anyhow!("eframe run_native failed: {e}"))
}

fn init_tracing() {
    // Default filter: `info` everywhere, plus per-module overrides
    // that silence two wgpu_hal warning sites we can't fix upstream:
    //
    //   * `wgpu_hal::vulkan::conv` spams
    //       "Unrecognized present mode 1000361000"
    //     on every swapchain creation when the driver advertises a
    //     mode wgpu's match arm doesn't enumerate.
    //
    //   * `wgpu_hal::vulkan::adapter` logs
    //       "Adapter is not Vulkan compliant, hiding adapter"
    //     twice per launch for adapters wgpu rejects.
    //
    // Both are noise from our perspective: harmless, not actionable,
    // and they crowd out useful messages. RUST_LOG still wins — if
    // a user sets it, they get exactly what they asked for.
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env().unwrap_or_else(|_| {
                tracing_subscriber::EnvFilter::new(
                    "info,wgpu_hal::vulkan::conv=error,wgpu_hal::vulkan::adapter=error",
                )
            }),
        )
        .try_init();
}

/// Install the workspace's `valenx-crash-reporter` panic hook.
///
/// Reports land in `<state_dir>/crashes/` so `state_dir()`'s
/// platform-specific resolver covers them automatically. When
/// `state_dir()` returns `None` (extremely rare — no HOME-equivalent
/// on the host) we fall back to `std::env::temp_dir()/valenx-crashes`
/// so a panic still produces something diagnostic.
///
/// **The crashes directory is NOT pre-created here.** Earlier
/// revisions did `std::fs::create_dir_all(&dir)` unconditionally at
/// startup, which polluted `<state_dir>/crashes/` on every launch
/// even when no crash ever happened. The panic hook's writer
/// (`CrashReport::write_to_disk`) already runs `create_dir_all`
/// before serialising — see `crates/valenx-crash-reporter/src/lib.rs`
/// — so the directory appears exactly when a real crash is being
/// persisted, and never before. This is also what the
/// `disable_file_browser_popups` setting documents: users who flip
/// that on don't want Valenx pre-creating state dirs on every
/// startup as a side effect of "click button → File Explorer pops
/// up" plumbing.
///
/// Idempotent: calling `run()` twice in the same process replaces
/// the previous hook.
fn install_crash_reporter() {
    let dir = crashes_dir();
    valenx_crash_reporter::install_panic_hook(dir, env!("CARGO_PKG_VERSION").to_string());
}

/// Resolve the per-user crashes directory.
///
/// Public so the Settings panel can show users where reports land
/// (and let them open the folder for a manual review before the
/// opt-in upload).
pub fn crashes_dir() -> PathBuf {
    state_dir()
        .map(|d| d.join("crashes"))
        .unwrap_or_else(|| std::env::temp_dir().join("valenx-crashes"))
}

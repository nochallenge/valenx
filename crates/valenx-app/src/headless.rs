//! Headless mode — run valenx batch tasks (compute, geometry export) with **no
//! window**, so an automation agent or CI can drive valenx on a machine with no
//! display. Selected by the `--headless` flag; the default stays the headed
//! interactive GUI.
//!
//! This is the user-requested "agents can do headless or headed work" split:
//! `valenx --headless <task>` runs a batch task and exits; `valenx` with no
//! flag opens the window as before.
//!
//! Honest scope: a small, growing set of batch tasks, not the full app. The
//! physics it reports is the same first-order **preliminary-design** model the
//! GUI uses (`valenx_astro`) — not flight-grade.

use std::path::Path;

use valenx_astro::{combust, solve_cycle, CycleInputs, Propellant};

/// Run the headless task named by `task_args` — the CLI tokens that follow
/// `--headless`. Defaults to `info` when no task is given.
///
/// # Errors
///
/// Returns an error for an unknown task, a missing output path, or a failed
/// geometry export.
pub fn run_headless(task_args: &[String]) -> anyhow::Result<()> {
    let task = task_args.first().map(String::as_str).unwrap_or("info");
    let arg1 = task_args.get(1).cloned();
    match task {
        "info" => print_info(),
        "cycle" => print_cycle(),
        "self-test" => return run_self_test(&task_args[1..]),
        "export-engine" => export_stl(&crate::rocket_mesh::detailed_engine_mesh(), arg1)?,
        "export-rocket" => export_stl(&crate::rocket_mesh::lv1_rocket_mesh(), arg1)?,
        "render-engine" => render_png(true, arg1)?,
        "render-rocket" => render_png(false, arg1)?,
        other => anyhow::bail!(
            "unknown headless task `{other}` — try: info | cycle | self-test | \
             export-engine <out.stl> | export-rocket <out.stl> | \
             render-engine <out.png> | render-rocket <out.png>"
        ),
    }
    Ok(())
}

/// Run valenx's baked-in **product self-test** and print the compact report
/// (one `id  PASS|FAIL|SKIP  <key value or reason>` line per product, then a
/// tally). Selected by `valenx --self-test [--group <G>] [--id <id>]` (the
/// `setup::run` arg-router maps the bare `--self-test` flag onto this `self-test`
/// headless task). Runs the [`crate::self_test`] registry entirely head-less —
/// each product check builds a fresh `ValenxApp::default()` and probes its panel
/// in a throwaway egui frame, so **no window and no `rfd` dialog** is ever
/// constructed. Exits non-zero (an `Err`) iff any product **FAILED** (skips are
/// not failures), so CI can gate on it.
fn run_self_test(args: &[String]) -> anyhow::Result<()> {
    let filter = crate::self_test::Filter::from_args(args);
    let report = crate::self_test::run_self_tests(&filter);
    // The report is the product; print it verbatim to stdout.
    print!("{}", report.render());
    if report.ok() {
        Ok(())
    } else {
        anyhow::bail!("self-test: {} product(s) FAILED", report.failed)
    }
}

/// Print the build identity, platform and the available headless tasks.
fn print_info() {
    println!("valenx {} — headless mode", env!("CARGO_PKG_VERSION"));
    println!("  platform : {}", std::env::consts::OS);
    println!(
        "  tasks    : info | cycle | export-engine <out.stl> | export-rocket <out.stl> | \
         render-engine <out.png> | render-rocket <out.png>"
    );
    println!("  headed   : run `valenx` with no --headless flag for the GUI");
}

/// Solve and print the Raptor-class methalox full-flow staged-combustion cycle
/// — a self-contained demonstration of headless compute.
fn print_cycle() {
    let inputs = CycleInputs::raptor_methalox();
    let cyc = solve_cycle(&inputs);
    let comb = combust(
        Propellant::Ch4Lox,
        inputs.mixture_ratio,
        inputs.chamber_pressure / 1.0e5,
    );
    println!("valenx headless — Raptor-class methalox FFSC cycle");
    println!(
        "  combustion  : Tc {:.0} K · gamma {:.3} · M {:.1} g/mol · c* {:.0} m/s",
        comb.chamber_temperature, comb.gamma, comb.molar_mass, comb.c_star
    );
    println!(
        "  chamber p   : {:.0} bar target · {:.0} bar cycle ceiling",
        inputs.chamber_pressure / 1.0e5,
        cyc.max_chamber_pressure / 1.0e5
    );
    println!(
        "  turbopumps  : ox {:.1} MW · fuel {:.1} MW (turbine inlet {:.0} K)",
        cyc.ox.turbine_power / 1.0e6,
        cyc.fuel.turbine_power / 1.0e6,
        inputs.ox.turbine_inlet_temperature
    );
    println!("  cycle closes: {}", if cyc.closes { "yes" } else { "no" });
    println!("  (first-order preliminary-design model — not flight-grade)");
}

/// Export a mesh to a binary STL at `path`.
fn export_stl(mesh: &valenx_mesh::Mesh, path: Option<String>) -> anyhow::Result<()> {
    let path = path.ok_or_else(|| anyhow::anyhow!("export needs an output path, e.g. out.stl"))?;
    valenx_mesh::write_stl_binary(mesh, Path::new(&path))
        .map_err(|e| anyhow::anyhow!("STL export failed: {e}"))?;
    let tris: usize = mesh
        .element_blocks
        .iter()
        .map(|b| b.connectivity.len() / 3)
        .sum();
    println!("wrote {path} ({tris} triangles)");
    Ok(())
}

/// Path-trace the engine or rocket headlessly and write it to a PNG — the same
/// `valenx-pathtrace` render the GUI's Render panel produces, with no window.
fn render_png(engine: bool, path: Option<String>) -> anyhow::Result<()> {
    let path = path.ok_or_else(|| anyhow::anyhow!("render needs an output path, e.g. out.png"))?;
    let (w, h, pixels) = if engine {
        crate::render_workbench::render_engine(480, 160, 6, 1.1)
    } else {
        crate::render_workbench::render_rocket(480, 160, 6, 1.1, false)
    }
    .map_err(|e| anyhow::anyhow!("render failed: {e}"))?;
    // Encode to an in-memory buffer, then publish atomically (sidecar →
    // rename) so a torn write can't leave a corrupt PNG — matches the
    // workspace's durable-write guard (no raw File::create in production).
    let mut buf: Vec<u8> = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut buf, w as u32, h as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header()?;
        writer.write_image_data(&pixels)?;
        writer.finish()?;
    }
    valenx_core::io_caps::atomic_write_bytes(std::path::Path::new(&path), &buf)
        .map_err(|e| anyhow::anyhow!("write {path}: {e}"))?;
    println!("wrote {path} ({w}x{h})");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_cycle_and_default_run_clean() {
        run_headless(&[]).unwrap(); // no task → info
        run_headless(&["info".to_string()]).unwrap();
        run_headless(&["cycle".to_string()]).unwrap();
    }

    #[test]
    fn export_engine_writes_a_binary_stl() {
        let path = std::env::temp_dir().join("valenx_headless_engine_test.stl");
        let p = path.to_string_lossy().into_owned();
        run_headless(&["export-engine".to_string(), p]).unwrap();
        let len = std::fs::metadata(&path).unwrap().len();
        // Binary STL = an 84-byte header + 50 bytes per triangle.
        assert!(len > 84, "STL too small ({len} bytes)");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn missing_path_and_unknown_task_error() {
        assert!(run_headless(&["export-engine".to_string()]).is_err());
        assert!(run_headless(&["render-engine".to_string()]).is_err());
        assert!(run_headless(&["frobnicate".to_string()]).is_err());
    }
}

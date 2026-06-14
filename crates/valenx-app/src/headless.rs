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
        "export-engine" => export_stl(&crate::rocket_mesh::detailed_engine_mesh(), arg1)?,
        "export-rocket" => export_stl(&crate::rocket_mesh::lv1_rocket_mesh(), arg1)?,
        other => anyhow::bail!(
            "unknown headless task `{other}` — try: info | cycle | \
             export-engine <out.stl> | export-rocket <out.stl>"
        ),
    }
    Ok(())
}

/// Print the build identity, platform and the available headless tasks.
fn print_info() {
    println!("valenx {} — headless mode", env!("CARGO_PKG_VERSION"));
    println!("  platform : {}", std::env::consts::OS);
    println!("  tasks    : info | cycle | export-engine <out.stl> | export-rocket <out.stl>");
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
        assert!(run_headless(&["frobnicate".to_string()]).is_err());
    }
}

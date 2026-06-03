//! End-to-end smoke test for the bundled `examples/dock/1iep_demo/`
//! case. Reads the case TOML, runs the native engine, and asserts the
//! output PDBQT has at least one MODEL block with a `VINA RESULT:`
//! score line.
//!
//! Why this test exists: Valenx has no CLI entry point (the binary
//! drives the GUI). The plan asked for a smoke run via the CLI; this
//! integration test is the equivalent — it drives `valenx_dock::dock()`
//! with the same inputs the GUI would read from `case.toml`, so any
//! regression that breaks the demo case will fail CI.
//!
//! Output is written to a temp dir, not back into the example
//! directory, so the test never mutates the checked-in example.
//!
//! Exhaustiveness is forced to the value from `case.toml` (2), which
//! keeps wall-clock under a minute on a typical dev machine.

use std::path::Path;

#[test]
fn smoke_run_1iep_demo_writes_poses() {
    // 1. Locate the bundled demo case (workspace-root-relative).
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let demo = manifest
        .parent() // -> crates/
        .and_then(|p| p.parent()) // -> workspace root
        .expect("workspace root above crates/")
        .join("examples/dock/1iep_demo");
    assert!(
        demo.join("case.toml").is_file(),
        "demo case.toml missing at {}",
        demo.display()
    );

    // 2. Parse the [bio.vina] block to know what to feed the engine.
    // We mirror the schema VinaInput::from_case_dir uses but only
    // pull the fields valenx_dock::dock needs (no need to take a
    // hard dep on the adapter crate just for TOML reading).
    let case_text = std::fs::read_to_string(demo.join("case.toml")).unwrap();
    let parsed: toml::Value = toml::from_str(&case_text).unwrap();
    let block = &parsed["bio"]["vina"];

    let receptor_path = demo.join(block["receptor"].as_str().unwrap());
    let ligand_path = demo.join(block["ligand"].as_str().unwrap());
    let receptor = std::fs::read_to_string(&receptor_path).unwrap();
    let ligand = std::fs::read_to_string(&ligand_path).unwrap();

    let center: [f64; 3] = block["center"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_float().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let size: [f64; 3] = block["size"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_float().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let exhaustiveness = block["exhaustiveness"].as_integer().unwrap() as u32;
    let num_modes = block["num_modes"].as_integer().unwrap() as u32;
    let energy_range = block["energy_range"].as_float().unwrap();
    assert_eq!(
        block["engine"].as_str().unwrap(),
        "native",
        "demo case should target the native engine"
    );

    // 3. Build the config and run.
    let cfg = valenx_dock::DockConfig {
        center: nalgebra::Vector3::new(center[0], center[1], center[2]),
        size: nalgebra::Vector3::new(size[0], size[1], size[2]),
        exhaustiveness,
        num_modes,
        energy_range,
        seed: 42,
        ..Default::default()
    };

    // Write to temp — never mutate the checked-in example tree.
    let out = std::env::temp_dir().join("valenx_dock_smoke_1iep.pdbqt");
    let _ = std::fs::remove_file(&out);
    let started = std::time::Instant::now();
    let poses =
        valenx_dock::dock(&receptor, &ligand, &cfg, &out, None).expect("native dock failed");
    eprintln!(
        "1iep demo: {} poses in {:.1}s (exhaustiveness={exhaustiveness})",
        poses.len(),
        started.elapsed().as_secs_f64()
    );
    assert!(!poses.is_empty(), "native dock returned no poses");

    // 4. Verify the output PDBQT has the expected shape.
    let written = std::fs::read_to_string(&out).expect("output PDBQT not written");
    assert!(
        written.contains("MODEL"),
        "output PDBQT has no MODEL block:\n{written}"
    );
    assert!(
        written.contains("VINA RESULT:"),
        "output PDBQT has no `VINA RESULT:` score line:\n{written}"
    );
    let _ = std::fs::remove_file(&out);
}

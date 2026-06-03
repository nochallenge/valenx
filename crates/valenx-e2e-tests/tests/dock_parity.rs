//! Parity check: native valenx-dock vs reference AutoDock Vina binary.
//! Skipped unless the `VINA_BIN` environment variable points at a
//! Vina 1.2 executable.

use std::process::Command;

#[test]
fn native_dock_top1_within_2_angstroms_of_reference() {
    let Ok(vina_bin) = std::env::var("VINA_BIN") else {
        eprintln!("skipping: VINA_BIN not set");
        return;
    };

    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/dock/1iep");
    let receptor = std::fs::read_to_string(fixture.join("receptor.pdbqt")).unwrap();
    let ligand = std::fs::read_to_string(fixture.join("ligand.pdbqt")).unwrap();
    let box_toml: toml::Value =
        toml::from_str(&std::fs::read_to_string(fixture.join("box.toml")).unwrap()).unwrap();
    let center: [f64; 3] = box_toml["center"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_float().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let size: [f64; 3] = box_toml["size"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_float().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    // 1. Run reference Vina.
    let tmp = std::env::temp_dir().join("valenx_dock_parity_ref.pdbqt");
    let status = Command::new(&vina_bin)
        .args([
            "--receptor",
            fixture.join("receptor.pdbqt").to_str().unwrap(),
            "--ligand",
            fixture.join("ligand.pdbqt").to_str().unwrap(),
            "--out",
            tmp.to_str().unwrap(),
            "--center_x",
            &center[0].to_string(),
            "--center_y",
            &center[1].to_string(),
            "--center_z",
            &center[2].to_string(),
            "--size_x",
            &size[0].to_string(),
            "--size_y",
            &size[1].to_string(),
            "--size_z",
            &size[2].to_string(),
            "--exhaustiveness",
            "8",
            "--num_modes",
            "9",
            "--seed",
            "42",
        ])
        .status()
        .expect("vina launch");
    assert!(status.success(), "reference vina returned non-zero");
    let reference_text = std::fs::read_to_string(&tmp).unwrap();
    let _ = std::fs::remove_file(&tmp);

    // 2. Run native.
    let cfg = valenx_dock::DockConfig {
        center: nalgebra::Vector3::new(center[0], center[1], center[2]),
        size: nalgebra::Vector3::new(size[0], size[1], size[2]),
        exhaustiveness: 8,
        num_modes: 9,
        seed: 42,
        ..Default::default()
    };
    let tmp_native = std::env::temp_dir().join("valenx_dock_parity_native.pdbqt");
    let poses = valenx_dock::dock(&receptor, &ligand, &cfg, &tmp_native, None).unwrap();
    assert!(!poses.is_empty(), "native dock returned no poses");
    let _ = std::fs::remove_file(&tmp_native);

    // 3. Compare top-1 RMSD. Parse reference's MODEL 1 atoms and
    // native's top pose atoms; compute heavy-atom RMSD.
    let ref_top = first_model_atoms(&reference_text);
    let lig = valenx_dock::ligand::Ligand::from_pdbqt(&ligand).unwrap();
    let native_top_positions = lig.apply_pose(&poses[0].0);
    // Reject silently-truncated comparisons: zip() would happily
    // run over min(ref, native, lig) atoms and report a misleading
    // RMSD if anything was lost in PDBQT round-trip.
    assert_eq!(
        ref_top.len(),
        native_top_positions.len(),
        "reference and native atom counts differ: ref={} native={}",
        ref_top.len(),
        native_top_positions.len(),
    );
    assert_eq!(
        native_top_positions.len(),
        lig.atoms.len(),
        "native positions vs ligand atoms length mismatch",
    );
    let mut sum_sq = 0.0;
    let mut n = 0;
    for ((rp, lp), atom) in ref_top
        .iter()
        .zip(native_top_positions.iter())
        .zip(lig.atoms.iter())
    {
        use valenx_dock::atom_type::Ad4AtomType;
        if matches!(atom.ad4_type, Ad4AtomType::H | Ad4AtomType::HD) {
            continue;
        }
        sum_sq += (rp - lp).norm_squared();
        n += 1;
    }
    let rmsd = (sum_sq / n as f64).sqrt();
    eprintln!("top-1 RMSD reference vs native = {rmsd:.3} Å");
    assert!(rmsd <= 2.0, "top-1 RMSD {rmsd:.3} Å exceeds 2.0 Å target");
}

fn first_model_atoms(text: &str) -> Vec<nalgebra::Vector3<f64>> {
    use valenx_bio::format::pdbqt::{parse, PdbqtRecord};
    let block = text.split("MODEL ").nth(1).unwrap_or("");
    let recs = parse(block).unwrap_or_default();
    recs.into_iter()
        .filter_map(|r| {
            if let PdbqtRecord::Atom(a) = r {
                Some(a.position)
            } else {
                None
            }
        })
        .collect()
}

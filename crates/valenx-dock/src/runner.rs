//! Top-level docking entry point.

use std::path::Path;

use crate::cluster::cluster_poses;
use crate::config::DockConfig;
use crate::error::DockError;
use crate::grid::GridBundle;
use crate::ligand::Ligand;
use crate::pose::Pose;
use crate::receptor::Receptor;
use crate::search::mc::exhaustiveness_search;

/// Progress callback signature. `(fraction, message)`; invoked only
/// from the calling thread between phases, so no Send/Sync bound is
/// required. Keep it cheap — it runs on the caller's thread while
/// the search workers are paused.
pub type ProgressFn<'a> = &'a (dyn Fn(f64, &str));

/// Run a full docking job: read receptor + ligand PDBQT, search,
/// cluster, write output PDBQT.
///
/// Returns the clustered, ranked poses (best first).
pub fn dock(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    config: &DockConfig,
    output_path: &Path,
    progress: Option<ProgressFn>,
) -> Result<Vec<(Pose, f64)>, DockError> {
    config.validate()?;
    let receptor = Receptor::from_pdbqt(receptor_pdbqt)?;
    let ligand = Ligand::from_pdbqt(ligand_pdbqt)?;
    if let Some(p) = progress {
        p(0.10, "building receptor grids");
    }
    let grids = GridBundle::build(
        &receptor,
        &ligand,
        config.grid_origin(),
        config.grid_spacing,
        config.grid_dims(),
    );
    if let Some(p) = progress {
        p(0.25, "performing search");
    }
    let raw = exhaustiveness_search(
        &ligand,
        &grids,
        config.center,
        config.size,
        config.exhaustiveness as usize,
        50, // inner ILS iterations per chain
        config.seed,
    );
    if let Some(p) = progress {
        p(0.85, "clustering poses");
    }
    let mut clustered = cluster_poses(&ligand, &raw, 2.0);
    // Cap by num_modes + energy_range.
    if let Some(best) = clustered.first().map(|p| p.1) {
        clustered.retain(|p| p.1 - best <= config.energy_range);
    }
    clustered.truncate(config.num_modes as usize);
    if let Some(p) = progress {
        p(0.95, "writing output");
    }
    write_output(&ligand, &clustered, output_path, ligand_pdbqt)?;
    if let Some(p) = progress {
        p(1.0, "done");
    }
    Ok(clustered)
}

/// `dock()` variant that emits structured events to a callback. The
/// callback runs on the calling thread (between phases); make it cheap.
pub fn dock_with_events<F: FnMut(&crate::events::DockEvent)>(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    config: &DockConfig,
    output_path: &Path,
    mut on_event: F,
) -> Result<Vec<(Pose, f64)>, DockError> {
    use crate::atom_type::Ad4AtomType;
    use crate::events::DockEvent;
    config.validate().map_err(|e| {
        on_event(&DockEvent::Error {
            code: e.code(),
            category: format!("{:?}", e.category()),
            message: format!("{e}"),
        });
        e
    })?;
    let receptor = Receptor::from_pdbqt(receptor_pdbqt).map_err(|e| {
        on_event(&DockEvent::Error {
            code: e.code(),
            category: format!("{:?}", e.category()),
            message: format!("{e}"),
        });
        e
    })?;
    let ligand = Ligand::from_pdbqt(ligand_pdbqt).map_err(|e| {
        on_event(&DockEvent::Error {
            code: e.code(),
            category: format!("{:?}", e.category()),
            message: format!("{e}"),
        });
        e
    })?;
    let heavy = ligand
        .atoms
        .iter()
        .filter(|a| !matches!(a.ad4_type, Ad4AtomType::H | Ad4AtomType::HD))
        .count();
    on_event(&DockEvent::Started {
        receptor: "<inline>".into(),
        ligand: "<inline>".into(),
        ligand_heavy_atoms: heavy,
        n_torsions: ligand.n_torsions(),
    });
    let grids = GridBundle::build(
        &receptor,
        &ligand,
        config.grid_origin(),
        config.grid_spacing,
        config.grid_dims(),
    );
    let total_voxels: usize = grids
        .grids
        .values()
        .map(|g| g.dims.0 * g.dims.1 * g.dims.2)
        .sum();
    on_event(&DockEvent::GridBuilt {
        n_grids: grids.grids.len(),
        total_voxels,
    });
    let t0 = std::time::Instant::now();
    let raw = crate::search::mc::exhaustiveness_search(
        &ligand,
        &grids,
        config.center,
        config.size,
        config.exhaustiveness as usize,
        50,
        config.seed,
    );
    on_event(&DockEvent::SearchProgress {
        fraction: 1.0,
        best_score: raw.first().map(|p| p.1).unwrap_or(0.0),
    });
    let mut clustered = crate::cluster::cluster_poses(&ligand, &raw, 2.0);
    if let Some(best) = clustered.first().map(|p| p.1) {
        clustered.retain(|p| p.1 - best <= config.energy_range);
    }
    clustered.truncate(config.num_modes as usize);
    let top_score = clustered.first().map(|p| p.1).unwrap_or(0.0);
    for (i, (pose, score)) in clustered.iter().enumerate() {
        let r = if i == 0 {
            0.0
        } else {
            crate::cluster::rmsd(&ligand, &clustered[0].0, pose)
        };
        on_event(&DockEvent::PoseFound {
            rank: i + 1,
            score: *score,
            rmsd_to_top: r,
        });
        let _ = top_score; // unused but useful in future refactors
        let _ = pose;
    }
    write_output(&ligand, &clustered, output_path, ligand_pdbqt)?;
    on_event(&DockEvent::Complete {
        n_poses: clustered.len(),
        wall_seconds: t0.elapsed().as_secs_f64(),
    });
    Ok(clustered)
}

pub(crate) fn write_output(
    ligand: &Ligand,
    poses: &[(Pose, f64)],
    output_path: &Path,
    ligand_pdbqt: &str,
) -> Result<(), DockError> {
    use valenx_bio::format::pdbqt::{parse, write_pose_ensemble, PdbqtRecord};
    // Reparse to keep ATOM records (with names, charges, ad4 strings)
    // in their original PDBQT order so the writer round-trips cleanly.
    let records = parse(ligand_pdbqt)?;
    let atoms: Vec<_> = records
        .into_iter()
        .filter_map(|r| {
            if let PdbqtRecord::Atom(a) = r {
                Some(a)
            } else {
                None
            }
        })
        .collect();
    let pose_atoms: Vec<(Vec<nalgebra::Vector3<f64>>, f64)> = poses
        .iter()
        .map(|(p, s)| (ligand.apply_pose(p), *s))
        .collect();
    let body = write_pose_ensemble(&atoms, &pose_atoms);
    atomic_write_pdbqt(output_path, &body)
}

/// Round-26 H1: write `body` to `output_path` via the atomic-rename
/// pattern so concurrent multi-MB writes never interleave bytes at
/// the destination. Factored out of [`write_output`] so the H1 test
/// can exercise the write logic directly with a synthetic multi-MB
/// body — the docking pipeline itself produces deterministic
/// output for the same input, making it awkward to construct "two
/// distinct large ensembles" via real `dock()` calls.
///
/// ## Round-27 STRUCTURAL consolidation
///
/// Now a thin wrapper around
/// [`valenx_core::io_caps::atomic_write_bytes`]. The canonical
/// helper provides all the crash-safety invariants this site used
/// to inline (unique `<pid>.<counter>` sidecar, `O_NOFOLLOW` /
/// `FILE_FLAG_OPEN_REPARSE_POINT`, fsync-before-rename, parent-dir
/// fsync on Unix) — and via consolidation now also picks up the
/// parent-fsync that was only landed in `valenx-app::state_paths`
/// in R26 M3 (closing M1's "dock parent-fsync sister gap"). All 4
/// atomic-write sites (state_paths, dock, crash-reporter,
/// render-bridge) delegate to the same code now, so any future
/// crash-safety fix lands once instead of being back-ported across
/// 4 copies.
///
/// The wrapper translates the canonical helper's `io::Error` into
/// `DockError::Io` so the dock-runner's error type contract stays
/// stable.
pub(crate) fn atomic_write_pdbqt(output_path: &Path, body: &str) -> Result<(), DockError> {
    valenx_core::io_caps::atomic_write_bytes(output_path, body.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dock_end_to_end_writes_output_with_at_least_one_model() {
        let receptor =
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
";
        let ligand = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        let tmp = std::env::temp_dir().join("valenx_dock_test_out.pdbqt");
        let cfg = DockConfig {
            center: nalgebra::Vector3::new(0.75, 0.0, 0.0),
            size: nalgebra::Vector3::new(6.0, 6.0, 6.0),
            exhaustiveness: 2,
            num_modes: 3,
            seed: 7,
            ..DockConfig::default()
        };
        let poses = dock(receptor, ligand, &cfg, &tmp, None).unwrap();
        assert!(!poses.is_empty(), "got no poses");
        let body = std::fs::read_to_string(&tmp).unwrap();
        assert!(body.contains("MODEL 1"));
        assert!(body.contains("VINA RESULT"));
        let _ = std::fs::remove_file(&tmp);
    }

    /// RED→GREEN (round-24 M4): re-running `dock` against the same
    /// output path must overwrite the previous file (pre-fix `fs::write`
    /// silently overwrote; post-fix the new `remove_file + create_new`
    /// pattern preserves the same overwrite semantic for legitimate
    /// repeated runs while refusing to follow a TOCTOU symlink). This
    /// is the regression-anchor for the overwrite contract — TOCTOU
    /// itself can't be reproduced reliably in a single-process test
    /// (Unix `unlink` is atomic so no window exists between unlink
    /// and create_new for the same thread).
    #[test]
    fn dock_overwrites_existing_output_post_m4() {
        let receptor =
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
";
        let ligand = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        let tmp = std::env::temp_dir().join("valenx_dock_m4_overwrite.pdbqt");
        // Pre-seed the output path with stale junk.
        std::fs::write(&tmp, b"STALE CONTENTS").unwrap();
        let cfg = DockConfig {
            center: nalgebra::Vector3::new(0.75, 0.0, 0.0),
            size: nalgebra::Vector3::new(6.0, 6.0, 6.0),
            exhaustiveness: 2,
            num_modes: 3,
            seed: 7,
            ..DockConfig::default()
        };
        let _ = dock(receptor, ligand, &cfg, &tmp, None).unwrap();
        let body = std::fs::read_to_string(&tmp).unwrap();
        // The stale string must be gone — overwrite preserved.
        assert!(!body.contains("STALE CONTENTS"));
        assert!(body.contains("MODEL 1"));
        let _ = std::fs::remove_file(&tmp);
    }

    /// RED→GREEN (round-25 H2): concurrent `dock` runs against the
    /// same output path must all succeed without an `AlreadyExists`
    /// error. Pre-fix (round-24 M4) used `remove_file + create_new`,
    /// which had a DoS race: when thread A and thread B both unlinked
    /// then create_new'd, B's create_new could land between A's
    /// unlink and A's create_new — B would succeed and A would error
    /// `AlreadyExists` (or vice versa). With the round-25 fix
    /// (`truncate(true)`) both opens proceed regardless of timing;
    /// the kernel serialises the writes. We assert ALL of N concurrent
    /// runs return Ok.
    #[test]
    fn dock_concurrent_writes_no_dos_race_round25_h2() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        let receptor =
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
";
        let ligand = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        let tmp = std::env::temp_dir().join(format!(
            "valenx_dock_h2_concurrent_{}.pdbqt",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // 4 concurrent writers is enough to surface the race
        // reproducibly without making the test slow — each dock run
        // takes ~50ms on a debug build.
        const N: usize = 4;
        let barrier = Arc::new(Barrier::new(N));
        let mut handles = Vec::with_capacity(N);
        for _ in 0..N {
            let receptor = receptor.to_string();
            let ligand = ligand.to_string();
            let tmp = tmp.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || {
                barrier.wait();
                let cfg = DockConfig {
                    center: nalgebra::Vector3::new(0.75, 0.0, 0.0),
                    size: nalgebra::Vector3::new(6.0, 6.0, 6.0),
                    exhaustiveness: 2,
                    num_modes: 3,
                    seed: 7,
                    ..DockConfig::default()
                };
                dock(&receptor, &ligand, &cfg, &tmp, None)
            }));
        }
        let mut ok = 0;
        let mut errs: Vec<String> = Vec::new();
        for h in handles {
            match h.join().unwrap() {
                Ok(_) => ok += 1,
                Err(e) => errs.push(format!("{e}")),
            }
        }
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(
            ok, N,
            "all {N} concurrent dock runs must succeed (got {ok}), errs: {errs:?}",
        );
    }

    /// RED→GREEN (round-26 H1): two concurrent writers each writing
    /// a distinct multi-MB body must end with one writer's content
    /// at the output path — NOT an interleaved mix of both.
    ///
    /// Pre-fix (round-25 H2) used
    /// `write(true).create(true).truncate(true) + O_NOFOLLOW` on
    /// the output path itself; two threads sharing that file
    /// descriptor saw their `write_all` syscalls interleaved by the
    /// kernel for any payload > PIPE_BUF (4 KiB on Linux). 5 MiB
    /// ensures the interleaving is observable. Post-fix each writer
    /// owns its own `<output>.tmp.<pid>.<counter>` sidecar (no
    /// shared handle), then `fs::rename` atomically promotes ONE
    /// sidecar over the target — the final file contains exactly
    /// one writer's content.
    #[test]
    fn atomic_write_pdbqt_concurrent_multi_mb_no_interleaving_round26_h1() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        // 5 MiB blocks of distinct fill bytes ('A' / 'B') so any
        // interleaved write that lands in the final file is
        // detectable as a non-contiguous run.
        const BODY_BYTES: usize = 5 * 1024 * 1024;
        let body_a = "A".repeat(BODY_BYTES);
        let body_b = "B".repeat(BODY_BYTES);
        let tmp = std::env::temp_dir().join(format!(
            "valenx_dock_h1_interleave_{}.pdbqt",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // Two threads, joined start via a Barrier so the writes
        // overlap in time.
        let barrier = Arc::new(Barrier::new(2));
        let body_a_clone = body_a.clone();
        let body_b_clone = body_b.clone();
        let tmp_a = tmp.clone();
        let tmp_b = tmp.clone();
        let bar_a = Arc::clone(&barrier);
        let bar_b = Arc::clone(&barrier);
        let ha = thread::spawn(move || {
            bar_a.wait();
            atomic_write_pdbqt(&tmp_a, &body_a_clone)
        });
        let hb = thread::spawn(move || {
            bar_b.wait();
            atomic_write_pdbqt(&tmp_b, &body_b_clone)
        });
        ha.join().unwrap().expect("writer A");
        hb.join().unwrap().expect("writer B");
        let final_body = std::fs::read_to_string(&tmp).expect("read final");
        // Best part of the contract: the final body is EXACTLY one
        // writer's content, never a mix. We check by comparing
        // against both candidates and asserting at least one matches.
        let matches_a = final_body == body_a;
        let matches_b = final_body == body_b;
        let _ = std::fs::remove_file(&tmp);
        // Also clean up any leaked sidecars (best-effort).
        if let Some(parent) = tmp.parent() {
            if let Ok(rd) = std::fs::read_dir(parent) {
                for entry in rd.flatten() {
                    let name = entry.file_name();
                    if let Some(s) = name.to_str() {
                        if s.starts_with(tmp.file_name().and_then(|n| n.to_str()).unwrap_or(""))
                            && s.contains(".tmp.")
                        {
                            let _ = std::fs::remove_file(entry.path());
                        }
                    }
                }
            }
        }
        assert!(
            matches_a || matches_b,
            "final body must equal exactly one writer's content — \
             interleaving detected (len={} expected={BODY_BYTES})",
            final_body.len(),
        );
    }
}

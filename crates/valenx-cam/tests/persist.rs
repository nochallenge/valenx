//! Persistence round-trip tests for [`valenx_cam::persist::CamFile`].

use nalgebra::Vector3;
use valenx_cam::operation::{Operation, PocketParams, PocketStrategy, ProfileParams};
use valenx_cam::persist::CamFile;
use valenx_cam::stock::Stock;
use valenx_cam::tool::{Tool, ToolKind};

#[test]
fn empty_round_trips() {
    let f = CamFile::new();
    let ron = f.to_ron().unwrap();
    let parsed = CamFile::from_ron(&ron).unwrap();
    assert_eq!(parsed.version, 2); // Phase 17 bumped to 2 (added fixture + setups)
    assert!(parsed.tools.is_empty());
    assert!(parsed.operations.is_empty());
}

#[test]
fn populated_round_trips() {
    let mut f = CamFile::new();
    f.stock = Stock::new(
        Vector3::new(-10.0, -10.0, 0.0),
        Vector3::new(50.0, 50.0, 20.0),
        "6061",
    )
    .unwrap();
    f.tools
        .push(Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap());
    f.tools
        .push(Tool::new(2, "Drill3", ToolKind::Drill, 3.0, 30.0, 2, "HSS").unwrap());
    f.operations
        .push(Operation::Profile(ProfileParams::default()));
    f.operations.push(Operation::Pocket(PocketParams {
        strategy: PocketStrategy::Spiral,
        step_over: 1.0,
        ..Default::default()
    }));
    let ron = f.to_ron().unwrap();
    let parsed = CamFile::from_ron(&ron).unwrap();
    assert_eq!(parsed.tools.len(), 2);
    assert_eq!(parsed.operations.len(), 2);
    assert_eq!(parsed.stock.material, "6061");
    assert!((parsed.stock.size.x - 50.0).abs() < 1e-9);
    if let Operation::Pocket(p) = &parsed.operations[1] {
        assert_eq!(p.strategy, PocketStrategy::Spiral);
        assert!((p.step_over - 1.0).abs() < 1e-9);
    } else {
        panic!("operation 1 should be Pocket");
    }
}

#[test]
fn file_round_trip() {
    let mut f = CamFile::new();
    f.tools
        .push(Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap());
    let tmp = std::env::temp_dir().join("valenx_cam_persist_test.ron");
    f.write_to(&tmp).unwrap();
    let parsed = CamFile::read_from(&tmp).unwrap();
    assert_eq!(parsed.tools.len(), 1);
    assert_eq!(parsed.tools[0].name, "EM6");
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn malformed_ron_errors() {
    let bad = CamFile::from_ron("not a ron file");
    assert!(bad.is_err());
    assert_eq!(bad.unwrap_err().code(), "cam.ron");
}

/// Round-28 H2 RED→GREEN (Unix) — pre-fix `CamFile::write_to` went
/// through `std::fs::write`, which follows leaf symlinks. An
/// attacker who pre-created the target as a symlink to a sentinel
/// file would have the CAM file's RON written through the symlink,
/// clobbering the sentinel. Post-fix the canonical
/// `valenx_core::io_caps::atomic_write_str` opens its sidecar
/// O_NOFOLLOW and renames it over the target — the original target
/// (a symlink) is replaced by a regular file, but the symlink's
/// target (the sentinel) is never written through.
#[cfg(unix)]
#[test]
fn write_to_does_not_clobber_symlink_target_round28_h2() {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let dir = std::env::temp_dir().join(format!("valenx-cam-r28-h2-{nanos}"));
    std::fs::create_dir_all(&dir).expect("create tmp dir");
    let sentinel = dir.join("sentinel.txt");
    std::fs::write(&sentinel, b"untouched\n").expect("write sentinel");
    let dst = dir.join("session.ron");
    std::os::unix::fs::symlink(&sentinel, &dst).expect("create symlink");

    let f = CamFile::new();
    f.write_to(&dst).expect("write_to over symlink");

    // Sentinel must still hold its original bytes — the canonical
    // helper never wrote through the symlink.
    let after = std::fs::read(&sentinel).expect("read sentinel");
    assert_eq!(after, b"untouched\n", "sentinel was clobbered");
    // The symlink itself has been replaced by a regular file
    // holding the RON.
    let dst_meta = std::fs::symlink_metadata(&dst).expect("metadata");
    assert!(!dst_meta.file_type().is_symlink(), "symlink survived");

    let _ = std::fs::remove_dir_all(&dir);
}

/// Round-28 H2 RED→GREEN — post-fix invariant: concurrent
/// `write_to` calls all return Ok, no `AlreadyExists` surfaces, and
/// the final file parses as valid `CamFile`. The structural anchor
/// is that every writer gets its own
/// `<basename>.tmp.<pid>.<counter>` sidecar; the rename then
/// atomically publishes a complete file from whichever thread wins.
/// Pre-fix `std::fs::write` happens to not raise `AlreadyExists`
/// because it opens with truncate (not `create_new`), so this test
/// passes both pre-fix and post-fix on Windows — the real bug is
/// the truncate-then-write window which can produce torn content
/// AND the leaf-symlink follow (covered by the Unix test above).
/// This concurrency test guards against a future regression that
/// re-introduces `create_new`-based sidecars without per-writer
/// namespacing.
#[test]
fn write_to_concurrent_calls_yield_valid_ron_round28_h2() {
    use std::sync::{Arc, Barrier};
    use valenx_cam::operation::{Operation, ProfileParams};
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let dir = std::env::temp_dir().join(format!("valenx-cam-r28-h2-conc-{nanos}"));
    std::fs::create_dir_all(&dir).expect("create tmp dir");
    let dst = Arc::new(dir.join("session.ron"));
    let n = 8usize;
    let barrier = Arc::new(Barrier::new(n));
    let mut handles = Vec::with_capacity(n);
    for i in 0..n {
        let dst = Arc::clone(&dst);
        let b = Arc::clone(&barrier);
        handles.push(std::thread::spawn(move || {
            let mut f = CamFile::new();
            for _ in 0..i {
                f.operations
                    .push(Operation::Profile(ProfileParams::default()));
            }
            b.wait();
            f.write_to(&dst)
        }));
    }
    for h in handles {
        let res = h.join().expect("thread did not panic");
        res.expect("concurrent write_to must not error");
    }
    let parsed = CamFile::read_from(&dst).expect("final file parses without torn RON");
    assert_eq!(parsed.version, 2);
    let _ = std::fs::remove_dir_all(&dir);
}

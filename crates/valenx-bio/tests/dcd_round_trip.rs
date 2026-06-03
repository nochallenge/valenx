//! End-to-end DCD reader tests against hand-synthesised byte
//! streams. The unit tests in `src/format/dcd.rs` cover the same
//! shapes; the duplicates here pin the integration-test path so a
//! breaking API change in the public `read()` surface trips the
//! integration suite as well.

use nalgebra::Vector3;
use valenx_bio::format::dcd::{self, DcdError};

/// Wrap a record body with the leading + trailing 4-byte LE length
/// prefixes Fortran unformatted I/O expects.
fn wrap(body: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(body.len() + 8);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(body);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out
}

fn build_header(nframes: u32) -> Vec<u8> {
    // 84-byte body: 4 magic + 9*i32 + f32 + i32 + 8*i32 + i32.
    let mut h = Vec::with_capacity(84);
    h.extend_from_slice(b"CORD");
    h.extend_from_slice(&(nframes as i32).to_le_bytes());
    for _ in 0..3 {
        h.extend_from_slice(&0i32.to_le_bytes());
    }
    for _ in 0..4 {
        h.extend_from_slice(&0i32.to_le_bytes());
    }
    h.extend_from_slice(&0i32.to_le_bytes());
    h.extend_from_slice(&0.001f32.to_le_bytes());
    h.extend_from_slice(&0i32.to_le_bytes());
    for _ in 0..8 {
        h.extend_from_slice(&0i32.to_le_bytes());
    }
    h.extend_from_slice(&24i32.to_le_bytes());
    debug_assert_eq!(h.len(), 84);
    h
}

fn build_titles() -> Vec<u8> {
    let mut t = Vec::with_capacity(4 + 80);
    t.extend_from_slice(&1i32.to_le_bytes());
    let title = b"integration-test DCD";
    let mut padded = [b' '; 80];
    padded[..title.len()].copy_from_slice(title);
    t.extend_from_slice(&padded);
    t
}

fn pack_f32s(values: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

#[test]
fn read_synthesised_two_frame_dcd() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&wrap(&build_header(2)));
    bytes.extend_from_slice(&wrap(&build_titles()));
    bytes.extend_from_slice(&wrap(&3i32.to_le_bytes()));
    // Frame 0 — three atoms along the x-axis.
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.0f32, 1.0, 2.0])));
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.0f32, 0.0, 0.0])));
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.0f32, 0.0, 0.0])));
    // Frame 1 — slightly displaced.
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.5f32, 1.5, 2.5])));
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.1f32, 0.1, 0.1])));
    bytes.extend_from_slice(&wrap(&pack_f32s(&[0.2f32, 0.2, 0.2])));

    let traj = dcd::read(&bytes, "integration").expect("parses");
    assert_eq!(traj.id, "integration");
    assert_eq!(traj.frame_count(), 2);
    assert_eq!(traj.atom_count(), Some(3));
    let f0 = traj.frame(0).expect("frame 0");
    assert_eq!(f0[0], Vector3::new(0.0, 0.0, 0.0));
    assert_eq!(f0[2], Vector3::new(2.0, 0.0, 0.0));
    let f1 = traj.frame(1).expect("frame 1");
    let eps = 1e-5;
    assert!((f1[0].x - 0.5).abs() < eps);
    assert!((f1[2].z - 0.2).abs() < eps);
}

#[test]
fn read_rejects_non_cord_magic() {
    let mut header = build_header(0);
    header[0..4].copy_from_slice(b"FAKE");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&wrap(&header));
    bytes.extend_from_slice(&wrap(&build_titles()));
    bytes.extend_from_slice(&wrap(&0i32.to_le_bytes()));

    let err = dcd::read(&bytes, "fake").expect_err("magic check");
    match err {
        DcdError::Magic { found, expected } => {
            assert_eq!(&found, b"FAKE");
            assert_eq!(expected, "CORD");
        }
        other => panic!("expected Magic, got {other:?}"),
    }
}

#[test]
fn read_handles_zero_frames() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&wrap(&build_header(0)));
    bytes.extend_from_slice(&wrap(&build_titles()));
    bytes.extend_from_slice(&wrap(&12i32.to_le_bytes()));

    let traj = dcd::read(&bytes, "zero").expect("zero-frame DCD parses");
    assert_eq!(traj.frame_count(), 0);
    assert!(traj.frames.is_empty());
}

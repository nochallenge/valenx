//! Trajectory container and writers — **roadmap feature 6**.
//!
//! A *trajectory* is the time series of coordinate frames an MD run
//! produces. This module provides:
//!
//! - [`Trajectory`] — an in-memory container of frames (each frame a
//!   `Vec<Vector3<f64>>`), all with the same atom count.
//! - A **binary** writer / reader, [`write_binary`] / [`read_binary`],
//!   in a compact DCD-class layout: a small magic-tagged header
//!   (atom count, frame count, time step) followed by the raw `f32`
//!   coordinates of every frame. Like CHARMM/NAMD DCD it stores
//!   single-precision coordinates separated by component
//!   (`x[0..N] y[0..N] z[0..N]`); it is a self-describing native
//!   format, not byte-compatible with the historical Fortran DCD
//!   record framing.
//! - A **framed-text** writer / reader, [`write_text`] / [`read_text`]
//!   — a human-readable `FRAME` / coordinate-line format for small
//!   trajectories and debugging.
//!
//! All integers are written little-endian by hand, so the format is
//! stable and dependency-free.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::system::System;

/// An in-memory trajectory: a sequence of coordinate frames.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trajectory {
    /// Number of atoms in every frame.
    natoms: usize,
    /// Time between consecutive frames (ps).
    dt: f64,
    /// The frames, in order; each is `natoms` positions (nm).
    frames: Vec<Vec<Vector3<f64>>>,
}

impl Trajectory {
    /// An empty trajectory for `natoms`-atom frames spaced `dt` ps
    /// apart.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `natoms` is zero or `dt` is not finite
    /// and positive.
    pub fn new(natoms: usize, dt: f64) -> Result<Self> {
        if natoms == 0 {
            return Err(MdError::invalid("natoms", "must be at least 1"));
        }
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        Ok(Trajectory {
            natoms,
            dt,
            frames: Vec::new(),
        })
    }

    /// Number of atoms per frame.
    pub fn natoms(&self) -> usize {
        self.natoms
    }

    /// Time between frames (ps).
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Number of stored frames.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether the trajectory has no frames.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// The frames as a slice.
    pub fn frames(&self) -> &[Vec<Vector3<f64>>] {
        &self.frames
    }

    /// Borrows frame `index`, if it exists.
    pub fn frame(&self, index: usize) -> Option<&[Vector3<f64>]> {
        self.frames.get(index).map(|f| f.as_slice())
    }

    /// Appends a frame.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the frame's atom count does
    /// not match the trajectory.
    pub fn push_frame(&mut self, frame: Vec<Vector3<f64>>) -> Result<()> {
        if frame.len() != self.natoms {
            return Err(MdError::dimension(format!(
                "frame has {} atoms, trajectory expects {}",
                frame.len(),
                self.natoms
            )));
        }
        self.frames.push(frame);
        Ok(())
    }

    /// Appends the current positions of a [`System`] as a frame.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] on an atom-count mismatch.
    pub fn push_system(&mut self, system: &System) -> Result<()> {
        self.push_frame(system.positions.clone())
    }
}

/// Magic bytes at the start of a binary trajectory file.
const BINARY_MAGIC: &[u8; 8] = b"VLXMDTRJ";

/// Serialises a trajectory to the compact binary (DCD-class) format.
pub fn write_binary(traj: &Trajectory) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(BINARY_MAGIC);
    // Header: version, natoms, nframes (u32 LE), dt (f64 LE).
    out.extend_from_slice(&1u32.to_le_bytes());
    out.extend_from_slice(&(traj.natoms as u32).to_le_bytes());
    out.extend_from_slice(&(traj.frames.len() as u32).to_le_bytes());
    out.extend_from_slice(&traj.dt.to_le_bytes());
    // Each frame: x[], then y[], then z[] as f32 LE (the DCD layout).
    for frame in &traj.frames {
        for axis in 0..3 {
            for p in frame {
                let v = p[axis] as f32;
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
    }
    out
}

/// Reads a binary (DCD-class) trajectory.
///
/// # Errors
/// [`MdError::Parse`] on a bad magic / version, a truncated header, or
/// a body shorter than the header promises.
pub fn read_binary(bytes: &[u8]) -> Result<Trajectory> {
    if bytes.len() < 8 + 4 + 4 + 4 + 8 {
        return Err(MdError::parse("dcd", "file is shorter than the header"));
    }
    if &bytes[0..8] != BINARY_MAGIC {
        return Err(MdError::parse("dcd", "bad magic bytes"));
    }
    let u32_at = |off: usize| -> u32 {
        u32::from_le_bytes([
            bytes[off],
            bytes[off + 1],
            bytes[off + 2],
            bytes[off + 3],
        ])
    };
    let version = u32_at(8);
    if version != 1 {
        return Err(MdError::parse(
            "dcd",
            format!("unsupported version {version}"),
        ));
    }
    let natoms = u32_at(12) as usize;
    let nframes = u32_at(16) as usize;
    let mut dt_bytes = [0u8; 8];
    dt_bytes.copy_from_slice(&bytes[20..28]);
    let dt = f64::from_le_bytes(dt_bytes);

    if natoms == 0 {
        return Err(MdError::parse("dcd", "header declares zero atoms"));
    }
    let mut traj = Trajectory::new(natoms, if dt > 0.0 { dt } else { 1.0 })?;

    let floats_per_frame = natoms * 3;
    let frame_bytes = floats_per_frame * 4;
    let expected = 28 + nframes * frame_bytes;
    if bytes.len() < expected {
        return Err(MdError::parse(
            "dcd",
            format!("body is {} bytes, expected {expected}", bytes.len()),
        ));
    }
    let mut offset = 28;
    for _ in 0..nframes {
        let mut xs = vec![0f32; natoms];
        let mut ys = vec![0f32; natoms];
        let mut zs = vec![0f32; natoms];
        for buf in [&mut xs, &mut ys, &mut zs] {
            for slot in buf.iter_mut() {
                let mut b = [0u8; 4];
                b.copy_from_slice(&bytes[offset..offset + 4]);
                *slot = f32::from_le_bytes(b);
                offset += 4;
            }
        }
        let frame: Vec<Vector3<f64>> = (0..natoms)
            .map(|i| Vector3::new(xs[i] as f64, ys[i] as f64, zs[i] as f64))
            .collect();
        traj.push_frame(frame)?;
    }
    Ok(traj)
}

/// Serialises a trajectory to the human-readable framed-text format.
pub fn write_text(traj: &Trajectory) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "# valenx-md trajectory  natoms={} nframes={} dt={}\n",
        traj.natoms,
        traj.frames.len(),
        traj.dt
    ));
    for (f, frame) in traj.frames.iter().enumerate() {
        out.push_str(&format!("FRAME {f}\n"));
        for p in frame {
            out.push_str(&format!("{:14.8} {:14.8} {:14.8}\n", p.x, p.y, p.z));
        }
    }
    out
}

/// Reads a framed-text trajectory.
///
/// # Errors
/// [`MdError::Parse`] on a missing header, a malformed `FRAME` marker,
/// a bad coordinate, or a frame of the wrong length.
pub fn read_text(text: &str) -> Result<Trajectory> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| MdError::parse("trajectory-text", "missing header line"))?;
    // Parse `natoms=` and `dt=` out of the header.
    let natoms = header_value(header, "natoms=")
        .and_then(|s| s.parse::<usize>().ok())
        .ok_or_else(|| MdError::parse("trajectory-text", "header missing natoms="))?;
    let dt = header_value(header, "dt=")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0);
    if natoms == 0 {
        return Err(MdError::parse("trajectory-text", "natoms is zero"));
    }
    let mut traj = Trajectory::new(natoms, if dt > 0.0 { dt } else { 1.0 })?;

    let mut current: Vec<Vector3<f64>> = Vec::new();
    let mut in_frame = false;
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("FRAME") {
            // Close the previous frame.
            if in_frame {
                finish_frame(&mut traj, &mut current, natoms)?;
            }
            let _ = rest; // the frame index is informational
            in_frame = true;
            current = Vec::with_capacity(natoms);
        } else {
            let mut fields = trimmed.split_whitespace();
            let parse = |s: Option<&str>| -> Result<f64> {
                s.ok_or_else(|| {
                    MdError::parse("trajectory-text", "short coordinate line")
                })?
                .parse::<f64>()
                .map_err(|_| MdError::parse("trajectory-text", "bad coordinate"))
            };
            let x = parse(fields.next())?;
            let y = parse(fields.next())?;
            let z = parse(fields.next())?;
            current.push(Vector3::new(x, y, z));
        }
    }
    if in_frame {
        finish_frame(&mut traj, &mut current, natoms)?;
    }
    Ok(traj)
}

/// Pushes a completed frame, checking its length.
fn finish_frame(
    traj: &mut Trajectory,
    current: &mut Vec<Vector3<f64>>,
    natoms: usize,
) -> Result<()> {
    if current.len() != natoms {
        return Err(MdError::parse(
            "trajectory-text",
            format!("a frame has {} atoms, expected {natoms}", current.len()),
        ));
    }
    traj.push_frame(std::mem::take(current))?;
    Ok(())
}

/// Extracts the token following `key` in a `key=value`-style header.
fn header_value<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    header
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix(key))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_trajectory() -> Trajectory {
        let mut traj = Trajectory::new(3, 0.002).unwrap();
        for f in 0..5 {
            let frame: Vec<Vector3<f64>> = (0..3)
                .map(|a| {
                    Vector3::new(
                        a as f64 + 0.1 * f as f64,
                        a as f64 * 2.0,
                        -(a as f64) + 0.05 * f as f64,
                    )
                })
                .collect();
            traj.push_frame(frame).unwrap();
        }
        traj
    }

    #[test]
    fn rejects_bad_construction() {
        assert!(Trajectory::new(0, 0.001).is_err());
        assert!(Trajectory::new(3, 0.0).is_err());
        assert!(Trajectory::new(3, -1.0).is_err());
    }

    #[test]
    fn push_frame_checks_atom_count() {
        let mut traj = Trajectory::new(3, 0.001).unwrap();
        assert!(traj.push_frame(vec![Vector3::zeros(); 2]).is_err());
        assert!(traj.push_frame(vec![Vector3::zeros(); 3]).is_ok());
        assert_eq!(traj.len(), 1);
    }

    #[test]
    fn binary_round_trip_is_lossless_to_f32() {
        let traj = sample_trajectory();
        let bytes = write_binary(&traj);
        let back = read_binary(&bytes).unwrap();
        assert_eq!(back.len(), traj.len());
        assert_eq!(back.natoms(), traj.natoms());
        assert!((back.dt() - traj.dt()).abs() < 1e-9);
        for (fa, fb) in back.frames().iter().zip(traj.frames()) {
            for (a, b) in fa.iter().zip(fb) {
                // f32 storage -> ~1e-6 relative tolerance.
                assert!((a - b).norm() < 1e-4, "{a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn text_round_trip_is_exact() {
        let traj = sample_trajectory();
        let text = write_text(&traj);
        let back = read_text(&text).unwrap();
        assert_eq!(back.len(), traj.len());
        for (fa, fb) in back.frames().iter().zip(traj.frames()) {
            for (a, b) in fa.iter().zip(fb) {
                assert!((a - b).norm() < 1e-7);
            }
        }
    }

    #[test]
    fn binary_rejects_corrupt_input() {
        assert!(read_binary(b"too short").is_err());
        let mut bytes = write_binary(&sample_trajectory());
        bytes[0] = b'X'; // break the magic
        assert!(read_binary(&bytes).is_err());
        // Truncate the body.
        let truncated = &write_binary(&sample_trajectory())[..40];
        assert!(read_binary(truncated).is_err());
    }

    #[test]
    fn text_rejects_corrupt_input() {
        assert!(read_text("").is_err());
        // Missing natoms.
        assert!(read_text("# trajectory dt=0.001\nFRAME 0\n").is_err());
        // Wrong-length frame.
        let bad = "# natoms=3 dt=0.001\nFRAME 0\n1 2 3\n4 5 6\n";
        assert!(read_text(bad).is_err());
    }

    #[test]
    fn empty_trajectory_round_trips() {
        let traj = Trajectory::new(4, 0.001).unwrap();
        let back = read_binary(&write_binary(&traj)).unwrap();
        assert!(back.is_empty());
        assert_eq!(back.natoms(), 4);
        let back_text = read_text(&write_text(&traj)).unwrap();
        assert!(back_text.is_empty());
    }
}

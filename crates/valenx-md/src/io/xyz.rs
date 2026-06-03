//! XYZ reader and writer — **roadmap feature 4**.
//!
//! The XYZ format is the simplest molecular file there is:
//!
//! ```text
//! <atom count>
//! <free-form comment line>
//! <element>  <x>  <y>  <z>
//! <element>  <x>  <y>  <z>
//! ...
//! ```
//!
//! Coordinates are Ångström by convention; this module converts to/
//! from the crate's nm. A file may hold a single frame or — for a
//! trajectory — many frames back to back, each with its own
//! count/comment header. [`read_xyz`] returns the first frame as a
//! [`System`]; [`read_xyz_frames`] returns every frame's coordinates.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::io::pdb::element_mass;
use crate::io::ANGSTROM_PER_NM;
use crate::system::{Atom, System, Topology};

/// Parses the first frame of an XYZ string into a [`System`] (no box,
/// no bonds).
///
/// # Errors
/// [`MdError::Parse`] on a missing / non-numeric atom count or a
/// malformed coordinate line.
pub fn read_xyz(text: &str) -> Result<System> {
    let mut lines = text.lines();
    let (topology, positions) = read_one_frame(&mut lines)?
        .ok_or_else(|| MdError::parse("xyz", "file is empty"))?;
    System::new(topology, positions)
}

/// Parses *every* frame of a multi-frame XYZ string, returning each
/// frame's coordinate list (nm).
///
/// The element list is taken from the first frame; all frames must
/// agree on the atom count.
///
/// # Errors
/// [`MdError::Parse`] on any malformed frame; [`MdError::DimensionMismatch`]
/// if frames disagree on atom count.
pub fn read_xyz_frames(text: &str) -> Result<Vec<Vec<Vector3<f64>>>> {
    let mut lines = text.lines();
    let mut frames = Vec::new();
    let mut expected: Option<usize> = None;
    while let Some((topology, positions)) = read_one_frame(&mut lines)? {
        match expected {
            None => expected = Some(topology.len()),
            Some(n) if n != topology.len() => {
                return Err(MdError::dimension(format!(
                    "XYZ frame has {} atoms, expected {n}",
                    topology.len()
                )));
            }
            _ => {}
        }
        frames.push(positions);
    }
    if frames.is_empty() {
        return Err(MdError::parse("xyz", "file is empty"));
    }
    Ok(frames)
}

/// Reads one frame from a line iterator. Returns `Ok(None)` at clean
/// end-of-input.
fn read_one_frame<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<Option<(Topology, Vec<Vector3<f64>>)>> {
    // Skip blank lines before the count.
    let count_line = loop {
        match lines.next() {
            None => return Ok(None),
            Some(l) if l.trim().is_empty() => continue,
            Some(l) => break l,
        }
    };
    let count: usize = count_line.trim().parse().map_err(|_| {
        MdError::parse("xyz", format!("expected an atom count, got `{count_line}`"))
    })?;
    // The comment line (may be blank, but must exist).
    let _comment = lines
        .next()
        .ok_or_else(|| MdError::parse("xyz", "missing comment line"))?;

    let mut topology = Topology::new();
    // Do NOT pre-size from `count`: it is caller-controlled, so
    // `Vec::with_capacity(count)` on an adversarial header (e.g.
    // `99999999999`) would attempt a huge allocation and abort the
    // process. Grow on demand — the per-atom `lines.next()` below
    // already bounds the loop by the real input length.
    let mut positions: Vec<Vector3<f64>> = Vec::new();
    for atom_index in 0..count {
        let line = lines.next().ok_or_else(|| {
            MdError::parse(
                "xyz",
                format!("expected {count} atoms, file ended after {atom_index}"),
            )
        })?;
        let mut fields = line.split_whitespace();
        let element = fields
            .next()
            .ok_or_else(|| MdError::parse("xyz", "empty atom line"))?
            .to_string();
        let parse = |s: Option<&str>, what: &str| -> Result<f64> {
            s.ok_or_else(|| MdError::parse("xyz", format!("missing {what} coordinate")))?
                .parse::<f64>()
                .map_err(|_| MdError::parse("xyz", format!("bad {what} coordinate")))
        };
        let x = parse(fields.next(), "x")?;
        let y = parse(fields.next(), "y")?;
        let z = parse(fields.next(), "z")?;
        let mass = element_mass(&element);
        let atom = Atom::new(element.clone(), mass, 0.0)
            .map_err(|e| MdError::parse("xyz", e.to_string()))?
            .with_element(element);
        topology.push_atom(atom);
        positions.push(Vector3::new(
            x / ANGSTROM_PER_NM,
            y / ANGSTROM_PER_NM,
            z / ANGSTROM_PER_NM,
        ));
    }
    Ok(Some((topology, positions)))
}

/// Serialises a [`System`] to a single-frame XYZ string.
///
/// `comment` becomes the second line; pass `""` for a blank comment.
pub fn write_xyz(system: &System, comment: &str) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}\n", system.len()));
    out.push_str(comment);
    out.push('\n');
    for (atom, pos) in system.topology.atoms.iter().zip(&system.positions) {
        let symbol = if atom.element.is_empty() {
            atom.type_name.as_str()
        } else {
            atom.element.as_str()
        };
        out.push_str(&format!(
            "{:<3} {:14.8} {:14.8} {:14.8}\n",
            symbol,
            pos.x * ANGSTROM_PER_NM,
            pos.y * ANGSTROM_PER_NM,
            pos.z * ANGSTROM_PER_NM,
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "\
3
water molecule
O    0.000   0.000   0.000
H    0.957   0.000   0.000
H   -0.240   0.927   0.000
";

    #[test]
    fn reads_a_single_frame() {
        let sys = read_xyz(WATER).unwrap();
        assert_eq!(sys.len(), 3);
        assert_eq!(sys.topology.atoms[0].element, "O");
        assert!((sys.topology.atoms[0].mass - 15.999).abs() < 1e-3);
        // 0.957 Å -> 0.0957 nm.
        assert!((sys.positions[1].x - 0.0957).abs() < 1e-6);
    }

    #[test]
    fn round_trip_preserves_coordinates() {
        let sys = read_xyz(WATER).unwrap();
        let text = write_xyz(&sys, "round trip");
        let back = read_xyz(&text).unwrap();
        assert_eq!(back.len(), 3);
        for (a, b) in back.positions.iter().zip(&sys.positions) {
            assert!((a - b).norm() < 1e-8);
        }
    }

    #[test]
    fn reads_a_multi_frame_trajectory() {
        let traj = format!("{WATER}{WATER}{WATER}");
        let frames = read_xyz_frames(&traj).unwrap();
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[0].len(), 3);
    }

    #[test]
    fn rejects_malformed_files() {
        assert!(read_xyz("").is_err());
        // Bad count.
        assert!(read_xyz("notanumber\ncomment\n").is_err());
        // Truncated.
        assert!(read_xyz("3\ncomment\nO 0 0 0\n").is_err());
        // Missing coordinate.
        assert!(read_xyz("1\nc\nO 0 0\n").is_err());
    }

    #[test]
    fn frame_count_mismatch_is_rejected() {
        let bad = "2\na\nO 0 0 0\nH 1 0 0\n1\nb\nO 0 0 0\n";
        assert!(read_xyz_frames(bad).is_err());
    }

    #[test]
    fn absurd_atom_count_is_rejected_without_a_huge_allocation() {
        // A header claiming ~1e11 atoms must not drive a multi-hundred-
        // gigabyte `Vec::with_capacity`; the truncation check rejects it.
        let bad = "99999999999\ncomment\nO 0 0 0\n";
        assert!(read_xyz(bad).is_err());
        assert!(read_xyz_frames(bad).is_err());
    }

    #[test]
    fn garbage_input_never_panics() {
        for bad in ["", "\0\0\0", "x", "1", "1\n", "-5\nc\n", "1\nc\nO a b c\n"] {
            let _ = read_xyz(bad);
            let _ = read_xyz_frames(bad);
        }
    }
}

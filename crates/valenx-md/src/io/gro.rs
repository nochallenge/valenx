//! GRO reader and writer — **roadmap feature 5**.
//!
//! The GROMACS `.gro` format is a fixed-column structure file that —
//! unlike PDB — natively carries **velocities** and works in
//! **nanometres**, the crate's own length unit. Layout:
//!
//! ```text
//! <free-form title line>
//! <atom count>
//! <residue records...>
//! <box-vector line>
//! ```
//!
//! Each residue record is:
//!
//! ```text
//! COLUMNS  FIELD
//!  1- 5    residue number
//!  6-10    residue name
//! 11-15    atom name
//! 16-20    atom serial
//! 21-28    x   (nm,  %8.3f)
//! 29-36    y
//! 37-44    z
//! 45-52    vx  (nm/ps, optional)
//! 53-60    vy
//! 61-68    vz
//! ```
//!
//! The box line is three floats (an orthorhombic box) or nine (a
//! triclinic box in GROMACS order: `xx yy zz xy xz yx yz zx zy`).
//!
//! GRO carries no masses; the reader assigns one from the atom name's
//! leading element via the shared [`crate::io::pdb::element_mass`]
//! table.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::io::pdb::element_mass;
use crate::pbc::SimBox;
use crate::system::{Atom, System, Topology};

/// Element guess from a GRO atom name's leading letters.
fn guess_element(atom_name: &str) -> String {
    let name = atom_name.trim();
    let upper = name.to_ascii_uppercase();
    for two in ["CL", "NA", "MG", "FE", "ZN", "CA", "BR"] {
        if upper.starts_with(two) {
            return two.to_string();
        }
    }
    name.chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase().to_string())
        .unwrap_or_default()
}

/// Slices a fixed-column field, returning `""` past the line end.
fn col(line: &str, start: usize, end: usize) -> &str {
    let bytes = line.as_bytes();
    if start >= bytes.len() {
        return "";
    }
    let end = end.min(bytes.len());
    std::str::from_utf8(&bytes[start..end]).unwrap_or("").trim()
}

/// Parses a GRO string into a [`System`], including velocities and the
/// periodic box if present.
///
/// # Errors
/// [`MdError::Parse`] on a missing count, a malformed coordinate, or a
/// box line that is neither 3 nor 9 numbers.
pub fn read_gro(text: &str) -> Result<System> {
    let mut lines = text.lines();
    // Title.
    let _title = lines
        .next()
        .ok_or_else(|| MdError::parse("gro", "missing title line"))?;
    // Atom count.
    let count_line = lines
        .next()
        .ok_or_else(|| MdError::parse("gro", "missing atom-count line"))?;
    let count: usize = count_line
        .trim()
        .parse()
        .map_err(|_| MdError::parse("gro", format!("bad atom count `{}`", count_line.trim())))?;

    let mut topology = Topology::new();
    // Do NOT pre-size from the caller-controlled `count`: an adversarial
    // header (e.g. `99999999999`) would make `Vec::with_capacity`
    // attempt a huge allocation and abort the process. Grow on demand —
    // the per-atom `lines.next()` already bounds the loop by the real
    // input length.
    let mut positions: Vec<Vector3<f64>> = Vec::new();
    let mut velocities: Vec<Vector3<f64>> = Vec::new();
    let mut any_velocity = false;

    for atom_index in 0..count {
        let line = lines.next().ok_or_else(|| {
            MdError::parse(
                "gro",
                format!("expected {count} atoms, file ended after {atom_index}"),
            )
        })?;
        let residue_id: i32 = col(line, 0, 5).parse().unwrap_or(0);
        let residue = col(line, 5, 10).to_string();
        let atom_name = col(line, 10, 15).to_string();
        let parse_coord = |s: &str, what: &str| -> Result<f64> {
            s.parse::<f64>().map_err(|_| {
                MdError::parse(
                    "gro",
                    format!("atom {}: bad {what} value `{s}`", atom_index + 1),
                )
            })
        };
        let x = parse_coord(col(line, 20, 28), "x")?;
        let y = parse_coord(col(line, 28, 36), "y")?;
        let z = parse_coord(col(line, 36, 44), "z")?;
        // Velocities are optional.
        let vx_field = col(line, 44, 52);
        let (vx, vy, vz) = if vx_field.is_empty() {
            (0.0, 0.0, 0.0)
        } else {
            any_velocity = true;
            (
                parse_coord(vx_field, "vx")?,
                parse_coord(col(line, 52, 60), "vy")?,
                parse_coord(col(line, 60, 68), "vz")?,
            )
        };
        let element = guess_element(&atom_name);
        let mass = element_mass(&element);
        let type_name = if element.is_empty() {
            atom_name.clone()
        } else {
            element.clone()
        };
        let atom = Atom::new(type_name, mass, 0.0)
            .map_err(|e| MdError::parse("gro", e.to_string()))?
            .with_element(element)
            .with_name(atom_name)
            .with_residue(residue, residue_id);
        topology.push_atom(atom);
        positions.push(Vector3::new(x, y, z));
        velocities.push(Vector3::new(vx, vy, vz));
    }

    // Box line.
    let cell = match lines.next() {
        None => SimBox::none(),
        Some(box_line) => parse_box(box_line)?,
    };

    let mut system = System::new(topology, positions)?.with_cell(cell);
    if any_velocity {
        system.set_velocities(velocities)?;
    }
    Ok(system)
}

/// Parses the GRO box-vector line.
fn parse_box(line: &str) -> Result<SimBox> {
    let nums: Vec<f64> = line
        .split_whitespace()
        .map(|s| s.parse::<f64>())
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| MdError::parse("gro", "non-numeric value on the box line"))?;
    match nums.len() {
        0 => Ok(SimBox::none()),
        3 => SimBox::orthorhombic(nums[0], nums[1], nums[2]),
        9 => {
            // GROMACS order: xx yy zz xy xz yx yz zx zy. The lattice
            // vectors are the rows (a = (xx,xy,xz), ...).
            let a = Vector3::new(nums[0], nums[3], nums[4]);
            let b = Vector3::new(nums[5], nums[1], nums[6]);
            let c = Vector3::new(nums[7], nums[8], nums[2]);
            SimBox::triclinic(a, b, c)
        }
        n => Err(MdError::parse(
            "gro",
            format!("box line has {n} numbers, expected 3 or 9"),
        )),
    }
}

/// Serialises a [`System`] to a GRO string.
///
/// Velocities are written if any atom has a non-zero velocity. The box
/// line is the 3-number orthorhombic form for an orthorhombic cell,
/// the 9-number form otherwise, and absent for a non-periodic system.
pub fn write_gro(system: &System, title: &str) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str(&format!("{}\n", system.len()));

    let write_velocities = system.velocities.iter().any(|v| v.norm_squared() > 0.0);
    for (i, (atom, pos)) in system
        .topology
        .atoms
        .iter()
        .zip(&system.positions)
        .enumerate()
    {
        let serial = (i + 1) % 100_000;
        let resid = atom.residue_id.rem_euclid(100_000);
        let name = if atom.name.is_empty() {
            atom.type_name.as_str()
        } else {
            atom.name.as_str()
        };
        out.push_str(&format!(
            "{:>5}{:<5}{:>5}{:>5}{:8.3}{:8.3}{:8.3}",
            resid,
            truncate(&atom.residue, 5),
            truncate(name, 5),
            serial,
            pos.x,
            pos.y,
            pos.z,
        ));
        if write_velocities {
            let v = system.velocities[i];
            out.push_str(&format!("{:8.4}{:8.4}{:8.4}", v.x, v.y, v.z));
        }
        out.push('\n');
    }

    // Box line.
    if system.cell.is_periodic() {
        if system.cell.is_orthorhombic() {
            let [a, b, c] = system.cell.edge_lengths();
            out.push_str(&format!("{a:10.5}{b:10.5}{c:10.5}\n"));
        } else {
            let h = system.cell.matrix();
            // Rows are the lattice vectors; emit GROMACS order.
            out.push_str(&format!(
                "{:10.5}{:10.5}{:10.5}{:10.5}{:10.5}{:10.5}{:10.5}{:10.5}{:10.5}\n",
                h.m11, h.m22, h.m33, h.m12, h.m13, h.m21, h.m23, h.m31, h.m32,
            ));
        }
    } else {
        out.push_str("   0.00000   0.00000   0.00000\n");
    }
    out
}

/// Truncates a string to at most `n` characters (GRO fields are fixed
/// width).
fn truncate(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
MD system
4
    1SOL     OW    1   0.230   0.628   0.113   0.123  -0.456   0.789
    1SOL    HW1    2   0.137   0.626   0.150  -0.100   0.200  -0.300
    1SOL    HW2    3   0.231   0.589   0.021   0.050  -0.050   0.150
    2NA      NA    4   1.500   1.500   1.500   0.000   0.000   0.000
   3.00000   3.00000   3.00000
";

    #[test]
    fn reads_atoms_positions_velocities_and_box() {
        let sys = read_gro(SAMPLE).unwrap();
        assert_eq!(sys.len(), 4);
        // GRO is already in nm.
        assert!((sys.positions[0].x - 0.230).abs() < 1e-6);
        // Velocities read.
        assert!((sys.velocities[0].x - 0.123).abs() < 1e-6);
        assert!((sys.velocities[0].y - (-0.456)).abs() < 1e-6);
        // Box.
        assert!(sys.cell.is_periodic());
        assert!((sys.cell.volume() - 27.0).abs() < 1e-6);
    }

    #[test]
    fn reads_residue_metadata() {
        let sys = read_gro(SAMPLE).unwrap();
        assert_eq!(sys.topology.atoms[0].residue, "SOL");
        assert_eq!(sys.topology.atoms[0].name, "OW");
        assert_eq!(sys.topology.atoms[3].name, "NA");
        // NA atom name -> sodium element + mass.
        assert!((sys.topology.atoms[3].mass - 22.990).abs() < 1e-3);
    }

    #[test]
    fn round_trip_preserves_positions_and_velocities() {
        let sys = read_gro(SAMPLE).unwrap();
        let text = write_gro(&sys, "round trip");
        let back = read_gro(&text).unwrap();
        assert_eq!(back.len(), sys.len());
        for (a, b) in back.positions.iter().zip(&sys.positions) {
            assert!((a - b).norm() < 1e-3);
        }
        for (a, b) in back.velocities.iter().zip(&sys.velocities) {
            assert!((a - b).norm() < 1e-3);
        }
        assert!((back.cell.volume() - sys.cell.volume()).abs() < 1e-3);
    }

    #[test]
    fn handles_a_file_without_velocities() {
        let no_vel = "\
no velocities
2
    1MOL      C    1   0.100   0.200   0.300
    1MOL      C    2   0.400   0.500   0.600
   2.00000   2.00000   2.00000
";
        let sys = read_gro(no_vel).unwrap();
        assert_eq!(sys.len(), 2);
        for v in &sys.velocities {
            assert!(v.norm() < 1e-12);
        }
    }

    #[test]
    fn triclinic_box_round_trips() {
        let tri = "\
triclinic
1
    1MOL      C    1   0.100   0.200   0.300
   3.00000   3.00000   3.00000   0.50000   0.00000   0.00000   0.50000   0.00000   0.00000
";
        let sys = read_gro(tri).unwrap();
        assert!(sys.cell.is_periodic());
        assert!(!sys.cell.is_orthorhombic());
        let text = write_gro(&sys, "tri rt");
        let back = read_gro(&text).unwrap();
        assert!((back.cell.volume() - sys.cell.volume()).abs() < 1e-3);
    }

    #[test]
    fn rejects_malformed_files() {
        assert!(read_gro("").is_err());
        assert!(read_gro("title\nnotanumber\n").is_err());
        // Truncated atom list.
        assert!(read_gro("title\n3\n    1MOL C 1 0.1 0.2 0.3\n").is_err());
    }

    #[test]
    fn absurd_atom_count_is_rejected_without_a_huge_allocation() {
        // ~1e11 declared atoms must not drive a huge `Vec::with_capacity`.
        let bad = "title\n99999999999\n    1MOL C 1 0.1 0.2 0.3\n";
        assert!(read_gro(bad).is_err());
    }

    #[test]
    fn garbage_input_never_panics() {
        for bad in ["", "t", "t\n", "t\n2\n", "t\n-1\nx\n", "\u{0}\u{0}"] {
            let _ = read_gro(bad);
        }
    }
}

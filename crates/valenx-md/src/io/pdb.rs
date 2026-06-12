//! PDB reader and writer — **roadmap feature 3**.
//!
//! The Protein Data Bank format is the lingua franca of molecular
//! structure. This module reads and writes the subset an MD engine
//! needs: `ATOM` / `HETATM` coordinate records and the `CRYST1` unit-
//! cell record.
//!
//! PDB is a **fixed-column** format — fields are identified by their
//! character positions, not by whitespace. The reader honours the
//! official columns:
//!
//! ```text
//! COLUMNS  FIELD
//!  1- 6    record name  ("ATOM  " / "HETATM")
//!  7-11    serial
//! 13-16    atom name
//! 18-20    residue name
//! 22       chain id
//! 23-26    residue sequence number
//! 31-38    x  (Å)
//! 39-46    y  (Å)
//! 47-54    z  (Å)
//! 77-78    element symbol
//! ```
//!
//! Coordinates are Ångström in the file and nanometre in the returned
//! [`System`]; the reader divides by 10 and the writer multiplies.
//!
//! Masses are not stored in a PDB file, so the reader assigns each
//! atom a mass from a built-in element table (see
//! [`crate::io::pdb::element_mass`]); a structure read this way is
//! ready for geometry analysis and, once a force field is attached,
//! for dynamics.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::io::ANGSTROM_PER_NM;
use crate::pbc::SimBox;
use crate::system::{Atom, System, Topology};

/// Returns a representative atomic mass (u) for an element symbol.
///
/// Covers the elements common in biomolecular structures; an unknown
/// symbol falls back to carbon's mass so a read never fails on an
/// exotic atom — the caller can override the mass afterwards.
pub fn element_mass(symbol: &str) -> f64 {
    match symbol.trim().to_ascii_uppercase().as_str() {
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "NA" => 22.990,
        "MG" => 24.305,
        "P" => 30.974,
        "S" => 32.06,
        "CL" => 35.45,
        "K" => 39.098,
        "CA" => 40.078,
        "FE" => 55.845,
        "ZN" => 65.38,
        "BR" => 79.904,
        "I" => 126.904,
        _ => 12.011,
    }
}

/// Best-effort element guess from a 4-character PDB atom name.
fn guess_element(atom_name: &str) -> String {
    let name = atom_name.trim();
    // Two-letter elements that show up in biomolecules.
    let upper = name.to_ascii_uppercase();
    for two in ["CL", "NA", "MG", "FE", "ZN", "CA", "BR"] {
        if upper.starts_with(two) {
            return two.to_string();
        }
    }
    // Otherwise the first alphabetic character is the element.
    name.chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase().to_string())
        .unwrap_or_default()
}

/// Slices a fixed-column field `[start, end)` (0-based) from a record,
/// returning `""` if the line is too short.
fn col(line: &str, start: usize, end: usize) -> &str {
    let bytes = line.as_bytes();
    if start >= bytes.len() {
        return "";
    }
    let end = end.min(bytes.len());
    // PDB is ASCII; slicing on byte indices is safe here.
    std::str::from_utf8(&bytes[start..end]).unwrap_or("").trim()
}

/// Parses a PDB string into a [`System`].
///
/// `ATOM` and `HETATM` records become atoms; a `CRYST1` record, if
/// present, becomes the periodic box (orthorhombic cells only — a
/// non-orthorhombic `CRYST1` yields a [`MdError::Parse`]). The
/// topology has no bonds (PDB `CONECT` records are not read in v1).
///
/// # Errors
/// [`MdError::Parse`] on a malformed coordinate or `CRYST1` record.
pub fn read_pdb(text: &str) -> Result<System> {
    let mut topology = Topology::new();
    let mut positions = Vec::new();
    let mut cell = SimBox::none();

    for (lineno, line) in text.lines().enumerate() {
        let record = col(line, 0, 6);
        if record == "ATOM" || record == "HETATM" {
            let name = col(line, 12, 16).to_string();
            let residue = col(line, 17, 20).to_string();
            let chain = col(line, 21, 22).to_string();
            let residue_id: i32 = col(line, 22, 26).parse().unwrap_or(0);
            let parse_coord = |s: &str, what: &str| -> Result<f64> {
                s.parse::<f64>().map_err(|_| {
                    MdError::parse(
                        "pdb",
                        format!("line {}: bad {what} coordinate `{s}`", lineno + 1),
                    )
                })
            };
            let x = parse_coord(col(line, 30, 38), "x")?;
            let y = parse_coord(col(line, 38, 46), "y")?;
            let z = parse_coord(col(line, 46, 54), "z")?;
            let mut element = col(line, 76, 78).to_string();
            if element.is_empty() {
                element = guess_element(&name);
            }
            let mass = element_mass(&element);
            let type_name = if element.is_empty() {
                name.clone()
            } else {
                element.clone()
            };
            let atom = Atom::new(type_name, mass, 0.0)
                .map_err(|e| MdError::parse("pdb", format!("line {}: {e}", lineno + 1)))?
                .with_element(element)
                .with_name(name)
                .with_residue(residue, residue_id)
                .with_chain(chain);
            topology.push_atom(atom);
            positions.push(Vector3::new(
                x / ANGSTROM_PER_NM,
                y / ANGSTROM_PER_NM,
                z / ANGSTROM_PER_NM,
            ));
        } else if record == "CRYST1" {
            let a: f64 = col(line, 6, 15).parse().unwrap_or(0.0);
            let b: f64 = col(line, 15, 24).parse().unwrap_or(0.0);
            let c: f64 = col(line, 24, 33).parse().unwrap_or(0.0);
            let alpha: f64 = col(line, 33, 40).parse().unwrap_or(90.0);
            let beta: f64 = col(line, 40, 47).parse().unwrap_or(90.0);
            let gamma: f64 = col(line, 47, 54).parse().unwrap_or(90.0);
            if a > 0.0 && b > 0.0 && c > 0.0 {
                let orthogonal = (alpha - 90.0).abs() < 1e-3
                    && (beta - 90.0).abs() < 1e-3
                    && (gamma - 90.0).abs() < 1e-3;
                if !orthogonal {
                    return Err(MdError::parse(
                        "pdb",
                        "non-orthorhombic CRYST1 cells are not supported in v1",
                    ));
                }
                cell = SimBox::orthorhombic(
                    a / ANGSTROM_PER_NM,
                    b / ANGSTROM_PER_NM,
                    c / ANGSTROM_PER_NM,
                )?;
            }
        }
        // TER / END / other records are ignored.
    }

    if topology.is_empty() {
        return Err(MdError::parse("pdb", "no ATOM / HETATM records found"));
    }
    Ok(System::new(topology, positions)?.with_cell(cell))
}

/// Serialises a [`System`] to a PDB string.
///
/// Writes a `CRYST1` record for a periodic orthorhombic box, then one
/// `ATOM` record per atom (coordinates converted nm → Å), then `END`.
pub fn write_pdb(system: &System) -> String {
    let mut out = String::new();
    // CRYST1 for a periodic box.
    if system.cell.is_periodic() {
        let [a, b, c] = system.cell.edge_lengths();
        out.push_str(&format!(
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} P 1           1\n",
            a * ANGSTROM_PER_NM,
            b * ANGSTROM_PER_NM,
            c * ANGSTROM_PER_NM,
            90.0,
            90.0,
            90.0
        ));
    }
    for (i, (atom, pos)) in system
        .topology
        .atoms
        .iter()
        .zip(&system.positions)
        .enumerate()
    {
        let serial = (i + 1) % 100_000;
        let resid = atom.residue_id.rem_euclid(10_000);
        let chain = atom.chain.chars().next().unwrap_or(' ');
        // Atom name is left-justified in columns 13-16 for >3 chars,
        // otherwise conventionally shifted; keep it simple and
        // left-justify within the 4-wide field.
        let name = if atom.name.is_empty() {
            atom.type_name.as_str()
        } else {
            atom.name.as_str()
        };
        out.push_str(&format!(
            "ATOM  {serial:>5} {name:<4} {res:>3} {chain}{resid:>4}    \
{x:8.3}{y:8.3}{z:8.3}{occ:6.2}{temp:6.2}          {elem:>2}\n",
            serial = serial,
            name = name,
            res = atom.residue,
            chain = chain,
            resid = resid,
            x = pos.x * ANGSTROM_PER_NM,
            y = pos.y * ANGSTROM_PER_NM,
            z = pos.z * ANGSTROM_PER_NM,
            occ = 1.0,
            temp = 0.0,
            elem = atom.element,
        ));
    }
    out.push_str("END\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C
ATOM      3  C   ALA A   1      13.140   6.099  -5.255  1.00  0.00           C
HETATM    4  O   HOH A   2      20.000  20.000  20.000  1.00  0.00           O
END
";

    #[test]
    fn reads_atoms_and_box() {
        let sys = read_pdb(SAMPLE).unwrap();
        assert_eq!(sys.len(), 4);
        assert!(sys.cell.is_periodic());
        // 30 Å -> 3 nm.
        assert!((sys.cell.edge_lengths()[0] - 3.0).abs() < 1e-6);
        // First atom is nitrogen.
        assert_eq!(sys.topology.atoms[0].element, "N");
        assert!((sys.topology.atoms[0].mass - 14.007).abs() < 1e-3);
        // Coordinates converted to nm.
        assert!((sys.positions[0].x - 1.1104).abs() < 1e-6);
    }

    #[test]
    fn reads_residue_and_chain_metadata() {
        let sys = read_pdb(SAMPLE).unwrap();
        assert_eq!(sys.topology.atoms[1].name, "CA");
        assert_eq!(sys.topology.atoms[1].residue, "ALA");
        assert_eq!(sys.topology.atoms[1].chain, "A");
        assert_eq!(sys.topology.atoms[3].residue, "HOH");
    }

    #[test]
    fn round_trip_preserves_coordinates() {
        let sys = read_pdb(SAMPLE).unwrap();
        let text = write_pdb(&sys);
        let back = read_pdb(&text).unwrap();
        assert_eq!(back.len(), sys.len());
        for (a, b) in back.positions.iter().zip(&sys.positions) {
            assert!((a - b).norm() < 1e-3, "{a:?} vs {b:?}");
        }
        assert!((back.cell.edge_lengths()[0] - sys.cell.edge_lengths()[0]).abs() < 1e-3);
    }

    #[test]
    fn rejects_empty_and_malformed() {
        assert!(read_pdb("").is_err());
        assert!(read_pdb("REMARK just a comment\n").is_err());
        let bad = "ATOM      1  N   ALA A   1      XX.XXX   6.134  -6.504\n";
        assert!(read_pdb(bad).is_err());
    }

    #[test]
    fn rejects_triclinic_cryst1() {
        let triclinic = "\
CRYST1   30.000   30.000   30.000  90.00  60.00  90.00 P 1           1
ATOM      1  N   ALA A   1      11.000   6.000  -6.000  1.00  0.00           N
END
";
        assert!(read_pdb(triclinic).is_err());
    }

    #[test]
    fn element_mass_table_is_reasonable() {
        assert!((element_mass("C") - 12.011).abs() < 1e-3);
        assert!((element_mass(" o ") - 15.999).abs() < 1e-3);
        // Unknown falls back to carbon, not a panic.
        assert!((element_mass("Xx") - 12.011).abs() < 1e-3);
    }
}

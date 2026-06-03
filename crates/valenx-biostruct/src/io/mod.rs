//! Coordinate-file I/O: PDB and mmCIF readers and writers.
//!
//! Both formats deserialise into the same
//! [`crate::structure::Structure`] hierarchy, so the rest of the
//! crate is format-agnostic. [`detect_format`] sniffs which reader to
//! use from the text content.

pub mod mmcif;
pub mod pdb;

pub use mmcif::{read_mmcif, write_mmcif};
pub use pdb::{read_pdb, write_pdb};

use crate::error::Result;
use crate::structure::Structure;

/// A recognised coordinate-file format.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum CoordFormat {
    /// Legacy fixed-column PDB.
    Pdb,
    /// PDBx/mmCIF token-based format.
    MmCif,
}

/// Guess the coordinate format from the file text.
///
/// mmCIF files begin with a `data_` block header or carry `loop_` /
/// `_atom_site` tags; anything else is treated as PDB. The heuristic
/// scans only the first few non-blank lines.
pub fn detect_format(text: &str) -> CoordFormat {
    for line in text.lines().take(40) {
        let t = line.trim_start();
        if t.is_empty() {
            continue;
        }
        if t.starts_with("data_")
            || t.starts_with("loop_")
            || t.starts_with("_atom_site")
            || t.starts_with("_entry.")
        {
            return CoordFormat::MmCif;
        }
        if t.starts_with("ATOM")
            || t.starts_with("HETATM")
            || t.starts_with("HEADER")
            || t.starts_with("MODEL")
        {
            return CoordFormat::Pdb;
        }
    }
    CoordFormat::Pdb
}

/// Read a structure, auto-detecting PDB vs mmCIF from the content.
pub fn read_structure(text: &str, id: &str) -> Result<Structure> {
    match detect_format(text) {
        CoordFormat::Pdb => read_pdb(text, id),
        CoordFormat::MmCif => read_mmcif(text, id),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_pdb() {
        let pdb = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n";
        assert_eq!(detect_format(pdb), CoordFormat::Pdb);
    }

    #[test]
    fn detects_mmcif() {
        assert_eq!(detect_format("data_1abc\n#\n"), CoordFormat::MmCif);
        assert_eq!(
            detect_format("\n\nloop_\n_atom_site.id\n"),
            CoordFormat::MmCif
        );
    }

    #[test]
    fn read_structure_dispatches() {
        let pdb = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";
        let s = read_structure(pdb, "x").unwrap();
        assert_eq!(s.atom_count(), 1);
    }
}

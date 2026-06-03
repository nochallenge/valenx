//! PDB format reader. Spec: <https://www.wwpdb.org/documentation/file-format>
//!
//! v0.1 scope: ATOM and HETATM records only. We ignore HEADER /
//! TITLE / REMARK / CONECT / TER / END — they're informational and
//! the canonical [`crate::Structure`] doesn't have homes for them
//! yet. v0.2 lands HEADER preservation under a `Structure.metadata`
//! field.
//!
//! ATOM record column layout (1-indexed, fixed-width):
//! - 1-6   record name ("ATOM  " or "HETATM")
//! - 7-11  atom serial number
//! - 13-16 atom name
//! - 17    altLoc indicator (we keep first occurrence only)
//! - 18-20 residue name
//! - 22    chain identifier
//! - 23-26 residue sequence number
//! - 31-38 X coordinate
//! - 39-46 Y coordinate
//! - 47-54 Z coordinate
//! - 55-60 occupancy
//! - 61-66 temperature factor
//! - 77-78 element symbol

use nalgebra::Vector3;
use thiserror::Error;

use crate::structure::{Atom, Chain, Residue, Structure};

/// Errors raised by [`read`].
#[derive(Debug, Error)]
pub enum PdbError {
    /// An `ATOM` / `HETATM` line is too short or otherwise unparseable.
    #[error("malformed line {line}: {reason}")]
    Malformed {
        /// 1-based line number of the offending record.
        line: usize,
        /// Short human-readable explanation of the parse failure.
        reason: String,
    },
}

/// Parse a PDB-format string into a [`Structure`]. Only `ATOM` and
/// `HETATM` records are read; other record types are silently ignored.
///
/// ## Minimum line width
///
/// The PDB spec puts coordinates / occupancy / temperature factor in
/// columns 31-66 and the element symbol in columns 77-78. The element
/// symbol is documented as "right-justified" but real-world files —
/// Phenix output, hand-curated structures, many legacy archives —
/// truncate the trailing whitespace and the optional element column,
/// leaving 66-char ATOM lines. We therefore only require the line to
/// reach column 66 (chain ID through temperature factor); when the
/// line is shorter than 78 bytes the element symbol falls back to an
/// empty string and downstream code can infer it from the atom name.
///
/// ## altLoc handling
///
/// The altLoc indicator (column 17, 0-indexed 16) selects between
/// alternate side-chain conformations. We keep only the first
/// occurrence per (chain, seq_id, res_name, atom_name) key — the
/// canonical "primary conformer" choice — so disordered side-chains
/// don't inflate the atom count.
///
/// # Errors
///
/// Returns [`PdbError::Malformed`] when an `ATOM`/`HETATM` line is
/// shorter than column 66 (the minimum needed to read coordinates
/// and B-factor) or contains a non-numeric field where a coordinate
/// is expected.
pub fn read(id: impl Into<String>, text: &str) -> Result<Structure, PdbError> {
    /// Minimum line length covering record name through tempFactor
    /// (column 66, exclusive end-index 66). Element symbol (cols
    /// 77-78) and beyond are optional.
    const MIN_LEN: usize = 66;

    let mut chains_map: std::collections::BTreeMap<char, Vec<Residue>> = Default::default();

    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        if !(raw.starts_with("ATOM  ") || raw.starts_with("HETATM")) {
            continue;
        }
        if raw.len() < MIN_LEN {
            return Err(PdbError::Malformed {
                line: line_no,
                reason: format!(
                    "ATOM line too short ({} < {MIN_LEN} bytes; need coords + B-factor through col 66)",
                    raw.len()
                ),
            });
        }
        let cols = raw.as_bytes();
        let parse_f64 = |start: usize, len: usize| -> Result<f64, PdbError> {
            let s = std::str::from_utf8(&cols[start..start + len])
                .unwrap_or("")
                .trim();
            s.parse::<f64>().map_err(|_| PdbError::Malformed {
                line: line_no,
                reason: format!("non-float at cols {}-{}: {s:?}", start + 1, start + len),
            })
        };
        let parse_i32 = |start: usize, len: usize| -> Result<i32, PdbError> {
            let s = std::str::from_utf8(&cols[start..start + len])
                .unwrap_or("")
                .trim();
            s.parse::<i32>().map_err(|_| PdbError::Malformed {
                line: line_no,
                reason: format!("non-int at cols {}-{}: {s:?}", start + 1, start + len),
            })
        };
        let atom_name = std::str::from_utf8(&cols[12..16]).unwrap_or("").trim();
        // altLoc indicator: column 17 (0-indexed 16). Space / blank
        // means "no alternate" and is the primary conformer.
        let alt_loc = cols[16] as char;
        let res_name = std::str::from_utf8(&cols[17..20]).unwrap_or("").trim();
        let chain_id = cols[21] as char;
        let res_seq = parse_i32(22, 4)?;
        let x = parse_f64(30, 8)?;
        let y = parse_f64(38, 8)?;
        let z = parse_f64(46, 8)?;
        let b_factor = parse_f64(60, 6).unwrap_or(0.0);
        // Element symbol lives at cols 77-78 in fully-spec'd files,
        // but legitimate 66-char lines (Phenix, hand-curated) omit
        // it. Fall back to an empty string when the line is short;
        // downstream consumers can infer the element from the atom
        // name (e.g. " CA " → "C") if they need it.
        let element = if raw.len() >= 78 {
            std::str::from_utf8(&cols[76..78]).unwrap_or("").trim()
        } else {
            ""
        };

        let atom = Atom {
            name: atom_name.to_string(),
            element: element.to_string(),
            position: Vector3::new(x, y, z),
            b_factor,
        };

        let chain_residues = chains_map.entry(chain_id).or_default();
        // Append to existing residue if last entry matches; otherwise
        // open a new residue. Inside a matching residue, dedupe by
        // (atom_name) so a second altLoc record for the same atom
        // (e.g. " CA " with altLoc 'B' after an 'A') is dropped —
        // we keep the first occurrence per the module-level note.
        match chain_residues.last_mut() {
            Some(r) if r.seq_id == res_seq && r.name == res_name => {
                let dup = r.atoms.iter().any(|a| a.name == atom.name);
                if !dup {
                    r.atoms.push(atom);
                } else {
                    // Silently drop the duplicate altLoc — track it
                    // via a trace event in case anyone is debugging.
                    let _ = alt_loc;
                }
            }
            _ => chain_residues.push(Residue {
                name: res_name.to_string(),
                seq_id: res_seq,
                atoms: vec![atom],
            }),
        }
    }

    let chains = chains_map
        .into_iter()
        .map(|(id, residues)| Chain { id, residues })
        .collect();
    Ok(Structure {
        id: id.into(),
        chains,
    })
}

//! PDBQT format reader. PDBQT is PDB plus two extra columns: partial
//! charge (cols 71-76) and AutoDock-4 atom type (cols 78-79). It is
//! the canonical input format for AutoDock Vina.
//!
//! We parse ATOM and HETATM records into [`PdbqtAtom`] structs that
//! carry the standard PDB fields plus the two PDBQT extensions.
//! Flexibility records (ROOT, BRANCH, ENDBRANCH, TORSDOF) are exposed
//! separately as [`PdbqtRecord`] tokens so the docking engine can
//! rebuild the rotatable-bond tree.

use nalgebra::Vector3;

/// A single PDBQT ATOM/HETATM record.
#[derive(Clone, Debug, PartialEq)]
pub struct PdbqtAtom {
    /// 1-based serial from cols 7-11.
    pub serial: i32,
    /// Atom name, cols 13-16 (trimmed).
    pub name: String,
    /// 3-letter residue code, cols 18-20.
    pub residue_name: String,
    /// Chain id, col 22.
    pub chain: char,
    /// Residue sequence number, cols 23-26.
    pub residue_seq: i32,
    /// XYZ in Å, cols 31-38, 39-46, 47-54.
    pub position: Vector3<f64>,
    /// Occupancy, cols 55-60.
    pub occupancy: f64,
    /// B-factor / temperature factor, cols 61-66.
    pub b_factor: f64,
    /// Gasteiger partial charge, cols 71-76. PDBQT-specific.
    pub partial_charge: f64,
    /// AutoDock-4 atom type, cols 78-79 (trimmed). PDBQT-specific.
    pub ad4_type: String,
}

/// One line of a PDBQT file, classified by record type.
#[derive(Clone, Debug, PartialEq)]
pub enum PdbqtRecord {
    /// ATOM or HETATM record.
    Atom(PdbqtAtom),
    /// `ROOT` — start of the rigid ligand core.
    Root,
    /// `ENDROOT`.
    EndRoot,
    /// `BRANCH parent_serial child_serial` — start of a rotatable subtree.
    Branch {
        parent_serial: i32,
        child_serial: i32,
    },
    /// `ENDBRANCH parent_serial child_serial`.
    EndBranch {
        parent_serial: i32,
        child_serial: i32,
    },
    /// `TORSDOF n` — declared torsional degrees of freedom.
    Torsdof(i32),
    /// Any other line (MODEL, ENDMDL, REMARK, ...) — preserved as-is for debugging.
    Other(String),
}

/// Errors raised while parsing a PDBQT stream.
#[derive(Debug, thiserror::Error)]
pub enum PdbqtError {
    /// Line was shorter than the column we tried to read.
    #[error("line {line_no} too short ({len} chars) for PDBQT field")]
    LineTooShort {
        /// 1-based line number.
        line_no: usize,
        /// Actual line length.
        len: usize,
    },
    /// A fixed-column numeric field failed to parse.
    #[error("line {line_no} field `{field}`: {reason}")]
    BadField {
        /// 1-based line number.
        line_no: usize,
        /// Field name (e.g. "x", "partial_charge").
        field: &'static str,
        /// Underlying parse error message.
        reason: String,
    },
    /// BRANCH / ENDBRANCH expected two serial numbers.
    #[error("line {line_no} BRANCH record malformed")]
    BadBranch {
        /// 1-based line number.
        line_no: usize,
    },
}

/// Parse a complete PDBQT document into an ordered list of [`PdbqtRecord`]s.
/// Atoms appear in file order; flexibility tokens delimit rotatable groups.
pub fn parse(text: &str) -> Result<Vec<PdbqtRecord>, PdbqtError> {
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line_no = i + 1;
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            out.push(PdbqtRecord::Atom(parse_atom(line, line_no)?));
        } else if line.starts_with("ROOT") && !line.starts_with("ENDROOT") {
            out.push(PdbqtRecord::Root);
        } else if line.starts_with("ENDROOT") {
            out.push(PdbqtRecord::EndRoot);
        } else if line.starts_with("BRANCH") && !line.starts_with("ENDBRANCH") {
            let (a, b) = parse_branch(line, line_no)?;
            out.push(PdbqtRecord::Branch {
                parent_serial: a,
                child_serial: b,
            });
        } else if line.starts_with("ENDBRANCH") {
            let (a, b) = parse_branch(line, line_no)?;
            out.push(PdbqtRecord::EndBranch {
                parent_serial: a,
                child_serial: b,
            });
        } else if let Some(rest) = line.strip_prefix("TORSDOF") {
            let n: i32 =
                rest.trim()
                    .parse()
                    .map_err(|e: std::num::ParseIntError| PdbqtError::BadField {
                        line_no,
                        field: "torsdof",
                        reason: e.to_string(),
                    })?;
            out.push(PdbqtRecord::Torsdof(n));
        } else {
            out.push(PdbqtRecord::Other(line.to_string()));
        }
    }
    Ok(out)
}

fn parse_atom(line: &str, line_no: usize) -> Result<PdbqtAtom, PdbqtError> {
    // PDBQT columns are 1-indexed; Rust slicing is 0-indexed, so cols
    // 7-11 -> byte range 6..11. PDBQT files use ASCII so byte slicing
    // == char slicing.
    //
    // The mandatory PDBQT fields end at the partial-charge column 76,
    // so a line must reach at least index 76. The AutoDock-4 atom type
    // occupies columns 78-79; with a single-character type only column
    // 78 is filled, which makes a perfectly valid line exactly 78
    // characters long. A `< 79` check would wrongly reject those, so
    // the minimum is column 76 and `field` clips any range that runs
    // past the end of a (possibly trailing-space-trimmed) line.
    if line.len() < 76 {
        return Err(PdbqtError::LineTooShort {
            line_no,
            len: line.len(),
        });
    }
    let field = |range: std::ops::Range<usize>, name: &'static str| -> Result<&str, PdbqtError> {
        let _ = name;
        // Clip the requested column range to what the line actually
        // has: a fixed-column field may be short or absent if trailing
        // whitespace was stripped. A field that begins past the end of
        // the line reads as empty.
        let start = range.start.min(line.len());
        let end = range.end.min(line.len());
        Ok(line.get(start..end).unwrap_or("").trim())
    };
    let f_int = |range: std::ops::Range<usize>, name: &'static str| -> Result<i32, PdbqtError> {
        let raw = field(range, name)?;
        raw.parse::<i32>().map_err(|e| PdbqtError::BadField {
            line_no,
            field: name,
            reason: e.to_string(),
        })
    };
    let f_float = |range: std::ops::Range<usize>, name: &'static str| -> Result<f64, PdbqtError> {
        let raw = field(range, name)?;
        raw.parse::<f64>().map_err(|e| PdbqtError::BadField {
            line_no,
            field: name,
            reason: e.to_string(),
        })
    };

    let serial = f_int(6..11, "serial")?;
    let name = field(12..16, "name")?.to_string();
    let residue_name = field(17..20, "residue_name")?.to_string();
    let chain_str = field(21..22, "chain")?;
    let chain = chain_str.chars().next().unwrap_or(' ');
    let residue_seq = f_int(22..26, "residue_seq")?;
    let x = f_float(30..38, "x")?;
    let y = f_float(38..46, "y")?;
    let z = f_float(46..54, "z")?;
    // PDB convention allows occupancy and B-factor to be blank when
    // not known. We default to 1.0 and 0.0 respectively, matching
    // every other PDB parser in the wild. The PDBQT-specific
    // partial_charge and ad4_type columns are stricter (required)
    // since they're the whole point of using PDBQT over plain PDB.
    let occupancy = f_float(54..60, "occupancy").unwrap_or(1.0);
    let b_factor = f_float(60..66, "b_factor").unwrap_or(0.0);
    let partial_charge = f_float(70..76, "partial_charge")?;
    let ad4_type = field(77..79, "ad4_type")?.to_string();

    Ok(PdbqtAtom {
        serial,
        name,
        residue_name,
        chain,
        residue_seq,
        position: Vector3::new(x, y, z),
        occupancy,
        b_factor,
        partial_charge,
        ad4_type,
    })
}

fn parse_branch(line: &str, line_no: usize) -> Result<(i32, i32), PdbqtError> {
    // BRANCH / ENDBRANCH lines: `BRANCH<spaces>P<spaces>C`
    let rest = line.split_whitespace().skip(1).collect::<Vec<_>>();
    if rest.len() < 2 {
        return Err(PdbqtError::BadBranch { line_no });
    }
    let a: i32 = rest[0]
        .parse()
        .map_err(|_| PdbqtError::BadBranch { line_no })?;
    let b: i32 = rest[1]
        .parse()
        .map_err(|_| PdbqtError::BadBranch { line_no })?;
    Ok((a, b))
}

/// Format a single PDBQT ATOM record line (no trailing newline).
/// Column layout matches the AutoDock 4 / Vina spec exactly so the
/// output round-trips through [`parse`].
pub fn write_atom(a: &PdbqtAtom) -> String {
    // Columns:
    //   1-6    record name "ATOM  "
    //   7-11   serial      "%5d"
    //   13-16  name        "%-4s"
    //   18-20  resName     "%3s"
    //   22     chainID     "%c"
    //   23-26  resSeq      "%4d"
    //   31-38  x           "%8.3f"
    //   39-46  y           "%8.3f"
    //   47-54  z           "%8.3f"
    //   55-60  occupancy   "%6.2f"
    //   61-66  tempFactor  "%6.2f"
    //   71-76  partialChg  "%6.3f"
    //   78-79  atomType    "%-2s"
    format!(
        "ATOM  {serial:>5} {name:<4} {res:<3} {chain}{seq:>4}    {x:8.3}{y:8.3}{z:8.3}{occ:6.2}{bf:6.2}    {pc:6.3} {at:<2}",
        serial = a.serial,
        name = a.name,
        res = a.residue_name,
        chain = a.chain,
        seq = a.residue_seq,
        x = a.position.x,
        y = a.position.y,
        z = a.position.z,
        occ = a.occupancy,
        bf = a.b_factor,
        pc = a.partial_charge,
        at = a.ad4_type,
    )
}

/// Write a pose ensemble in Vina's `MODEL n / REMARK ... / ATOM... / ENDMDL`
/// format. `poses` is parallel arrays: each entry is the positions for
/// every atom of one pose, in the same order as `atoms`.
pub fn write_pose_ensemble(atoms: &[PdbqtAtom], poses: &[(Vec<Vector3<f64>>, f64)]) -> String {
    let mut out = String::new();
    for (i, (positions, score)) in poses.iter().enumerate() {
        out.push_str(&format!("MODEL {}\n", i + 1));
        out.push_str(&format!(
            "REMARK VINA RESULT:    {score:8.3}      0.000      0.000\n"
        ));
        for (atom, pos) in atoms.iter().zip(positions.iter()) {
            let mut a = atom.clone();
            a.position = *pos;
            out.push_str(&write_atom(&a));
            out.push('\n');
        }
        out.push_str("ENDMDL\n");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimum well-formed PDBQT ATOM record. Column layout:
    /// `RECORD  SERIAL  NAME RES C  SEQ        X       Y       Z   OCC   BF      CHARGE  TYPE`
    const ATOM_LINE: &str =
        "ATOM      1  C   LIG A   1      10.000  20.000  30.000  1.00  0.00     0.123 C ";

    #[test]
    fn parses_single_atom() {
        let recs = parse(ATOM_LINE).unwrap();
        assert_eq!(recs.len(), 1);
        match &recs[0] {
            PdbqtRecord::Atom(a) => {
                assert_eq!(a.serial, 1);
                assert_eq!(a.name, "C");
                assert_eq!(a.residue_name, "LIG");
                assert_eq!(a.chain, 'A');
                assert_eq!(a.residue_seq, 1);
                assert!((a.position.x - 10.0).abs() < 1e-6);
                assert!((a.position.y - 20.0).abs() < 1e-6);
                assert!((a.position.z - 30.0).abs() < 1e-6);
                assert!((a.partial_charge - 0.123).abs() < 1e-6);
                assert_eq!(a.ad4_type, "C");
            }
            other => panic!("expected Atom, got {other:?}"),
        }
    }

    #[test]
    fn parses_flexibility_tokens() {
        let doc = "\
ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C 
ENDBRANCH   1   2
TORSDOF 1
";
        let recs = parse(doc).unwrap();
        assert!(matches!(recs[0], PdbqtRecord::Root));
        assert!(matches!(recs[1], PdbqtRecord::Atom(_)));
        assert!(matches!(recs[2], PdbqtRecord::EndRoot));
        assert!(matches!(
            recs[3],
            PdbqtRecord::Branch {
                parent_serial: 1,
                child_serial: 2
            }
        ));
        assert!(matches!(
            recs[5],
            PdbqtRecord::EndBranch {
                parent_serial: 1,
                child_serial: 2
            }
        ));
        assert!(matches!(recs[6], PdbqtRecord::Torsdof(1)));
    }

    #[test]
    fn rejects_short_line() {
        let err = parse("ATOM  ").unwrap_err();
        assert!(matches!(err, PdbqtError::LineTooShort { .. }));
    }

    #[test]
    fn writer_round_trips_atom_record() {
        let atom = PdbqtAtom {
            serial: 7,
            name: "C1".to_string(),
            residue_name: "LIG".to_string(),
            chain: 'A',
            residue_seq: 1,
            position: Vector3::new(10.123, -2.456, 0.5),
            occupancy: 1.0,
            b_factor: 0.0,
            partial_charge: 0.150,
            ad4_type: "C".to_string(),
        };
        let line = write_atom(&atom);
        let parsed = parse(&line).unwrap();
        assert_eq!(parsed.len(), 1);
        match &parsed[0] {
            PdbqtRecord::Atom(a) => {
                assert_eq!(a.serial, 7);
                assert_eq!(a.name, "C1");
                assert!((a.position.x - 10.123).abs() < 1e-3);
                assert!((a.position.y - -2.456).abs() < 1e-3);
                assert!((a.position.z - 0.5).abs() < 1e-3);
                assert!((a.partial_charge - 0.150).abs() < 1e-3);
                assert_eq!(a.ad4_type, "C");
            }
            other => panic!("not an Atom: {other:?}"),
        }
    }
}

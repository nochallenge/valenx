//! MDL MOL (V2000) and SDF readers and writers.
//!
//! The MOL connection table is the workhorse interchange format: a
//! header block, a counts line, an atom block (element + 2D/3D
//! coordinates), a bond block (atom pair + bond type) and an `M  END`
//! terminator. An SD file (SDF) is a stream of MOL records, each
//! followed by `> <tag>` property blocks and a `$$$$` delimiter.
//!
//! [`read_mol`] / [`write_mol`] handle one V2000 record;
//! [`read_sdf`] / [`write_sdf`] handle a multi-record SD file,
//! round-tripping the per-molecule [`Molecule::properties`].
//!
//! **v1 simplifications:** only the V2000 dialect is parsed (V3000's
//! tagged block format is not); the atom-block stereo-parity and
//! charge columns are read into the model, but the assorted `M  CHG`
//! / `M  ISO` / `M  RAD` property lines are parsed for charge and
//! isotope while rarer property lines are skipped; query / R-group
//! features of the format are out of scope.

use crate::error::{CheminfError, Result};
use crate::molecule::{Atom, Bond, BondOrder, BondStereo, Molecule};

/// Parse a single V2000 MOL record into a [`Molecule`].
pub fn read_mol(text: &str) -> Result<Molecule> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 4 {
        return Err(CheminfError::parse(
            "mol",
            "truncated header — need ≥4 lines",
        ));
    }
    // line 0: molecule name; 1: program/timestamp; 2: comment; 3: counts
    let mut mol = Molecule::new();
    mol.name = lines[0].trim().to_string();

    let counts = lines[3];
    let n_atoms: usize = counts
        .get(0..3)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| CheminfError::parse("mol", "bad atom count in counts line"))?;
    let n_bonds: usize = counts
        .get(3..6)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| CheminfError::parse("mol", "bad bond count in counts line"))?;

    let atom_start = 4;
    if lines.len() < atom_start + n_atoms + n_bonds {
        return Err(CheminfError::parse(
            "mol",
            "file shorter than declared atom + bond counts",
        ));
    }

    let mut coords: Vec<[f64; 3]> = Vec::with_capacity(n_atoms);
    for i in 0..n_atoms {
        let line = lines[atom_start + i];
        // V2000 atom line is fixed-width: x(10) y(10) z(10) ' ' sym(3)...
        let x: f64 = parse_f(line, 0, 10)?;
        let y: f64 = parse_f(line, 10, 20)?;
        let z: f64 = parse_f(line, 20, 30)?;
        let sym = line
            .get(31..34)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .ok_or_else(|| CheminfError::parse("mol", format!("missing element on atom {i}")))?;
        let element = crate::element::by_symbol(sym).ok_or_else(|| {
            CheminfError::parse("mol", format!("unknown element `{sym}` on atom {i}"))
        })?;
        let mut atom = Atom::new(element.number);
        // legacy charge field (column 36-39): 0=none,1..7 mapped charges
        if let Some(c) = line.get(36..39).and_then(|s| s.trim().parse::<i32>().ok()) {
            atom.formal_charge = legacy_charge(c);
        }
        mol.add_atom(atom);
        coords.push([x, y, z]);
    }

    let bond_start = atom_start + n_atoms;
    for i in 0..n_bonds {
        let line = lines[bond_start + i];
        let a: usize = parse_u(line, 0, 3)?;
        let b: usize = parse_u(line, 3, 6)?;
        let ty: u32 = parse_u(line, 6, 9)? as u32;
        if a == 0 || b == 0 || a > n_atoms || b > n_atoms {
            return Err(CheminfError::parse(
                "mol",
                format!("bond {i} references an out-of-range atom"),
            ));
        }
        let order = match ty {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        // bond stereo column (10-12): 1 = wedge up, 6 = wedge down
        let stereo = match parse_u(line, 9, 12).unwrap_or(0) {
            1 => BondStereo::Up,
            6 => BondStereo::Down,
            _ => BondStereo::None,
        };
        mol.add_bond(Bond {
            a: a - 1,
            b: b - 1,
            order,
            aromatic: order == BondOrder::Aromatic,
            stereo,
        });
    }

    // property lines: M  CHG and M  ISO override the legacy fields
    for line in &lines[bond_start + n_bonds..] {
        if line.starts_with("M  END") {
            break;
        }
        if let Some(rest) = line.strip_prefix("M  CHG") {
            apply_chg_iso(&mut mol, rest, true);
        } else if let Some(rest) = line.strip_prefix("M  ISO") {
            apply_chg_iso(&mut mol, rest, false);
        }
    }

    // only keep coordinates if at least one is non-zero (a 2D/3D block)
    if coords.iter().any(|c| c.iter().any(|v| v.abs() > 1e-9)) {
        let any_z = coords.iter().any(|c| c[2].abs() > 1e-6);
        mol.coords = coords;
        mol.coords_3d = any_z;
    }

    crate::perceive::hydrogen::recompute_implicit_hydrogens(&mut mol);
    mol.validate()?;
    Ok(mol)
}

/// Apply an `M  CHG` (charge) or `M  ISO` (isotope) property line.
fn apply_chg_iso(mol: &mut Molecule, rest: &str, is_charge: bool) {
    let nums: Vec<i32> = rest
        .split_whitespace()
        .filter_map(|t| t.parse().ok())
        .collect();
    if nums.is_empty() {
        return;
    }
    let pairs = nums[0] as usize;
    for k in 0..pairs {
        let idx = 1 + k * 2;
        if idx + 1 >= nums.len() {
            break;
        }
        let atom = (nums[idx] - 1) as usize;
        let value = nums[idx + 1];
        if let Some(a) = mol.atoms.get_mut(atom) {
            if is_charge {
                a.formal_charge = value.clamp(-9, 9) as i8;
            } else {
                a.isotope = u16::try_from(value).ok();
            }
        }
    }
}

/// Translate the legacy single-digit charge code into a formal charge.
fn legacy_charge(code: i32) -> i8 {
    match code {
        1 => 3,
        2 => 2,
        3 => 1,
        5 => -1,
        6 => -2,
        7 => -3,
        _ => 0,
    }
}

fn parse_f(line: &str, lo: usize, hi: usize) -> Result<f64> {
    line.get(lo..hi)
        .map(|s| s.trim())
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| CheminfError::parse("mol", format!("bad number in columns {lo}..{hi}")))
}

fn parse_u(line: &str, lo: usize, hi: usize) -> Result<usize> {
    line.get(lo..hi)
        .map(|s| s.trim())
        .and_then(|s| {
            if s.is_empty() {
                Some(0)
            } else {
                s.parse().ok()
            }
        })
        .ok_or_else(|| CheminfError::parse("mol", format!("bad integer in columns {lo}..{hi}")))
}

/// Write a [`Molecule`] as a V2000 MOL record (no trailing `$$$$`).
pub fn write_mol(mol: &Molecule) -> String {
    let mut out = String::new();
    // header block
    out.push_str(&mol.name);
    out.push('\n');
    out.push_str("  valenx-cheminf\n");
    out.push('\n');
    // counts line
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));
    // atom block
    for (i, a) in mol.atoms.iter().enumerate() {
        let c = mol.coords.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
        let charge_code = match a.formal_charge {
            3 => 1,
            2 => 2,
            1 => 3,
            -1 => 5,
            -2 => 6,
            -3 => 7,
            _ => 0,
        };
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}  0  0  0  0  0  0  0  0  0  0\n",
            c[0],
            c[1],
            c[2],
            a.symbol(),
            0, // mass-difference field, superseded by M ISO
            charge_code,
        ));
    }
    // bond block
    for b in &mol.bonds {
        let ty = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 1,
            BondOrder::Aromatic => 4,
        };
        let st = match b.stereo {
            BondStereo::Up => 1,
            BondStereo::Down => 6,
            _ => 0,
        };
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0\n",
            b.a + 1,
            b.b + 1,
            ty,
            st
        ));
    }
    // M  CHG / M  ISO for atoms that carry them (more reliable than the
    // legacy column for |charge| > 3)
    let charged: Vec<(usize, i8)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.formal_charge != 0)
        .map(|(i, a)| (i, a.formal_charge))
        .collect();
    if !charged.is_empty() {
        for chunk in charged.chunks(8) {
            out.push_str(&format!("M  CHG{:>3}", chunk.len()));
            for (i, c) in chunk {
                out.push_str(&format!("{:>4}{:>4}", i + 1, c));
            }
            out.push('\n');
        }
    }
    let isos: Vec<(usize, u16)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| a.isotope.map(|iso| (i, iso)))
        .collect();
    if !isos.is_empty() {
        for chunk in isos.chunks(8) {
            out.push_str(&format!("M  ISO{:>3}", chunk.len()));
            for (i, iso) in chunk {
                out.push_str(&format!("{:>4}{:>4}", i + 1, iso));
            }
            out.push('\n');
        }
    }
    out.push_str("M  END\n");
    out
}

/// Parse a multi-record SD file. Each record's `> <tag>` blocks become
/// entries in that molecule's [`Molecule::properties`].
pub fn read_sdf(text: &str) -> Result<Vec<Molecule>> {
    let mut molecules = Vec::new();
    for (idx, record) in text.split("$$$$").enumerate() {
        // Strip ONLY the single newline that terminates the preceding
        // `$$$$` line — not every leading blank line. A V2000 MOL
        // record's first line is the molecule *name*, which is
        // legitimately blank for an unnamed molecule; trimming all
        // leading newlines would delete that blank name line and shift
        // the counts line, so the record would no longer parse. The
        // first record has no preceding `$$$$`, so nothing is stripped.
        let record = if idx == 0 {
            record
        } else {
            record
                .strip_prefix("\r\n")
                .or_else(|| record.strip_prefix('\n'))
                .unwrap_or(record)
        };
        let trimmed = record.trim_end_matches(['\r', '\n']);
        if trimmed.trim().is_empty() {
            continue;
        }
        // the MOL block ends at `M  END`; properties follow
        let end_pos = trimmed.find("M  END");
        let (mol_block, prop_block) = match end_pos {
            Some(p) => (&trimmed[..p + 6], &trimmed[p + 6..]),
            None => (trimmed, ""),
        };
        let mut mol = read_mol(mol_block)
            .map_err(|e| CheminfError::parse("sdf", format!("record {} — {e}", idx + 1)))?;
        parse_sdf_properties(prop_block, &mut mol);
        molecules.push(mol);
    }
    if molecules.is_empty() {
        return Err(CheminfError::parse("sdf", "no records found"));
    }
    Ok(molecules)
}

/// Parse the `> <tag>\nvalue` property lines after a MOL block.
fn parse_sdf_properties(block: &str, mol: &mut Molecule) {
    let lines: Vec<&str> = block.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i].trim();
        if let Some(tag) = extract_tag(line) {
            // value lines run until a blank line
            let mut value = String::new();
            i += 1;
            while i < lines.len() && !lines[i].trim().is_empty() {
                if !value.is_empty() {
                    value.push('\n');
                }
                value.push_str(lines[i]);
                i += 1;
            }
            mol.properties.push((tag, value));
        }
        i += 1;
    }
}

/// Extract the tag from a `> <Name>` or `>  <Name>` property header.
fn extract_tag(line: &str) -> Option<String> {
    if !line.starts_with('>') {
        return None;
    }
    let lt = line.find('<')?;
    let gt = line[lt + 1..].find('>')?;
    Some(line[lt + 1..lt + 1 + gt].to_string())
}

/// Write a slice of molecules as an SD file (`$$$$`-delimited records).
pub fn write_sdf(molecules: &[Molecule]) -> String {
    let mut out = String::new();
    for mol in molecules {
        out.push_str(&write_mol(mol));
        for (key, value) in &mol.properties {
            out.push_str(&format!(">  <{key}>\n{value}\n\n"));
        }
        out.push_str("$$$$\n");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn write_then_read_round_trips() {
        let mut m = parse_smiles("CCO").unwrap();
        m.name = "ethanol".to_string();
        // give it coordinates so the atom block is meaningful
        m.coords = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.2, 1.0, 0.0]];
        let text = write_mol(&m);
        let back = read_mol(&text).unwrap();
        assert_eq!(back.atom_count(), 3);
        assert_eq!(back.bond_count(), 2);
        assert_eq!(back.name, "ethanol");
        assert_eq!(back.atoms[2].atomic_number, 8);
    }

    #[test]
    fn charge_round_trips() {
        let m = parse_smiles("[NH4+]").unwrap();
        let text = write_mol(&m);
        let back = read_mol(&text).unwrap();
        assert_eq!(back.atoms[0].formal_charge, 1);
    }

    #[test]
    fn sdf_multi_record_with_properties() {
        let mut a = parse_smiles("CCO").unwrap();
        a.name = "ethanol".to_string();
        a.properties.push(("ID".to_string(), "1".to_string()));
        let mut b = parse_smiles("c1ccccc1").unwrap();
        b.name = "benzene".to_string();
        b.properties.push(("ID".to_string(), "2".to_string()));

        let sdf = write_sdf(&[a, b]);
        assert_eq!(sdf.matches("$$$$").count(), 2);

        let records = read_sdf(&sdf).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "ethanol");
        assert_eq!(
            records[0].properties.iter().find(|(k, _)| k == "ID"),
            Some(&("ID".to_string(), "1".to_string()))
        );
        assert_eq!(records[1].atom_count(), 6);
    }

    #[test]
    fn rejects_truncated() {
        assert!(read_mol("too\nshort\n").is_err());
        assert!(read_sdf("").is_err());
    }

    #[test]
    fn isotope_round_trips() {
        let m = parse_smiles("[13CH4]").unwrap();
        let text = write_mol(&m);
        assert!(text.contains("M  ISO"));
        let back = read_mol(&text).unwrap();
        assert_eq!(back.atoms[0].isotope, Some(13));
    }
}

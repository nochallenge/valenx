//! PDB-format reader and writer.
//!
//! Handles the column-oriented PDB record set most analysis needs:
//! `HEADER` / `TITLE`, `ATOM` / `HETATM`, `MODEL` / `ENDMDL`,
//! `TER`, `SEQRES`, `HELIX`, `SHEET`, `SSBOND`, `CONECT`, and the
//! `REMARK 350` `BIOMT` operators.
//!
//! ## Scope of this v1
//!
//! The PDB format is fixed-column. The reader slices by column index
//! (PDB v3.3 spec) and tolerates short lines by clamping. It does
//! **not** interpret `ANISOU`, `MODRES`, `LINK`,
//! `SCALE` / `ORIGX` / `CRYST1` transforms, or the full `REMARK`
//! taxonomy beyond `350` — those records are skipped without error.
//! The writer emits a self-consistent subset (`HEADER`, `TITLE`,
//! `SSBOND`, `ATOM` / `HETATM`, `TER`, `MODEL` / `ENDMDL`, `CONECT`,
//! `END`).

use crate::error::{BiostructError, Result};
use crate::structure::{
    Atom, Chain, Disulfide, Model, Residue, Structure, SymmetryOperator,
};
use nalgebra::Point3;

/// Parse a PDB-format string into a [`Structure`].
///
/// `id` becomes the structure id when the file carries no usable
/// `HEADER` id.
pub fn read_pdb(text: &str, id: &str) -> Result<Structure> {
    let mut structure = Structure::new(id);
    structure.models.clear();

    // The "current" model being filled. A file with no explicit
    // MODEL record gets an implicit model 1.
    let mut current_model = Model::new(1);
    let mut seen_model_record = false;
    let mut seqres: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    // REMARK 350 BIOMT accumulator: serial -> 3 rows being filled.
    let mut biomt: std::collections::HashMap<i32, [[f64; 4]; 3]> =
        std::collections::HashMap::new();

    for (lineno0, raw) in text.lines().enumerate() {
        let lineno = lineno0 + 1;
        let record = record_name(raw);
        match record.as_str() {
            "HEADER" => {
                // id_code occupies columns 63-66.
                let id_field = slice(raw, 62, 66).trim();
                if !id_field.is_empty() && structure.id == id {
                    structure.id = id_field.to_string();
                }
            }
            "TITLE" => {
                let chunk = slice(raw, 10, 80).trim_end();
                if !structure.title.is_empty() {
                    structure.title.push(' ');
                }
                structure.title.push_str(chunk.trim());
            }
            "SEQRES" => {
                let chain = slice(raw, 11, 12).trim().to_string();
                let body = slice(raw, 19, 70);
                let entry = seqres.entry(chain).or_default();
                for tlc in body.split_whitespace() {
                    entry.push(crate::structure::residue_one_letter(tlc));
                }
            }
            "HELIX" => {
                // init chain col 20, init seq 22-25; end chain 32,
                // end seq 34-37.
                let chain = slice(raw, 19, 20).trim().to_string();
                if let (Ok(s), Ok(e)) = (
                    slice(raw, 21, 25).trim().parse::<i32>(),
                    slice(raw, 33, 37).trim().parse::<i32>(),
                ) {
                    structure.helix_records.push((chain, s, e));
                }
            }
            "SHEET" => {
                // init chain col 22, init seq 23-26; end chain 33,
                // end seq 34-37.
                let chain = slice(raw, 21, 22).trim().to_string();
                if let (Ok(s), Ok(e)) = (
                    slice(raw, 22, 26).trim().parse::<i32>(),
                    slice(raw, 33, 37).trim().parse::<i32>(),
                ) {
                    structure.sheet_records.push((chain, s, e));
                }
            }
            "SSBOND" => {
                if let Some(ss) = parse_ssbond(raw) {
                    structure.disulfides.push(ss);
                }
            }
            "CONECT" => {
                if let Some(rec) = parse_conect(raw) {
                    structure.conect.push(rec);
                }
            }
            "REMARK" => {
                parse_remark_350_biomt(raw, &mut biomt);
            }
            "MODEL" => {
                if seen_model_record && !current_model.chains.is_empty() {
                    structure.models.push(current_model);
                }
                let serial = slice(raw, 10, 14).trim().parse::<i32>().unwrap_or(
                    structure.models.len() as i32 + 1,
                );
                current_model = Model::new(serial);
                seen_model_record = true;
            }
            "ENDMDL" => {
                structure.models.push(std::mem::replace(
                    &mut current_model,
                    Model::new(0),
                ));
            }
            "ATOM" | "HETATM" => {
                let (chain_id, residue, atom) = parse_atom_record(raw, lineno)?;
                push_atom(&mut current_model, &chain_id, residue, atom);
            }
            "TER" => { /* chain terminator — structural Vec order already segments */ }
            _ => { /* ANISOU, SCALE, LINK, MODRES, … — skipped */ }
        }
    }

    // Flush a trailing implicit / unterminated model.
    if !current_model.chains.is_empty() {
        structure.models.push(current_model);
    }
    if structure.models.is_empty() {
        return Err(BiostructError::parse(
            "pdb",
            0,
            "no ATOM/HETATM records found",
        ));
    }

    // Attach SEQRES sequences to chains.
    for model in &mut structure.models {
        for chain in &mut model.chains {
            if let Some(seq) = seqres.get(&chain.id) {
                chain.seqres = Some(seq.clone());
            }
        }
    }

    // Finalise BIOMT operators (only rows with all 3 filled).
    let mut serials: Vec<i32> = biomt.keys().copied().collect();
    serials.sort_unstable();
    for serial in serials {
        let rows = &biomt[&serial];
        structure.assembly_operators.push(SymmetryOperator {
            serial,
            rotation: [
                [rows[0][0], rows[0][1], rows[0][2]],
                [rows[1][0], rows[1][1], rows[1][2]],
                [rows[2][0], rows[2][1], rows[2][2]],
            ],
            translation: [rows[0][3], rows[1][3], rows[2][3]],
        });
    }

    Ok(structure)
}

/// Serialise a [`Structure`] back into PDB-format text.
///
/// Emits `HEADER`, `TITLE`, `SSBOND`, then per model the `ATOM` /
/// `HETATM` records (wrapped in `MODEL` / `ENDMDL` when there is more
/// than one model), a `TER` after each polymer chain, the `CONECT`
/// records, and a final `END`.
pub fn write_pdb(structure: &Structure) -> String {
    let mut out = String::new();
    // HEADER: classification blank, date blank, id in 63-66.
    out.push_str(&format!(
        "HEADER    {:<40}{:>4}{:<14}\n",
        "", structure.id, ""
    ));
    if !structure.title.is_empty() {
        for (i, chunk) in wrap_title(&structure.title).into_iter().enumerate() {
            if i == 0 {
                out.push_str(&format!("TITLE     {chunk}\n"));
            } else {
                out.push_str(&format!("TITLE   {:>2} {chunk}\n", i + 1));
            }
        }
    }
    // SSBOND disulfide records.
    for (i, ss) in structure.disulfides.iter().enumerate() {
        let (ca, sa, ia) = &ss.partner_a;
        let (cb, sb, ib) = &ss.partner_b;
        out.push_str(&format!(
            "SSBOND {serial:>3} CYS {ca:>1} {sa:>4}{ia:>1}   CYS {cb:>1} {sb:>4}{ib:>1}\n",
            serial = (i + 1).min(999),
            ca = ca.chars().next().unwrap_or('A'),
            sa = sa,
            ia = ia,
            cb = cb.chars().next().unwrap_or('A'),
            sb = sb,
            ib = ib,
        ));
    }

    let multi = structure.models.len() > 1;
    let mut serial: i32 = 1;
    for model in &structure.models {
        if multi {
            out.push_str(&format!("MODEL     {:>4}\n", model.serial));
        }
        for chain in &model.chains {
            let mut last_polymer: Option<&Residue> = None;
            for residue in &chain.residues {
                for atom in &residue.atoms {
                    out.push_str(&format_atom_record(
                        serial, atom, residue, &chain.id,
                    ));
                    out.push('\n');
                    serial += 1;
                }
                if residue.is_amino_acid() || residue.is_nucleotide() {
                    last_polymer = Some(residue);
                }
            }
            if let Some(r) = last_polymer {
                out.push_str(&format!(
                    "TER   {:>5}      {:>3} {:>1}{:>4}{}\n",
                    serial,
                    r.name,
                    chain.id,
                    r.seq_num,
                    r.ins_code,
                ));
                serial += 1;
            }
        }
        if multi {
            out.push_str("ENDMDL\n");
        }
    }
    // CONECT records — up to four partner serials per line.
    for (base, partners) in &structure.conect {
        for chunk in partners.chunks(4) {
            out.push_str(&format!("CONECT{:>5}", base.min(&99999)));
            for p in chunk {
                out.push_str(&format!("{:>5}", p.min(&99999)));
            }
            out.push('\n');
        }
    }
    out.push_str("END\n");
    out
}

// --- record-level helpers --------------------------------------------

/// Extract the 6-character record name, padded / trimmed.
fn record_name(line: &str) -> String {
    slice(line, 0, 6).trim().to_ascii_uppercase()
}

/// Slice a line by 0-based half-open column range `[start, end)`,
/// clamping to the line length. PDB columns in the spec are 1-based
/// and inclusive — convert with `start = col1-1`, `end = col2`.
fn slice(line: &str, start: usize, end: usize) -> &str {
    let bytes = line.as_bytes();
    let s = start.min(bytes.len());
    let e = end.min(bytes.len());
    if s >= e {
        return "";
    }
    // PDB files are ASCII; slicing on byte indices is safe here, but
    // guard against a stray multibyte char by falling back to "".
    line.get(s..e).unwrap_or("")
}

/// Parse one `ATOM` / `HETATM` record into `(chain_id, residue
/// shell, atom)`. The residue shell carries the residue identity but
/// no atoms; the caller merges atoms into the matching residue.
fn parse_atom_record(line: &str, lineno: usize) -> Result<(String, Residue, Atom)> {
    let hetatm = record_name(line) == "HETATM";
    let serial = slice(line, 6, 11).trim().parse::<i32>().unwrap_or(0);
    // Keep the raw 4-column name field: its leading-space alignment is
    // the cue that distinguishes a single-letter element (" CA " = Cα)
    // from a two-letter one ("CA  " = calcium). The stored atom name is
    // the trimmed form.
    let name_field = slice(line, 12, 16).to_string();
    let name = name_field.trim().to_string();
    let alt_loc = slice(line, 16, 17).chars().next().unwrap_or(' ');
    let res_name = slice(line, 17, 20).trim().to_ascii_uppercase();
    let chain_id = {
        let c = slice(line, 21, 22).trim();
        if c.is_empty() { "A".to_string() } else { c.to_string() }
    };
    let seq_num = slice(line, 22, 26)
        .trim()
        .parse::<i32>()
        .map_err(|_| {
            BiostructError::parse("pdb", lineno, "non-integer residue sequence number")
        })?;
    let ins_code = slice(line, 26, 27).chars().next().unwrap_or(' ');
    let x = parse_coord(line, 30, 38, lineno, "x")?;
    let y = parse_coord(line, 38, 46, lineno, "y")?;
    let z = parse_coord(line, 46, 54, lineno, "z")?;
    let occupancy = slice(line, 54, 60).trim().parse::<f64>().unwrap_or(1.0);
    let b_factor = slice(line, 60, 66).trim().parse::<f64>().unwrap_or(0.0);
    let element_field = slice(line, 76, 78).trim().to_ascii_uppercase();
    let element = if element_field.is_empty() {
        guess_element(&name_field)
    } else {
        element_field
    };
    let charge = parse_charge(slice(line, 78, 80));

    let mut residue = Residue::new(res_name, seq_num);
    residue.ins_code = ins_code;
    residue.hetatm = hetatm;

    let atom = Atom {
        serial,
        name,
        alt_loc,
        element,
        coord: Point3::new(x, y, z),
        occupancy,
        b_factor,
        charge,
    };
    Ok((chain_id, residue, atom))
}

/// Parse a fixed-width coordinate field, erroring on a non-number.
fn parse_coord(line: &str, s: usize, e: usize, lineno: usize, axis: &str) -> Result<f64> {
    slice(line, s, e)
        .trim()
        .parse::<f64>()
        .map_err(|_| BiostructError::parse("pdb", lineno, format!("bad {axis} coordinate")))
}

/// Parse the 2-character PDB charge field (`"2+"`, `"1-"`, `""`).
fn parse_charge(field: &str) -> i32 {
    let f = field.trim();
    if f.len() < 2 {
        return 0;
    }
    let digit = f.chars().next().and_then(|c| c.to_digit(10));
    let sign = f.chars().nth(1);
    match (digit, sign) {
        (Some(d), Some('+')) => d as i32,
        (Some(d), Some('-')) => -(d as i32),
        _ => 0,
    }
}

/// Guess an element symbol from an atom name when the element column
/// is blank.
///
/// The PDB atom-name field (columns 13-16) carries the disambiguation
/// that a context-free string does not: a **single-letter** element is
/// right-justified so that the element letter sits in column 14,
/// leaving column 13 blank (`" CA "` is carbon-α). A genuine
/// **two-letter** element fills columns 13-14 (`"CA  "` is calcium,
/// `"FE  "` is iron). So when `name` is passed with that leading-space
/// alignment preserved, this distinguishes Cα from a calcium ion
/// correctly.
///
/// If `name` arrives already trimmed (no column information — e.g. from
/// mmCIF), the leading-space cue is gone; the heuristic then takes the
/// leading alphabetic run, treating a polymer-style name (a biopolymer
/// element C/N/O/S/P/H followed by a Greek-position letter) as that
/// single-letter element, since Cα-class atoms vastly outnumber the
/// ambiguous metals in real structures.
pub fn guess_element(name: &str) -> String {
    // A leading space in the *untrimmed* name is the column-13 cue for
    // a single-letter element.
    let single_letter_aligned = name.starts_with(' ') && !name.trim().is_empty();

    let trimmed = name.trim();
    let alpha: String = trimmed.chars().filter(|c| c.is_ascii_alphabetic()).collect();
    if alpha.is_empty() {
        return "C".to_string();
    }
    // Hydrogen atom names often start with a digit then H.
    if trimmed.chars().next().map(|c| c.is_ascii_digit()) == Some(true)
        && alpha.starts_with('H')
    {
        return "H".to_string();
    }
    let upper = alpha.to_ascii_uppercase();

    // Column-aligned single-letter element: the leading space settles
    // it — take the first letter only.
    if single_letter_aligned {
        return upper[0..1].to_string();
    }

    // Recognised two-letter elements that appear in PDB files.
    const TWO: &[&str] = &[
        "FE", "ZN", "MG", "MN", "CA", "NA", "CL", "CU", "NI", "CO", "SE", "BR", "CD", "HG",
        "PT", "AU", "AG", "MO", "AS",
    ];
    if upper.len() >= 2 && TWO.contains(&&upper[0..2]) {
        // No column cue. A polymer atom name is a biopolymer element
        // (C/N/O/S/P/H) followed by a Greek-position letter
        // (α/β/γ/δ/ε/ζ/η) — read those as the single-letter element.
        let first = &upper[0..1];
        let second = upper.as_bytes()[1];
        let biopolymer = matches!(first, "C" | "N" | "O" | "S" | "P" | "H");
        let greek_position = matches!(second, b'A' | b'B' | b'G' | b'D' | b'E' | b'Z' | b'H');
        if biopolymer && greek_position {
            return first.to_string();
        }
        return upper[0..2].to_string();
    }
    upper[0..1].to_string()
}

/// Merge a parsed atom into `model`, creating the chain / residue
/// when first seen. Residues sharing `(seq_num, ins_code, name)`
/// within a chain accumulate atoms.
fn push_atom(model: &mut Model, chain_id: &str, residue: Residue, atom: Atom) {
    let chain = match model.chains.iter_mut().find(|c| c.id == chain_id) {
        Some(c) => c,
        None => {
            model.chains.push(Chain::new(chain_id));
            model.chains.last_mut().unwrap()
        }
    };
    // A residue matches the most-recently-added one with the same id;
    // PDB files keep a residue's atoms contiguous, so checking the
    // tail is correct and O(1).
    let matches_tail = chain
        .residues
        .last()
        .map(|r| {
            r.seq_num == residue.seq_num
                && r.ins_code == residue.ins_code
                && r.name == residue.name
        })
        .unwrap_or(false);
    if matches_tail {
        chain.residues.last_mut().unwrap().atoms.push(atom);
    } else {
        let mut r = residue;
        r.atoms.push(atom);
        chain.residues.push(r);
    }
}

/// Accumulate a `REMARK 350   BIOMT{1,2,3}` line into `biomt`.
fn parse_remark_350_biomt(line: &str, biomt: &mut std::collections::HashMap<i32, [[f64; 4]; 3]>) {
    // REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
    let body = slice(line, 7, 80);
    let toks: Vec<&str> = body.split_whitespace().collect();
    if toks.len() < 7 || toks[0] != "350" {
        return;
    }
    let row = match toks[1] {
        "BIOMT1" => 0usize,
        "BIOMT2" => 1usize,
        "BIOMT3" => 2usize,
        _ => return,
    };
    let serial = match toks[2].parse::<i32>() {
        Ok(s) => s,
        Err(_) => return,
    };
    let vals: Vec<f64> = toks[3..7].iter().filter_map(|t| t.parse().ok()).collect();
    if vals.len() != 4 {
        return;
    }
    let entry = biomt.entry(serial).or_insert([[0.0; 4]; 3]);
    entry[row] = [vals[0], vals[1], vals[2], vals[3]];
}

/// Parse a `SSBOND` disulfide-bond record.
///
/// PDB column layout: partner-1 chain at col 16, seq 18-21, iCode 22;
/// partner-2 chain at col 30, seq 32-35, iCode 36.
fn parse_ssbond(line: &str) -> Option<Disulfide> {
    let chain_a = slice(line, 15, 16).trim();
    let seq_a = slice(line, 17, 21).trim().parse::<i32>().ok()?;
    let ins_a = slice(line, 21, 22).chars().next().unwrap_or(' ');
    let chain_b = slice(line, 29, 30).trim();
    let seq_b = slice(line, 31, 35).trim().parse::<i32>().ok()?;
    let ins_b = slice(line, 35, 36).chars().next().unwrap_or(' ');
    Some(Disulfide {
        partner_a: (
            if chain_a.is_empty() { "A" } else { chain_a }.to_string(),
            seq_a,
            ins_a,
        ),
        partner_b: (
            if chain_b.is_empty() { "A" } else { chain_b }.to_string(),
            seq_b,
            ins_b,
        ),
    })
}

/// Parse a `CONECT` record into `(atom_serial, [partner serials])`.
///
/// PDB `CONECT` carries the base atom serial in columns 7-11 and up to
/// four bonded-partner serials in the following 5-column fields. This
/// reader accepts any number of whitespace-separated serials so it
/// also handles non-strict files.
fn parse_conect(line: &str) -> Option<(i32, Vec<i32>)> {
    let mut serials = slice(line, 6, 80)
        .split_whitespace()
        .filter_map(|t| t.parse::<i32>().ok());
    let base = serials.next()?;
    let partners: Vec<i32> = serials.collect();
    if partners.is_empty() {
        return None;
    }
    Some((base, partners))
}

/// Format a single `ATOM` / `HETATM` record (no trailing newline).
fn format_atom_record(serial: i32, atom: &Atom, residue: &Residue, chain_id: &str) -> String {
    let record = if residue.hetatm { "HETATM" } else { "ATOM  " };
    // Atom name alignment: 4-char-or-longer names start at col 13,
    // shorter element-1 names are offset by one for the classic
    // " CA " look.
    let name = if atom.name.len() >= 4 || atom.element.len() == 2 {
        format!("{:<4}", atom.name)
    } else {
        format!(" {:<3}", atom.name)
    };
    let charge_field = match atom.charge {
        0 => "  ".to_string(),
        c if c > 0 => format!("{}+", c.min(9)),
        c => format!("{}-", (-c).min(9)),
    };
    format!(
        "{record}{serial:>5} {name}{alt:>1}{resn:>3} {chain:>1}{seq:>4}{ins:>1}   \
         {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{b:>6.2}          {elem:>2}{charge}",
        record = record,
        serial = serial.min(99999),
        name = name,
        alt = atom.alt_loc,
        resn = residue.name,
        chain = chain_id.chars().next().unwrap_or('A'),
        seq = residue.seq_num,
        ins = residue.ins_code,
        x = atom.coord.x,
        y = atom.coord.y,
        z = atom.coord.z,
        occ = atom.occupancy,
        b = atom.b_factor,
        elem = atom.element,
        charge = charge_field,
    )
}

/// Wrap a long title into <=70-char `TITLE`-continuation chunks.
fn wrap_title(title: &str) -> Vec<String> {
    let mut chunks = Vec::new();
    let mut current = String::new();
    for word in title.split_whitespace() {
        if current.len() + word.len() + 1 > 69 && !current.is_empty() {
            chunks.push(std::mem::take(&mut current));
        }
        if !current.is_empty() {
            current.push(' ');
        }
        current.push_str(word);
    }
    if !current.is_empty() {
        chunks.push(current);
    }
    if chunks.is_empty() {
        chunks.push(String::new());
    }
    chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINI_PDB: &str = "\
HEADER    HYDROLASE                               01-JAN-00   1ABC
TITLE     A SMALL TEST STRUCTURE
SEQRES   1 A    3  ALA GLY SER
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00 21.00           C
ATOM      3  C   ALA A   1      13.085   6.564  -5.180  1.00 19.50           C
ATOM      4  O   ALA A   1      13.339   7.742  -5.430  1.00 18.00           O
ATOM      5  CB  ALA A   1      10.768   6.882  -4.205  0.50 22.00           C
ATOM      6  N   GLY A   2      13.965   5.628  -4.926  1.00 17.00           N
ATOM      7  CA  GLY A   2      15.398   5.927  -4.911  1.00 16.00           C
ATOM      8  C   GLY A   2      16.180   4.741  -4.378  1.00 15.50           C
HETATM    9 ZN    ZN A 101      20.000  20.000  20.000  1.00 30.00          ZN2+
END
";

    #[test]
    fn reads_atoms_and_residues() {
        let s = read_pdb(MINI_PDB, "fallback").unwrap();
        assert_eq!(s.id, "1ABC");
        assert_eq!(s.title, "A SMALL TEST STRUCTURE");
        let m = s.first_model();
        assert_eq!(m.chains.len(), 1);
        let chain = m.chain("A").unwrap();
        // ALA, GLY, ZN
        assert_eq!(chain.residues.len(), 3);
        assert_eq!(chain.residues[0].atoms.len(), 5);
        assert_eq!(chain.residues[0].ca().unwrap().b_factor, 21.0);
    }

    #[test]
    fn reads_hetatm_and_charge() {
        let s = read_pdb(MINI_PDB, "x").unwrap();
        let zn = s.first_model().chain("A").unwrap().residue(101, ' ').unwrap();
        assert!(zn.hetatm);
        assert_eq!(zn.atoms[0].element, "ZN");
        assert_eq!(zn.atoms[0].charge, 2);
    }

    #[test]
    fn reads_seqres() {
        let s = read_pdb(MINI_PDB, "x").unwrap();
        assert_eq!(
            s.first_model().chain("A").unwrap().seqres.as_deref(),
            Some("AGS")
        );
    }

    #[test]
    fn element_guessing() {
        assert_eq!(guess_element("CA"), "C"); // calpha, single-letter C
        assert_eq!(guess_element("FE"), "FE");
        assert_eq!(guess_element("1HB"), "H");
        assert_eq!(guess_element("OD1"), "O");
    }

    #[test]
    fn round_trips_coordinates() {
        let s = read_pdb(MINI_PDB, "x").unwrap();
        let text = write_pdb(&s);
        let s2 = read_pdb(&text, "x").unwrap();
        let a1 = s.first_model().chain("A").unwrap().residues[0].atoms[0].coord;
        let a2 = s2.first_model().chain("A").unwrap().residues[0].atoms[0].coord;
        assert!((a1 - a2).norm() < 1e-3, "coords drifted: {a1} vs {a2}");
        assert_eq!(s2.atom_count(), s.atom_count());
    }

    #[test]
    fn multi_model_round_trip() {
        let two = "\
MODEL        1
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ENDMDL
END
";
        let s = read_pdb(two, "x").unwrap();
        assert_eq!(s.models.len(), 2);
        let text = write_pdb(&s);
        assert!(text.contains("MODEL"));
        let s2 = read_pdb(&text, "x").unwrap();
        assert_eq!(s2.models.len(), 2);
    }

    #[test]
    fn empty_file_errors() {
        assert!(read_pdb("REMARK just a remark\n", "x").is_err());
    }

    #[test]
    fn parses_and_writes_ssbond() {
        let pdb = "\
SSBOND   1 CYS A    6    CYS A  127
ATOM      1  SG  CYS A   6       0.000   0.000   0.000  1.00  0.00           S
ATOM      2  SG  CYS A 127       2.050   0.000   0.000  1.00  0.00           S
END
";
        let s = read_pdb(pdb, "x").unwrap();
        assert_eq!(s.disulfides.len(), 1);
        assert_eq!(s.disulfides[0].partner_a, ("A".to_string(), 6, ' '));
        assert_eq!(s.disulfides[0].partner_b, ("A".to_string(), 127, ' '));
        // Round-trip: the writer emits SSBOND, the reader parses it back.
        let text = write_pdb(&s);
        assert!(text.contains("SSBOND"));
        let s2 = read_pdb(&text, "x").unwrap();
        assert_eq!(s2.disulfides, s.disulfides);
    }

    #[test]
    fn parses_and_writes_conect() {
        let pdb = "\
HETATM    1  C1  LIG A 200       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  LIG A 200       1.300   0.000   0.000  1.00  0.00           O
CONECT    1    2
CONECT    2    1
END
";
        let s = read_pdb(pdb, "x").unwrap();
        assert_eq!(s.conect.len(), 2);
        assert_eq!(s.conect[0], (1, vec![2]));
        let text = write_pdb(&s);
        assert!(text.contains("CONECT"));
        let s2 = read_pdb(&text, "x").unwrap();
        assert_eq!(s2.conect, s.conect);
    }

    #[test]
    fn parses_biomt_operators() {
        let pdb = "\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        5.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
END
";
        let s = read_pdb(pdb, "x").unwrap();
        assert_eq!(s.assembly_operators.len(), 2);
        assert!(s.assembly_operators[0].is_identity());
        assert_eq!(s.assembly_operators[1].translation[0], 5.0);
    }
}

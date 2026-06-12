//! mmCIF-format reader and writer.
//!
//! mmCIF (the PDBx/mmCIF dictionary) is the modern archival format.
//! Unlike PDB it is **token-based**, not column-based: data lives in
//! `loop_` tables keyed by named tags.
//!
//! This v1 reads the categories an analysis pipeline needs:
//!
//! - `_atom_site` — the coordinate loop (the only mandatory one).
//! - `_struct.title` and `_entry.id` — header metadata.
//! - `_entity_poly_seq` — the full polymer sequence (the SEQRES
//!   equivalent), when present.
//! - `_pdbx_struct_assembly_gen` / `_pdbx_struct_oper_list` — the
//!   biological-assembly operators.
//!
//! The reader implements a real CIF tokeniser (whitespace-separated
//! values, `'…'` and `"…"` quoting, and `;…;` multi-line text
//! blocks). It does **not** evaluate `save_` frames, `data_` block
//! cross-references, or dictionary types. The writer emits a minimal
//! single-`data_` block with `_entry`, `_struct` and an `_atom_site`
//! loop — enough to round-trip coordinates.

use crate::error::{BiostructError, Result};
use crate::structure::{Atom, Chain, Model, Residue, Structure, SymmetryOperator};
use nalgebra::Point3;
use std::collections::HashMap;

/// Parse an mmCIF-format string into a [`Structure`].
pub fn read_mmcif(text: &str, id: &str) -> Result<Structure> {
    let tokens = tokenize(text);
    let mut structure = Structure::new(id);
    structure.models.clear();

    // Single-pass cursor over the token stream.
    let mut i = 0;
    let mut title: Option<String> = None;
    let mut entry_id: Option<String> = None;
    let mut atom_loop: Option<LoopTable> = None;
    let mut oper_loop: Option<LoopTable> = None;
    let mut entity_seq: Vec<(String, String)> = Vec::new(); // (asym_id, mon_id)

    while i < tokens.len() {
        match &tokens[i] {
            Token::Tag(tag) => {
                // A bare `key value` pair (not in a loop).
                if i + 1 < tokens.len() {
                    if let Token::Value(v) = &tokens[i + 1] {
                        match tag.as_str() {
                            "_struct.title" => title = Some(v.clone()),
                            "_entry.id" => entry_id = Some(v.clone()),
                            _ => {}
                        }
                        i += 2;
                        continue;
                    }
                }
                i += 1;
            }
            Token::Loop => {
                let (table, next) = parse_loop(&tokens, i + 1)?;
                i = next;
                if table.has_category("_atom_site") {
                    atom_loop = Some(table);
                } else if table.has_category("_pdbx_struct_oper_list") {
                    oper_loop = Some(table);
                } else if table.has_category("_entity_poly_seq") {
                    for row in &table.rows {
                        let asym = table
                            .get(row, "_entity_poly_seq.entity_id")
                            .unwrap_or("?")
                            .to_string();
                        let mon = table
                            .get(row, "_entity_poly_seq.mon_id")
                            .unwrap_or("UNK")
                            .to_string();
                        entity_seq.push((asym, mon));
                    }
                }
            }
            _ => i += 1,
        }
    }

    if let Some(id_v) = entry_id {
        if structure.id == id {
            structure.id = id_v;
        }
    }
    if let Some(t) = title {
        structure.title = t;
    }

    let atom_loop =
        atom_loop.ok_or_else(|| BiostructError::parse("mmcif", 0, "no _atom_site loop found"))?;
    build_models_from_atom_site(&atom_loop, &mut structure)?;

    if let Some(loops) = oper_loop {
        for row in &loops.rows {
            if let Some(op) = parse_oper_row(&loops, row) {
                structure.assembly_operators.push(op);
            }
        }
    }

    // Attach entity_poly_seq sequences. mmCIF keys by entity, but a
    // v1 maps entity -> one-letter run and assigns to chains whose
    // label matches; absent a robust entity->asym map we attach by
    // chain order when counts line up.
    if !entity_seq.is_empty() {
        let mut per_entity: HashMap<String, String> = HashMap::new();
        for (ent, mon) in &entity_seq {
            per_entity
                .entry(ent.clone())
                .or_default()
                .push(crate::structure::residue_one_letter(mon));
        }
        // If exactly one entity, give its sequence to every polymer
        // chain; otherwise leave seqres unset (observed_sequence
        // still works).
        if per_entity.len() == 1 {
            let seq = per_entity.values().next().cloned();
            for m in &mut structure.models {
                for c in &mut m.chains {
                    if c.polymer_residues().len() == seq.as_ref().map_or(0, |s| s.len()) {
                        c.seqres = seq.clone();
                    }
                }
            }
        }
    }

    if structure.models.is_empty() {
        return Err(BiostructError::parse(
            "mmcif",
            0,
            "_atom_site loop produced no models",
        ));
    }
    Ok(structure)
}

/// Serialise a [`Structure`] into a minimal mmCIF `data_` block.
pub fn write_mmcif(structure: &Structure) -> String {
    let mut out = String::new();
    out.push_str(&format!("data_{}\n#\n", sanitize_block(&structure.id)));
    out.push_str(&format!("_entry.id   {}\n#\n", quote_value(&structure.id)));
    if !structure.title.is_empty() {
        out.push_str(&format!(
            "_struct.title   {}\n#\n",
            quote_value(&structure.title)
        ));
    }
    out.push_str("loop_\n");
    for tag in ATOM_SITE_TAGS {
        out.push_str(&format!("{tag}\n"));
    }
    let mut serial = 1;
    for model in &structure.models {
        for chain in &model.chains {
            for residue in &chain.residues {
                for atom in &residue.atoms {
                    let group = if residue.hetatm { "HETATM" } else { "ATOM" };
                    out.push_str(&format!(
                        "{group} {serial} {elem} {atom_name} {alt} {resn} {chain} \
                         {seq} {ins} {x:.3} {y:.3} {z:.3} {occ:.2} {b:.2} {model}\n",
                        group = group,
                        serial = serial,
                        elem = quote_value(&atom.element),
                        atom_name = quote_value(&atom.name),
                        alt = if atom.alt_loc == ' ' {
                            ".".to_string()
                        } else {
                            atom.alt_loc.to_string()
                        },
                        resn = quote_value(&residue.name),
                        chain = quote_value(&chain.id),
                        seq = residue.seq_num,
                        ins = if residue.ins_code == ' ' {
                            "?".to_string()
                        } else {
                            residue.ins_code.to_string()
                        },
                        x = atom.coord.x,
                        y = atom.coord.y,
                        z = atom.coord.z,
                        occ = atom.occupancy,
                        b = atom.b_factor,
                        model = model.serial,
                    ));
                    serial += 1;
                }
            }
        }
    }
    out.push_str("#\n");
    out
}

/// The `_atom_site` column tags the writer emits, in order.
const ATOM_SITE_TAGS: &[&str] = &[
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_PDB_model_num",
];

// --- CIF tokeniser ---------------------------------------------------

/// A single CIF lexical token.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// A `loop_` keyword.
    Loop,
    /// A tag beginning with `_`.
    Tag(String),
    /// A data value (already unquoted).
    Value(String),
    /// A `data_…` block header (ignored after the first).
    DataBlock,
}

/// Tokenise CIF text: handles `#` comments, `'…'` / `"…"` quoting and
/// `;…;` multi-line text fields.
fn tokenize(text: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut lines = text.lines().peekable();
    while let Some(line) = lines.next() {
        // A `;` in column 1 opens a multi-line text block.
        if let Some(rest) = line.strip_prefix(';') {
            let mut buf = rest.to_string();
            for inner in lines.by_ref() {
                if inner.starts_with(';') {
                    break;
                }
                buf.push('\n');
                buf.push_str(inner);
            }
            tokens.push(Token::Value(buf.trim().to_string()));
            continue;
        }
        tokenize_line(line, &mut tokens);
    }
    tokens
}

/// Tokenise a single (non-`;`-block) line.
fn tokenize_line(line: &str, tokens: &mut Vec<Token>) {
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        if c == '#' {
            break; // comment to end of line
        }
        if c == '\'' || c == '"' {
            let quote = c;
            i += 1;
            let mut buf = String::new();
            while i < chars.len() {
                // A closing quote must be followed by whitespace or
                // end-of-line per the CIF spec.
                if chars[i] == quote && (i + 1 >= chars.len() || chars[i + 1].is_whitespace()) {
                    i += 1;
                    break;
                }
                buf.push(chars[i]);
                i += 1;
            }
            tokens.push(Token::Value(buf));
            continue;
        }
        // Bare token.
        let start = i;
        while i < chars.len() && !chars[i].is_whitespace() {
            i += 1;
        }
        let word: String = chars[start..i].iter().collect();
        let lower = word.to_ascii_lowercase();
        if lower == "loop_" {
            tokens.push(Token::Loop);
        } else if lower.starts_with("data_") {
            tokens.push(Token::DataBlock);
        } else if word.starts_with('_') {
            tokens.push(Token::Tag(word));
        } else {
            tokens.push(Token::Value(word));
        }
    }
}

/// A parsed `loop_` table: a list of column tags and the row values.
struct LoopTable {
    tags: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl LoopTable {
    /// Whether any column tag belongs to the named category
    /// (`"_atom_site"`, …).
    fn has_category(&self, category: &str) -> bool {
        let prefix = format!("{category}.");
        self.tags.iter().any(|t| t.starts_with(&prefix))
    }

    /// Column index of an exact tag.
    fn col(&self, tag: &str) -> Option<usize> {
        self.tags.iter().position(|t| t == tag)
    }

    /// Fetch a cell of `row` by tag. Returns `None` for unknown tags
    /// and treats CIF nulls `.` / `?` as `None`.
    fn get<'a>(&self, row: &'a [String], tag: &str) -> Option<&'a str> {
        let idx = self.col(tag)?;
        let v = row.get(idx)?.as_str();
        if v == "." || v == "?" {
            None
        } else {
            Some(v)
        }
    }
}

/// Parse a `loop_` table starting at token index `start` (just after
/// the `loop_` keyword). Returns the table and the index of the first
/// token after it.
fn parse_loop(tokens: &[Token], start: usize) -> Result<(LoopTable, usize)> {
    let mut tags = Vec::new();
    let mut i = start;
    while i < tokens.len() {
        if let Token::Tag(t) = &tokens[i] {
            tags.push(t.clone());
            i += 1;
        } else {
            break;
        }
    }
    if tags.is_empty() {
        return Err(BiostructError::parse("mmcif", 0, "loop_ with no tags"));
    }
    let ncol = tags.len();
    let mut rows = Vec::new();
    let mut current: Vec<String> = Vec::with_capacity(ncol);
    while i < tokens.len() {
        match &tokens[i] {
            Token::Value(v) => {
                current.push(v.clone());
                if current.len() == ncol {
                    rows.push(std::mem::replace(&mut current, Vec::with_capacity(ncol)));
                }
                i += 1;
            }
            // Any non-value token ends the loop.
            _ => break,
        }
    }
    Ok((LoopTable { tags, rows }, i))
}

/// Build coordinate models from a parsed `_atom_site` loop.
fn build_models_from_atom_site(table: &LoopTable, structure: &mut Structure) -> Result<()> {
    // Resolve column indices once. mmCIF files vary in which optional
    // columns they carry; only the coordinates are mandatory.
    let need = |tag: &str| -> Result<usize> {
        table
            .col(tag)
            .ok_or_else(|| BiostructError::parse("mmcif", 0, format!("missing {tag} column")))
    };
    let c_x = need("_atom_site.Cartn_x")?;
    let c_y = need("_atom_site.Cartn_y")?;
    let c_z = need("_atom_site.Cartn_z")?;

    // mmCIF auth_* columns mirror PDB numbering; prefer them when
    // present so the residue ids match a PDB of the same entry.
    let c_atom = table
        .col("_atom_site.label_atom_id")
        .or_else(|| table.col("_atom_site.auth_atom_id"));
    let c_comp = table
        .col("_atom_site.auth_comp_id")
        .or_else(|| table.col("_atom_site.label_comp_id"));
    let c_asym = table
        .col("_atom_site.auth_asym_id")
        .or_else(|| table.col("_atom_site.label_asym_id"));
    let c_seq = table
        .col("_atom_site.auth_seq_id")
        .or_else(|| table.col("_atom_site.label_seq_id"));
    let c_group = table.col("_atom_site.group_PDB");
    let c_alt = table.col("_atom_site.label_alt_id");
    let c_elem = table.col("_atom_site.type_symbol");
    let c_occ = table.col("_atom_site.occupancy");
    let c_b = table.col("_atom_site.B_iso_or_equiv");
    let c_ins = table.col("_atom_site.pdbx_PDB_ins_code");
    let c_model = table.col("_atom_site.pdbx_PDB_model_num");
    let c_serial = table.col("_atom_site.id");

    // model_serial -> Model
    let mut models: Vec<Model> = Vec::new();
    let cell = |row: &[String], idx: Option<usize>| -> String {
        idx.and_then(|i| row.get(i)).cloned().unwrap_or_default()
    };

    for (rownum, row) in table.rows.iter().enumerate() {
        let x: f64 = row.get(c_x).and_then(|v| v.parse().ok()).ok_or_else(|| {
            BiostructError::parse("mmcif", 0, format!("bad x in row {}", rownum + 1))
        })?;
        let y: f64 = row.get(c_y).and_then(|v| v.parse().ok()).ok_or_else(|| {
            BiostructError::parse("mmcif", 0, format!("bad y in row {}", rownum + 1))
        })?;
        let z: f64 = row.get(c_z).and_then(|v| v.parse().ok()).ok_or_else(|| {
            BiostructError::parse("mmcif", 0, format!("bad z in row {}", rownum + 1))
        })?;

        let model_serial: i32 = cell(row, c_model).parse().unwrap_or(1);
        let group = cell(row, c_group);
        let hetatm = group.eq_ignore_ascii_case("HETATM");
        let atom_name = {
            let n = cell(row, c_atom);
            // mmCIF atom names may carry double-quote padding already
            // stripped by the tokeniser.
            n.trim().to_string()
        };
        let alt_raw = cell(row, c_alt);
        let alt_loc = if alt_raw == "." || alt_raw == "?" || alt_raw.is_empty() {
            ' '
        } else {
            alt_raw.chars().next().unwrap_or(' ')
        };
        let res_name = cell(row, c_comp).trim().to_ascii_uppercase();
        let chain_id = {
            let c = cell(row, c_asym);
            if c.is_empty() {
                "A".to_string()
            } else {
                c
            }
        };
        let seq_num: i32 = cell(row, c_seq).parse().unwrap_or(0);
        let ins_raw = cell(row, c_ins);
        let ins_code = if ins_raw == "." || ins_raw == "?" || ins_raw.is_empty() {
            ' '
        } else {
            ins_raw.chars().next().unwrap_or(' ')
        };
        let element = {
            let e = cell(row, c_elem).trim().to_ascii_uppercase();
            if e.is_empty() {
                crate::io::pdb::guess_element(&atom_name)
            } else {
                e
            }
        };
        let occupancy: f64 = cell(row, c_occ).parse().unwrap_or(1.0);
        let b_factor: f64 = cell(row, c_b).parse().unwrap_or(0.0);
        let serial: i32 = cell(row, c_serial).parse().unwrap_or(rownum as i32 + 1);

        let model = match models.iter_mut().find(|m| m.serial == model_serial) {
            Some(m) => m,
            None => {
                models.push(Model::new(model_serial));
                models.last_mut().unwrap()
            }
        };
        let chain = match model.chains.iter_mut().find(|c| c.id == chain_id) {
            Some(c) => c,
            None => {
                model.chains.push(Chain::new(&chain_id));
                model.chains.last_mut().unwrap()
            }
        };
        let atom = Atom {
            serial,
            name: atom_name,
            alt_loc,
            element,
            coord: Point3::new(x, y, z),
            occupancy,
            b_factor,
            charge: 0,
        };
        let matches_tail = chain
            .residues
            .last()
            .map(|r| r.seq_num == seq_num && r.ins_code == ins_code && r.name == res_name)
            .unwrap_or(false);
        if matches_tail {
            chain.residues.last_mut().unwrap().atoms.push(atom);
        } else {
            let mut r = Residue::new(&res_name, seq_num);
            r.ins_code = ins_code;
            r.hetatm = hetatm;
            r.atoms.push(atom);
            chain.residues.push(r);
        }
    }

    models.sort_by_key(|m| m.serial);
    structure.models = models;
    Ok(())
}

/// Parse one `_pdbx_struct_oper_list` row into a [`SymmetryOperator`].
fn parse_oper_row(table: &LoopTable, row: &[String]) -> Option<SymmetryOperator> {
    let serial: i32 = table
        .get(row, "_pdbx_struct_oper_list.id")
        .and_then(|v| v.parse().ok())
        .unwrap_or(1);
    let mut rotation = [[0.0; 3]; 3];
    let mut translation = [0.0; 3];
    for (i, &ti) in [1, 2, 3].iter().enumerate() {
        for (j, &tj) in [1, 2, 3].iter().enumerate() {
            let tag = format!("_pdbx_struct_oper_list.matrix[{ti}][{tj}]");
            rotation[i][j] = table.get(row, &tag).and_then(|v| v.parse().ok())?;
        }
        let vt = format!("_pdbx_struct_oper_list.vector[{ti}]");
        translation[i] = table
            .get(row, &vt)
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0);
    }
    Some(SymmetryOperator {
        serial,
        rotation,
        translation,
    })
}

// --- writer helpers --------------------------------------------------

/// Quote a CIF value if it contains whitespace or is empty.
fn quote_value(v: &str) -> String {
    if v.is_empty() {
        return "?".to_string();
    }
    if v.chars().any(|c| c.is_whitespace()) || v.starts_with('_') || v.starts_with('#') {
        if v.contains('\'') {
            format!("\"{v}\"")
        } else {
            format!("'{v}'")
        }
    } else {
        v.to_string()
    }
}

/// Make a CIF `data_` block name (no whitespace).
fn sanitize_block(id: &str) -> String {
    let s: String = id
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect();
    if s.is_empty() {
        "structure".to_string()
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINI_CIF: &str = "\
data_1ABC
#
_entry.id   1ABC
#
_struct.title   'A small test structure'
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N   . ALA A 1 ? 11.104 6.134 -6.504 1.00 20.00 1
ATOM   2  C  CA  . ALA A 1 ? 11.639 6.071 -5.147 1.00 21.00 1
ATOM   3  C  C   . ALA A 1 ? 13.085 6.564 -5.180 1.00 19.50 1
ATOM   4  N  N   . GLY A 2 ? 13.965 5.628 -4.926 1.00 17.00 1
ATOM   5  C  CA  . GLY A 2 ? 15.398 5.927 -4.911 1.00 16.00 1
HETATM 6  ZN ZN  . ZN  A 3 ? 20.000 20.000 20.000 1.00 30.00 1
#
";

    #[test]
    fn tokenizer_handles_quotes() {
        let toks = tokenize("_struct.title   'A small test structure'\n");
        assert_eq!(toks.len(), 2);
        assert_eq!(toks[1], Token::Value("A small test structure".to_string()));
    }

    #[test]
    fn tokenizer_handles_multiline_text() {
        let toks = tokenize("_x\n;line one\nline two\n;\n");
        assert!(matches!(&toks[1], Token::Value(v) if v.contains("line one")));
    }

    #[test]
    fn reads_atom_site_loop() {
        let s = read_mmcif(MINI_CIF, "fb").unwrap();
        assert_eq!(s.id, "1ABC");
        assert_eq!(s.title, "A small test structure");
        let m = s.first_model();
        let chain = m.chain("A").unwrap();
        assert_eq!(chain.residues.len(), 3); // ALA GLY ZN
        assert_eq!(chain.residues[0].atoms.len(), 3);
        assert_eq!(chain.residues[0].ca().unwrap().b_factor, 21.0);
    }

    #[test]
    fn reads_hetatm() {
        let s = read_mmcif(MINI_CIF, "x").unwrap();
        let zn = s.first_model().chain("A").unwrap().residue(3, ' ').unwrap();
        assert!(zn.hetatm);
        assert_eq!(zn.atoms[0].element, "ZN");
    }

    #[test]
    fn round_trips_coordinates() {
        let s = read_mmcif(MINI_CIF, "x").unwrap();
        let text = write_mmcif(&s);
        let s2 = read_mmcif(&text, "x").unwrap();
        assert_eq!(s2.atom_count(), s.atom_count());
        let a1 = s.first_model().atoms().next().unwrap().coord;
        let a2 = s2.first_model().atoms().next().unwrap().coord;
        assert!((a1 - a2).norm() < 1e-3);
    }

    #[test]
    fn reads_oper_list() {
        let cif = "\
data_x
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 1 0 0 0 1 0 0 0
2 -1 0 0 0 -1 0 0 0 1 5 0 0
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA ALA A 1 0 0 0
#
";
        let s = read_mmcif(cif, "x").unwrap();
        assert_eq!(s.assembly_operators.len(), 2);
        assert_eq!(s.assembly_operators[1].translation[0], 5.0);
    }

    #[test]
    fn missing_atom_site_errors() {
        assert!(read_mmcif("data_x\n_entry.id  x\n", "x").is_err());
    }
}

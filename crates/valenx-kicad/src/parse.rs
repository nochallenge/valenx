//! Minimal `.kicad_pcb` S-expression reader. v1 covers outline +
//! drills + components; nets/tracks/zones defer to Phase 42.5.

use std::path::Path;

use nalgebra::Vector3;

use crate::board::{Component, KicadBoard};
use crate::error::KicadError;

/// Read a `.kicad_pcb` file. v1 minimal — outline + drills +
/// components, no nets / tracks / zones.
///
/// Round-23 sweep: bounded at [`valenx_core::io_caps::MAX_KICAD_FILE_BYTES`]
/// (1 GiB) — production PCBs with thousands of footprints cross
/// 100 MiB; 1 GiB matches the DXF cap.
pub fn import_kicad_pcb(path: impl AsRef<Path>) -> Result<KicadBoard, KicadError> {
    let text = valenx_core::io_caps::read_capped_to_string(
        path.as_ref(),
        valenx_core::io_caps::MAX_KICAD_FILE_BYTES as usize,
    )?;
    from_str(&text)
}

/// Parse a `.kicad_pcb` string. Hand-rolled S-expression tokenizer
/// because we only need a tiny subset of the surface (outline +
/// drills + footprints).
pub fn from_str(text: &str) -> Result<KicadBoard, KicadError> {
    let tokens = tokenize(text);
    let (root, _) = parse_sexpr(&tokens, 0, 0)?;
    let Sexpr::List(top) = &root else {
        return Err(KicadError::Parse("root must be a list".into()));
    };
    if top.first().map(name_of) != Some("kicad_pcb") {
        return Err(KicadError::Parse(format!(
            "expected kicad_pcb root, got {:?}",
            top.first().map(name_of)
        )));
    }

    let mut board = KicadBoard::new_default();

    for child in top.iter().skip(1) {
        if let Sexpr::List(items) = child {
            match items.first().map(name_of) {
                Some("general") => parse_general(items, &mut board),
                Some("gr_line") => parse_gr_line(items, &mut board),
                Some("via") => parse_via_as_drill(items, &mut board),
                Some("footprint") | Some("module") => {
                    parse_footprint(items, &mut board);
                }
                _ => {}
            }
        }
    }

    Ok(board)
}

fn parse_general(items: &[Sexpr], board: &mut KicadBoard) {
    for it in items {
        if let Sexpr::List(inner) = it {
            if name_of_first(inner) == "thickness" {
                if let Some(Sexpr::Atom(s)) = inner.get(1) {
                    if let Ok(v) = s.parse::<f64>() {
                        board.thickness_mm = v;
                    }
                }
            }
        }
    }
}

fn parse_gr_line(items: &[Sexpr], board: &mut KicadBoard) {
    // gr_line (start x y) (end x y) (layer Edge.Cuts) ...
    let mut layer = String::new();
    let mut start: Option<[f64; 2]> = None;
    let mut end: Option<[f64; 2]> = None;
    for it in items.iter().skip(1) {
        if let Sexpr::List(inner) = it {
            match name_of_first(inner) {
                "start" => start = parse_xy(inner),
                "end" => end = parse_xy(inner),
                "layer" => {
                    if let Some(Sexpr::Atom(s)) = inner.get(1) {
                        layer = s.clone();
                    }
                }
                _ => {}
            }
        }
    }
    if layer == "Edge.Cuts" {
        if let (Some(s), Some(e)) = (start, end) {
            // Build the outline by chaining endpoints. We push start
            // when the outline is empty, then keep extending with the
            // newest end point.
            if board.outline.is_empty() {
                board.outline.push(s);
            }
            board.outline.push(e);
        }
    }
}

fn parse_via_as_drill(items: &[Sexpr], board: &mut KicadBoard) {
    let mut pos = Vector3::zeros();
    let mut drill = 0.6;
    for it in items.iter().skip(1) {
        if let Sexpr::List(inner) = it {
            match name_of_first(inner) {
                "at" => {
                    if let Some(xy) = parse_xy(inner) {
                        pos = Vector3::new(xy[0], xy[1], 0.0);
                    }
                }
                "drill" => {
                    if let Some(Sexpr::Atom(s)) = inner.get(1) {
                        if let Ok(v) = s.parse::<f64>() {
                            drill = v;
                        }
                    }
                }
                _ => {}
            }
        }
    }
    board.drill_holes.push((pos, drill));
}

fn parse_footprint(items: &[Sexpr], board: &mut KicadBoard) {
    let footprint_name = items
        .get(1)
        .and_then(|s| match s {
            Sexpr::Atom(a) => Some(a.clone()),
            _ => None,
        })
        .unwrap_or_default();
    let mut at = Vector3::zeros();
    let mut rot = 0.0;
    let mut reference = String::new();
    let mut model_3d_path: Option<String> = None;
    for it in items.iter().skip(2) {
        if let Sexpr::List(inner) = it {
            match name_of_first(inner) {
                "at" => {
                    if let Some(xy) = parse_xy(inner) {
                        at = Vector3::new(xy[0], xy[1], 0.0);
                    }
                    if let Some(Sexpr::Atom(s)) = inner.get(3) {
                        if let Ok(v) = s.parse::<f64>() {
                            rot = v;
                        }
                    }
                }
                "fp_text" => {
                    // (fp_text reference R1 ...)
                    if let Some(Sexpr::Atom(kind)) = inner.get(1) {
                        if kind == "reference" {
                            if let Some(Sexpr::Atom(s)) = inner.get(2) {
                                reference = s.clone();
                            }
                        }
                    }
                }
                "model" => {
                    if let Some(Sexpr::Atom(s)) = inner.get(1) {
                        model_3d_path = Some(s.clone());
                    }
                }
                _ => {}
            }
        }
    }
    board.components.push(Component {
        ref_designator: reference,
        footprint_name,
        position: at,
        rotation_deg: rot,
        model_3d_path,
    });
}

fn parse_xy(list: &[Sexpr]) -> Option<[f64; 2]> {
    let x = match list.get(1)? {
        Sexpr::Atom(s) => s.parse::<f64>().ok()?,
        _ => return None,
    };
    let y = match list.get(2)? {
        Sexpr::Atom(s) => s.parse::<f64>().ok()?,
        _ => return None,
    };
    Some([x, y])
}

// ---- S-expression mini-parser ----

#[derive(Clone, Debug)]
enum Sexpr {
    Atom(String),
    List(Vec<Sexpr>),
}

fn name_of(s: &Sexpr) -> &str {
    match s {
        Sexpr::Atom(a) => a.as_str(),
        Sexpr::List(_) => "",
    }
}

/// The head-symbol name of a list's first element, or `""` if the list is
/// empty. Used so a malformed empty `()` can't panic `inner[0]`.
fn name_of_first(items: &[Sexpr]) -> &str {
    items.first().map(name_of).unwrap_or("")
}

fn tokenize(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_str = false;
    for c in text.chars() {
        if in_str {
            if c == '"' {
                out.push(cur.clone());
                cur.clear();
                in_str = false;
            } else {
                cur.push(c);
            }
            continue;
        }
        match c {
            '(' | ')' => {
                if !cur.is_empty() {
                    out.push(cur.clone());
                    cur.clear();
                }
                out.push(c.to_string());
            }
            '"' => in_str = true,
            c if c.is_whitespace() => {
                if !cur.is_empty() {
                    out.push(cur.clone());
                    cur.clear();
                }
            }
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

fn parse_sexpr(tokens: &[String], idx: usize, depth: usize) -> Result<(Sexpr, usize), KicadError> {
    // Bound recursion: a deeply nested `(((…` would otherwise overflow the
    // stack (the 1 GiB file cap allows millions of nesting levels).
    const MAX_DEPTH: usize = 512;
    if depth > MAX_DEPTH {
        return Err(KicadError::Parse("s-expression nested too deeply".into()));
    }
    if idx >= tokens.len() {
        return Err(KicadError::Parse("unexpected EOF".into()));
    }
    let t = &tokens[idx];
    if t == "(" {
        let mut items = Vec::new();
        let mut i = idx + 1;
        while i < tokens.len() && tokens[i] != ")" {
            let (inner, next) = parse_sexpr(tokens, i, depth + 1)?;
            items.push(inner);
            i = next;
        }
        if i >= tokens.len() {
            return Err(KicadError::Parse("unmatched (".into()));
        }
        Ok((Sexpr::List(items), i + 1))
    } else {
        Ok((Sexpr::Atom(t.clone()), idx + 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_minimal_pcb() {
        let text = "(kicad_pcb (general (thickness 1.6)) \
            (gr_line (start 0 0) (end 100 0) (layer Edge.Cuts)) \
            (gr_line (start 100 0) (end 100 80) (layer Edge.Cuts)) \
            (gr_line (start 100 80) (end 0 80) (layer Edge.Cuts)) \
            (gr_line (start 0 80) (end 0 0) (layer Edge.Cuts)) \
            (via (at 5 5) (drill 0.6)) \
            (footprint \"Resistor_SMD:R_0805\" (at 50 40 90) \
              (fp_text reference R1) \
              (model \"R_0805.step\"))) ";
        let board = from_str(text).expect("parse");
        assert!((board.thickness_mm - 1.6).abs() < 1e-9);
        assert!(board.outline.len() >= 4);
        assert_eq!(board.drill_holes.len(), 1);
        assert_eq!(board.components.len(), 1);
        assert_eq!(board.components[0].ref_designator, "R1");
        assert!(board.components[0].rotation_deg == 90.0);
    }

    #[test]
    fn empty_section_list_does_not_panic() {
        // A malformed empty `()` inside a section must be skipped, not panic
        // `inner[0]`. (Previously `name_of(&inner[0])` indexed an empty Vec.)
        assert!(from_str("(kicad_pcb (general ()))").is_ok());
        assert!(from_str("(kicad_pcb (gr_line ()))").is_ok());
    }

    #[test]
    fn deeply_nested_input_errors_not_stack_overflow() {
        // A pathologically nested s-expression must hit the depth cap and
        // return an error rather than overflow the stack.
        assert!(from_str(&"(".repeat(5000)).is_err());
    }

    #[test]
    fn rejects_non_kicad_root() {
        let text = "(other (x 1))";
        assert!(matches!(from_str(text), Err(KicadError::Parse(_))));
    }
}

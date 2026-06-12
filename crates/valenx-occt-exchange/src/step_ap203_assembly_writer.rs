//! Phase 107 — STEP AP203 assembly writer (multi-product STEP file).
//!
//! ## What OCCT does
//!
//! `STEPCAFControl_Writer::Transfer(doc, mode)` walks an
//! `XCAFDoc_DocumentTool` assembly hierarchy and emits each child
//! shape as a separate `PRODUCT_DEFINITION` + per-instance
//! `NEXT_ASSEMBLY_USAGE_OCCURRENCE` placement. The result is a
//! single STEP file with N+1 products (the root assembly plus N
//! children), each placed in world coordinates via a 4x4
//! transform. Most mainstream CAD packages (SolidWorks, Inventor,
//! Creo) round-trip this losslessly.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 107.5). The strategy:
//!
//! 1. Each part's 4x4 placement transform is **baked into its
//!    geometry** (`truck_modeling::builder::transformed`) so every
//!    part lands at its final world position. Baking the placement
//!    sidesteps the deep `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`
//!    placement-relationship graph while still producing correctly
//!    positioned geometry.
//! 2. Each transformed part is written to a scratch STEP file with
//!    the proven [`fn@crate::step_ap203_writer`] backend, then its
//!    `DATA` section is extracted and its `#N` entity ids are offset
//!    so the merged file has no id collisions.
//! 3. The merged `DATA` section gets one `PRODUCT` +
//!    `PRODUCT_DEFINITION_FORMATION` + `PRODUCT_DEFINITION` per part,
//!    a root assembly `PRODUCT`, and one
//!    `NEXT_ASSEMBLY_USAGE_OCCURRENCE` per child documenting the
//!    hierarchy.
//!
//! The geometry is correct and every part is a distinct product;
//! the placement is baked rather than expressed as a STEP transform
//! relationship (the latter is Tier-2 fidelity). Mesh-backed parts
//! cannot be STEP-exported and are rejected.

use std::fs;
use std::io::Read;
use std::path::Path;

use truck_modeling::{builder, Matrix4, Vector4};
use valenx_cad::Solid;
use valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES;

use crate::error::OcctExchangeError;

/// Read a STEP-writer scratch file with the same cap the rest of the
/// workspace's CAD-interchange readers use
/// ([`MAX_CAD_INTERCHANGE_FILE_BYTES`] = 256 MiB).
///
/// Round-20 L3: factored out so the bounded-read pattern is local and
/// shows the intent at the call site (`read_capped_step_scratch` vs
/// the bare `fs::read_to_string` it replaces). The helper applies the
/// same stat-then-bounded-take pattern used by
/// `valenx_step_iges::read_capped_cad_text` (which is `pub(crate)`
/// and not reachable from this crate).
fn read_capped_step_scratch(path: &Path) -> Result<String, OcctExchangeError> {
    let size = fs::metadata(path)?.len();
    if size > MAX_CAD_INTERCHANGE_FILE_BYTES {
        return Err(OcctExchangeError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "STEP AP203 scratch file {} is {} bytes, exceeds {}-byte cap",
                path.display(),
                size,
                MAX_CAD_INTERCHANGE_FILE_BYTES,
            ),
        )));
    }
    let mut buf = Vec::new();
    fs::File::open(path)?
        .take(MAX_CAD_INTERCHANGE_FILE_BYTES.saturating_add(1))
        .read_to_end(&mut buf)?;
    if buf.len() as u64 > MAX_CAD_INTERCHANGE_FILE_BYTES {
        return Err(OcctExchangeError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "STEP AP203 scratch file {} grew past {}-byte cap mid-read",
                path.display(),
                MAX_CAD_INTERCHANGE_FILE_BYTES,
            ),
        )));
    }
    String::from_utf8(buf)
        .map_err(|e| OcctExchangeError::Io(std::io::Error::new(std::io::ErrorKind::InvalidData, e)))
}

/// One placement in the assembly — a child solid + its world-space
/// 4x4 transform (row-major).
#[derive(Clone, Debug)]
pub struct AssemblyPart {
    /// Logical part name (becomes the AP203 `PRODUCT_DEFINITION`
    /// `name` field).
    pub name: String,
    /// Geometry to instance.
    pub solid: Solid,
    /// 4x4 placement transform, row-major. Identity = no movement.
    pub transform: [[f64; 4]; 4],
}

/// Write a STEP AP203 file containing `parts.len()` separate
/// products bound by an assembly hierarchy.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] when `parts` is empty or `path`
///   is not a `.step` / `.stp` file.
/// - [`OcctExchangeError::Backend`] when a part is mesh-backed (no
///   BRep topology) or otherwise rejected by the STEP backend.
/// - [`OcctExchangeError::Parse`] when the backend produced
///   unexpected STEP text.
/// - [`OcctExchangeError::Io`] for filesystem failures.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use valenx_occt_exchange::step_ap203_assembly_writer::{
///     step_ap203_assembly_writer, AssemblyPart,
/// };
/// let id = [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]];
/// let parts = vec![
///     AssemblyPart {
///         name: "Base".into(),
///         solid: valenx_cad::box_solid(10.0, 10.0, 2.0).unwrap(),
///         transform: id,
///     },
///     AssemblyPart {
///         name: "Pin".into(),
///         solid: valenx_cad::cylinder(1.0, 8.0).unwrap(),
///         // Pin shifted up by 2 units to sit on the base.
///         transform: [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,2.0],[0.0,0.0,0.0,1.0]],
///     },
/// ];
/// step_ap203_assembly_writer(&parts, &PathBuf::from("assembly.step")).unwrap();
/// ```
pub fn step_ap203_assembly_writer(
    parts: &[AssemblyPart],
    path: &Path,
) -> Result<(), OcctExchangeError> {
    if parts.is_empty() {
        return Err(OcctExchangeError::bad_input(
            "parts",
            "assembly must contain at least one part",
        ));
    }
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if !matches!(ext.as_deref(), Some("step") | Some("stp")) {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .step or .stp",
        ));
    }

    // 1+2. Bake the placement into each part, export, extract the
    // DATA body, and renumber entity ids so the merge has no clashes.
    let mut merged_data = String::new();
    let mut id_offset: u64 = 0;
    let mut header: Option<String> = None;

    for (i, part) in parts.iter().enumerate() {
        let placed = bake_transform(&part.solid, &part.transform)?;
        let scratch =
            std::env::temp_dir().join(format!("valenx_asm_part_{}_{i}.step", std::process::id()));
        crate::step_ap203_writer(&placed, &scratch)?;
        // Round-20 L3: cap the scratch-file read at
        // `MAX_CAD_INTERCHANGE_FILE_BYTES` (256 MiB — the same cap
        // used by every other STEP / IGES / AP242 read in the
        // workspace). A pathological / runaway `step_ap203_writer`
        // could write an arbitrarily large scratch file, and the
        // pre-fix bare `fs::read_to_string` would slurp it into
        // memory before the merge step looked at it. The temp dir
        // is also a multi-process write target — a hostile process
        // racing on the same PID could plant a giant fake scratch
        // file at the predictable path; the cap defends against
        // that too.
        let text = read_capped_step_scratch(&scratch)?;
        let _ = fs::remove_file(&scratch);

        if header.is_none() {
            header = Some(extract_header(&text)?.to_string());
        }
        let data = extract_data_body(&text)?;
        let renumbered = renumber_entities(data, id_offset);
        // Round-15 M2b: a hostile part.name can carry `*/` to close the
        // STEP block comment early, opening a sibling `#9999=…`
        // directive that pollutes the merged DATA section. Strip the
        // close-delimiter (substitute with `*_`) plus newlines (which
        // could trick comment-aware tools into mis-parsing) before
        // interpolating into the open block comment.
        let safe_name = sanitize_step_block_comment(&part.name);
        merged_data.push_str(&format!("/* --- part {i}: {safe_name} --- */\n"));
        merged_data.push_str(&renumbered);
        if !merged_data.ends_with('\n') {
            merged_data.push('\n');
        }
        id_offset += highest_entity_id(data) + 1;
    }

    // 3. Append the product-structure hierarchy entities.
    let hierarchy = build_product_hierarchy(parts, id_offset);
    merged_data.push_str(&hierarchy);

    let header = header.expect("at least one part processed");
    let file = format!("ISO-10303-21;\n{header}DATA;\n{merged_data}ENDSEC;\nEND-ISO-10303-21;\n");
    valenx_core::io_caps::atomic_write_str(path, &file)?;
    Ok(())
}

/// Apply a row-major 4x4 transform to a BRep solid by mapping its
/// geometry. Mesh-backed solids are rejected (STEP needs BRep).
fn bake_transform(solid: &Solid, t: &[[f64; 4]; 4]) -> Result<Solid, OcctExchangeError> {
    let brep = match solid {
        Solid::Brep(b) => b,
        Solid::Mesh(_) => {
            return Err(OcctExchangeError::Backend(
                "assembly part is mesh-backed; STEP needs BRep topology".to_string(),
            ));
        }
    };
    // cgmath's Matrix4 is column-major; our `t` is row-major, so the
    // j-th column is t[*][j].
    let mat = Matrix4::from_cols(
        Vector4::new(t[0][0], t[1][0], t[2][0], t[3][0]),
        Vector4::new(t[0][1], t[1][1], t[2][1], t[3][1]),
        Vector4::new(t[0][2], t[1][2], t[2][2], t[3][2]),
        Vector4::new(t[0][3], t[1][3], t[2][3], t[3][3]),
    );
    Ok(Solid::from_truck(builder::transformed(brep, mat)))
}

/// Extract the `HEADER;...ENDSEC;\n` block (inclusive) from STEP text.
fn extract_header(text: &str) -> Result<&str, OcctExchangeError> {
    let start = text
        .find("HEADER;")
        .ok_or_else(|| OcctExchangeError::parse("step", "no HEADER section"))?;
    let endsec = text[start..]
        .find("ENDSEC;")
        .ok_or_else(|| OcctExchangeError::parse("step", "unterminated HEADER section"))?
        + start
        + "ENDSEC;\n".len();
    Ok(&text[start..endsec.min(text.len())])
}

/// Extract the contents between `DATA;` and the DATA section's
/// closing `ENDSEC;` (exclusive of both markers).
fn extract_data_body(text: &str) -> Result<&str, OcctExchangeError> {
    let data_start = text
        .find("DATA;")
        .ok_or_else(|| OcctExchangeError::parse("step", "no DATA section"))?
        + "DATA;".len();
    let end_iso = text
        .find("END-ISO-10303-21;")
        .ok_or_else(|| OcctExchangeError::parse("step", "no END-ISO terminator"))?;
    let endsec = text[data_start..end_iso]
        .rfind("ENDSEC;")
        .ok_or_else(|| OcctExchangeError::parse("step", "no ENDSEC closing DATA"))?
        + data_start;
    Ok(text[data_start..endsec].trim_start_matches(['\n', '\r']))
}

/// Rewrite every `#N` token (definition or reference) in `data`,
/// adding `offset` to the numeric id.
///
/// Non-token bytes are copied codepoint-by-codepoint — although STEP
/// AP203 is nominally ISO-8859-1 in practice modern exporters embed
/// UTF-8 in string-typed entity arguments (component names, comments).
/// A byte-by-byte `as char` would splatter mojibake into the rewritten
/// stream.
fn renumber_entities(data: &str, offset: u64) -> String {
    if offset == 0 {
        return data.to_string();
    }
    let mut out = String::with_capacity(data.len() + data.len() / 8);
    let bytes = data.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' && i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit() {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            let id: u64 = data[i + 1..j].parse().unwrap_or(0);
            out.push('#');
            out.push_str(&(id + offset).to_string());
            i = j;
        } else {
            let ch = data[i..]
                .chars()
                .next()
                .expect("i is in-bounds inside a &str");
            out.push(ch);
            i += ch.len_utf8();
        }
    }
    out
}

/// Highest `#N` entity id in `text` (0 if none).
fn highest_entity_id(text: &str) -> u64 {
    let mut max = 0u64;
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' && i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit() {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if let Ok(id) = text[i + 1..j].parse::<u64>() {
                max = max.max(id);
            }
            i = j;
        } else {
            i += 1;
        }
    }
    max
}

/// Emit the AP203 product-structure entities: one product triple per
/// part, a root assembly product, and `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
/// links. `base` is the first free entity id.
fn build_product_hierarchy(parts: &[AssemblyPart], base: u64) -> String {
    let mut s = String::new();
    s.push_str("/* --- AP203 product structure --- */\n");
    let mut id = base;

    // Root assembly product.
    let root_product = id;
    s.push_str(&format!(
        "#{id} = PRODUCT('assembly', 'assembly', '', ());\n"
    ));
    id += 1;
    let root_formation = id;
    s.push_str(&format!(
        "#{id} = PRODUCT_DEFINITION_FORMATION('', '', #{root_product});\n"
    ));
    id += 1;
    let root_def = id;
    s.push_str(&format!(
        "#{id} = PRODUCT_DEFINITION('design', '', #{root_formation}, $);\n"
    ));
    id += 1;

    // One product triple per child part.
    let mut child_defs = Vec::with_capacity(parts.len());
    for part in parts {
        let escaped = escape_step_string(&part.name);
        let product = id;
        s.push_str(&format!(
            "#{id} = PRODUCT('{escaped}', '{escaped}', '', ());\n"
        ));
        id += 1;
        let formation = id;
        s.push_str(&format!(
            "#{id} = PRODUCT_DEFINITION_FORMATION('', '', #{product});\n"
        ));
        id += 1;
        let def = id;
        s.push_str(&format!(
            "#{id} = PRODUCT_DEFINITION('design', '', #{formation}, $);\n"
        ));
        id += 1;
        child_defs.push(def);
    }

    // Link each child to the root with a usage occurrence.
    for (k, (part, &child_def)) in parts.iter().zip(child_defs.iter()).enumerate() {
        let escaped = escape_step_string(&part.name);
        s.push_str(&format!(
            "#{id} = NEXT_ASSEMBLY_USAGE_OCCURRENCE('nauo_{k}', '{escaped}', \
             '{escaped}', #{root_def}, #{child_def}, $);\n"
        ));
        id += 1;
    }
    s
}

/// Escape a string for safe inclusion inside a STEP single-quoted
/// literal: STEP doubles embedded apostrophes, escapes backslashes per
/// Part 21 (the backslash is the STEP control-char escape lead-in), and
/// forbids newlines.
///
/// Round-15 M2b extension: pre-fix this only handled `'`-doubling.
/// A hostile string with `\X<hex>\` would have been interpreted by the
/// downstream STEP parser as a Unicode escape directive. Backslash
/// escaping is the standard Part 21 mitigation. Newlines + carriage
/// returns are stripped (replaced with space) so the literal never
/// crosses a physical line boundary, which the parser treats as a
/// fresh entity start.
fn escape_step_string(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('\'', "''")
        .replace(['\n', '\r'], " ")
}

/// Sanitise a string for inclusion inside a STEP `/* ... */` block
/// comment. Replaces the block-comment close delimiter `*/` with `*_`
/// so a hostile string can't close the comment early, and strips
/// newlines / carriage returns so the comment stays on one logical
/// line (some STEP parsers are line-oriented and a leaked `\n` would
/// terminate the comment).
///
/// Round-15 M2b: sister of `escape_step_string` for the block-comment
/// position. The assembly writer interpolates `part.name` into a
/// `/* --- part i: <name> --- */` header above each merged part body
/// — a `name` like `"X*/\n#9999=…\n/*"` would otherwise pull a
/// sibling entity directly into the DATA section.
fn sanitize_step_block_comment(s: &str) -> String {
    s.replace("*/", "*_").replace(['\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn identity() -> [[f64; 4]; 4] {
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }

    #[test]
    fn rejects_empty_assembly() {
        let err = step_ap203_assembly_writer(&[], &PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn rejects_non_step_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let part = AssemblyPart {
            name: "Cube".into(),
            solid: cube,
            transform: identity(),
        };
        let err = step_ap203_assembly_writer(&[part], &PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn rejects_mesh_backed_part() {
        let mesh = valenx_mesh::Mesh::new("m");
        let part = AssemblyPart {
            name: "MeshPart".into(),
            solid: valenx_cad::Solid::from_mesh(mesh),
            transform: identity(),
        };
        let err = step_ap203_assembly_writer(&[part], &PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.backend");
    }

    #[test]
    fn renumber_offsets_every_id_token() {
        let data = "#1 = A(#2, #3);\n#2 = B();\n#3 = C(#1);\n";
        let out = renumber_entities(data, 100);
        assert!(out.contains("#101 = A(#102, #103)"));
        assert!(out.contains("#102 = B()"));
        assert!(out.contains("#103 = C(#101)"));
        // Zero offset is a no-op.
        assert_eq!(renumber_entities(data, 0), data);
    }

    #[test]
    fn highest_entity_id_scans_max() {
        assert_eq!(highest_entity_id("#5 = X(#99); #12 = Y();"), 99);
        assert_eq!(highest_entity_id("no ids"), 0);
    }

    #[test]
    fn extract_data_body_strips_markers() {
        let text = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('S'));\nENDSEC;\n\
                    DATA;\n#1 = POINT();\nENDSEC;\nEND-ISO-10303-21;\n";
        let body = extract_data_body(text).unwrap();
        assert!(body.contains("#1 = POINT()"));
        assert!(!body.contains("DATA;"));
        assert!(!body.contains("ENDSEC"));
    }

    #[test]
    fn escape_step_string_doubles_apostrophes() {
        assert_eq!(escape_step_string("O'Brien"), "O''Brien");
        assert_eq!(escape_step_string("line\nbreak"), "line break");
    }

    #[test]
    fn build_hierarchy_emits_products_and_links() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let parts = vec![
            AssemblyPart {
                name: "A".into(),
                solid: cube.clone(),
                transform: identity(),
            },
            AssemblyPart {
                name: "B".into(),
                solid: cube,
                transform: identity(),
            },
        ];
        let h = build_product_hierarchy(&parts, 1000);
        // Root + 2 children = 3 PRODUCT entities.
        assert_eq!(h.matches("= PRODUCT(").count(), 3);
        // 3 PRODUCT_DEFINITION entities.
        assert_eq!(h.matches("PRODUCT_DEFINITION(").count(), 3);
        // 2 usage occurrences for the 2 children.
        assert_eq!(h.matches("NEXT_ASSEMBLY_USAGE_OCCURRENCE(").count(), 2);
    }

    // -----------------------------------------------------------------
    // Round-15 M2b RED→GREEN: STEP AP203 block-comment + Part-21 string
    // injection via part.name. Sister of the round-14 IFC `ifc_str` fix.
    // -----------------------------------------------------------------

    #[test]
    fn part_name_block_comment_close_does_not_inject_entities() {
        // Pre-fix payload: a hostile part.name carries `*/` to close
        // the `/* --- part i: <name> --- */` block comment, then a
        // sibling `#9999=PRODUCT(...)` entity that pollutes the
        // merged STEP DATA section, then `/*` to swallow the trailer.
        // `escape_step_string` and the block-comment writer must both
        // neutralise the payload — verify by checking the merged
        // hierarchy text never contains `*/` inside a comment or a
        // stray `#9999=` directive.
        let parts = vec![AssemblyPart {
            name: "X*/\n#9999=PRODUCT('PWNED',$);\n/*".into(),
            solid: valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap(),
            transform: identity(),
        }];
        let hierarchy = build_product_hierarchy(&parts, 100);
        // The block-comment writer's payload is `/* --- part i: ... ---
        // */\n`; it must contain exactly one `*/` (the legitimate
        // close), not three (legit + payload's own close + payload's
        // re-open didn't terminate). The hierarchy text in this test
        // doesn't include the block-comment header from `step_ap203_
        // assembly_writer`, so check `escape_step_string` directly
        // through the hierarchy.
        let pwned_count = hierarchy.matches("PWNED").count();
        // The hostile string ends up sanitised either to a no-PWNED
        // form or kept as a literal arg inside a single-quoted string.
        // Either way no sibling #9999 entity opens.
        let _ = pwned_count; // not the primary assertion
                             // No fresh #9999= directive at column 0 (or after whitespace).
        assert!(
            !hierarchy.contains("\n#9999="),
            "part.name injection opened sibling #9999= entity: {hierarchy}"
        );
        // No literal newline inside any single-quoted STEP literal.
        // Walk the lines and check each contains a balanced number of
        // `'` chars (zero or even).
        for line in hierarchy.lines() {
            let quote_count = line.matches('\'').count();
            assert_eq!(
                quote_count % 2,
                0,
                "line has unbalanced STEP single-quotes (newline leaked into literal): {line}"
            );
        }
    }

    #[test]
    fn escape_step_string_strips_newlines_and_escapes_backslash() {
        // Round-15 M2b: extend `escape_step_string` to cover backslash
        // (Part-21 control-char escape) and strip embedded newlines.
        // Pre-fix only handled `'` doubling.
        // Newlines: replaced with space (already handled).
        assert!(!escape_step_string("a\nb").contains('\n'));
        assert!(!escape_step_string("a\rb").contains('\r'));
        // Backslash: must be escaped to `\\\\` so the STEP parser
        // doesn't interpret the next char as a control sequence.
        let escaped = escape_step_string("a\\b");
        // Either the backslash is stripped or doubled — both are safe;
        // the unsafe form is a single backslash. Pin no-single-backslash.
        let single_backslashes = escaped
            .chars()
            .enumerate()
            .filter(|(i, c)| {
                *c == '\\'
                    && escaped.chars().nth(i + 1).is_none_or(|next| next != '\\')
                    && (*i == 0 || escaped.chars().nth(i - 1).is_none_or(|prev| prev != '\\'))
            })
            .count();
        assert_eq!(
            single_backslashes, 0,
            "backslash must be escaped (Part-21 control char): {escaped}"
        );
        // Apostrophe doubling preserved (existing behaviour).
        assert_eq!(escape_step_string("O'Brien"), "O''Brien");
    }

    #[test]
    fn block_comment_close_in_part_name_does_not_leak_data() {
        // Direct test of the block-comment writer in the assembly
        // writer's loop body. The `*/` payload in part.name must not
        // close the comment prematurely; verify by writing a real
        // 2-part assembly with one hostile name and inspecting the
        // merged file for any `#9999=` directive landed outside a
        // proper entity definition.
        let safe_part = AssemblyPart {
            name: "OK".into(),
            solid: valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap(),
            transform: identity(),
        };
        let hostile_part = AssemblyPart {
            name: "Bad*/\n#9999=PRODUCT('PWNED',$);\n/*X".into(),
            solid: valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap(),
            transform: identity(),
        };
        let out =
            std::env::temp_dir().join(format!("valenx_asm_round15_{}.step", std::process::id()));
        step_ap203_assembly_writer(&[safe_part, hostile_part], &out).expect("assembly write");
        let text = fs::read_to_string(&out).expect("read back");
        let _ = fs::remove_file(&out);
        // No fresh #9999= directive opened at column 0 (a real STEP
        // entity definition). The hostile #9999= substring still
        // appears inside the single-quoted PRODUCT literal — that's
        // benign data, not a directive.
        assert!(
            !text
                .lines()
                .any(|l| l.starts_with("#9999=") || l.starts_with("#9999 =")),
            "part.name injection opened sibling #9999= directive at column 0: {text}"
        );
        // The merged hierarchy block ends with /* --- AP203 product
        // structure --- */ and each `/* --- part N: ... --- */` block
        // comment from the merger phase. Each block comment must close
        // exactly once on its opening line — count `*/` in the part-
        // header lines.
        for line in text.lines() {
            if line.starts_with("/* --- part ") {
                let close_count = line.matches("*/").count();
                assert_eq!(
                    close_count, 1,
                    "part-block comment line has != 1 close delimiter: {line}"
                );
            }
        }
    }

    #[test]
    fn writes_two_part_assembly_to_disk() {
        // End-to-end: two real BRep parts -> one merged STEP file.
        let base = valenx_cad::box_solid(10.0, 10.0, 2.0).unwrap();
        let pin = valenx_cad::cylinder(1.0, 8.0).unwrap();
        let parts = vec![
            AssemblyPart {
                name: "Base".into(),
                solid: base,
                transform: identity(),
            },
            AssemblyPart {
                name: "Pin".into(),
                solid: pin,
                transform: [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 2.0],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        ];
        let out = std::env::temp_dir().join(format!("valenx_asm_test_{}.step", std::process::id()));
        step_ap203_assembly_writer(&parts, &out).expect("assembly write");
        let text = fs::read_to_string(&out).expect("read back");
        let _ = fs::remove_file(&out);

        // The merged file has a valid STEP envelope.
        assert!(text.starts_with("ISO-10303-21;"));
        assert!(text.trim_end().ends_with("END-ISO-10303-21;"));
        // Both parts contributed geometry.
        assert!(text.contains("/* --- part 0: Base --- */"));
        assert!(text.contains("/* --- part 1: Pin --- */"));
        // Product structure present: the AP203 hierarchy contributes
        // a root + one PRODUCT per child, by name...
        assert!(text.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"));
        assert!(text.contains("PRODUCT('assembly', 'assembly'"));
        assert!(text.contains("PRODUCT('Base', 'Base'"));
        assert!(text.contains("PRODUCT('Pin', 'Pin'"));
        // ...and each merged part body keeps its own truck-stepio
        // PRODUCT envelope, so the file has the 3 hierarchy products
        // plus one per part body (5 for a 2-part assembly).
        assert_eq!(
            text.matches("= PRODUCT(").count(),
            3 + parts.len(),
            "3 hierarchy products + one PRODUCT per merged part body"
        );
        // Entity ids are unique across the merged DATA section.
        let mut ids: Vec<u64> = Vec::new();
        for line in text.lines() {
            let t = line.trim_start();
            if let Some(rest) = t.strip_prefix('#') {
                let d: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
                if let Ok(id) = d.parse::<u64>() {
                    ids.push(id);
                }
            }
        }
        let mut sorted = ids.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), ids.len(), "entity ids must be unique");
    }

    /// Round-20 L3 RED→GREEN: the scratch-file reader rejects a
    /// file larger than the shared CAD-interchange cap. Pre-fix the
    /// scratch read was a bare `fs::read_to_string` — a pathological
    /// writer (or hostile process racing on the predictable
    /// `valenx_asm_part_<pid>_<i>.step` path in `/tmp`) could plant
    /// a multi-GB fake scratch file and OOM the merge step.
    #[test]
    fn read_capped_step_scratch_rejects_oversize_file() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx_asm_r20l3_oversize_{}.step",
            std::process::id()
        ));
        let mut f = std::fs::File::create(&tmp).unwrap();
        // 257 MiB — past the 256 MiB MAX_CAD_INTERCHANGE_FILE_BYTES.
        // Sparse via set_len so the test doesn't actually write 257 MiB.
        f.set_len(MAX_CAD_INTERCHANGE_FILE_BYTES + 1).unwrap();
        f.write_all(b"x").unwrap();
        drop(f);
        let err = read_capped_step_scratch(&tmp)
            .expect_err("round-20 L3: oversized scratch must reject as IO error");
        match err {
            OcctExchangeError::Io(io) => {
                assert_eq!(
                    io.kind(),
                    std::io::ErrorKind::InvalidData,
                    "over-cap scratch file must surface as InvalidData",
                );
                let msg = io.to_string();
                assert!(
                    msg.contains("cap") || msg.contains("exceed"),
                    "expected cap-mention, got: {msg}"
                );
            }
            other => panic!("expected OcctExchangeError::Io, got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}

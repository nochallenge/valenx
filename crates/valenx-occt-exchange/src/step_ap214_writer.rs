//! Phase 103 — STEP AP214 (ISO 10303-214, Core Data for Automotive
//! Mechanical Design Processes) writer.
//!
//! ## What OCCT does
//!
//! `STEPControl_Writer` with `write.step.schema` = `"AP214IS"` emits
//! a STEP file conforming to AP214's automotive-flavoured product
//! data schema. AP214 is a superset of AP203 — same geometric
//! coverage, plus PDM authority blocks (`PERSON_AND_ORGANIZATION`,
//! `APPROVAL_*`, `SECURITY_CLASSIFICATION_LEVEL`), per-face colour
//! attribution, and assembly hierarchy. It's still in active use at
//! European automotive OEMs (VW Group, BMW, PSA).
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 103.5). The geometry is written
//! by the proven [`fn@crate::step_ap203_writer`] backend
//! (truck-stepio). The resulting text is then post-processed so the
//! file is a genuine AP214 file rather than an AP203 look-alike:
//!
//! 1. The `FILE_SCHEMA` line is rewritten to declare the AP214
//!    schema name (`AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }`) —
//!    AP214-aware readers key off exactly this string.
//! 2. Every caller-supplied colour is emitted as a real `COLOUR_RGB`
//!    entity appended to the `DATA` section, with fresh entity ids
//!    that do not collide with truck-stepio's. The colours are
//!    recoverable by [`fn@crate::step_ap214_reader`] (and any STEP
//!    reader's colour scan).
//!
//! The colours are **not** wired into the full
//! `STYLED_ITEM → PRESENTATION_STYLE_ASSIGNMENT → SURFACE_STYLE_USAGE`
//! graph that binds a colour to a specific `ADVANCED_FACE` — that
//! deep per-face binding is Tier-2 fidelity work. The colour data is
//! present and round-trips; it just is not yet face-anchored.

use std::path::Path;

use valenx_cad::Solid;
use valenx_step_iges::ap242::Ap242Color;

use crate::error::OcctExchangeError;

/// AP214 schema identifier string. AP214-aware readers branch on
/// exactly this token in the `FILE_SCHEMA` entity.
const AP214_SCHEMA: &str = "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }";

/// Write `solid` to `path` as ISO 10303-214 STEP text, with per-solid
/// colour attribution.
///
/// The geometry is written via the AP203 backend, then the file is
/// upgraded in place to declare the AP214 schema and to carry the
/// supplied colours as `COLOUR_RGB` entities.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if `path` is not `.step` / `.stp`.
/// - [`OcctExchangeError::Backend`] when the AP203 backend refuses
///   the solid (mesh-backed, empty boundary).
/// - [`OcctExchangeError::Io`] for filesystem failures.
/// - [`OcctExchangeError::Parse`] if the backend produced a STEP file
///   without the expected `FILE_SCHEMA` / `ENDSEC` markers.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use valenx_occt_exchange::step_ap214_writer;
/// let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
/// step_ap214_writer(&cube, &[], &PathBuf::from("cube_ap214.step")).unwrap();
/// ```
pub fn step_ap214_writer(
    solid: &Solid,
    colors: &[Ap242Color],
    path: &Path,
) -> Result<(), OcctExchangeError> {
    // 1. Write the geometry with the proven AP203 backend.
    crate::step_ap203_writer(solid, path)?;

    // 2. Re-read the file and upgrade it to AP214.
    //
    // Round-21 M1: cap the re-read at MAX_CAD_INTERCHANGE_FILE_BYTES
    // (256 MiB). Pre-fix this was a bare `fs::read_to_string` — even
    // though the file was just written by us, a hostile filesystem
    // (or a concurrent writer racing under the same path) could
    // present a multi-GB file before this read.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES as usize,
    )?;
    let upgraded = upgrade_to_ap214(&text, colors)?;
    valenx_core::io_caps::atomic_write_str(path, &upgraded)?;
    Ok(())
}

/// Rewrite STEP text so it declares the AP214 schema and carries the
/// supplied colours as `COLOUR_RGB` entities. Pure function — kept
/// filesystem-free so it is unit-testable.
fn upgrade_to_ap214(text: &str, colors: &[Ap242Color]) -> Result<String, OcctExchangeError> {
    // --- schema rewrite ---
    // truck-stepio emits `FILE_SCHEMA(('ISO-10303-042'));`. Replace
    // the schema token with the AP214 identifier.
    let schema_start = text.find("FILE_SCHEMA").ok_or_else(|| {
        OcctExchangeError::parse("step header", "no FILE_SCHEMA entity in backend output")
    })?;
    let schema_end = text[schema_start..].find(';').ok_or_else(|| {
        OcctExchangeError::parse("step header", "unterminated FILE_SCHEMA entity")
    })? + schema_start
        + 1;
    let mut out = String::with_capacity(text.len() + colors.len() * 96);
    out.push_str(&text[..schema_start]);
    out.push_str(&format!("FILE_SCHEMA(('{AP214_SCHEMA}'));"));
    let body = &text[schema_end..];

    // --- colour injection ---
    if colors.is_empty() {
        out.push_str(body);
        return Ok(out);
    }

    // Find the closing `ENDSEC;` of the DATA section. truck-stepio
    // writes `...ENDSEC;\nEND-ISO-10303-21;`. The colour entities go
    // just before that ENDSEC.
    let end_iso = body.find("END-ISO-10303-21;").ok_or_else(|| {
        OcctExchangeError::parse("step data", "no END-ISO-10303-21 terminator")
    })?;
    // The ENDSEC that closes DATA is the last ENDSEC before END-ISO.
    let data_endsec = body[..end_iso].rfind("ENDSEC;").ok_or_else(|| {
        OcctExchangeError::parse("step data", "no ENDSEC closing the DATA section")
    })?;

    out.push_str(&body[..data_endsec]);
    // Fresh entity ids start above truck-stepio's highest id.
    let mut next_id = highest_entity_id(body) + 1;
    for c in colors {
        out.push_str(&format!(
            "/* valenx-ap214 colour: owner={} */\n\
             #{next_id} = COLOUR_RGB('{}', {:.6}, {:.6}, {:.6});\n",
            c.owner, c.owner, c.r, c.g, c.b,
        ));
        next_id += 1;
    }
    out.push_str(&body[data_endsec..]);
    Ok(out)
}

/// Highest `#N` entity id appearing in the STEP text. Returns 0 when
/// the text has no entities.
fn highest_entity_id(text: &str) -> u64 {
    let mut max = 0u64;
    for line in text.lines() {
        let trimmed = line.trim_start();
        let Some(rest) = trimmed.strip_prefix('#') else {
            continue;
        };
        let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        if let Ok(id) = digits.parse::<u64>() {
            max = max.max(id);
        }
    }
    max
}

#[cfg(test)]
mod tests {
    use super::*;

    fn color(owner: &str, r: f32, g: f32, b: f32) -> Ap242Color {
        Ap242Color {
            owner: owner.to_string(),
            r,
            g,
            b,
        }
    }

    #[test]
    fn upgrade_rewrites_schema() {
        let input = "ISO-10303-21;\nHEADER;\n\
                     FILE_SCHEMA(('ISO-10303-042'));\nENDSEC;\n\
                     DATA;\n#1 = CARTESIAN_POINT('', (0.0,0.0,0.0));\n\
                     ENDSEC;\nEND-ISO-10303-21;\n";
        let out = upgrade_to_ap214(input, &[]).unwrap();
        assert!(out.contains("AUTOMOTIVE_DESIGN"));
        assert!(!out.contains("ISO-10303-042"));
        // Geometry untouched.
        assert!(out.contains("#1 = CARTESIAN_POINT"));
    }

    #[test]
    fn upgrade_injects_colour_entities() {
        let input = "ISO-10303-21;\nHEADER;\n\
                     FILE_SCHEMA(('ISO-10303-042'));\nENDSEC;\n\
                     DATA;\n#7 = CARTESIAN_POINT('', (0.0,0.0,0.0));\n\
                     ENDSEC;\nEND-ISO-10303-21;\n";
        let out = upgrade_to_ap214(
            input,
            &[color("solid", 1.0, 0.5, 0.25), color("face_2", 0.0, 0.0, 1.0)],
        )
        .unwrap();
        assert!(out.contains("COLOUR_RGB('solid', 1.000000, 0.500000, 0.250000)"));
        assert!(out.contains("COLOUR_RGB('face_2', 0.000000, 0.000000, 1.000000)"));
        // New ids must not collide with #7.
        assert!(out.contains("#8 = COLOUR_RGB"));
        assert!(out.contains("#9 = COLOUR_RGB"));
        // Colours sit inside the DATA section, before the terminator.
        let endsec = out.rfind("ENDSEC;").unwrap();
        let colour = out.find("COLOUR_RGB").unwrap();
        assert!(colour < endsec, "colours must precede the closing ENDSEC");
    }

    #[test]
    fn highest_entity_id_finds_max() {
        let text = "#1 = A;\n#42 = B;\n#7 = C;\nnot an entity\n";
        assert_eq!(highest_entity_id(text), 42);
        assert_eq!(highest_entity_id("no entities here"), 0);
    }

    #[test]
    fn upgrade_rejects_text_without_schema() {
        let err = upgrade_to_ap214("garbage with no header", &[]).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.parse");
    }

    #[test]
    fn rejects_non_step_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = step_ap214_writer(&cube, &[], std::path::Path::new("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}

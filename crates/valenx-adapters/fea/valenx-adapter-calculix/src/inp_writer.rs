//! Generate a CalculiX `.inp` Abaqus-flavoured input deck from a
//! canonical [`Mesh`] plus a [`LinearStaticInput`].
//!
//! The writer emits a single self-contained `.inp` — no `*INCLUDE`
//! pulling in external mesh files. Users who already have hand-
//! written `.inp` decks drive `ccx` directly and don't need this.

use std::fmt::Write;
use std::fs;
use std::io;
use std::path::Path;

use valenx_core::adapter_helpers::sanitize_structured_identifier;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::case_input::{Boundary, LinearStaticInput, Load, OutputField};

/// Generate the entire `.inp` as a single string.
pub fn generate(mesh: &Mesh, input: &LinearStaticInput, name: &str) -> String {
    let mut out = String::with_capacity(4096);

    writeln!(out, "*HEADING").unwrap();
    // Round-15 M2a: project name and material.name come straight from
    // user-supplied data (case.toml / project shell). Pre-fix they
    // interpolated raw into a CCX card argument, so a payload like
    // `"demo\n*BOUNDARY\n NALL, 1, 3, 0.0\n"` would open a sibling
    // *BOUNDARY card that silently re-fixes DOFs at run time. Sister
    // of the Elmer SIF writer fix.
    let safe_name = sanitize_structured_identifier(name);
    writeln!(out, " valenx/{safe_name} ({:?})", input.analysis).unwrap();
    writeln!(out).unwrap();

    write_nodes(&mut out, mesh);
    writeln!(out).unwrap();
    write_element_blocks(&mut out, mesh);
    writeln!(out).unwrap();

    write_node_sets(&mut out, mesh);
    writeln!(out).unwrap();

    write_material(&mut out, input);
    writeln!(out).unwrap();

    // Assign the material to every 3D element block via a SOLID
    // SECTION. 1D / shell elements would want SHELL / BEAM sections
    // when those land.
    write_solid_sections(&mut out, mesh, &input.material.name);
    writeln!(out).unwrap();

    write_boundaries(&mut out, &input.boundaries);
    writeln!(out).unwrap();

    write_step(&mut out, input);

    out
}

fn write_nodes(out: &mut String, mesh: &Mesh) {
    writeln!(out, "*NODE, NSET=NALL").unwrap();
    for (i, n) in mesh.nodes.iter().enumerate() {
        // CalculiX node IDs are 1-based.
        writeln!(out, " {}, {:.10e}, {:.10e}, {:.10e}", i + 1, n.x, n.y, n.z).unwrap();
    }
}

fn write_element_blocks(out: &mut String, mesh: &Mesh) {
    let mut element_tag: u32 = 0;
    for (block_idx, block) in mesh.element_blocks.iter().enumerate() {
        let Some(ccx_type) = canonical_to_ccx(block.element_type) else {
            // Unknown element type — comment it out rather than emit
            // a broken block.
            writeln!(
                out,
                "** skipped {} elements of type {:?} — no CalculiX mapping",
                block.count(),
                block.element_type
            )
            .unwrap();
            element_tag += block.count() as u32;
            continue;
        };
        let elset = format!("EBLOCK_{block_idx}");
        writeln!(out, "*ELEMENT, TYPE={ccx_type}, ELSET={elset}").unwrap();
        let npe = block.element_type.nodes_per_element();
        for i in 0..block.count() {
            element_tag += 1;
            let start = i * npe;
            let nodes = &block.connectivity[start..start + npe];
            // CalculiX expects 1-based node IDs separated by commas.
            write!(out, " {element_tag}").unwrap();
            for &n in nodes {
                write!(out, ", {}", n + 1).unwrap();
            }
            writeln!(out).unwrap();
        }
    }
}

fn write_node_sets(out: &mut String, mesh: &Mesh) {
    for group in &mesh.boundaries {
        if group.element_indices.is_empty() {
            continue;
        }
        // Translate the boundary group's element indices into node
        // ids. For `Nodes`-kind groups the `element_indices` field
        // is reused as direct node indices; for `Faces` we walk the
        // element blocks to gather unique node ids. MVP: treat
        // everything as node ids and let users bring the right
        // groups. Faces-to-nodes lowering lands with a richer
        // `BoundaryGroup` type.
        let nset = sanitize_nset_name(&group.name);
        writeln!(out, "*NSET, NSET={nset}").unwrap();
        // Eight per line is CalculiX convention.
        let mut written = 0;
        for &id in &group.element_indices {
            // 1-based
            if written > 0 && written % 8 == 0 {
                writeln!(out).unwrap();
            } else if written > 0 {
                write!(out, ", ").unwrap();
            }
            write!(out, "{}", id + 1).unwrap();
            written += 1;
        }
        if written > 0 {
            writeln!(out).unwrap();
        }
    }
}

fn write_material(out: &mut String, input: &LinearStaticInput) {
    // Round-15 M2a: see `generate()`. material.name flows from
    // case.toml unchecked.
    let mat_name = sanitize_structured_identifier(&input.material.name);
    writeln!(out, "*MATERIAL, NAME={mat_name}").unwrap();
    writeln!(out, "*ELASTIC").unwrap();
    writeln!(
        out,
        " {:.6e}, {:.6}",
        input.material.youngs_modulus, input.material.poissons_ratio
    )
    .unwrap();
    if let Some(rho) = input.material.density {
        writeln!(out, "*DENSITY").unwrap();
        writeln!(out, " {rho:.6e}").unwrap();
    }
}

fn write_solid_sections(out: &mut String, mesh: &Mesh, material: &str) {
    // Round-15 M2a sister site: material name interpolates into the
    // MATERIAL= field of *SOLID SECTION, so the same payload that
    // attacks *MATERIAL, NAME= works here too.
    let material = sanitize_structured_identifier(material);
    for (block_idx, block) in mesh.element_blocks.iter().enumerate() {
        if block.element_type.dim() != 3 {
            continue;
        }
        writeln!(
            out,
            "*SOLID SECTION, ELSET=EBLOCK_{block_idx}, MATERIAL={material}"
        )
        .unwrap();
    }
}

fn write_boundaries(out: &mut String, bcs: &[Boundary]) {
    if bcs.is_empty() {
        return;
    }
    writeln!(out, "*BOUNDARY").unwrap();
    for bc in bcs {
        writeln!(
            out,
            " {}, {}, {}, {:.6e}",
            sanitize_nset_name(&bc.nset),
            bc.dof_start,
            bc.dof_end,
            bc.value
        )
        .unwrap();
    }
}

fn write_step(out: &mut String, input: &LinearStaticInput) {
    // NLGEOM is only meaningful for the structural analyses
    // (*STATIC and *DYNAMIC). Thermal + Modal silently ignore the
    // flag — CalculiX would error on `*STEP, NLGEOM` followed by
    // `*HEAT TRANSFER` because heat conduction is inherently
    // small-deformation in CCX's formulation.
    let nlgeom_applies = matches!(
        input.analysis,
        crate::case_input::AnalysisKind::LinearStatic
            | crate::case_input::AnalysisKind::LinearDynamic
    );
    if input.step.nlgeom && nlgeom_applies {
        writeln!(out, "*STEP, NLGEOM").unwrap();
    } else {
        writeln!(out, "*STEP").unwrap();
    }
    writeln!(out, "{}", input.analysis.ccx_card()).unwrap();
    // Every analysis that consumes a `<delta_t>, <total_t>` data row
    // gets one — *STATIC + *DYNAMIC + transient *HEAT TRANSFER.
    // *FREQUENCY and steady *HEAT TRANSFER skip the line; CalculiX
    // would otherwise read it as an eigen request count or a property
    // it can't parse, depending on the card.
    if input.analysis.needs_increment_line() {
        // Under NLGEOM the data row optionally extends to a 4-tuple:
        // <inc_init>, <time_total>, <inc_min>, <inc_max>. Emit the
        // longer form when the user specifies bounds; otherwise the
        // 2-tuple is enough and CalculiX picks defaults.
        match (input.step.inc_min, input.step.inc_max) {
            (Some(lo), Some(hi)) if input.step.nlgeom && nlgeom_applies => {
                writeln!(
                    out,
                    " {:.6e}, {:.6e}, {:.6e}, {:.6e}",
                    input.step.time_increment, input.step.time_total, lo, hi
                )
                .unwrap();
            }
            _ => {
                writeln!(
                    out,
                    " {:.6e}, {:.6e}",
                    input.step.time_increment, input.step.time_total
                )
                .unwrap();
            }
        }
    }

    write_loads(out, &input.loads);
    write_output_requests(out, &input.step.output_fields);

    writeln!(out, "*END STEP").unwrap();
}

fn write_loads(out: &mut String, loads: &[Load]) {
    if loads.is_empty() {
        return;
    }
    writeln!(out, "*CLOAD").unwrap();
    for load in loads {
        writeln!(
            out,
            " {}, {}, {:.6e}",
            sanitize_nset_name(&load.nset),
            load.dof,
            load.force
        )
        .unwrap();
    }
}

fn write_output_requests(out: &mut String, fields: &[OutputField]) {
    let node_fields: Vec<&OutputField> = fields
        .iter()
        .filter(|f| matches!(f, OutputField::U | OutputField::Rf | OutputField::Nt))
        .collect();
    let element_fields: Vec<&OutputField> = fields
        .iter()
        .filter(|f| matches!(f, OutputField::S))
        .collect();

    if !node_fields.is_empty() {
        writeln!(out, "*NODE PRINT, NSET=NALL").unwrap();
        let codes: Vec<&str> = node_fields.iter().map(|f| f.ccx_code()).collect();
        writeln!(out, " {}", codes.join(", ")).unwrap();
        writeln!(out, "*NODE FILE").unwrap();
        writeln!(out, " {}", codes.join(", ")).unwrap();
    }
    if !element_fields.is_empty() {
        writeln!(out, "*EL PRINT, ELSET=EALL").unwrap();
        let codes: Vec<&str> = element_fields.iter().map(|f| f.ccx_code()).collect();
        writeln!(out, " {}", codes.join(", ")).unwrap();
        writeln!(out, "*EL FILE").unwrap();
        writeln!(out, " {}", codes.join(", ")).unwrap();
    }
}

/// Map canonical `ElementType` to CalculiX / Abaqus element code.
/// Unknown types return `None` and the element block is skipped
/// with a `**` comment in the output.
fn canonical_to_ccx(et: ElementType) -> Option<&'static str> {
    match et {
        ElementType::Tet4 => Some("C3D4"),
        ElementType::Tet10 => Some("C3D10"),
        ElementType::Hex8 => Some("C3D8"),
        ElementType::Hex20 => Some("C3D20"),
        ElementType::Prism6 => Some("C3D6"),
        ElementType::Pyr5 => Some("C3D5"),
        ElementType::Tri3 => Some("CPS3"),
        ElementType::Tri6 => Some("CPS6"),
        ElementType::Quad4 => Some("CPS4"),
        ElementType::Line2 => None, // beams need section props first
    }
}

fn sanitize_nset_name(name: &str) -> String {
    // CalculiX allows letters, digits, underscore, hyphen; all upper.
    let mut s = String::with_capacity(name.len());
    for c in name.chars() {
        if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
            s.push(c.to_ascii_uppercase());
        } else {
            s.push('_');
        }
    }
    if s.is_empty() {
        "UNNAMED".to_string()
    } else {
        s
    }
}

/// Write the generated deck into a file, creating parents.
pub fn write_to_file(
    mesh: &Mesh,
    input: &LinearStaticInput,
    name: &str,
    path: &Path,
) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    valenx_core::io_caps::atomic_write_str(path, &generate(mesh, input, name))
}

// Suppress the unused-helper hint on `write_element_blocks`' loop
// variable when the compiler is being chatty.
#[allow(dead_code)]
fn _ref_element_block(_b: &ElementBlock) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::case_input::{AnalysisKind, LinearStaticInput, Material, OutputField, Step};
    use nalgebra::Vector3;
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    fn one_tet_mesh() -> Mesh {
        let mut m = Mesh::new("cantilever");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        m.recompute_stats();
        m
    }

    fn sample_input() -> LinearStaticInput {
        LinearStaticInput {
            analysis: crate::case_input::AnalysisKind::LinearStatic,
            material: Material {
                name: "steel".into(),
                youngs_modulus: 210e9,
                poissons_ratio: 0.3,
                density: Some(7850.0),
            },
            mesh_source: std::path::PathBuf::from("mesh.canonical.json"),
            boundaries: vec![Boundary {
                nset: "fixed".into(),
                dof_start: 1,
                dof_end: 3,
                value: 0.0,
            }],
            loads: vec![Load {
                nset: "tip".into(),
                dof: 2,
                force: -1000.0,
            }],
            step: Step {
                time_total: 1.0,
                time_increment: 1.0,
                output_fields: vec![OutputField::U, OutputField::S, OutputField::Rf],
                nlgeom: false,
                inc_min: None,
                inc_max: None,
            },
        }
    }

    #[test]
    fn inp_has_expected_sections() {
        let text = generate(&one_tet_mesh(), &sample_input(), "cantilever");
        assert!(text.contains("*HEADING"));
        assert!(text.contains("*NODE, NSET=NALL"));
        assert!(text.contains("*ELEMENT, TYPE=C3D4"));
        assert!(text.contains("*MATERIAL, NAME=steel"));
        assert!(text.contains("*ELASTIC"));
        assert!(text.contains("*DENSITY"));
        assert!(text.contains("*SOLID SECTION"));
        assert!(text.contains("*BOUNDARY"));
        assert!(text.contains("*STEP"));
        assert!(text.contains("*STATIC"));
        assert!(text.contains("*CLOAD"));
        assert!(text.contains("*NODE PRINT"));
        assert!(text.contains("*END STEP"));
    }

    #[test]
    fn nodes_are_one_based() {
        let text = generate(&one_tet_mesh(), &sample_input(), "c");
        assert!(text.contains("\n 1, "));
        assert!(text.contains("\n 4, "));
        // 0-index would appear as " 0, " at a line start, which is
        // never valid in CalculiX.
        assert!(!text.contains("\n 0, "));
    }

    #[test]
    fn nset_names_are_sanitised() {
        assert_eq!(sanitize_nset_name("fixed"), "FIXED");
        assert_eq!(sanitize_nset_name("top face"), "TOP_FACE");
        assert_eq!(sanitize_nset_name("inlet.outlet"), "INLET_OUTLET");
        assert_eq!(sanitize_nset_name(""), "UNNAMED");
    }

    #[test]
    fn thermal_step_writes_heat_card() {
        let mut input = sample_input();
        input.analysis = crate::case_input::AnalysisKind::Thermal;
        input.step.output_fields = vec![OutputField::Nt];
        let text = generate(&one_tet_mesh(), &input, "heat");
        assert!(text.contains("*HEAT TRANSFER, STEADY STATE"));
        assert!(text.contains("*NODE PRINT"));
        assert!(text.contains("NT"));
        // Steady analyses don't emit the time-increment data row.
        // The lines under *HEAT TRANSFER, STEADY STATE should be the
        // load + output cards directly, not a `<delta_t>, <total_t>`.
        let after_heat = text
            .split("*HEAT TRANSFER, STEADY STATE")
            .nth(1)
            .expect("split on heat card");
        let next_line = after_heat
            .lines()
            .find(|l| !l.trim().is_empty())
            .unwrap_or("");
        assert!(
            next_line.starts_with("*"),
            "steady thermal must not emit a numeric data row, got: {next_line}"
        );
    }

    #[test]
    fn dynamic_step_writes_dynamic_card_and_increment_line() {
        let mut input = sample_input();
        input.analysis = crate::case_input::AnalysisKind::LinearDynamic;
        input.step.time_total = 0.01;
        input.step.time_increment = 1e-4;
        let text = generate(&one_tet_mesh(), &input, "drop-test");
        assert!(text.contains("*DYNAMIC"), "got: {text}");
        // The data row right after *DYNAMIC carries delta_t, total_t.
        // Format is scientific (`{:.6e}`) so we look for both values.
        assert!(
            text.contains("1.000000e-4") || text.contains("1.000000e-04"),
            "delta_t (1e-4) missing from inp: {text}"
        );
        assert!(
            text.contains("1.000000e-2") || text.contains("1.000000e-02"),
            "time_total (1e-2) missing from inp: {text}"
        );
        // Static card must NOT also appear — we only emit the picked
        // analysis kind.
        assert!(
            !text.contains("*STATIC\n"),
            "static card leaked into dynamic step"
        );
    }

    #[test]
    fn thermal_transient_writes_heat_card_without_steady_qualifier() {
        let mut input = sample_input();
        input.analysis = crate::case_input::AnalysisKind::ThermalTransient;
        input.step.output_fields = vec![OutputField::Nt];
        input.step.time_total = 60.0;
        input.step.time_increment = 0.1;
        let text = generate(&one_tet_mesh(), &input, "cooldown");
        assert!(text.contains("*HEAT TRANSFER\n"), "got: {text}");
        assert!(
            !text.contains("*HEAT TRANSFER, STEADY STATE"),
            "transient thermal must not emit the steady-state qualifier"
        );
        assert!(text.contains("1.000000e-1") || text.contains("1.000000e-01"));
        assert!(text.contains("6.000000e1") || text.contains("6.000000e01"));
    }

    #[test]
    fn unknown_element_type_is_skipped_not_broken() {
        let mut m = one_tet_mesh();
        let mut line = ElementBlock::new(ElementType::Line2);
        line.connectivity = vec![0, 1];
        m.element_blocks.push(line);
        let text = generate(&m, &sample_input(), "c");
        assert!(text.contains("skipped"));
        // The valid tet block still landed.
        assert!(text.contains("*ELEMENT, TYPE=C3D4"));
    }

    // -----------------------------------------------------------------
    // NLGEOM nonlinear-static
    // -----------------------------------------------------------------

    #[test]
    fn nlgeom_flag_emits_step_with_nlgeom_keyword() {
        let mut input = sample_input();
        input.step.nlgeom = true;
        let text = generate(&one_tet_mesh(), &input, "c");
        assert!(text.contains("*STEP, NLGEOM"), "got:\n{text}");
        // The bare *STEP without keyword must not appear when NLGEOM
        // is set — that would imply we wrote two step cards.
        assert!(!text.contains("*STEP\n"), "got:\n{text}");
    }

    #[test]
    fn nlgeom_unset_keeps_plain_step_card_for_backward_compat() {
        // Existing linear-static cases (default `nlgeom = false`)
        // must keep emitting the plain `*STEP` card so prior runs
        // stay byte-identical.
        let text = generate(&one_tet_mesh(), &sample_input(), "c");
        assert!(text.contains("*STEP\n"), "got:\n{text}");
        assert!(!text.contains("*STEP, NLGEOM"), "got:\n{text}");
    }

    #[test]
    fn nlgeom_emits_4tuple_data_row_when_inc_bounds_set() {
        // Under NLGEOM with both inc_min + inc_max specified, the
        // increment-control data row extends from
        //   <delta_t>, <total_t>
        // to
        //   <delta_t>, <total_t>, <inc_min>, <inc_max>
        let mut input = sample_input();
        input.step.nlgeom = true;
        input.step.time_increment = 0.1;
        input.step.time_total = 1.0;
        input.step.inc_min = Some(1e-6);
        input.step.inc_max = Some(0.2);
        let text = generate(&one_tet_mesh(), &input, "c");
        // 4 floats separated by commas.
        let row = text
            .lines()
            .find(|line| line.contains("1.000000e-1, 1.000000e0"))
            .expect("expected control row in output");
        let comma_count = row.chars().filter(|&c| c == ',').count();
        assert_eq!(comma_count, 3, "row should have 4 fields (3 commas): {row}");
    }

    #[test]
    fn nlgeom_falls_back_to_2tuple_when_only_one_bound_given() {
        // CalculiX expects either zero or both bounds; emitting just
        // inc_min would silently mis-parse. Stay conservative: drop
        // back to the 2-tuple data row.
        let mut input = sample_input();
        input.step.nlgeom = true;
        input.step.inc_min = Some(1e-6);
        input.step.inc_max = None;
        let text = generate(&one_tet_mesh(), &input, "c");
        let row = text
            .lines()
            .find(|line| line.contains("1.000000e0, 1.000000e0"))
            .expect("expected control row");
        let comma_count = row.chars().filter(|&c| c == ',').count();
        assert_eq!(comma_count, 1, "row should have 2 fields (1 comma): {row}");
    }

    #[test]
    fn nlgeom_is_silently_ignored_for_thermal_analyses() {
        // Heat-transfer cases shouldn't see NLGEOM — CalculiX would
        // refuse the `*HEAT TRANSFER` card paired with a `NLGEOM`
        // step keyword.
        let mut input = sample_input();
        input.analysis = AnalysisKind::Thermal;
        input.step.nlgeom = true;
        input.step.output_fields = vec![OutputField::Nt];
        let text = generate(&one_tet_mesh(), &input, "c");
        assert!(!text.contains("NLGEOM"), "got:\n{text}");
        assert!(text.contains("*HEAT TRANSFER, STEADY STATE"));
    }

    #[test]
    fn nlgeom_applies_to_dynamic_structural_analyses() {
        // Geometric nonlinearity is meaningful for *DYNAMIC too —
        // large rotations during a drop test are the canonical
        // example.
        let mut input = sample_input();
        input.analysis = AnalysisKind::LinearDynamic;
        input.step.nlgeom = true;
        let text = generate(&one_tet_mesh(), &input, "c");
        assert!(text.contains("*STEP, NLGEOM"), "got:\n{text}");
        assert!(text.contains("*DYNAMIC"));
    }

    // -----------------------------------------------------------------
    // Round-15 M2a RED→GREEN: INP identifier-injection via material.name
    // (*MATERIAL, NAME=...) and project name (*HEADING). Sister of
    // the SIF writer fix in valenx-adapter-elmer.
    // -----------------------------------------------------------------

    #[test]
    fn inp_material_name_injection_attempt_is_neutralised() {
        // Pre-fix payload: a hostile material.name lets case.toml
        // inject a sibling *MATERIAL or *BOUNDARY card. The pre-fix
        // writer interpolated material.name raw into
        // `*MATERIAL, NAME={}`, so a `"X\n*MATERIAL, NAME=X2\n..."`
        // payload would land as a fresh CCX card at column 0.
        // After sanitisation the payload collapses into a single
        // ugly material name on one line — no fresh cards open.
        let mut input = sample_input();
        // Strip the sample_input boundaries so the only *BOUNDARY in
        // the file would be one the injection opens (sample_input
        // ships a real fixed BC; we need a clean canvas).
        input.boundaries.clear();
        input.material.name = "X\n*MATERIAL, NAME=X2\n*BOUNDARY\n NALL, 1, 3, 0.0".into();
        let text = generate(&one_tet_mesh(), &input, "c");
        // No fresh *BOUNDARY card opened at column 0 by the injection
        // (no real boundaries in this case so any *BOUNDARY would be
        // proof of injection).
        assert!(
            !text.lines().any(|l| l == "*BOUNDARY"),
            "material.name injection opened sibling *BOUNDARY card: {text}"
        );
        // Exactly one *MATERIAL card despite the payload trying to
        // open a second one.
        assert_eq!(
            text.lines().filter(|l| l.starts_with("*MATERIAL")).count(),
            1,
            "material.name injection cloned the *MATERIAL card: {text}"
        );
    }

    #[test]
    fn inp_heading_name_injection_attempt_is_neutralised() {
        // *HEADING / project name injection: a payload with
        // `\n*BOUNDARY\n...` would otherwise emit a sibling card right
        // after the heading line. The sanitiser strips newlines so
        // the *HEADING block stays single-line.
        let project_name =
            "demo\n*BOUNDARY\n NALL, 1, 3, 0.0\n*PIRATED";
        let text = generate(&one_tet_mesh(), &sample_input(), project_name);
        // Walk the lines and assert the *HEADING block is exactly one
        // header line + one content line.
        let mut iter = text.lines();
        let h = iter.next().unwrap();
        assert_eq!(h, "*HEADING", "first line must be *HEADING");
        let content = iter.next().unwrap();
        assert!(
            !content.starts_with('*'),
            "heading content must not start a fresh CCX card: {content}"
        );
        // No *PIRATED card landed anywhere.
        assert!(
            !text.contains("\n*PIRATED"),
            "project-name injection leaked sibling card: {text}"
        );
    }
}

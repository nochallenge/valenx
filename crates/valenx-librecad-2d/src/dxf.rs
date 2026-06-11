//! DXF reader + writer.
//!
//! We support the AutoCAD R12 ASCII DXF format — group-code lines
//! followed by value lines. Supported entities: LINE, CIRCLE, ARC,
//! LWPOLYLINE, SPLINE, TEXT, MTEXT, DIMENSION, HATCH, INSERT. Layers
//! are read from the LAYER section; blocks from BLOCKS.
//!
//! The reader is intentionally permissive — unknown group codes are
//! skipped — so DXF files produced by other tools but containing the
//! supported entity set still round-trip cleanly.

use crate::drawing::{Block, Drawing2D, Entity2D, Layer};
use crate::error::LibreCadError;

/// Read every DXF section we know about; ignore the rest.
///
/// Round-23 sweep: bounded at [`valenx_core::io_caps::MAX_DXF_FILE_BYTES`]
/// (1 GiB) — production multi-layer DXFs can cross 100 MiB; 1 GiB
/// refuses pathological inputs that would OOM `String::from_utf8`.
pub fn read_full(path: &str) -> Result<Drawing2D, LibreCadError> {
    let text = valenx_core::io_caps::read_capped_to_string(
        std::path::Path::new(path),
        valenx_core::io_caps::MAX_DXF_FILE_BYTES as usize,
    )
    .map_err(|e| LibreCadError::Io(e.to_string()))?;
    parse(&text)
}

/// Parse a DXF string. Used by `read_full` and the round-trip tests.
pub fn parse(text: &str) -> Result<Drawing2D, LibreCadError> {
    let pairs = lex(text)?;
    let mut drawing = Drawing2D::new();
    let mut i = 0;
    while i < pairs.len() {
        let (code, val) = &pairs[i];
        if *code == 0 && val == "SECTION" {
            // Section header: next is a `2` code with the section name.
            if i + 1 < pairs.len() && pairs[i + 1].0 == 2 {
                let section = &pairs[i + 1].1;
                match section.as_str() {
                    "TABLES" => i = read_tables(&pairs, i + 2, &mut drawing)?,
                    "BLOCKS" => i = read_blocks(&pairs, i + 2, &mut drawing)?,
                    "ENTITIES" => i = read_entities(&pairs, i + 2, &mut drawing)?,
                    _ => i = skip_to_endsec(&pairs, i + 2),
                }
                continue;
            }
        }
        i += 1;
    }
    Ok(drawing)
}

fn skip_to_endsec(pairs: &[(i32, String)], start: usize) -> usize {
    let mut i = start;
    while i < pairs.len() {
        if pairs[i].0 == 0 && pairs[i].1 == "ENDSEC" {
            return i + 1;
        }
        i += 1;
    }
    i
}

fn read_tables(
    pairs: &[(i32, String)],
    start: usize,
    d: &mut Drawing2D,
) -> Result<usize, LibreCadError> {
    let mut i = start;
    while i < pairs.len() {
        if pairs[i].0 == 0 && pairs[i].1 == "ENDSEC" {
            return Ok(i + 1);
        }
        if pairs[i].0 == 0 && pairs[i].1 == "LAYER" {
            // Read until next 0-code or ENDTAB.
            let (layer, next) = read_layer(pairs, i + 1)?;
            if !d.layers.iter().any(|l| l.name == layer.name) {
                d.layers.push(layer);
            }
            i = next;
            continue;
        }
        i += 1;
    }
    Ok(i)
}

fn read_layer(pairs: &[(i32, String)], start: usize) -> Result<(Layer, usize), LibreCadError> {
    let mut name = String::new();
    let mut color = 7u8;
    let mut linetype = "CONTINUOUS".to_string();
    let mut i = start;
    while i < pairs.len() {
        let (code, val) = &pairs[i];
        if *code == 0 {
            break;
        }
        match *code {
            2 => name = val.clone(),
            62 => {
                color = val.parse::<i32>().unwrap_or(7).clamp(0, 255) as u8;
            }
            6 => linetype = val.clone(),
            _ => {}
        }
        i += 1;
    }
    Ok((
        Layer {
            name,
            color,
            linetype,
            visible: true,
        },
        i,
    ))
}

fn read_blocks(
    pairs: &[(i32, String)],
    start: usize,
    d: &mut Drawing2D,
) -> Result<usize, LibreCadError> {
    let mut i = start;
    while i < pairs.len() {
        if pairs[i].0 == 0 && pairs[i].1 == "ENDSEC" {
            return Ok(i + 1);
        }
        if pairs[i].0 == 0 && pairs[i].1 == "BLOCK" {
            let (block, next) = read_block(pairs, i + 1)?;
            d.blocks.push(block);
            i = next;
            continue;
        }
        i += 1;
    }
    Ok(i)
}

fn read_block(pairs: &[(i32, String)], start: usize) -> Result<(Block, usize), LibreCadError> {
    let mut name = String::new();
    let mut origin = [0.0_f64; 2];
    let mut entities = Vec::new();
    let mut i = start;
    while i < pairs.len() {
        let (code, val) = &pairs[i];
        if *code == 0 && val == "ENDBLK" {
            // Consume any 5 / 100 codes following ENDBLK, then break.
            i += 1;
            while i < pairs.len() && pairs[i].0 != 0 {
                i += 1;
            }
            break;
        }
        if *code == 0 && val != "BLOCK" {
            let kind = val.clone();
            let (entity, next) = read_entity(&kind, pairs, i + 1)?;
            if let Some(e) = entity {
                entities.push(e);
            }
            i = next;
            continue;
        }
        match *code {
            2 => name = val.clone(),
            10 => origin[0] = val.parse().unwrap_or(0.0),
            20 => origin[1] = val.parse().unwrap_or(0.0),
            _ => {}
        }
        i += 1;
    }
    Ok((
        Block {
            name,
            origin,
            entities,
        },
        i,
    ))
}

fn read_entities(
    pairs: &[(i32, String)],
    start: usize,
    d: &mut Drawing2D,
) -> Result<usize, LibreCadError> {
    let mut i = start;
    while i < pairs.len() {
        if pairs[i].0 == 0 && pairs[i].1 == "ENDSEC" {
            return Ok(i + 1);
        }
        if pairs[i].0 == 0 {
            let kind = pairs[i].1.clone();
            let (entity, next) = read_entity(&kind, pairs, i + 1)?;
            if let Some(e) = entity {
                d.entities.push(e);
            }
            i = next;
            continue;
        }
        i += 1;
    }
    Ok(i)
}

fn read_entity(
    kind: &str,
    pairs: &[(i32, String)],
    start: usize,
) -> Result<(Option<Entity2D>, usize), LibreCadError> {
    let mut i = start;
    let mut layer = "0".to_string();
    let mut a = [0.0_f64; 2];
    let mut b = [0.0_f64; 2];
    let mut centre = [0.0_f64; 2];
    let mut text_pos = [0.0_f64; 2];
    let mut radius = 0.0_f64;
    let mut start_angle = 0.0_f64;
    let mut end_angle = 0.0_f64;
    let mut text = String::new();
    let mut height = 1.0_f64;
    let mut width = 1.0_f64;
    let mut block_name = String::new();
    let mut scale = 1.0_f64;
    let mut rotation = 0.0_f64;
    let mut vertices: Vec<[f64; 2]> = Vec::new();
    let mut closed = false;
    let mut current_vtx_x = 0.0_f64;
    let mut current_vtx_y = 0.0_f64;
    let mut have_vtx_x = false;
    let mut have_vtx_y = false;
    let mut degree = 3u8;
    let mut pattern = "SOLID".to_string();
    while i < pairs.len() && pairs[i].0 != 0 {
        let (code, val) = &pairs[i];
        match *code {
            8 => layer = val.clone(),
            10 => {
                a[0] = val.parse().unwrap_or(0.0);
                centre[0] = a[0];
                current_vtx_x = a[0];
                have_vtx_x = true;
            }
            20 => {
                a[1] = val.parse().unwrap_or(0.0);
                centre[1] = a[1];
                current_vtx_y = a[1];
                have_vtx_y = true;
            }
            11 => b[0] = val.parse().unwrap_or(0.0),
            21 => b[1] = val.parse().unwrap_or(0.0),
            // 13/23 carry a DIMENSION's text midpoint (`text_pos`); kept
            // distinct from a's 10/20 so the second point can't clobber it.
            13 => text_pos[0] = val.parse().unwrap_or(0.0),
            23 => text_pos[1] = val.parse().unwrap_or(0.0),
            40 => {
                radius = val.parse().unwrap_or(0.0);
                height = radius;
                scale = radius;
            }
            41 => width = val.parse().unwrap_or(width),
            50 => {
                start_angle = val.parse().unwrap_or(0.0);
                rotation = start_angle;
            }
            51 => end_angle = val.parse().unwrap_or(0.0),
            1 => text = val.clone(),
            2 => block_name = val.clone(),
            70 => {
                let flag: i32 = val.parse().unwrap_or(0);
                closed = (flag & 1) != 0;
                degree = (flag.unsigned_abs() as u8).clamp(1, 7);
            }
            71 => {
                pattern = val.clone();
            }
            _ => {}
        }
        // LWPOLYLINE flushes a vertex whenever it has both X and Y.
        if kind == "LWPOLYLINE" && have_vtx_x && have_vtx_y && (*code == 20 || *code == 21) {
            vertices.push([current_vtx_x, current_vtx_y]);
            have_vtx_x = false;
            have_vtx_y = false;
        }
        i += 1;
    }
    let entity = match kind {
        "LINE" => Some(Entity2D::Line { layer, a, b }),
        "CIRCLE" => Some(Entity2D::Circle {
            layer,
            centre,
            radius,
        }),
        "ARC" => Some(Entity2D::Arc {
            layer,
            centre,
            radius,
            start_angle_deg: start_angle,
            end_angle_deg: end_angle,
        }),
        "LWPOLYLINE" => Some(Entity2D::Polyline {
            layer,
            vertices,
            closed,
        }),
        "SPLINE" => Some(Entity2D::Spline {
            layer,
            control_points: Vec::new(), // v1: no per-CP parsing yet
            degree,
        }),
        "HATCH" => Some(Entity2D::Hatch {
            layer,
            boundary: Vec::new(), // v1: no boundary parsing yet
            pattern,
        }),
        "TEXT" => Some(Entity2D::Text {
            layer,
            position: a,
            height,
            text,
        }),
        "MTEXT" => Some(Entity2D::MText {
            layer,
            position: a,
            height,
            width,
            text,
        }),
        "DIMENSION" => Some(Entity2D::Dimension {
            layer,
            a,
            b,
            text_pos,
            text,
        }),
        "INSERT" => Some(Entity2D::Insert {
            layer,
            block: block_name,
            position: a,
            scale,
            rotation_deg: rotation,
        }),
        _ => None,
    };
    Ok((entity, i))
}

fn lex(text: &str) -> Result<Vec<(i32, String)>, LibreCadError> {
    let mut out = Vec::new();
    let mut lines = text.lines().enumerate();
    while let Some((idx, code_line)) = lines.next() {
        let code_line = code_line.trim();
        if code_line.is_empty() {
            continue;
        }
        let code: i32 = code_line.parse().map_err(|_| LibreCadError::DxfParse {
            line: idx + 1,
            message: format!("expected integer group code, got `{code_line}`"),
        })?;
        let (vidx, value_line) = lines.next().ok_or_else(|| LibreCadError::DxfParse {
            line: idx + 1,
            message: "unexpected EOF after group code".into(),
        })?;
        let _ = vidx;
        out.push((code, value_line.trim().to_string()));
    }
    Ok(out)
}

/// Write a [`Drawing2D`] back out as DXF.
pub fn write_full(d: &Drawing2D, path: &str) -> Result<(), LibreCadError> {
    let text = serialise(d);
    valenx_core::io_caps::atomic_write_str(std::path::Path::new(path), &text)
        .map_err(|e| LibreCadError::Io(e.to_string()))?;
    Ok(())
}

/// Serialise a [`Drawing2D`] to DXF text. Used by both `write_full`
/// and round-trip tests.
pub fn serialise(d: &Drawing2D) -> String {
    let mut s = String::new();
    s.push_str("0\nSECTION\n2\nTABLES\n");
    for layer in &d.layers {
        s.push_str("0\nLAYER\n");
        push_pair(&mut s, 2, &layer.name);
        push_pair_int(&mut s, 62, layer.color as i32);
        push_pair(&mut s, 6, &layer.linetype);
    }
    s.push_str("0\nENDSEC\n");
    if !d.blocks.is_empty() {
        s.push_str("0\nSECTION\n2\nBLOCKS\n");
        for block in &d.blocks {
            s.push_str("0\nBLOCK\n");
            push_pair(&mut s, 2, &block.name);
            push_pair_f(&mut s, 10, block.origin[0]);
            push_pair_f(&mut s, 20, block.origin[1]);
            for entity in &block.entities {
                serialise_entity(&mut s, entity);
            }
            s.push_str("0\nENDBLK\n");
        }
        s.push_str("0\nENDSEC\n");
    }
    s.push_str("0\nSECTION\n2\nENTITIES\n");
    for entity in &d.entities {
        serialise_entity(&mut s, entity);
    }
    s.push_str("0\nENDSEC\n0\nEOF\n");
    s
}

fn serialise_entity(s: &mut String, e: &Entity2D) {
    s.push_str(&format!("0\n{}\n", e.kind()));
    push_pair(s, 8, e.layer());
    match e {
        Entity2D::Line { a, b, .. } => {
            push_pair_f(s, 10, a[0]);
            push_pair_f(s, 20, a[1]);
            push_pair_f(s, 11, b[0]);
            push_pair_f(s, 21, b[1]);
        }
        Entity2D::Circle { centre, radius, .. } => {
            push_pair_f(s, 10, centre[0]);
            push_pair_f(s, 20, centre[1]);
            push_pair_f(s, 40, *radius);
        }
        Entity2D::Arc {
            centre,
            radius,
            start_angle_deg,
            end_angle_deg,
            ..
        } => {
            push_pair_f(s, 10, centre[0]);
            push_pair_f(s, 20, centre[1]);
            push_pair_f(s, 40, *radius);
            push_pair_f(s, 50, *start_angle_deg);
            push_pair_f(s, 51, *end_angle_deg);
        }
        Entity2D::Polyline {
            vertices, closed, ..
        } => {
            push_pair_int(s, 70, if *closed { 1 } else { 0 });
            for v in vertices {
                push_pair_f(s, 10, v[0]);
                push_pair_f(s, 20, v[1]);
            }
        }
        Entity2D::Spline { degree, .. } => {
            push_pair_int(s, 70, *degree as i32);
        }
        Entity2D::Hatch { pattern, .. } => {
            push_pair(s, 71, pattern);
        }
        Entity2D::Text {
            position,
            height,
            text,
            ..
        } => {
            push_pair_f(s, 10, position[0]);
            push_pair_f(s, 20, position[1]);
            push_pair_f(s, 40, *height);
            push_pair(s, 1, text);
        }
        Entity2D::MText {
            position,
            height,
            width,
            text,
            ..
        } => {
            push_pair_f(s, 10, position[0]);
            push_pair_f(s, 20, position[1]);
            push_pair_f(s, 40, *height);
            push_pair_f(s, 41, *width);
            push_pair(s, 1, text);
        }
        Entity2D::Dimension {
            a, b, text_pos, text, ..
        } => {
            push_pair_f(s, 10, a[0]);
            push_pair_f(s, 20, a[1]);
            push_pair_f(s, 11, b[0]);
            push_pair_f(s, 21, b[1]);
            push_pair_f(s, 13, text_pos[0]);
            push_pair_f(s, 23, text_pos[1]);
            push_pair(s, 1, text);
        }
        Entity2D::Insert {
            block,
            position,
            scale,
            rotation_deg,
            ..
        } => {
            push_pair(s, 2, block);
            push_pair_f(s, 10, position[0]);
            push_pair_f(s, 20, position[1]);
            push_pair_f(s, 40, *scale);
            push_pair_f(s, 50, *rotation_deg);
        }
    }
}

fn push_pair(s: &mut String, code: i32, val: &str) {
    s.push_str(&format!("{code}\n{val}\n"));
}

fn push_pair_f(s: &mut String, code: i32, val: f64) {
    s.push_str(&format!("{code}\n{val}\n"));
}

fn push_pair_int(s: &mut String, code: i32, val: i32) {
    s.push_str(&format!("{code}\n{val}\n"));
}

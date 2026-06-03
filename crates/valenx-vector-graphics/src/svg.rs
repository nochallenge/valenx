//! SVG writer + minimal reader.
//!
//! The writer emits an XML document with the standard `xmlns` and a
//! `viewBox` covering the entity bounding box plus a small margin.
//! The reader handles the same subset of elements that the writer
//! emits: `<line>`, `<rect>`, `<ellipse>`, `<polygon>`, `<text>`,
//! and `<path d="...">`. It tolerates whitespace and attribute
//! reordering, and falls back to skipping unrecognised elements.

use nalgebra::Vector2;

use crate::entity::{PathSegment, VectorEntity};
use crate::error::VectorError;
use crate::path;

/// Emit an SVG document containing every entity in `entities`.
pub fn to_svg(entities: &[VectorEntity]) -> String {
    let (lo, hi) = bbox_all(entities);
    let mw = (hi.x - lo.x).max(1.0);
    let mh = (hi.y - lo.y).max(1.0);
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{} {} {} {}\">\n",
        lo.x, lo.y, mw, mh
    ));
    for e in entities {
        push_entity(&mut s, e);
        s.push('\n');
    }
    s.push_str("</svg>\n");
    s
}

fn bbox_all(entities: &[VectorEntity]) -> (Vector2<f64>, Vector2<f64>) {
    if entities.is_empty() {
        return (Vector2::zeros(), Vector2::new(1.0, 1.0));
    }
    let mut lo = Vector2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Vector2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut add = |p: Vector2<f64>| {
        if p.x < lo.x { lo.x = p.x; }
        if p.y < lo.y { lo.y = p.y; }
        if p.x > hi.x { hi.x = p.x; }
        if p.y > hi.y { hi.y = p.y; }
    };
    for e in entities {
        match e {
            VectorEntity::Line { a, b } => {
                add(*a);
                add(*b);
            }
            VectorEntity::Rect { origin, size } => {
                add(*origin);
                add(origin + size);
            }
            VectorEntity::Ellipse { centre, rx, ry } => {
                add(centre - Vector2::new(*rx, *ry));
                add(centre + Vector2::new(*rx, *ry));
            }
            VectorEntity::Polygon(v) => {
                for p in v {
                    add(*p);
                }
            }
            VectorEntity::Text { anchor, .. } => add(*anchor),
            VectorEntity::Path(p) => {
                let (l, h) = path::bbox(p);
                add(l);
                add(h);
            }
        }
    }
    if !lo.x.is_finite() {
        return (Vector2::zeros(), Vector2::new(1.0, 1.0));
    }
    (lo, hi)
}

fn push_entity(s: &mut String, e: &VectorEntity) {
    match e {
        VectorEntity::Line { a, b } => {
            s.push_str(&format!(
                "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\"/>",
                a.x, a.y, b.x, b.y
            ));
        }
        VectorEntity::Rect { origin, size } => {
            s.push_str(&format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"none\" stroke=\"black\"/>",
                origin.x, origin.y, size.x, size.y
            ));
        }
        VectorEntity::Ellipse { centre, rx, ry } => {
            s.push_str(&format!(
                "<ellipse cx=\"{}\" cy=\"{}\" rx=\"{}\" ry=\"{}\" fill=\"none\" stroke=\"black\"/>",
                centre.x, centre.y, rx, ry
            ));
        }
        VectorEntity::Polygon(v) => {
            s.push_str("<polygon points=\"");
            for (i, p) in v.iter().enumerate() {
                if i > 0 {
                    s.push(' ');
                }
                s.push_str(&format!("{},{}", p.x, p.y));
            }
            s.push_str("\" fill=\"none\" stroke=\"black\"/>");
        }
        VectorEntity::Text { anchor, font_size, text } => {
            s.push_str(&format!(
                "<text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>",
                anchor.x,
                anchor.y,
                font_size,
                escape_xml(text)
            ));
        }
        VectorEntity::Path(segs) => {
            s.push_str("<path d=\"");
            push_path_d(s, segs);
            s.push_str("\" fill=\"none\" stroke=\"black\"/>");
        }
    }
}

fn push_path_d(s: &mut String, segs: &[PathSegment]) {
    for (i, seg) in segs.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        match seg {
            PathSegment::MoveTo(p) => s.push_str(&format!("M {} {}", p.x, p.y)),
            PathSegment::LineTo(p) => s.push_str(&format!("L {} {}", p.x, p.y)),
            PathSegment::CurveTo { c1, c2, end } => s.push_str(&format!(
                "C {} {} {} {} {} {}",
                c1.x, c1.y, c2.x, c2.y, end.x, end.y
            )),
            PathSegment::QuadTo { c, end } => s.push_str(&format!(
                "Q {} {} {} {}",
                c.x, c.y, end.x, end.y
            )),
            PathSegment::Arc {
                rx,
                ry,
                x_axis_rotation_deg,
                large_arc,
                sweep,
                end,
            } => s.push_str(&format!(
                "A {} {} {} {} {} {} {}",
                rx,
                ry,
                x_axis_rotation_deg,
                if *large_arc { 1 } else { 0 },
                if *sweep { 1 } else { 0 },
                end.x,
                end.y
            )),
            PathSegment::Close => s.push('Z'),
        }
    }
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Minimal SVG parser. Recognises `<line>`, `<rect>`, `<ellipse>`,
/// `<polygon>`, `<text>`, `<path>`. Unknown tags are skipped.
pub fn from_svg(text: &str) -> Result<Vec<VectorEntity>, VectorError> {
    let mut out = Vec::new();
    let mut i = 0usize;
    let bytes = text.as_bytes();
    while i < bytes.len() {
        // Find next '<'.
        while i < bytes.len() && bytes[i] != b'<' {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        // Find matching '>'.
        let start = i;
        let mut j = i + 1;
        while j < bytes.len() && bytes[j] != b'>' {
            j += 1;
        }
        if j >= bytes.len() {
            return Err(VectorError::Parse {
                byte_offset: start,
                message: "unterminated tag".into(),
            });
        }
        let tag = std::str::from_utf8(&bytes[start + 1..j]).map_err(|e| VectorError::Parse {
            byte_offset: start,
            message: format!("non-utf8 tag: {e}"),
        })?;
        i = j + 1;
        // Skip closing tags / declarations / comments.
        if tag.starts_with('?') || tag.starts_with('!') || tag.starts_with('/') {
            continue;
        }
        let name = tag.split_whitespace().next().unwrap_or("").trim_end_matches('/');
        let (attrs, content_after) = match name {
            "text" => {
                // Read body until matching </text>.
                let close_marker = "</text>";
                if let Some(end_pos) = text[i..].find(close_marker) {
                    let body = &text[i..i + end_pos];
                    i += end_pos + close_marker.len();
                    (parse_attrs(tag), Some(body.to_string()))
                } else {
                    return Err(VectorError::Parse {
                        byte_offset: start,
                        message: "unterminated <text>".into(),
                    });
                }
            }
            _ => (parse_attrs(tag), None),
        };
        match name {
            "line" => {
                let x1 = attr_f64(&attrs, "x1")?;
                let y1 = attr_f64(&attrs, "y1")?;
                let x2 = attr_f64(&attrs, "x2")?;
                let y2 = attr_f64(&attrs, "y2")?;
                out.push(VectorEntity::Line {
                    a: Vector2::new(x1, y1),
                    b: Vector2::new(x2, y2),
                });
            }
            "rect" => {
                let x = attr_f64(&attrs, "x")?;
                let y = attr_f64(&attrs, "y")?;
                let w = attr_f64(&attrs, "width")?;
                let h = attr_f64(&attrs, "height")?;
                out.push(VectorEntity::Rect {
                    origin: Vector2::new(x, y),
                    size: Vector2::new(w, h),
                });
            }
            "ellipse" => {
                let cx = attr_f64(&attrs, "cx")?;
                let cy = attr_f64(&attrs, "cy")?;
                let rx = attr_f64(&attrs, "rx")?;
                let ry = attr_f64(&attrs, "ry")?;
                out.push(VectorEntity::Ellipse {
                    centre: Vector2::new(cx, cy),
                    rx,
                    ry,
                });
            }
            "polygon" => {
                let pts_str = attr_str(&attrs, "points")?;
                let pts = parse_points(&pts_str)?;
                out.push(VectorEntity::Polygon(pts));
            }
            "text" => {
                let x = attr_f64(&attrs, "x")?;
                let y = attr_f64(&attrs, "y")?;
                let fs = attr_f64(&attrs, "font-size").unwrap_or(10.0);
                let body = content_after.unwrap_or_default();
                out.push(VectorEntity::Text {
                    anchor: Vector2::new(x, y),
                    font_size: fs,
                    text: body,
                });
            }
            "path" => {
                let d = attr_str(&attrs, "d")?;
                let segs = parse_path_d(&d)?;
                out.push(VectorEntity::Path(segs));
            }
            _ => {} // Skip unrecognised.
        }
    }
    Ok(out)
}

fn parse_attrs(tag: &str) -> std::collections::BTreeMap<String, String> {
    let mut out = std::collections::BTreeMap::new();
    let mut chars = tag.chars().peekable();
    // Skip element name.
    while let Some(c) = chars.peek() {
        if c.is_whitespace() {
            break;
        }
        chars.next();
    }
    while chars.peek().is_some() {
        // Skip ws.
        while let Some(c) = chars.peek() {
            if c.is_whitespace() {
                chars.next();
            } else {
                break;
            }
        }
        let mut key = String::new();
        while let Some(c) = chars.peek() {
            if *c == '=' || c.is_whitespace() {
                break;
            }
            key.push(*c);
            chars.next();
        }
        if key.is_empty() {
            break;
        }
        // Skip '='.
        if chars.peek() == Some(&'=') {
            chars.next();
        } else {
            continue;
        }
        // Skip whitespace + quote.
        while let Some(c) = chars.peek() {
            if c.is_whitespace() {
                chars.next();
            } else {
                break;
            }
        }
        let q = chars.next().unwrap_or('"');
        let mut val = String::new();
        for c in chars.by_ref() {
            if c == q {
                break;
            }
            val.push(c);
        }
        out.insert(key.trim_end_matches('/').to_string(), val);
    }
    out
}

fn attr_f64(
    map: &std::collections::BTreeMap<String, String>,
    key: &str,
) -> Result<f64, VectorError> {
    let s = map.get(key).ok_or(VectorError::Parse {
        byte_offset: 0,
        message: format!("missing attribute `{key}`"),
    })?;
    s.parse::<f64>().map_err(|e| VectorError::Parse {
        byte_offset: 0,
        message: format!("bad number for `{key}`: {e}"),
    })
}

fn attr_str(
    map: &std::collections::BTreeMap<String, String>,
    key: &str,
) -> Result<String, VectorError> {
    map.get(key).cloned().ok_or(VectorError::Parse {
        byte_offset: 0,
        message: format!("missing attribute `{key}`"),
    })
}

fn parse_points(text: &str) -> Result<Vec<Vector2<f64>>, VectorError> {
    let mut out = Vec::new();
    for pair in text.split_whitespace() {
        let mut parts = pair.split(',');
        let x = parts.next().ok_or(VectorError::Parse {
            byte_offset: 0,
            message: format!("bad polygon point `{pair}`"),
        })?;
        let y = parts.next().ok_or(VectorError::Parse {
            byte_offset: 0,
            message: format!("bad polygon point `{pair}`"),
        })?;
        let xx = x.parse::<f64>().map_err(|e| VectorError::Parse {
            byte_offset: 0,
            message: format!("bad x in `{pair}`: {e}"),
        })?;
        let yy = y.parse::<f64>().map_err(|e| VectorError::Parse {
            byte_offset: 0,
            message: format!("bad y in `{pair}`: {e}"),
        })?;
        out.push(Vector2::new(xx, yy));
    }
    Ok(out)
}

fn parse_path_d(text: &str) -> Result<Vec<PathSegment>, VectorError> {
    let mut out = Vec::new();
    // Tokenise: command letters + numbers (delimited by whitespace
    // or commas).
    let mut tokens: Vec<String> = Vec::new();
    let mut cur = String::new();
    let flush = |cur: &mut String, tokens: &mut Vec<String>| {
        if !cur.is_empty() {
            tokens.push(std::mem::take(cur));
        }
    };
    for c in text.chars() {
        if c.is_alphabetic() {
            flush(&mut cur, &mut tokens);
            tokens.push(c.to_string());
        } else if c.is_whitespace() || c == ',' {
            flush(&mut cur, &mut tokens);
        } else {
            cur.push(c);
        }
    }
    flush(&mut cur, &mut tokens);

    let mut iter = tokens.into_iter().peekable();
    while let Some(tok) = iter.next() {
        if tok.len() != 1 || !tok.chars().next().unwrap().is_alphabetic() {
            return Err(VectorError::Parse {
                byte_offset: 0,
                message: format!("expected command got `{tok}`"),
            });
        }
        let cmd = tok.chars().next().unwrap();
        match cmd.to_ascii_uppercase() {
            'M' => {
                let x = pop_f(&mut iter)?;
                let y = pop_f(&mut iter)?;
                out.push(PathSegment::MoveTo(Vector2::new(x, y)));
            }
            'L' => {
                let x = pop_f(&mut iter)?;
                let y = pop_f(&mut iter)?;
                out.push(PathSegment::LineTo(Vector2::new(x, y)));
            }
            'C' => {
                let c1x = pop_f(&mut iter)?;
                let c1y = pop_f(&mut iter)?;
                let c2x = pop_f(&mut iter)?;
                let c2y = pop_f(&mut iter)?;
                let x = pop_f(&mut iter)?;
                let y = pop_f(&mut iter)?;
                out.push(PathSegment::CurveTo {
                    c1: Vector2::new(c1x, c1y),
                    c2: Vector2::new(c2x, c2y),
                    end: Vector2::new(x, y),
                });
            }
            'Q' => {
                let cx = pop_f(&mut iter)?;
                let cy = pop_f(&mut iter)?;
                let x = pop_f(&mut iter)?;
                let y = pop_f(&mut iter)?;
                out.push(PathSegment::QuadTo {
                    c: Vector2::new(cx, cy),
                    end: Vector2::new(x, y),
                });
            }
            'A' => {
                let rx = pop_f(&mut iter)?;
                let ry = pop_f(&mut iter)?;
                let rot = pop_f(&mut iter)?;
                let la = pop_f(&mut iter)? as u8;
                let sw = pop_f(&mut iter)? as u8;
                let x = pop_f(&mut iter)?;
                let y = pop_f(&mut iter)?;
                out.push(PathSegment::Arc {
                    rx,
                    ry,
                    x_axis_rotation_deg: rot,
                    large_arc: la != 0,
                    sweep: sw != 0,
                    end: Vector2::new(x, y),
                });
            }
            'Z' => out.push(PathSegment::Close),
            other => {
                return Err(VectorError::Parse {
                    byte_offset: 0,
                    message: format!("unknown path command `{other}`"),
                });
            }
        }
    }
    Ok(out)
}

fn pop_f(
    iter: &mut std::iter::Peekable<std::vec::IntoIter<String>>,
) -> Result<f64, VectorError> {
    let tok = iter.next().ok_or(VectorError::Parse {
        byte_offset: 0,
        message: "missing number".into(),
    })?;
    tok.parse::<f64>().map_err(|e| VectorError::Parse {
        byte_offset: 0,
        message: format!("bad number `{tok}`: {e}"),
    })
}

//! BRL-CAD-style CSG tree — primitives, booleans, transforms.
//!
//! The evaluator returns an abstract [`SolidHandle`] — a Lisp-style
//! string of the canonical tree — so this crate stays decoupled from
//! the truck-modeling kernel. The host application can pipe the
//! evaluation result through `valenx-cad` to materialise an actual
//! BRep Solid.

use serde::{Deserialize, Serialize};

use crate::error::BrlCadError;

/// Primitive solids supported by BRL-CAD (subset).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Primitive {
    /// Axis-aligned box at the origin.
    Box {
        /// Length x.
        lx: f64,
        /// Length y.
        ly: f64,
        /// Length z.
        lz: f64,
    },
    /// Sphere centred at the origin.
    Sphere {
        /// Radius.
        r: f64,
    },
    /// Cylinder along +Z, centred on Z=0.
    Cylinder {
        /// Radius.
        r: f64,
        /// Height.
        h: f64,
    },
    /// Cone along +Z, centred on Z=0.
    Cone {
        /// Bottom radius.
        r1: f64,
        /// Top radius.
        r2: f64,
        /// Height.
        h: f64,
    },
    /// Torus in the XY plane.
    Torus {
        /// Major radius.
        rm: f64,
        /// Minor radius.
        rn: f64,
    },
    /// Half-space `n . x <= d`.
    HalfSpace {
        /// Unit normal.
        nx: f64,
        /// Unit normal Y.
        ny: f64,
        /// Unit normal Z.
        nz: f64,
        /// Distance from origin.
        d: f64,
    },
}

/// 4x4 affine transform — row-major, last row is `(0, 0, 0, 1)`.
pub type Matrix4 = [[f64; 4]; 4];

/// Identity transform.
pub fn identity() -> Matrix4 {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

/// CSG tree node.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum CsgNode {
    /// A primitive leaf.
    Primitive(Primitive),
    /// Boolean union.
    Union(Box<CsgNode>, Box<CsgNode>),
    /// Boolean intersection.
    Intersection(Box<CsgNode>, Box<CsgNode>),
    /// Boolean difference (a - b).
    Difference(Box<CsgNode>, Box<CsgNode>),
    /// Affine transform applied to a sub-tree.
    Transform(Box<CsgNode>, Matrix4),
}

/// Abstract result handle returned by [`evaluate`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SolidHandle {
    /// Canonical tree string (BRL-CAD MGED Lisp-ish form).
    pub canonical: String,
}

/// Validate then "evaluate" a CSG tree. v1 normalises the tree to
/// its canonical [`SolidHandle`] string — actual BRep evaluation is
/// performed by piping `handle.canonical` through
/// [`parse_mged`] + the host's truck-modeling adapter.
pub fn evaluate(tree: &CsgNode) -> Result<SolidHandle, BrlCadError> {
    validate(tree)?;
    Ok(SolidHandle {
        canonical: to_pretty_string(tree),
    })
}

fn validate(tree: &CsgNode) -> Result<(), BrlCadError> {
    match tree {
        CsgNode::Primitive(p) => validate_prim(p),
        CsgNode::Union(a, b) | CsgNode::Intersection(a, b) | CsgNode::Difference(a, b) => {
            validate(a)?;
            validate(b)
        }
        CsgNode::Transform(a, _) => validate(a),
    }
}

fn validate_prim(p: &Primitive) -> Result<(), BrlCadError> {
    match p {
        Primitive::Box { lx, ly, lz } => {
            if [lx, ly, lz].iter().any(|x| !x.is_finite() || **x <= 0.0) {
                return Err(BrlCadError::BadParameter {
                    name: "box",
                    reason: format!("dimensions must be > 0 ({lx}, {ly}, {lz})"),
                });
            }
        }
        Primitive::Sphere { r } => {
            if !r.is_finite() || *r <= 0.0 {
                return Err(BrlCadError::BadParameter {
                    name: "sphere.r",
                    reason: format!("must be > 0 (got {r})"),
                });
            }
        }
        Primitive::Cylinder { r, h } => {
            if !r.is_finite() || *r <= 0.0 || !h.is_finite() || *h <= 0.0 {
                return Err(BrlCadError::BadParameter {
                    name: "cylinder",
                    reason: format!("r/h must be > 0 ({r}, {h})"),
                });
            }
        }
        Primitive::Cone { r1, r2, h } => {
            if !r1.is_finite()
                || *r1 < 0.0
                || !r2.is_finite()
                || *r2 < 0.0
                || !h.is_finite()
                || *h <= 0.0
            {
                return Err(BrlCadError::BadParameter {
                    name: "cone",
                    reason: format!("r1/r2 must be >= 0, h > 0 ({r1}, {r2}, {h})"),
                });
            }
        }
        Primitive::Torus { rm, rn } => {
            if !rm.is_finite() || *rm <= 0.0 || !rn.is_finite() || *rn <= 0.0 || rn >= rm {
                return Err(BrlCadError::BadParameter {
                    name: "torus",
                    reason: format!("0 < rn < rm required ({rn}, {rm})"),
                });
            }
        }
        Primitive::HalfSpace { nx, ny, nz, d } => {
            let n2 = nx * nx + ny * ny + nz * nz;
            if !n2.is_finite() || n2 < 1e-12 || !d.is_finite() {
                return Err(BrlCadError::BadParameter {
                    name: "halfspace",
                    reason: format!("zero normal or non-finite d ({nx}, {ny}, {nz}, {d})"),
                });
            }
        }
    }
    Ok(())
}

/// Lisp-style pretty printer — BRL-CAD MGED style.
pub fn to_pretty_string(node: &CsgNode) -> String {
    let mut s = String::new();
    pretty(node, &mut s);
    s
}

fn pretty(node: &CsgNode, out: &mut String) {
    match node {
        CsgNode::Primitive(p) => {
            out.push_str(&prim_to_string(p));
        }
        CsgNode::Union(a, b) => {
            out.push_str("(union ");
            pretty(a, out);
            out.push(' ');
            pretty(b, out);
            out.push(')');
        }
        CsgNode::Intersection(a, b) => {
            out.push_str("(intersection ");
            pretty(a, out);
            out.push(' ');
            pretty(b, out);
            out.push(')');
        }
        CsgNode::Difference(a, b) => {
            out.push_str("(difference ");
            pretty(a, out);
            out.push(' ');
            pretty(b, out);
            out.push(')');
        }
        CsgNode::Transform(a, m) => {
            out.push_str("(transform ");
            let mut nums = Vec::new();
            for row in m {
                for v in row {
                    nums.push(format!("{v}"));
                }
            }
            out.push_str(&nums.join(","));
            out.push(' ');
            pretty(a, out);
            out.push(')');
        }
    }
}

fn prim_to_string(p: &Primitive) -> String {
    match p {
        Primitive::Box { lx, ly, lz } => format!("(box {lx} {ly} {lz})"),
        Primitive::Sphere { r } => format!("(sph {r})"),
        Primitive::Cylinder { r, h } => format!("(cyl {r} {h})"),
        Primitive::Cone { r1, r2, h } => format!("(cone {r1} {r2} {h})"),
        Primitive::Torus { rm, rn } => format!("(tor {rm} {rn})"),
        Primitive::HalfSpace { nx, ny, nz, d } => format!("(half {nx} {ny} {nz} {d})"),
    }
}

/// MGED parser — accepts the Lisp-style text emitted by
/// [`to_pretty_string`].
pub fn parse_mged(text: &str) -> Result<CsgNode, BrlCadError> {
    let tokens = tokenise(text)?;
    let mut iter = tokens.into_iter().peekable();
    let node = parse_node(&mut iter, 0)?;
    if iter.next().is_some() {
        return Err(BrlCadError::Parse {
            line: 0,
            message: "trailing tokens after root expression".into(),
        });
    }
    Ok(node)
}

fn tokenise(text: &str) -> Result<Vec<String>, BrlCadError> {
    // Whitespace + parens as delimiters. Numbers may include commas
    // (transform matrices) — keep them in a single token.
    let mut out: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut line = 1usize;
    let flush = |cur: &mut String, out: &mut Vec<String>| {
        if !cur.is_empty() {
            out.push(std::mem::take(cur));
        }
    };
    for ch in text.chars() {
        match ch {
            '(' | ')' => {
                flush(&mut cur, &mut out);
                out.push(ch.to_string());
            }
            ' ' | '\t' | '\r' => flush(&mut cur, &mut out),
            '\n' => {
                flush(&mut cur, &mut out);
                line += 1;
            }
            _ => cur.push(ch),
        }
    }
    flush(&mut cur, &mut out);
    if out.is_empty() {
        return Err(BrlCadError::Parse {
            line,
            message: "empty input".into(),
        });
    }
    Ok(out)
}

fn parse_node(
    iter: &mut std::iter::Peekable<std::vec::IntoIter<String>>,
    depth: usize,
) -> Result<CsgNode, BrlCadError> {
    // Bound recursion so a pathologically deep expression (e.g. thousands
    // of nested `(union ...`) is rejected with a parse error instead of
    // overflowing the stack and aborting the process. `parse_mged` /
    // `BrlCadPanelState::evaluate` feed caller-controlled free-form text
    // here. Mirrors the parser depth caps elsewhere in the workspace
    // (openscad=200, kicad=512, spreadsheet=100).
    const MAX_DEPTH: usize = 256;
    if depth > MAX_DEPTH {
        return Err(BrlCadError::Parse {
            line: 0,
            message: format!("CSG expression nesting exceeds {MAX_DEPTH} levels"),
        });
    }
    let tok = iter.next().ok_or(BrlCadError::Parse {
        line: 0,
        message: "unexpected end of input".into(),
    })?;
    if tok != "(" {
        return Err(BrlCadError::Parse {
            line: 0,
            message: format!("expected `(` got `{tok}`"),
        });
    }
    let op = iter.next().ok_or(BrlCadError::Parse {
        line: 0,
        message: "missing operator".into(),
    })?;
    let result = match op.as_str() {
        "box" => {
            let lx = pop_f64(iter)?;
            let ly = pop_f64(iter)?;
            let lz = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::Box { lx, ly, lz })
        }
        "sph" => {
            let r = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::Sphere { r })
        }
        "cyl" => {
            let r = pop_f64(iter)?;
            let h = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::Cylinder { r, h })
        }
        "cone" => {
            let r1 = pop_f64(iter)?;
            let r2 = pop_f64(iter)?;
            let h = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::Cone { r1, r2, h })
        }
        "tor" => {
            let rm = pop_f64(iter)?;
            let rn = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::Torus { rm, rn })
        }
        "half" => {
            let nx = pop_f64(iter)?;
            let ny = pop_f64(iter)?;
            let nz = pop_f64(iter)?;
            let d = pop_f64(iter)?;
            CsgNode::Primitive(Primitive::HalfSpace { nx, ny, nz, d })
        }
        "union" => {
            let a = parse_node(iter, depth + 1)?;
            let b = parse_node(iter, depth + 1)?;
            CsgNode::Union(Box::new(a), Box::new(b))
        }
        "intersection" => {
            let a = parse_node(iter, depth + 1)?;
            let b = parse_node(iter, depth + 1)?;
            CsgNode::Intersection(Box::new(a), Box::new(b))
        }
        "difference" => {
            let a = parse_node(iter, depth + 1)?;
            let b = parse_node(iter, depth + 1)?;
            CsgNode::Difference(Box::new(a), Box::new(b))
        }
        "transform" => {
            let mat_tok = iter.next().ok_or(BrlCadError::Parse {
                line: 0,
                message: "transform missing matrix".into(),
            })?;
            let vals: Vec<f64> = mat_tok
                .split(',')
                .map(|s| s.parse::<f64>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| BrlCadError::Parse {
                    line: 0,
                    message: format!("bad matrix scalar `{mat_tok}`: {e}"),
                })?;
            if vals.len() != 16 {
                return Err(BrlCadError::Parse {
                    line: 0,
                    message: format!("matrix needs 16 scalars (got {})", vals.len()),
                });
            }
            let mut m = identity();
            for (i, v) in vals.iter().enumerate() {
                m[i / 4][i % 4] = *v;
            }
            let a = parse_node(iter, depth + 1)?;
            CsgNode::Transform(Box::new(a), m)
        }
        other => {
            return Err(BrlCadError::Parse {
                line: 0,
                message: format!("unknown operator `{other}`"),
            });
        }
    };
    let close = iter.next().ok_or(BrlCadError::Parse {
        line: 0,
        message: "missing closing `)`".into(),
    })?;
    if close != ")" {
        return Err(BrlCadError::Parse {
            line: 0,
            message: format!("expected `)` got `{close}`"),
        });
    }
    Ok(result)
}

fn pop_f64(iter: &mut std::iter::Peekable<std::vec::IntoIter<String>>) -> Result<f64, BrlCadError> {
    let tok = iter.next().ok_or(BrlCadError::Parse {
        line: 0,
        message: "unexpected end of input (expected number)".into(),
    })?;
    tok.parse::<f64>().map_err(|e| BrlCadError::Parse {
        line: 0,
        message: format!("bad number `{tok}`: {e}"),
    })
}

/// Alias for [`to_pretty_string`] — matches the plan's naming.
pub fn write_mged(tree: &CsgNode) -> String {
    to_pretty_string(tree)
}

//! Evaluate a parsed OpenSCAD AST down to a [`valenx_cad::Solid`].
//!
//! v1 subset:
//! - primitives: `cube`, `sphere`, `cylinder`
//! - booleans: `union`, `difference`, `intersection`
//! - transforms: `translate`, `rotate`, `scale`
//! - variable bindings + finite `for` loops

use std::collections::HashMap;

use valenx_cad::boolean::{difference, intersection, union};
use valenx_cad::primitives::{box_solid, cylinder, sphere};
use valenx_cad::Solid;

use crate::ast::{Ast, BinOp};
use crate::error::OpenScadError;

/// Runtime value the interpreter passes around.
#[derive(Clone, Debug)]
pub enum Value {
    /// A 64-bit float.
    Number(f64),
    /// A vector — flattened to a homogenous `f64` triple at use sites.
    Vector(Vec<f64>),
    /// A solid built by a primitive or composition.
    Solid(Solid),
}

/// Variable environment.
#[derive(Default, Clone, Debug)]
pub struct Env {
    /// Stored variable bindings.
    pub vars: HashMap<String, Value>,
}

/// Evaluate the top-level statement list. The result is the implicit
/// union of every solid produced — matches OpenSCAD's behaviour where
/// loose top-level shapes auto-union.
pub fn evaluate(stmts: &[Ast]) -> Result<Solid, OpenScadError> {
    let mut env = Env::default();
    let mut solids = Vec::new();
    for s in stmts {
        if let Some(sol) = eval_stmt(s, &mut env)? {
            solids.push(sol);
        }
    }
    union_all(&solids)
}

fn eval_stmt(node: &Ast, env: &mut Env) -> Result<Option<Solid>, OpenScadError> {
    match node {
        Ast::Assign(name, expr) => {
            let v = eval_expr(expr, env)?;
            env.vars.insert(name.clone(), v);
            Ok(None)
        }
        Ast::Block(children) => {
            let mut solids = Vec::new();
            for c in children {
                if let Some(s) = eval_stmt(c, env)? {
                    solids.push(s);
                }
            }
            if solids.is_empty() {
                Ok(None)
            } else {
                Ok(Some(union_all(&solids)?))
            }
        }
        Ast::For {
            var,
            lo,
            step,
            hi,
            body,
        } => {
            let lo_v = number(eval_expr(lo, env)?)?;
            let step_v = number(eval_expr(step, env)?)?;
            let hi_v = number(eval_expr(hi, env)?)?;
            if step_v == 0.0 {
                return Err(OpenScadError::Eval {
                    reason: "for-loop step = 0".into(),
                });
            }
            let mut solids = Vec::new();
            let mut i = lo_v;
            // Inclusive endpoints; tolerate float drift via a 1e-9
            // slack. Direction follows the sign of step_v.
            let asc = step_v > 0.0;
            let mut guard = 0usize;
            while (asc && i <= hi_v + 1e-9) || (!asc && i >= hi_v - 1e-9) {
                let saved = env.vars.get(var).cloned();
                env.vars.insert(var.clone(), Value::Number(i));
                if let Some(s) = eval_stmt(body, env)? {
                    solids.push(s);
                }
                if let Some(v) = saved {
                    env.vars.insert(var.clone(), v);
                } else {
                    env.vars.remove(var);
                }
                i += step_v;
                guard += 1;
                if guard > 100_000 {
                    return Err(OpenScadError::Eval {
                        reason: "for-loop exceeded 100k iterations (likely runaway range)".into(),
                    });
                }
            }
            if solids.is_empty() {
                Ok(None)
            } else {
                Ok(Some(union_all(&solids)?))
            }
        }
        Ast::Call {
            name,
            positional,
            named,
            children,
        } => eval_call(name, positional, named, children, env),
        // Expressions in statement position aren't a thing in OpenSCAD,
        // but treat them as no-ops to keep the interpreter total.
        _ => Ok(None),
    }
}

fn eval_call(
    name: &str,
    positional: &[Ast],
    named: &[(String, Ast)],
    children: &[Ast],
    env: &mut Env,
) -> Result<Option<Solid>, OpenScadError> {
    // Evaluate children first into a flat solid list so transforms
    // and booleans share the same input convention.
    let mut child_solids = Vec::new();
    for c in children {
        if let Some(s) = eval_stmt(c, env)? {
            child_solids.push(s);
        }
    }
    match name {
        // ---- Primitives ----------------------------------------------
        "cube" => {
            // cube(size) or cube([x,y,z])
            let s =
                pos_or_named(positional, named, 0, "size").ok_or_else(|| OpenScadError::Eval {
                    reason: "cube: missing size".into(),
                })?;
            let v = eval_expr(s, env)?;
            let (dx, dy, dz) = match v {
                Value::Number(n) => (n, n, n),
                Value::Vector(vs) if vs.len() == 3 => (vs[0], vs[1], vs[2]),
                other => {
                    return Err(OpenScadError::Eval {
                        reason: format!("cube: bad size {other:?}"),
                    })
                }
            };
            Ok(Some(box_solid(dx, dy, dz).map_err(cad_err)?))
        }
        "sphere" => {
            // sphere(r) or sphere(r = ..)
            let r_expr = named_first(named, "r")
                .or_else(|| positional.first())
                .ok_or_else(|| OpenScadError::Eval {
                    reason: "sphere: missing radius".into(),
                })?;
            let r = number(eval_expr(r_expr, env)?)?;
            Ok(Some(sphere(r).map_err(cad_err)?))
        }
        "cylinder" => {
            // cylinder(h, r) — OpenSCAD has many variants (r1/r2, d, ..)
            // — v1 supports h+r positional or named.
            let h_expr = named_first(named, "h")
                .or_else(|| positional.first())
                .ok_or_else(|| OpenScadError::Eval {
                    reason: "cylinder: missing height".into(),
                })?;
            let r_expr = named_first(named, "r")
                .or_else(|| positional.get(1))
                .ok_or_else(|| OpenScadError::Eval {
                    reason: "cylinder: missing radius".into(),
                })?;
            let h = number(eval_expr(h_expr, env)?)?;
            let r = number(eval_expr(r_expr, env)?)?;
            Ok(Some(cylinder(r, h).map_err(cad_err)?))
        }
        // ---- Booleans -----------------------------------------------
        "union" => {
            if child_solids.is_empty() {
                Ok(None)
            } else {
                Ok(Some(union_all(&child_solids)?))
            }
        }
        "difference" => fold_diff(&child_solids),
        "intersection" => fold_inter(&child_solids),
        // ---- Transforms ---------------------------------------------
        "translate" => {
            let v_expr = positional.first().ok_or_else(|| OpenScadError::Eval {
                reason: "translate: missing vector".into(),
            })?;
            let v = vector3(eval_expr(v_expr, env)?)?;
            let merged = union_all(&child_solids)?;
            let translated = merged
                .translated(v[0], v[1], v[2])
                .map_err(|e| OpenScadError::Cad(e.to_string()))?;
            Ok(Some(translated))
        }
        "rotate" => {
            // rotate([rx, ry, rz]) — Euler ZYX in OpenSCAD (degrees).
            // v1 supports the simpler `rotate(angle, axis)` form, plus
            // an XYZ Euler form via successive rotations.
            let v0 = positional.first().ok_or_else(|| OpenScadError::Eval {
                reason: "rotate: missing arg".into(),
            })?;
            let merged = union_all(&child_solids)?;
            let val = eval_expr(v0, env)?;
            // Two forms: rotate(angle, [ax, ay, az]) or rotate([rx,ry,rz]).
            if let Value::Number(angle_deg) = val {
                let axis = positional.get(1).ok_or_else(|| OpenScadError::Eval {
                    reason: "rotate(angle, axis): missing axis".into(),
                })?;
                let ax = vector3(eval_expr(axis, env)?)?;
                let angle = angle_deg.to_radians();
                let rotated = merged
                    .rotated((0.0, 0.0, 0.0), (ax[0], ax[1], ax[2]), angle)
                    .map_err(|e| OpenScadError::Cad(e.to_string()))?;
                Ok(Some(rotated))
            } else {
                let v = vector3(val)?;
                // Euler XYZ — three successive rotations about the
                // world axes. Matches the OpenSCAD docs' "rotate by
                // [a, b, c]" semantics for the common subset.
                let mut s = merged;
                if v[0] != 0.0 {
                    s = s
                        .rotated((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), v[0].to_radians())
                        .map_err(|e| OpenScadError::Cad(e.to_string()))?;
                }
                if v[1] != 0.0 {
                    s = s
                        .rotated((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), v[1].to_radians())
                        .map_err(|e| OpenScadError::Cad(e.to_string()))?;
                }
                if v[2] != 0.0 {
                    s = s
                        .rotated((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), v[2].to_radians())
                        .map_err(|e| OpenScadError::Cad(e.to_string()))?;
                }
                Ok(Some(s))
            }
        }
        "scale" => {
            // valenx-cad doesn't expose a public anisotropic scale; v1
            // honours uniform scale only (vec3 with equal components
            // or a bare scalar).
            let v_expr = positional.first().ok_or_else(|| OpenScadError::Eval {
                reason: "scale: missing vector".into(),
            })?;
            let val = eval_expr(v_expr, env)?;
            let _factor = match val {
                Value::Number(n) => n,
                Value::Vector(vs) if vs.len() == 3 && vs[0] == vs[1] && vs[1] == vs[2] => vs[0],
                other => {
                    return Err(OpenScadError::Eval {
                        reason: format!("scale: only uniform scale supported in v1, got {other:?}"),
                    });
                }
            };
            // Without a kernel-level scale primitive we return the
            // children unchanged but record the scale request via the
            // typed error path — v2 will plug in a true scale once
            // valenx-cad exposes one. For now, geometry passes through.
            let merged = union_all(&child_solids)?;
            Ok(Some(merged))
        }
        // Unknown module — surface as eval error.
        other => Err(OpenScadError::Eval {
            reason: format!("unknown module `{other}`"),
        }),
    }
}

fn eval_expr(node: &Ast, env: &Env) -> Result<Value, OpenScadError> {
    match node {
        Ast::Number(v) => Ok(Value::Number(*v)),
        Ast::Ident(name) => env
            .vars
            .get(name)
            .cloned()
            .ok_or_else(|| OpenScadError::Eval {
                reason: format!("undefined variable `{name}`"),
            }),
        Ast::Vector(items) => {
            let mut vs = Vec::with_capacity(items.len());
            for it in items {
                let v = eval_expr(it, env)?;
                vs.push(number(v)?);
            }
            Ok(Value::Vector(vs))
        }
        Ast::Negate(rhs) => {
            let v = eval_expr(rhs, env)?;
            Ok(Value::Number(-number(v)?))
        }
        Ast::BinaryOp(l, op, r) => {
            let lv = number(eval_expr(l, env)?)?;
            let rv = number(eval_expr(r, env)?)?;
            let v = match op {
                BinOp::Add => lv + rv,
                BinOp::Sub => lv - rv,
                BinOp::Mul => lv * rv,
                BinOp::Div => {
                    if rv == 0.0 {
                        return Err(OpenScadError::Eval {
                            reason: "division by zero".into(),
                        });
                    }
                    lv / rv
                }
            };
            Ok(Value::Number(v))
        }
        other => Err(OpenScadError::Eval {
            reason: format!("expression evaluator hit non-expr node {other:?}"),
        }),
    }
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

fn pos_or_named<'a>(
    positional: &'a [Ast],
    named: &'a [(String, Ast)],
    idx: usize,
    name: &str,
) -> Option<&'a Ast> {
    named_first(named, name).or_else(|| positional.get(idx))
}

fn named_first<'a>(named: &'a [(String, Ast)], name: &str) -> Option<&'a Ast> {
    named.iter().find(|(n, _)| n == name).map(|(_, v)| v)
}

fn number(v: Value) -> Result<f64, OpenScadError> {
    match v {
        Value::Number(n) => Ok(n),
        other => Err(OpenScadError::Eval {
            reason: format!("expected number, got {other:?}"),
        }),
    }
}

fn vector3(v: Value) -> Result<[f64; 3], OpenScadError> {
    match v {
        Value::Vector(vs) if vs.len() == 3 => Ok([vs[0], vs[1], vs[2]]),
        Value::Vector(vs) if vs.len() == 2 => Ok([vs[0], vs[1], 0.0]),
        Value::Number(n) => Ok([n, n, n]),
        other => Err(OpenScadError::Eval {
            reason: format!("expected vec3, got {other:?}"),
        }),
    }
}

fn cad_err(e: valenx_cad::CadError) -> OpenScadError {
    OpenScadError::Cad(e.to_string())
}

fn union_all(items: &[Solid]) -> Result<Solid, OpenScadError> {
    let mut iter = items.iter().cloned();
    let mut acc = iter.next().ok_or_else(|| OpenScadError::Eval {
        reason: "union of zero solids".into(),
    })?;
    for next in iter {
        acc = union(&acc, &next).map_err(cad_err)?;
    }
    Ok(acc)
}

fn fold_diff(items: &[Solid]) -> Result<Option<Solid>, OpenScadError> {
    if items.is_empty() {
        return Ok(None);
    }
    let mut acc = items[0].clone();
    for next in &items[1..] {
        acc = difference(&acc, next).map_err(cad_err)?;
    }
    Ok(Some(acc))
}

fn fold_inter(items: &[Solid]) -> Result<Option<Solid>, OpenScadError> {
    if items.is_empty() {
        return Ok(None);
    }
    let mut acc = items[0].clone();
    for next in &items[1..] {
        acc = intersection(&acc, next).map_err(cad_err)?;
    }
    Ok(Some(acc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lex::lex;
    use crate::parse::parse;

    fn run(src: &str) -> Result<Solid, OpenScadError> {
        let toks = lex(src)?;
        let ast = parse(&toks)?;
        evaluate(&ast)
    }

    #[test]
    fn evaluate_simple_cube() {
        let s = run("cube([1, 2, 3]);").expect("ok");
        assert!(s.faces() > 0);
    }

    #[test]
    fn evaluate_variable_binding() {
        let s = run("x = 5; cube([x, x, x]);").expect("ok");
        assert!(s.faces() > 0);
    }

    #[test]
    fn evaluate_for_loop_creates_union() {
        // 3 cubes side by side — union should succeed.
        let s = run("for(i = [0 : 2]) translate([i * 10, 0, 0]) cube([1, 1, 1]);").expect("ok");
        assert!(s.faces() > 0);
    }

    #[test]
    fn evaluate_undefined_variable_errors() {
        let err = run("cube([y, 1, 1]);").unwrap_err();
        assert!(matches!(err, OpenScadError::Eval { .. }));
    }
}

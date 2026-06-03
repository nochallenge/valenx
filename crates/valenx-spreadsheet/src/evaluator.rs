//! Recursive formula evaluator.
//!
//! Walks an [`Expr`] AST against a [`Spreadsheet`] context, resolving
//! cell refs to their numeric value (re-evaluating their formulas
//! lazily as it goes) and dispatching function calls to a built-in
//! function table.
//!
//! Booleans are encoded as `1.0` for true and `0.0` for false to keep
//! the value space uniform (spreadsheet land has no separate boolean
//! type). The comparison helpers (`gt`, `lt`, `eq`, `and`, `or`,
//! `not`) all follow that convention.
//!
//! ## Built-in functions
//!
//! Trig (radians): `sin`, `cos`, `tan`, `asin`, `acos`, `atan`,
//! `atan2(y, x)`.
//! Power / log: `sqrt`, `pow(b, e)`, `exp`, `ln`, `log10`.
//! Sign / round: `abs`, `floor`, `ceil`, `round`.
//! Aggregates: `min(a, b)`, `max(a, b)`.
//! Other: `mod(a, b)`.
//! Logic: `if(cond, then, else)`, `gt(a, b)`, `lt(a, b)`,
//! `eq(a, b)`, `and(a, b)`, `or(a, b)`, `not(x)`.
//!
//! ## Constants
//!
//! `pi` (`std::f64::consts::PI`), `e` (`std::f64::consts::E`).

use std::collections::HashSet;

use crate::cell::{Cell, CellRef};
use crate::error::SpreadsheetError;
use crate::formula::{BinOp, Expr, UnOp};
use crate::parser;
use crate::sheet::Spreadsheet;

/// Round-3 fix: maximum recursion depth for `evaluate_expr`. A deeply
/// nested expression like `=1+(1+(1+(1+(...))))` or a wide call tree
/// would otherwise blow the OS stack before the per-cell cycle
/// detector kicked in. 100 is plenty for any human-authored formula
/// (Excel itself caps at 64) and small enough to keep the stack frame
/// well under a megabyte even with debuginfo.
pub const MAX_EVAL_DEPTH: usize = 100;

/// Walk `expr` against `ctx`, returning the numeric result.
///
/// Top-level helper that internally bootstraps an [`Evaluator`] with
/// an empty visited set. Re-entrant calls (cell refs to formulas to
/// cell refs to ...) go through [`Evaluator::evaluate_expr`] with the
/// shared set.
///
/// # Errors
///
/// - [`SpreadsheetError::EvaluationError`] for unknown
///   names/functions, division by zero, etc.
/// - [`SpreadsheetError::CircularReference`] when a cell's formula
///   transitively depends on itself.
/// - [`SpreadsheetError::EvaluationError`] when nesting exceeds
///   [`MAX_EVAL_DEPTH`].
pub fn evaluate(expr: &Expr, ctx: &Spreadsheet) -> Result<f64, SpreadsheetError> {
    let mut visited = HashSet::new();
    Evaluator {
        ctx,
        visited: &mut visited,
        depth: 0,
    }
    .evaluate_expr(expr)
}

/// Per-evaluation context: borrows the [`Spreadsheet`] and tracks the
/// visited cell set for circular-reference detection.
pub struct Evaluator<'a> {
    /// Workbook to resolve [`CellRef`] expressions against.
    pub ctx: &'a Spreadsheet,
    /// Cells we're currently in the middle of evaluating - used to
    /// detect cycles before recursing.
    pub visited: &'a mut HashSet<CellRef>,
    /// Current recursion depth — incremented on every `evaluate_expr`
    /// call so we can short-circuit before the OS stack overflows.
    pub depth: usize,
}

impl<'a> Evaluator<'a> {
    /// Evaluate one expression.
    pub fn evaluate_expr(&mut self, expr: &Expr) -> Result<f64, SpreadsheetError> {
        if self.depth >= MAX_EVAL_DEPTH {
            return Err(SpreadsheetError::EvaluationError(format!(
                "expression nesting exceeds maximum depth {MAX_EVAL_DEPTH} \
                 — refusing to recurse further (stack overflow guard)"
            )));
        }
        self.depth += 1;
        let result = self.evaluate_expr_inner(expr);
        self.depth -= 1;
        result
    }

    fn evaluate_expr_inner(&mut self, expr: &Expr) -> Result<f64, SpreadsheetError> {
        match expr {
            Expr::Number(n) => Ok(*n),
            Expr::Ref(r) => self.resolve_ref(r),
            Expr::Name(name) => apply_constant(name),
            Expr::Call(name, args) => {
                let mut values = Vec::with_capacity(args.len());
                for a in args {
                    values.push(self.evaluate_expr(a)?);
                }
                apply_function(name, &values)
            }
            Expr::Binary(op, lhs, rhs) => {
                let l = self.evaluate_expr(lhs)?;
                let r = self.evaluate_expr(rhs)?;
                apply_binary(*op, l, r)
            }
            Expr::Unary(op, rhs) => {
                let r = self.evaluate_expr(rhs)?;
                Ok(apply_unary(*op, r))
            }
        }
    }

    fn resolve_ref(&mut self, r: &CellRef) -> Result<f64, SpreadsheetError> {
        if !self.visited.insert(r.clone()) {
            return Err(SpreadsheetError::CircularReference(r.to_string()));
        }
        let cell = self.ctx.cell(r).clone();
        let result = match cell {
            Cell::Empty => Ok(0.0),
            Cell::Number(n) => Ok(n),
            Cell::Text(s) => Err(SpreadsheetError::EvaluationError(format!(
                "cell `{r}` contains text `{s}` - expected a number"
            ))),
            Cell::Formula(src) => {
                let parsed = parser::parse(&src).map_err(|e| match e {
                    SpreadsheetError::ParseError {
                        input,
                        position,
                        reason,
                    } => SpreadsheetError::EvaluationError(format!(
                        "cell `{r}` formula `{input}`: parse error at {position}: {reason}"
                    )),
                    other => other,
                })?;
                self.evaluate_expr(&parsed)
            }
        };
        // We're done with this cell; remove it from the visited set so
        // sibling refs to the same cell (legitimate, not cyclic) don't
        // false-positive.
        self.visited.remove(r);
        result
    }
}

fn apply_binary(op: BinOp, l: f64, r: f64) -> Result<f64, SpreadsheetError> {
    Ok(match op {
        BinOp::Add => l + r,
        BinOp::Sub => l - r,
        BinOp::Mul => l * r,
        BinOp::Div => {
            if r == 0.0 {
                return Err(SpreadsheetError::EvaluationError(
                    "division by zero".to_string(),
                ));
            }
            l / r
        }
        BinOp::Pow => l.powf(r),
    })
}

fn apply_unary(op: UnOp, r: f64) -> f64 {
    match op {
        UnOp::Neg => -r,
    }
}

fn apply_constant(name: &str) -> Result<f64, SpreadsheetError> {
    match name.to_ascii_lowercase().as_str() {
        "pi" => Ok(std::f64::consts::PI),
        "e" => Ok(std::f64::consts::E),
        "true" => Ok(1.0),
        "false" => Ok(0.0),
        other => Err(SpreadsheetError::EvaluationError(format!(
            "unknown name `{other}`"
        ))),
    }
}

fn expect_arity(name: &str, args: &[f64], n: usize) -> Result<(), SpreadsheetError> {
    if args.len() != n {
        return Err(SpreadsheetError::EvaluationError(format!(
            "function `{name}` expects {n} argument(s), got {}",
            args.len()
        )));
    }
    Ok(())
}

fn apply_function(name: &str, args: &[f64]) -> Result<f64, SpreadsheetError> {
    let lname = name.to_ascii_lowercase();
    Ok(match lname.as_str() {
        // Trig (radians).
        "sin" => {
            expect_arity(name, args, 1)?;
            args[0].sin()
        }
        "cos" => {
            expect_arity(name, args, 1)?;
            args[0].cos()
        }
        "tan" => {
            expect_arity(name, args, 1)?;
            args[0].tan()
        }
        "asin" => {
            expect_arity(name, args, 1)?;
            args[0].asin()
        }
        "acos" => {
            expect_arity(name, args, 1)?;
            args[0].acos()
        }
        "atan" => {
            expect_arity(name, args, 1)?;
            args[0].atan()
        }
        "atan2" => {
            expect_arity(name, args, 2)?;
            args[0].atan2(args[1])
        }
        // Power / log.
        "sqrt" => {
            expect_arity(name, args, 1)?;
            args[0].sqrt()
        }
        "pow" => {
            expect_arity(name, args, 2)?;
            args[0].powf(args[1])
        }
        "exp" => {
            expect_arity(name, args, 1)?;
            args[0].exp()
        }
        "ln" => {
            expect_arity(name, args, 1)?;
            args[0].ln()
        }
        "log10" => {
            expect_arity(name, args, 1)?;
            args[0].log10()
        }
        // Sign / round.
        "abs" => {
            expect_arity(name, args, 1)?;
            args[0].abs()
        }
        "floor" => {
            expect_arity(name, args, 1)?;
            args[0].floor()
        }
        "ceil" => {
            expect_arity(name, args, 1)?;
            args[0].ceil()
        }
        "round" => {
            expect_arity(name, args, 1)?;
            args[0].round()
        }
        // Aggregates.
        "min" => {
            expect_arity(name, args, 2)?;
            args[0].min(args[1])
        }
        "max" => {
            expect_arity(name, args, 2)?;
            args[0].max(args[1])
        }
        "mod" => {
            expect_arity(name, args, 2)?;
            if args[1] == 0.0 {
                return Err(SpreadsheetError::EvaluationError("mod by zero".into()));
            }
            args[0].rem_euclid(args[1])
        }
        // Logic.
        "if" => {
            expect_arity(name, args, 3)?;
            if args[0] != 0.0 {
                args[1]
            } else {
                args[2]
            }
        }
        "gt" => {
            expect_arity(name, args, 2)?;
            if args[0] > args[1] {
                1.0
            } else {
                0.0
            }
        }
        "lt" => {
            expect_arity(name, args, 2)?;
            if args[0] < args[1] {
                1.0
            } else {
                0.0
            }
        }
        "eq" => {
            expect_arity(name, args, 2)?;
            if (args[0] - args[1]).abs() < 1e-12 {
                1.0
            } else {
                0.0
            }
        }
        "and" => {
            expect_arity(name, args, 2)?;
            if args[0] != 0.0 && args[1] != 0.0 {
                1.0
            } else {
                0.0
            }
        }
        "or" => {
            expect_arity(name, args, 2)?;
            if args[0] != 0.0 || args[1] != 0.0 {
                1.0
            } else {
                0.0
            }
        }
        "not" => {
            expect_arity(name, args, 1)?;
            if args[0] == 0.0 {
                1.0
            } else {
                0.0
            }
        }
        other => {
            return Err(SpreadsheetError::EvaluationError(format!(
                "unknown function `{other}`"
            )))
        }
    })
}

impl Spreadsheet {
    /// Evaluate the cell at `r`, returning the numeric result.
    ///
    /// Convenience wrapper around [`evaluate`] that bootstraps the
    /// visited-set with `r` itself so a single cell that refers back
    /// to its own coordinate (`=A1` in cell A1) is caught as a cycle.
    ///
    /// # Errors
    ///
    /// See [`evaluate`].
    pub fn evaluate_cell(&self, r: &CellRef) -> Result<f64, SpreadsheetError> {
        let mut visited = HashSet::new();
        let cell = self.cell(r).clone();
        match cell {
            Cell::Empty => Ok(0.0),
            Cell::Number(n) => Ok(n),
            Cell::Text(s) => Err(SpreadsheetError::EvaluationError(format!(
                "cell `{r}` contains text `{s}` - expected a number"
            ))),
            Cell::Formula(src) => {
                // Bootstrap with this cell already in the visited set so
                // direct self-reference (`=A1` in cell A1) is rejected.
                visited.insert(r.clone());
                let parsed = parser::parse(&src)?;
                let mut ev = Evaluator {
                    ctx: self,
                    visited: &mut visited,
                    depth: 0,
                };
                ev.evaluate_expr(&parsed)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ss_with(values: &[(&str, Cell)]) -> Spreadsheet {
        let mut ss = Spreadsheet::new();
        ss.add_sheet("S");
        for (a1, cell) in values {
            let r = CellRef::from_a1("S", a1).unwrap();
            ss.set_cell(&r, cell.clone()).unwrap();
        }
        ss
    }

    fn run_expr(s: &str, ss: &Spreadsheet) -> Result<f64, SpreadsheetError> {
        let expr = parser::parse(s)?;
        evaluate(&expr, ss)
    }

    #[test]
    fn const_literal() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("42", &ss).unwrap(), 42.0);
    }

    #[test]
    fn addition_subtraction() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("1+2", &ss).unwrap(), 3.0);
        assert_eq!(run_expr("5-3", &ss).unwrap(), 2.0);
    }

    #[test]
    fn multiplication_division() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("3*4", &ss).unwrap(), 12.0);
        assert_eq!(run_expr("10/4", &ss).unwrap(), 2.5);
    }

    #[test]
    fn precedence() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("2+3*4", &ss).unwrap(), 14.0);
        assert_eq!(run_expr("(2+3)*4", &ss).unwrap(), 20.0);
    }

    #[test]
    fn unary_negation() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("-3+1", &ss).unwrap(), -2.0);
        assert_eq!(run_expr("--3", &ss).unwrap(), 3.0);
    }

    #[test]
    fn power_right_associative() {
        let ss = Spreadsheet::new();
        // 2^3^2 = 2^9 = 512 (right-assoc)
        assert_eq!(run_expr("2^3^2", &ss).unwrap(), 512.0);
    }

    #[test]
    fn division_by_zero_errors() {
        let ss = Spreadsheet::new();
        let e = run_expr("1/0", &ss).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.evaluation_error");
    }

    #[test]
    fn cell_ref_resolves() {
        let ss = ss_with(&[("A1", Cell::Number(5.0))]);
        assert_eq!(run_expr("S.A1 + 2", &ss).unwrap(), 7.0);
    }

    #[test]
    fn chained_cell_refs() {
        // A1 = 5, A2 = A1 * 2, A3 = A2 + 1.
        let ss = ss_with(&[
            ("A1", Cell::Number(5.0)),
            ("A2", Cell::Formula("S.A1 * 2".into())),
            ("A3", Cell::Formula("S.A2 + 1".into())),
        ]);
        let r = CellRef::from_a1("S", "A3").unwrap();
        assert_eq!(ss.evaluate_cell(&r).unwrap(), 11.0);
    }

    #[test]
    fn empty_cell_is_zero() {
        let ss = ss_with(&[]);
        let r = CellRef::from_a1("S", "A1").unwrap();
        assert_eq!(ss.evaluate_cell(&r).unwrap(), 0.0);
    }

    #[test]
    fn text_cell_errors() {
        let ss = ss_with(&[("A1", Cell::Text("hi".into()))]);
        let r = CellRef::from_a1("S", "A1").unwrap();
        let e = ss.evaluate_cell(&r).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.evaluation_error");
    }

    #[test]
    fn circular_self_reference_errors() {
        let ss = ss_with(&[("A1", Cell::Formula("S.A1 + 1".into()))]);
        let r = CellRef::from_a1("S", "A1").unwrap();
        let e = ss.evaluate_cell(&r).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.circular_reference");
    }

    #[test]
    fn circular_indirect_reference_errors() {
        // A1 -> A2 -> A1.
        let ss = ss_with(&[
            ("A1", Cell::Formula("S.A2 + 1".into())),
            ("A2", Cell::Formula("S.A1 + 1".into())),
        ]);
        let r = CellRef::from_a1("S", "A1").unwrap();
        let e = ss.evaluate_cell(&r).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.circular_reference");
    }

    #[test]
    fn diamond_reference_is_not_circular() {
        // A1 = 1, A2 = A1, A3 = A1, A4 = A2 + A3 -> 2.
        // Cell A1 is visited twice but on disjoint stack frames; the
        // visited-set must allow that.
        let ss = ss_with(&[
            ("A1", Cell::Number(1.0)),
            ("A2", Cell::Formula("S.A1".into())),
            ("A3", Cell::Formula("S.A1".into())),
            ("A4", Cell::Formula("S.A2 + S.A3".into())),
        ]);
        let r = CellRef::from_a1("S", "A4").unwrap();
        assert_eq!(ss.evaluate_cell(&r).unwrap(), 2.0);
    }

    #[test]
    fn pi_constant() {
        let ss = Spreadsheet::new();
        let v = run_expr("pi", &ss).unwrap();
        assert!((v - std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn e_constant() {
        let ss = Spreadsheet::new();
        let v = run_expr("e", &ss).unwrap();
        assert!((v - std::f64::consts::E).abs() < 1e-12);
    }

    #[test]
    fn sin_pi_over_2() {
        let ss = Spreadsheet::new();
        let v = run_expr("sin(pi/2)", &ss).unwrap();
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sqrt_and_pow() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("sqrt(16)", &ss).unwrap(), 4.0);
        assert_eq!(run_expr("pow(2, 10)", &ss).unwrap(), 1024.0);
    }

    #[test]
    fn min_max_abs_round() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("min(3, 5)", &ss).unwrap(), 3.0);
        assert_eq!(run_expr("max(3, 5)", &ss).unwrap(), 5.0);
        assert_eq!(run_expr("abs(-7)", &ss).unwrap(), 7.0);
        assert_eq!(run_expr("round(2.6)", &ss).unwrap(), 3.0);
        assert_eq!(run_expr("floor(2.9)", &ss).unwrap(), 2.0);
        assert_eq!(run_expr("ceil(2.1)", &ss).unwrap(), 3.0);
    }

    #[test]
    fn if_function() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("if(1, 10, 20)", &ss).unwrap(), 10.0);
        assert_eq!(run_expr("if(0, 10, 20)", &ss).unwrap(), 20.0);
        // if + comparison
        assert_eq!(run_expr("if(gt(5, 3), 100, 200)", &ss).unwrap(), 100.0);
    }

    #[test]
    fn logic_helpers() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("and(1, 0)", &ss).unwrap(), 0.0);
        assert_eq!(run_expr("and(1, 1)", &ss).unwrap(), 1.0);
        assert_eq!(run_expr("or(0, 1)", &ss).unwrap(), 1.0);
        assert_eq!(run_expr("not(0)", &ss).unwrap(), 1.0);
        assert_eq!(run_expr("not(5)", &ss).unwrap(), 0.0);
        assert_eq!(run_expr("eq(2, 2)", &ss).unwrap(), 1.0);
    }

    #[test]
    fn mod_function() {
        let ss = Spreadsheet::new();
        assert_eq!(run_expr("mod(10, 3)", &ss).unwrap(), 1.0);
    }

    #[test]
    fn unknown_function_errors() {
        let ss = Spreadsheet::new();
        let e = run_expr("nope(1)", &ss).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.evaluation_error");
    }

    #[test]
    fn unknown_name_errors() {
        let ss = Spreadsheet::new();
        let e = run_expr("alpha", &ss).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.evaluation_error");
    }

    #[test]
    fn case_insensitive_names() {
        let ss = Spreadsheet::new();
        let v = run_expr("SIN(PI/2)", &ss).unwrap();
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn wrong_arity_errors() {
        let ss = Spreadsheet::new();
        let e = run_expr("sin(1, 2)", &ss).unwrap_err();
        assert_eq!(e.code(), "spreadsheet.evaluation_error");
    }

    /// Round-3 fix: a deeply-nested expression tree like
    /// `1*(1*(1*(1*(1*(...)))))` produces a Binary expression chain
    /// that the recursive evaluator would otherwise walk all the way
    /// down — blowing the OS stack before any cycle detector kicked
    /// in. The depth cap (MAX_EVAL_DEPTH = 100) bounds recursion
    /// well below the OS stack size.
    #[test]
    fn depth_cap_rejects_pathological_nesting() {
        let ss = Spreadsheet::new();
        // Build `1*(1*(1*(...)))` 200 levels deep — twice the cap.
        let mut expr = String::from("1");
        for _ in 0..200 {
            expr = format!("1*({expr})");
        }
        let err = run_expr(&expr, &ss).expect_err("deep nesting must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("nesting") || msg.contains("depth"),
            "msg: {msg}"
        );
    }

    /// And the inverse: an expression nested within the cap should
    /// still evaluate cleanly.
    #[test]
    fn depth_cap_does_not_reject_normal_nesting() {
        let ss = Spreadsheet::new();
        // 40 levels deep — well under the cap of 100.
        let mut expr = String::from("1");
        for _ in 0..40 {
            expr = format!("1*({expr})");
        }
        assert_eq!(run_expr(&expr, &ss).unwrap(), 1.0);
    }
}

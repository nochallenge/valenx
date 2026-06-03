//! OpenSCAD abstract syntax tree.

use serde::{Deserialize, Serialize};

/// Binary operator (`+ - * /`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum BinOp {
    /// Add.
    Add,
    /// Subtract.
    Sub,
    /// Multiply.
    Mul,
    /// Divide.
    Div,
}

/// AST node.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Ast {
    /// Numeric literal.
    Number(f64),
    /// Identifier reference (variable lookup).
    Ident(String),
    /// `[a, b, c]` vector literal.
    Vector(Vec<Ast>),
    /// `a op b` arithmetic.
    BinaryOp(Box<Ast>, BinOp, Box<Ast>),
    /// `-a` unary negate.
    Negate(Box<Ast>),
    /// Function / module call. Examples:
    ///
    /// - `cube([1,2,3])`
    /// - `translate([x,0,0]) cube()`
    /// - `translate([x,0,0]) { cube(); sphere(); }`
    Call {
        /// Callee name.
        name: String,
        /// Positional args.
        positional: Vec<Ast>,
        /// Named args (`r = 5`).
        named: Vec<(String, Ast)>,
        /// Child nodes (for transforms / boolean blocks).
        children: Vec<Ast>,
    },
    /// `x = expr;` variable binding statement.
    Assign(String, Box<Ast>),
    /// `for(i = [lo : hi]) child` finite range loop. `step` defaults to
    /// 1; OpenSCAD's `[lo : step : hi]` form is unified here.
    For {
        /// Loop variable name.
        var: String,
        /// Start value.
        lo: Box<Ast>,
        /// Step value (default 1.0).
        step: Box<Ast>,
        /// End value (inclusive).
        hi: Box<Ast>,
        /// Body.
        body: Box<Ast>,
    },
    /// Bare `{ ... }` block. Acts like an implicit union — evaluated
    /// as the union of its children.
    Block(Vec<Ast>),
}

//! Formula AST + token stream produced by the [`crate::parser`].
//!
//! Tokens live here so the lexer and the parser share a single
//! definition; the AST [`Expr`] is what the [`crate::evaluator`]
//! traverses to produce a numeric result.
//!
//! Operator precedence (lowest → highest):
//!
//! ```text
//! +/-           (additive)
//! *, /, mod     (multiplicative)
//! ^             (exponent, right-associative)
//! unary -, +    (prefix)
//! function call, parens, atom
//! ```

use crate::cell::CellRef;

/// Lexer token. Emitted by [`crate::parser::tokenize`].
#[derive(Clone, Debug, PartialEq)]
pub enum Token {
    /// Numeric literal — already parsed to `f64` by the lexer.
    Number(f64),
    /// Bare identifier (may be a function name, a named constant, or
    /// the sheet-name portion of a cell reference).
    Ident(String),
    /// `+`
    Plus,
    /// `-`
    Minus,
    /// `*`
    Star,
    /// `/`
    Slash,
    /// `^`
    Caret,
    /// `(`
    LParen,
    /// `)`
    RParen,
    /// `,`
    Comma,
    /// `.`
    Dot,
    /// End-of-stream sentinel.
    Eof,
}

impl Token {
    /// Short label for error messages (`"Number"`, `"+"`, etc.).
    pub fn label(&self) -> &'static str {
        match self {
            Token::Number(_) => "number",
            Token::Ident(_) => "identifier",
            Token::Plus => "`+`",
            Token::Minus => "`-`",
            Token::Star => "`*`",
            Token::Slash => "`/`",
            Token::Caret => "`^`",
            Token::LParen => "`(`",
            Token::RParen => "`)`",
            Token::Comma => "`,`",
            Token::Dot => "`.`",
            Token::Eof => "end of input",
        }
    }
}

/// Binary operators recognised by the parser.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum BinOp {
    /// `+`
    Add,
    /// `-`
    Sub,
    /// `*`
    Mul,
    /// `/`
    Div,
    /// `^`
    Pow,
}

/// Unary operators recognised by the parser.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum UnOp {
    /// `-x`
    Neg,
}

/// Formula AST node — what [`crate::parser::parse`] produces and
/// what [`crate::evaluator::evaluate`] consumes.
#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    /// Numeric literal (`3.14`).
    Number(f64),
    /// Reference to another cell (`Sheet.A1`).
    Ref(CellRef),
    /// Bare named constant (`pi`, `e`) — the evaluator looks it up in
    /// its constants table. Unknown names error.
    Name(String),
    /// Function call (`sin(x)`, `if(cond, then, else)`).
    Call(String, Vec<Expr>),
    /// `lhs op rhs`.
    Binary(BinOp, Box<Expr>, Box<Expr>),
    /// `op rhs`.
    Unary(UnOp, Box<Expr>),
}

impl Expr {
    /// Convenience constructor for binary nodes — avoids the `Box::new`
    /// boilerplate at every call site.
    pub fn binary(op: BinOp, lhs: Expr, rhs: Expr) -> Self {
        Expr::Binary(op, Box::new(lhs), Box::new(rhs))
    }

    /// Convenience constructor for unary nodes.
    pub fn unary(op: UnOp, rhs: Expr) -> Self {
        Expr::Unary(op, Box::new(rhs))
    }
}

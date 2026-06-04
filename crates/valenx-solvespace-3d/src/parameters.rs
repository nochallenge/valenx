//! Named parameters with expressions — the "Change Parameters" layer.
//!
//! The constraint zoo takes literal `target: f64` dimensions. A real
//! parametric modeller (Fusion's *Parameters*, SolidWorks' equations) instead
//! lets you name a dimension and drive it with an **expression** that may
//! reference other parameters: `width = 50`, `height = width * 1.5`,
//! `bolt_circle = pi * d`. Editing one value re-drives every dependent.
//!
//! This module is that layer: a [`ParameterTable`] of `name → expression`, a
//! small arithmetic interpreter (numbers, `+ - * /`, parentheses, unary minus,
//! the functions `sqrt/sin/cos/tan/abs`, and the constants `pi`/`e`), parameter
//! references resolved recursively, and **cycle detection**. Feed the resolved
//! value into a constraint's `target` to drive geometry from a parameter.
//!
//! The interpreter computes over `f64` only — it parses and folds arithmetic,
//! never executes arbitrary code.

use std::collections::BTreeMap;

/// Error from parsing or resolving a parameter expression.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParamError {
    /// The expression could not be parsed.
    Parse(String),
    /// A referenced parameter does not exist.
    Undefined(String),
    /// A parameter (transitively) references itself.
    Cycle(String),
}

impl std::fmt::Display for ParamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamError::Parse(m) => write!(f, "parse error: {m}"),
            ParamError::Undefined(n) => write!(f, "undefined parameter '{n}'"),
            ParamError::Cycle(n) => write!(f, "cyclic parameter reference involving '{n}'"),
        }
    }
}

impl std::error::Error for ParamError {}

/// A table of named parameters, each defined by an expression string.
#[derive(Clone, Debug, Default)]
pub struct ParameterTable {
    params: BTreeMap<String, String>,
}

impl ParameterTable {
    /// Empty table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Define or replace parameter `name` with `expr` (an expression that may
    /// reference other parameters).
    pub fn set(&mut self, name: &str, expr: &str) {
        self.params.insert(name.to_string(), expr.to_string());
    }

    /// The expression string defining `name`, if any.
    pub fn expr(&self, name: &str) -> Option<&str> {
        self.params.get(name).map(|s| s.as_str())
    }

    /// Parameter names, sorted.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.params.keys().map(|s| s.as_str())
    }

    /// Resolve the numeric value of parameter `name`.
    pub fn value(&self, name: &str) -> Result<f64, ParamError> {
        let src = self
            .params
            .get(name)
            .ok_or_else(|| ParamError::Undefined(name.to_string()))?;
        let ast = parse(src)?;
        let mut visiting = vec![name.to_string()];
        compute_ast(&ast, self, &mut visiting)
    }

    /// Compute an arbitrary expression in this table's context (its parameter
    /// references resolve against this table).
    pub fn compute(&self, expr: &str) -> Result<f64, ParamError> {
        let ast = parse(expr)?;
        let mut visiting = Vec::new();
        compute_ast(&ast, self, &mut visiting)
    }
}

// --- expression engine: tokenize → parse to AST → fold to a number ----------

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Num(f64),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
}

fn tokenize(s: &str) -> Result<Vec<Tok>, ParamError> {
    let chars: Vec<char> = s.chars().collect();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        match c {
            '+' => toks.push(Tok::Plus),
            '-' => toks.push(Tok::Minus),
            '*' => toks.push(Tok::Star),
            '/' => toks.push(Tok::Slash),
            '(' => toks.push(Tok::LParen),
            ')' => toks.push(Tok::RParen),
            _ if c.is_ascii_digit() || c == '.' => {
                let start = i;
                while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                    i += 1;
                }
                let text: String = chars[start..i].iter().collect();
                let n = text
                    .parse::<f64>()
                    .map_err(|_| ParamError::Parse(format!("bad number '{text}'")))?;
                toks.push(Tok::Num(n));
                continue;
            }
            _ if c.is_ascii_alphabetic() || c == '_' => {
                let start = i;
                while i < chars.len() && (chars[i].is_ascii_alphanumeric() || chars[i] == '_') {
                    i += 1;
                }
                toks.push(Tok::Ident(chars[start..i].iter().collect()));
                continue;
            }
            _ => return Err(ParamError::Parse(format!("unexpected character '{c}'"))),
        }
        i += 1;
    }
    Ok(toks)
}

#[derive(Debug, Clone)]
enum Ast {
    Num(f64),
    Ident(String),
    Neg(Box<Ast>),
    Add(Box<Ast>, Box<Ast>),
    Sub(Box<Ast>, Box<Ast>),
    Mul(Box<Ast>, Box<Ast>),
    Div(Box<Ast>, Box<Ast>),
    Call(String, Box<Ast>),
}

struct Parser {
    toks: Vec<Tok>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos)
    }
    fn advance(&mut self) -> Option<Tok> {
        let t = self.toks.get(self.pos).cloned();
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn parse_expr(&mut self) -> Result<Ast, ParamError> {
        let mut lhs = self.parse_term()?;
        while let Some(t) = self.peek() {
            match t {
                Tok::Plus => {
                    self.pos += 1;
                    lhs = Ast::Add(Box::new(lhs), Box::new(self.parse_term()?));
                }
                Tok::Minus => {
                    self.pos += 1;
                    lhs = Ast::Sub(Box::new(lhs), Box::new(self.parse_term()?));
                }
                _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_term(&mut self) -> Result<Ast, ParamError> {
        let mut lhs = self.parse_factor()?;
        while let Some(t) = self.peek() {
            match t {
                Tok::Star => {
                    self.pos += 1;
                    lhs = Ast::Mul(Box::new(lhs), Box::new(self.parse_factor()?));
                }
                Tok::Slash => {
                    self.pos += 1;
                    lhs = Ast::Div(Box::new(lhs), Box::new(self.parse_factor()?));
                }
                _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_factor(&mut self) -> Result<Ast, ParamError> {
        match self.advance() {
            Some(Tok::Num(n)) => Ok(Ast::Num(n)),
            Some(Tok::Minus) => Ok(Ast::Neg(Box::new(self.parse_factor()?))),
            Some(Tok::Plus) => self.parse_factor(),
            Some(Tok::LParen) => {
                let e = self.parse_expr()?;
                match self.advance() {
                    Some(Tok::RParen) => Ok(e),
                    _ => Err(ParamError::Parse("expected ')'".into())),
                }
            }
            Some(Tok::Ident(id)) => {
                if matches!(self.peek(), Some(Tok::LParen)) {
                    self.pos += 1; // consume '('
                    let arg = self.parse_expr()?;
                    match self.advance() {
                        Some(Tok::RParen) => Ok(Ast::Call(id, Box::new(arg))),
                        _ => Err(ParamError::Parse("expected ')' after argument".into())),
                    }
                } else {
                    Ok(Ast::Ident(id))
                }
            }
            other => Err(ParamError::Parse(format!("unexpected token: {other:?}"))),
        }
    }
}

fn parse(expr: &str) -> Result<Ast, ParamError> {
    let mut p = Parser { toks: tokenize(expr)?, pos: 0 };
    let ast = p.parse_expr()?;
    if p.pos != p.toks.len() {
        return Err(ParamError::Parse("trailing tokens after expression".into()));
    }
    Ok(ast)
}

fn apply_fn(name: &str, x: f64) -> Result<f64, ParamError> {
    Ok(match name {
        "sqrt" => x.sqrt(),
        "sin" => x.sin(),
        "cos" => x.cos(),
        "tan" => x.tan(),
        "abs" => x.abs(),
        _ => return Err(ParamError::Parse(format!("unknown function '{name}'"))),
    })
}

fn compute_ast(ast: &Ast, table: &ParameterTable, visiting: &mut Vec<String>) -> Result<f64, ParamError> {
    match ast {
        Ast::Num(n) => Ok(*n),
        Ast::Neg(a) => Ok(-compute_ast(a, table, visiting)?),
        Ast::Add(a, b) => Ok(compute_ast(a, table, visiting)? + compute_ast(b, table, visiting)?),
        Ast::Sub(a, b) => Ok(compute_ast(a, table, visiting)? - compute_ast(b, table, visiting)?),
        Ast::Mul(a, b) => Ok(compute_ast(a, table, visiting)? * compute_ast(b, table, visiting)?),
        Ast::Div(a, b) => Ok(compute_ast(a, table, visiting)? / compute_ast(b, table, visiting)?),
        Ast::Call(name, arg) => {
            let x = compute_ast(arg, table, visiting)?;
            apply_fn(name, x)
        }
        Ast::Ident(id) => {
            match id.as_str() {
                "pi" => return Ok(std::f64::consts::PI),
                "e" => return Ok(std::f64::consts::E),
                _ => {}
            }
            if visiting.iter().any(|n| n == id) {
                return Err(ParamError::Cycle(id.clone()));
            }
            let src = table
                .params
                .get(id)
                .ok_or_else(|| ParamError::Undefined(id.clone()))?;
            let sub = parse(src)?;
            visiting.push(id.clone());
            let v = compute_ast(&sub, table, visiting)?;
            visiting.pop();
            Ok(v)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn arithmetic_and_precedence() {
        let t = ParameterTable::new();
        assert!(approx(t.compute("5").unwrap(), 5.0));
        assert!(approx(t.compute("2 + 3 * 4").unwrap(), 14.0));
        assert!(approx(t.compute("(2 + 3) * 4").unwrap(), 20.0));
        assert!(approx(t.compute("-3 + 1").unwrap(), -2.0));
        assert!(approx(t.compute("10 / 4").unwrap(), 2.5));
    }

    #[test]
    fn functions_and_constants() {
        let t = ParameterTable::new();
        assert!(approx(t.compute("sqrt(16)").unwrap(), 4.0));
        assert!(approx(t.compute("cos(0)").unwrap(), 1.0));
        assert!(approx(t.compute("abs(0 - 7)").unwrap(), 7.0));
        assert!(approx(t.compute("pi").unwrap(), std::f64::consts::PI));
    }

    #[test]
    fn parameter_references_resolve_recursively() {
        let mut t = ParameterTable::new();
        t.set("width", "50");
        t.set("height", "width * 1.5");
        t.set("area", "width * height");
        assert!(approx(t.value("width").unwrap(), 50.0));
        assert!(approx(t.value("height").unwrap(), 75.0));
        assert!(approx(t.value("area").unwrap(), 3750.0));
        assert!(approx(t.compute("height / 2 + 1").unwrap(), 38.5));
    }

    #[test]
    fn undefined_and_cyclic_references_error() {
        let mut t = ParameterTable::new();
        assert_eq!(t.compute("nope"), Err(ParamError::Undefined("nope".into())));
        t.set("a", "b + 1");
        t.set("b", "a + 1");
        assert!(matches!(t.value("a"), Err(ParamError::Cycle(_))));
    }
}

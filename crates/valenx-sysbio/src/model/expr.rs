//! Tiny arithmetic / boolean expression AST — feature 31.
//!
//! [`Expr`] is a small typed AST over species amounts, global parameter
//! values, simulation time and the usual arithmetic / comparison /
//! boolean operators. It is the lingua franca of the three pieces of
//! commercial-depth machinery added in this pass:
//!
//! - SBML L3-style **event triggers** — boolean expressions watched
//!   each integrator step (`Expr` produces a number; the trigger is
//!   "fired" when that number transitions from `<= 0` to `> 0`).
//! - SBML L3-style **event assignments** — numeric expressions assigned
//!   to species / parameters when an event fires.
//! - SBML L3-style **assignment** and **rate** rules — numeric
//!   expressions enforced every integrator output (assignment) or
//!   folded into the ODE right-hand side (rate).
//!
//! Variables are addressed by *index* (the same indexing scheme the
//! rate laws use), so the hot evaluation loop is allocation- and
//! lookup-free. A symbolic [`Expr::Var`] / [`Expr::Param`] is the
//! source-level form; the model builder is responsible for
//! re-indexing if it ever re-orders species or parameters.
//!
//! ## Scope
//!
//! The AST covers what real SBML rules need in practice: `+`, `-`,
//! `*`, `/`, `^`, unary `-`, comparison (`<`, `<=`, `>`, `>=`, `=`,
//! `!=`), boolean `and` / `or` / `not`, and the `exp`, `ln`, `log10`,
//! `sqrt`, `abs`, `floor`, `ceil`, `min`, `max`, `pow` MathML
//! built-ins. A general MathML evaluator (csymbol delays, piecewise
//! lambdas, user-defined functions) is out of scope — the same v1 line
//! the existing `<kineticLaw>` annotation draws.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

/// Symbolic arithmetic / boolean expression computable against a
/// state vector, a parameter vector and the simulation time.
///
/// Numeric results are `f64`. Boolean results are encoded as `1.0`
/// (true) / `0.0` (false) - the unified representation lets a single
/// [`value`](Expr::value) signature serve both trigger and assignment
/// expressions, and lets a trigger be computed as a *signed*
/// crossing function (positive => true, non-positive => false).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Expr {
    /// A numeric literal.
    Const(f64),
    /// The simulation time `t`.
    Time,
    /// A species amount by index into the species vector.
    Var(usize),
    /// A global parameter value by index into the parameter vector.
    Param(usize),
    /// Negation `-a`.
    Neg(Box<Expr>),
    /// Addition `a + b`.
    Add(Box<Expr>, Box<Expr>),
    /// Subtraction `a - b`.
    Sub(Box<Expr>, Box<Expr>),
    /// Multiplication `a * b`.
    Mul(Box<Expr>, Box<Expr>),
    /// Division `a / b`. A divide-by-zero yields `0.0` rather than a
    /// `NaN`, matching the defensive style of the existing rate-law
    /// evaluator.
    Div(Box<Expr>, Box<Expr>),
    /// Power `a^b`.
    Pow(Box<Expr>, Box<Expr>),
    /// `a < b` (encoded `1.0` / `0.0`).
    Lt(Box<Expr>, Box<Expr>),
    /// `a <= b`.
    Le(Box<Expr>, Box<Expr>),
    /// `a > b`.
    Gt(Box<Expr>, Box<Expr>),
    /// `a >= b`.
    Ge(Box<Expr>, Box<Expr>),
    /// `a == b` (numeric equality up to a tiny tolerance).
    Eq(Box<Expr>, Box<Expr>),
    /// `a != b`.
    Ne(Box<Expr>, Box<Expr>),
    /// Boolean `and`.
    And(Box<Expr>, Box<Expr>),
    /// Boolean `or`.
    Or(Box<Expr>, Box<Expr>),
    /// Boolean `not`.
    Not(Box<Expr>),
    /// MathML `min`.
    Min(Box<Expr>, Box<Expr>),
    /// MathML `max`.
    Max(Box<Expr>, Box<Expr>),
    /// `exp(a)`.
    Exp(Box<Expr>),
    /// `ln(a)`.
    Ln(Box<Expr>),
    /// `log10(a)`.
    Log10(Box<Expr>),
    /// `sqrt(a)`.
    Sqrt(Box<Expr>),
    /// `|a|`.
    Abs(Box<Expr>),
    /// `floor(a)`.
    Floor(Box<Expr>),
    /// `ceil(a)`.
    Ceil(Box<Expr>),
}

// The `add` / `sub` / `mul` / `div` / `neg` constructors below have the
// same names as `std::ops` trait methods but are AST builders, not
// arithmetic. Implementing the traits would force `Expr + Expr` to
// allocate an unrelated boxed node and would surprise callers - the
// explicit factory methods are exactly the right shape for the trigger
// / rule expression builders. Suppress the lint here.
#[allow(clippy::should_implement_trait)]
impl Expr {
    /// Constant convenience.
    pub fn k(v: f64) -> Expr {
        Expr::Const(v)
    }
    /// Species reference convenience.
    pub fn var(i: usize) -> Expr {
        Expr::Var(i)
    }
    /// Parameter reference convenience.
    pub fn param(i: usize) -> Expr {
        Expr::Param(i)
    }
    /// Build `Add` with boxing.
    pub fn add(a: Expr, b: Expr) -> Expr {
        Expr::Add(Box::new(a), Box::new(b))
    }
    /// Build `Sub` with boxing.
    pub fn sub(a: Expr, b: Expr) -> Expr {
        Expr::Sub(Box::new(a), Box::new(b))
    }
    /// Build `Mul` with boxing.
    pub fn mul(a: Expr, b: Expr) -> Expr {
        Expr::Mul(Box::new(a), Box::new(b))
    }
    /// Build `Div` with boxing.
    pub fn div(a: Expr, b: Expr) -> Expr {
        Expr::Div(Box::new(a), Box::new(b))
    }
    /// Build `Pow` with boxing.
    pub fn pow(a: Expr, b: Expr) -> Expr {
        Expr::Pow(Box::new(a), Box::new(b))
    }
    /// Build `Neg` with boxing.
    pub fn neg(a: Expr) -> Expr {
        Expr::Neg(Box::new(a))
    }
    /// Build a comparison `a < b`.
    pub fn lt(a: Expr, b: Expr) -> Expr {
        Expr::Lt(Box::new(a), Box::new(b))
    }
    /// Build `a > b`.
    pub fn gt(a: Expr, b: Expr) -> Expr {
        Expr::Gt(Box::new(a), Box::new(b))
    }
    /// Build `a <= b`.
    pub fn le(a: Expr, b: Expr) -> Expr {
        Expr::Le(Box::new(a), Box::new(b))
    }
    /// Build `a >= b`.
    pub fn ge(a: Expr, b: Expr) -> Expr {
        Expr::Ge(Box::new(a), Box::new(b))
    }
    /// Build `a && b`.
    pub fn and(a: Expr, b: Expr) -> Expr {
        Expr::And(Box::new(a), Box::new(b))
    }
    /// Build `a || b`.
    pub fn or(a: Expr, b: Expr) -> Expr {
        Expr::Or(Box::new(a), Box::new(b))
    }
    /// Build `!a`.
    pub fn not_(a: Expr) -> Expr {
        Expr::Not(Box::new(a))
    }
    /// Build `min(a, b)`.
    pub fn min_(a: Expr, b: Expr) -> Expr {
        Expr::Min(Box::new(a), Box::new(b))
    }
    /// Build `max(a, b)`.
    pub fn max_(a: Expr, b: Expr) -> Expr {
        Expr::Max(Box::new(a), Box::new(b))
    }

    /// Compute the expression value at species amounts `y`, parameter
    /// values `p` and simulation time `t`.
    ///
    /// Out-of-range indices contribute a `0.0` (defensive - a validated
    /// model never produces them, but a hand-built tester might). A
    /// boolean operand is converted via `> 0.5` so the truth table
    /// matches the encoding the comparison ops emit.
    pub fn value(&self, y: &[f64], p: &[f64], t: f64) -> f64 {
        let truthy = |x: f64| -> bool { x > 0.5 };
        let pack = |b: bool| -> f64 {
            if b {
                1.0
            } else {
                0.0
            }
        };
        match self {
            Expr::Const(v) => *v,
            Expr::Time => t,
            Expr::Var(i) => y.get(*i).copied().unwrap_or(0.0),
            Expr::Param(i) => p.get(*i).copied().unwrap_or(0.0),
            Expr::Neg(a) => -a.value(y, p, t),
            Expr::Add(a, b1) => a.value(y, p, t) + b1.value(y, p, t),
            Expr::Sub(a, b1) => a.value(y, p, t) - b1.value(y, p, t),
            Expr::Mul(a, b1) => a.value(y, p, t) * b1.value(y, p, t),
            Expr::Div(a, b1) => {
                let d = b1.value(y, p, t);
                if d.abs() < 1e-300 {
                    0.0
                } else {
                    a.value(y, p, t) / d
                }
            }
            Expr::Pow(a, b1) => a.value(y, p, t).powf(b1.value(y, p, t)),
            Expr::Lt(a, b1) => pack(a.value(y, p, t) < b1.value(y, p, t)),
            Expr::Le(a, b1) => pack(a.value(y, p, t) <= b1.value(y, p, t)),
            Expr::Gt(a, b1) => pack(a.value(y, p, t) > b1.value(y, p, t)),
            Expr::Ge(a, b1) => pack(a.value(y, p, t) >= b1.value(y, p, t)),
            Expr::Eq(a, b1) => pack((a.value(y, p, t) - b1.value(y, p, t)).abs() < 1e-12),
            Expr::Ne(a, b1) => pack((a.value(y, p, t) - b1.value(y, p, t)).abs() >= 1e-12),
            Expr::And(a, b1) => pack(truthy(a.value(y, p, t)) && truthy(b1.value(y, p, t))),
            Expr::Or(a, b1) => pack(truthy(a.value(y, p, t)) || truthy(b1.value(y, p, t))),
            Expr::Not(a) => pack(!truthy(a.value(y, p, t))),
            Expr::Min(a, b1) => a.value(y, p, t).min(b1.value(y, p, t)),
            Expr::Max(a, b1) => a.value(y, p, t).max(b1.value(y, p, t)),
            Expr::Exp(a) => a.value(y, p, t).exp(),
            Expr::Ln(a) => {
                let v = a.value(y, p, t);
                if v <= 0.0 {
                    f64::NEG_INFINITY
                } else {
                    v.ln()
                }
            }
            Expr::Log10(a) => {
                let v = a.value(y, p, t);
                if v <= 0.0 {
                    f64::NEG_INFINITY
                } else {
                    v.log10()
                }
            }
            Expr::Sqrt(a) => a.value(y, p, t).max(0.0).sqrt(),
            Expr::Abs(a) => a.value(y, p, t).abs(),
            Expr::Floor(a) => a.value(y, p, t).floor(),
            Expr::Ceil(a) => a.value(y, p, t).ceil(),
        }
    }

    /// Whether this expression reads species `i`.
    pub fn reads_var(&self, i: usize) -> bool {
        self.var_deps().contains(&i)
    }

    /// Whether this expression reads parameter `i`.
    pub fn reads_param(&self, i: usize) -> bool {
        self.param_deps().contains(&i)
    }

    /// The set of species indices read anywhere in the expression.
    pub fn var_deps(&self) -> BTreeSet<usize> {
        let mut out = BTreeSet::new();
        self.collect_var_deps(&mut out);
        out
    }

    /// The set of parameter indices read anywhere in the expression.
    pub fn param_deps(&self) -> BTreeSet<usize> {
        let mut out = BTreeSet::new();
        self.collect_param_deps(&mut out);
        out
    }

    fn collect_var_deps(&self, into: &mut BTreeSet<usize>) {
        match self {
            Expr::Var(i) => {
                into.insert(*i);
            }
            Expr::Const(_) | Expr::Time | Expr::Param(_) => {}
            Expr::Neg(a)
            | Expr::Not(a)
            | Expr::Exp(a)
            | Expr::Ln(a)
            | Expr::Log10(a)
            | Expr::Sqrt(a)
            | Expr::Abs(a)
            | Expr::Floor(a)
            | Expr::Ceil(a) => a.collect_var_deps(into),
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Pow(a, b)
            | Expr::Lt(a, b)
            | Expr::Le(a, b)
            | Expr::Gt(a, b)
            | Expr::Ge(a, b)
            | Expr::Eq(a, b)
            | Expr::Ne(a, b)
            | Expr::And(a, b)
            | Expr::Or(a, b)
            | Expr::Min(a, b)
            | Expr::Max(a, b) => {
                a.collect_var_deps(into);
                b.collect_var_deps(into);
            }
        }
    }

    fn collect_param_deps(&self, into: &mut BTreeSet<usize>) {
        match self {
            Expr::Param(i) => {
                into.insert(*i);
            }
            Expr::Const(_) | Expr::Time | Expr::Var(_) => {}
            Expr::Neg(a)
            | Expr::Not(a)
            | Expr::Exp(a)
            | Expr::Ln(a)
            | Expr::Log10(a)
            | Expr::Sqrt(a)
            | Expr::Abs(a)
            | Expr::Floor(a)
            | Expr::Ceil(a) => a.collect_param_deps(into),
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Pow(a, b)
            | Expr::Lt(a, b)
            | Expr::Le(a, b)
            | Expr::Gt(a, b)
            | Expr::Ge(a, b)
            | Expr::Eq(a, b)
            | Expr::Ne(a, b)
            | Expr::And(a, b)
            | Expr::Or(a, b)
            | Expr::Min(a, b)
            | Expr::Max(a, b) => {
                a.collect_param_deps(into);
                b.collect_param_deps(into);
            }
        }
    }

    /// Serialise the expression to a compact ASCII form. The SBML
    /// writer embeds the result inside a `sbml:expr` attribute; the
    /// reader round-trips it via [`Expr::parse`]. This is the same
    /// "annotation as ground truth" pattern the rate-law writer uses.
    ///
    /// Variables become `s<i>`, parameters become `p<i>`, time is
    /// `t`. Operators use parenthesised infix notation so the parser
    /// can stay flat and predictable.
    pub fn to_string_compact(&self) -> String {
        match self {
            Expr::Const(v) => format!("{v}"),
            Expr::Time => "t".to_string(),
            Expr::Var(i) => format!("s{i}"),
            Expr::Param(i) => format!("p{i}"),
            Expr::Neg(a) => format!("(-{})", a.to_string_compact()),
            Expr::Add(a, b) => format!("({}+{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Sub(a, b) => format!("({}-{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Mul(a, b) => format!("({}*{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Div(a, b) => format!("({}/{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Pow(a, b) => format!("({}^{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Lt(a, b) => format!("({}<{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Le(a, b) => format!("({}<={})", a.to_string_compact(), b.to_string_compact()),
            Expr::Gt(a, b) => format!("({}>{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Ge(a, b) => format!("({}>={})", a.to_string_compact(), b.to_string_compact()),
            Expr::Eq(a, b) => format!("({}=={})", a.to_string_compact(), b.to_string_compact()),
            Expr::Ne(a, b) => format!("({}!={})", a.to_string_compact(), b.to_string_compact()),
            Expr::And(a, b) => format!("({}&{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Or(a, b) => format!("({}|{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Not(a) => format!("(!{})", a.to_string_compact()),
            Expr::Min(a, b) => format!("min({},{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Max(a, b) => format!("max({},{})", a.to_string_compact(), b.to_string_compact()),
            Expr::Exp(a) => format!("exp({})", a.to_string_compact()),
            Expr::Ln(a) => format!("ln({})", a.to_string_compact()),
            Expr::Log10(a) => format!("log10({})", a.to_string_compact()),
            Expr::Sqrt(a) => format!("sqrt({})", a.to_string_compact()),
            Expr::Abs(a) => format!("abs({})", a.to_string_compact()),
            Expr::Floor(a) => format!("floor({})", a.to_string_compact()),
            Expr::Ceil(a) => format!("ceil({})", a.to_string_compact()),
        }
    }

    /// Parse the compact ASCII form back into an `Expr`.
    ///
    /// The grammar accepts the exact output of
    /// [`Expr::to_string_compact`] plus a handful of writability niceties:
    /// optional whitespace around tokens, and the symbolic forms `&&`
    /// / `||` for `&` / `|`. Returns `None` on any structural error.
    pub fn parse(s: &str) -> Option<Expr> {
        let mut p = Parser {
            chars: s.chars().filter(|c| !c.is_whitespace()).collect(),
            pos: 0,
        };
        let e = p.parse_or()?;
        if p.pos < p.chars.len() {
            return None;
        }
        Some(e)
    }
}

// --- compact-ASCII parser --------------------------------------------

struct Parser {
    chars: Vec<char>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }
    fn eat(&mut self, c: char) -> bool {
        if self.peek() == Some(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }
    fn eat_str(&mut self, s: &str) -> bool {
        let chars: Vec<char> = s.chars().collect();
        if self.pos + chars.len() > self.chars.len() {
            return false;
        }
        for (i, &c) in chars.iter().enumerate() {
            if self.chars[self.pos + i] != c {
                return false;
            }
        }
        self.pos += chars.len();
        true
    }

    fn parse_or(&mut self) -> Option<Expr> {
        let mut left = self.parse_and()?;
        loop {
            if self.eat_str("||") || self.eat('|') {
                let right = self.parse_and()?;
                left = Expr::or(left, right);
            } else {
                break;
            }
        }
        Some(left)
    }
    fn parse_and(&mut self) -> Option<Expr> {
        let mut left = self.parse_cmp()?;
        loop {
            if self.eat_str("&&") || self.eat('&') {
                let right = self.parse_cmp()?;
                left = Expr::and(left, right);
            } else {
                break;
            }
        }
        Some(left)
    }
    fn parse_cmp(&mut self) -> Option<Expr> {
        let left = self.parse_add()?;
        if self.eat_str("<=") {
            let right = self.parse_add()?;
            return Some(Expr::le(left, right));
        }
        if self.eat_str(">=") {
            let right = self.parse_add()?;
            return Some(Expr::ge(left, right));
        }
        if self.eat_str("==") {
            let right = self.parse_add()?;
            return Some(Expr::Eq(Box::new(left), Box::new(right)));
        }
        if self.eat_str("!=") {
            let right = self.parse_add()?;
            return Some(Expr::Ne(Box::new(left), Box::new(right)));
        }
        if self.eat('<') {
            let right = self.parse_add()?;
            return Some(Expr::lt(left, right));
        }
        if self.eat('>') {
            let right = self.parse_add()?;
            return Some(Expr::gt(left, right));
        }
        Some(left)
    }
    fn parse_add(&mut self) -> Option<Expr> {
        let mut left = self.parse_mul()?;
        loop {
            if self.eat('+') {
                let right = self.parse_mul()?;
                left = Expr::add(left, right);
            } else if self.eat('-') {
                let right = self.parse_mul()?;
                left = Expr::sub(left, right);
            } else {
                break;
            }
        }
        Some(left)
    }
    fn parse_mul(&mut self) -> Option<Expr> {
        let mut left = self.parse_pow()?;
        loop {
            if self.eat('*') {
                let right = self.parse_pow()?;
                left = Expr::mul(left, right);
            } else if self.eat('/') {
                let right = self.parse_pow()?;
                left = Expr::div(left, right);
            } else {
                break;
            }
        }
        Some(left)
    }
    fn parse_pow(&mut self) -> Option<Expr> {
        let base = self.parse_unary()?;
        if self.eat('^') {
            let exp = self.parse_unary()?;
            return Some(Expr::pow(base, exp));
        }
        Some(base)
    }
    fn parse_unary(&mut self) -> Option<Expr> {
        if self.eat('-') {
            let a = self.parse_unary()?;
            return Some(Expr::neg(a));
        }
        if self.eat('!') {
            let a = self.parse_unary()?;
            return Some(Expr::not_(a));
        }
        self.parse_primary()
    }
    fn parse_primary(&mut self) -> Option<Expr> {
        if self.eat('(') {
            let e = self.parse_or()?;
            if !self.eat(')') {
                return None;
            }
            return Some(e);
        }
        let c = self.peek()?;
        if c.is_alphabetic() || c == '_' {
            let start = self.pos;
            while let Some(ch) = self.peek() {
                if ch.is_alphanumeric() || ch == '_' {
                    self.pos += 1;
                } else {
                    break;
                }
            }
            let name: String = self.chars[start..self.pos].iter().collect();
            if name == "t" {
                return Some(Expr::Time);
            }
            if let Some(rest) = name.strip_prefix('s') {
                if let Ok(i) = rest.parse::<usize>() {
                    return Some(Expr::Var(i));
                }
            }
            if let Some(rest) = name.strip_prefix('p') {
                if let Ok(i) = rest.parse::<usize>() {
                    return Some(Expr::Param(i));
                }
            }
            if !self.eat('(') {
                return None;
            }
            let a = self.parse_or()?;
            let b = if self.eat(',') { Some(self.parse_or()?) } else { None };
            if !self.eat(')') {
                return None;
            }
            return Some(match (name.as_str(), b) {
                ("min", Some(b)) => Expr::min_(a, b),
                ("max", Some(b)) => Expr::max_(a, b),
                ("pow", Some(b)) => Expr::pow(a, b),
                ("exp", None) => Expr::Exp(Box::new(a)),
                ("ln", None) => Expr::Ln(Box::new(a)),
                ("log10", None) => Expr::Log10(Box::new(a)),
                ("sqrt", None) => Expr::Sqrt(Box::new(a)),
                ("abs", None) => Expr::Abs(Box::new(a)),
                ("floor", None) => Expr::Floor(Box::new(a)),
                ("ceil", None) => Expr::Ceil(Box::new(a)),
                _ => return None,
            });
        }
        let start = self.pos;
        while let Some(ch) = self.peek() {
            if ch.is_ascii_digit() || ch == '.' || ch == 'e' || ch == 'E' || ch == '+' || ch == '-' {
                if ch == '+' || ch == '-' {
                    if start < self.pos {
                        let prev = self.chars[self.pos - 1];
                        if prev != 'e' && prev != 'E' {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                self.pos += 1;
            } else {
                break;
            }
        }
        if start == self.pos {
            return None;
        }
        let lit: String = self.chars[start..self.pos].iter().collect();
        lit.parse::<f64>().ok().map(Expr::Const)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_value() {
        let e = Expr::mul(Expr::add(Expr::k(2.0), Expr::k(3.0)), Expr::var(0));
        assert_eq!(e.value(&[4.0], &[], 0.0), 20.0);
    }

    #[test]
    fn comparison_value() {
        let e = Expr::gt(Expr::var(0), Expr::param(0));
        assert_eq!(e.value(&[2.0], &[1.0], 0.0), 1.0);
        assert_eq!(e.value(&[0.5], &[1.0], 0.0), 0.0);
    }

    #[test]
    fn boolean_ops() {
        let e = Expr::and(
            Expr::gt(Expr::var(0), Expr::k(1.0)),
            Expr::lt(Expr::var(1), Expr::k(5.0)),
        );
        assert_eq!(e.value(&[2.0, 3.0], &[], 0.0), 1.0);
        assert_eq!(e.value(&[2.0, 9.0], &[], 0.0), 0.0);
    }

    #[test]
    fn time_dependency() {
        let e = Expr::ge(Expr::Time, Expr::k(5.0));
        assert_eq!(e.value(&[], &[], 4.999), 0.0);
        assert_eq!(e.value(&[], &[], 5.0), 1.0);
        assert_eq!(e.value(&[], &[], 6.0), 1.0);
    }

    #[test]
    fn divide_by_zero_is_zero() {
        let e = Expr::div(Expr::k(1.0), Expr::k(0.0));
        assert_eq!(e.value(&[], &[], 0.0), 0.0);
    }

    #[test]
    fn dependency_tracking() {
        let e = Expr::add(
            Expr::var(0),
            Expr::mul(Expr::var(2), Expr::param(1)),
        );
        let vd = e.var_deps();
        assert!(vd.contains(&0));
        assert!(vd.contains(&2));
        assert!(!vd.contains(&1));
        let pd = e.param_deps();
        assert!(pd.contains(&1));
        assert!(!pd.contains(&0));
    }

    #[test]
    fn round_trip_parse() {
        let cases = vec![
            Expr::add(Expr::var(0), Expr::param(1)),
            Expr::mul(Expr::var(0), Expr::k(2.5)),
            Expr::gt(Expr::var(0), Expr::k(5.0)),
            Expr::and(
                Expr::gt(Expr::var(0), Expr::k(1.0)),
                Expr::lt(Expr::var(1), Expr::Time),
            ),
            Expr::min_(Expr::var(0), Expr::k(1.0)),
            Expr::Exp(Box::new(Expr::var(0))),
            Expr::sub(Expr::k(1.0), Expr::div(Expr::var(0), Expr::var(1))),
            Expr::pow(Expr::var(0), Expr::k(3.0)),
            Expr::neg(Expr::var(0)),
        ];
        for e in cases {
            let s = e.to_string_compact();
            let back = Expr::parse(&s)
                .unwrap_or_else(|| panic!("could not parse `{s}`"));
            let y = vec![2.0, 3.0];
            let p = vec![1.0, 4.0];
            for &t in &[0.0_f64, 1.5, 7.0] {
                let lhs = e.value(&y, &p, t);
                let rhs = back.value(&y, &p, t);
                assert!(
                    (lhs - rhs).abs() < 1e-9,
                    "{s}: {lhs} vs {rhs} at t={t}",
                );
            }
        }
    }

    #[test]
    fn parse_handles_whitespace_and_double_ops() {
        let e = Expr::parse("(s0 > 1) && (s1 < 5)").unwrap();
        assert_eq!(e.value(&[2.0, 3.0], &[], 0.0), 1.0);
        assert_eq!(e.value(&[2.0, 9.0], &[], 0.0), 0.0);
    }

    #[test]
    fn parse_rejects_garbage() {
        assert!(Expr::parse("((s0+").is_none());
        assert!(Expr::parse("foo(s0)").is_none());
        assert!(Expr::parse("").is_none());
    }
}

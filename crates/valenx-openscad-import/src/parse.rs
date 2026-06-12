//! Recursive-descent parser for the OpenSCAD subset.

use crate::ast::{Ast, BinOp};
use crate::error::OpenScadError;
use crate::lex::Token;

/// Positional + named argument list returned by the parser.
type ArgList = (Vec<Ast>, Vec<(String, Ast)>);

/// Parse a token stream into a list of top-level statements.
///
/// Each statement becomes one [`Ast`] node — assignments, calls, or
/// blocks. The interpreter folds the list into a single expression by
/// implicit unioning of the resulting solids.
pub fn parse(tokens: &[Token]) -> Result<Vec<Ast>, OpenScadError> {
    let mut p = Parser {
        toks: tokens,
        i: 0,
        depth: 0,
    };
    let mut out = Vec::new();
    while p.peek().is_some() {
        let stmt = p.parse_statement()?;
        out.push(stmt);
    }
    Ok(out)
}

/// Cap the parser's recursion depth so a pathologically nested source
/// (deeply nested `(`/`[` expressions or `{}` blocks) returns an error
/// instead of overflowing the stack.
const MAX_PARSE_DEPTH: usize = 200;

struct Parser<'a> {
    toks: &'a [Token],
    i: usize,
    /// Current recursion depth, bounded by `MAX_PARSE_DEPTH`.
    depth: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&Token> {
        self.toks.get(self.i)
    }

    fn bump(&mut self) -> Option<&Token> {
        let t = self.toks.get(self.i);
        if t.is_some() {
            self.i += 1;
        }
        t
    }

    fn expect(&mut self, want: &Token) -> Result<(), OpenScadError> {
        if let Some(t) = self.peek() {
            if std::mem::discriminant(t) == std::mem::discriminant(want) {
                self.i += 1;
                return Ok(());
            }
            return Err(OpenScadError::Parse {
                reason: format!("expected {want:?}, got {t:?}"),
            });
        }
        Err(OpenScadError::Parse {
            reason: format!("expected {want:?}, got EOF"),
        })
    }

    fn parse_statement(&mut self) -> Result<Ast, OpenScadError> {
        if self.depth >= MAX_PARSE_DEPTH {
            return Err(OpenScadError::Parse {
                reason: format!("statement nesting exceeds the {MAX_PARSE_DEPTH}-deep cap"),
            });
        }
        self.depth += 1;
        let r = self.parse_statement_impl();
        self.depth -= 1;
        r
    }

    fn parse_statement_impl(&mut self) -> Result<Ast, OpenScadError> {
        // `ident = expr ;` is an assignment.
        if let (Some(Token::Ident(name)), Some(Token::Eq)) =
            (self.toks.get(self.i), self.toks.get(self.i + 1))
        {
            let name = name.clone();
            self.i += 2; // ident =
            let value = self.parse_expr()?;
            self.expect(&Token::Semi)?;
            return Ok(Ast::Assign(name, Box::new(value)));
        }
        // `for(var = [lo : hi]) body`
        if let Some(Token::Ident(s)) = self.peek() {
            if s == "for" {
                self.i += 1; // for
                self.expect(&Token::LParen)?;
                let var = match self.bump() {
                    Some(Token::Ident(v)) => v.clone(),
                    other => {
                        return Err(OpenScadError::Parse {
                            reason: format!("expected loop var ident, got {other:?}"),
                        })
                    }
                };
                self.expect(&Token::Eq)?;
                // Range literal: `[lo : hi]` or `[lo : step : hi]`.
                self.expect(&Token::LBracket)?;
                let lo = self.parse_expr()?;
                self.expect(&Token::Colon)?;
                let mid = self.parse_expr()?;
                let (step, hi) = if let Some(Token::Colon) = self.peek() {
                    self.i += 1; // colon
                    let hi = self.parse_expr()?;
                    (mid, hi)
                } else {
                    (Ast::Number(1.0), mid)
                };
                self.expect(&Token::RBracket)?;
                self.expect(&Token::RParen)?;
                let body = self.parse_statement()?;
                return Ok(Ast::For {
                    var,
                    lo: Box::new(lo),
                    step: Box::new(step),
                    hi: Box::new(hi),
                    body: Box::new(body),
                });
            }
        }
        // Bare `{ ... }` block.
        if matches!(self.peek(), Some(Token::LBrace)) {
            return self.parse_block();
        }
        // Otherwise: call / transform statement.
        self.parse_call_stmt()
    }

    fn parse_block(&mut self) -> Result<Ast, OpenScadError> {
        self.expect(&Token::LBrace)?;
        let mut children = Vec::new();
        while !matches!(self.peek(), Some(Token::RBrace) | None) {
            children.push(self.parse_statement()?);
        }
        self.expect(&Token::RBrace)?;
        Ok(Ast::Block(children))
    }

    fn parse_call_stmt(&mut self) -> Result<Ast, OpenScadError> {
        let name = match self.bump() {
            Some(Token::Ident(n)) => n.clone(),
            other => {
                return Err(OpenScadError::Parse {
                    reason: format!("expected ident at start of statement, got {other:?}"),
                })
            }
        };
        self.expect(&Token::LParen)?;
        let (positional, named) = self.parse_arglist()?;
        self.expect(&Token::RParen)?;
        // Either followed by `;` (leaf call) or by another statement /
        // block (transform call carrying children).
        let children = match self.peek() {
            Some(Token::Semi) => {
                self.i += 1;
                Vec::new()
            }
            Some(Token::LBrace) => {
                let blk = self.parse_block()?;
                if let Ast::Block(ch) = blk {
                    ch
                } else {
                    unreachable!()
                }
            }
            _ => {
                let child = self.parse_statement()?;
                vec![child]
            }
        };
        Ok(Ast::Call {
            name,
            positional,
            named,
            children,
        })
    }

    fn parse_arglist(&mut self) -> Result<ArgList, OpenScadError> {
        let mut positional = Vec::new();
        let mut named = Vec::new();
        if matches!(self.peek(), Some(Token::RParen)) {
            return Ok((positional, named));
        }
        loop {
            // Lookahead for `ident = expr` named arg.
            if let (Some(Token::Ident(n)), Some(Token::Eq)) =
                (self.toks.get(self.i), self.toks.get(self.i + 1))
            {
                let name = n.clone();
                self.i += 2;
                let val = self.parse_expr()?;
                named.push((name, val));
            } else {
                let e = self.parse_expr()?;
                positional.push(e);
            }
            match self.peek() {
                Some(Token::Comma) => {
                    self.i += 1;
                }
                _ => break,
            }
        }
        Ok((positional, named))
    }

    /// Expression parser — Pratt-style with two precedence levels:
    /// `+ -` (low) then `* /` (high).
    fn parse_expr(&mut self) -> Result<Ast, OpenScadError> {
        if self.depth >= MAX_PARSE_DEPTH {
            return Err(OpenScadError::Parse {
                reason: format!("expression nesting exceeds the {MAX_PARSE_DEPTH}-deep cap"),
            });
        }
        self.depth += 1;
        let r = self.parse_expr_impl();
        self.depth -= 1;
        r
    }

    fn parse_expr_impl(&mut self) -> Result<Ast, OpenScadError> {
        let mut lhs = self.parse_mul()?;
        while let Some(t) = self.peek() {
            let op = match t {
                Token::Plus => BinOp::Add,
                Token::Minus => BinOp::Sub,
                _ => break,
            };
            self.i += 1;
            let rhs = self.parse_mul()?;
            lhs = Ast::BinaryOp(Box::new(lhs), op, Box::new(rhs));
        }
        Ok(lhs)
    }

    fn parse_mul(&mut self) -> Result<Ast, OpenScadError> {
        let mut lhs = self.parse_unary()?;
        while let Some(t) = self.peek() {
            let op = match t {
                Token::Star => BinOp::Mul,
                Token::Slash => BinOp::Div,
                _ => break,
            };
            self.i += 1;
            let rhs = self.parse_unary()?;
            lhs = Ast::BinaryOp(Box::new(lhs), op, Box::new(rhs));
        }
        Ok(lhs)
    }

    fn parse_unary(&mut self) -> Result<Ast, OpenScadError> {
        if matches!(self.peek(), Some(Token::Minus)) {
            self.i += 1;
            let rhs = self.parse_atom()?;
            return Ok(Ast::Negate(Box::new(rhs)));
        }
        self.parse_atom()
    }

    fn parse_atom(&mut self) -> Result<Ast, OpenScadError> {
        match self.bump().cloned() {
            Some(Token::Number(v)) => Ok(Ast::Number(v)),
            Some(Token::Ident(s)) => {
                // Could be a function call inside an expression
                // (e.g. `sin(x)`), but for the v1 subset we only need
                // bare identifier references — extending to in-expr
                // calls is a parser-only patch.
                Ok(Ast::Ident(s))
            }
            Some(Token::LBracket) => {
                let mut items = Vec::new();
                if !matches!(self.peek(), Some(Token::RBracket)) {
                    loop {
                        items.push(self.parse_expr()?);
                        match self.peek() {
                            Some(Token::Comma) => {
                                self.i += 1;
                            }
                            _ => break,
                        }
                    }
                }
                self.expect(&Token::RBracket)?;
                Ok(Ast::Vector(items))
            }
            Some(Token::LParen) => {
                let e = self.parse_expr()?;
                self.expect(&Token::RParen)?;
                Ok(e)
            }
            other => Err(OpenScadError::Parse {
                reason: format!("expected expression atom, got {other:?}"),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lex::lex;

    #[test]
    fn parse_cube_call() {
        let toks = lex("cube([1, 2, 3]);").expect("lex");
        let ast = parse(&toks).expect("parse");
        assert_eq!(ast.len(), 1);
        if let Ast::Call {
            name, positional, ..
        } = &ast[0]
        {
            assert_eq!(name, "cube");
            assert_eq!(positional.len(), 1);
        } else {
            panic!("expected call");
        }
    }

    #[test]
    fn deeply_nested_expression_errors_not_stack_overflow() {
        // ~300 nested parens exceeds MAX_PARSE_DEPTH; the parser must return an
        // error rather than overflow the stack on a malicious .scad source.
        let src = format!("x={}1{};", "(".repeat(300), ")".repeat(300));
        let toks = lex(&src).expect("lex");
        assert!(parse(&toks).is_err());
    }

    #[test]
    fn parse_translated_block() {
        let toks = lex("translate([10, 0, 0]) { cube([1, 1, 1]); sphere(r = 2); }").expect("lex");
        let ast = parse(&toks).expect("parse");
        assert_eq!(ast.len(), 1);
        if let Ast::Call { name, children, .. } = &ast[0] {
            assert_eq!(name, "translate");
            assert_eq!(children.len(), 2);
        } else {
            panic!("expected outer call");
        }
    }

    #[test]
    fn parse_for_loop() {
        let toks = lex("for(i = [0 : 2]) translate([i * 5, 0, 0]) cube([1,1,1]);").expect("lex");
        let ast = parse(&toks).expect("parse");
        assert_eq!(ast.len(), 1);
        assert!(matches!(ast[0], Ast::For { .. }));
    }

    #[test]
    fn parse_assignment() {
        let toks = lex("x = 5;").expect("lex");
        let ast = parse(&toks).expect("parse");
        assert_eq!(ast.len(), 1);
        if let Ast::Assign(name, _) = &ast[0] {
            assert_eq!(name, "x");
        } else {
            panic!();
        }
    }
}

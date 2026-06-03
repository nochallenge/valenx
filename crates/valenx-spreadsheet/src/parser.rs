//! Hand-rolled lexer + recursive-descent parser for formulas.
//!
//! Grammar (lowest → highest precedence):
//!
//! ```text
//! expr    := term (('+' | '-') term)*
//! term    := factor (('*' | '/') factor)*
//! factor  := power ('^' power)*           -- right-associative
//! power   := '-' power | primary
//! primary := NUMBER
//!          | NAME ('.' (NAME | NUMBER))?  -- cell ref or named constant
//!          | NAME '(' args ')'             -- function call
//!          | '(' expr ')'
//! args    := expr (',' expr)*
//! ```
//!
//! The lexer is character-by-character to keep dependencies to zero.
//! Whitespace is skipped between tokens. Identifiers accept the usual
//! ASCII alphanumeric + `_`. Numbers accept `123`, `1.5`, `1.5e10`,
//! `1.5e-3`.

use crate::cell::CellRef;
use crate::error::SpreadsheetError;
use crate::formula::{BinOp, Expr, Token, UnOp};

/// Public entry — parse a formula source string into an [`Expr`].
///
/// A leading `=` is allowed and stripped (matches the convention where
/// the user types `=A1+B2` in a cell).
///
/// # Errors
///
/// Returns [`SpreadsheetError::ParseError`] with position info when
/// the lexer or parser hits an unexpected character / token.
pub fn parse(input: &str) -> Result<Expr, SpreadsheetError> {
    let src = input.strip_prefix('=').unwrap_or(input);
    let tokens = tokenize(src)?;
    let mut p = Parser::new(input, tokens);
    let expr = p.parse_expr()?;
    p.expect_eof()?;
    Ok(expr)
}

/// Lex a formula source string into a flat token vector. Each token
/// also carries the character offset where it begins in the source —
/// the parser uses that for error reporting.
pub fn tokenize(src: &str) -> Result<Vec<TokenAt>, SpreadsheetError> {
    let chars: Vec<char> = src.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        let start = i;
        match c {
            '+' => {
                tokens.push(TokenAt {
                    tok: Token::Plus,
                    pos: start,
                });
                i += 1;
            }
            '-' => {
                tokens.push(TokenAt {
                    tok: Token::Minus,
                    pos: start,
                });
                i += 1;
            }
            '*' => {
                tokens.push(TokenAt {
                    tok: Token::Star,
                    pos: start,
                });
                i += 1;
            }
            '/' => {
                tokens.push(TokenAt {
                    tok: Token::Slash,
                    pos: start,
                });
                i += 1;
            }
            '^' => {
                tokens.push(TokenAt {
                    tok: Token::Caret,
                    pos: start,
                });
                i += 1;
            }
            '(' => {
                tokens.push(TokenAt {
                    tok: Token::LParen,
                    pos: start,
                });
                i += 1;
            }
            ')' => {
                tokens.push(TokenAt {
                    tok: Token::RParen,
                    pos: start,
                });
                i += 1;
            }
            ',' => {
                tokens.push(TokenAt {
                    tok: Token::Comma,
                    pos: start,
                });
                i += 1;
            }
            '.' => {
                tokens.push(TokenAt {
                    tok: Token::Dot,
                    pos: start,
                });
                i += 1;
            }
            d if d.is_ascii_digit()
                || (d == '.' && i + 1 < chars.len() && chars[i + 1].is_ascii_digit()) =>
            {
                // Number — digits, optional dot, optional exponent.
                // Note: a bare `.` is handled above as a Dot token; a
                // dot followed by a digit starts a number like `.5`.
                let mut j = i;
                let mut seen_dot = false;
                let mut seen_e = false;
                while j < chars.len() {
                    let ch = chars[j];
                    if ch.is_ascii_digit() {
                        j += 1;
                    } else if ch == '.' && !seen_dot && !seen_e {
                        seen_dot = true;
                        j += 1;
                    } else if (ch == 'e' || ch == 'E') && !seen_e {
                        seen_e = true;
                        j += 1;
                        if j < chars.len() && (chars[j] == '+' || chars[j] == '-') {
                            j += 1;
                        }
                    } else {
                        break;
                    }
                }
                let text: String = chars[i..j].iter().collect();
                let n: f64 = text.parse().map_err(|_| SpreadsheetError::ParseError {
                    input: src.to_string(),
                    position: start,
                    reason: format!("could not parse number `{text}`"),
                })?;
                tokens.push(TokenAt {
                    tok: Token::Number(n),
                    pos: start,
                });
                i = j;
            }
            a if a.is_ascii_alphabetic() || a == '_' => {
                // Identifier — alphanumeric + `_`.
                let mut j = i;
                while j < chars.len() && (chars[j].is_ascii_alphanumeric() || chars[j] == '_') {
                    j += 1;
                }
                let name: String = chars[i..j].iter().collect();
                tokens.push(TokenAt {
                    tok: Token::Ident(name),
                    pos: start,
                });
                i = j;
            }
            other => {
                return Err(SpreadsheetError::ParseError {
                    input: src.to_string(),
                    position: start,
                    reason: format!("unexpected character `{other}`"),
                });
            }
        }
    }
    tokens.push(TokenAt {
        tok: Token::Eof,
        pos: chars.len(),
    });
    Ok(tokens)
}

/// A token plus its starting character offset in the source. Used by
/// the parser's error messages.
#[derive(Clone, Debug)]
pub struct TokenAt {
    /// The token payload.
    pub tok: Token,
    /// Zero-based character offset where the token starts.
    pub pos: usize,
}

/// Cap the parser's recursion depth so a pathologically nested
/// expression can't stack-overflow the host. Round-8 wired this into
/// the LParen arm of [`parse_primary`]; round-12 extended it to
/// [`parse_factor`]'s `^` recursion and [`parse_power`]'s unary `-`
/// recursion, since those two paths recurse without ever consuming a
/// paren. All three recursive entry points now share this cap. 100
/// levels is far past anything a human would type into a cell and
/// well below the default thread stack size — the cap fires long
/// before the OS does.
const MAX_PARSE_DEPTH: usize = 100;

struct Parser<'a> {
    /// Original (un-stripped) input — retained so error messages can
    /// echo what the user typed.
    original: &'a str,
    tokens: Vec<TokenAt>,
    pos: usize,
    /// Recursion-depth counter. Round-8 wired the cap into
    /// [`parse_primary`]'s LParen arm; round-12 extended it to
    /// [`parse_factor`] (right-associative `^`) and [`parse_power`]
    /// (unary `-`), because those two arms also recurse without ever
    /// consuming an LParen — so the LParen-only check could be
    /// bypassed by chains like `=-----...x` or `=2^2^2^...^2`. All
    /// three recursive entry points now bump-and-check uniformly.
    depth: usize,
}

impl<'a> Parser<'a> {
    fn new(original: &'a str, tokens: Vec<TokenAt>) -> Self {
        Self {
            original,
            tokens,
            pos: 0,
            depth: 0,
        }
    }

    fn peek(&self) -> &TokenAt {
        &self.tokens[self.pos]
    }

    fn bump(&mut self) -> TokenAt {
        let t = self.tokens[self.pos].clone();
        if !matches!(t.tok, Token::Eof) {
            self.pos += 1;
        }
        t
    }

    fn err(&self, pos: usize, reason: impl Into<String>) -> SpreadsheetError {
        SpreadsheetError::ParseError {
            input: self.original.to_string(),
            position: pos,
            reason: reason.into(),
        }
    }

    fn expect_eof(&self) -> Result<(), SpreadsheetError> {
        let t = self.peek();
        if !matches!(t.tok, Token::Eof) {
            return Err(self.err(
                t.pos,
                format!("unexpected trailing {} after expression", t.tok.label()),
            ));
        }
        Ok(())
    }

    fn parse_expr(&mut self) -> Result<Expr, SpreadsheetError> {
        let mut lhs = self.parse_term()?;
        loop {
            match self.peek().tok {
                Token::Plus => {
                    self.bump();
                    let rhs = self.parse_term()?;
                    lhs = Expr::binary(BinOp::Add, lhs, rhs);
                }
                Token::Minus => {
                    self.bump();
                    let rhs = self.parse_term()?;
                    lhs = Expr::binary(BinOp::Sub, lhs, rhs);
                }
                _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_term(&mut self) -> Result<Expr, SpreadsheetError> {
        let mut lhs = self.parse_factor()?;
        loop {
            match self.peek().tok {
                Token::Star => {
                    self.bump();
                    let rhs = self.parse_factor()?;
                    lhs = Expr::binary(BinOp::Mul, lhs, rhs);
                }
                Token::Slash => {
                    self.bump();
                    let rhs = self.parse_factor()?;
                    lhs = Expr::binary(BinOp::Div, lhs, rhs);
                }
                _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_factor(&mut self) -> Result<Expr, SpreadsheetError> {
        // factor := power ('^' power)*  — right-associative.
        let lhs = self.parse_power()?;
        if matches!(self.peek().tok, Token::Caret) {
            // Round-12: cap recursion depth here. `^` is right-
            // associative so we recurse into `parse_factor` on the
            // RHS — a long chain `2^2^2^...^2` would otherwise blow
            // the stack despite the LParen-arm cap, because no
            // parens are ever consumed.
            let caret = self.peek().clone();
            if self.depth >= MAX_PARSE_DEPTH {
                return Err(self.err(
                    caret.pos,
                    format!(
                        "expression nesting exceeds the {MAX_PARSE_DEPTH}-deep cap"
                    ),
                ));
            }
            self.bump();
            self.depth += 1;
            let rhs = self.parse_factor();
            self.depth -= 1;
            let rhs = rhs?;
            return Ok(Expr::binary(BinOp::Pow, lhs, rhs));
        }
        Ok(lhs)
    }

    fn parse_power(&mut self) -> Result<Expr, SpreadsheetError> {
        if matches!(self.peek().tok, Token::Minus) {
            // Round-12: cap recursion depth here. Unary `-` recurses
            // into `parse_power` — `=-----...x` would otherwise blow
            // the stack despite the LParen-arm cap, because no
            // parens are ever consumed.
            let minus = self.peek().clone();
            if self.depth >= MAX_PARSE_DEPTH {
                return Err(self.err(
                    minus.pos,
                    format!(
                        "expression nesting exceeds the {MAX_PARSE_DEPTH}-deep cap"
                    ),
                ));
            }
            self.bump();
            self.depth += 1;
            let rhs = self.parse_power();
            self.depth -= 1;
            let rhs = rhs?;
            return Ok(Expr::unary(UnOp::Neg, rhs));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Expr, SpreadsheetError> {
        let t = self.bump();
        match t.tok {
            Token::Number(n) => Ok(Expr::Number(n)),
            Token::LParen => {
                // Round-8: cap recursion depth here — LParen is the
                // only arm in `parse_primary` that re-enters
                // `parse_expr` unbounded. Without this cap a long
                // string of `((((` characters would recurse until
                // the OS stack guard kills the thread.
                if self.depth >= MAX_PARSE_DEPTH {
                    return Err(self.err(
                        t.pos,
                        format!(
                            "expression nesting exceeds the {MAX_PARSE_DEPTH}-deep cap"
                        ),
                    ));
                }
                self.depth += 1;
                let e = self.parse_expr();
                self.depth -= 1;
                let e = e?;
                let close = self.bump();
                if !matches!(close.tok, Token::RParen) {
                    return Err(self.err(
                        close.pos,
                        format!("expected `)`, got {}", close.tok.label()),
                    ));
                }
                Ok(e)
            }
            Token::Ident(name) => {
                // 3 possibilities:
                //   1) NAME '.' (NAME | NUMBER)       — cell ref `Sheet.A1`
                //   2) NAME '(' args ')'              — function call
                //   3) NAME                            — bare named constant
                if matches!(self.peek().tok, Token::Dot) {
                    self.bump();
                    let after = self.bump();
                    let a1 = match after.tok {
                        Token::Ident(s) => s,
                        Token::Number(n) => {
                            // Allow `Sheet.A1` to parse where `A1` was
                            // tokenised as `A` then `1` — except the
                            // lexer eats `A1` as a single ident because
                            // identifiers accept alphanumeric. So the
                            // only path here is e.g. `Sheet.123` which
                            // should error.
                            return Err(self.err(
                                after.pos,
                                format!("expected cell coordinate after `.`, got number `{n}`"),
                            ));
                        }
                        other => {
                            return Err(self.err(
                                after.pos,
                                format!(
                                    "expected cell coordinate after `.`, got {}",
                                    other.label()
                                ),
                            ));
                        }
                    };
                    let r = CellRef::from_a1(name, &a1)
                        .map_err(|e| self.err(t.pos, format!("bad cell ref: {e}")))?;
                    return Ok(Expr::Ref(r));
                }
                if matches!(self.peek().tok, Token::LParen) {
                    self.bump();
                    let mut args = Vec::new();
                    if !matches!(self.peek().tok, Token::RParen) {
                        args.push(self.parse_expr()?);
                        while matches!(self.peek().tok, Token::Comma) {
                            self.bump();
                            args.push(self.parse_expr()?);
                        }
                    }
                    let close = self.bump();
                    if !matches!(close.tok, Token::RParen) {
                        return Err(self.err(
                            close.pos,
                            format!(
                                "expected `)` or `,` in call to `{name}`, got {}",
                                close.tok.label()
                            ),
                        ));
                    }
                    return Ok(Expr::Call(name, args));
                }
                Ok(Expr::Name(name))
            }
            Token::Eof => Err(self.err(
                t.pos,
                "unexpected end of input — expected expression".to_string(),
            )),
            other => Err(self.err(
                t.pos,
                format!("unexpected {} — expected expression", other.label()),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_ok(s: &str) -> Expr {
        match parse(s) {
            Ok(e) => e,
            Err(err) => panic!("parse({s:?}) failed: {err:?}"),
        }
    }

    #[test]
    fn tokenize_number() {
        let toks = tokenize("3.14").unwrap();
        assert_eq!(toks.len(), 2); // [Number, Eof]
        assert!(matches!(toks[0].tok, Token::Number(_)));
    }

    #[test]
    fn tokenize_skips_whitespace() {
        let toks = tokenize("  1  +  2  ").unwrap();
        let kinds: Vec<_> = toks.iter().map(|t| t.tok.label()).collect();
        assert_eq!(kinds, vec!["number", "`+`", "number", "end of input"]);
    }

    #[test]
    fn tokenize_compound() {
        let toks = tokenize("Sheet1.A1 + sin(pi)").unwrap();
        let kinds: Vec<_> = toks.iter().map(|t| t.tok.label()).collect();
        assert_eq!(
            kinds,
            vec![
                "identifier", // Sheet1
                "`.`",
                "identifier", // A1
                "`+`",
                "identifier", // sin
                "`(`",
                "identifier", // pi
                "`)`",
                "end of input",
            ]
        );
    }

    #[test]
    fn tokenize_scientific_number() {
        let toks = tokenize("1.5e-3").unwrap();
        assert!(matches!(toks[0].tok, Token::Number(n) if (n - 1.5e-3).abs() < 1e-15));
    }

    #[test]
    fn tokenize_unknown_char_errors() {
        let e = tokenize("1 @ 2").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.parse_error");
    }

    #[test]
    fn parse_atom_number() {
        assert_eq!(parse_ok("42"), Expr::Number(42.0));
    }

    #[test]
    fn parse_addition() {
        assert_eq!(
            parse_ok("1+2"),
            Expr::binary(BinOp::Add, Expr::Number(1.0), Expr::Number(2.0))
        );
    }

    #[test]
    fn parse_precedence() {
        // 3*4+5 → (3*4) + 5 → Binary(Add, Binary(Mul, 3, 4), 5)
        let e = parse_ok("3*4+5");
        let expected = Expr::binary(
            BinOp::Add,
            Expr::binary(BinOp::Mul, Expr::Number(3.0), Expr::Number(4.0)),
            Expr::Number(5.0),
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_parens_override() {
        // (1+2)*3 → Binary(Mul, Binary(Add, 1, 2), 3)
        let e = parse_ok("(1+2)*3");
        let expected = Expr::binary(
            BinOp::Mul,
            Expr::binary(BinOp::Add, Expr::Number(1.0), Expr::Number(2.0)),
            Expr::Number(3.0),
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_function_call() {
        // sin(pi/2) → Call("sin", [Binary(Div, Name("pi"), 2)])
        let e = parse_ok("sin(pi/2)");
        let expected = Expr::Call(
            "sin".into(),
            vec![Expr::binary(
                BinOp::Div,
                Expr::Name("pi".into()),
                Expr::Number(2.0),
            )],
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_cell_ref() {
        let e = parse_ok("Sheet1.A1 + 2");
        let expected = Expr::binary(
            BinOp::Add,
            Expr::Ref(CellRef::parse("Sheet1.A1").unwrap()),
            Expr::Number(2.0),
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_unary_minus() {
        // -3+1 = -2
        let e = parse_ok("-3+1");
        let expected = Expr::binary(
            BinOp::Add,
            Expr::unary(UnOp::Neg, Expr::Number(3.0)),
            Expr::Number(1.0),
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_power_right_associative() {
        // 2^3^2 → 2^(3^2)
        let e = parse_ok("2^3^2");
        let expected = Expr::binary(
            BinOp::Pow,
            Expr::Number(2.0),
            Expr::binary(BinOp::Pow, Expr::Number(3.0), Expr::Number(2.0)),
        );
        assert_eq!(e, expected);
    }

    #[test]
    fn parse_leading_eq_is_stripped() {
        let e = parse_ok("=1+2");
        assert!(matches!(e, Expr::Binary(BinOp::Add, _, _)));
    }

    #[test]
    fn parse_call_no_args() {
        let e = parse_ok("now()");
        assert_eq!(e, Expr::Call("now".into(), vec![]));
    }

    #[test]
    fn parse_call_three_args() {
        let e = parse_ok("if(1, 2, 3)");
        assert_eq!(
            e,
            Expr::Call(
                "if".into(),
                vec![Expr::Number(1.0), Expr::Number(2.0), Expr::Number(3.0)]
            )
        );
    }

    #[test]
    fn parse_unclosed_paren_errors() {
        let err = parse("(1+2").unwrap_err();
        assert_eq!(err.code(), "spreadsheet.parse_error");
    }

    #[test]
    fn parse_trailing_garbage_errors() {
        let err = parse("1+2 garbage").unwrap_err();
        assert_eq!(err.code(), "spreadsheet.parse_error");
    }

    #[test]
    fn parse_empty_errors() {
        let err = parse("").unwrap_err();
        assert_eq!(err.code(), "spreadsheet.parse_error");
    }

    #[test]
    fn parse_error_reports_position() {
        // `1 @ 2` — `@` is at char 2.
        let err = parse("1 @ 2").unwrap_err();
        match err {
            SpreadsheetError::ParseError { position, .. } => assert_eq!(position, 2),
            other => panic!("expected ParseError, got {other:?}"),
        }
    }

    /// Round-12 M2: a chain of 200 leading `-` signs (unary minus)
    /// drives `parse_power -> parse_power` recursion without ever
    /// touching an LParen, so it would bypass the round-8 cap and
    /// stack-overflow the host. The depth bump in `parse_power`
    /// rejects the input with a structured error.
    #[test]
    fn parse_rejects_pathological_unary_minus_chain() {
        let depth = 200;
        let mut s = String::with_capacity(depth + 1);
        for _ in 0..depth {
            s.push('-');
        }
        s.push('1');
        let err = parse(&s).unwrap_err();
        match err {
            SpreadsheetError::ParseError { reason, .. } => {
                assert!(
                    reason.contains("nesting") || reason.contains(&MAX_PARSE_DEPTH.to_string()),
                    "expected nesting-depth message, got: {reason}"
                );
            }
            other => panic!("expected ParseError, got {other:?}"),
        }
    }

    /// Round-12 M2: a chain of 200 right-associative `^` operators
    /// drives `parse_factor -> parse_factor` recursion without ever
    /// touching an LParen, so it would bypass the round-8 cap and
    /// stack-overflow the host. The depth bump in `parse_factor`
    /// rejects the input with a structured error.
    #[test]
    fn parse_rejects_pathological_caret_chain() {
        let depth = 200;
        let mut s = String::with_capacity(depth * 3 + 1);
        s.push('2');
        for _ in 0..depth {
            s.push_str("^2");
        }
        let err = parse(&s).unwrap_err();
        match err {
            SpreadsheetError::ParseError { reason, .. } => {
                assert!(
                    reason.contains("nesting") || reason.contains(&MAX_PARSE_DEPTH.to_string()),
                    "expected nesting-depth message, got: {reason}"
                );
            }
            other => panic!("expected ParseError, got {other:?}"),
        }
    }

    #[test]
    fn parse_rejects_pathological_paren_nesting() {
        // Round-8 RED→GREEN: pre-fix, a 2000-deep nested expression
        // would recurse through `parse_expr` -> `parse_primary` and
        // overflow the default thread stack, crashing the host. The
        // MAX_PARSE_DEPTH cap rejects the input with a structured
        // error before recursion gets close to the OS limit.
        let depth = 2_000;
        let mut s = String::with_capacity(depth * 2 + 1);
        for _ in 0..depth {
            s.push('(');
        }
        s.push('1');
        for _ in 0..depth {
            s.push(')');
        }
        let err = parse(&s).unwrap_err();
        match err {
            SpreadsheetError::ParseError { reason, .. } => {
                assert!(
                    reason.contains("nesting") || reason.contains(&MAX_PARSE_DEPTH.to_string()),
                    "expected nesting-depth message, got: {reason}"
                );
            }
            other => panic!("expected ParseError, got {other:?}"),
        }
    }
}

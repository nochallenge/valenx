//! OpenSCAD lexer.

use crate::error::OpenScadError;

/// A token from an OpenSCAD source file.
#[derive(Clone, Debug, PartialEq)]
pub enum Token {
    /// Identifier or keyword.
    Ident(String),
    /// Numeric literal (always parsed as `f64`).
    Number(f64),
    /// String literal contents (without quotes).
    Str(String),
    /// `{`
    LBrace,
    /// `}`
    RBrace,
    /// `(`
    LParen,
    /// `)`
    RParen,
    /// `[`
    LBracket,
    /// `]`
    RBracket,
    /// `,`
    Comma,
    /// `;`
    Semi,
    /// `=`
    Eq,
    /// `+`
    Plus,
    /// `-`
    Minus,
    /// `*`
    Star,
    /// `/`
    Slash,
    /// `<`
    Lt,
    /// `>`
    Gt,
    /// `!`
    Bang,
    /// `:`
    Colon,
}

/// True when `s` is an OpenSCAD reserved word.
pub fn is_keyword(s: &str) -> bool {
    matches!(
        s,
        "module"
            | "function"
            | "if"
            | "else"
            | "for"
            | "intersection_for"
            | "assign"
            | "let"
            | "cube"
            | "sphere"
            | "cylinder"
            | "union"
            | "difference"
            | "intersection"
            | "translate"
            | "rotate"
            | "scale"
            | "mirror"
            | "linear_extrude"
            | "rotate_extrude"
    )
}

/// Tokenise an entire OpenSCAD source string.
///
/// Whitespace and `//` / `/* */` comments are stripped. Returns the
/// flat token list in source order — the parser owns position tracking.
pub fn lex(src: &str) -> Result<Vec<Token>, OpenScadError> {
    let mut out = Vec::new();
    let bytes = src.as_bytes();
    let n = bytes.len();
    let mut i = 0;
    while i < n {
        let c = bytes[i];
        // Whitespace.
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        // Single-line comment.
        if c == b'/' && i + 1 < n && bytes[i + 1] == b'/' {
            while i < n && bytes[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        // Block comment.
        if c == b'/' && i + 1 < n && bytes[i + 1] == b'*' {
            i += 2;
            while i + 1 < n && !(bytes[i] == b'*' && bytes[i + 1] == b'/') {
                i += 1;
            }
            i = (i + 2).min(n);
            continue;
        }
        // Single-character punctuation.
        let single = match c {
            b'{' => Some(Token::LBrace),
            b'}' => Some(Token::RBrace),
            b'(' => Some(Token::LParen),
            b')' => Some(Token::RParen),
            b'[' => Some(Token::LBracket),
            b']' => Some(Token::RBracket),
            b',' => Some(Token::Comma),
            b';' => Some(Token::Semi),
            b'=' => Some(Token::Eq),
            b'+' => Some(Token::Plus),
            b'*' => Some(Token::Star),
            b'/' => Some(Token::Slash),
            b'<' => Some(Token::Lt),
            b'>' => Some(Token::Gt),
            b'!' => Some(Token::Bang),
            b':' => Some(Token::Colon),
            _ => None,
        };
        if let Some(tok) = single {
            out.push(tok);
            i += 1;
            continue;
        }
        // `-` could be unary or binary — let the parser sort it out.
        if c == b'-' {
            out.push(Token::Minus);
            i += 1;
            continue;
        }
        // Number: integer / decimal / scientific. Always returns f64.
        if c.is_ascii_digit() || (c == b'.' && i + 1 < n && bytes[i + 1].is_ascii_digit()) {
            let start = i;
            while i < n
                && (bytes[i].is_ascii_digit()
                    || bytes[i] == b'.'
                    || bytes[i] == b'e'
                    || bytes[i] == b'E'
                    || ((bytes[i] == b'+' || bytes[i] == b'-')
                        && i > start
                        && (bytes[i - 1] == b'e' || bytes[i - 1] == b'E')))
            {
                i += 1;
            }
            let s = std::str::from_utf8(&bytes[start..i]).map_err(|_| OpenScadError::Lex {
                pos: start,
                reason: "non-utf8 in number".into(),
            })?;
            let v: f64 = s.parse().map_err(|_| OpenScadError::Lex {
                pos: start,
                reason: format!("bad number `{s}`"),
            })?;
            out.push(Token::Number(v));
            continue;
        }
        // String literal — double-quoted, no escape processing beyond
        // `\"` / `\\` (sufficient for the basic shape constructors).
        if c == b'"' {
            i += 1;
            let mut buf = String::new();
            while i < n && bytes[i] != b'"' {
                if bytes[i] == b'\\' && i + 1 < n {
                    match bytes[i + 1] {
                        b'"' => buf.push('"'),
                        b'\\' => buf.push('\\'),
                        b'n' => buf.push('\n'),
                        b't' => buf.push('\t'),
                        // Unknown escape — keep the escaped byte
                        // verbatim. `other` is a single ASCII byte
                        // (anything multibyte would have been the
                        // continuation of a UTF-8 codepoint starting
                        // BEFORE the backslash, which can't happen
                        // because `\` itself is single-byte).
                        other => buf.push(other as char),
                    }
                    i += 2;
                } else {
                    // Not an escape — copy the next full UTF-8 codepoint
                    // verbatim. `bytes[i] as char` would otherwise emit
                    // mojibake for every multibyte char in the string
                    // literal (e.g. Greek letters in comments-as-strings).
                    let ch = src[i..]
                        .chars()
                        .next()
                        .expect("i is in-bounds inside a &str");
                    buf.push(ch);
                    i += ch.len_utf8();
                }
            }
            if i >= n {
                return Err(OpenScadError::Lex {
                    pos: i,
                    reason: "unterminated string".into(),
                });
            }
            i += 1; // closing quote
            out.push(Token::Str(buf));
            continue;
        }
        // Identifier — starts with letter/underscore, then alnum/`_`.
        if c.is_ascii_alphabetic() || c == b'_' || c == b'$' {
            let start = i;
            i += 1;
            while i < n
                && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_' || bytes[i] == b'$')
            {
                i += 1;
            }
            let s = std::str::from_utf8(&bytes[start..i]).map_err(|_| OpenScadError::Lex {
                pos: start,
                reason: "non-utf8 ident".into(),
            })?;
            out.push(Token::Ident(s.to_string()));
            continue;
        }
        return Err(OpenScadError::Lex {
            pos: i,
            reason: format!("unexpected byte 0x{c:02x}"),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lex_simple_cube_call() {
        let toks = lex("cube([1, 2, 3]);").expect("ok");
        assert_eq!(
            toks,
            vec![
                Token::Ident("cube".into()),
                Token::LParen,
                Token::LBracket,
                Token::Number(1.0),
                Token::Comma,
                Token::Number(2.0),
                Token::Comma,
                Token::Number(3.0),
                Token::RBracket,
                Token::RParen,
                Token::Semi,
            ]
        );
    }

    #[test]
    fn lex_skips_comments() {
        let toks = lex("// hi\nfoo /* mid */ bar").expect("ok");
        assert_eq!(
            toks,
            vec![Token::Ident("foo".into()), Token::Ident("bar".into())]
        );
    }

    #[test]
    fn lex_handles_negative_and_scientific() {
        let toks = lex("-1.5e2").expect("ok");
        assert_eq!(toks, vec![Token::Minus, Token::Number(150.0)]);
    }

    #[test]
    fn is_keyword_matches_subset() {
        assert!(is_keyword("cube"));
        assert!(is_keyword("union"));
        assert!(is_keyword("for"));
        assert!(!is_keyword("foo"));
    }
}

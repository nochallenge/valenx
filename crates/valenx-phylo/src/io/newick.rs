//! Newick tree format — reader and writer.
//!
//! Newick encodes a rooted tree as a nested parenthesised expression:
//!
//! ```text
//! ((A:0.1,B:0.2):0.3,C:0.5);
//! ```
//!
//! Each `(...)` is an internal node; a bare name is a leaf; `:value`
//! after either is the branch length leading to its parent; a name
//! after a `)` labels the internal node (a clade name, or — by a
//! widespread convention — a bootstrap support value). The whole
//! expression ends with `;`.
//!
//! The reader here is a hand-written recursive-descent parser. It
//! handles: nested parentheses, quoted labels (`'A name'`),
//! underscores-as-spaces in unquoted labels (the historic convention),
//! branch lengths in decimal or scientific notation, internal-node
//! labels, and a trailing semicolon. Comments in `[...]` are stripped.
//! It does **not** implement the extended New Hampshire (NHX)
//! `[&&NHX:...]` annotation grammar — NHX comment blocks are discarded.

use crate::error::{PhyloError, Result};
use crate::tree::{Node, NodeId, Tree};

/// Maximum `(...)` nesting depth the Newick parser will descend before
/// rejecting the input.
///
/// `parse_subtree` recurses once per opening parenthesis, so without a
/// bound a deeply-nested input (`((((…))))`) overflows the call stack
/// and *aborts the process* — a stack overflow is not a catchable
/// panic. Real phylogenies are rarely more than ~100 levels deep (a
/// fully-pectinate "ladder" tree of N taxa is N-1 deep, and even a
/// thousand-taxon ladder is an unusual thing to feed a viewer), so a
/// 1 000-level cap rejects only pathological / malicious nesting while
/// leaving every realistic tree untouched.
///
/// The cap must fire *before* the live recursion frames exhaust the
/// stack, or the guard is useless. Measured in a debug build, this
/// parser overflows an 8 MiB stack (the default main-thread size on
/// Windows/macOS and the common Linux default) somewhere between 6 000
/// and 8 000 nested frames. 1 000 keeps a ~6× margin below that floor
/// even in debug — and release frames are smaller still — so the guard
/// is guaranteed to return an error rather than let the stack blow.
const MAX_NEWICK_DEPTH: usize = 1_000;

/// Parses a Newick string into a [`Tree`].
///
/// The returned tree is flagged [`rooted`](Tree::rooted)` = true`
/// (Newick is a rooted format); callers working with unrooted data can
/// clear the flag.
///
/// # Errors
/// [`PhyloError::Parse`] on unbalanced parentheses, an empty input, a
/// malformed branch length, or trailing junk after the `;`.
pub fn read_newick(input: &str) -> Result<Tree> {
    let cleaned = strip_comments(input);
    let trimmed = cleaned.trim();
    if trimmed.is_empty() {
        return Err(PhyloError::parse("newick", "empty input"));
    }
    let body = trimmed.strip_suffix(';').unwrap_or(trimmed);
    let mut parser = Parser {
        chars: body.chars().collect(),
        pos: 0,
        tree: Tree::building(),
    };
    let root = parser.parse_subtree(0)?;
    parser.skip_ws();
    if parser.pos != parser.chars.len() {
        return Err(PhyloError::parse(
            "newick",
            format!("unexpected trailing input at position {}", parser.pos),
        ));
    }
    parser
        .tree
        .finish_building(root, true)
        .map_err(|e| PhyloError::parse("newick", e.to_string()))
}

/// Removes `[...]` comment blocks (including NHX annotations).
fn strip_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth = 0u32;
    for c in s.chars() {
        match c {
            '[' => depth += 1,
            ']' if depth > 0 => depth -= 1,
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

struct Parser {
    chars: Vec<char>,
    pos: usize,
    tree: Tree,
}

impl Parser {
    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(c) if c.is_whitespace()) {
            self.pos += 1;
        }
    }

    /// Parses one subtree (leaf or `(...)` clade) and returns its id.
    ///
    /// `depth` is the current `(...)` nesting level; it is checked
    /// against [`MAX_NEWICK_DEPTH`] on entering a clade so a maliciously
    /// deep input is rejected with an error rather than overflowing the
    /// stack.
    fn parse_subtree(&mut self, depth: usize) -> Result<NodeId> {
        self.skip_ws();
        let id = if self.peek() == Some('(') {
            if depth >= MAX_NEWICK_DEPTH {
                return Err(PhyloError::parse(
                    "newick",
                    format!("nesting too deep (max {MAX_NEWICK_DEPTH})"),
                ));
            }
            self.pos += 1; // consume '('
            let mut children = Vec::new();
            loop {
                let child = self.parse_subtree(depth + 1)?;
                children.push(child);
                self.skip_ws();
                match self.peek() {
                    Some(',') => {
                        self.pos += 1;
                    }
                    Some(')') => {
                        self.pos += 1;
                        break;
                    }
                    other => {
                        return Err(PhyloError::parse(
                            "newick",
                            format!("expected ',' or ')', found {other:?}"),
                        ));
                    }
                }
            }
            let node = self.tree.push_node(Node {
                label: None,
                branch_length: None,
                parent: None,
                children: children.clone(),
            });
            for &c in &children {
                self.tree.node_mut(c).parent = Some(node);
            }
            node
        } else {
            // A leaf — created with no label yet; label read below.
            self.tree.push_node(Node {
                label: None,
                branch_length: None,
                parent: None,
                children: Vec::new(),
            })
        };

        // Optional label (leaf name or internal clade label).
        let label = self.read_label();
        if let Some(l) = label {
            if !l.is_empty() {
                self.tree.node_mut(id).label = Some(l);
            }
        }

        // Optional `:branch_length`.
        self.skip_ws();
        if self.peek() == Some(':') {
            self.pos += 1;
            let bl = self.read_number()?;
            self.tree.node_mut(id).branch_length = Some(bl);
        }
        Ok(id)
    }

    /// Reads a node label: a single-quoted string, or an unquoted run
    /// up to the next structural character. Returns `None` if no label
    /// is present.
    fn read_label(&mut self) -> Option<String> {
        self.skip_ws();
        match self.peek() {
            Some('\'') => {
                self.pos += 1;
                let mut s = String::new();
                while let Some(c) = self.peek() {
                    self.pos += 1;
                    if c == '\'' {
                        // A doubled '' is an escaped single quote.
                        if self.peek() == Some('\'') {
                            self.pos += 1;
                            s.push('\'');
                            continue;
                        }
                        break;
                    }
                    s.push(c);
                }
                Some(s)
            }
            Some(c) if !is_structural(c) => {
                let mut s = String::new();
                while let Some(c) = self.peek() {
                    if is_structural(c) {
                        break;
                    }
                    // Historic convention: unquoted underscore == space.
                    s.push(if c == '_' { ' ' } else { c });
                    self.pos += 1;
                }
                Some(s.trim().to_string())
            }
            _ => None,
        }
    }

    /// Reads a decimal / scientific-notation number.
    fn read_number(&mut self) -> Result<f64> {
        self.skip_ws();
        let start = self.pos;
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E') {
                self.pos += 1;
            } else {
                break;
            }
        }
        let slice: String = self.chars[start..self.pos].iter().collect();
        slice
            .parse::<f64>()
            .map_err(|_| PhyloError::parse("newick", format!("bad branch length `{slice}`")))
    }
}

/// `true` for characters that terminate an unquoted label.
fn is_structural(c: char) -> bool {
    matches!(c, '(' | ')' | ',' | ':' | ';' | '\'')
}

/// Whether an unquoted label needs single-quoting on output.
fn needs_quotes(label: &str) -> bool {
    label.is_empty()
        || label
            .chars()
            .any(|c| is_structural(c) || c.is_whitespace() || c == '[' || c == ']')
}

/// Serialises a [`Tree`] to a one-line Newick string ending with `;`.
///
/// Branch lengths are written with [`f64`]'s default formatting when
/// present and omitted when `None`. Internal-node labels are written
/// after the closing parenthesis. Labels containing structural
/// characters or whitespace are single-quoted.
pub fn write_newick(tree: &Tree) -> String {
    let mut out = String::new();
    write_subtree(tree, tree.root(), &mut out);
    out.push(';');
    out
}

fn write_subtree(tree: &Tree, id: NodeId, out: &mut String) {
    let node = tree.node(id);
    if node.is_internal() {
        out.push('(');
        for (i, &c) in node.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            write_subtree(tree, c, out);
        }
        out.push(')');
    }
    if let Some(label) = &node.label {
        if needs_quotes(label) {
            out.push('\'');
            out.push_str(&label.replace('\'', "''"));
            out.push('\'');
        } else {
            out.push_str(label);
        }
    }
    if let Some(bl) = node.branch_length {
        out.push(':');
        out.push_str(&format_branch_length(bl));
    }
}

/// Formats a branch length compactly: an integer-valued length prints
/// without a trailing `.0`, everything else uses `{}`.
fn format_branch_length(bl: f64) -> String {
    if bl.fract() == 0.0 && bl.is_finite() && bl.abs() < 1e15 {
        format!("{}", bl as i64)
    } else {
        format!("{bl}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_tree() {
        let t = read_newick("((A:0.1,B:0.2):0.3,C:0.5);").unwrap();
        assert_eq!(t.leaf_count(), 3);
        assert!(t.rooted);
        let a = t.find("A").unwrap();
        assert!((t.node(a).branch_length.unwrap() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn round_trip_preserves_topology_and_lengths() {
        let src = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        let t = read_newick(src).unwrap();
        let out = write_newick(&t);
        let t2 = read_newick(&out).unwrap();
        assert_eq!(t.leaf_labels(), t2.leaf_labels());
        assert_eq!(t.leaf_count(), t2.leaf_count());
        // Patristic A..D survives the round trip.
        let d1 = t.patristic_distance(t.find("A").unwrap(), t.find("D").unwrap());
        let d2 = t2.patristic_distance(t2.find("A").unwrap(), t2.find("D").unwrap());
        assert!((d1 - d2).abs() < 1e-12);
    }

    #[test]
    fn handles_internal_labels_and_no_branch_lengths() {
        let t = read_newick("((A,B)clade1,C);").unwrap();
        assert_eq!(t.leaf_count(), 3);
        let internal = t.find("clade1").unwrap();
        assert!(t.node(internal).is_internal());
    }

    #[test]
    fn quoted_labels_with_spaces_round_trip() {
        let t = read_newick("('Homo sapiens':1.0,'Pan troglodytes':1.0);").unwrap();
        assert!(t.find("Homo sapiens").is_some());
        let out = write_newick(&t);
        assert!(out.contains('\''), "quoted label lost: {out}");
        let t2 = read_newick(&out).unwrap();
        assert!(t2.find("Homo sapiens").is_some());
    }

    #[test]
    fn strips_comments() {
        let t = read_newick("(A[a comment],B);").unwrap();
        assert_eq!(t.leaf_count(), 2);
        assert!(t.find("A").is_some());
    }

    #[test]
    fn rejects_unbalanced_parentheses() {
        assert!(read_newick("((A,B);").is_err());
        assert!(read_newick("").is_err());
    }

    #[test]
    fn scientific_notation_branch_length() {
        let t = read_newick("(A:1.5e-3,B:2E-2);").unwrap();
        let a = t.find("A").unwrap();
        assert!((t.node(a).branch_length.unwrap() - 0.0015).abs() < 1e-12);
    }

    #[test]
    fn deeply_nested_input_is_rejected_not_a_stack_overflow() {
        // A pathological / malicious Newick string of nothing but
        // opening parentheses. `parse_subtree` recurses once per `(`, so
        // an unbounded parser overflows the stack on a deeply-nested
        // input and *aborts the process* (a stack overflow is NOT a
        // catchable panic). `cap + 1` opening parens drive the recursion
        // exactly one level past the guard, which must return a clean
        // `Err` reporting the nesting is too deep. Because the cap is far
        // below the stack-overflow floor (see `MAX_NEWICK_DEPTH`), these
        // ~1 001 frames fit comfortably on the default test stack, so the
        // guard's own return path never overflows.
        let depth = MAX_NEWICK_DEPTH + 1;
        let pathological: String = std::iter::repeat_n('(', depth).collect();
        let err = read_newick(&pathological).expect_err("deep nesting must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("nesting too deep"),
            "expected a depth-limit error, got: {msg}"
        );
    }
}

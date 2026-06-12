//! Structure-selection mini-language (v1).
//!
//! A compact selection grammar in the spirit of PyMOL / ChimeraX /
//! VMD selections. A selection is a boolean expression over atom
//! predicates joined by `and` / `or` and negated with `not`:
//!
//! ```text
//! chain A and resi 10-50 and name CA
//! (chain A or chain B) and not hetatm
//! element FE or resn HOH
//! backbone and chain A
//! ```
//!
//! ## Supported predicates
//!
//! | Token | Matches |
//! |---|---|
//! | `chain <id>[,<id>…]` | atoms in those chains |
//! | `resi <n>[-<m>][,…]` | residues by sequence number / range |
//! | `resn <name>[,…]` | residues by name (`ALA`, `HOH`, …) |
//! | `name <atom>[,…]` | atoms by name (`CA`, `N`, …) |
//! | `element <sym>[,…]` | atoms by element symbol |
//! | `hetatm` | `HETATM` atoms |
//! | `protein` / `nucleic` / `water` | polymer-class predicates |
//! | `backbone` / `sidechain` | protein backbone / non-backbone |
//! | `all` / `none` | everything / nothing |
//!
//! ## Scope of this v1
//!
//! No distance operators (`within`, `around`), no `byres` expansion,
//! no secondary-structure predicates. The parser is a hand-rolled
//! recursive-descent over a whitespace tokeniser with `()` grouping;
//! precedence is `not` > `and` > `or`.

use crate::error::{BiostructError, Result};
use crate::structure::{Atom, Model, Residue, ResidueKind};

/// A parsed, reusable atom-selection expression.
#[derive(Debug, Clone, PartialEq)]
pub struct Selection {
    expr: Expr,
}

/// One reference to a selected atom within a [`Model`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct AtomRef {
    /// Index into `Model::chains`.
    pub chain: usize,
    /// Index into `Chain::residues`.
    pub residue: usize,
    /// Index into `Residue::atoms`.
    pub atom: usize,
}

impl Selection {
    /// Parse a selection string into a reusable [`Selection`].
    pub fn parse(query: &str) -> Result<Selection> {
        let tokens = lex(query);
        if tokens.is_empty() {
            return Err(BiostructError::invalid_selection(query, "empty selection"));
        }
        let mut parser = Parser {
            tokens: &tokens,
            pos: 0,
            query,
        };
        let expr = parser.parse_or()?;
        if parser.pos != tokens.len() {
            return Err(BiostructError::invalid_selection(
                query,
                format!("unexpected token `{}`", tokens[parser.pos]),
            ));
        }
        Ok(Selection { expr })
    }

    /// Whether a single atom (with its residue / chain context)
    /// satisfies the selection.
    pub fn matches(&self, atom: &Atom, residue: &Residue, chain_id: &str) -> bool {
        test_expr(&self.expr, atom, residue, chain_id)
    }

    /// Collect every [`AtomRef`] in `model` that satisfies the
    /// selection.
    pub fn select(&self, model: &Model) -> Vec<AtomRef> {
        let mut out = Vec::new();
        for (ci, chain) in model.chains.iter().enumerate() {
            for (ri, residue) in chain.residues.iter().enumerate() {
                for (ai, atom) in residue.atoms.iter().enumerate() {
                    if self.matches(atom, residue, &chain.id) {
                        out.push(AtomRef {
                            chain: ci,
                            residue: ri,
                            atom: ai,
                        });
                    }
                }
            }
        }
        out
    }

    /// Resolve the selection to borrowed `&Atom`s.
    pub fn select_atoms<'m>(&self, model: &'m Model) -> Vec<&'m Atom> {
        self.select(model)
            .into_iter()
            .map(|r| &model.chains[r.chain].residues[r.residue].atoms[r.atom])
            .collect()
    }
}

/// One-shot convenience: parse `query` and resolve it against `model`.
pub fn select(model: &Model, query: &str) -> Result<Vec<AtomRef>> {
    Ok(Selection::parse(query)?.select(model))
}

// --- expression tree -------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Expr {
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Not(Box<Expr>),
    Pred(Pred),
}

#[derive(Debug, Clone, PartialEq)]
enum Pred {
    Chain(Vec<String>),
    /// Inclusive `(low, high)` residue-number ranges.
    Resi(Vec<(i32, i32)>),
    Resn(Vec<String>),
    Name(Vec<String>),
    Element(Vec<String>),
    Hetatm,
    Protein,
    Nucleic,
    Water,
    Backbone,
    Sidechain,
    All,
    None,
}

/// Protein-backbone atom names.
const BACKBONE: &[&str] = &["N", "CA", "C", "O", "OXT"];

/// Recursively test an atom against a parsed selection expression.
fn test_expr(expr: &Expr, atom: &Atom, residue: &Residue, chain_id: &str) -> bool {
    match expr {
        Expr::And(a, b) => {
            test_expr(a, atom, residue, chain_id) && test_expr(b, atom, residue, chain_id)
        }
        Expr::Or(a, b) => {
            test_expr(a, atom, residue, chain_id) || test_expr(b, atom, residue, chain_id)
        }
        Expr::Not(a) => !test_expr(a, atom, residue, chain_id),
        Expr::Pred(p) => test_pred(p, atom, residue, chain_id),
    }
}

/// Test an atom against a single leaf predicate.
fn test_pred(p: &Pred, atom: &Atom, residue: &Residue, chain_id: &str) -> bool {
    match p {
        Pred::Chain(ids) => ids.iter().any(|c| c.eq_ignore_ascii_case(chain_id)),
        Pred::Resi(ranges) => ranges
            .iter()
            .any(|(lo, hi)| residue.seq_num >= *lo && residue.seq_num <= *hi),
        Pred::Resn(names) => names.iter().any(|n| n.eq_ignore_ascii_case(&residue.name)),
        Pred::Name(names) => names.iter().any(|n| n.eq_ignore_ascii_case(&atom.name)),
        Pred::Element(syms) => syms.iter().any(|s| s.eq_ignore_ascii_case(&atom.element)),
        Pred::Hetatm => residue.hetatm,
        Pred::Protein => residue.kind() == ResidueKind::AminoAcid,
        Pred::Nucleic => {
            matches!(residue.kind(), ResidueKind::Dna | ResidueKind::Rna)
        }
        Pred::Water => residue.kind() == ResidueKind::Water,
        Pred::Backbone => {
            residue.kind() == ResidueKind::AminoAcid && BACKBONE.contains(&atom.name.as_str())
        }
        Pred::Sidechain => {
            residue.kind() == ResidueKind::AminoAcid
                && !BACKBONE.contains(&atom.name.as_str())
                && !atom.is_hydrogen()
        }
        Pred::All => true,
        Pred::None => false,
    }
}

// --- lexer -----------------------------------------------------------

/// Split a selection string into tokens, treating `(` and `)` as
/// standalone tokens.
fn lex(query: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current = String::new();
    for c in query.chars() {
        match c {
            '(' | ')' => {
                if !current.is_empty() {
                    tokens.push(std::mem::take(&mut current));
                }
                tokens.push(c.to_string());
            }
            c if c.is_whitespace() => {
                if !current.is_empty() {
                    tokens.push(std::mem::take(&mut current));
                }
            }
            c => current.push(c),
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

// --- recursive-descent parser ----------------------------------------

struct Parser<'a> {
    tokens: &'a [String],
    pos: usize,
    query: &'a str,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&str> {
        self.tokens.get(self.pos).map(|s| s.as_str())
    }

    fn bump(&mut self) -> Option<&str> {
        let t = self.tokens.get(self.pos).map(|s| s.as_str());
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn err(&self, reason: impl Into<String>) -> BiostructError {
        BiostructError::invalid_selection(self.query, reason)
    }

    /// `or` has the lowest precedence.
    fn parse_or(&mut self) -> Result<Expr> {
        let mut left = self.parse_and()?;
        while matches!(
            self.peek().map(|s| s.to_ascii_lowercase()).as_deref(),
            Some("or")
        ) {
            self.bump();
            let right = self.parse_and()?;
            left = Expr::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    /// `and` binds tighter than `or`.
    fn parse_and(&mut self) -> Result<Expr> {
        let mut left = self.parse_not()?;
        while matches!(
            self.peek().map(|s| s.to_ascii_lowercase()).as_deref(),
            Some("and")
        ) {
            self.bump();
            let right = self.parse_not()?;
            left = Expr::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    /// `not` binds tightest, before primaries.
    fn parse_not(&mut self) -> Result<Expr> {
        if matches!(
            self.peek().map(|s| s.to_ascii_lowercase()).as_deref(),
            Some("not")
        ) {
            self.bump();
            let inner = self.parse_not()?;
            return Ok(Expr::Not(Box::new(inner)));
        }
        self.parse_primary()
    }

    /// A parenthesised group or a single predicate.
    fn parse_primary(&mut self) -> Result<Expr> {
        match self.peek() {
            Some("(") => {
                self.bump();
                let inner = self.parse_or()?;
                match self.bump() {
                    Some(")") => Ok(inner),
                    _ => Err(self.err("unbalanced `(`")),
                }
            }
            Some(_) => self.parse_pred(),
            None => Err(self.err("expected a predicate")),
        }
    }

    /// A keyword predicate, consuming its argument tokens.
    fn parse_pred(&mut self) -> Result<Expr> {
        let kw = match self.bump() {
            Some(t) => t.to_ascii_lowercase(),
            None => return Err(self.err("expected a predicate keyword")),
        };
        let pred = match kw.as_str() {
            "all" => Pred::All,
            "none" => Pred::None,
            "hetatm" | "hetero" => Pred::Hetatm,
            "protein" | "polymer.protein" => Pred::Protein,
            "nucleic" | "polymer.nucleic" => Pred::Nucleic,
            "water" | "waters" | "solvent" => Pred::Water,
            "backbone" | "mainchain" => Pred::Backbone,
            "sidechain" | "sc" => Pred::Sidechain,
            "chain" => Pred::Chain(self.argument_list("chain")?),
            "resn" | "resname" => Pred::Resn(self.argument_list("resn")?),
            "name" => Pred::Name(self.argument_list("name")?),
            "element" | "elem" => Pred::Element(self.argument_list("element")?),
            "resi" | "resid" | "resseq" => {
                let raw = self.argument_list("resi")?;
                let mut ranges = Vec::new();
                for item in raw {
                    ranges.push(parse_range(&item).map_err(|e| self.err(e))?);
                }
                Pred::Resi(ranges)
            }
            other => {
                return Err(self.err(format!("unknown selection keyword `{other}`")));
            }
        };
        Ok(Expr::Pred(pred))
    }

    /// Consume the next token as a comma-separated argument list.
    fn argument_list(&mut self, kw: &str) -> Result<Vec<String>> {
        let tok = match self.bump() {
            Some(t) => t.to_string(),
            None => return Err(self.err(format!("`{kw}` needs an argument"))),
        };
        let items: Vec<String> = tok
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect();
        if items.is_empty() {
            return Err(self.err(format!("`{kw}` argument is empty")));
        }
        Ok(items)
    }
}

/// Parse a single `n` or `n-m` residue range item.
fn parse_range(item: &str) -> std::result::Result<(i32, i32), String> {
    // A leading '-' on a negative number complicates `n-m` splitting;
    // split on the first '-' that is not at position 0.
    if let Some(dash) = item.char_indices().skip(1).find(|(_, c)| *c == '-') {
        let (a, b) = item.split_at(dash.0);
        let b = &b[1..];
        let lo: i32 = a
            .trim()
            .parse()
            .map_err(|_| format!("bad range start in `{item}`"))?;
        let hi: i32 = b
            .trim()
            .parse()
            .map_err(|_| format!("bad range end in `{item}`"))?;
        if lo > hi {
            return Err(format!("range start > end in `{item}`"));
        }
        Ok((lo, hi))
    } else {
        let n: i32 = item
            .trim()
            .parse()
            .map_err(|_| format!("bad residue number `{item}`"))?;
        Ok((n, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain, Residue};
    use nalgebra::Point3;

    fn demo_model() -> Model {
        let mut m = Model::new(1);
        let mut a = Chain::new("A");
        for seq in 1..=3 {
            let mut r = Residue::new("ALA", seq);
            r.atoms
                .push(Atom::new("N", "N", Point3::new(0.0, 0.0, 0.0)));
            r.atoms
                .push(Atom::new("CA", "C", Point3::new(1.0, 0.0, 0.0)));
            r.atoms
                .push(Atom::new("CB", "C", Point3::new(2.0, 0.0, 0.0)));
            a.residues.push(r);
        }
        let mut b = Chain::new("B");
        let mut w = Residue::new("HOH", 100);
        w.hetatm = true;
        w.atoms
            .push(Atom::new("O", "O", Point3::new(9.0, 9.0, 9.0)));
        b.residues.push(w);
        m.chains.push(a);
        m.chains.push(b);
        m
    }

    #[test]
    fn chain_predicate() {
        let m = demo_model();
        let hits = select(&m, "chain A").unwrap();
        assert_eq!(hits.len(), 9); // 3 residues x 3 atoms
    }

    #[test]
    fn and_combination() {
        let m = demo_model();
        let hits = select(&m, "chain A and name CA").unwrap();
        assert_eq!(hits.len(), 3);
    }

    #[test]
    fn residue_range() {
        let m = demo_model();
        let hits = select(&m, "chain A and resi 2-3 and name CA").unwrap();
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn not_and_or() {
        let m = demo_model();
        let hits = select(&m, "not hetatm").unwrap();
        assert_eq!(hits.len(), 9);
        let hits = select(&m, "hetatm or name CA").unwrap();
        assert_eq!(hits.len(), 1 + 3);
    }

    #[test]
    fn parens_group() {
        let m = demo_model();
        let hits = select(&m, "(chain A or chain B) and not name N").unwrap();
        // chain A: CA, CB per residue = 6; chain B water O = 1
        assert_eq!(hits.len(), 7);
    }

    #[test]
    fn backbone_sidechain() {
        let m = demo_model();
        let bb = select(&m, "backbone").unwrap();
        assert_eq!(bb.len(), 3 * 2); // N, CA per residue
        let sc = select(&m, "sidechain").unwrap();
        assert_eq!(sc.len(), 3); // CB per residue
    }

    #[test]
    fn water_and_element() {
        let m = demo_model();
        assert_eq!(select(&m, "water").unwrap().len(), 1);
        assert_eq!(select(&m, "element C").unwrap().len(), 6);
        assert_eq!(select(&m, "all").unwrap().len(), 10);
        assert_eq!(select(&m, "none").unwrap().len(), 0);
    }

    #[test]
    fn bad_selections_error() {
        assert!(Selection::parse("").is_err());
        assert!(Selection::parse("chain").is_err());
        assert!(Selection::parse("frobnicate A").is_err());
        assert!(Selection::parse("(chain A").is_err());
        assert!(Selection::parse("resi 10-5").is_err());
    }

    #[test]
    fn comma_argument_lists() {
        let m = demo_model();
        let hits = select(&m, "name N,CA").unwrap();
        assert_eq!(hits.len(), 3 * 2);
    }
}

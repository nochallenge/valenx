//! SMARTS substructure queries — a query parser plus VF2 matching.
//!
//! A [`SmartsPattern`] is a small query graph: each query atom is a
//! list of [`AtomExpr`] primitives ANDed together, each query bond a
//! list of [`BondExpr`] primitives. [`SmartsPattern::parse`] reads the
//! SMARTS string; [`SmartsPattern::find_matches`] runs the VF2
//! subgraph-isomorphism algorithm against a target [`Molecule`].
//!
//! Supported SMARTS primitives:
//!
//! - element by symbol (`C`, `Cl`, `O`) and aromatic / aliphatic forms
//!   (`c`, `n`); the wildcard `*` and `a` / `A` (any aromatic / any
//!   aliphatic);
//! - bracket primitives `[#6]` (atomic number), `[+]` / `[-]` /
//!   `[+2]` (charge), `[H2]` (total-H count), `[D3]` (degree),
//!   `[X4]` (connectivity), `[R]` / `[R0]` (ring membership),
//!   `[r6]` (smallest ring size), `[A]` / `[a]`;
//! - logical `,` (OR) inside a bracket; primitives juxtaposed are
//!   ANDed; a leading `!` negates a primitive;
//! - bonds `-`, `=`, `#`, `:`, `~` (any), `@` (ring bond), and the
//!   implicit "single or aromatic" bond;
//! - branches, ring-closure digits — the same structural grammar as
//!   SMILES;
//! - the `.` component separator — `[C:1].[O:2]` is a two-component
//!   query graph that matches two disconnected fragments, which is what
//!   a multi-reactant reaction template needs.
//!
//! **v1 simplifications:** recursive SMARTS `$(...)`, explicit
//! component-level grouping parentheses `(...)` and the `;`
//! low-precedence AND are not parsed. The supported subset covers the
//! functional-group and pharmacophore queries that drive descriptors,
//! standardisation and reaction matching in this crate.

use crate::error::{CheminfError, Result};
use crate::molecule::{BondOrder, Molecule};

/// One atom-matching primitive. A query atom matches iff *all* of its
/// primitives match (with `Or` handling alternation internally).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomExpr {
    /// Matches any atom (`*`).
    Any,
    /// Matches any aromatic atom (`a`).
    AnyAromatic,
    /// Matches any aliphatic (non-aromatic) atom (`A`).
    AnyAliphatic,
    /// Matches a specific element, with an aromatic-state requirement.
    /// `aromatic = Some(true)` → lowercase `c`; `Some(false)` →
    /// uppercase `C`; `None` → `[#6]`, either state.
    Element {
        /// Atomic number.
        z: u8,
        /// Required aromatic state, or `None` for don't-care.
        aromatic: Option<bool>,
    },
    /// Formal charge equals `i32`.
    Charge(i32),
    /// Total hydrogen count equals `u8`.
    TotalH(u8),
    /// Heavy-atom degree (explicit connections) equals `u8`.
    Degree(u8),
    /// Total connectivity (degree + implicit H) equals `u8` — `X`.
    Connectivity(u8),
    /// Ring-membership count equals `u8` (`R0` = acyclic, `R` ≥ 1).
    RingCount(Option<u8>),
    /// Member of a ring of exactly this size (`r6`).
    RingSize(u8),
    /// Logical OR of sub-expressions (the `,` operator).
    Or(Vec<AtomExpr>),
    /// Logical AND of sub-expressions (juxtaposition).
    And(Vec<AtomExpr>),
    /// Negation of a sub-expression (`!`).
    Not(Box<AtomExpr>),
}

/// One bond-matching primitive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondExpr {
    /// Single bond (`-`).
    Single,
    /// Double bond (`=`).
    Double,
    /// Triple bond (`#`).
    Triple,
    /// Aromatic bond (`:`).
    Aromatic,
    /// Any bond (`~`).
    Any,
    /// Any ring bond (`@`).
    Ring,
    /// The default SMARTS bond — single or aromatic.
    SingleOrAromatic,
}

/// A query atom: a primitive expression plus the atom-map number for
/// reaction transforms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryAtom {
    /// Conjunction of matching primitives.
    pub expr: AtomExpr,
    /// Atom-map number (`:1`), `0` if unmapped.
    pub map: u32,
}

/// A query bond between two query atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryBond {
    /// First query-atom index.
    pub a: usize,
    /// Second query-atom index.
    pub b: usize,
    /// The bond-matching expression.
    pub expr: BondExpr,
}

/// A compiled SMARTS query graph.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct SmartsPattern {
    /// Query atoms.
    pub atoms: Vec<QueryAtom>,
    /// Query bonds.
    pub bonds: Vec<QueryBond>,
}

impl SmartsPattern {
    /// Parse a SMARTS string into a query graph.
    pub fn parse(input: &str) -> Result<Self> {
        let s = input.trim();
        if s.is_empty() {
            return Err(CheminfError::parse("smarts", "empty pattern"));
        }
        SmartsParser::new(s).run()
    }

    /// Number of query atoms.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Find every match of this pattern in `target`. Each match is a
    /// `Vec` mapping query-atom index → target-atom index. All
    /// occurrences are returned (including symmetry-equivalent ones).
    pub fn find_matches(&self, target: &Molecule) -> Vec<Vec<usize>> {
        if self.atoms.is_empty() {
            return Vec::new();
        }
        let ctx = MatchContext::new(self, target);
        let mut results = Vec::new();
        let mut mapping = vec![usize::MAX; self.atoms.len()];
        let mut used = vec![false; target.atoms.len()];
        ctx.vf2(0, &mut mapping, &mut used, &mut results, false);
        results
    }

    /// Return the first match, or `None` — cheaper when only a yes/no
    /// answer or one hit is needed.
    pub fn find_first(&self, target: &Molecule) -> Option<Vec<usize>> {
        if self.atoms.is_empty() {
            return None;
        }
        let ctx = MatchContext::new(self, target);
        let mut results = Vec::new();
        let mut mapping = vec![usize::MAX; self.atoms.len()];
        let mut used = vec![false; target.atoms.len()];
        ctx.vf2(0, &mut mapping, &mut used, &mut results, true);
        results.into_iter().next()
    }

    /// `true` if the pattern occurs at least once in `target`.
    pub fn matches(&self, target: &Molecule) -> bool {
        self.find_first(target).is_some()
    }

    /// Count distinct matches whose *atom set* differs — collapses the
    /// symmetry-equivalent orderings of one occurrence.
    pub fn count_unique(&self, target: &Molecule) -> usize {
        let mut seen: Vec<Vec<usize>> = Vec::new();
        for m in self.find_matches(target) {
            let mut s = m.clone();
            s.sort_unstable();
            if !seen.contains(&s) {
                seen.push(s);
            }
        }
        seen.len()
    }
}

// --- VF2 matcher ------------------------------------------------------

/// Pre-computed context for one (pattern, target) matching run: the
/// pattern's atom order (DFS so each new query atom links to a
/// previously-matched one) and adjacency tables.
struct MatchContext<'a> {
    pattern: &'a SmartsPattern,
    target: &'a Molecule,
    /// Query-atom indices in DFS visitation order.
    order: Vec<usize>,
    /// SSSR ring analysis of the target, computed once.
    ring_info: crate::perceive::rings::RingInfo,
}

impl<'a> MatchContext<'a> {
    fn new(pattern: &'a SmartsPattern, target: &'a Molecule) -> Self {
        // DFS over the query graph from atom 0 so each step adds an
        // atom adjacent to one already placed → cheap pruning.
        let n = pattern.atoms.len();
        let mut order = Vec::with_capacity(n);
        let mut seen = vec![false; n];
        for start in 0..n {
            if seen[start] {
                continue;
            }
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(u) = stack.pop() {
                order.push(u);
                let mut nbrs: Vec<usize> = pattern
                    .bonds
                    .iter()
                    .filter_map(|qb| {
                        if qb.a == u {
                            Some(qb.b)
                        } else if qb.b == u {
                            Some(qb.a)
                        } else {
                            None
                        }
                    })
                    .collect();
                nbrs.sort_unstable();
                for v in nbrs {
                    if !seen[v] {
                        seen[v] = true;
                        stack.push(v);
                    }
                }
            }
        }
        MatchContext {
            pattern,
            target,
            order,
            ring_info: crate::perceive::rings::sssr(target),
        }
    }

    /// Recursive VF2 search. `depth` is the position in `self.order`;
    /// `mapping[q]` is the target atom assigned to query atom `q`.
    fn vf2(
        &self,
        depth: usize,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
        first_only: bool,
    ) {
        if first_only && !results.is_empty() {
            return;
        }
        if depth == self.order.len() {
            results.push(mapping.to_vec());
            return;
        }
        let q = self.order[depth];
        for t in 0..self.target.atoms.len() {
            if used[t] {
                continue;
            }
            if !self.atom_matches(q, t) {
                continue;
            }
            if !self.bonds_consistent(q, t, mapping) {
                continue;
            }
            mapping[q] = t;
            used[t] = true;
            self.vf2(depth + 1, mapping, used, results, first_only);
            mapping[q] = usize::MAX;
            used[t] = false;
            if first_only && !results.is_empty() {
                return;
            }
        }
    }

    /// Every already-placed query neighbour of `q` must be joined to
    /// `t` by a target bond satisfying the query bond expression.
    fn bonds_consistent(&self, q: usize, t: usize, mapping: &[usize]) -> bool {
        for qb in &self.pattern.bonds {
            let other_q = if qb.a == q {
                qb.b
            } else if qb.b == q {
                qb.a
            } else {
                continue;
            };
            let mapped = mapping[other_q];
            if mapped == usize::MAX {
                continue; // neighbour not placed yet
            }
            match self.target.bond_between(t, mapped) {
                None => return false,
                Some(bi) => {
                    let in_ring = self.ring_info.bond_in_ring(bi);
                    if !bond_expr_matches(qb.expr, &self.target.bonds[bi], in_ring) {
                        return false;
                    }
                }
            }
        }
        true
    }

    fn atom_matches(&self, q: usize, t: usize) -> bool {
        atom_expr_matches(
            &self.pattern.atoms[q].expr,
            self.target,
            t,
            &self.ring_info,
        )
    }
}

/// Evaluate an [`AtomExpr`] against target atom `t`, using the
/// pre-computed `rings` for ring-membership primitives.
fn atom_expr_matches(
    expr: &AtomExpr,
    mol: &Molecule,
    t: usize,
    rings: &crate::perceive::rings::RingInfo,
) -> bool {
    let a = &mol.atoms[t];
    match expr {
        AtomExpr::Any => true,
        AtomExpr::AnyAromatic => a.aromatic,
        AtomExpr::AnyAliphatic => !a.aromatic && !a.is_dummy(),
        AtomExpr::Element { z, aromatic } => {
            a.atomic_number == *z
                && match aromatic {
                    Some(want) => a.aromatic == *want,
                    None => true,
                }
        }
        AtomExpr::Charge(c) => i32::from(a.formal_charge) == *c,
        AtomExpr::TotalH(h) => a.total_h() == *h,
        AtomExpr::Degree(d) => mol.degree(t) as u8 == *d,
        AtomExpr::Connectivity(x) => (mol.degree(t) + usize::from(a.total_h())) as u8 == *x,
        AtomExpr::RingCount(want) => {
            let count = rings.atom_ring_count.get(t).copied().unwrap_or(0) as u8;
            match want {
                Some(n) => count == *n,
                None => count >= 1,
            }
        }
        AtomExpr::RingSize(size) => rings
            .rings
            .iter()
            .any(|r| r.size() as u8 == *size && r.contains_atom(t)),
        AtomExpr::Or(subs) => subs.iter().any(|e| atom_expr_matches(e, mol, t, rings)),
        AtomExpr::And(subs) => subs.iter().all(|e| atom_expr_matches(e, mol, t, rings)),
        AtomExpr::Not(inner) => !atom_expr_matches(inner, mol, t, rings),
    }
}

/// Evaluate a [`BondExpr`] against a concrete target bond; `in_ring`
/// says whether that bond is in an SSSR ring (for the `@` primitive).
fn bond_expr_matches(expr: BondExpr, bond: &crate::molecule::Bond, in_ring: bool) -> bool {
    match expr {
        BondExpr::Any => true,
        BondExpr::Single => bond.order == BondOrder::Single,
        BondExpr::Double => bond.order == BondOrder::Double,
        BondExpr::Triple => bond.order == BondOrder::Triple,
        BondExpr::Aromatic => bond.order == BondOrder::Aromatic || bond.aromatic,
        BondExpr::Ring => in_ring,
        BondExpr::SingleOrAromatic => {
            matches!(bond.order, BondOrder::Single | BondOrder::Aromatic) || bond.aromatic
        }
    }
}

// --- SMARTS parser ----------------------------------------------------

struct SmartsParser<'a> {
    chars: Vec<char>,
    pos: usize,
    src: &'a str,
    pat: SmartsPattern,
    branch: Vec<usize>,
    prev: Option<usize>,
    pending_bond: Option<BondExpr>,
    rings: std::collections::HashMap<u8, (usize, Option<BondExpr>)>,
}

impl<'a> SmartsParser<'a> {
    fn new(src: &'a str) -> Self {
        SmartsParser {
            chars: src.chars().collect(),
            pos: 0,
            src,
            pat: SmartsPattern::default(),
            branch: Vec::new(),
            prev: None,
            pending_bond: None,
            rings: std::collections::HashMap::new(),
        }
    }

    fn err(&self, d: impl Into<String>) -> CheminfError {
        CheminfError::parse(
            "smarts",
            format!("{} (at {} of `{}`)", d.into(), self.pos, self.src),
        )
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn run(mut self) -> Result<SmartsPattern> {
        while let Some(c) = self.peek() {
            match c {
                '(' => {
                    let anchor = self.prev.ok_or_else(|| self.err("'(' with no atom"))?;
                    self.branch.push(anchor);
                    self.pos += 1;
                }
                ')' => {
                    self.prev =
                        Some(self.branch.pop().ok_or_else(|| self.err("unbalanced ')'"))?);
                    self.pos += 1;
                }
                '-' | '=' | '#' | ':' | '~' | '@' => {
                    self.read_bond(c)?;
                }
                '[' => {
                    let atom = self.read_bracket()?;
                    self.connect(atom);
                }
                '0'..='9' | '%' => self.read_ring()?,
                '.' => {
                    // Component separator: the atom after a `.` starts a
                    // new disconnected query component, with no bond to
                    // the previous atom. The VF2 matcher already handles
                    // a multi-component query graph (it seeds a DFS from
                    // every unvisited atom), so a `.`-joined pattern such
                    // as `[C:1].[O:2]` matches two distinct fragments —
                    // exactly what a multi-reactant reaction template
                    // needs.
                    if self.prev.is_none() {
                        return Err(self.err("'.' with no preceding atom"));
                    }
                    if self.pending_bond.is_some() {
                        return Err(self.err("bond symbol immediately before '.'"));
                    }
                    if !self.branch.is_empty() {
                        return Err(self.err("'.' inside a branch"));
                    }
                    self.prev = None;
                    self.pos += 1;
                }
                'A' | 'a' | '*' => {
                    let expr = match c {
                        'A' => AtomExpr::AnyAliphatic,
                        'a' => AtomExpr::AnyAromatic,
                        '*' => AtomExpr::Any,
                        _ => unreachable!(),
                    };
                    self.pos += 1;
                    let atom = self.push_atom(expr, 0);
                    self.connect(atom);
                }
                c if c.is_ascii_alphabetic() => {
                    let atom = self.read_organic()?;
                    self.connect(atom);
                }
                ' ' | '\t' | '\n' | '\r' => break,
                other => return Err(self.err(format!("unexpected '{other}'"))),
            }
        }
        if !self.branch.is_empty() {
            return Err(self.err("unclosed branch"));
        }
        if let Some((&d, _)) = self.rings.iter().next() {
            return Err(self.err(format!("unclosed ring digit {d}")));
        }
        if self.pat.atoms.is_empty() {
            return Err(self.err("no query atoms"));
        }
        Ok(self.pat)
    }

    fn push_atom(&mut self, expr: AtomExpr, map: u32) -> usize {
        self.pat.atoms.push(QueryAtom { expr, map });
        self.pat.atoms.len() - 1
    }

    fn connect(&mut self, atom: usize) {
        if let Some(prev) = self.prev {
            let expr = self.pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic);
            self.pat.bonds.push(QueryBond {
                a: prev,
                b: atom,
                expr,
            });
        } else {
            self.pending_bond = None;
        }
        self.prev = Some(atom);
    }

    fn read_bond(&mut self, c: char) -> Result<()> {
        let expr = match c {
            '-' => BondExpr::Single,
            '=' => BondExpr::Double,
            '#' => BondExpr::Triple,
            ':' => BondExpr::Aromatic,
            '~' => BondExpr::Any,
            '@' => BondExpr::Ring,
            _ => unreachable!(),
        };
        if self.pending_bond.is_some() {
            return Err(self.err("two bond symbols"));
        }
        self.pending_bond = Some(expr);
        self.pos += 1;
        Ok(())
    }

    fn read_organic(&mut self) -> Result<usize> {
        let c = self.chars[self.pos];
        // two-letter elements first
        let two: String = self.chars[self.pos..].iter().take(2).collect();
        if c == 'C' && two == "Cl" {
            self.pos += 2;
            let atom = self.push_atom(
                AtomExpr::Element {
                    z: 17,
                    aromatic: Some(false),
                },
                0,
            );
            return Ok(atom);
        }
        if c == 'B' && two == "Br" {
            self.pos += 2;
            let atom = self.push_atom(
                AtomExpr::Element {
                    z: 35,
                    aromatic: Some(false),
                },
                0,
            );
            return Ok(atom);
        }
        let aromatic = c.is_ascii_lowercase();
        let sym = c.to_ascii_uppercase().to_string();
        let z = crate::element::by_symbol(&sym)
            .map(|e| e.number)
            .ok_or_else(|| self.err(format!("unknown organic atom `{c}`")))?;
        self.pos += 1;
        let atom = self.push_atom(
            AtomExpr::Element {
                z,
                aromatic: Some(aromatic),
            },
            0,
        );
        Ok(atom)
    }

    /// Parse a `[ ... ]` bracket query atom — a sequence of primitives
    /// combined with `,` (OR) and juxtaposition (AND), each optionally
    /// `!`-negated.
    fn read_bracket(&mut self) -> Result<usize> {
        debug_assert_eq!(self.chars[self.pos], '[');
        self.pos += 1;
        let mut or_groups: Vec<AtomExpr> = Vec::new();
        let mut and_group: Vec<AtomExpr> = Vec::new();
        let mut map = 0u32;

        loop {
            match self.peek() {
                None => return Err(self.err("unterminated bracket atom")),
                Some(']') => {
                    self.pos += 1;
                    break;
                }
                Some(',') => {
                    self.pos += 1;
                    or_groups.push(fold_and(std::mem::take(&mut and_group)));
                }
                Some(':') => {
                    self.pos += 1;
                    let mut n = String::new();
                    while let Some(d @ '0'..='9') = self.peek() {
                        n.push(d);
                        self.pos += 1;
                    }
                    map = n.parse().unwrap_or(0);
                }
                Some(_) => {
                    let prim = self.read_primitive()?;
                    and_group.push(prim);
                }
            }
        }
        or_groups.push(fold_and(and_group));
        let expr = if or_groups.len() == 1 {
            or_groups.pop().unwrap()
        } else {
            AtomExpr::Or(or_groups)
        };
        Ok(self.push_atom(expr, map))
    }

    /// One bracket primitive, possibly `!`-negated.
    fn read_primitive(&mut self) -> Result<AtomExpr> {
        let mut negate = false;
        while self.peek() == Some('!') {
            negate = !negate;
            self.pos += 1;
        }
        let base = self.read_primitive_inner()?;
        Ok(if negate {
            AtomExpr::Not(Box::new(base))
        } else {
            base
        })
    }

    fn read_primitive_inner(&mut self) -> Result<AtomExpr> {
        let c = self.peek().ok_or_else(|| self.err("primitive expected"))?;
        match c {
            '*' => {
                self.pos += 1;
                Ok(AtomExpr::Any)
            }
            'a' => {
                self.pos += 1;
                Ok(AtomExpr::AnyAromatic)
            }
            'A' => {
                self.pos += 1;
                Ok(AtomExpr::AnyAliphatic)
            }
            '#' => {
                self.pos += 1;
                let z = self.read_uint()? as u8;
                Ok(AtomExpr::Element { z, aromatic: None })
            }
            '+' => {
                self.pos += 1;
                let mut n = 1i32;
                while self.peek() == Some('+') {
                    n += 1;
                    self.pos += 1;
                }
                if n == 1 {
                    if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                        n = d.to_digit(10).unwrap() as i32;
                        self.pos += 1;
                    }
                }
                Ok(AtomExpr::Charge(n))
            }
            '-' => {
                self.pos += 1;
                let mut n = 1i32;
                while self.peek() == Some('-') {
                    n += 1;
                    self.pos += 1;
                }
                if n == 1 {
                    if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                        n = d.to_digit(10).unwrap() as i32;
                        self.pos += 1;
                    }
                }
                Ok(AtomExpr::Charge(-n))
            }
            'H' => {
                self.pos += 1;
                let n = self.read_uint_default(1)?;
                Ok(AtomExpr::TotalH(n as u8))
            }
            'D' => {
                self.pos += 1;
                let n = self.read_uint_default(1)?;
                Ok(AtomExpr::Degree(n as u8))
            }
            'X' => {
                self.pos += 1;
                let n = self.read_uint_default(1)?;
                Ok(AtomExpr::Connectivity(n as u8))
            }
            'R' => {
                self.pos += 1;
                if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                    self.pos += 1;
                    Ok(AtomExpr::RingCount(Some(d.to_digit(10).unwrap() as u8)))
                } else {
                    Ok(AtomExpr::RingCount(None))
                }
            }
            'r' => {
                self.pos += 1;
                let n = self.read_uint_default(3)?;
                Ok(AtomExpr::RingSize(n as u8))
            }
            c if c.is_ascii_uppercase() => {
                // an element symbol, maybe two letters
                let start = self.pos;
                self.pos += 1;
                if let Some(c2) = self.peek() {
                    if c2.is_ascii_lowercase() {
                        let two: String = self.chars[start..start + 2].iter().collect();
                        if crate::element::by_symbol(&two).is_some() {
                            self.pos += 1;
                        }
                    }
                }
                let sym: String = self.chars[start..self.pos].iter().collect();
                let z = crate::element::by_symbol(&sym)
                    .map(|e| e.number)
                    .ok_or_else(|| self.err(format!("unknown element `{sym}`")))?;
                Ok(AtomExpr::Element {
                    z,
                    aromatic: Some(false),
                })
            }
            c if c.is_ascii_lowercase() => {
                // aromatic element
                let sym = c.to_ascii_uppercase().to_string();
                let z = crate::element::by_symbol(&sym)
                    .map(|e| e.number)
                    .ok_or_else(|| self.err(format!("unknown aromatic element `{c}`")))?;
                self.pos += 1;
                Ok(AtomExpr::Element {
                    z,
                    aromatic: Some(true),
                })
            }
            other => Err(self.err(format!("unsupported SMARTS primitive '{other}'"))),
        }
    }

    fn read_uint(&mut self) -> Result<u32> {
        let mut n = String::new();
        while let Some(d @ '0'..='9') = self.peek() {
            n.push(d);
            self.pos += 1;
        }
        n.parse()
            .map_err(|_| self.err("expected a non-negative integer"))
    }

    fn read_uint_default(&mut self, default: u32) -> Result<u32> {
        if self.peek().map(|c| c.is_ascii_digit()).unwrap_or(false) {
            self.read_uint()
        } else {
            Ok(default)
        }
    }

    fn read_ring(&mut self) -> Result<()> {
        let label: u8 = if self.peek() == Some('%') {
            self.pos += 1;
            let d1 = self
                .peek()
                .filter(|c| c.is_ascii_digit())
                .ok_or_else(|| self.err("'%' needs two digits"))?;
            self.pos += 1;
            let d2 = self
                .peek()
                .filter(|c| c.is_ascii_digit())
                .ok_or_else(|| self.err("'%' needs two digits"))?;
            self.pos += 1;
            (d1.to_digit(10).unwrap() * 10 + d2.to_digit(10).unwrap()) as u8
        } else {
            let d = self.chars[self.pos];
            self.pos += 1;
            d.to_digit(10).unwrap() as u8
        };
        let atom = self
            .prev
            .ok_or_else(|| self.err("ring digit with no atom"))?;
        let pending = self.pending_bond.take();
        if let Some((open_atom, open_bond)) = self.rings.remove(&label) {
            let expr = pending.or(open_bond).unwrap_or(BondExpr::SingleOrAromatic);
            self.pat.bonds.push(QueryBond {
                a: open_atom,
                b: atom,
                expr,
            });
        } else {
            self.rings.insert(label, (atom, pending));
        }
        Ok(())
    }
}

/// Collapse a list of AND-ed primitives into a single [`AtomExpr`].
fn fold_and(mut group: Vec<AtomExpr>) -> AtomExpr {
    match group.len() {
        0 => AtomExpr::Any,
        1 => group.pop().unwrap(),
        _ => AtomExpr::And(group),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn parse_simple_pattern() {
        let p = SmartsPattern::parse("CCO").unwrap();
        assert_eq!(p.atom_count(), 3);
        assert_eq!(p.bonds.len(), 2);
    }

    #[test]
    fn element_match() {
        let mol = mol_from_smiles("CCO").unwrap();
        let p = SmartsPattern::parse("O").unwrap();
        assert!(p.matches(&mol));
        assert_eq!(p.count_unique(&mol), 1);

        let p = SmartsPattern::parse("N").unwrap();
        assert!(!p.matches(&mol));
    }

    #[test]
    fn carbonyl_pattern() {
        let acetone = mol_from_smiles("CC(=O)C").unwrap();
        let p = SmartsPattern::parse("C=O").unwrap();
        assert!(p.matches(&acetone));
        assert_eq!(p.count_unique(&acetone), 1);

        let ethanol = mol_from_smiles("CCO").unwrap();
        assert!(!p.matches(&ethanol));
    }

    #[test]
    fn aromatic_query() {
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let p = SmartsPattern::parse("c").unwrap();
        // 6 aromatic carbons
        assert_eq!(p.count_unique(&benzene), 6);

        let cyclohexane = mol_from_smiles("C1CCCCC1").unwrap();
        assert_eq!(p.count_unique(&cyclohexane), 0);
    }

    #[test]
    fn bracket_primitives() {
        let mol = mol_from_smiles("CC(=O)O").unwrap(); // acetic acid
        // carboxyl carbon: a carbon with degree 3
        let p = SmartsPattern::parse("[#6][D3]").unwrap();
        assert!(p.matches(&mol));

        // charged atom query on an ammonium
        let amm = mol_from_smiles("[NH4+]").unwrap();
        let p = SmartsPattern::parse("[N+]").unwrap();
        assert!(p.matches(&amm));
    }

    #[test]
    fn ring_membership_query() {
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let p = SmartsPattern::parse("[R]").unwrap();
        assert_eq!(p.count_unique(&benzene), 6);

        let p = SmartsPattern::parse("[r6]").unwrap();
        assert_eq!(p.count_unique(&benzene), 6);

        let hexane = mol_from_smiles("CCCCCC").unwrap();
        let p = SmartsPattern::parse("[R]").unwrap();
        assert_eq!(p.count_unique(&hexane), 0);
    }

    #[test]
    fn or_and_not() {
        let mol = mol_from_smiles("CCN").unwrap();
        // nitrogen OR oxygen
        let p = SmartsPattern::parse("[N,O]").unwrap();
        assert!(p.matches(&mol));
        // not carbon
        let p = SmartsPattern::parse("[!#6]").unwrap();
        assert_eq!(p.count_unique(&mol), 1);
    }

    #[test]
    fn rejects_bad_smarts() {
        assert!(SmartsPattern::parse("").is_err());
        assert!(SmartsPattern::parse("C(").is_err());
        assert!(SmartsPattern::parse("[Zz]").is_err());
    }
}

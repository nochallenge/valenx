//! SMILES parser — string → [`Molecule`].
//!
//! Covers the structural core of the Daylight / OpenSMILES grammar:
//!
//! - the organic-subset atoms `B C N O P S F Cl Br I` and the aromatic
//!   lowercase forms `b c n o p s`;
//! - bracket atoms `[<isotope><symbol><chirality><Hn><charge><:map>]`
//!   for everything else (any element, explicit H count, formal
//!   charge, isotope, atom-map number, `@` / `@@` chirality);
//! - bonds `- = # : / \` and the implicit single / aromatic bond;
//! - branches `( … )` to arbitrary depth;
//! - ring-closure digits `1`–`9` and the two-digit `%nn` form, with
//!   an optional bond symbol on either the opening or closing digit;
//! - the dot `.` fragment separator.
//!
//! After the graph is built the parser fills valence-based implicit
//! hydrogens for organic-subset atoms (bracket atoms keep exactly the
//! hydrogens written) and runs a light aromatic-bond fix-up so a
//! lowercase ring gets [`BondOrder::Aromatic`] ring bonds.
//!
//! **v1 simplifications:** extended tetrahedral / allene / square-planar
//! / trigonal-bipyramidal chirality classes (`@TH`, `@AL`, `@SP`,
//! `@TB`, `@OH`) are not parsed — only plain `@` / `@@`; the wildcard
//! `*` is accepted as a dummy atom. Reaction SMILES (`>`) is handled by
//! [`crate::reaction`], not here.

use crate::element::{self, AROMATIC_ELEMENTS, ORGANIC_SUBSET};
use crate::error::{CheminfError, Result};
use crate::molecule::{Atom, Bond, BondOrder, BondStereo, Chirality, Molecule};

/// Parse a SMILES string into a [`Molecule`].
///
/// Implicit hydrogens are filled for organic-subset atoms; bracket
/// atoms keep exactly the hydrogens written between the brackets.
pub fn parse_smiles(input: &str) -> Result<Molecule> {
    let s = input.trim();
    if s.is_empty() {
        return Err(CheminfError::parse("smiles", "empty input"));
    }
    let mut p = Parser::new(s);
    p.run()?;
    let bracket_atoms = std::mem::take(&mut p.bracket_atoms);
    let mut mol = p.finish()?;
    fill_implicit_hydrogens(&mut mol);
    // A bracket atom's hydrogen count is exactly what the brackets
    // carried (the `explicit_h`); it gets NO valence-based implicit
    // hydrogens. `fill_implicit_hydrogens` cannot tell an organic
    // bracket atom written with no `Hn` (e.g. `[O-]`, 0 H) from an
    // unbracketed one, so any implicit H it added to a bracket atom is
    // cleared back to zero here.
    for &i in &bracket_atoms {
        if i < mol.atoms.len() {
            mol.atoms[i].implicit_h = 0;
        }
    }
    fixup_aromatic_bonds(&mut mol);
    mol.validate()?;
    Ok(mol)
}

/// State carried while a ring-closure digit is open: the atom that
/// opened it and any bond symbol attached to the opening digit.
#[derive(Clone)]
struct RingOpen {
    atom: usize,
    order: Option<BondOrder>,
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    src: &'a str,
    mol: Molecule,
    /// Atom-index stack — `top` is "where the next atom bonds to".
    branch: Vec<usize>,
    /// Pending bond symbol read before the next atom.
    pending_bond: Option<BondOrder>,
    pending_stereo: BondStereo,
    /// Open ring closures keyed by digit (1..=99).
    rings: std::collections::HashMap<u8, RingOpen>,
    /// Index of the most recently added atom, or `None` at start /
    /// after a dot.
    prev: Option<usize>,
    /// Indices of atoms written in `[ … ]` brackets. A bracket atom's
    /// hydrogen count is *exactly* what was written (0 if no `Hn`), so
    /// these atoms must be excluded from valence-based implicit-H
    /// filling — `[O-]` is the bare oxide ion, not hydroxide `[OH-]`.
    bracket_atoms: Vec<usize>,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Self {
        Parser {
            chars: src.chars().collect(),
            pos: 0,
            src,
            mol: Molecule::new(),
            branch: Vec::new(),
            pending_bond: None,
            pending_stereo: BondStereo::None,
            rings: std::collections::HashMap::new(),
            prev: None,
            bracket_atoms: Vec::new(),
        }
    }

    fn err(&self, detail: impl Into<String>) -> CheminfError {
        let d = detail.into();
        CheminfError::parse(
            "smiles",
            format!("{d} (at position {} of `{}`)", self.pos, self.src),
        )
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn run(&mut self) -> Result<()> {
        while let Some(c) = self.peek() {
            match c {
                '(' => {
                    let anchor = self
                        .prev
                        .ok_or_else(|| self.err("'(' with no preceding atom"))?;
                    self.branch.push(anchor);
                    self.pos += 1;
                }
                ')' => {
                    self.prev = self
                        .branch
                        .pop()
                        .ok_or_else(|| self.err("unbalanced ')'"))?
                        .into();
                    self.pos += 1;
                }
                '.' => {
                    self.prev = None;
                    self.pending_bond = None;
                    self.pending_stereo = BondStereo::None;
                    self.pos += 1;
                }
                '-' | '=' | '#' | '$' | ':' | '/' | '\\' => {
                    self.read_bond_symbol(c)?;
                }
                '[' => {
                    let atom = self.read_bracket_atom()?;
                    self.bracket_atoms.push(atom);
                    self.connect(atom)?;
                }
                '0'..='9' | '%' => {
                    self.read_ring_closure()?;
                }
                'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' => {
                    let atom = self.read_organic_atom()?;
                    self.connect(atom)?;
                }
                'b' | 'c' | 'n' | 'o' | 'p' | 's' => {
                    let atom = self.read_organic_aromatic()?;
                    self.connect(atom)?;
                }
                '*' => {
                    let atom = self.mol.add_atom(Atom::new(0));
                    self.pos += 1;
                    self.connect(atom)?;
                }
                ' ' | '\t' | '\n' | '\r' => {
                    // Trailing whitespace / a name field ends the SMILES.
                    break;
                }
                other => return Err(self.err(format!("unexpected character '{other}'"))),
            }
        }
        Ok(())
    }

    fn read_bond_symbol(&mut self, c: char) -> Result<()> {
        if self.pending_bond.is_some() || self.pending_stereo != BondStereo::None {
            return Err(self.err("two bond symbols in a row"));
        }
        match c {
            '-' => self.pending_bond = Some(BondOrder::Single),
            '=' => self.pending_bond = Some(BondOrder::Double),
            '#' => self.pending_bond = Some(BondOrder::Triple),
            '$' => self.pending_bond = Some(BondOrder::Quadruple),
            ':' => self.pending_bond = Some(BondOrder::Aromatic),
            '/' => {
                self.pending_bond = Some(BondOrder::Single);
                self.pending_stereo = BondStereo::Up;
            }
            '\\' => {
                self.pending_bond = Some(BondOrder::Single);
                self.pending_stereo = BondStereo::Down;
            }
            _ => unreachable!(),
        }
        self.pos += 1;
        Ok(())
    }

    /// Bond an atom that is *not* a ring digit to the current `prev`.
    fn connect(&mut self, atom: usize) -> Result<()> {
        if let Some(prev) = self.prev {
            let order = self.pending_bond.take().unwrap_or_else(|| {
                if self.mol.atoms[prev].aromatic && self.mol.atoms[atom].aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            });
            self.mol.add_bond(Bond {
                a: prev,
                b: atom,
                order,
                aromatic: order == BondOrder::Aromatic,
                stereo: std::mem::replace(&mut self.pending_stereo, BondStereo::None),
            });
        } else {
            self.pending_bond = None;
            self.pending_stereo = BondStereo::None;
        }
        self.prev = Some(atom);
        Ok(())
    }

    fn read_organic_atom(&mut self) -> Result<usize> {
        // Two-letter halogens take precedence over one-letter atoms.
        let c = self.chars[self.pos];
        let two: String = self.chars[self.pos..].iter().take(2).collect();
        let (sym, len) = match (c, two.as_str()) {
            ('C', "Cl") => ("Cl", 2),
            ('B', "Br") => ("Br", 2),
            _ => (
                match c {
                    'B' => "B",
                    'C' => "C",
                    'N' => "N",
                    'O' => "O",
                    'P' => "P",
                    'S' => "S",
                    'F' => "F",
                    'I' => "I",
                    _ => return Err(self.err("not an organic-subset atom")),
                },
                1,
            ),
        };
        let z = element::by_symbol(sym).map(|e| e.number).unwrap();
        self.pos += len;
        Ok(self.mol.add_atom(Atom::new(z)))
    }

    fn read_organic_aromatic(&mut self) -> Result<usize> {
        let c = self.chars[self.pos];
        let sym = match c {
            'b' => "B",
            'c' => "C",
            'n' => "N",
            'o' => "O",
            'p' => "P",
            's' => "S",
            _ => return Err(self.err("not an aromatic organic atom")),
        };
        let z = element::by_symbol(sym).map(|e| e.number).unwrap();
        if !AROMATIC_ELEMENTS.contains(&z) {
            return Err(self.err("element cannot be aromatic"));
        }
        self.pos += 1;
        let mut atom = Atom::new(z);
        atom.aromatic = true;
        Ok(self.mol.add_atom(atom))
    }

    /// Parse a `[ … ]` bracket atom.
    fn read_bracket_atom(&mut self) -> Result<usize> {
        debug_assert_eq!(self.chars[self.pos], '[');
        self.pos += 1;
        let mut atom = Atom::new(0);

        // optional isotope
        let mut iso = String::new();
        while let Some(d @ '0'..='9') = self.peek() {
            iso.push(d);
            self.pos += 1;
        }
        if !iso.is_empty() {
            atom.isotope = iso.parse::<u16>().ok();
        }

        // element symbol — uppercase (+ optional lowercase) or aromatic
        let sym_start = self.pos;
        let aromatic;
        let c = self
            .peek()
            .ok_or_else(|| self.err("unterminated bracket atom"))?;
        if c == '*' {
            self.pos += 1;
            atom.atomic_number = 0;
            aromatic = false;
        } else if c.is_ascii_uppercase() {
            self.pos += 1;
            if let Some(c2) = self.peek() {
                if c2.is_ascii_lowercase() {
                    let two: String = self.chars[sym_start..sym_start + 2].iter().collect();
                    if element::by_symbol(&two).is_some() {
                        self.pos += 1;
                    }
                }
            }
            let sym: String = self.chars[sym_start..self.pos].iter().collect();
            atom.atomic_number = element::by_symbol(&sym)
                .map(|e| e.number)
                .ok_or_else(|| self.err(format!("unknown element `{sym}`")))?;
            aromatic = false;
        } else if c.is_ascii_lowercase() {
            // aromatic element: c, n, o, p, s, b, se, as, te
            self.pos += 1;
            if let Some(c2) = self.peek() {
                if c2.is_ascii_lowercase() {
                    let two: String = [c.to_ascii_uppercase(), c2].iter().collect::<String>();
                    if element::by_symbol(&two).is_some() {
                        self.pos += 1;
                    }
                }
            }
            let sym: String = self.chars[sym_start..self.pos].iter().collect::<String>();
            let cap: String = sym
                .char_indices()
                .map(|(i, ch)| if i == 0 { ch.to_ascii_uppercase() } else { ch })
                .collect();
            atom.atomic_number = element::by_symbol(&cap)
                .map(|e| e.number)
                .ok_or_else(|| self.err(format!("unknown aromatic element `{sym}`")))?;
            aromatic = true;
        } else {
            return Err(self.err("expected an element symbol in bracket atom"));
        }
        atom.aromatic = aromatic;

        // optional chirality — plain `@` / `@@` only; the extended
        // tetrahedral / allene / square-planar / TB / OH classes
        // (`@TH1`, `@AL2`, `@SP1`, `@TB5`, `@OH3`) are not parsed in v1
        if self.peek() == Some('@') {
            self.pos += 1;
            if self.peek() == Some('@') {
                self.pos += 1;
                atom.chirality = Chirality::Cw;
            } else {
                atom.chirality = Chirality::Ccw;
            }
        }

        // optional explicit hydrogen count
        if self.peek() == Some('H') {
            self.pos += 1;
            let mut n = String::new();
            while let Some(d @ '0'..='9') = self.peek() {
                n.push(d);
                self.pos += 1;
            }
            atom.explicit_h = if n.is_empty() {
                1
            } else {
                n.parse::<u8>().unwrap_or(1)
            };
        }

        // optional formal charge
        match self.peek() {
            Some('+') => {
                self.pos += 1;
                atom.formal_charge = self.read_charge_magnitude(1)?;
            }
            Some('-') => {
                self.pos += 1;
                atom.formal_charge = self.read_charge_magnitude(-1)?;
            }
            _ => {}
        }

        // optional atom-map number
        if self.peek() == Some(':') {
            self.pos += 1;
            let mut n = String::new();
            while let Some(d @ '0'..='9') = self.peek() {
                n.push(d);
                self.pos += 1;
            }
            atom.map_number = n.parse::<u32>().unwrap_or(0);
        }

        if self.peek() != Some(']') {
            return Err(self.err("unterminated bracket atom — expected ']'"));
        }
        self.pos += 1;
        Ok(self.mol.add_atom(atom))
    }

    /// After a `+` or `-`, read either a run of the same sign (`++`),
    /// or an explicit magnitude (`+2`), or nothing (lone `+`).
    fn read_charge_magnitude(&mut self, sign: i8) -> Result<i8> {
        // run of identical signs
        let run_char = if sign > 0 { '+' } else { '-' };
        let mut count: i32 = 1;
        while self.peek() == Some(run_char) {
            count += 1;
            self.pos += 1;
        }
        if count == 1 {
            // maybe an explicit digit
            let mut n = String::new();
            while let Some(d @ '0'..='9') = self.peek() {
                n.push(d);
                self.pos += 1;
            }
            if !n.is_empty() {
                count = n
                    .parse::<i32>()
                    .map_err(|_| self.err("bad charge magnitude"))?;
            }
        }
        let signed = i32::from(sign) * count;
        i8::try_from(signed).map_err(|_| self.err("charge out of range"))
    }

    /// Parse a ring-closure digit or `%nn` two-digit closure.
    fn read_ring_closure(&mut self) -> Result<()> {
        let label: u8 = if self.peek() == Some('%') {
            self.pos += 1;
            let d1 = self
                .peek()
                .filter(|c| c.is_ascii_digit())
                .ok_or_else(|| self.err("'%' must be followed by two digits"))?;
            self.pos += 1;
            let d2 = self
                .peek()
                .filter(|c| c.is_ascii_digit())
                .ok_or_else(|| self.err("'%' must be followed by two digits"))?;
            self.pos += 1;
            (d1.to_digit(10).unwrap() * 10 + d2.to_digit(10).unwrap()) as u8
        } else {
            let d = self.chars[self.pos];
            self.pos += 1;
            d.to_digit(10).unwrap() as u8
        };

        let atom = self
            .prev
            .ok_or_else(|| self.err("ring-closure digit with no preceding atom"))?;
        let pending = self.pending_bond.take();
        self.pending_stereo = BondStereo::None;

        if let Some(open) = self.rings.remove(&label) {
            // Closing an existing ring — wire the bond.
            if open.atom == atom {
                return Err(self.err("ring closure bonds an atom to itself"));
            }
            let order = pending.or(open.order).unwrap_or_else(|| {
                if self.mol.atoms[open.atom].aromatic && self.mol.atoms[atom].aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            });
            if self.mol.bond_between(open.atom, atom).is_none() {
                self.mol.add_bond(Bond {
                    a: open.atom,
                    b: atom,
                    order,
                    aromatic: order == BondOrder::Aromatic,
                    stereo: BondStereo::None,
                });
            }
        } else {
            self.rings.insert(
                label,
                RingOpen {
                    atom,
                    order: pending,
                },
            );
        }
        Ok(())
    }

    fn finish(self) -> Result<Molecule> {
        if !self.branch.is_empty() {
            return Err(CheminfError::parse(
                "smiles",
                format!("{} unclosed '(' branch(es)", self.branch.len()),
            ));
        }
        if let Some((&d, _)) = self.rings.iter().next() {
            return Err(CheminfError::parse(
                "smiles",
                format!("unclosed ring-bond digit {d}"),
            ));
        }
        if self.mol.is_empty() {
            return Err(CheminfError::parse("smiles", "no atoms parsed"));
        }
        Ok(self.mol)
    }
}

/// Fill valence-based implicit hydrogens for organic-subset atoms.
///
/// Bracket atoms (any atom whose hydrogens were written explicitly, or
/// any non-organic-subset element) keep exactly what they declared.
/// Organic-subset atoms get `H = standard_valence - bond_order_sum`,
/// clamped at zero, with the formal charge shifting the target valence
/// (e.g. `[NH4+]` style behaviour for an unbracketed `N` is *not*
/// triggered — only bracketed atoms can carry an explicit count, so a
/// charged organic-subset atom is rare and handled by charge offset).
pub fn fill_implicit_hydrogens(mol: &mut Molecule) {
    // Mark which atoms came from a bracket (explicit H already set, or
    // an element outside the organic subset). We approximate: an atom
    // is "organic-subset implicit" iff its element is in ORGANIC_SUBSET
    // and it carries no explicit H. Bracket organic atoms with H0 will
    // therefore also be filled — acceptable for v1 since `[CH0]` is
    // unusual; callers needing exact control build the molecule
    // directly.
    for i in 0..mol.atoms.len() {
        let a = &mol.atoms[i];
        if a.is_dummy() || a.is_hydrogen() {
            continue;
        }
        if !ORGANIC_SUBSET.contains(&a.atomic_number) {
            continue;
        }
        if a.explicit_h > 0 {
            continue;
        }
        let z = a.atomic_number;
        let charge = a.formal_charge;
        let aromatic = a.aromatic;
        let bond_sum = mol.explicit_valence(i);
        let valences = element::standard_valences(z);
        if valences.is_empty() {
            continue;
        }
        // Aromatic atoms: count one of the ring's π electrons as part
        // of the σ framework — an aromatic carbon with two ring
        // neighbours wants one H. We treat aromatic bond sum (1.5+1.5)
        // = 3 and a target of 4 → 1 H. That falls out naturally from
        // BondOrder::Aromatic = 1.5, so no special case is needed
        // beyond rounding.
        let need = pick_valence(valences, bond_sum, charge, z);
        let filled = (need - bond_sum).round();
        let h = if filled > 0.0 { filled as u8 } else { 0 };
        mol.atoms[i].implicit_h = h;
        // an aromatic atom with no ring H still keeps aromatic flag
        let _ = aromatic;
    }
}

/// Choose the standard valence an atom is aiming for, given its current
/// bond-order sum and formal charge. Picks the smallest tabulated
/// valence ≥ the bond sum (after the charge offset); falls back to the
/// largest if the atom is already hypervalent.
fn pick_valence(valences: &[u8], bond_sum: f64, charge: i8, z: u8) -> f64 {
    // Charge offset: a cation of an electron-rich atom (N+, O+) gains
    // a bond slot; an anion loses one. For carbon a + or - removes a
    // valence. We model the common organic cases.
    let offset = match z {
        7 | 15 => i32::from(charge),            // N+, P+ → +1 valence
        8 | 16 => i32::from(charge),            // O+, S+ → +1; O- → -1
        6 => -i32::from(charge.unsigned_abs()), // C+, C- both lose a slot
        _ => i32::from(charge),
    };
    for &v in valences {
        let target = f64::from(v) + f64::from(offset);
        if target >= bond_sum - 1e-6 {
            return target.max(0.0);
        }
    }
    let last = f64::from(*valences.last().unwrap()) + f64::from(offset);
    last.max(bond_sum)
}

/// Ensure every bond inside an all-aromatic ring run carries
/// [`BondOrder::Aromatic`]. The parser already sets aromatic bonds for
/// adjacent lowercase atoms, but a ring-closure bond between two
/// aromatic atoms written far apart may have come through as `Single`;
/// this promotes those.
fn fixup_aromatic_bonds(mol: &mut Molecule) {
    for bi in 0..mol.bonds.len() {
        let (a, b) = (mol.bonds[bi].a, mol.bonds[bi].b);
        if mol.atoms[a].aromatic
            && mol.atoms[b].aromatic
            && mol.bonds[bi].order == BondOrder::Single
        {
            // Only promote if both atoms truly sit in a ring — a
            // biaryl single bond joins two aromatic atoms but is not
            // itself aromatic. Use a quick ring test: the bond is in a
            // ring iff a path a→b avoiding this bond exists.
            if path_exists_excluding(mol, a, b, bi) {
                mol.bonds[bi].order = BondOrder::Aromatic;
                mol.bonds[bi].aromatic = true;
            }
        }
    }
}

/// True if a bond-path connects `from` and `to` without using bond
/// `excl`. Used to tell a ring bond from a chain / biaryl bond.
fn path_exists_excluding(mol: &Molecule, from: usize, to: usize, excl: usize) -> bool {
    let mut seen = vec![false; mol.atoms.len()];
    let mut stack = vec![from];
    seen[from] = true;
    while let Some(u) = stack.pop() {
        for (bi, bond) in mol.bonds.iter().enumerate() {
            if bi == excl {
                continue;
            }
            if let Some(v) = bond.other(u) {
                if v == to {
                    return true;
                }
                if !seen[v] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ethanol() {
        let m = parse_smiles("CCO").unwrap();
        assert_eq!(m.atom_count(), 3);
        assert_eq!(m.bond_count(), 2);
        // CH3-CH2-OH : implicit H = 3, 2, 1
        assert_eq!(m.atoms[0].implicit_h, 3);
        assert_eq!(m.atoms[1].implicit_h, 2);
        assert_eq!(m.atoms[2].implicit_h, 1);
    }

    #[test]
    fn branches_and_rings() {
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(benzene.atom_count(), 6);
        assert_eq!(benzene.bond_count(), 6);
        assert!(benzene.atoms.iter().all(|a| a.aromatic));
        assert!(benzene.bonds.iter().all(|b| b.order == BondOrder::Aromatic));
        // each aromatic CH carries one implicit H
        assert!(benzene.atoms.iter().all(|a| a.implicit_h == 1));

        let isobutane = parse_smiles("CC(C)C").unwrap();
        assert_eq!(isobutane.atom_count(), 4);
        assert_eq!(isobutane.degree(1), 3);
    }

    #[test]
    fn bracket_atoms() {
        let m = parse_smiles("[NH4+]").unwrap();
        assert_eq!(m.atoms[0].atomic_number, 7);
        assert_eq!(m.atoms[0].explicit_h, 4);
        assert_eq!(m.atoms[0].formal_charge, 1);

        let m = parse_smiles("[13CH4]").unwrap();
        assert_eq!(m.atoms[0].isotope, Some(13));
        assert_eq!(m.atoms[0].explicit_h, 4);

        let m = parse_smiles("[O-]").unwrap();
        assert_eq!(m.atoms[0].formal_charge, -1);

        let m = parse_smiles("[Fe+3]").unwrap();
        assert_eq!(m.atoms[0].atomic_number, 26);
        assert_eq!(m.atoms[0].formal_charge, 3);
    }

    #[test]
    fn bond_orders_and_fragments() {
        let m = parse_smiles("C=C").unwrap();
        assert_eq!(m.bonds[0].order, BondOrder::Double);
        assert_eq!(m.atoms[0].implicit_h, 2);

        let m = parse_smiles("C#N").unwrap();
        assert_eq!(m.bonds[0].order, BondOrder::Triple);

        let salt = parse_smiles("[Na+].[Cl-]").unwrap();
        assert_eq!(salt.atom_count(), 2);
        assert_eq!(salt.bond_count(), 0);
        assert_eq!(salt.component_count(), 2);
    }

    #[test]
    fn chirality_and_maps() {
        let m = parse_smiles("[C@H](N)(O)C").unwrap();
        assert_eq!(m.atoms[0].chirality, Chirality::Ccw);
        let m = parse_smiles("[C@@H](N)(O)C").unwrap();
        assert_eq!(m.atoms[0].chirality, Chirality::Cw);
        let m = parse_smiles("[CH3:1]").unwrap();
        assert_eq!(m.atoms[0].map_number, 1);
    }

    #[test]
    fn two_digit_ring_closure() {
        let m = parse_smiles("C%10CCCCC%10").unwrap();
        assert_eq!(m.atom_count(), 6);
        assert_eq!(m.bond_count(), 6);
    }

    #[test]
    fn rejects_malformed() {
        assert!(parse_smiles("").is_err());
        assert!(parse_smiles("C(").is_err());
        assert!(parse_smiles("C)").is_err());
        assert!(parse_smiles("C1CC").is_err()); // unclosed ring
        assert!(parse_smiles("[Xx]").is_err()); // unknown element
        assert!(parse_smiles("C==C").is_err()); // double bond symbol
    }

    #[test]
    fn pyridine_aromatic_nitrogen() {
        let m = parse_smiles("c1ccncc1").unwrap();
        assert_eq!(m.atom_count(), 6);
        let n = m.atoms.iter().find(|a| a.atomic_number == 7).unwrap();
        assert!(n.aromatic);
        // aromatic N with two ring bonds carries no H
        assert_eq!(n.implicit_h, 0);
    }
}

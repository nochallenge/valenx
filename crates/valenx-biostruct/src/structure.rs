//! The macromolecular-structure hierarchy:
//! [`Structure`] > [`Model`] > [`Chain`] > [`Residue`] > [`Atom`].
//!
//! This mirrors the Biopython `Structure / Model / Chain / Residue /
//! Atom` (SMCRA) tree and the mmCIF data model. Containers are plain
//! `Vec`s — for the structure sizes a desktop tool handles (tens of
//! thousands of atoms) linear scans are fast enough and keep the model
//! trivially `serde`-serialisable and clonable.
//!
//! ## Identity
//!
//! - A [`Residue`] is keyed by `(seq_num, ins_code)` — the PDB
//!   residue sequence number plus the single-character insertion
//!   code. Two residues `100` and `100A` are distinct.
//! - An [`Atom`] carries an `alt_loc` character. The empty `' '`
//!   alt-loc means "no alternate". Most analysis takes the
//!   [highest-occupancy](Residue::primary_atoms) alt-loc per name.
//! - `hetatm` distinguishes `HETATM` records (ligands, ions, waters,
//!   modified residues) from polymer `ATOM` records.

use nalgebra::Point3;
use serde::{Deserialize, Serialize};

/// A single atom: an element, a name, a 3-D coordinate and the
/// crystallographic occupancy / B-factor / alt-loc metadata.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// PDB atom serial number (column 7-11). Informational only — the
    /// crate never keys on it.
    pub serial: i32,
    /// Atom name, whitespace-trimmed (`"CA"`, `"N"`, `"OD1"`, `"C1'"`).
    pub name: String,
    /// Alternate-location indicator. `' '` means "no alternate".
    pub alt_loc: char,
    /// Element symbol, upper-cased (`"C"`, `"N"`, `"FE"`). Derived
    /// from the PDB element column when present, else guessed from the
    /// atom name.
    pub element: String,
    /// Cartesian coordinate in ångström.
    pub coord: Point3<f64>,
    /// Crystallographic occupancy in `[0, 1]`.
    pub occupancy: f64,
    /// Isotropic B-factor (temperature factor), ångström².
    pub b_factor: f64,
    /// Formal charge as written in the PDB charge column (column
    /// 79-80), e.g. `+1`. `0` when blank.
    pub charge: i32,
}

impl Atom {
    /// Construct an atom with default occupancy `1.0`, B-factor `0.0`,
    /// no alt-loc and zero charge.
    pub fn new(name: impl Into<String>, element: impl Into<String>, coord: Point3<f64>) -> Self {
        Atom {
            serial: 0,
            name: name.into(),
            alt_loc: ' ',
            element: element.into(),
            coord,
            occupancy: 1.0,
            b_factor: 0.0,
            charge: 0,
        }
    }

    /// Euclidean distance to another atom, in ångström.
    pub fn distance(&self, other: &Atom) -> f64 {
        (self.coord - other.coord).norm()
    }

    /// Whether this atom is a hydrogen (element `H` or `D`).
    pub fn is_hydrogen(&self) -> bool {
        matches!(self.element.as_str(), "H" | "D")
    }
}

/// Polymer type a [`Residue`] belongs to, inferred from its name.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum ResidueKind {
    /// One of the 20 standard amino acids (or a recognised modified
    /// residue such as `MSE`).
    AminoAcid,
    /// A DNA nucleotide (`DA`, `DC`, `DG`, `DT`, `DI`).
    Dna,
    /// An RNA nucleotide (`A`, `C`, `G`, `U`, `I`).
    Rna,
    /// A water molecule (`HOH`, `WAT`, `DOD`).
    Water,
    /// Anything else — a ligand, ion or unrecognised group.
    Other,
}

/// A residue: a named group of atoms keyed by sequence number and
/// insertion code.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Residue {
    /// Residue name, trimmed and upper-cased (`"ALA"`, `"DA"`,
    /// `"HOH"`).
    pub name: String,
    /// PDB residue sequence number (column 23-26).
    pub seq_num: i32,
    /// Insertion code (column 27). `' '` when absent.
    pub ins_code: char,
    /// `true` for `HETATM` groups, `false` for polymer `ATOM` records.
    pub hetatm: bool,
    /// Member atoms, in file order.
    pub atoms: Vec<Atom>,
}

impl Residue {
    /// Construct an empty residue.
    pub fn new(name: impl Into<String>, seq_num: i32) -> Self {
        Residue {
            name: name.into(),
            seq_num,
            ins_code: ' ',
            hetatm: false,
            atoms: Vec::new(),
        }
    }

    /// The `(seq_num, ins_code)` identity tuple.
    pub fn id(&self) -> (i32, char) {
        (self.seq_num, self.ins_code)
    }

    /// Classify this residue's polymer type from its name.
    pub fn kind(&self) -> ResidueKind {
        classify_residue(&self.name)
    }

    /// Whether this residue is one of the 20 standard amino acids or a
    /// recognised modified residue.
    pub fn is_amino_acid(&self) -> bool {
        self.kind() == ResidueKind::AminoAcid
    }

    /// Whether this residue is a DNA or RNA nucleotide.
    pub fn is_nucleotide(&self) -> bool {
        matches!(self.kind(), ResidueKind::Dna | ResidueKind::Rna)
    }

    /// First atom with the given (trimmed) name, ignoring alt-loc.
    pub fn atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name == name)
    }

    /// Highest-occupancy atom with the given name. When several
    /// alt-locs share a name this returns the dominant conformer; ties
    /// resolve to the first encountered.
    pub fn primary_atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().filter(|a| a.name == name).max_by(|a, b| {
            a.occupancy
                .partial_cmp(&b.occupancy)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// One representative atom per distinct atom name (the
    /// highest-occupancy alt-loc of each). Returned in first-seen
    /// order. This is the alt-loc-collapsed view most geometry uses.
    pub fn primary_atoms(&self) -> Vec<&Atom> {
        let mut seen: Vec<&str> = Vec::new();
        let mut out: Vec<&Atom> = Vec::new();
        for atom in &self.atoms {
            if seen.contains(&atom.name.as_str()) {
                continue;
            }
            seen.push(&atom.name);
            if let Some(best) = self.primary_atom(&atom.name) {
                out.push(best);
            }
        }
        out
    }

    /// The alpha-carbon (`CA`) for amino acids, highest-occupancy.
    pub fn ca(&self) -> Option<&Atom> {
        self.primary_atom("CA")
    }

    /// Mean coordinate of all atoms (unweighted geometric centre).
    pub fn centroid(&self) -> Option<Point3<f64>> {
        if self.atoms.is_empty() {
            return None;
        }
        let mut acc = nalgebra::Vector3::zeros();
        for a in &self.atoms {
            acc += a.coord.coords;
        }
        Some(Point3::from(acc / self.atoms.len() as f64))
    }
}

/// A chain: an ordered list of residues sharing a one- (or multi-)
/// character chain identifier.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Chain {
    /// Chain identifier (`"A"`, `"B"`; mmCIF allows multi-character).
    pub id: String,
    /// Residues in file / sequence order.
    pub residues: Vec<Residue>,
    /// `SEQRES`-declared one-letter sequence, when the file carried
    /// `SEQRES` records. This is the *full* construct sequence
    /// including residues missing from the coordinates.
    pub seqres: Option<String>,
}

impl Chain {
    /// Construct an empty chain.
    pub fn new(id: impl Into<String>) -> Self {
        Chain {
            id: id.into(),
            residues: Vec::new(),
            seqres: None,
        }
    }

    /// Find a residue by its `(seq_num, ins_code)` identity.
    pub fn residue(&self, seq_num: i32, ins_code: char) -> Option<&Residue> {
        self.residues
            .iter()
            .find(|r| r.seq_num == seq_num && r.ins_code == ins_code)
    }

    /// All polymer (amino-acid or nucleotide) residues, in order.
    pub fn polymer_residues(&self) -> Vec<&Residue> {
        self.residues
            .iter()
            .filter(|r| r.is_amino_acid() || r.is_nucleotide())
            .collect()
    }

    /// One-letter sequence of the *observed* residues (those present
    /// in the coordinates), in chain order. Amino acids map to the
    /// 20-letter code, nucleotides to `ACGTU`, unknowns to `X`.
    pub fn observed_sequence(&self) -> String {
        self.residues
            .iter()
            .filter(|r| r.is_amino_acid() || r.is_nucleotide())
            .map(|r| residue_one_letter(&r.name))
            .collect()
    }
}

/// A model: one full coordinate set. Crystal structures have a single
/// model; NMR ensembles carry many (`MODEL` / `ENDMDL` records).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Model {
    /// Model serial number (PDB `MODEL` record). `1` for the default.
    pub serial: i32,
    /// Chains in file order.
    pub chains: Vec<Chain>,
}

impl Model {
    /// Construct an empty model.
    pub fn new(serial: i32) -> Self {
        Model {
            serial,
            chains: Vec::new(),
        }
    }

    /// Find a chain by identifier.
    pub fn chain(&self, id: &str) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    /// Mutable chain lookup.
    pub fn chain_mut(&mut self, id: &str) -> Option<&mut Chain> {
        self.chains.iter_mut().find(|c| c.id == id)
    }

    /// Every residue across every chain, in order.
    pub fn residues(&self) -> impl Iterator<Item = &Residue> {
        self.chains.iter().flat_map(|c| c.residues.iter())
    }

    /// Every atom across every chain / residue, in order.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.chains
            .iter()
            .flat_map(|c| c.residues.iter())
            .flat_map(|r| r.atoms.iter())
    }

    /// Total atom count.
    pub fn atom_count(&self) -> usize {
        self.atoms().count()
    }
}

/// A `BIOMT`-style symmetry operator: a 3×3 rotation plus a 3-vector
/// translation, applied as `x' = R·x + t`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SymmetryOperator {
    /// Operator serial number as written in the `REMARK 350` block.
    pub serial: i32,
    /// 3×3 rotation matrix, row-major as `[[r00,r01,r02], …]`.
    pub rotation: [[f64; 3]; 3],
    /// Translation vector in ångström.
    pub translation: [f64; 3],
}

impl SymmetryOperator {
    /// The identity operator (no rotation, no translation).
    pub fn identity() -> Self {
        SymmetryOperator {
            serial: 1,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0, 0.0, 0.0],
        }
    }

    /// Apply the operator to a point.
    pub fn apply(&self, p: &Point3<f64>) -> Point3<f64> {
        let r = &self.rotation;
        let x = r[0][0] * p.x + r[0][1] * p.y + r[0][2] * p.z + self.translation[0];
        let y = r[1][0] * p.x + r[1][1] * p.y + r[1][2] * p.z + self.translation[1];
        let z = r[2][0] * p.x + r[2][1] * p.y + r[2][2] * p.z + self.translation[2];
        Point3::new(x, y, z)
    }

    /// Whether this operator is (numerically) the identity.
    pub fn is_identity(&self) -> bool {
        let r = &self.rotation;
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for i in 0..3 {
            for j in 0..3 {
                if (r[i][j] - id[i][j]).abs() > 1e-6 {
                    return false;
                }
            }
            if self.translation[i].abs() > 1e-6 {
                return false;
            }
        }
        true
    }
}

/// A disulfide bond between two cysteine residues, as parsed from a
/// PDB `SSBOND` record. Each partner is `(chain, seq_num, ins_code)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Disulfide {
    /// First cysteine, `(chain, seq_num, ins_code)`.
    pub partner_a: (String, i32, char),
    /// Second cysteine, `(chain, seq_num, ins_code)`.
    pub partner_b: (String, i32, char),
}

/// A complete macromolecular structure: an identifier, one or more
/// coordinate [`Model`]s and the header-derived metadata.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Structure {
    /// 4-character PDB id or any user label.
    pub id: String,
    /// Free-text title (PDB `TITLE`, mmCIF `_struct.title`).
    pub title: String,
    /// Coordinate models. Always at least one after a successful
    /// parse.
    pub models: Vec<Model>,
    /// `REMARK 350` biological-assembly operators, when present.
    pub assembly_operators: Vec<SymmetryOperator>,
    /// `HELIX`-record secondary-structure spans, as
    /// `(chain, start_seq, end_seq)`.
    pub helix_records: Vec<(String, i32, i32)>,
    /// `SHEET`-record strand spans, as `(chain, start_seq, end_seq)`.
    pub sheet_records: Vec<(String, i32, i32)>,
    /// `SSBOND`-record disulfide bonds, when present.
    pub disulfides: Vec<Disulfide>,
    /// `CONECT`-record explicit connectivity, as `(serial, [partner
    /// serials])` — atom serial numbers, not indices. Used mainly for
    /// hetero-group bonds the PDB writer round-trips.
    pub conect: Vec<(i32, Vec<i32>)>,
}

impl Structure {
    /// Construct an empty structure with a single empty model.
    pub fn new(id: impl Into<String>) -> Self {
        Structure {
            id: id.into(),
            title: String::new(),
            models: vec![Model::new(1)],
            assembly_operators: Vec::new(),
            helix_records: Vec::new(),
            sheet_records: Vec::new(),
            disulfides: Vec::new(),
            conect: Vec::new(),
        }
    }

    /// The first model — the one geometry analyses use by default.
    /// Panics only if `models` is empty, which a successful parse
    /// never produces.
    pub fn first_model(&self) -> &Model {
        &self.models[0]
    }

    /// Mutable first model.
    pub fn first_model_mut(&mut self) -> &mut Model {
        &mut self.models[0]
    }

    /// Total atom count across the first model.
    pub fn atom_count(&self) -> usize {
        self.models.first().map(|m| m.atom_count()).unwrap_or(0)
    }

    /// Validate the hierarchy: at least one model, each model has at
    /// least one chain, each chain has at least one residue, each
    /// residue has at least one atom.
    pub fn validate(&self) -> crate::error::Result<()> {
        use crate::error::BiostructError;
        if self.models.is_empty() {
            return Err(BiostructError::invalid_structure("structure has no models"));
        }
        for m in &self.models {
            if m.chains.is_empty() {
                return Err(BiostructError::invalid_structure(format!(
                    "model {} has no chains",
                    m.serial
                )));
            }
            for c in &m.chains {
                if c.residues.is_empty() {
                    return Err(BiostructError::invalid_structure(format!(
                        "chain {} has no residues",
                        c.id
                    )));
                }
                for r in &c.residues {
                    if r.atoms.is_empty() {
                        return Err(BiostructError::invalid_structure(format!(
                            "residue {} {} has no atoms",
                            r.name, r.seq_num
                        )));
                    }
                }
            }
        }
        Ok(())
    }
}

// --- Residue classification & one-letter codes -----------------------

/// The 20 standard amino-acid three-letter codes plus a handful of
/// common modified residues that should still classify as amino acids.
const AMINO_ACIDS: &[(&str, char)] = &[
    ("ALA", 'A'),
    ("ARG", 'R'),
    ("ASN", 'N'),
    ("ASP", 'D'),
    ("CYS", 'C'),
    ("GLN", 'Q'),
    ("GLU", 'E'),
    ("GLY", 'G'),
    ("HIS", 'H'),
    ("ILE", 'I'),
    ("LEU", 'L'),
    ("LYS", 'K'),
    ("MET", 'M'),
    ("PHE", 'F'),
    ("PRO", 'P'),
    ("SER", 'S'),
    ("THR", 'T'),
    ("TRP", 'W'),
    ("TYR", 'Y'),
    ("VAL", 'V'),
    // common modified residues — map to their parent
    ("MSE", 'M'), // selenomethionine
    ("SEC", 'U'), // selenocysteine
    ("PYL", 'O'), // pyrrolysine
    ("HYP", 'P'), // hydroxyproline
    ("HSD", 'H'), // CHARMM histidine protonation states
    ("HSE", 'H'),
    ("HSP", 'H'),
    ("CSO", 'C'),
    ("PTR", 'Y'),
    ("SEP", 'S'),
    ("TPO", 'T'),
];

const DNA_BASES: &[(&str, char)] = &[
    ("DA", 'A'),
    ("DC", 'C'),
    ("DG", 'G'),
    ("DT", 'T'),
    ("DI", 'I'),
    ("DU", 'U'),
];

const RNA_BASES: &[(&str, char)] = &[
    ("A", 'A'),
    ("C", 'C'),
    ("G", 'G'),
    ("U", 'U'),
    ("I", 'I'),
    ("RA", 'A'),
    ("RC", 'C'),
    ("RG", 'G'),
    ("RU", 'U'),
];

const WATER_NAMES: &[&str] = &["HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"];

/// Classify a residue name into a [`ResidueKind`].
pub fn classify_residue(name: &str) -> ResidueKind {
    let n = name.trim().to_ascii_uppercase();
    if WATER_NAMES.contains(&n.as_str()) {
        return ResidueKind::Water;
    }
    if AMINO_ACIDS.iter().any(|(c, _)| *c == n) {
        return ResidueKind::AminoAcid;
    }
    if DNA_BASES.iter().any(|(c, _)| *c == n) {
        return ResidueKind::Dna;
    }
    if RNA_BASES.iter().any(|(c, _)| *c == n) {
        return ResidueKind::Rna;
    }
    ResidueKind::Other
}

/// One-letter code for a residue name. Amino acids map to the
/// 20-letter alphabet, nucleotides to `ACGTU`, everything else to
/// `X`.
pub fn residue_one_letter(name: &str) -> char {
    let n = name.trim().to_ascii_uppercase();
    if let Some((_, c)) = AMINO_ACIDS.iter().find(|(code, _)| *code == n) {
        return *c;
    }
    if let Some((_, c)) = DNA_BASES.iter().find(|(code, _)| *code == n) {
        return *c;
    }
    if let Some((_, c)) = RNA_BASES.iter().find(|(code, _)| *code == n) {
        return *c;
    }
    'X'
}

/// Whether `name` is a recognised nucleotide three-letter code.
pub fn is_nucleotide_name(name: &str) -> bool {
    matches!(classify_residue(name), ResidueKind::Dna | ResidueKind::Rna)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aa_residue(name: &str, seq: i32) -> Residue {
        let mut r = Residue::new(name, seq);
        r.atoms
            .push(Atom::new("N", "N", Point3::new(0.0, 0.0, 0.0)));
        r.atoms
            .push(Atom::new("CA", "C", Point3::new(1.5, 0.0, 0.0)));
        r.atoms
            .push(Atom::new("C", "C", Point3::new(2.5, 1.0, 0.0)));
        r
    }

    #[test]
    fn classification() {
        assert_eq!(classify_residue("ALA"), ResidueKind::AminoAcid);
        assert_eq!(classify_residue(" mse "), ResidueKind::AminoAcid);
        assert_eq!(classify_residue("DA"), ResidueKind::Dna);
        assert_eq!(classify_residue("U"), ResidueKind::Rna);
        assert_eq!(classify_residue("HOH"), ResidueKind::Water);
        assert_eq!(classify_residue("ZN"), ResidueKind::Other);
    }

    #[test]
    fn one_letter_codes() {
        assert_eq!(residue_one_letter("TRP"), 'W');
        assert_eq!(residue_one_letter("MSE"), 'M');
        assert_eq!(residue_one_letter("DG"), 'G');
        assert_eq!(residue_one_letter("LIG"), 'X');
    }

    #[test]
    fn residue_atom_lookup() {
        let r = aa_residue("ALA", 1);
        assert!(r.ca().is_some());
        assert_eq!(r.atom("N").unwrap().element, "N");
        assert!(r.atom("CB").is_none());
        assert_eq!(r.id(), (1, ' '));
    }

    #[test]
    fn primary_atom_picks_highest_occupancy() {
        let mut r = Residue::new("SER", 5);
        let mut lo = Atom::new("OG", "O", Point3::new(0.0, 0.0, 0.0));
        lo.alt_loc = 'A';
        lo.occupancy = 0.3;
        let mut hi = Atom::new("OG", "O", Point3::new(1.0, 0.0, 0.0));
        hi.alt_loc = 'B';
        hi.occupancy = 0.7;
        r.atoms.push(lo);
        r.atoms.push(hi);
        let p = r.primary_atom("OG").unwrap();
        assert!((p.occupancy - 0.7).abs() < 1e-9);
        assert_eq!(r.primary_atoms().len(), 1);
    }

    #[test]
    fn chain_observed_sequence() {
        let mut c = Chain::new("A");
        c.residues.push(aa_residue("ALA", 1));
        c.residues.push(aa_residue("GLY", 2));
        c.residues.push(aa_residue("TRP", 3));
        assert_eq!(c.observed_sequence(), "AGW");
        assert_eq!(c.polymer_residues().len(), 3);
    }

    #[test]
    fn structure_validation() {
        let mut s = Structure::new("TEST");
        assert!(s.validate().is_err()); // empty model
        let mut c = Chain::new("A");
        c.residues.push(aa_residue("ALA", 1));
        s.first_model_mut().chains.push(c);
        assert!(s.validate().is_ok());
        assert_eq!(s.atom_count(), 3);
    }

    #[test]
    fn symmetry_operator_identity() {
        let op = SymmetryOperator::identity();
        assert!(op.is_identity());
        let p = Point3::new(3.0, -2.0, 7.0);
        let q = op.apply(&p);
        assert!((p - q).norm() < 1e-12);
    }

    #[test]
    fn symmetry_operator_translation() {
        let mut op = SymmetryOperator::identity();
        op.translation = [1.0, 2.0, 3.0];
        assert!(!op.is_identity());
        let q = op.apply(&Point3::origin());
        assert!((q - Point3::new(1.0, 2.0, 3.0)).norm() < 1e-12);
    }
}

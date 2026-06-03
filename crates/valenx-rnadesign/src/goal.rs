//! Features 1–2 — the design goal and the design constraints.
//!
//! Every synthetic-RNA design begins with two declarations:
//!
//! - a [`DesignGoal`] — *what* the RNA must do. Three families:
//!   - a **structural RNA** ([`DesignGoal::Structural`]) — fold to a
//!     given target dot-bracket (an aptamer, ribozyme, riboswitch,
//!     engineered hairpin or a custom structure);
//!   - a **coding mRNA** ([`DesignGoal::Coding`]) — translate a given
//!     protein sequence, assembled into a five-part therapeutic mRNA;
//!   - a **hybrid** ([`DesignGoal::Hybrid`]) — a coding RNA that must
//!     *also* present required structural elements (a structured 5′UTR
//!     motif, a 3′-end hairpin).
//! - a [`DesignConstraints`] set — the *boundaries* every candidate
//!   must respect: GC range, length window, forbidden motifs, forbidden
//!   restriction sites, required / forbidden subsequences, the host
//!   organism (for codon usage) and the modified-nucleoside scheme.
//!
//! Both types validate themselves ([`DesignGoal::validate`],
//! [`DesignConstraints::validate`]) so a malformed goal is caught at
//! the front of the workflow rather than deep inside the optimiser.
//!
//! ## v1 scope — honest framing
//!
//! A [`DesignGoal`] is a *specification*, not a promise. The downstream
//! workflow produces a sequence predicted — by a nearest-neighbor
//! energy model — to satisfy the goal; it is a strong in-silico
//! candidate, never a guarantee of in-vivo behaviour. See the crate
//! root for the full disclaimer.

use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use valenx_genediting::mrna::codon::ExpressionHost;
use valenx_genediting::mrna::tailcap::MrnaUseCase;
use valenx_genediting::mrna::uridine::ModifiedNucleoside;
use valenx_rnastruct::Structure;

/// The functional family a structural-RNA design targets.
///
/// Purely descriptive — it labels the intent and tunes a few defaults
/// (a ribozyme wants a stable catalytic core, an aptamer a defined
/// binding pocket) but the folding target is always the dot-bracket the
/// caller supplies.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum StructuralClass {
    /// An aptamer — an RNA that folds into a ligand-binding pocket.
    Aptamer,
    /// A ribozyme — a catalytically active RNA.
    Ribozyme,
    /// A riboswitch — a regulatory RNA with a ligand-sensing aptamer
    /// domain (a two-state design; see [`crate::design::riboswitch`]).
    Riboswitch,
    /// An engineered hairpin / stem-loop element.
    Hairpin,
    /// A custom structure with no standard functional label.
    Custom,
}

impl StructuralClass {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            StructuralClass::Aptamer => "aptamer",
            StructuralClass::Ribozyme => "ribozyme",
            StructuralClass::Riboswitch => "riboswitch",
            StructuralClass::Hairpin => "hairpin",
            StructuralClass::Custom => "custom structure",
        }
    }
}

/// A required structural element a [`DesignGoal::Hybrid`] design must
/// present somewhere in the transcript.
///
/// The element is given as a dot-bracket fragment plus a coarse
/// placement hint. The hybrid designer seats the element and checks the
/// assembled transcript folds it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralElement {
    /// A short dot-bracket fragment the element must fold into.
    pub dot_bracket: String,
    /// Where in the transcript the element belongs.
    pub placement: ElementPlacement,
    /// A human-readable label for the element (`"5' hairpin"`, …).
    pub label: String,
}

/// Coarse placement of a [`StructuralElement`] within a transcript.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementPlacement {
    /// In the 5′ untranslated region.
    FivePrimeUtr,
    /// In the 3′ untranslated region.
    ThreePrimeUtr,
}

/// The design goal — what the RNA must do (feature 1).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum DesignGoal {
    /// Design a structural RNA that folds to `target`.
    Structural {
        /// The target secondary structure, as a dot-bracket-derived
        /// [`Structure`].
        target: Structure,
        /// The functional family (labels intent, tunes defaults).
        class: StructuralClass,
    },
    /// Design a coding mRNA that translates `protein`, assembled as a
    /// five-part therapeutic construct.
    Coding {
        /// The protein to encode, as one-letter amino-acid codes
        /// (a trailing `*` stop is allowed but not required).
        protein: Vec<u8>,
        /// The intended therapeutic use (drives poly-A / cap choice).
        use_case: MrnaUseCase,
    },
    /// Design a coding RNA that *also* presents required structural
    /// elements.
    Hybrid {
        /// The protein to encode.
        protein: Vec<u8>,
        /// The intended therapeutic use.
        use_case: MrnaUseCase,
        /// Structural elements the assembled transcript must fold.
        elements: Vec<StructuralElement>,
    },
}

impl DesignGoal {
    /// A structural-RNA goal from a target dot-bracket string.
    ///
    /// # Errors
    /// [`RnaDesignError::Goal`] if the dot-bracket fails to parse or is
    /// pseudoknotted (the inverse folder is pseudoknot-free).
    pub fn structural(dot_bracket: &str, class: StructuralClass) -> Result<Self> {
        let target = Structure::from_dot_bracket(dot_bracket)
            .map_err(|e| RnaDesignError::goal("target", e.to_string()))?;
        if target.is_empty() {
            return Err(RnaDesignError::goal("target", "target structure is empty"));
        }
        if target.has_pseudoknot() {
            return Err(RnaDesignError::goal(
                "target",
                "target is pseudoknotted — the inverse folder is pseudoknot-free",
            ));
        }
        Ok(DesignGoal::Structural { target, class })
    }

    /// A coding-mRNA goal encoding `protein`.
    pub fn coding(protein: impl Into<Vec<u8>>, use_case: MrnaUseCase) -> Self {
        DesignGoal::Coding {
            protein: protein.into(),
            use_case,
        }
    }

    /// A short human-readable name for the goal kind.
    pub fn kind_name(&self) -> &'static str {
        match self {
            DesignGoal::Structural { .. } => "structural RNA",
            DesignGoal::Coding { .. } => "coding mRNA",
            DesignGoal::Hybrid { .. } => "hybrid coding RNA with structural elements",
        }
    }

    /// `true` when the goal involves a coding sequence (a `Coding` or
    /// `Hybrid` goal) — the designs that run codon optimisation.
    pub fn is_coding(&self) -> bool {
        matches!(self, DesignGoal::Coding { .. } | DesignGoal::Hybrid { .. })
    }

    /// Validates the goal — a non-empty, non-pseudoknotted target for a
    /// structural goal; a non-empty protein for a coding / hybrid goal.
    ///
    /// # Errors
    /// [`RnaDesignError::Goal`] describing the first problem found.
    pub fn validate(&self) -> Result<()> {
        match self {
            DesignGoal::Structural { target, .. } => {
                if target.is_empty() {
                    return Err(RnaDesignError::goal("target", "target structure is empty"));
                }
                if target.has_pseudoknot() {
                    return Err(RnaDesignError::goal(
                        "target",
                        "target is pseudoknotted — unreachable by the MFE folder",
                    ));
                }
            }
            DesignGoal::Coding { protein, .. } => validate_protein(protein)?,
            DesignGoal::Hybrid {
                protein, elements, ..
            } => {
                validate_protein(protein)?;
                if elements.is_empty() {
                    return Err(RnaDesignError::goal(
                        "elements",
                        "a hybrid goal must list at least one structural element",
                    ));
                }
                for el in elements {
                    Structure::from_dot_bracket(&el.dot_bracket).map_err(|e| {
                        RnaDesignError::goal("element", format!("`{}`: {e}", el.label))
                    })?;
                }
            }
        }
        Ok(())
    }
}

/// Validates a one-letter protein sequence: non-empty, every residue a
/// legal amino-acid code (or a single trailing stop `*`).
fn validate_protein(protein: &[u8]) -> Result<()> {
    if protein.is_empty() {
        return Err(RnaDesignError::goal("protein", "protein sequence is empty"));
    }
    for (i, &b) in protein.iter().enumerate() {
        let u = b.to_ascii_uppercase();
        let is_aa = u.is_ascii_uppercase() && u != b'B' && u != b'J' && u != b'O' && u != b'U'
            && u != b'X' && u != b'Z';
        let is_stop = u == b'*' && i + 1 == protein.len();
        if !is_aa && !is_stop {
            return Err(RnaDesignError::goal(
                "protein",
                format!("`{}` at position {i} is not a standard amino-acid code", b as char),
            ));
        }
    }
    Ok(())
}

/// The design constraints — the boundaries every candidate must respect
/// (feature 2).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignConstraints {
    /// Minimum acceptable GC fraction in `[0, 1]`.
    pub gc_min: f64,
    /// Maximum acceptable GC fraction in `[0, 1]`.
    pub gc_max: f64,
    /// Minimum acceptable sequence length, nt (`0` = no lower bound).
    pub length_min: usize,
    /// Maximum acceptable sequence length, nt (`0` = no upper bound).
    pub length_max: usize,
    /// Sequence motifs that must not appear anywhere (RNA or DNA
    /// alphabet — both are matched after transcription to RNA).
    pub forbidden_motifs: Vec<String>,
    /// Names of restriction enzymes (from the `valenx-bioseq` database)
    /// whose recognition sites must not appear in the DNA template.
    pub forbidden_restriction_sites: Vec<String>,
    /// Subsequences that must appear at least once.
    pub required_subsequences: Vec<String>,
    /// Subsequences that must not appear (a superset use of
    /// `forbidden_motifs` for caller-supplied named sequences).
    pub forbidden_subsequences: Vec<String>,
    /// The expression host — drives codon usage for a coding design.
    pub host: ExpressionHost,
    /// The modified-nucleoside scheme the physical RNA will be made
    /// with (`m1Ψ` for a therapeutic construct).
    pub nucleoside: ModifiedNucleoside,
    /// Maximum acceptable homopolymer-run length (a synthesizability
    /// limit — long single-base runs are hard to synthesise / transcribe).
    pub max_homopolymer: usize,
}

impl Default for DesignConstraints {
    /// Sensible therapeutic-RNA defaults: GC 40–60 %, no length bounds,
    /// no forbidden motifs, human host, `m1Ψ`, homopolymer runs capped
    /// at 6 nt.
    fn default() -> Self {
        DesignConstraints {
            gc_min: 0.40,
            gc_max: 0.60,
            length_min: 0,
            length_max: 0,
            forbidden_motifs: Vec::new(),
            forbidden_restriction_sites: Vec::new(),
            required_subsequences: Vec::new(),
            forbidden_subsequences: Vec::new(),
            host: ExpressionHost::Human,
            nucleoside: ModifiedNucleoside::N1MethylPseudouridine,
            max_homopolymer: 6,
        }
    }
}

impl DesignConstraints {
    /// A fresh constraint set with the therapeutic defaults
    /// ([`DesignConstraints::default`]).
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the GC fraction window.
    pub fn with_gc_range(mut self, gc_min: f64, gc_max: f64) -> Self {
        self.gc_min = gc_min;
        self.gc_max = gc_max;
        self
    }

    /// Sets the length window (`0` for either bound = unbounded).
    pub fn with_length_range(mut self, length_min: usize, length_max: usize) -> Self {
        self.length_min = length_min;
        self.length_max = length_max;
        self
    }

    /// Adds a forbidden sequence motif.
    pub fn forbid_motif(mut self, motif: impl Into<String>) -> Self {
        self.forbidden_motifs.push(motif.into());
        self
    }

    /// Adds a forbidden restriction-enzyme site by enzyme name.
    pub fn forbid_restriction_site(mut self, enzyme: impl Into<String>) -> Self {
        self.forbidden_restriction_sites.push(enzyme.into());
        self
    }

    /// Sets the expression host.
    pub fn with_host(mut self, host: ExpressionHost) -> Self {
        self.host = host;
        self
    }

    /// `true` if `length` falls inside the configured window (an
    /// unbounded side always passes).
    pub fn length_ok(&self, length: usize) -> bool {
        (self.length_min == 0 || length >= self.length_min)
            && (self.length_max == 0 || length <= self.length_max)
    }

    /// `true` if `gc` falls inside the configured GC window.
    pub fn gc_ok(&self, gc: f64) -> bool {
        gc >= self.gc_min - 1e-9 && gc <= self.gc_max + 1e-9
    }

    /// Validates the constraint set — a well-ordered GC range inside
    /// `[0, 1]`, a well-ordered length window, a positive homopolymer
    /// cap.
    ///
    /// # Errors
    /// [`RnaDesignError::Goal`] describing the first problem found.
    pub fn validate(&self) -> Result<()> {
        if !(0.0..=1.0).contains(&self.gc_min) || !(0.0..=1.0).contains(&self.gc_max) {
            return Err(RnaDesignError::goal(
                "gc_range",
                "GC fractions must lie in [0, 1]",
            ));
        }
        if self.gc_min > self.gc_max {
            return Err(RnaDesignError::goal(
                "gc_range",
                format!("gc_min {:.2} exceeds gc_max {:.2}", self.gc_min, self.gc_max),
            ));
        }
        if self.length_min != 0 && self.length_max != 0 && self.length_min > self.length_max {
            return Err(RnaDesignError::goal(
                "length",
                format!(
                    "length_min {} exceeds length_max {}",
                    self.length_min, self.length_max
                ),
            ));
        }
        if self.max_homopolymer == 0 {
            return Err(RnaDesignError::goal(
                "constraints",
                "max_homopolymer must be at least 1",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_goal_from_dot_bracket() {
        let g = DesignGoal::structural("((((....))))", StructuralClass::Hairpin).unwrap();
        assert!(g.validate().is_ok());
        assert_eq!(g.kind_name(), "structural RNA");
        assert!(!g.is_coding());
    }

    #[test]
    fn structural_goal_rejects_pseudoknot() {
        let err = DesignGoal::structural("((..[[..))..]]", StructuralClass::Custom).unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn structural_goal_rejects_garbage_dot_bracket() {
        assert!(DesignGoal::structural("((.x)", StructuralClass::Custom).is_err());
        assert!(DesignGoal::structural("", StructuralClass::Custom).is_err());
    }

    #[test]
    fn coding_goal_validates_protein() {
        let g = DesignGoal::coding(b"MKVLA".to_vec(), MrnaUseCase::Vaccine);
        assert!(g.validate().is_ok());
        assert!(g.is_coding());
    }

    #[test]
    fn coding_goal_rejects_bad_protein() {
        // 'X' is an ambiguity code, not a standard residue.
        let g = DesignGoal::coding(b"MKXLA".to_vec(), MrnaUseCase::Vaccine);
        assert!(g.validate().is_err());
        // Empty protein.
        let g = DesignGoal::coding(Vec::new(), MrnaUseCase::Vaccine);
        assert!(g.validate().is_err());
    }

    #[test]
    fn coding_goal_allows_trailing_stop() {
        let g = DesignGoal::coding(b"MKVL*".to_vec(), MrnaUseCase::Vaccine);
        assert!(g.validate().is_ok());
        // A stop in the middle is rejected.
        let g = DesignGoal::coding(b"MK*VL".to_vec(), MrnaUseCase::Vaccine);
        assert!(g.validate().is_err());
    }

    #[test]
    fn hybrid_goal_needs_elements() {
        let g = DesignGoal::Hybrid {
            protein: b"MKVL".to_vec(),
            use_case: MrnaUseCase::Vaccine,
            elements: Vec::new(),
        };
        assert!(g.validate().is_err());
        let g = DesignGoal::Hybrid {
            protein: b"MKVL".to_vec(),
            use_case: MrnaUseCase::Vaccine,
            elements: vec![StructuralElement {
                dot_bracket: "((((....))))".to_string(),
                placement: ElementPlacement::FivePrimeUtr,
                label: "5' hairpin".to_string(),
            }],
        };
        assert!(g.validate().is_ok());
    }

    #[test]
    fn constraints_default_is_valid() {
        assert!(DesignConstraints::default().validate().is_ok());
    }

    #[test]
    fn constraints_reject_inverted_gc_range() {
        let c = DesignConstraints::new().with_gc_range(0.7, 0.3);
        assert!(c.validate().is_err());
        let c = DesignConstraints::new().with_gc_range(-0.1, 0.5);
        assert!(c.validate().is_err());
    }

    #[test]
    fn constraints_reject_inverted_length_range() {
        let c = DesignConstraints::new().with_length_range(200, 100);
        assert!(c.validate().is_err());
    }

    #[test]
    fn length_and_gc_predicates() {
        let c = DesignConstraints::new()
            .with_gc_range(0.4, 0.6)
            .with_length_range(50, 150);
        assert!(c.length_ok(100));
        assert!(!c.length_ok(40));
        assert!(!c.length_ok(200));
        assert!(c.gc_ok(0.5));
        assert!(!c.gc_ok(0.7));
        // Unbounded length always passes.
        let unbounded = DesignConstraints::new();
        assert!(unbounded.length_ok(1));
        assert!(unbounded.length_ok(100_000));
    }

    #[test]
    fn builders_accumulate() {
        let c = DesignConstraints::new()
            .forbid_motif("AAUAAA")
            .forbid_restriction_site("EcoRI")
            .with_host(ExpressionHost::EColi);
        assert_eq!(c.forbidden_motifs, vec!["AAUAAA"]);
        assert_eq!(c.forbidden_restriction_sites, vec!["EcoRI"]);
        assert_eq!(c.host, ExpressionHost::EColi);
    }
}

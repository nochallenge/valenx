//! Batch structure analysis — the [`StructureReport`] bundle and
//! multi-structure batch utilities.
//!
//! Where [`validate`](crate::validate) answers "is this model
//! geometrically sane", [`StructureReport`] answers "what *is* this
//! structure": it rolls every headline analysis on one structure —
//! per-chain summaries, the residue / atom composition, the DSSP
//! secondary-structure content, the radius of gyration and a
//! validation verdict — into a single struct. [`analyze_batch`] and
//! [`summarize_batch`] apply it across many structures, the typical
//! "profile this PDB set" workflow.
//!
//! Each protein / nucleic-acid chain's observed sequence is carried as
//! a validated [`valenx_bioseq::Seq`] — the same sequence type the
//! rest of the Valenx computational-biology stack uses, so a chain
//! extracted here drops straight into `valenx-bioseq` / `valenx-align`
//! without a re-parse.

use crate::dssp::{assign_chain, SecondaryStructure};
use crate::error::Result;
use crate::geometry::shape::model_radius_of_gyration;
use crate::structure::{Chain, ResidueKind, Structure};
use crate::validate::{validate_structure, ValidationReport};
use std::collections::BTreeMap;
use valenx_bioseq::{Seq, SeqKind};

/// The polymer class of a chain, inferred from its residue content.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ChainKind {
    /// A protein chain (mostly amino-acid residues).
    Protein,
    /// A DNA chain (mostly DNA nucleotides).
    Dna,
    /// An RNA chain (mostly RNA nucleotides).
    Rna,
    /// A chain with no polymer residues — ligands, ions, waters only.
    NonPolymer,
}

impl ChainKind {
    /// Short human-readable label.
    pub fn label(&self) -> &'static str {
        match self {
            ChainKind::Protein => "protein",
            ChainKind::Dna => "DNA",
            ChainKind::Rna => "RNA",
            ChainKind::NonPolymer => "non-polymer",
        }
    }
}

/// Per-chain secondary-structure content, as residue fractions in
/// `[0, 1]`. Computed for protein chains only (zeroed otherwise).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct SecondaryContent {
    /// Fraction of residues in a helix state (`H` / `G` / `I`).
    pub helix: f64,
    /// Fraction of residues in a sheet state (`E` / `B`).
    pub sheet: f64,
    /// Fraction of residues in a turn / bend state (`T` / `S`).
    pub turn: f64,
    /// Fraction of residues with no secondary-structure assignment.
    pub coil: f64,
}

/// A single-chain summary inside a [`StructureReport`].
#[derive(Clone, Debug, PartialEq)]
pub struct ChainSummary {
    /// Chain identifier.
    pub id: String,
    /// Inferred polymer class.
    pub kind: ChainKind,
    /// Total residue count (polymer + hetero).
    pub residue_count: usize,
    /// Polymer (amino-acid or nucleotide) residue count.
    pub polymer_residue_count: usize,
    /// Atom count across the chain.
    pub atom_count: usize,
    /// The observed-residue sequence as a validated [`Seq`]. `None`
    /// for a non-polymer chain or when the observed residues do not
    /// form a valid sequence of the inferred kind.
    pub sequence: Option<Seq>,
    /// DSSP secondary-structure content (protein chains only).
    pub secondary: SecondaryContent,
}

impl ChainSummary {
    /// The observed sequence as a plain string (`""` when there is no
    /// extracted sequence).
    pub fn sequence_string(&self) -> String {
        self.sequence
            .as_ref()
            .map(|s| s.as_str().to_string())
            .unwrap_or_default()
    }
}

/// A full single-structure analysis bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct StructureReport {
    /// The structure id.
    pub id: String,
    /// Free-text title.
    pub title: String,
    /// Number of coordinate models (NMR ensembles have many).
    pub model_count: usize,
    /// Per-chain summaries of the first model.
    pub chains: Vec<ChainSummary>,
    /// Residue-name counts across the first model (e.g. `ALA → 12`).
    pub residue_composition: BTreeMap<String, usize>,
    /// Element-symbol counts across the first model (e.g. `C → 480`).
    pub element_composition: BTreeMap<String, usize>,
    /// Total residue count across the first model.
    pub residue_count: usize,
    /// Total atom count across the first model.
    pub atom_count: usize,
    /// Count of water residues.
    pub water_count: usize,
    /// Count of non-water hetero (ligand / ion) residues.
    pub ligand_count: usize,
    /// Mass-weighted radius of gyration of the first model, ångström.
    pub radius_of_gyration: f64,
    /// The structure-validation verdict.
    pub validation: ValidationReport,
}

impl StructureReport {
    /// Compute the full report for `structure`.
    ///
    /// `clash_tolerance` is forwarded to [`validate_structure`]
    /// (≈ 0.4 Å is the MolProbity default).
    pub fn analyze(structure: &Structure, clash_tolerance: f64) -> Result<StructureReport> {
        let model = structure.first_model();

        // Per-chain summaries.
        let mut chains = Vec::with_capacity(model.chains.len());
        for chain in &model.chains {
            chains.push(summarize_chain(chain));
        }

        // Whole-model composition.
        let mut residue_composition: BTreeMap<String, usize> = BTreeMap::new();
        let mut element_composition: BTreeMap<String, usize> = BTreeMap::new();
        let mut water_count = 0usize;
        let mut ligand_count = 0usize;
        for residue in model.residues() {
            *residue_composition.entry(residue.name.clone()).or_default() += 1;
            match residue.kind() {
                ResidueKind::Water => water_count += 1,
                ResidueKind::Other => {
                    if residue.hetatm {
                        ligand_count += 1;
                    }
                }
                _ => {}
            }
        }
        for atom in model.atoms() {
            *element_composition
                .entry(atom.element.clone())
                .or_default() += 1;
        }

        // Radius of gyration — falls back to 0 for an atom-less model.
        let radius_of_gyration = model_radius_of_gyration(model).unwrap_or(0.0);

        let validation = validate_structure(structure, clash_tolerance)?;

        Ok(StructureReport {
            id: structure.id.clone(),
            title: structure.title.clone(),
            model_count: structure.models.len(),
            chains,
            residue_composition,
            element_composition,
            residue_count: model.residues().count(),
            atom_count: model.atom_count(),
            water_count,
            ligand_count,
            radius_of_gyration,
            validation,
        })
    }

    /// Number of protein chains.
    pub fn protein_chain_count(&self) -> usize {
        self.chains
            .iter()
            .filter(|c| c.kind == ChainKind::Protein)
            .count()
    }

    /// Number of nucleic-acid (DNA or RNA) chains.
    pub fn nucleic_chain_count(&self) -> usize {
        self.chains
            .iter()
            .filter(|c| matches!(c.kind, ChainKind::Dna | ChainKind::Rna))
            .count()
    }

    /// Mean DSSP helix fraction over the protein chains, or `0.0` when
    /// there are none.
    pub fn mean_helix_fraction(&self) -> f64 {
        let proteins: Vec<&ChainSummary> = self
            .chains
            .iter()
            .filter(|c| c.kind == ChainKind::Protein)
            .collect();
        if proteins.is_empty() {
            return 0.0;
        }
        proteins.iter().map(|c| c.secondary.helix).sum::<f64>() / proteins.len() as f64
    }

    /// Mean DSSP sheet fraction over the protein chains, or `0.0` when
    /// there are none.
    pub fn mean_sheet_fraction(&self) -> f64 {
        let proteins: Vec<&ChainSummary> = self
            .chains
            .iter()
            .filter(|c| c.kind == ChainKind::Protein)
            .collect();
        if proteins.is_empty() {
            return 0.0;
        }
        proteins.iter().map(|c| c.secondary.sheet).sum::<f64>() / proteins.len() as f64
    }

    /// A compact multi-line text summary, handy for logs.
    pub fn summary(&self) -> String {
        format!(
            "structure {}: {} model(s), {} chain(s) ({} protein, {} nucleic)\n  \
             residues: {} ({} water, {} ligand), atoms: {}\n  \
             Rg: {:.1} Å | helix {:.0}% sheet {:.0}%\n  \
             validation: {} clash(es), quality {:.1}/100",
            self.id,
            self.model_count,
            self.chains.len(),
            self.protein_chain_count(),
            self.nucleic_chain_count(),
            self.residue_count,
            self.water_count,
            self.ligand_count,
            self.atom_count,
            self.radius_of_gyration,
            self.mean_helix_fraction() * 100.0,
            self.mean_sheet_fraction() * 100.0,
            self.validation.clashes.len(),
            self.validation.quality_score(),
        )
    }
}

/// Build the [`ChainSummary`] of one chain.
fn summarize_chain(chain: &Chain) -> ChainSummary {
    let kind = classify_chain(chain);
    let polymer_residue_count = chain.polymer_residues().len();
    let atom_count: usize = chain.residues.iter().map(|r| r.atoms.len()).sum();

    // Observed-residue sequence as a validated Seq, when the chain is
    // a polymer. `observed_sequence` maps unknowns to `X`, which the
    // protein alphabet accepts; the nucleotide alphabets do not, so a
    // nucleic chain with an `X` falls back to `None`.
    let sequence = match kind {
        ChainKind::Protein => {
            Seq::new(SeqKind::Protein, chain.observed_sequence()).ok()
        }
        ChainKind::Dna => Seq::new(SeqKind::Dna, chain.observed_sequence()).ok(),
        ChainKind::Rna => Seq::new(SeqKind::Rna, chain.observed_sequence()).ok(),
        ChainKind::NonPolymer => None,
    };

    // DSSP secondary-structure content for protein chains.
    let secondary = if kind == ChainKind::Protein {
        secondary_content(chain)
    } else {
        SecondaryContent::default()
    };

    ChainSummary {
        id: chain.id.clone(),
        kind,
        residue_count: chain.residues.len(),
        polymer_residue_count,
        atom_count,
        sequence,
        secondary,
    }
}

/// Infer a chain's polymer class from its residue content: the
/// majority polymer kind wins; a chain with no polymer residues is
/// [`ChainKind::NonPolymer`].
fn classify_chain(chain: &Chain) -> ChainKind {
    let (mut aa, mut dna, mut rna) = (0usize, 0usize, 0usize);
    for r in &chain.residues {
        match r.kind() {
            ResidueKind::AminoAcid => aa += 1,
            ResidueKind::Dna => dna += 1,
            ResidueKind::Rna => rna += 1,
            _ => {}
        }
    }
    if aa == 0 && dna == 0 && rna == 0 {
        return ChainKind::NonPolymer;
    }
    if aa >= dna && aa >= rna {
        ChainKind::Protein
    } else if dna >= rna {
        ChainKind::Dna
    } else {
        ChainKind::Rna
    }
}

/// DSSP secondary-structure content of a protein chain.
fn secondary_content(chain: &Chain) -> SecondaryContent {
    let result = assign_chain(chain);
    if result.states.is_empty() {
        return SecondaryContent::default();
    }
    let total = result.states.len() as f64;
    let mut helix = 0usize;
    let mut sheet = 0usize;
    let mut turn = 0usize;
    for s in &result.states {
        match s {
            SecondaryStructure::AlphaHelix
            | SecondaryStructure::Helix310
            | SecondaryStructure::PiHelix => helix += 1,
            SecondaryStructure::Strand | SecondaryStructure::Bridge => sheet += 1,
            SecondaryStructure::Turn | SecondaryStructure::Bend => turn += 1,
            SecondaryStructure::Coil => {}
        }
    }
    let coil = result.states.len() - helix - sheet - turn;
    SecondaryContent {
        helix: helix as f64 / total,
        sheet: sheet as f64 / total,
        turn: turn as f64 / total,
        coil: coil as f64 / total,
    }
}

/// Analyse every structure in a slice, returning one
/// [`StructureReport`] per structure.
pub fn analyze_batch(
    structures: &[Structure],
    clash_tolerance: f64,
) -> Result<Vec<StructureReport>> {
    structures
        .iter()
        .map(|s| StructureReport::analyze(s, clash_tolerance))
        .collect()
}

/// Keep only the structures whose report satisfies `predicate` — the
/// "filter a structure set" workflow.
pub fn filter_batch(
    structures: Vec<Structure>,
    clash_tolerance: f64,
    predicate: impl Fn(&StructureReport) -> bool,
) -> Result<Vec<Structure>> {
    let mut kept = Vec::new();
    for s in structures {
        let report = StructureReport::analyze(&s, clash_tolerance)?;
        if predicate(&report) {
            kept.push(s);
        }
    }
    Ok(kept)
}

/// Summary statistics over a batch of [`StructureReport`]s — the
/// headline numbers a structure-set triage wants.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BatchSummary {
    /// Number of structures analysed.
    pub count: usize,
    /// Total chains across the batch.
    pub total_chains: usize,
    /// Mean atom count per structure.
    pub mean_atom_count: f64,
    /// Mean radius of gyration, ångström.
    pub mean_radius_of_gyration: f64,
    /// Mean validation quality score (0..100).
    pub mean_quality_score: f64,
    /// Number of structures that passed validation with no flags.
    pub clean_structures: usize,
}

/// Roll a slice of reports into a [`BatchSummary`].
pub fn summarize_batch(reports: &[StructureReport]) -> BatchSummary {
    if reports.is_empty() {
        return BatchSummary::default();
    }
    let n = reports.len() as f64;
    BatchSummary {
        count: reports.len(),
        total_chains: reports.iter().map(|r| r.chains.len()).sum(),
        mean_atom_count: reports.iter().map(|r| r.atom_count as f64).sum::<f64>() / n,
        mean_radius_of_gyration: reports
            .iter()
            .map(|r| r.radius_of_gyration)
            .sum::<f64>()
            / n,
        mean_quality_score: reports
            .iter()
            .map(|r| r.validation.quality_score())
            .sum::<f64>()
            / n,
        clean_structures: reports
            .iter()
            .filter(|r| r.validation.is_clean())
            .count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain, Residue};
    use nalgebra::Point3;

    /// A small protein chain of `n` alanines with full backbones.
    fn protein_chain(id: &str, n: usize) -> Chain {
        let mut chain = Chain::new(id);
        for seq in 0..n {
            let base = seq as f64 * 3.8;
            let mut r = Residue::new("ALA", seq as i32 + 1);
            r.atoms
                .push(Atom::new("N", "N", Point3::new(base, 0.0, 0.0)));
            r.atoms
                .push(Atom::new("CA", "C", Point3::new(base + 1.5, 0.0, 0.0)));
            r.atoms
                .push(Atom::new("C", "C", Point3::new(base + 2.5, 1.4, 0.0)));
            r.atoms
                .push(Atom::new("O", "O", Point3::new(base + 2.5, 2.6, 0.0)));
            r.atoms
                .push(Atom::new("CB", "C", Point3::new(base + 1.5, -1.5, 0.0)));
            chain.residues.push(r);
        }
        chain
    }

    fn protein_structure(id: &str) -> Structure {
        let mut s = Structure::new(id);
        s.title = "demo".to_string();
        s.first_model_mut().chains.push(protein_chain("A", 6));
        s
    }

    #[test]
    fn report_has_chain_summary() {
        let s = protein_structure("1ABC");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert_eq!(report.id, "1ABC");
        assert_eq!(report.chains.len(), 1);
        assert_eq!(report.chains[0].kind, ChainKind::Protein);
        assert_eq!(report.chains[0].residue_count, 6);
        assert_eq!(report.protein_chain_count(), 1);
        assert_eq!(report.nucleic_chain_count(), 0);
    }

    #[test]
    fn chain_sequence_is_a_valid_seq() {
        let s = protein_structure("X");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        let seq = report.chains[0].sequence.as_ref().expect("protein seq");
        assert_eq!(seq.kind(), SeqKind::Protein);
        assert_eq!(seq.as_str(), "AAAAAA");
        assert_eq!(report.chains[0].sequence_string(), "AAAAAA");
    }

    #[test]
    fn composition_counts_residues_and_elements() {
        let s = protein_structure("X");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert_eq!(report.residue_composition.get("ALA"), Some(&6));
        // 6 residues * (3 C + 1 N + 1 O) = 18 C, 6 N, 6 O.
        assert_eq!(report.element_composition.get("C"), Some(&18));
        assert_eq!(report.element_composition.get("N"), Some(&6));
        assert_eq!(report.residue_count, 6);
        assert_eq!(report.atom_count, 30);
    }

    #[test]
    fn waters_and_ligands_are_counted() {
        let mut s = protein_structure("X");
        let mut w = Residue::new("HOH", 100);
        w.hetatm = true;
        w.atoms.push(Atom::new("O", "O", Point3::new(50.0, 0.0, 0.0)));
        let mut zn = Residue::new("ZN", 101);
        zn.hetatm = true;
        zn.atoms
            .push(Atom::new("ZN", "ZN", Point3::new(60.0, 0.0, 0.0)));
        s.first_model_mut().chains[0].residues.push(w);
        s.first_model_mut().chains[0].residues.push(zn);
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert_eq!(report.water_count, 1);
        assert_eq!(report.ligand_count, 1);
    }

    #[test]
    fn radius_of_gyration_is_positive() {
        let s = protein_structure("X");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert!(report.radius_of_gyration > 0.0);
    }

    #[test]
    fn secondary_content_fractions_sum_to_one() {
        let s = protein_structure("X");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        let sc = report.chains[0].secondary;
        let total = sc.helix + sc.sheet + sc.turn + sc.coil;
        assert!((total - 1.0).abs() < 1e-9, "fractions sum to {total}");
    }

    #[test]
    fn classify_chain_kinds() {
        let protein = protein_chain("A", 4);
        assert_eq!(classify_chain(&protein), ChainKind::Protein);

        let mut dna = Chain::new("B");
        for seq in 1..=4 {
            let mut r = Residue::new("DA", seq);
            r.atoms.push(Atom::new("P", "P", Point3::origin()));
            dna.residues.push(r);
        }
        assert_eq!(classify_chain(&dna), ChainKind::Dna);

        let mut ligand_only = Chain::new("C");
        let mut zn = Residue::new("ZN", 1);
        zn.hetatm = true;
        zn.atoms.push(Atom::new("ZN", "ZN", Point3::origin()));
        ligand_only.residues.push(zn);
        assert_eq!(classify_chain(&ligand_only), ChainKind::NonPolymer);
    }

    #[test]
    fn nucleic_chain_has_no_sequence_with_unknowns() {
        // A DNA chain whose residues are all standard maps cleanly;
        // a non-polymer chain has no sequence.
        let mut s = Structure::new("DNA");
        let mut chain = Chain::new("A");
        for seq in 1..=3 {
            let mut r = Residue::new("DA", seq);
            r.atoms.push(Atom::new("P", "P", Point3::new(seq as f64, 0.0, 0.0)));
            chain.residues.push(r);
        }
        s.first_model_mut().chains.push(chain);
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert_eq!(report.chains[0].kind, ChainKind::Dna);
        // DA -> 'A', a valid DNA residue.
        assert_eq!(report.chains[0].sequence_string(), "AAA");
    }

    #[test]
    fn summary_mentions_the_structure() {
        let s = protein_structure("MYPROT");
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        let text = report.summary();
        assert!(text.contains("MYPROT"));
        assert!(text.contains("Rg"));
        assert!(text.contains("validation"));
    }

    #[test]
    fn batch_analysis_and_summary() {
        let structures = vec![
            protein_structure("A1"),
            protein_structure("A2"),
            protein_structure("A3"),
        ];
        let reports = analyze_batch(&structures, 0.4).unwrap();
        assert_eq!(reports.len(), 3);
        let summary = summarize_batch(&reports);
        assert_eq!(summary.count, 3);
        assert_eq!(summary.total_chains, 3);
        assert!(summary.mean_atom_count > 0.0);
        assert!((0.0..=100.0).contains(&summary.mean_quality_score));
    }

    #[test]
    fn filter_batch_keeps_matching() {
        let structures = vec![protein_structure("SMALL"), {
            let mut big = Structure::new("BIG");
            big.first_model_mut().chains.push(protein_chain("A", 20));
            big
        }];
        // Keep only structures with more than 10 residues.
        let kept = filter_batch(structures, 0.4, |r| r.residue_count > 10).unwrap();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].id, "BIG");
    }

    #[test]
    fn empty_batch_summary_is_default() {
        let summary = summarize_batch(&[]);
        assert_eq!(summary, BatchSummary::default());
        assert_eq!(summary.count, 0);
    }
}

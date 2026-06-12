//! Curated reference-gene lists for the commercial-depth safety screen.
//!
//! Production gene-editing design tools (Benchling, Synthego, IDT,
//! CRISPOR) cross-reference every cut site against three well-known
//! community catalogues to grade safety:
//!
//! - **Essential genes** — genes whose loss-of-function is lethal to
//!   the cell (or strongly fitness-impairing). Cutting an essential
//!   gene with an NHEJ-prone nuclease is the classical "kill the cell
//!   you wanted to engineer" failure mode. Source convention: the
//!   `OGEE` and `DepMap` essential-gene panels.
//! - **Cancer-driver genes** — genes whose perturbation has been
//!   causally linked to cancer transformation. Cutting them off-target
//!   raises an oncogenicity flag. Source convention: the Sanger
//!   Cancer Gene Census tier-1 list.
//! - **Safe-harbor loci** — well-characterised genomic locations
//!   (AAVS1, hROSA26, CCR5) where stable integration is tolerated
//!   without disrupting essential or oncogenic function. An edit at a
//!   safe-harbor locus is *low* risk and should be reported as such.
//!
//! ## v1 scope
//!
//! Every list here is a **representative published subset** —
//! ~100-150 entries of the ~1500-2000 the full curated databases
//! contain. The lists are encoded inline; the safety aggregator
//! consumes them and reports the per-edit risk verdict. The names are
//! drawn from the standard community catalogues (HGNC symbols), but
//! the full catalogues live behind versioned web databases and are
//! the caller's job to plug in for a production install. The on-disk
//! gene-symbol comparison is case-insensitive; the lists carry no
//! genomic coordinates because the catalogues' coordinates are
//! reference-build-dependent and the safety screen here works at the
//! **symbol level** — exactly the granularity an in-silico design
//! tool needs.

use serde::{Deserialize, Serialize};

/// A reference gene-symbol list with a kind label.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReferenceGeneList {
    /// What kind of list this is.
    pub kind: ReferenceListKind,
    /// HGNC-style gene symbols in the list, uppercased and sorted.
    pub symbols: Vec<String>,
}

/// Tag for what a [`ReferenceGeneList`] represents.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ReferenceListKind {
    /// Essential genes whose disruption is fitness-lethal.
    Essential,
    /// Cancer-driver / oncogene-or-tumour-suppressor genes.
    CancerDriver,
    /// Safe-harbor loci (low integration risk).
    SafeHarbor,
}

impl ReferenceGeneList {
    /// `true` when `symbol` (case-insensitive) is in this list.
    pub fn contains(&self, symbol: &str) -> bool {
        let u = symbol.to_ascii_uppercase();
        self.symbols.binary_search(&u).is_ok()
    }

    /// Number of symbols in this list.
    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    /// `true` when empty (never produced; clippy hygiene).
    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }
}

/// The curated essential-gene list (representative subset, ~110
/// symbols spanning the major essential-gene families: ribosomal
/// proteins, RNA polymerases, splicing, replication, the proteasome,
/// the spliceosome, tRNA synthetases, the cell cycle).
///
/// Sources of the qualitative panel: OGEE v3 (Chen et al. 2020) /
/// DepMap (Tsherniak et al. 2017) / MaGeCK common essentials. Symbols
/// are HGNC; the subset is curated for *recognisability* (every entry
/// is one of the textbook essential-gene families) and for *coverage*
/// of the safety-screen smoke test (any cut near a ribosomal protein
/// should raise a warning).
pub fn essential_genes() -> ReferenceGeneList {
    let symbols: &[&str] = &[
        // RNA Pol II + transcription core.
        "POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E", "POLR2F", "POLR2G", "POLR2H", "POLR2I",
        "POLR2J", "POLR2K", "POLR2L", // DNA polymerases + replication.
        "POLA1", "POLA2", "POLD1", "POLD2", "POLE", "POLE2", "PCNA", "RFC1", "RFC2", "RFC3",
        "RFC4", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
        // Ribosomal proteins — large subunit (RPL).
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL9", "RPL10", "RPL11", "RPL12", "RPL13",
        "RPL14", "RPL15", "RPL17", "RPL18", "RPL19", "RPL21", "RPL22", "RPL23", "RPL24", "RPL26",
        "RPL27", "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL36", "RPL37",
        "RPL38", "RPL39", "RPLP0", "RPLP1", "RPLP2",
        // Ribosomal proteins — small subunit (RPS).
        "RPS2", "RPS3", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12",
        "RPS13", "RPS14", "RPS15", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23",
        "RPS24", "RPS25", "RPS26", "RPS27", "RPS28", "RPS29", "RPSA",
        // Spliceosome core.
        "PRPF8", "PRPF19", "PRPF31", "PRPF38A", "SNRNP200", "SNRPA1", "SNRPB", "SNRPD1", "SNRPD2",
        "SNRPD3", "SNRPE", "SNRPF", "SNRPG", // tRNA synthetases (a sample).
        "AARS1", "DARS1", "EARS2", "FARSA", "GARS1", "HARS1", "IARS1", "KARS1", "LARS1",
        // Proteasome core.
        "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2",
        // Cell cycle / housekeeping master regulators that are
        // unambiguously essential in human cells.
        "CDK1", "CDK7", "CDK9", "MYC", "CCNB1", "CCNH",
    ];
    let mut v: Vec<String> = symbols.iter().map(|s| s.to_ascii_uppercase()).collect();
    v.sort();
    v.dedup();
    ReferenceGeneList {
        kind: ReferenceListKind::Essential,
        symbols: v,
    }
}

/// The curated cancer-driver list (representative subset, ~110
/// symbols spanning the major Sanger Cancer Gene Census tier-1
/// driver families: receptor tyrosine kinases, RAS-MAPK, PI3K-AKT,
/// the chromatin / transcription factor panel, tumour suppressors,
/// and the well-studied fusion partners).
///
/// Source convention: Sanger Cancer Gene Census tier-1
/// (Sondka et al. 2018 / 2024 updates), one of the most-cited
/// oncogene reference catalogues. The subset is curated for
/// *recognisability* (every entry is a textbook cancer driver) and
/// for *coverage* of the smoke-test set (TP53, KRAS, MYC, EGFR all
/// included so the safety screen unambiguously fires on them).
pub fn cancer_driver_genes() -> ReferenceGeneList {
    let symbols: &[&str] = &[
        // Receptor tyrosine kinases.
        "EGFR", "ERBB2", "ERBB3", "ERBB4", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "MET", "ALK", "ROS1",
        "RET", "KIT", "PDGFRA", "PDGFRB", // RAS-MAPK pathway.
        "KRAS", "NRAS", "HRAS", "BRAF", "MAP2K1", "MAP2K2", "NF1", // PI3K-AKT-mTOR.
        "PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "AKT3", "MTOR", "PTEN", "TSC1", "TSC2",
        "STK11", // Cell cycle + DNA damage response.
        "TP53", "RB1", "CDKN2A", "CDK4", "CDK6", "CCND1", "CCNE1", "MDM2", "MDM4", "ATM", "ATR",
        "CHEK2", "BRCA1", "BRCA2", "PALB2", // Transcription factors + chromatin.
        "MYC", "MYCN", "MYCL", "EZH2", "KDM6A", "KMT2A", "KMT2C", "KMT2D", "CREBBP", "EP300",
        "ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", // Wnt + Hedgehog + Notch.
        "APC", "CTNNB1", "AXIN1", "AXIN2", "PTCH1", "SMO", "GLI1", "NOTCH1", "NOTCH2", "FBXW7",
        // Translation initiation + ribosome biogenesis.
        "DICER1", "DROSHA", // Apoptosis + survival.
        "BCL2", "BCL6", "BCL10", "MCL1", "BIRC3", "CASP8",
        // Hematological-malignancy core.
        "FLT3", "NPM1", "DNMT3A", "TET2", "IDH1", "IDH2", "WT1", "RUNX1", "JAK2", "JAK3", "CALR",
        "MPL", "CBL", "ASXL1", "SF3B1", "SRSF2", "U2AF1", "ETV6", "GATA2", "PAX5", "IKZF1",
        // Mismatch repair + replication-fidelity tumour suppressors.
        "MLH1", "MSH2", "MSH6", "PMS2", // VHL + renal/pheo drivers.
        "VHL", "SDHA", "SDHB", "SDHC", "SDHD",
    ];
    let mut v: Vec<String> = symbols.iter().map(|s| s.to_ascii_uppercase()).collect();
    v.sort();
    v.dedup();
    ReferenceGeneList {
        kind: ReferenceListKind::CancerDriver,
        symbols: v,
    }
}

/// The curated safe-harbor loci.
///
/// "Safe harbor" is a well-characterised location where stable
/// transgene integration is well tolerated — no disruption of a
/// nearby essential or oncogenic gene, predictable expression. The
/// three textbook human safe harbors are:
///
/// - **AAVS1 / PPP1R12C intron 1** — the classical adeno-associated-
///   virus integration site on chr19.
/// - **hROSA26 / THUMPD3** — the human ortholog of the mouse ROSA26
///   safe harbor.
/// - **CCR5** — disruption-tolerated, the CCR5-Δ32 natural allele is
///   common and benign for non-HIV-related function.
///
/// Source convention: Sadelain et al. 2012, Papapetrou & Schambach
/// 2016 / IDT / Benchling safe-harbor panels.
pub fn safe_harbor_loci() -> ReferenceGeneList {
    let symbols: &[&str] = &[
        "AAVS1", "PPP1R12C", // host gene for AAVS1
        "HROSA26", "ROSA26", "THUMPD3", "CCR5",
    ];
    let mut v: Vec<String> = symbols.iter().map(|s| s.to_ascii_uppercase()).collect();
    v.sort();
    v.dedup();
    ReferenceGeneList {
        kind: ReferenceListKind::SafeHarbor,
        symbols: v,
    }
}

/// A bundle of every curated list — the convenience entry point for
/// the safety aggregator.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReferenceGeneDatabase {
    /// The essential-gene list.
    pub essential: ReferenceGeneList,
    /// The cancer-driver list.
    pub cancer_drivers: ReferenceGeneList,
    /// The safe-harbor list.
    pub safe_harbors: ReferenceGeneList,
}

impl ReferenceGeneDatabase {
    /// Builds the default curated database — every list inline,
    /// always available without I/O.
    pub fn curated() -> Self {
        ReferenceGeneDatabase {
            essential: essential_genes(),
            cancer_drivers: cancer_driver_genes(),
            safe_harbors: safe_harbor_loci(),
        }
    }

    /// `true` when `symbol` is an essential gene (case-insensitive).
    pub fn is_essential(&self, symbol: &str) -> bool {
        self.essential.contains(symbol)
    }

    /// `true` when `symbol` is a cancer-driver gene.
    pub fn is_cancer_driver(&self, symbol: &str) -> bool {
        self.cancer_drivers.contains(symbol)
    }

    /// `true` when `symbol` is a safe-harbor locus.
    pub fn is_safe_harbor(&self, symbol: &str) -> bool {
        self.safe_harbors.contains(symbol)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn essential_list_includes_ribosomal_genes() {
        let g = essential_genes();
        assert!(g.contains("RPS6"));
        assert!(g.contains("rpl10")); // case-insensitive
        assert!(g.contains("POLR2A"));
        assert!(!g.contains("DEFINITELY_NOT_A_GENE"));
    }

    #[test]
    fn cancer_driver_list_includes_textbook_drivers() {
        let g = cancer_driver_genes();
        assert!(g.contains("TP53"));
        assert!(g.contains("KRAS"));
        assert!(g.contains("MYC"));
        assert!(g.contains("EGFR"));
        assert!(g.contains("BRCA1"));
    }

    #[test]
    fn safe_harbor_includes_the_three_classics() {
        let g = safe_harbor_loci();
        assert!(g.contains("AAVS1"));
        assert!(g.contains("CCR5"));
        // hROSA26 + its alias both present.
        assert!(g.contains("HROSA26"));
        assert!(g.contains("ROSA26"));
    }

    #[test]
    fn lists_are_sorted_for_binary_search() {
        for g in [essential_genes(), cancer_driver_genes(), safe_harbor_loci()] {
            let mut sorted = g.symbols.clone();
            sorted.sort();
            assert_eq!(g.symbols, sorted, "{:?} list not sorted", g.kind);
            // No duplicates.
            let mut dedup = sorted.clone();
            dedup.dedup();
            assert_eq!(
                dedup.len(),
                sorted.len(),
                "{:?} list has duplicates",
                g.kind
            );
        }
    }

    #[test]
    fn database_convenience_accessors() {
        let db = ReferenceGeneDatabase::curated();
        assert!(db.is_essential("RPS6"));
        assert!(db.is_cancer_driver("TP53"));
        assert!(db.is_safe_harbor("AAVS1"));
        // A gene that is a driver but not essential nor a harbor.
        assert!(!db.is_essential("TP53"));
        assert!(!db.is_safe_harbor("TP53"));
    }

    #[test]
    fn lists_have_meaningful_size() {
        // Sanity floor: enough to give the safety screen real coverage
        // (the smoke-test expectations target ribosomal / cancer
        // panels of ~100 symbols).
        assert!(essential_genes().len() >= 80);
        assert!(cancer_driver_genes().len() >= 80);
        assert!(safe_harbor_loci().len() >= 3);
    }
}

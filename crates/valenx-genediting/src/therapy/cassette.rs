//! Feature 24 — gene-therapy expression-cassette design.
//!
//! A gene therapy delivers a *transgene* packaged inside a viral
//! vector. The DNA payload — the **expression cassette** — has, 5′→3′:
//!
//! - a **promoter** (and optionally an enhancer) driving transcription;
//! - the **transgene** coding sequence;
//! - a **polyadenylation signal**;
//! - optional **regulatory elements** — a WPRE to boost expression,
//!   ITRs / LTRs that the vector requires.
//!
//! The hard constraint is **payload size**. A recombinant AAV packages
//! only ~4.7 kb between its ITRs; a lentivirus tolerates a much larger
//! ~8–9 kb between its LTRs. A cassette that overruns its vector
//! simply will not package.
//!
//! This module assembles a cassette layout, sums its size, and checks
//! it against the chosen vector's capacity ([`design_cassette`]).
//!
//! ## v1 scope
//!
//! The cassette designer works with *sized, named* elements — it
//! tracks lengths and the part order and validates the transgene CDS;
//! it does not assemble the literal vector backbone sequence. The
//! capacity limits are the well-established consensus figures (AAV
//! ~4.7 kb, lentivirus ~8.5 kb usable).

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgt, transcribe};
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;

/// A viral delivery vector for a gene-therapy cassette.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DeliveryVector {
    /// Recombinant adeno-associated virus — small, non-integrating,
    /// the dominant in-vivo vector. ~4.7 kb ITR-to-ITR capacity.
    Aav,
    /// Lentivirus — integrating, large cargo, the workhorse for
    /// ex-vivo cell engineering. ~8.5 kb usable capacity.
    Lentivirus,
    /// Adenovirus — large cargo, episomal, strongly immunogenic.
    /// ~8 kb in a first-generation vector (much more gutless).
    Adenovirus,
}

impl DeliveryVector {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            DeliveryVector::Aav => "AAV",
            DeliveryVector::Lentivirus => "lentivirus",
            DeliveryVector::Adenovirus => "adenovirus",
        }
    }

    /// The usable packaging capacity of the vector in base pairs — the
    /// space available for the expression cassette (the consensus
    /// figure for the cassette payload, not the whole genome).
    pub fn capacity_bp(self) -> usize {
        match self {
            DeliveryVector::Aav => 4700,
            DeliveryVector::Lentivirus => 8500,
            DeliveryVector::Adenovirus => 8000,
        }
    }

    /// `true` when the vector integrates into the host genome
    /// (lentivirus). Integration carries an insertional-mutagenesis
    /// consideration the safety module surfaces.
    pub fn integrates(self) -> bool {
        matches!(self, DeliveryVector::Lentivirus)
    }
}

/// A promoter choice for the cassette.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Promoter {
    /// CMV — a strong, ubiquitous viral promoter; can silence over
    /// time in vivo.
    Cmv,
    /// EF-1α — a strong, ubiquitous mammalian promoter; resists
    /// silencing better than CMV.
    Ef1Alpha,
    /// CAG — a strong synthetic hybrid (CMV enhancer + chicken
    /// β-actin promoter).
    Cag,
    /// A small synthetic / tissue-specific promoter — chosen when the
    /// cassette is size-constrained or expression must be restricted.
    SyntheticSmall,
}

impl Promoter {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            Promoter::Cmv => "CMV",
            Promoter::Ef1Alpha => "EF-1alpha",
            Promoter::Cag => "CAG",
            Promoter::SyntheticSmall => "synthetic small promoter",
        }
    }

    /// A representative length of the promoter element in base pairs.
    pub fn length_bp(self) -> usize {
        match self {
            Promoter::Cmv => 600,
            Promoter::Ef1Alpha => 1200,
            Promoter::Cag => 1700,
            Promoter::SyntheticSmall => 250,
        }
    }
}

/// A polyadenylation-signal choice for the cassette.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PolyASignal {
    /// SV40 late polyA — compact (~130 bp).
    Sv40,
    /// Bovine growth-hormone polyA — strong (~225 bp).
    Bgh,
    /// A short synthetic polyA signal (~50 bp) for tight size budgets.
    SyntheticShort,
}

impl PolyASignal {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            PolyASignal::Sv40 => "SV40 late polyA",
            PolyASignal::Bgh => "bGH polyA",
            PolyASignal::SyntheticShort => "synthetic short polyA",
        }
    }

    /// A representative length in base pairs.
    pub fn length_bp(self) -> usize {
        match self {
            PolyASignal::Sv40 => 130,
            PolyASignal::Bgh => 225,
            PolyASignal::SyntheticShort => 50,
        }
    }
}

/// A request to design a gene-therapy expression cassette.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CassetteRequest {
    /// The transgene coding sequence (DNA or RNA; `ATG`…stop, length
    /// divisible by three).
    pub transgene: Vec<u8>,
    /// The delivery vector — sets the capacity limit.
    pub vector: DeliveryVector,
    /// The promoter.
    pub promoter: Promoter,
    /// The polyA signal.
    pub poly_a: PolyASignal,
    /// Include a WPRE (woodchuck post-transcriptional regulatory
    /// element, ~600 bp) to boost expression. Costs payload space.
    pub include_wpre: bool,
}

impl CassetteRequest {
    /// A request with sensible defaults: EF-1α promoter, bGH polyA, no
    /// WPRE.
    pub fn new(transgene: impl Into<Vec<u8>>, vector: DeliveryVector) -> Self {
        CassetteRequest {
            transgene: transgene.into(),
            vector,
            promoter: Promoter::Ef1Alpha,
            poly_a: PolyASignal::Bgh,
            include_wpre: false,
        }
    }
}

/// One sized element of an assembled cassette.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CassetteElement {
    /// Element name (`"promoter (EF-1alpha)"`, `"transgene"`, …).
    pub name: String,
    /// Element length in base pairs.
    pub length_bp: usize,
}

/// A designed gene-therapy expression cassette.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ExpressionCassette {
    /// The cassette elements, 5′→3′.
    pub elements: Vec<CassetteElement>,
    /// Total cassette size in base pairs.
    pub total_bp: usize,
    /// The delivery vector and its usable capacity.
    pub vector: DeliveryVector,
    /// `true` when the cassette fits inside the vector's capacity.
    pub fits_vector: bool,
    /// Spare capacity in base pairs (`0` when the cassette overruns).
    pub headroom_bp: usize,
    /// Human-readable notes — fit status, recommendations.
    pub notes: Vec<String>,
}

impl ExpressionCassette {
    /// The transgene length in base pairs (the `"transgene"` element).
    pub fn transgene_bp(&self) -> usize {
        self.elements
            .iter()
            .find(|e| e.name == "transgene")
            .map(|e| e.length_bp)
            .unwrap_or(0)
    }
}

/// The representative size of the AAV ITR pair / lentiviral LTR pair
/// the cassette must sit between (counted against the payload).
fn vector_overhead_bp(vector: DeliveryVector) -> (usize, &'static str) {
    match vector {
        DeliveryVector::Aav => (290, "AAV ITRs (2 x ~145 bp)"),
        DeliveryVector::Lentivirus => (1400, "lentiviral LTRs + psi packaging signal"),
        DeliveryVector::Adenovirus => (200, "adenoviral ITRs + packaging signal"),
    }
}

/// The WPRE element length in base pairs.
const WPRE_BP: usize = 600;

/// Designs a gene-therapy expression cassette (feature 24).
///
/// Validates the transgene as a coding sequence, assembles the cassette
/// (vector overhead + promoter + transgene + polyA + optional WPRE),
/// sums the size and checks it against the vector's packaging capacity.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] if the transgene is empty, not a
/// multiple of three, non-ACGT/U, or not a proper CDS.
pub fn design_cassette(req: &CassetteRequest) -> Result<ExpressionCassette> {
    // Validate the transgene CDS.
    let rna = transcribe(&req.transgene);
    if rna.is_empty() || rna.len() % 3 != 0 {
        return Err(GeneditingError::invalid_target(
            "cds",
            "transgene must be a non-empty CDS of length divisible by 3",
        ));
    }
    // The transgene may be given as DNA or RNA; validate as DNA.
    let dna: Vec<u8> = rna
        .iter()
        .map(|&b| if b == b'U' { b'T' } else { b })
        .collect();
    if !is_acgt(&dna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "transgene must be A/C/G/T(/U)",
        ));
    }
    let code = GeneticCode::standard();
    if !code.is_start_codon(&[dna[0], dna[1], dna[2]]) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "transgene does not begin with a start codon",
        ));
    }
    let n = dna.len();
    if !code.is_stop_codon(&[dna[n - 3], dna[n - 2], dna[n - 1]]) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "transgene does not end with a stop codon",
        ));
    }

    // Assemble the cassette elements.
    let (overhead, overhead_name) = vector_overhead_bp(req.vector);
    let mut elements = vec![
        CassetteElement {
            name: format!("vector overhead ({overhead_name})"),
            length_bp: overhead,
        },
        CassetteElement {
            name: format!("promoter ({})", req.promoter.name()),
            length_bp: req.promoter.length_bp(),
        },
        CassetteElement {
            name: "transgene".to_string(),
            length_bp: dna.len(),
        },
    ];
    if req.include_wpre {
        elements.push(CassetteElement {
            name: "WPRE".to_string(),
            length_bp: WPRE_BP,
        });
    }
    elements.push(CassetteElement {
        name: format!("polyA ({})", req.poly_a.name()),
        length_bp: req.poly_a.length_bp(),
    });

    let total_bp: usize = elements.iter().map(|e| e.length_bp).sum();
    let capacity = req.vector.capacity_bp();
    let fits = total_bp <= capacity;
    let headroom = capacity.saturating_sub(total_bp);

    let mut notes = Vec::new();
    if fits {
        notes.push(format!(
            "Cassette is {total_bp} bp — fits the {} ~{} bp capacity with {} bp headroom.",
            req.vector.name(),
            capacity,
            headroom,
        ));
    } else {
        notes.push(format!(
            "Cassette is {total_bp} bp — OVER the {} ~{} bp capacity by {} bp; \
             it will not package.",
            req.vector.name(),
            capacity,
            total_bp - capacity,
        ));
        // Concrete size-reduction advice.
        if req.promoter != Promoter::SyntheticSmall {
            notes.push(
                "Consider a smaller synthetic promoter to recover payload space.".to_string(),
            );
        }
        if req.include_wpre {
            notes.push("Dropping the WPRE recovers ~600 bp.".to_string());
        }
        if req.vector == DeliveryVector::Aav {
            notes.push(
                "If the transgene is intrinsically too large for AAV, a \
                 dual-AAV split or a lentiviral vector is the usual route."
                    .to_string(),
            );
        }
    }
    if req.vector.integrates() {
        notes.push(
            "Lentivirus integrates into the host genome — see the safety \
             screen's insertional-mutagenesis note."
                .to_string(),
        );
    }

    Ok(ExpressionCassette {
        elements,
        total_bp,
        vector: req.vector,
        fits_vector: fits,
        headroom_bp: headroom,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_transgene() -> Vec<u8> {
        // ATG + a handful of codons + stop = 21 bp.
        b"ATGGCCCTGGAAGAAGCCTAA".to_vec()
    }

    /// A transgene of a requested length (in codons, start + stop
    /// included). Builds `ATG (GCC)* TAA`.
    fn transgene_of_codons(sense_codons: usize) -> Vec<u8> {
        let mut v = b"ATG".to_vec();
        for _ in 0..sense_codons {
            v.extend_from_slice(b"GCC");
        }
        v.extend_from_slice(b"TAA");
        v
    }

    #[test]
    fn small_cassette_fits_aav() {
        let req = CassetteRequest::new(small_transgene(), DeliveryVector::Aav);
        let c = design_cassette(&req).unwrap();
        assert!(c.fits_vector);
        assert!(c.headroom_bp > 0);
        assert!(c.notes[0].contains("fits"));
    }

    #[test]
    fn cassette_orders_elements() {
        let req = CassetteRequest::new(small_transgene(), DeliveryVector::Aav);
        let c = design_cassette(&req).unwrap();
        // First element is the vector overhead, last is the polyA.
        assert!(c.elements[0].name.contains("overhead"));
        assert!(c.elements.last().unwrap().name.contains("polyA"));
        assert!(c.elements.iter().any(|e| e.name == "transgene"));
    }

    #[test]
    fn oversized_transgene_overruns_aav() {
        // A ~5 kb transgene cannot fit AAV's 4.7 kb.
        let big = transgene_of_codons(1700); // ~5106 bp
        let req = CassetteRequest::new(big, DeliveryVector::Aav);
        let c = design_cassette(&req).unwrap();
        assert!(!c.fits_vector);
        assert_eq!(c.headroom_bp, 0);
        assert!(c.notes.iter().any(|n| n.contains("OVER")));
        // Advice to switch vectors is present.
        assert!(c
            .notes
            .iter()
            .any(|n| n.contains("dual-AAV") || n.contains("lentiviral")));
    }

    #[test]
    fn same_transgene_fits_lentivirus() {
        let big = transgene_of_codons(1700);
        let req = CassetteRequest::new(big, DeliveryVector::Lentivirus);
        let c = design_cassette(&req).unwrap();
        assert!(c.fits_vector, "lentivirus has the headroom AAV lacks");
    }

    #[test]
    fn wpre_consumes_payload_space() {
        let mut req = CassetteRequest::new(small_transgene(), DeliveryVector::Aav);
        let without = design_cassette(&req).unwrap();
        req.include_wpre = true;
        let with = design_cassette(&req).unwrap();
        assert_eq!(with.total_bp, without.total_bp + WPRE_BP);
    }

    #[test]
    fn lentivirus_notes_integration() {
        let req = CassetteRequest::new(small_transgene(), DeliveryVector::Lentivirus);
        let c = design_cassette(&req).unwrap();
        assert!(c.notes.iter().any(|n| n.contains("integrates")));
    }

    #[test]
    fn rejects_transgene_without_start_codon() {
        let req = CassetteRequest::new(b"GCCGCCGCCTAA".to_vec(), DeliveryVector::Aav);
        assert_eq!(
            design_cassette(&req).unwrap_err().code(),
            "genediting.invalid_target"
        );
    }

    #[test]
    fn rejects_transgene_not_multiple_of_three() {
        let req = CassetteRequest::new(b"ATGGCC".to_vec(), DeliveryVector::Aav);
        // ATGGCC is /3 but has no stop — still rejected.
        assert!(design_cassette(&req).is_err());
        let req2 = CassetteRequest::new(b"ATGGC".to_vec(), DeliveryVector::Aav);
        assert!(design_cassette(&req2).is_err());
    }

    #[test]
    fn capacity_ordering_is_sane() {
        assert!(DeliveryVector::Aav.capacity_bp() < DeliveryVector::Lentivirus.capacity_bp());
        assert!(DeliveryVector::Lentivirus.integrates());
        assert!(!DeliveryVector::Aav.integrates());
    }

    #[test]
    fn transgene_bp_helper() {
        let req = CassetteRequest::new(small_transgene(), DeliveryVector::Aav);
        let c = design_cassette(&req).unwrap();
        assert_eq!(c.transgene_bp(), small_transgene().len());
    }
}

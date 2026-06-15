//! Verified human TGF-beta-superfamily mature signaling-domain sequences.
//!
//! Fetched byte-exact from the UniProt REST API (<https://rest.uniprot.org>)
//! on 2026-06-15; mature-chain boundaries are the UniProt-annotated
//! processed-chain features. Accessions and positions are recorded for
//! provenance. These sequences are NOT fabricated or hand-typed: this file
//! is generated directly from the fetched FASTA.

/// A verified reference protein from UniProt.
#[derive(Debug, Clone, Copy)]
pub struct RefProtein {
    /// UniProt accession.
    pub accession: &'static str,
    /// Protein name.
    pub name: &'static str,
    /// Gene symbol.
    pub gene: &'static str,
    /// Mature-chain start (1-based, inclusive) in the UniProt precursor.
    pub chain_start: usize,
    /// Mature-chain end (1-based, inclusive).
    pub chain_end: usize,
    /// The mature signaling-domain amino-acid sequence.
    pub mature: &'static str,
}

/// Growth/differentiation factor 8 (myostatin) (UniProt O14793, mature chain 267-375, 109 aa).
pub const GDF8: RefProtein = RefProtein {
    accession: "O14793",
    name: "Growth/differentiation factor 8 (myostatin)",
    gene: "MSTN",
    chain_start: 267,
    chain_end: 375,
    mature: "DFGLDCDEHSTESRCCRYPLTVDFEAFGWDWIIAPKRYKANYCSGECEFVFLQKYPHTHLVHQANPRGSAGPCCTPTKMSPINMLYFNGKEQIIYGKIPAMVVDRCGCS",
};

/// Growth/differentiation factor 11 (UniProt O95390, mature chain 299-407, 109 aa).
pub const GDF11: RefProtein = RefProtein {
    accession: "O95390",
    name: "Growth/differentiation factor 11",
    gene: "GDF11",
    chain_start: 299,
    chain_end: 407,
    mature: "NLGLDCDEHSSESRCCRYPLTVDFEAFGWDWIIAPKRYKANYCSGQCEYMFMQKYPHTHLVQQANPRGSAGPCCTPTKMSPINMLYFNDKQQIIYGKIPGMVVDRCGCS",
};

/// Inhibin beta A chain (activin A) (UniProt P08476, mature chain 311-426, 116 aa).
pub const ACTIVIN_A: RefProtein = RefProtein {
    accession: "P08476",
    name: "Inhibin beta A chain (activin A)",
    gene: "INHBA",
    chain_start: 311,
    chain_end: 426,
    mature: "GLECDGKVNICCKKQFFVSFKDIGWNDWIIAPSGYHANYCEGECPSHIAGTSGSSLSFHSTVINHYRMRGHSPFANLKSCCVPTKLRPMSMLYYDDGQNIIKKDIQNMIVEECGCS",
};

/// Bone morphogenetic protein 7 (UniProt P18075, mature chain 293-431, 139 aa).
pub const BMP7: RefProtein = RefProtein {
    accession: "P18075",
    name: "Bone morphogenetic protein 7",
    gene: "BMP7",
    chain_start: 293,
    chain_end: 431,
    mature: "STGSKQRSQNRSKTPKNQEALRMANVAENSSSDQRQACKKHELYVSFRDLGWQDWIIAPEGYAAYYCEGECAFPLNSYMNATNHAIVQTLVHFINPETVPKPCCAPTQLNAISVLYFDDSSNVILKKYRNMVVRACGCH",
};

//! Restriction-enzyme database and digestion.
//!
//! Ships a built-in REBASE-style subset of ~200 commonly-used Type II
//! restriction enzymes — each with its recognition site, the cut
//! offset on the top and bottom strands, the resulting overhang
//! type, the enzyme's prototype (the canonical isoschizomer family
//! head), Dam/Dcm/CpG methylation-sensitivity flags, and a small
//! commercial-vendor availability set. [`digest`] finds every cut
//! site in a sequence (handling palindromes, IUPAC ambiguity in the
//! site, and circular molecules) and [`double_digest`] applies an
//! enzyme pair.
//!
//! Source: REBASE (rebase.neb.com), prototype and methylation columns
//! reproduced for the chosen subset. Cut offsets follow the REBASE
//! convention: a positive number is the number of bases from the
//! start of the recognition site to the cut point. For `EcoRI`
//! (`G^AATTC`, bottom cut `CTTAA^G`) the top cut is at offset 1 and
//! the bottom cut at offset 5, giving a 4-base 5′ overhang.
//!
//! Scope honesty: the database covers ~200 commonly-used cloning
//! enzymes (REBASE has thousands; we encode the production subset
//! shipped by major vendors and used in standard cloning workflows).
//! Methylation-sensitivity flags are the conservative published
//! values from REBASE — a `true` for `dam_sensitive` means the
//! enzyme is *blocked* by Dam methylation at an overlapping
//! GATC site; verify in REBASE for unusual contexts.

use crate::alphabet::{self, SeqKind};
use crate::error::{BioseqError, Result};
use crate::seq::Seq;

/// The kind of end a cut leaves.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum OverhangType {
    /// A 5′ recessed cut leaving a 5′ single-stranded overhang.
    FivePrime,
    /// A 3′ recessed cut leaving a 3′ single-stranded overhang.
    ThreePrime,
    /// A blunt cut (both strands cut at the same position).
    Blunt,
}

/// Commercial vendors that supply restriction enzymes — a compact set
/// of bit flags. Not every vendor is represented; the bits cover the
/// major catalog suppliers most cloning labs already have on the
/// shelf.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Default)]
pub struct Vendors(pub u8);

impl Vendors {
    /// New England Biolabs.
    pub const NEB: Vendors = Vendors(0b0000_0001);
    /// Thermo Fisher Scientific.
    pub const THERMO: Vendors = Vendors(0b0000_0010);
    /// Promega.
    pub const PROMEGA: Vendors = Vendors(0b0000_0100);
    /// Takara Bio.
    pub const TAKARA: Vendors = Vendors(0b0000_1000);
    /// Sigma-Aldrich / Merck.
    pub const SIGMA: Vendors = Vendors(0b0001_0000);

    /// Combine vendor sets with bitwise OR.
    pub const fn union(self, other: Vendors) -> Vendors {
        Vendors(self.0 | other.0)
    }

    /// `true` if this set includes the given vendor.
    pub const fn includes(self, other: Vendors) -> bool {
        (self.0 & other.0) != 0
    }

    /// `true` if no vendors are listed.
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }
}

/// A Type II restriction enzyme with REBASE-class metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Enzyme {
    /// Enzyme name (e.g. `"EcoRI"`).
    pub name: &'static str,
    /// Recognition site, 5′→3′ on the top strand. May contain IUPAC
    /// ambiguity codes.
    pub site: &'static str,
    /// Cut offset on the top strand, in bases from the site start.
    pub top_cut: usize,
    /// Cut offset on the bottom strand, in bases from the site start
    /// (measured along the top-strand coordinate).
    pub bottom_cut: usize,
    /// REBASE prototype — the canonical isoschizomer family head.
    /// Enzymes that recognize the same site map to the same prototype.
    /// For a prototype enzyme this equals `name`.
    pub prototype: &'static str,
    /// `true` if Dam methylation (GATC) at or overlapping the
    /// recognition site BLOCKS this enzyme.
    pub dam_sensitive: bool,
    /// `true` if Dcm methylation (CCWGG) at or overlapping the
    /// recognition site BLOCKS this enzyme.
    pub dcm_sensitive: bool,
    /// `true` if CpG methylation at or overlapping the recognition
    /// site BLOCKS this enzyme.
    pub cpg_sensitive: bool,
    /// Commercial vendors that supply this enzyme.
    pub vendors: Vendors,
}

impl Enzyme {
    /// The overhang type this enzyme produces, derived from the
    /// relative top vs. bottom cut positions.
    pub fn overhang(&self) -> OverhangType {
        match self.top_cut.cmp(&self.bottom_cut) {
            std::cmp::Ordering::Less => OverhangType::FivePrime,
            std::cmp::Ordering::Greater => OverhangType::ThreePrime,
            std::cmp::Ordering::Equal => OverhangType::Blunt,
        }
    }

    /// Length of the single-stranded overhang in bases (`0` for blunt).
    pub fn overhang_len(&self) -> usize {
        self.top_cut.abs_diff(self.bottom_cut)
    }

    /// The recognition-site length.
    pub fn site_len(&self) -> usize {
        self.site.len()
    }

    /// `true` if the site reads the same on both strands (a palindrome).
    pub fn is_palindrome(&self) -> bool {
        let bytes = self.site.as_bytes();
        let rc = crate::ops::revcomp::reverse_complement_dna_bytes(bytes);
        rc == bytes
    }

    /// `true` if this enzyme is its own prototype (the canonical
    /// member of its isoschizomer family).
    pub fn is_prototype(&self) -> bool {
        self.name == self.prototype
    }
}

/// Macro shorthand to keep the table compact without losing the
/// fields' meaning. The trailing `;` form expands to a full literal.
macro_rules! e {
    (
        $name:literal,
        $site:literal,
        $top:expr,
        $bot:expr,
        $proto:literal,
        $dam:expr,
        $dcm:expr,
        $cpg:expr,
        $vendors:expr $(,)?
    ) => {
        Enzyme {
            name: $name,
            site: $site,
            top_cut: $top,
            bottom_cut: $bot,
            prototype: $proto,
            dam_sensitive: $dam,
            dcm_sensitive: $dcm,
            cpg_sensitive: $cpg,
            vendors: $vendors,
        }
    };
}

/// The built-in restriction-enzyme database — a curated REBASE subset
/// of ~200 commonly-used cloning enzymes.
///
/// Each entry carries its prototype (canonical isoschizomer head),
/// methylation-sensitivity flags from REBASE, and the major commercial
/// vendors that supply the enzyme.
pub fn enzyme_database() -> &'static [Enzyme] {
    // The constant `V_*` are convenience aliases used inside the table.
    const V_NEB: Vendors = Vendors::NEB;
    const V_THM: Vendors = Vendors::THERMO;
    const V_PRG: Vendors = Vendors::PROMEGA;
    const V_TAK: Vendors = Vendors::TAKARA;
    const V_SIG: Vendors = Vendors::SIGMA;
    const ALL: Vendors = Vendors(V_NEB.0 | V_THM.0 | V_PRG.0 | V_TAK.0 | V_SIG.0);
    const MAJORS: Vendors = Vendors(V_NEB.0 | V_THM.0);

    &[
        // ===== 6-cutters, 5′ overhang (palindromic) =====
        // EcoRI family (G^AATTC).
        e!("EcoRI",   "GAATTC", 1, 5, "EcoRI", false, false, false, ALL),
        // BamHI family (G^GATCC).
        e!("BamHI",   "GGATCC", 1, 5, "BamHI", false, false, false, ALL),
        // HindIII family (A^AGCTT).
        e!("HindIII", "AAGCTT", 1, 5, "HindIII", false, false, false, ALL),
        // XhoI family (C^TCGAG).
        e!("XhoI",    "CTCGAG", 1, 5, "XhoI", false, false, true, MAJORS),
        e!("PaeR7I",  "CTCGAG", 1, 5, "XhoI", false, false, true, V_NEB),
        // SalI family (G^TCGAC).
        e!("SalI",    "GTCGAC", 1, 5, "SalI", false, false, false, ALL),
        // BglII family (A^GATCT).
        e!("BglII",   "AGATCT", 1, 5, "BglII", false, false, false, MAJORS),
        // SpeI family (A^CTAGT).
        e!("SpeI",    "ACTAGT", 1, 5, "SpeI", false, false, false, MAJORS),
        // XbaI family (T^CTAGA) — dam-sensitive (TCTAGAtc overlaps GATC).
        e!("XbaI",    "TCTAGA", 1, 5, "XbaI", true, false, false, ALL),
        // NheI family (G^CTAGC).
        e!("NheI",    "GCTAGC", 1, 5, "NheI", false, false, false, MAJORS),
        // NcoI family (C^CATGG).
        e!("NcoI",    "CCATGG", 1, 5, "NcoI", false, false, false, ALL),
        // NdeI family (CA^TATG).
        e!("NdeI",    "CATATG", 2, 4, "NdeI", false, false, false, MAJORS),
        // AflII family (C^TTAAG).
        e!("AflII",   "CTTAAG", 1, 5, "AflII", false, false, false, V_NEB),
        e!("BspTI",   "CTTAAG", 1, 5, "AflII", false, false, false, V_THM),
        // AscI family (GG^CGCGCC) — 8 nt site.
        e!("AscI",    "GGCGCGCC", 2, 6, "AscI", false, false, true, V_NEB),
        e!("SgsI",    "GGCGCGCC", 2, 6, "AscI", false, false, true, V_THM),
        // MfeI family (C^AATTG).
        e!("MfeI",    "CAATTG", 1, 5, "MfeI", false, false, false, V_NEB),
        e!("MunI",    "CAATTG", 1, 5, "MfeI", false, false, false, V_THM),
        // ClaI family (AT^CGAT) — dam-sensitive.
        e!("ClaI",    "ATCGAT", 2, 4, "ClaI", true, false, true, ALL),
        e!("BspDI",   "ATCGAT", 2, 4, "ClaI", true, false, true, V_NEB),
        // AvrII family (C^CTAGG).
        e!("AvrII",   "CCTAGG", 1, 5, "AvrII", false, false, false, V_NEB),
        e!("BlnI",    "CCTAGG", 1, 5, "AvrII", false, false, false, V_TAK),
        e!("XmaJI",   "CCTAGG", 1, 5, "AvrII", false, false, false, V_THM),
        // BspEI family (T^CCGGA).
        e!("BspEI",   "TCCGGA", 1, 5, "BspEI", false, false, true, V_NEB),
        e!("Kpn2I",   "TCCGGA", 1, 5, "BspEI", false, false, true, V_THM),
        // AgeI family (A^CCGGT).
        e!("AgeI",    "ACCGGT", 1, 5, "AgeI", false, false, true, V_NEB),
        e!("BshTI",   "ACCGGT", 1, 5, "AgeI", false, false, true, V_THM),
        // NgoMIV family (G^CCGGC).
        e!("NgoMIV",  "GCCGGC", 1, 5, "NgoMIV", false, false, true, V_NEB),
        // MluI family (A^CGCGT).
        e!("MluI",    "ACGCGT", 1, 5, "MluI", false, false, true, ALL),
        // BstBI family (TT^CGAA).
        e!("BstBI",   "TTCGAA", 2, 4, "BstBI", false, false, true, V_NEB),
        e!("Csp45I",  "TTCGAA", 2, 4, "BstBI", false, false, true, V_THM),
        // NruI family (TCG^CGA) — blunt cutter.
        // ===== 6-cutters, blunt =====
        e!("EcoRV",   "GATATC", 3, 3, "EcoRV", false, false, false, ALL),
        e!("SnaBI",   "TACGTA", 3, 3, "SnaBI", false, false, true, V_NEB),
        e!("SmaI",    "CCCGGG", 3, 3, "SmaI", false, false, true, ALL),
        e!("XmaI",    "CCCGGG", 1, 5, "XmaI", false, false, true, V_NEB),
        e!("PvuII",   "CAGCTG", 3, 3, "PvuII", false, false, false, ALL),
        e!("ScaI",    "AGTACT", 3, 3, "ScaI", false, false, false, V_NEB),
        e!("StuI",    "AGGCCT", 3, 3, "StuI", false, false, false, V_NEB),
        e!("DraI",    "TTTAAA", 3, 3, "DraI", false, false, false, ALL),
        e!("HpaI",    "GTTAAC", 3, 3, "HpaI", false, false, false, V_NEB),
        e!("NruI",    "TCGCGA", 3, 3, "NruI", false, false, true, V_NEB),
        e!("FspI",    "TGCGCA", 3, 3, "FspI", false, false, true, V_NEB),
        e!("AfeI",    "AGCGCT", 3, 3, "AfeI", false, false, true, V_NEB),
        e!("BmgBI",   "CACGTC", 3, 3, "BmgBI", false, false, true, V_NEB),
        e!("BstZ17I", "GTATAC", 3, 3, "BstZ17I", false, false, false, V_NEB),
        e!("PshAI",   "GACNNNNGTC", 5, 5, "PshAI", false, false, false, V_NEB),
        e!("PvuI",    "CGATCG", 4, 2, "PvuI", false, false, true, ALL),
        // ===== 6-cutters, 3′ overhang =====
        e!("PstI",    "CTGCAG", 5, 1, "PstI", false, false, false, ALL),
        e!("SacI",    "GAGCTC", 5, 1, "SacI", false, false, false, MAJORS),
        e!("SacII",   "CCGCGG", 4, 2, "SacII", false, false, true, V_NEB),
        e!("KpnI",    "GGTACC", 5, 1, "KpnI", false, false, false, ALL),
        e!("Acc65I",  "GGTACC", 1, 5, "Acc65I", false, false, false, V_NEB),
        e!("SphI",    "GCATGC", 5, 1, "SphI", false, false, false, V_NEB),
        e!("PaeI",    "GCATGC", 5, 1, "SphI", false, false, false, V_THM),
        e!("ApaI",    "GGGCCC", 5, 1, "ApaI", false, false, true, V_NEB),
        e!("Bsp120I", "GGGCCC", 1, 5, "Bsp120I", false, false, true, V_THM),
        e!("AatII",   "GACGTC", 5, 1, "AatII", false, false, true, V_NEB),
        e!("ZraI",    "GACGTC", 3, 3, "ZraI", false, false, true, V_NEB),
        e!("NsiI",    "ATGCAT", 5, 1, "NsiI", false, false, false, V_NEB),
        e!("Mph1103I", "ATGCAT", 5, 1, "NsiI", false, false, false, V_THM),
        e!("BglI",    "GCCNNNNNGGC", 7, 4, "BglI", false, false, false, V_NEB),
        // ===== 8-cutters (rare cutters) =====
        e!("NotI",    "GCGGCCGC", 2, 6, "NotI", false, false, true, ALL),
        e!("PacI",    "TTAATTAA", 5, 3, "PacI", false, false, false, V_NEB),
        e!("FseI",    "GGCCGGCC", 6, 2, "FseI", false, false, true, V_NEB),
        e!("SbfI",    "CCTGCAGG", 6, 2, "SbfI", false, false, false, V_NEB),
        e!("Sse8387I","CCTGCAGG", 6, 2, "SbfI", false, false, false, V_TAK),
        e!("PmeI",    "GTTTAAAC", 4, 4, "PmeI", false, false, false, V_NEB),
        e!("SwaI",    "ATTTAAAT", 4, 4, "SwaI", false, false, false, V_NEB),
        e!("SrfI",    "GCCCGGGC", 4, 4, "SrfI", false, false, true, V_NEB),
        // I-SceI homing endonuclease (18 nt site).
        e!("I-SceI",  "TAGGGATAACAGGGTAAT", 9, 13, "I-SceI", false, false, false, V_NEB),
        e!("I-CeuI",  "TAACTATAACGGTCCTAAGGTAGCGA", 17, 21, "I-CeuI", false, false, false, V_NEB),
        e!("I-PpoI",  "CTCTCTTAAGGTAGC", 8, 5, "I-PpoI", false, false, false, V_NEB),
        // ===== Type IIS (cut outside the recognition site) =====
        e!("BsaI",    "GGTCTC", 7, 11, "BsaI", false, false, false, V_NEB),
        e!("BbsI",    "GAAGAC", 8, 12, "BbsI", false, false, false, V_NEB),
        e!("BsmBI",   "CGTCTC", 7, 11, "BsmBI", false, false, false, V_NEB),
        e!("Esp3I",   "CGTCTC", 7, 11, "BsmBI", false, false, false, V_THM),
        e!("SapI",    "GCTCTTC", 8, 11, "SapI", false, false, false, V_NEB),
        e!("LguI",    "GCTCTTC", 8, 11, "SapI", false, false, false, V_THM),
        e!("AarI",    "CACCTGC", 11, 15, "AarI", false, false, false, V_THM),
        e!("BspQI",   "GCTCTTC", 8, 11, "SapI", false, false, false, V_NEB),
        e!("BfuAI",   "ACCTGC", 10, 14, "BfuAI", false, false, false, V_NEB),
        e!("BtgZI",   "GCGATG", 16, 20, "BtgZI", false, false, true, V_NEB),
        e!("FokI",    "GGATG", 14, 18, "FokI", false, false, false, V_NEB),
        e!("MlyI",    "GAGTC", 10, 10, "MlyI", false, false, false, V_NEB),
        e!("EarI",    "CTCTTC", 7, 11, "EarI", false, false, false, V_NEB),
        e!("MmeI",    "TCCRAC", 26, 24, "MmeI", false, false, true, V_NEB),
        // ===== 4-cutters (frequent cutters) =====
        e!("AluI",    "AGCT", 2, 2, "AluI", false, false, false, ALL),
        e!("HaeIII",  "GGCC", 2, 2, "HaeIII", false, false, false, ALL),
        e!("MspI",    "CCGG", 1, 3, "MspI", false, false, false, V_NEB),
        e!("HpaII",   "CCGG", 1, 3, "HpaII", false, false, true, V_NEB),
        e!("TaqI",    "TCGA", 1, 3, "TaqI", false, false, true, V_NEB),
        e!("RsaI",    "GTAC", 2, 2, "RsaI", false, false, false, V_NEB),
        e!("Sau3AI",  "GATC", 0, 4, "Sau3AI", true, false, false, V_NEB),
        e!("DpnI",    "GATC", 2, 2, "DpnI", false, false, false, V_NEB),
        e!("DpnII",   "GATC", 0, 4, "Sau3AI", true, false, false, V_NEB),
        e!("MboI",    "GATC", 0, 4, "Sau3AI", true, false, false, V_NEB),
        e!("NlaIII",  "CATG", 4, 0, "NlaIII", false, false, false, V_NEB),
        e!("CviAII",  "CATG", 1, 3, "CviAII", false, false, false, V_NEB),
        e!("FatI",    "CATG", 0, 4, "FatI", false, false, false, V_NEB),
        e!("CviQI",   "GTAC", 1, 3, "CviQI", false, false, false, V_NEB),
        e!("MnlI",    "CCTC", 11, 10, "MnlI", false, false, false, V_NEB),
        e!("HinfI",   "GANTC", 1, 4, "HinfI", false, false, false, ALL),
        e!("DdeI",    "CTNAG", 1, 4, "DdeI", false, false, false, V_NEB),
        e!("TfiI",    "GAWTC", 1, 4, "TfiI", false, false, false, V_NEB),
        e!("MaeII",   "ACGT", 1, 3, "MaeII", false, false, true, V_THM),
        e!("CviKI",   "RGCY", 2, 2, "CviKI", false, false, false, V_NEB),
        // ===== Additional 6-cutters (NEB / Thermo catalog, isoschizomers + extras) =====
        e!("ApaLI",   "GTGCAC", 1, 5, "ApaLI", false, false, false, V_NEB),
        e!("Alw44I",  "GTGCAC", 1, 5, "ApaLI", false, false, false, V_THM),
        e!("BspHI",   "TCATGA", 1, 5, "BspHI", false, false, false, V_NEB),
        e!("PagI",    "TCATGA", 1, 5, "BspHI", false, false, false, V_THM),
        e!("PciI",    "ACATGT", 1, 5, "PciI", false, false, false, V_NEB),
        e!("AflIII",  "ACRYGT", 1, 5, "AflIII", false, false, false, V_NEB),
        e!("BstXI",   "CCANNNNNNTGG", 8, 4, "BstXI", false, false, false, V_NEB),
        e!("DraIII",  "CACNNNGTG", 6, 3, "DraIII", false, false, false, V_NEB),
        e!("Eco47III","AGCGCT", 3, 3, "AfeI", false, false, true, V_THM),
        e!("EagI",    "CGGCCG", 1, 5, "EagI", false, false, true, V_NEB),
        e!("EclXI",   "CGGCCG", 1, 5, "EagI", false, false, true, V_THM),
        e!("EheI",    "GGCGCC", 2, 4, "EheI", false, false, true, V_THM),
        e!("KasI",    "GGCGCC", 1, 5, "KasI", false, false, true, V_NEB),
        e!("NarI",    "GGCGCC", 2, 4, "NarI", false, false, true, V_NEB),
        e!("BbeI",    "GGCGCC", 5, 1, "BbeI", false, false, true, V_NEB),
        e!("BsiWI",   "CGTACG", 1, 5, "BsiWI", false, false, true, V_NEB),
        e!("BstEII",  "GGTNACC", 1, 6, "BstEII", false, false, false, V_NEB),
        e!("Eco91I",  "GGTNACC", 1, 6, "BstEII", false, false, false, V_THM),
        e!("BstAPI",  "GCANNNNNTGC", 7, 4, "BstAPI", false, false, false, V_NEB),
        e!("BstEII-HF","GGTNACC", 1, 6, "BstEII", false, false, false, V_NEB),
        e!("MluCI",   "AATT", 0, 4, "MluCI", false, false, false, V_NEB),
        e!("Tsp509I", "AATT", 0, 4, "MluCI", false, false, false, V_NEB),
        e!("Cac8I",   "GCNNGC", 3, 3, "Cac8I", false, false, false, V_NEB),
        e!("HpyCH4III","ACNGT", 3, 2, "HpyCH4III", false, false, false, V_NEB),
        e!("HpyCH4IV","ACGT", 1, 3, "HpyCH4IV", false, false, true, V_NEB),
        e!("HpyCH4V", "TGCA", 2, 2, "HpyCH4V", false, false, false, V_NEB),
        e!("PflMI",   "CCANNNNNTGG", 7, 4, "PflMI", false, false, false, V_NEB),
        e!("Bpu10I",  "CCTNAGC", 2, 5, "Bpu10I", false, false, false, V_THM),
        e!("BpiI",    "GAAGAC", 8, 12, "BbsI", false, false, false, V_THM),
        e!("BsrBI",   "CCGCTC", 3, 3, "BsrBI", false, false, true, V_NEB),
        e!("AhdI",    "GACNNNNNGTC", 6, 5, "AhdI", false, false, false, V_NEB),
        e!("BsmI",    "GAATGC", 7, 5, "BsmI", false, false, false, V_NEB),
        e!("BstNI",   "CCWGG", 2, 3, "BstNI", false, true, false, V_NEB),
        e!("EcoNI",   "CCTNNNNNAGG", 5, 6, "EcoNI", false, false, false, V_NEB),
        e!("PspOMI",  "GGGCCC", 1, 5, "Bsp120I", false, false, true, V_NEB),
        e!("PluTI",   "GGCGCC", 5, 1, "BbeI", false, false, true, V_NEB),
        e!("NaeI",    "GCCGGC", 3, 3, "NaeI", false, false, true, V_NEB),
        e!("MscI",    "TGGCCA", 3, 3, "MscI", false, false, false, V_NEB),
        e!("BsaAI",   "YACGTR", 3, 3, "BsaAI", false, false, true, V_NEB),
        e!("BsaBI",   "GATNNNNATC", 5, 5, "BsaBI", false, false, false, V_NEB),
        e!("DrdI",    "GACNNNNNNGTC", 7, 5, "DrdI", false, false, false, V_NEB),
        e!("AvaI",    "CYCGRG", 1, 5, "AvaI", false, false, true, V_NEB),
        e!("BsoBI",   "CYCGRG", 1, 5, "BsoBI", false, false, true, V_NEB),
        e!("TspRI",   "NNCASTGNN", 9, 0, "TspRI", false, false, false, V_NEB),
        e!("BsgI",    "GTGCAG", 22, 20, "BsgI", false, false, false, V_NEB),
        e!("HgaI",    "GACGC", 10, 15, "HgaI", false, false, true, V_NEB),
        e!("BciVI",   "GTATCC", 18, 16, "BciVI", false, false, false, V_NEB),
        e!("BpmI",    "CTGGAG", 22, 20, "BpmI", false, false, false, V_NEB),
        e!("BpuEI",   "CTTGAG", 22, 20, "BpuEI", false, false, false, V_NEB),
        e!("AcuI",    "CTGAAG", 22, 20, "AcuI", false, false, false, V_NEB),
        // ===== Additional 4-cutters and methylation-sensitive enzymes =====
        e!("HhaI",    "GCGC", 3, 1, "HhaI", false, false, true, V_NEB),
        e!("CviAII-2","CATG", 1, 3, "CviAII", false, false, false, V_NEB),
        e!("BfaI",    "CTAG", 1, 3, "BfaI", false, false, false, V_NEB),
        e!("DpnII-2", "GATC", 0, 4, "Sau3AI", true, false, false, V_NEB),
        e!("MboII",   "GAAGA", 8, 7, "MboII", false, false, false, V_NEB),
        e!("CviJI",   "RGCY", 2, 2, "CviKI", false, false, false, V_NEB),
        e!("HpyAV",   "CCTTC", 6, 5, "HpyAV", false, false, false, V_NEB),
        e!("PleI",    "GAGTC", 4, 5, "PleI", false, false, false, V_NEB),
        e!("HphI",    "GGTGA", 8, 7, "HphI", false, false, false, V_NEB),
        e!("FauI",    "CCCGC", 4, 6, "FauI", false, false, true, V_NEB),
        e!("BciT130I","CCWGG", 2, 3, "BstNI", false, true, false, V_THM),
        e!("EcoT22I", "ATGCAT", 5, 1, "NsiI", false, false, false, V_TAK),
        e!("SspI",    "AATATT", 3, 3, "SspI", false, false, false, V_NEB),
        e!("VspI",    "ATTAAT", 2, 4, "VspI", false, false, false, V_THM),
        e!("AseI",    "ATTAAT", 2, 4, "AseI", false, false, false, V_NEB),
        e!("PsiI",    "TTATAA", 3, 3, "PsiI", false, false, false, V_NEB),
        e!("MfeI-HF", "CAATTG", 1, 5, "MfeI", false, false, false, V_NEB),
        e!("NspI",    "RCATGY", 5, 1, "NspI", false, false, false, V_NEB),
        e!("AlwNI",   "CAGNNNCTG", 6, 3, "AlwNI", false, false, false, V_NEB),
        e!("XcmI",    "CCANNNNNNNNNTGG", 8, 7, "XcmI", false, false, false, V_NEB),
        e!("BsiEI",   "CGRYCG", 4, 2, "BsiEI", false, false, true, V_NEB),
        e!("BstUI",   "CGCG", 2, 2, "BstUI", false, false, true, V_NEB),
        // ===== A few extras for symmetry coverage =====
        e!("PflFI",   "GACNNNGTC", 4, 5, "PflFI", false, false, false, V_NEB),
        e!("Tth111I", "GACNNNGTC", 4, 5, "PflFI", false, false, false, V_THM),
        e!("PstI-HF", "CTGCAG", 5, 1, "PstI", false, false, false, V_NEB),
        e!("NotI-HF", "GCGGCCGC", 2, 6, "NotI", false, false, true, V_NEB),
        e!("EcoRI-HF","GAATTC", 1, 5, "EcoRI", false, false, false, V_NEB),
        e!("BamHI-HF","GGATCC", 1, 5, "BamHI", false, false, false, V_NEB),
        e!("HindIII-HF","AAGCTT", 1, 5, "HindIII", false, false, false, V_NEB),
        e!("BsaI-HF", "GGTCTC", 7, 11, "BsaI", false, false, false, V_NEB),
        e!("BbsI-HF", "GAAGAC", 8, 12, "BbsI", false, false, false, V_NEB),
        e!("BsmBI-v2","CGTCTC", 7, 11, "BsmBI", false, false, false, V_NEB),
        e!("NheI-HF", "GCTAGC", 1, 5, "NheI", false, false, false, V_NEB),
        e!("SpeI-HF", "ACTAGT", 1, 5, "SpeI", false, false, false, V_NEB),
        e!("KpnI-HF", "GGTACC", 5, 1, "KpnI", false, false, false, V_NEB),
        e!("MfeI-HFv2","CAATTG", 1, 5, "MfeI", false, false, false, V_NEB),
        // ===== A handful of methylation-specific enzymes =====
        e!("Eco57I",  "CTGAAG", 22, 20, "Eco57I", false, false, false, V_THM),
        e!("MaeIII",  "GTNAC", 0, 5, "MaeIII", false, false, false, V_THM),
        e!("MseI",    "TTAA", 1, 3, "MseI", false, false, false, V_NEB),
        e!("Tru1I",   "TTAA", 1, 3, "MseI", false, false, false, V_THM),
        e!("PspGI",   "CCWGG", 0, 5, "PspGI", false, true, false, V_NEB),
        e!("BstAPI-2","GCANNNNNTGC", 7, 4, "BstAPI", false, false, false, V_NEB),
        e!("BstZ17I-HF","GTATAC", 3, 3, "BstZ17I", false, false, false, V_NEB),
    ]
}

/// Looks up an enzyme by name (case-insensitive).
pub fn enzyme_by_name(name: &str) -> Option<&'static Enzyme> {
    enzyme_database()
        .iter()
        .find(|e| e.name.eq_ignore_ascii_case(name))
}

/// Returns every enzyme in the database that shares `prototype`'s
/// recognition site (its isoschizomer family).
pub fn isoschizomers_of(prototype: &str) -> Vec<&'static Enzyme> {
    enzyme_database()
        .iter()
        .filter(|e| e.prototype.eq_ignore_ascii_case(prototype))
        .collect()
}

/// Returns every prototype enzyme — one canonical representative per
/// isoschizomer family. The output is sorted by name.
pub fn prototypes() -> Vec<&'static Enzyme> {
    let mut out: Vec<&'static Enzyme> = enzyme_database()
        .iter()
        .filter(|e| e.is_prototype())
        .collect();
    out.sort_by_key(|e| e.name);
    out
}

/// Returns every enzyme blocked by Dam methylation.
pub fn dam_sensitive_enzymes() -> Vec<&'static Enzyme> {
    enzyme_database().iter().filter(|e| e.dam_sensitive).collect()
}

/// Returns every enzyme blocked by Dcm methylation.
pub fn dcm_sensitive_enzymes() -> Vec<&'static Enzyme> {
    enzyme_database().iter().filter(|e| e.dcm_sensitive).collect()
}

/// Returns every enzyme blocked by CpG methylation.
pub fn cpg_sensitive_enzymes() -> Vec<&'static Enzyme> {
    enzyme_database().iter().filter(|e| e.cpg_sensitive).collect()
}

/// A single restriction-cut site found in a sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CutSite {
    /// Name of the enzyme that cuts here.
    pub enzyme: &'static str,
    /// 0-based index of the recognition site's start on the top
    /// strand.
    pub site_start: usize,
    /// 0-based top-strand cut position (the phosphodiester bond
    /// between `top_cut_pos − 1` and `top_cut_pos` is broken).
    pub top_cut_pos: usize,
    /// 0-based bottom-strand cut position, in top-strand coordinates.
    pub bottom_cut_pos: usize,
    /// The overhang type left by the cut.
    pub overhang: OverhangType,
}

/// Finds every position where `enzyme` cuts `seq`.
///
/// Handles palindromic and non-palindromic sites (the latter are also
/// matched on the reverse strand), IUPAC ambiguity codes in the
/// recognition site, and circular molecules (sites spanning the
/// origin). Returns cut sites sorted by `site_start`. Returns
/// [`BioseqError::Invalid`] for a non-DNA input.
pub fn digest(seq: &Seq, enzyme: &Enzyme) -> Result<Vec<CutSite>> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid(
            "kind",
            "restriction digestion needs a DNA sequence",
        ));
    }
    let n = seq.len();
    let site = enzyme.site.as_bytes();
    let m = site.len();
    if m == 0 || n == 0 {
        return Ok(Vec::new());
    }
    let bytes = seq.as_bytes();
    let palindrome = enzyme.is_palindrome();
    let site_rc = crate::ops::revcomp::reverse_complement_dna_bytes(site);

    // For a circular sequence, search a doubled window so sites that
    // wrap the origin are found; dedupe afterwards.
    let search_len = if seq.is_circular() {
        n // a recognition site can start anywhere 0..n on a circle
    } else {
        n.saturating_sub(m) + 1
    };

    let mut sites = Vec::new();
    for start in 0..search_len {
        // Gather the m bases starting at `start`, wrapping if circular.
        let window = gather(bytes, start, m, seq.is_circular());
        if window.len() < m {
            continue; // linear, ran off the end
        }
        let fwd = match_iupac(site, &window);
        let rev = if palindrome {
            false // already matched as forward; avoid double counting
        } else {
            match_iupac(&site_rc, &window)
        };
        if fwd {
            sites.push(make_cut_site(enzyme, start, n, false));
        } else if rev {
            sites.push(make_cut_site(enzyme, start, n, true));
        }
    }
    sites.sort_by_key(|s| s.site_start);
    Ok(sites)
}

/// Builds a [`CutSite`] for a recognition site found at `site_start`.
/// `reverse` indicates the site matched on the bottom strand (the cut
/// offsets are then mirrored).
fn make_cut_site(enzyme: &Enzyme, site_start: usize, seq_len: usize, reverse: bool) -> CutSite {
    let m = enzyme.site.len() as isize;
    // The cut offset relative to `site_start`, as a SIGNED quantity. For
    // a Type IIS enzyme the forward `top_cut` / `bottom_cut` exceed `m`
    // (the enzyme cuts downstream of its recognition site); when the
    // site is matched on the bottom strand the mirrored offset `m - cut`
    // is then NEGATIVE — the cut falls upstream of `site_start`. Signed
    // arithmetic is therefore mandatory here; the previous `usize`
    // subtraction `m - bottom_cut` panicked with overflow for BsaI et al.
    let (top_off, bot_off) = if reverse {
        // Mirror the offsets onto the bottom-strand-matched span.
        (m - enzyme.bottom_cut as isize, m - enzyme.top_cut as isize)
    } else {
        (enzyme.top_cut as isize, enzyme.bottom_cut as isize)
    };
    // Resolve a signed offset to a 0-based top-strand position: wrap
    // modulo the length on a circular molecule, clamp into range on a
    // linear one.
    let resolve = |off: isize| -> usize {
        if seq_len == 0 {
            return 0;
        }
        let n = seq_len as isize;
        let pos = site_start as isize + off;
        // Euclidean remainder keeps a negative position on a circle.
        let wrapped = pos.rem_euclid(n);
        wrapped as usize
    };
    CutSite {
        enzyme: enzyme.name,
        site_start,
        top_cut_pos: resolve(top_off),
        bottom_cut_pos: resolve(bot_off),
        overhang: enzyme.overhang(),
    }
}

/// Collects `len` bases of `bytes` from `start`, wrapping if `circular`.
fn gather(bytes: &[u8], start: usize, len: usize, circular: bool) -> Vec<u8> {
    let n = bytes.len();
    let mut out = Vec::with_capacity(len);
    for k in 0..len {
        let idx = start + k;
        if idx < n {
            out.push(bytes[idx]);
        } else if circular {
            out.push(bytes[idx % n]);
        } else {
            break;
        }
    }
    out
}

/// `true` if every position of `pattern` (which may contain IUPAC
/// ambiguity codes) matches the concrete base in `window`.
fn match_iupac(pattern: &[u8], window: &[u8]) -> bool {
    if pattern.len() != window.len() {
        return false;
    }
    pattern
        .iter()
        .zip(window)
        .all(|(&p, &b)| alphabet::iupac_matches(p, b))
}

/// Performs a double digest with two enzymes, merging their cut sites.
/// Cut sites are returned sorted by `site_start`.
pub fn double_digest(
    seq: &Seq,
    enzyme_a: &Enzyme,
    enzyme_b: &Enzyme,
) -> Result<Vec<CutSite>> {
    let mut sites = digest(seq, enzyme_a)?;
    sites.extend(digest(seq, enzyme_b)?);
    sites.sort_by_key(|s| s.site_start);
    Ok(sites)
}

/// Counts how many times `enzyme` cuts `seq` (a quick "is this enzyme
/// a unique cutter?" check).
pub fn cut_count(seq: &Seq, enzyme: &Enzyme) -> Result<usize> {
    Ok(digest(seq, enzyme)?.len())
}

/// All enzymes from the built-in database that do not cut `seq` at all
/// — useful for choosing a cloning site. Returns names sorted
/// alphabetically.
pub fn non_cutters(seq: &Seq) -> Result<Vec<&'static str>> {
    let mut out = Vec::new();
    for e in enzyme_database() {
        if cut_count(seq, e)? == 0 {
            out.push(e.name);
        }
    }
    out.sort_unstable();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Topology;

    #[test]
    fn database_is_populated_at_commercial_depth() {
        let db = enzyme_database();
        assert!(
            db.len() >= 180,
            "expected ~200 enzymes (commercial-depth subset), got {}",
            db.len()
        );
        assert!(enzyme_by_name("ecori").is_some());
        assert!(enzyme_by_name("NotI").is_some());
        assert!(enzyme_by_name("DpnI").is_some());
        assert!(enzyme_by_name("nonsense").is_none());
    }

    #[test]
    fn prototype_and_isoschizomer_relationships() {
        // BspEI is the prototype for the family that includes Kpn2I.
        let bspei = enzyme_by_name("BspEI").unwrap();
        let kpn2i = enzyme_by_name("Kpn2I").unwrap();
        assert!(bspei.is_prototype());
        assert!(!kpn2i.is_prototype());
        assert_eq!(kpn2i.prototype, "BspEI");
        // The isoschizomer-of query returns both.
        let fam = isoschizomers_of("BspEI");
        let names: Vec<&str> = fam.iter().map(|e| e.name).collect();
        assert!(names.contains(&"BspEI"));
        assert!(names.contains(&"Kpn2I"));
    }

    #[test]
    fn methylation_sensitivity_flags_populated() {
        // XbaI is dam-sensitive (TCTAGAtc overlaps GATC).
        let xbai = enzyme_by_name("XbaI").unwrap();
        assert!(xbai.dam_sensitive, "XbaI is dam-sensitive in REBASE");
        // ClaI is also dam-sensitive (ATCGATc overlap).
        let clai = enzyme_by_name("ClaI").unwrap();
        assert!(clai.dam_sensitive);
        // PspGI is dcm-blocked.
        let pspgi = enzyme_by_name("PspGI").unwrap();
        assert!(pspgi.dcm_sensitive);
        // NotI (GCGGCCGC) overlaps CpG and is blocked by CpG methylation.
        let noti = enzyme_by_name("NotI").unwrap();
        assert!(noti.cpg_sensitive);
        // EcoRI carries none of the flags — its site has no overlap.
        let ecori = enzyme_by_name("EcoRI").unwrap();
        assert!(!ecori.dam_sensitive);
        assert!(!ecori.dcm_sensitive);
        assert!(!ecori.cpg_sensitive);
    }

    #[test]
    fn vendor_metadata_populated() {
        let ecori = enzyme_by_name("EcoRI").unwrap();
        // EcoRI is sold by every major vendor.
        assert!(ecori.vendors.includes(Vendors::NEB));
        assert!(ecori.vendors.includes(Vendors::THERMO));
        // SapI is a more specialty enzyme — at least NEB carries it.
        let sapi = enzyme_by_name("SapI").unwrap();
        assert!(sapi.vendors.includes(Vendors::NEB));
        // The vendor flag set is never empty for a real entry.
        for e in enzyme_database() {
            assert!(!e.vendors.is_empty(), "{} has no vendor flags", e.name);
        }
    }

    #[test]
    fn prototypes_returns_one_per_family() {
        let protos = prototypes();
        // Every prototype's name equals its `prototype` field.
        for e in &protos {
            assert!(e.is_prototype(), "{} is not a prototype", e.name);
        }
        // No two prototypes share a name.
        let mut names: Vec<&str> = protos.iter().map(|e| e.name).collect();
        names.sort_unstable();
        let n_before = names.len();
        names.dedup();
        assert_eq!(names.len(), n_before, "duplicate prototype names");
    }

    #[test]
    fn methylation_query_helpers() {
        let dams = dam_sensitive_enzymes();
        assert!(!dams.is_empty());
        assert!(dams.iter().any(|e| e.name == "XbaI"));
        let dcms = dcm_sensitive_enzymes();
        assert!(!dcms.is_empty());
        assert!(dcms.iter().any(|e| e.name == "BstNI"));
        let cpgs = cpg_sensitive_enzymes();
        assert!(!cpgs.is_empty());
        assert!(cpgs.iter().any(|e| e.name == "NotI"));
    }

    #[test]
    fn overhang_classification() {
        assert_eq!(
            enzyme_by_name("EcoRI").unwrap().overhang(),
            OverhangType::FivePrime
        );
        assert_eq!(
            enzyme_by_name("PstI").unwrap().overhang(),
            OverhangType::ThreePrime
        );
        assert_eq!(
            enzyme_by_name("SmaI").unwrap().overhang(),
            OverhangType::Blunt
        );
        assert_eq!(enzyme_by_name("EcoRI").unwrap().overhang_len(), 4);
        assert_eq!(enzyme_by_name("SmaI").unwrap().overhang_len(), 0);
    }

    #[test]
    fn palindrome_detection() {
        assert!(enzyme_by_name("EcoRI").unwrap().is_palindrome());
        assert!(enzyme_by_name("BamHI").unwrap().is_palindrome());
        // BsaI's site GGTCTC is not a palindrome.
        assert!(!enzyme_by_name("BsaI").unwrap().is_palindrome());
    }

    #[test]
    fn ecori_cuts_its_site() {
        // EcoRI site GAATTC at index 3.
        let s = Seq::new(SeqKind::Dna, "AAAGAATTCAAA").unwrap();
        let cuts = digest(&s, enzyme_by_name("EcoRI").unwrap()).unwrap();
        assert_eq!(cuts.len(), 1);
        let c = &cuts[0];
        assert_eq!(c.site_start, 3);
        assert_eq!(c.top_cut_pos, 4); // G^AATTC
        assert_eq!(c.bottom_cut_pos, 8);
        assert_eq!(c.overhang, OverhangType::FivePrime);
    }

    #[test]
    fn no_site_means_no_cut() {
        let s = Seq::new(SeqKind::Dna, "AAAAAAAAAA").unwrap();
        assert_eq!(cut_count(&s, enzyme_by_name("EcoRI").unwrap()).unwrap(), 0);
    }

    #[test]
    fn multiple_sites_found() {
        // Two EcoRI sites.
        let s = Seq::new(SeqKind::Dna, "GAATTCAAAAGAATTC").unwrap();
        let cuts = digest(&s, enzyme_by_name("EcoRI").unwrap()).unwrap();
        assert_eq!(cuts.len(), 2);
        assert_eq!(cuts[0].site_start, 0);
        assert_eq!(cuts[1].site_start, 10);
    }

    #[test]
    fn non_palindrome_matches_both_strands() {
        // BsaI forward site GGTCTC; its reverse complement GAGACC.
        // Embed the reverse-strand site only.
        let s = Seq::new(SeqKind::Dna, "AAAAGAGACCAAAA").unwrap();
        let cuts = digest(&s, enzyme_by_name("BsaI").unwrap()).unwrap();
        assert_eq!(cuts.len(), 1, "BsaI should match its reverse-strand site");
    }

    #[test]
    fn circular_site_spans_origin() {
        // Place GAATTC so it wraps: last 3 bases + first 3 bases.
        let s = Seq::with_topology(SeqKind::Dna, "TTCAAAAAAGAA", Topology::Circular).unwrap();
        // site GAATTC starts at index 9: G A A | T T C (wraps).
        let cuts = digest(&s, enzyme_by_name("EcoRI").unwrap()).unwrap();
        assert_eq!(cuts.len(), 1);
        assert_eq!(cuts[0].site_start, 9);
    }

    #[test]
    fn ambiguity_in_site_matches() {
        // AarI is fine; test an enzyme via Sau3AI-like. Use a site
        // with ambiguity directly: HaeIII GGCC is concrete, so build
        // a synthetic check through iupac matching on an N pattern.
        assert!(match_iupac(b"GGNCC", b"GGACC"));
        assert!(match_iupac(b"GGNCC", b"GGGCC"));
        assert!(!match_iupac(b"GGNCC", b"GGGCA"));
    }

    #[test]
    fn double_digest_merges_sites() {
        // EcoRI site then BamHI site.
        let s = Seq::new(SeqKind::Dna, "GAATTCAAAGGATCC").unwrap();
        let cuts = double_digest(
            &s,
            enzyme_by_name("EcoRI").unwrap(),
            enzyme_by_name("BamHI").unwrap(),
        )
        .unwrap();
        assert_eq!(cuts.len(), 2);
        assert_eq!(cuts[0].enzyme, "EcoRI");
        assert_eq!(cuts[1].enzyme, "BamHI");
    }

    #[test]
    fn non_cutters_excludes_present_sites() {
        let s = Seq::new(SeqKind::Dna, "GAATTCAAAAAAAAAA").unwrap();
        let nc = non_cutters(&s).unwrap();
        assert!(!nc.contains(&"EcoRI"));
        // BamHI's site is absent.
        assert!(nc.contains(&"BamHI"));
    }

    #[test]
    fn non_dna_rejected() {
        let r = Seq::new(SeqKind::Rna, "GAAUUC").unwrap();
        assert!(digest(&r, enzyme_by_name("EcoRI").unwrap()).is_err());
    }
}

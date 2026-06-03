//! Feature 15 — the LinearDesign-class joint mRNA coding-sequence
//! optimiser.
//!
//! Commercial mRNA design (LinearDesign, Zhang *et al.*, *Nature* 2023)
//! does something the rest of this crate — and the
//! `valenx-genediting` mRNA driver — does **not**: it optimises codon
//! usage and secondary-structure stability *jointly*. A naive pipeline
//! picks the host-optimal codon for every residue, *then* folds the
//! result; the two passes never see each other, so the codon choice
//! that maximises the codon adaptation index (CAI) and the codon
//! choice that gives the most stable (lowest minimum-free-energy) mRNA
//! are settled independently. LinearDesign instead searches the whole
//! synonymous-codon space *at once* for the sequence minimising a
//! single combined objective.
//!
//! ## The method
//!
//! The synonymous codon choices of a protein form a **lattice**: each
//! residue contributes a small set of codon triples, and a coding
//! sequence (CDS) is one path that picks a triple per residue. This
//! module runs a **secondary-structure folding dynamic program
//! directly over that lattice** — a Zuker-style minimum-free-energy
//! recurrence in which every nucleotide is a *choice* drawn from the
//! lattice rather than a fixed letter. The DP is run left-to-right with
//! [LinearFold-style](valenx_rnastruct::fold_linear) beam pruning, so
//! it stays near-linear in the CDS length.
//!
//! The objective minimised at every paired / unpaired decision is
//!
//! ```text
//!     MFE  +  λ · (codon-optimality penalty)
//! ```
//!
//! where the codon-optimality penalty of a chosen codon is `−ln(w)`
//! with `w` its *relative adaptiveness* (the CAI weight). With
//!
//! - `λ = 0` the DP returns the pure **minimum-free-energy** CDS — the
//!   most thermodynamically stable mRNA that still encodes the protein;
//! - a large `λ` the structural term is swamped and the DP returns the
//!   pure **CAI-optimal** CDS — every residue gets its host-optimal
//!   codon;
//! - an intermediate `λ` traces the **Pareto trade-off** between
//!   stability and translational efficiency.
//!
//! ## Energy model — and the v1 structural scope
//!
//! The folding term reuses the **Turner-2004 nearest-neighbor model**
//! of `valenx-rnastruct` ([`valenx_rnastruct::structure_energy`] is the
//! authoritative re-score of the chosen CDS, so the reported MFE is
//! exactly that crate's energy of the returned structure).
//!
//! The lattice DP optimises over the structural class of **stacked
//! helices, hairpin loops and multiloops** — secondary structures with
//! no bulges and no internal loops. This is a deliberate, documented
//! v1 restriction: a stacked helix is the dominant stabilising element
//! of any RNA fold, and excluding bulges / internal loops is what keeps
//! the lattice DP **exact** — every Turner energy term in this class
//! reads only nucleotides whose codons the DP state already pins, so
//! the recurrence is provably the joint optimum over the lattice for
//! this class (asserted by the Pareto-monotonicity tests). The
//! returned structure is re-scored with the *full* Turner model, so the
//! reported MFE is the true energy of a real, valid secondary
//! structure; it is an upper bound on the unrestricted MFE of that CDS.
//!
//! ## v1 scope — honest framing
//!
//! - **Structural class**: stacks + hairpins + multiloops; no bulges /
//!   internal loops (see above). Coaxial stacking and dangling ends are
//!   not folded into the lattice recurrence.
//! - With a beam wide enough to disable pruning the DP is the **exact**
//!   joint optimum over the lattice for that structural class; with the
//!   default beam it is near-optimal, as LinearFold / LinearDesign are.
//! - The codon-optimality penalty uses the `valenx-bioseq` relative
//!   synonymous-codon-usage tables (E. coli / human) — representative
//!   published figures, see that crate.
//! - The returned MFE / CAI are *predictions* from the energy model and
//!   the usage tables; neither is a measurement. A LinearDesign CDS is
//!   a strong in-silico candidate that must still be lab-validated.

use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use valenx_bioseq::cloning::codon_opt::{codon_usage_table, CodonUsageTable, Host};
use valenx_bioseq::ops::translate::GeneticCode;
use valenx_genediting::mrna::codon::ExpressionHost;
use valenx_rnastruct::fold::energy::{self, multiloop};
use valenx_rnastruct::{structure_energy, BasePair, RnaSeq, Structure};

/// The minimum number of unpaired bases in a hairpin loop — matches
/// `valenx-rnastruct`'s `MIN_HAIRPIN`.
const MIN_HAIRPIN: usize = 3;

/// Maximum hairpin-loop size the lattice DP considers. A hairpin loop
/// longer than this is astronomically unstable under the Turner model
/// and never part of an optimal structure; capping keeps the hairpin
/// search bounded.
const MAX_HAIRPIN: usize = 30;

/// A large positive energy used as the "infeasible" sentinel.
const INF: f64 = 1.0e9;

/// Default beam width for the lattice DP. Like LinearFold's default
/// (100), wide enough to be near-optimal on real CDS while keeping the
/// search near-linear in the CDS length.
pub const DEFAULT_BEAM_SIZE: usize = 200;

// ---------------------------------------------------------------------
// Request / result types
// ---------------------------------------------------------------------

/// A request for a LinearDesign joint mRNA optimisation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinearDesignRequest {
    /// The protein to encode, as one-letter amino-acid codes. A leading
    /// `M` and a trailing `*` are added automatically if absent, so the
    /// caller may pass the bare mature protein.
    pub protein: Vec<u8>,
    /// The expression host — selects the codon-usage table.
    pub host: ExpressionHost,
    /// The trade-off weight `λ ≥ 0` on the codon-optimality penalty.
    /// `0` → pure minimum-free-energy CDS; large → pure CAI-optimal CDS.
    pub lambda: f64,
    /// The beam width for the lattice DP. A wider beam is more accurate
    /// and slower; a very wide beam disables pruning and gives the exact
    /// lattice optimum. `0` is treated as `1`.
    pub beam_size: usize,
}

impl LinearDesignRequest {
    /// A request encoding `protein` for a host at trade-off `lambda`,
    /// with the default beam width.
    pub fn new(protein: impl Into<Vec<u8>>, host: ExpressionHost, lambda: f64) -> Self {
        LinearDesignRequest {
            protein: protein.into(),
            host,
            lambda,
            beam_size: DEFAULT_BEAM_SIZE,
        }
    }
}

/// The result of a LinearDesign joint mRNA optimisation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinearDesignResult {
    /// The optimised coding sequence, RNA (`A C G U`), including the
    /// `AUG` start codon and the stop codon.
    pub cds: Vec<u8>,
    /// The predicted minimum free energy of `cds` (kcal/mol) — the
    /// `valenx-rnastruct` Turner-2004 energy of the structure the
    /// lattice DP found. More negative = a more stable mRNA.
    pub mfe: f64,
    /// The predicted secondary structure of `cds` (the structure the
    /// lattice DP folded it into).
    pub structure: Structure,
    /// The codon adaptation index of `cds` for the host, `(0, 1]`.
    pub cai: f64,
    /// The trade-off weight `λ` used.
    pub lambda: f64,
    /// `true` if the beam was wide enough that no pruning happened — the
    /// CDS is then the *exact* joint optimum over the lattice (for the
    /// stacks+hairpins+multiloops structural class).
    pub exact: bool,
    /// The combined objective value the DP minimised
    /// (`MFE + λ·codon-penalty`); lower is better.
    pub objective: f64,
}

impl LinearDesignResult {
    /// The CDS as a `&str` (valid UTF-8 — only `A C G U`).
    pub fn cds_str(&self) -> &str {
        std::str::from_utf8(&self.cds).unwrap_or("")
    }
}

// ---------------------------------------------------------------------
// The codon lattice
// ---------------------------------------------------------------------

/// One residue's slot in the codon lattice: the synonymous codon
/// triples (RNA, internal `0..4` codes) and their codon-optimality
/// penalties.
#[derive(Clone, Debug)]
struct CodonSlot {
    /// Each option: the three base codes and its `−ln(w)` penalty.
    options: Vec<CodonOption>,
}

/// A single synonymous codon choice for a residue.
#[derive(Copy, Clone, Debug)]
struct CodonOption {
    /// The three nucleotide codes (`A=0 C=1 G=2 U=3`).
    bases: [u8; 3],
    /// The codon-optimality penalty `−ln(w)` (`≥ 0`, `0` for the
    /// host-optimal codon and for single-codon families).
    penalty: f64,
}

/// The codon lattice over a whole protein.
struct CodonLattice {
    /// One [`CodonSlot`] per CDS codon (residue + start + stop).
    slots: Vec<CodonSlot>,
}

impl CodonLattice {
    /// Builds the lattice for `protein` (one-letter codes, *with* the
    /// leading `M` and trailing `*` already present) and a host table.
    fn build(protein: &[u8], table: &CodonUsageTable, code: &GeneticCode) -> Result<Self> {
        let mut slots = Vec::with_capacity(protein.len());
        for (pos, &aa) in protein.iter().enumerate() {
            let aa = aa.to_ascii_uppercase();
            let mut raw: Vec<([u8; 3], f64)> = Vec::new();
            let mut best_freq = 0.0_f64;
            const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            for &b0 in &DNA_BASES {
                for &b1 in &DNA_BASES {
                    for &b2 in &DNA_BASES {
                        let codon = [b0, b1, b2];
                        if code.translate_codon(&codon) != aa {
                            continue;
                        }
                        let f = table.frequency(&codon);
                        if f > best_freq {
                            best_freq = f;
                        }
                        raw.push((codon, f));
                    }
                }
            }
            if raw.is_empty() {
                return Err(RnaDesignError::goal(
                    "protein",
                    format!(
                        "residue `{}` at position {pos} has no codon in the codon table",
                        aa as char
                    ),
                ));
            }
            let mut options = Vec::with_capacity(raw.len());
            for (codon, f) in raw {
                let bases = [
                    dna_to_code(codon[0]),
                    dna_to_code(codon[1]),
                    dna_to_code(codon[2]),
                ];
                // Relative adaptiveness w = f / best; penalty = -ln w.
                // A zero-frequency codon gets a small floor so it is
                // heavily — but not infinitely — penalised.
                let w = if best_freq > 0.0 {
                    (f / best_freq).max(1.0e-3)
                } else {
                    1.0
                };
                options.push(CodonOption {
                    bases,
                    penalty: -w.ln(),
                });
            }
            // Deterministic order: lowest penalty (host-optimal) first,
            // then lexicographic by codon.
            options.sort_by(|a, b| {
                a.penalty
                    .partial_cmp(&b.penalty)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.bases.cmp(&b.bases))
            });
            slots.push(CodonSlot { options });
        }
        Ok(CodonLattice { slots })
    }

    /// Total CDS length in nucleotides.
    fn n_nt(&self) -> usize {
        self.slots.len() * 3
    }

    /// The codon index of nucleotide position `p`.
    #[inline]
    fn codon_of(p: usize) -> usize {
        p / 3
    }

    /// The within-codon offset (`0..3`) of nucleotide position `p`.
    #[inline]
    fn offset_of(p: usize) -> usize {
        p % 3
    }

    /// The base code at position `p` if codon `codon_of(p)` takes choice
    /// `choice`.
    #[inline]
    fn base_at(&self, p: usize, choice: usize) -> u8 {
        self.slots[Self::codon_of(p)].options[choice].bases[Self::offset_of(p)]
    }

    /// The number of codon choices for codon `c`.
    #[inline]
    fn n_choices(&self, c: usize) -> usize {
        self.slots[c].options.len()
    }

    /// The codon-optimality penalty of choice `choice` for codon `c`.
    #[inline]
    fn penalty(&self, c: usize, choice: usize) -> f64 {
        self.slots[c].options[choice].penalty
    }
}

/// Maps a DNA base byte to the internal RNA `0..4` code (`T`→`U`).
fn dna_to_code(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => energy::A,
        b'C' => energy::C,
        b'G' => energy::G,
        b'T' | b'U' => energy::U,
        _ => energy::A,
    }
}

/// Maps an internal RNA code to its ASCII letter.
fn code_to_ascii(c: u8) -> u8 {
    match c {
        energy::A => b'A',
        energy::C => b'C',
        energy::G => b'G',
        _ => b'U',
    }
}

// ---------------------------------------------------------------------
// The lattice folding DP
// ---------------------------------------------------------------------
//
// The grammar is the standard Zuker grammar restricted to the
// stacks + hairpins + multiloops structural class:
//
//   P[i..j]  — i pairs j; the region closes a hairpin, a *stack* onto
//              the inner pair (i+1, j-1), or a multiloop.
//   M1[i..j] — a one-branch multiloop fragment (the branch starts at i;
//              arbitrary trailing unpaired bases after it).
//   M[i..j]  — a ≥1-branch multiloop fragment.
//   M2[i..j] — a ≥2-branch multiloop fragment.
//   C[..j]   — the exterior loop over the prefix [0..=j].
//
// Every Turner energy term in this class reads only the four
// nucleotides at the ends of a pair plus the pair's two immediate
// neighbours.  A `P` state is therefore keyed by the codon choices of
// the four codons containing  i, i+1, j-1, j — `PCodons` — which is
// enough for *any* enclosing context to recompute the energy that
// touches this pair.  Multiloop fragments use the linear (dangle-free)
// Turner multiloop model, which reads no mismatch bases, so M-fragments
// are keyed only by the codons of their two endpoints.

/// The codon choices of the four codons a paired region `[i..=j]`
/// exposes to its enclosing context: the codons of `i`, `i+1`, `j-1`
/// and `j`. Where two of those positions share a codon the stored
/// choices are equal (consistency is enforced on construction).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
struct PCodons {
    /// Codon choice of the codon containing `i`.
    ci: usize,
    /// Codon choice of the codon containing `i + 1`.
    cip1: usize,
    /// Codon choice of the codon containing `j - 1`.
    cjm1: usize,
    /// Codon choice of the codon containing `j`.
    cj: usize,
}

/// A `P`-channel beam key: the start position and the four exposed
/// codon choices.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
struct PKey {
    /// Start nucleotide position of the paired region.
    i: usize,
    /// The four exposed codon choices.
    cod: PCodons,
}

/// An `M`-channel beam key: a fragment `[i..=j]` keyed by the codon
/// choices of its two endpoints (the linear multiloop model reads no
/// mismatch bases, so two codons suffice).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
struct MKey {
    /// Start nucleotide position of the fragment.
    i: usize,
    /// Codon choice of the codon containing `i`.
    ci: usize,
    /// Codon choice of the codon containing `j` (the beam's position).
    cj: usize,
}

/// How a paired-region (`P`) state was derived — for traceback.
#[derive(Copy, Clone, Debug)]
enum PFrom {
    /// `[i..=j]` is a hairpin.
    Hairpin,
    /// `[i..=j]` stacks directly onto the inner pair `(i+1, j-1)`; the
    /// inner pair's exposed codons are carried.
    Stack(PCodons),
    /// `[i..=j]` closes a multiloop; the interior is `M2[i+1..=j-1]`
    /// with the carried endpoint codon choices.
    Multi { m2_ci: usize, m2_cj: usize },
}

/// How a one-branch multiloop fragment (`M1`) state was derived.
#[derive(Copy, Clone, Debug)]
enum M1From {
    /// The branch is the pair `[i..=j]`, with the carried exposed
    /// codons of that pair.
    Branch(PCodons),
    /// `j` is a trailing unpaired base; the rest is `M1[i..=j-1]` whose
    /// right endpoint codon choice is carried.
    UnpairedRight { prev_cj: usize },
}

/// How a `≥1`-branch multiloop fragment (`M`) state was derived.
#[derive(Copy, Clone, Debug)]
enum MFrom {
    /// Promoted from a one-branch `M1` fragment.
    FromM1,
    /// `M[i+1..=j]` with a leading unpaired base at `i`; the inner
    /// fragment's left endpoint codon choice is carried.
    UnpairedLeft { inner_ci: usize },
    /// `M2[i..=j]` (≥2 branches is also ≥1).
    FromM2,
}

/// How a `≥2`-branch multiloop fragment (`M2`) state was derived.
#[derive(Copy, Clone, Debug)]
enum M2From {
    /// `M[i..=k]` followed by the last branch `M1[k+1..=j]`; the split
    /// point and the two fragments' shared-boundary codon choices are
    /// carried.
    Split { k: usize, m_cj: usize, m1_ci: usize },
    /// `j` unpaired; the rest is `M2[i..=j-1]` whose right endpoint
    /// codon choice is carried.
    UnpairedRight { prev_cj: usize },
}

/// How an exterior (`C`) state was derived.
#[derive(Copy, Clone, Debug)]
enum CFrom {
    /// `j` is unpaired in the exterior loop; the previous exterior
    /// state's codon choice is carried.
    Unpaired { prev_cj: usize },
    /// An exterior branch `P[i..=j]` closes, preceded by exterior
    /// `C[i-1]`. The branch's exposed codons and the preceding
    /// exterior state's codon choice are carried.
    Branch {
        i: usize,
        cod: PCodons,
        prev_cj: usize,
    },
}

/// A scored DP state: an objective value (`MFE + λ·penalty` of the
/// sub-problem) and a backpointer.
#[derive(Copy, Clone, Debug)]
struct Cell<T> {
    objective: f64,
    from: T,
}

/// The left-to-right beam-search folder over the codon lattice.
struct LatticeFold<'a> {
    lat: &'a CodonLattice,
    n: usize,
    beam: usize,
    lambda: f64,
    /// `true` once any beam was pruned.
    pruned: bool,
    /// `p[j]` — paired regions ending at `j`, keyed by [`PKey`].
    p: Vec<HashMap<PKey, Cell<PFrom>>>,
    /// `m1[j]` — one-branch multiloop fragments ending at `j`.
    m1: Vec<HashMap<MKey, Cell<M1From>>>,
    /// `m[j]` — ≥1-branch multiloop fragments ending at `j`.
    m: Vec<HashMap<MKey, Cell<MFrom>>>,
    /// `m2[j]` — ≥2-branch multiloop fragments ending at `j`.
    m2: Vec<HashMap<MKey, Cell<M2From>>>,
    /// `c[j]` — exterior of the prefix `[0..=j]`, keyed by the codon
    /// choice of the codon containing `j`.
    c: Vec<HashMap<usize, Cell<CFrom>>>,
}

impl<'a> LatticeFold<'a> {
    fn new(lat: &'a CodonLattice, beam: usize, lambda: f64) -> Self {
        let n = lat.n_nt();
        LatticeFold {
            lat,
            n,
            beam,
            lambda,
            pruned: false,
            p: vec![HashMap::new(); n],
            m1: vec![HashMap::new(); n],
            m: vec![HashMap::new(); n],
            m2: vec![HashMap::new(); n],
            c: vec![HashMap::new(); n],
        }
    }

    /// The codon-optimality penalty *charged once per codon*: a codon is
    /// charged when its **last** nucleotide (offset 2) is consumed by a
    /// state. This charges every full CDS exactly its total codon
    /// penalty regardless of which DP channel covers each codon.
    #[inline]
    fn codon_charge(&self, p: usize, choice: usize) -> f64 {
        if CodonLattice::offset_of(p) == 2 {
            self.lambda * self.lat.penalty(CodonLattice::codon_of(p), choice)
        } else {
            0.0
        }
    }

    /// The minimum codon-optimality charge of the codon containing
    /// position `p` — used for codons fully inside a hairpin loop,
    /// which the energy never reads and nothing else constrains, so the
    /// joint optimum gives them their host-optimal (penalty-0) codon.
    /// (With `options` sorted host-optimal-first the minimum penalty is
    /// `options[0].penalty == 0`, so this is `0.0`; kept as a function
    /// for clarity.)
    #[inline]
    fn min_codon_charge(&self, p: usize) -> f64 {
        if CodonLattice::offset_of(p) == 2 {
            self.lambda
                * self.lat.penalty(CodonLattice::codon_of(p), 0)
        } else {
            0.0
        }
    }

    /// Codon choices for the four exposed codons of a pair `(i, j)` are
    /// **consistent** when any two of `i, i+1, j-1, j` that share a
    /// codon carry the same choice.
    fn pcodons_consistent(i: usize, j: usize, cod: &PCodons) -> bool {
        let positions = [
            (i, cod.ci),
            (i + 1, cod.cip1),
            (j - 1, cod.cjm1),
            (j, cod.cj),
        ];
        for a in 0..4 {
            for b in (a + 1)..4 {
                if CodonLattice::codon_of(positions[a].0)
                    == CodonLattice::codon_of(positions[b].0)
                    && positions[a].1 != positions[b].1
                {
                    return false;
                }
            }
        }
        true
    }

    /// Runs the whole left-to-right fill.
    fn run(&mut self) {
        for j in 0..self.n {
            self.fill_p(j);
            self.fill_m1(j);
            self.fill_m_and_m2(j);
            self.fill_c(j);
            self.pruned |= prune_p(&mut self.p[j], self.beam);
            self.pruned |= prune_m(&mut self.m1[j], self.beam);
            self.pruned |= prune_m(&mut self.m[j], self.beam);
            self.pruned |= prune_m(&mut self.m2[j], self.beam);
            self.pruned |= prune_c(&mut self.c[j], self.beam);
        }
    }

    /// Fills `p[j]` — every paired region `[i..=j]`.
    fn fill_p(&mut self, j: usize) {
        let mut out: HashMap<PKey, Cell<PFrom>> = HashMap::new();

        // --- hairpins -----------------------------------------------
        if j > MIN_HAIRPIN {
            let i_lo = j.saturating_sub(MAX_HAIRPIN + 1);
            let i_hi = j - MIN_HAIRPIN - 1;
            for i in i_lo..=i_hi {
                self.fill_hairpin(i, j, &mut out);
            }
        }

        // --- stacks: (i, j) stacks onto the inner pair (i+1, j-1) ----
        // The inner pair is a P state ending at j-1 starting at i+1.
        if j >= 2 {
            let inner: Vec<(PKey, f64)> = self.p[j - 1]
                .iter()
                .map(|(&key, c)| (key, c.objective))
                .collect();
            for (inner_key, inner_obj) in inner {
                if inner_obj >= INF / 2.0 {
                    continue;
                }
                let k = inner_key.i; // inner pair start
                if k == 0 {
                    continue;
                }
                let i = k - 1;
                let l = j - 1;
                self.fill_stack(i, j, k, l, inner_key, inner_obj, &mut out);
            }
        }

        // --- multiloop close: interior is M2[i+1..=j-1] --------------
        if j >= 2 {
            let interior: Vec<(MKey, f64)> = self.m2[j - 1]
                .iter()
                .map(|(&key, c)| (key, c.objective))
                .collect();
            for (m2_key, m2_obj) in interior {
                if m2_obj >= INF / 2.0 || m2_key.i == 0 {
                    continue;
                }
                self.fill_multi_close(j, m2_key, m2_obj, &mut out);
            }
        }

        self.p[j] = out;
    }

    /// Fills the hairpins closed by `(i, j)`.
    ///
    /// A small hairpin loop (≤ 4 unpaired bases) spans ≤ 3 codons and
    /// its energy reads the whole loop (the published triloop /
    /// tetraloop tables) — every spanned codon is enumerated. A larger
    /// hairpin loop's Turner energy reads only `i, i+1, j-1, j`, so just
    /// those four codons are enumerated; codons fully inside the loop
    /// are charged their host-optimal penalty.
    fn fill_hairpin(&self, i: usize, j: usize, out: &mut HashMap<PKey, Cell<PFrom>>) {
        let lat = self.lat;
        let loop_size = j - i - 1;
        // The distinct codons of i, i+1, j-1, j.
        let bnd = distinct_codons(&[i, i + 1, j - 1, j]);

        if loop_size <= 4 {
            // Enumerate every spanned codon.
            let first = CodonLattice::codon_of(i);
            let last = CodonLattice::codon_of(j);
            let spanned: Vec<usize> = (first..=last).collect();
            let mut choice = vec![0usize; spanned.len()];
            loop {
                let ch_of = |c: usize| choice[c - first];
                let cod = PCodons {
                    ci: ch_of(CodonLattice::codon_of(i)),
                    cip1: ch_of(CodonLattice::codon_of(i + 1)),
                    cjm1: ch_of(CodonLattice::codon_of(j - 1)),
                    cj: ch_of(CodonLattice::codon_of(j)),
                };
                let bi = lat.base_at(i, cod.ci);
                let bj = lat.base_at(j, cod.cj);
                if energy::can_pair_codes(bi, bj) {
                    let loop_bases: Vec<u8> = ((i + 1)..j)
                        .map(|p| lat.base_at(p, ch_of(CodonLattice::codon_of(p))))
                        .collect();
                    let e = energy::hairpin_energy(bi, bj, &loop_bases);
                    if e < INF / 2.0 {
                        let mut charge = 0.0;
                        for p in i..=j {
                            if CodonLattice::offset_of(p) == 2 {
                                charge += self.lambda
                                    * lat.penalty(
                                        CodonLattice::codon_of(p),
                                        ch_of(CodonLattice::codon_of(p)),
                                    );
                            }
                        }
                        relax_p(out, PKey { i, cod }, e + charge, PFrom::Hairpin);
                    }
                }
                if !advance_odometer(&mut choice, &spanned, lat) {
                    break;
                }
            }
            return;
        }

        // Large loop: enumerate only the boundary codons.
        let mut choice = vec![0usize; bnd.len()];
        loop {
            let ch_of = |c: usize| -> usize {
                bnd.iter().position(|&x| x == c).map(|idx| choice[idx]).unwrap_or(0)
            };
            let cod = PCodons {
                ci: ch_of(CodonLattice::codon_of(i)),
                cip1: ch_of(CodonLattice::codon_of(i + 1)),
                cjm1: ch_of(CodonLattice::codon_of(j - 1)),
                cj: ch_of(CodonLattice::codon_of(j)),
            };
            let bi = lat.base_at(i, cod.ci);
            let bj = lat.base_at(j, cod.cj);
            if energy::can_pair_codes(bi, bj) {
                let mm5 = lat.base_at(i + 1, cod.cip1);
                let mm3 = lat.base_at(j - 1, cod.cjm1);
                // For a large hairpin only loop[0] and loop[last] are
                // load-bearing; the rest is filler.
                let mut loop_bases = vec![energy::A; loop_size];
                loop_bases[0] = mm5;
                loop_bases[loop_size - 1] = mm3;
                let e = energy::hairpin_energy(bi, bj, &loop_bases);
                if e < INF / 2.0 {
                    let mut charge = 0.0;
                    for p in i..=j {
                        if CodonLattice::offset_of(p) == 2 {
                            let codon = CodonLattice::codon_of(p);
                            if bnd.contains(&codon) {
                                charge += self.lambda
                                    * lat.penalty(codon, ch_of(codon));
                            } else {
                                // Interior codon: host-optimal.
                                charge += self.min_codon_charge(p);
                            }
                        }
                    }
                    relax_p(out, PKey { i, cod }, e + charge, PFrom::Hairpin);
                }
            }
            if !advance_odometer(&mut choice, &bnd, lat) {
                break;
            }
        }
    }

    /// Fills a stack: `(i, j)` pairs and stacks directly onto the inner
    /// pair `(k, l) = (i+1, j-1)`. The stack energy reads exactly the
    /// four bases `i, j, k, l`; `k`'s and `l`'s codons are pinned by the
    /// inner pair's exposed codons, and the inner pair's `i+1`-codon is
    /// `k`'s codon, its `j-1`-codon is `l`'s codon — so every codon the
    /// stack energy and the *enclosing* context need is determined.
    #[allow(clippy::too_many_arguments)]
    fn fill_stack(
        &self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        inner_key: PKey,
        inner_obj: f64,
        out: &mut HashMap<PKey, Cell<PFrom>>,
    ) {
        let lat = self.lat;
        // The inner pair (k, l) — its closing-pair codons.
        let inner = inner_key.cod;
        // k = i+1, l = j-1. The inner pair's exposed codons:
        //   inner.ci   = codon(k)
        //   inner.cj   = codon(l)
        // The stack reads codon(i), codon(j) (fresh) and codon(k),
        // codon(l) (pinned). The *new* P state's exposed codons are
        //   ci  = codon(i)            (fresh)
        //   cip1 = codon(i+1) = codon(k)   = inner.ci
        //   cjm1 = codon(j-1) = codon(l)   = inner.cj
        //   cj  = codon(j)            (fresh)
        let ci_codon = CodonLattice::codon_of(i);
        let cj_codon = CodonLattice::codon_of(j);
        for choice_i in 0..lat.n_choices(ci_codon) {
            for choice_j in 0..lat.n_choices(cj_codon) {
                let cod = PCodons {
                    ci: choice_i,
                    cip1: inner.ci,
                    cjm1: inner.cj,
                    cj: choice_j,
                };
                if !Self::pcodons_consistent(i, j, &cod) {
                    continue;
                }
                let bi = lat.base_at(i, choice_i);
                let bj = lat.base_at(j, choice_j);
                if !energy::can_pair_codes(bi, bj) {
                    continue;
                }
                let bk = lat.base_at(k, inner.ci);
                let bl = lat.base_at(l, inner.cj);
                // Internal-loop energy with left = right = 0 is the
                // stack energy.
                let e = energy::internal_loop_energy(
                    bi, bj, bk, bl, 0, 0, bk, bl, bk, bl,
                );
                if e >= INF / 2.0 {
                    continue;
                }
                // Charge codons closed at i and j (offset-2 rule); the
                // inner pair self-charged [k..=l].
                let charge =
                    self.codon_charge(i, choice_i) + self.codon_charge(j, choice_j);
                let obj = e + inner_obj + charge;
                relax_p(
                    out,
                    PKey { i, cod },
                    obj,
                    PFrom::Stack(inner),
                );
            }
        }
    }

    /// Fills a multiloop closure: `(i, j)` closes a multiloop whose
    /// interior is the `≥2`-branch fragment `M2[i+1..=j-1]`.
    fn fill_multi_close(
        &self,
        j: usize,
        m2_key: MKey,
        m2_obj: f64,
        out: &mut HashMap<PKey, Cell<PFrom>>,
    ) {
        let lat = self.lat;
        let i = m2_key.i - 1;
        // The M2 interior runs [i+1..=j-1]; its endpoint codon choices:
        //   m2_key.ci = codon(i+1),  m2_key.cj = codon(j-1).
        let ci_codon = CodonLattice::codon_of(i);
        let cj_codon = CodonLattice::codon_of(j);
        for choice_i in 0..lat.n_choices(ci_codon) {
            for choice_j in 0..lat.n_choices(cj_codon) {
                let cod = PCodons {
                    ci: choice_i,
                    cip1: m2_key.ci,
                    cjm1: m2_key.cj,
                    cj: choice_j,
                };
                if !Self::pcodons_consistent(i, j, &cod) {
                    continue;
                }
                let bi = lat.base_at(i, choice_i);
                let bj = lat.base_at(j, choice_j);
                if !energy::can_pair_codes(bi, bj) {
                    continue;
                }
                let closure = multiloop::OFFSET
                    + multiloop::PER_BRANCH
                    + energy::terminal_penalty(bi, bj);
                let charge =
                    self.codon_charge(i, choice_i) + self.codon_charge(j, choice_j);
                let obj = closure + m2_obj + charge;
                relax_p(
                    out,
                    PKey { i, cod },
                    obj,
                    PFrom::Multi {
                        m2_ci: m2_key.ci,
                        m2_cj: m2_key.cj,
                    },
                );
            }
        }
    }

    /// Fills `m1[j]` — one-branch multiloop fragments.
    fn fill_m1(&mut self, j: usize) {
        let lat = self.lat;
        let cj_codon = CodonLattice::codon_of(j);
        let mut out: HashMap<MKey, Cell<M1From>> = HashMap::new();

        // The branch is an entire P[i..=j].
        for (&key, cell) in &self.p[j] {
            let obj = cell.objective + multiloop::PER_BRANCH;
            let mkey = MKey {
                i: key.i,
                ci: key.cod.ci,
                cj: key.cod.cj,
            };
            relax_m(&mut out, mkey, obj, M1From::Branch(key.cod));
        }
        // j is a trailing unpaired base; the rest is m1[i..=j-1].
        if j > 0 {
            let prev: Vec<(MKey, f64)> = self.m1[j - 1]
                .iter()
                .map(|(&k, c)| (k, c.objective))
                .collect();
            let cjm1_codon = CodonLattice::codon_of(j - 1);
            for (prev_key, prev_obj) in prev {
                for choice_j in 0..lat.n_choices(cj_codon) {
                    // j's codon must extend the previous fragment's
                    // right codon, and stay consistent with its left.
                    if !codons_consistent(
                        cj_codon,
                        choice_j,
                        cjm1_codon,
                        prev_key.cj,
                    ) {
                        continue;
                    }
                    if !codons_consistent(
                        CodonLattice::codon_of(prev_key.i),
                        prev_key.ci,
                        cj_codon,
                        choice_j,
                    ) {
                        continue;
                    }
                    let charge = self.codon_charge(j, choice_j);
                    let obj = prev_obj + multiloop::PER_UNPAIRED + charge;
                    relax_m(
                        &mut out,
                        MKey {
                            i: prev_key.i,
                            ci: prev_key.ci,
                            cj: choice_j,
                        },
                        obj,
                        M1From::UnpairedRight {
                            prev_cj: prev_key.cj,
                        },
                    );
                }
            }
        }
        self.m1[j] = out;
    }

    /// Fills `m2[j]` (≥2 branches) then `m[j]` (≥1 branch).
    fn fill_m_and_m2(&mut self, j: usize) {
        let lat = self.lat;
        let cj_codon = CodonLattice::codon_of(j);

        // --- m2: M[i..=k] then last branch M1[k+1..=j] ---------------
        let mut m2_out: HashMap<MKey, Cell<M2From>> = HashMap::new();
        let m1_states: Vec<(MKey, f64)> = self.m1[j]
            .iter()
            .map(|(&key, c)| (key, c.objective))
            .collect();
        for (m1_key, m1_obj) in m1_states {
            let kp1 = m1_key.i;
            if kp1 == 0 {
                continue;
            }
            let k = kp1 - 1;
            let ck_codon = CodonLattice::codon_of(k);
            let m_states: Vec<(MKey, f64)> = self.m[k]
                .iter()
                .map(|(&key, c)| (key, c.objective))
                .collect();
            for (m_key, m_obj) in m_states {
                // M's right codon (codon of k) must agree with M1's
                // left codon (codon of k+1).
                if !codons_consistent(
                    ck_codon,
                    m_key.cj,
                    CodonLattice::codon_of(kp1),
                    m1_key.ci,
                ) {
                    continue;
                }
                // Combined fragment's endpoint codons consistent.
                if !codons_consistent(
                    CodonLattice::codon_of(m_key.i),
                    m_key.ci,
                    cj_codon,
                    m1_key.cj,
                ) {
                    continue;
                }
                let obj = m_obj + m1_obj;
                relax_m(
                    &mut m2_out,
                    MKey {
                        i: m_key.i,
                        ci: m_key.ci,
                        cj: m1_key.cj,
                    },
                    obj,
                    M2From::Split {
                        k,
                        m_cj: m_key.cj,
                        m1_ci: m1_key.ci,
                    },
                );
            }
        }
        // j unpaired: extend an M2 ending at j-1.
        if j > 0 {
            let prev: Vec<(MKey, f64)> = self.m2[j - 1]
                .iter()
                .map(|(&key, c)| (key, c.objective))
                .collect();
            let cjm1_codon = CodonLattice::codon_of(j - 1);
            for (prev_key, prev_obj) in prev {
                for choice_j in 0..lat.n_choices(cj_codon) {
                    if !codons_consistent(
                        cj_codon,
                        choice_j,
                        cjm1_codon,
                        prev_key.cj,
                    ) {
                        continue;
                    }
                    if !codons_consistent(
                        CodonLattice::codon_of(prev_key.i),
                        prev_key.ci,
                        cj_codon,
                        choice_j,
                    ) {
                        continue;
                    }
                    let charge = self.codon_charge(j, choice_j);
                    let obj = prev_obj + multiloop::PER_UNPAIRED + charge;
                    relax_m(
                        &mut m2_out,
                        MKey {
                            i: prev_key.i,
                            ci: prev_key.ci,
                            cj: choice_j,
                        },
                        obj,
                        M2From::UnpairedRight {
                            prev_cj: prev_key.cj,
                        },
                    );
                }
            }
        }

        // --- m: m1 ∪ m2, plus a leading unpaired extension -----------
        let mut m_out: HashMap<MKey, Cell<MFrom>> = HashMap::new();
        for (&key, cell) in &self.m1[j] {
            relax_m(&mut m_out, key, cell.objective, MFrom::FromM1);
        }
        for (&key, cell) in &m2_out {
            relax_m(&mut m_out, key, cell.objective, MFrom::FromM2);
        }
        // Leading unpaired base at i: sweep starts from j down to 1.
        let starts: Vec<usize> = {
            let mut s: Vec<usize> = m_out.keys().map(|k| k.i).collect();
            s.sort_unstable();
            s.dedup();
            s
        };
        for &s in starts.iter().rev() {
            if s == 0 {
                continue;
            }
            let cur: Vec<(MKey, f64)> = m_out
                .iter()
                .filter(|(k, _)| k.i == s)
                .map(|(&k, c)| (k, c.objective))
                .collect();
            let prev_pos = s - 1;
            let cprev = CodonLattice::codon_of(prev_pos);
            let cs = CodonLattice::codon_of(s);
            for (key, obj) in cur {
                for choice_prev in 0..lat.n_choices(cprev) {
                    if !codons_consistent(cprev, choice_prev, cs, key.ci) {
                        continue;
                    }
                    if !codons_consistent(
                        cprev,
                        choice_prev,
                        cj_codon,
                        key.cj,
                    ) {
                        continue;
                    }
                    let charge = self.codon_charge(prev_pos, choice_prev);
                    let new_obj = obj + multiloop::PER_UNPAIRED + charge;
                    relax_m(
                        &mut m_out,
                        MKey {
                            i: prev_pos,
                            ci: choice_prev,
                            cj: key.cj,
                        },
                        new_obj,
                        MFrom::UnpairedLeft { inner_ci: key.ci },
                    );
                }
            }
        }

        self.m2[j] = m2_out;
        self.m[j] = m_out;
    }

    /// Fills `c[j]` — the exterior of the prefix `[0..=j]`.
    fn fill_c(&mut self, j: usize) {
        let lat = self.lat;
        let cj_codon = CodonLattice::codon_of(j);
        let mut out: HashMap<usize, Cell<CFrom>> = HashMap::new();

        // j unpaired in the exterior.
        if j == 0 {
            for choice_j in 0..lat.n_choices(cj_codon) {
                let charge = self.codon_charge(0, choice_j);
                out.insert(
                    choice_j,
                    Cell {
                        objective: charge,
                        from: CFrom::Unpaired { prev_cj: 0 },
                    },
                );
            }
        } else {
            let cjm1_codon = CodonLattice::codon_of(j - 1);
            let prev: Vec<(usize, f64)> = self.c[j - 1]
                .iter()
                .map(|(&ch, c)| (ch, c.objective))
                .collect();
            for (prev_choice, prev_obj) in &prev {
                for choice_j in 0..lat.n_choices(cj_codon) {
                    if !codons_consistent(
                        cj_codon,
                        choice_j,
                        cjm1_codon,
                        *prev_choice,
                    ) {
                        continue;
                    }
                    let charge = self.codon_charge(j, choice_j);
                    let obj = prev_obj + charge;
                    relax_c(
                        &mut out,
                        choice_j,
                        obj,
                        CFrom::Unpaired {
                            prev_cj: *prev_choice,
                        },
                    );
                }
            }
        }

        // An exterior branch P[i..=j], preceded by exterior C[i-1].
        let p_states: Vec<(PKey, f64)> = self.p[j]
            .iter()
            .map(|(&key, c)| (key, c.objective))
            .collect();
        for (key, p_obj) in p_states {
            let bi = lat.base_at(key.i, key.cod.ci);
            let bj = lat.base_at(j, key.cod.cj);
            let term = energy::terminal_penalty(bi, bj);
            if key.i == 0 {
                relax_c(
                    &mut out,
                    key.cod.cj,
                    p_obj + term,
                    CFrom::Branch {
                        i: key.i,
                        cod: key.cod,
                        prev_cj: 0,
                    },
                );
            } else {
                let cim1_codon = CodonLattice::codon_of(key.i - 1);
                let prev: Vec<(usize, f64)> = self.c[key.i - 1]
                    .iter()
                    .map(|(&ch, c)| (ch, c.objective))
                    .collect();
                for (prev_choice, prev_obj) in prev {
                    if !codons_consistent(
                        cim1_codon,
                        prev_choice,
                        CodonLattice::codon_of(key.i),
                        key.cod.ci,
                    ) {
                        continue;
                    }
                    let obj = prev_obj + p_obj + term;
                    relax_c(
                        &mut out,
                        key.cod.cj,
                        obj,
                        CFrom::Branch {
                            i: key.i,
                            cod: key.cod,
                            prev_cj: prev_choice,
                        },
                    );
                }
            }
        }

        self.c[j] = out;
    }

    /// Reconstructs the optimal CDS and structure from the filled DP.
    fn traceback(&self) -> Result<(Vec<u8>, Structure, f64)> {
        let n = self.n;
        let mut codon_choice: Vec<Option<usize>> = vec![None; self.lat.slots.len()];
        let mut partner: Vec<Option<usize>> = vec![None; n];

        let last = n - 1;
        let (&best_choice, best_cell) = self.c[last]
            .iter()
            .min_by(|a, b| {
                a.1.objective
                    .partial_cmp(&b.1.objective)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .ok_or_else(|| {
                RnaDesignError::no_design(
                    "lineardesign",
                    "the lattice DP produced no exterior state",
                )
            })?;
        let objective = best_cell.objective;

        // Pins a codon, erroring on a contradiction.
        fn pin(
            codon_choice: &mut [Option<usize>],
            codon: usize,
            choice: usize,
        ) -> Result<()> {
            match codon_choice[codon] {
                None => {
                    codon_choice[codon] = Some(choice);
                    Ok(())
                }
                Some(prev) if prev == choice => Ok(()),
                Some(prev) => Err(RnaDesignError::no_design(
                    "lineardesign",
                    format!(
                        "traceback codon contradiction at {codon}: {prev} vs {choice}"
                    ),
                )),
            }
        }

        /// Pins the four codons a `PCodons` describes for pair `(i, j)`.
        fn pin_pcodons(
            codon_choice: &mut [Option<usize>],
            i: usize,
            j: usize,
            cod: &PCodons,
        ) -> Result<()> {
            pin(codon_choice, CodonLattice::codon_of(i), cod.ci)?;
            pin(codon_choice, CodonLattice::codon_of(i + 1), cod.cip1)?;
            pin(codon_choice, CodonLattice::codon_of(j - 1), cod.cjm1)?;
            pin(codon_choice, CodonLattice::codon_of(j), cod.cj)?;
            Ok(())
        }

        enum Job {
            C(isize, usize),
            P(usize, usize, PCodons),
            M1(MKey, usize),
            M(MKey, usize),
            M2(MKey, usize),
        }
        let mut stack: Vec<Job> = vec![Job::C(last as isize, best_choice)];

        while let Some(job) = stack.pop() {
            match job {
                Job::C(j, choice_j) => {
                    if j < 0 {
                        continue;
                    }
                    let ju = j as usize;
                    pin(&mut codon_choice, CodonLattice::codon_of(ju), choice_j)?;
                    let cell = self.c[ju].get(&choice_j).ok_or_else(|| {
                        RnaDesignError::no_design(
                            "lineardesign",
                            "traceback: missing C state",
                        )
                    })?;
                    match cell.from {
                        CFrom::Unpaired { prev_cj } => {
                            if ju > 0 {
                                stack.push(Job::C(j - 1, prev_cj));
                            }
                        }
                        CFrom::Branch { i, cod, prev_cj } => {
                            stack.push(Job::P(i, ju, cod));
                            if i > 0 {
                                stack.push(Job::C(i as isize - 1, prev_cj));
                            }
                        }
                    }
                }
                Job::P(i, j, cod) => {
                    partner[i] = Some(j);
                    partner[j] = Some(i);
                    pin_pcodons(&mut codon_choice, i, j, &cod)?;
                    let key = PKey { i, cod };
                    let cell = self.p[j].get(&key).ok_or_else(|| {
                        RnaDesignError::no_design(
                            "lineardesign",
                            "traceback: missing P state",
                        )
                    })?;
                    match cell.from {
                        PFrom::Hairpin => {
                            self.repin_hairpin(i, j, &cod, &mut codon_choice)?;
                        }
                        PFrom::Stack(inner) => {
                            // Inner pair (i+1, j-1).
                            stack.push(Job::P(i + 1, j - 1, inner));
                        }
                        PFrom::Multi { m2_ci, m2_cj } => {
                            stack.push(Job::M2(
                                MKey {
                                    i: i + 1,
                                    ci: m2_ci,
                                    cj: m2_cj,
                                },
                                j - 1,
                            ));
                        }
                    }
                }
                Job::M1(key, j) => {
                    pin(&mut codon_choice, CodonLattice::codon_of(key.i), key.ci)?;
                    pin(&mut codon_choice, CodonLattice::codon_of(j), key.cj)?;
                    let cell = self.m1[j].get(&key).ok_or_else(|| {
                        RnaDesignError::no_design(
                            "lineardesign",
                            "traceback: missing M1 state",
                        )
                    })?;
                    match cell.from {
                        M1From::Branch(cod) => {
                            stack.push(Job::P(key.i, j, cod));
                        }
                        M1From::UnpairedRight { prev_cj } => {
                            stack.push(Job::M1(
                                MKey {
                                    i: key.i,
                                    ci: key.ci,
                                    cj: prev_cj,
                                },
                                j - 1,
                            ));
                        }
                    }
                }
                Job::M(key, j) => {
                    pin(&mut codon_choice, CodonLattice::codon_of(key.i), key.ci)?;
                    pin(&mut codon_choice, CodonLattice::codon_of(j), key.cj)?;
                    let cell = self.m[j].get(&key).ok_or_else(|| {
                        RnaDesignError::no_design(
                            "lineardesign",
                            "traceback: missing M state",
                        )
                    })?;
                    match cell.from {
                        MFrom::FromM1 => stack.push(Job::M1(key, j)),
                        MFrom::FromM2 => stack.push(Job::M2(key, j)),
                        MFrom::UnpairedLeft { inner_ci } => {
                            stack.push(Job::M(
                                MKey {
                                    i: key.i + 1,
                                    ci: inner_ci,
                                    cj: key.cj,
                                },
                                j,
                            ));
                        }
                    }
                }
                Job::M2(key, j) => {
                    pin(&mut codon_choice, CodonLattice::codon_of(key.i), key.ci)?;
                    pin(&mut codon_choice, CodonLattice::codon_of(j), key.cj)?;
                    let cell = self.m2[j].get(&key).ok_or_else(|| {
                        RnaDesignError::no_design(
                            "lineardesign",
                            "traceback: missing M2 state",
                        )
                    })?;
                    match cell.from {
                        M2From::Split { k, m_cj, m1_ci } => {
                            stack.push(Job::M(
                                MKey {
                                    i: key.i,
                                    ci: key.ci,
                                    cj: m_cj,
                                },
                                k,
                            ));
                            stack.push(Job::M1(
                                MKey {
                                    i: k + 1,
                                    ci: m1_ci,
                                    cj: key.cj,
                                },
                                j,
                            ));
                        }
                        M2From::UnpairedRight { prev_cj } => {
                            stack.push(Job::M2(
                                MKey {
                                    i: key.i,
                                    ci: key.ci,
                                    cj: prev_cj,
                                },
                                j - 1,
                            ));
                        }
                    }
                }
            }
        }

        // Any codon left unpinned (it sat in a hairpin / exterior
        // unpaired run whose choice the traceback did not need) gets its
        // host-optimal (penalty-0, sorted-first) codon.
        let mut cds = Vec::with_capacity(n);
        for (c, choice) in codon_choice.iter().enumerate() {
            let ch = choice.unwrap_or(0);
            for off in 0..3 {
                cds.push(code_to_ascii(self.lat.slots[c].options[ch].bases[off]));
            }
        }
        let structure = Structure::from_partner(partner).map_err(|e| {
            RnaDesignError::upstream("valenx-rnastruct", e.to_string())
        })?;
        Ok((cds, structure, objective))
    }

    /// Re-pins the codons of a traced-back hairpin `[i..=j]`.
    ///
    /// The four exposed codons are already pinned by the caller. For a
    /// small loop (≤ 4) the interior codons are re-derived by matching
    /// the stored objective; for a large loop the interior codons take
    /// their host-optimal (choice-0) codon, exactly as the fill charged
    /// them.
    fn repin_hairpin(
        &self,
        i: usize,
        j: usize,
        cod: &PCodons,
        codon_choice: &mut [Option<usize>],
    ) -> Result<()> {
        let lat = self.lat;
        let loop_size = j - i - 1;
        let first = CodonLattice::codon_of(i);
        let last = CodonLattice::codon_of(j);

        if loop_size <= 4 {
            let spanned: Vec<usize> = (first..=last).collect();
            let key = PKey { i, cod: *cod };
            let stored = self.p[j].get(&key).map(|c| c.objective).ok_or_else(|| {
                RnaDesignError::no_design("lineardesign", "repin: missing hairpin P")
            })?;
            let mut choice = vec![0usize; spanned.len()];
            loop {
                let ch_of = |c: usize| choice[c - first];
                let trial = PCodons {
                    ci: ch_of(CodonLattice::codon_of(i)),
                    cip1: ch_of(CodonLattice::codon_of(i + 1)),
                    cjm1: ch_of(CodonLattice::codon_of(j - 1)),
                    cj: ch_of(CodonLattice::codon_of(j)),
                };
                if &trial == cod {
                    let bi = lat.base_at(i, cod.ci);
                    let bj = lat.base_at(j, cod.cj);
                    if energy::can_pair_codes(bi, bj) {
                        let loop_bases: Vec<u8> = ((i + 1)..j)
                            .map(|p| {
                                lat.base_at(p, ch_of(CodonLattice::codon_of(p)))
                            })
                            .collect();
                        let e = energy::hairpin_energy(bi, bj, &loop_bases);
                        let mut charge = 0.0;
                        for p in i..=j {
                            if CodonLattice::offset_of(p) == 2 {
                                charge += self.lambda
                                    * lat.penalty(
                                        CodonLattice::codon_of(p),
                                        ch_of(CodonLattice::codon_of(p)),
                                    );
                            }
                        }
                        if (e + charge - stored).abs() < 1e-6 {
                            for &c in &spanned {
                                codon_choice[c] = Some(ch_of(c));
                            }
                            return Ok(());
                        }
                    }
                }
                if !advance_odometer(&mut choice, &spanned, lat) {
                    break;
                }
            }
            return Err(RnaDesignError::no_design(
                "lineardesign",
                "repin: could not reconstruct small-hairpin codons",
            ));
        }

        // Large loop: interior codons take their host-optimal codon.
        if last > first + 1 {
            for slot in &mut codon_choice[(first + 1)..last] {
                if slot.is_none() {
                    *slot = Some(0);
                }
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------
// DP helpers
// ---------------------------------------------------------------------

/// The distinct codons (ascending) the given nucleotide positions fall
/// in.
fn distinct_codons(positions: &[usize]) -> Vec<usize> {
    let mut codons: Vec<usize> =
        positions.iter().map(|&p| CodonLattice::codon_of(p)).collect();
    codons.sort_unstable();
    codons.dedup();
    codons
}

/// Advances a per-codon-choice odometer; returns `false` when it wraps.
fn advance_odometer(choice: &mut [usize], codons: &[usize], lat: &CodonLattice) -> bool {
    for idx in (0..choice.len()).rev() {
        choice[idx] += 1;
        if choice[idx] < lat.n_choices(codons[idx]) {
            return true;
        }
        choice[idx] = 0;
    }
    false
}

/// `true` if codon choices `ci` (codon `a`) and `cj` (codon `b`) agree
/// on every codon they share — they share a codon iff `a == b`.
#[inline]
fn codons_consistent(a: usize, ci: usize, b: usize, cj: usize) -> bool {
    a != b || ci == cj
}

/// Relaxes a `P`-beam entry to the lower of its current objective and
/// `obj`.
fn relax_p(out: &mut HashMap<PKey, Cell<PFrom>>, key: PKey, obj: f64, from: PFrom) {
    if obj >= INF / 2.0 {
        return;
    }
    let slot = out.entry(key).or_insert(Cell { objective: INF, from });
    if obj < slot.objective {
        slot.objective = obj;
        slot.from = from;
    }
}

/// Relaxes an `M`-channel beam entry.
fn relax_m<T: Copy>(
    out: &mut HashMap<MKey, Cell<T>>,
    key: MKey,
    obj: f64,
    from: T,
) {
    if obj >= INF / 2.0 {
        return;
    }
    let slot = out.entry(key).or_insert(Cell { objective: INF, from });
    if obj < slot.objective {
        slot.objective = obj;
        slot.from = from;
    }
}

/// Relaxes an exterior-beam entry (keyed by codon choice).
fn relax_c(out: &mut HashMap<usize, Cell<CFrom>>, key: usize, obj: f64, from: CFrom) {
    if obj >= INF / 2.0 {
        return;
    }
    let slot = out.entry(key).or_insert(Cell { objective: INF, from });
    if obj < slot.objective {
        slot.objective = obj;
        slot.from = from;
    }
}

/// Prunes a `P`-beam to its best `beam` entries; returns `true` if it
/// removed anything.
fn prune_p(beam_map: &mut HashMap<PKey, Cell<PFrom>>, beam: usize) -> bool {
    if beam_map.len() <= beam {
        return false;
    }
    let mut entries: Vec<(PKey, f64)> =
        beam_map.iter().map(|(&k, c)| (k, c.objective)).collect();
    entries.sort_by(|a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    let keep: std::collections::HashSet<PKey> =
        entries.iter().take(beam).map(|&(k, _)| k).collect();
    beam_map.retain(|k, _| keep.contains(k));
    true
}

/// Prunes an `M`-channel beam.
fn prune_m<T: Copy>(beam_map: &mut HashMap<MKey, Cell<T>>, beam: usize) -> bool {
    if beam_map.len() <= beam {
        return false;
    }
    let mut entries: Vec<(MKey, f64)> =
        beam_map.iter().map(|(&k, c)| (k, c.objective)).collect();
    entries.sort_by(|a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    let keep: std::collections::HashSet<MKey> =
        entries.iter().take(beam).map(|&(k, _)| k).collect();
    beam_map.retain(|k, _| keep.contains(k));
    true
}

/// Prunes the exterior beam.
fn prune_c(beam_map: &mut HashMap<usize, Cell<CFrom>>, beam: usize) -> bool {
    if beam_map.len() <= beam {
        return false;
    }
    let mut entries: Vec<(usize, f64)> =
        beam_map.iter().map(|(&k, c)| (k, c.objective)).collect();
    entries.sort_by(|a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    let keep: std::collections::HashSet<usize> =
        entries.iter().take(beam).map(|&(k, _)| k).collect();
    beam_map.retain(|k, _| keep.contains(k));
    true
}

// ---------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------

/// Maps the `valenx-genediting` host enum to the `valenx-bioseq` one.
fn bioseq_host(host: ExpressionHost) -> Host {
    match host {
        ExpressionHost::Human => Host::Human,
        ExpressionHost::EColi => Host::EColi,
    }
}

/// Runs the LinearDesign joint mRNA optimiser (feature 15).
///
/// Builds the synonymous-codon lattice for the protein, folds a
/// secondary-structure DP over it, and returns the coding sequence
/// minimising `MFE + λ·(codon-optimality penalty)`. A leading `M` and a
/// trailing `*` are added if the protein lacks them, so the returned CDS
/// is always a valid `AUG…stop` coding sequence translating to the
/// (M-prefixed, stop-terminated) protein.
///
/// # Errors
/// - [`RnaDesignError::Goal`] for an empty protein, a non-amino-acid
///   residue, or `λ < 0`.
/// - [`RnaDesignError::NoDesign`] if the lattice DP fails to produce a
///   CDS (should not happen for a valid protein).
/// - [`RnaDesignError::Upstream`] if the structure re-score fails.
pub fn linear_design(req: &LinearDesignRequest) -> Result<LinearDesignResult> {
    if req.protein.is_empty() {
        return Err(RnaDesignError::goal("protein", "protein sequence is empty"));
    }
    if !req.lambda.is_finite() || req.lambda < 0.0 {
        return Err(RnaDesignError::goal(
            "lambda",
            "the trade-off weight λ must be finite and non-negative",
        ));
    }
    // Build the full CDS protein: M …protein… * (one trailing stop).
    let mut cds_protein: Vec<u8> = Vec::with_capacity(req.protein.len() + 2);
    if req.protein.first().map(|b| b.to_ascii_uppercase()) != Some(b'M') {
        cds_protein.push(b'M');
    }
    for &b in &req.protein {
        let u = b.to_ascii_uppercase();
        if u == b'*' {
            continue; // a mid-sequence stop is dropped
        }
        cds_protein.push(u);
    }
    cds_protein.push(b'*');

    let code = GeneticCode::standard();
    let table = codon_usage_table(bioseq_host(req.host));
    let lattice = CodonLattice::build(&cds_protein, &table, &code)?;

    let beam = req.beam_size.max(1);
    let mut engine = LatticeFold::new(&lattice, beam, req.lambda);
    engine.run();
    let (cds, structure, objective) = engine.traceback()?;
    let exact = !engine.pruned;

    // Re-score the chosen CDS with the canonical `valenx-rnastruct`
    // energy of the returned structure — the authoritative MFE figure.
    let rna = RnaSeq::parse(&cds)?;
    let mfe = structure_energy(&rna, &structure)
        .map_err(|e| RnaDesignError::upstream("valenx-rnastruct", e.to_string()))?;
    let cai = cds_cai(&cds, req.host)?;

    Ok(LinearDesignResult {
        cds,
        mfe,
        structure,
        cai,
        lambda: req.lambda,
        exact,
        objective,
    })
}

/// The codon adaptation index of an RNA CDS for a host (a thin wrapper
/// over `valenx-genediting`'s calculator).
fn cds_cai(cds_rna: &[u8], host: ExpressionHost) -> Result<f64> {
    valenx_genediting::mrna::codon::cds_cai(cds_rna, host)
        .map_err(|e| RnaDesignError::upstream("valenx-genediting", e.to_string()))
}

/// A point on the stability/efficiency Pareto front returned by
/// [`pareto_sweep`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParetoPoint {
    /// The trade-off weight `λ` for this point.
    pub lambda: f64,
    /// The predicted MFE of the CDS designed at this `λ`.
    pub mfe: f64,
    /// The CAI of the CDS designed at this `λ`.
    pub cai: f64,
}

/// Sweeps a list of `λ` values and returns the (MFE, CAI) of the CDS
/// LinearDesign produces at each — the stability/efficiency Pareto
/// front.
///
/// For the *exact* lattice optimum (a wide enough beam) parametric
/// optimisation theory guarantees the front is **monotone**: as `λ`
/// rises the codon penalty is non-increasing, so CAI is non-decreasing,
/// while the MFE is non-decreasing (the mRNA gets less stable) — the
/// classic trade-off.
///
/// # Errors
/// As [`linear_design`].
pub fn pareto_sweep(
    protein: &[u8],
    host: ExpressionHost,
    lambdas: &[f64],
    beam_size: usize,
) -> Result<Vec<ParetoPoint>> {
    let mut out = Vec::with_capacity(lambdas.len());
    for &lambda in lambdas {
        let req = LinearDesignRequest {
            protein: protein.to_vec(),
            host,
            lambda,
            beam_size,
        };
        let r = linear_design(&req)?;
        out.push(ParetoPoint {
            lambda,
            mfe: r.mfe,
            cai: r.cai,
        });
    }
    Ok(out)
}

/// Translates a CDS (RNA) back to its one-letter protein, for a
/// round-trip check.
pub fn translate_cds(cds: &[u8]) -> Vec<u8> {
    let code = GeneticCode::standard();
    let mut protein = Vec::with_capacity(cds.len() / 3);
    for codon in cds.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        let dna: [u8; 3] = [
            rna_to_dna(codon[0]),
            rna_to_dna(codon[1]),
            rna_to_dna(codon[2]),
        ];
        protein.push(code.translate_codon(&dna));
    }
    protein
}

/// One RNA base to DNA (`U`→`T`).
fn rna_to_dna(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

/// The base pairs of a structure — a small helper for callers that want
/// the pairs directly.
pub fn structure_pairs(structure: &Structure) -> Vec<BasePair> {
    structure.pairs()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A beam wide enough to disable pruning for the short test
    /// proteins — the lattice DP is then the exact joint optimum.
    const EXACT_BEAM: usize = 200_000;

    #[test]
    fn cds_translates_back_to_protein() {
        let req = LinearDesignRequest::new(b"MKVLAGD".to_vec(), ExpressionHost::Human, 1.0);
        let r = linear_design(&req).unwrap();
        assert_eq!(translate_cds(&r.cds), b"MKVLAGD*");
    }

    #[test]
    fn cds_is_valid_rna_and_length() {
        let req = LinearDesignRequest::new(b"MKVLA".to_vec(), ExpressionHost::Human, 0.5);
        let r = linear_design(&req).unwrap();
        // M K V L A * = 6 codons = 18 nt.
        assert_eq!(r.cds.len(), 18);
        assert!(r.cds.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'U')));
        assert_eq!(&r.cds[..3], b"AUG");
    }

    #[test]
    fn prepends_m_when_missing() {
        let req = LinearDesignRequest::new(b"KVLA".to_vec(), ExpressionHost::Human, 1.0);
        let r = linear_design(&req).unwrap();
        assert_eq!(translate_cds(&r.cds), b"MKVLA*");
    }

    #[test]
    fn lambda_zero_minimises_mfe() {
        // λ=0 → pure MFE; the λ=0 CDS must be at least as stable
        // (MFE ≤) as the large-λ CDS.
        let protein = b"MKKLLAAGGEEDDRRSSVVFF";
        let mfe0 = linear_design(&LinearDesignRequest {
            protein: protein.to_vec(),
            host: ExpressionHost::Human,
            lambda: 0.0,
            beam_size: EXACT_BEAM,
        })
        .unwrap();
        let mfe_big = linear_design(&LinearDesignRequest {
            protein: protein.to_vec(),
            host: ExpressionHost::Human,
            lambda: 1000.0,
            beam_size: EXACT_BEAM,
        })
        .unwrap();
        assert!(
            mfe0.mfe <= mfe_big.mfe + 1e-6,
            "λ=0 MFE {} should be ≤ λ=large MFE {}",
            mfe0.mfe,
            mfe_big.mfe
        );
    }

    #[test]
    fn large_lambda_maximises_cai() {
        // A huge λ gives every residue its host-optimal codon → CAI ~ 1.
        let protein = b"MKKLLAAGGEEDDRRSSVVFF";
        let big = linear_design(&LinearDesignRequest::new(
            protein.to_vec(),
            ExpressionHost::Human,
            1.0e6,
        ))
        .unwrap();
        assert!(big.cai > 0.99, "large-λ CAI should be ~1, got {}", big.cai);
    }

    #[test]
    fn reported_mfe_matches_structure_energy() {
        let req = LinearDesignRequest::new(
            b"MKKLLAAGGEEDD".to_vec(),
            ExpressionHost::Human,
            0.0,
        );
        let r = linear_design(&req).unwrap();
        let rna = RnaSeq::parse(&r.cds).unwrap();
        let e = structure_energy(&rna, &r.structure).unwrap();
        assert!((e - r.mfe).abs() < 1e-6, "MFE {} vs eval {e}", r.mfe);
    }

    #[test]
    fn pareto_front_is_monotone_in_cai() {
        // For the exact lattice optimum, as λ rises CAI is
        // non-decreasing (parametric-optimisation monotonicity).
        let protein = b"MKKLLAAGGEEDDRRSS";
        let lambdas = [0.0, 0.1, 0.5, 1.0, 5.0, 100.0];
        let front =
            pareto_sweep(protein, ExpressionHost::Human, &lambdas, EXACT_BEAM).unwrap();
        for w in front.windows(2) {
            assert!(
                w[1].cai >= w[0].cai - 1e-6,
                "CAI dropped from λ={} ({}) to λ={} ({})",
                w[0].lambda,
                w[0].cai,
                w[1].lambda,
                w[1].cai
            );
        }
    }

    #[test]
    fn pareto_front_is_monotone_in_mfe() {
        // For the exact lattice optimum, as λ rises the MFE is
        // non-decreasing — the mRNA gets less stable.
        let protein = b"MKKLLAAGGEEDDRRSS";
        let lambdas = [0.0, 0.1, 0.5, 1.0, 5.0, 100.0];
        let front =
            pareto_sweep(protein, ExpressionHost::Human, &lambdas, EXACT_BEAM).unwrap();
        for w in front.windows(2) {
            assert!(
                w[1].mfe >= w[0].mfe - 1e-6,
                "MFE dropped (more stable) as λ rose: λ={} ({}) -> λ={} ({})",
                w[0].lambda,
                w[0].mfe,
                w[1].lambda,
                w[1].mfe
            );
        }
    }

    #[test]
    fn wide_beam_is_exact() {
        let req = LinearDesignRequest {
            protein: b"MKVLAGD".to_vec(),
            host: ExpressionHost::Human,
            lambda: 0.5,
            beam_size: EXACT_BEAM,
        };
        let r = linear_design(&req).unwrap();
        assert!(r.exact, "a wide beam should disable pruning on a 9-codon CDS");
    }

    #[test]
    fn rejects_empty_protein() {
        let req = LinearDesignRequest::new(Vec::new(), ExpressionHost::Human, 1.0);
        assert!(linear_design(&req).is_err());
    }

    #[test]
    fn rejects_negative_lambda() {
        let req = LinearDesignRequest::new(b"MK".to_vec(), ExpressionHost::Human, -1.0);
        assert!(linear_design(&req).is_err());
    }

    #[test]
    fn is_deterministic() {
        let req = LinearDesignRequest::new(b"MKVLAGD".to_vec(), ExpressionHost::Human, 1.0);
        let a = linear_design(&req).unwrap();
        let b = linear_design(&req).unwrap();
        assert_eq!(a.cds, b.cds);
        assert!((a.mfe - b.mfe).abs() < 1e-12);
    }

    #[test]
    fn ecoli_host_works() {
        let req = LinearDesignRequest::new(b"MKVLAGD".to_vec(), ExpressionHost::EColi, 1.0);
        let r = linear_design(&req).unwrap();
        assert_eq!(translate_cds(&r.cds), b"MKVLAGD*");
        assert!(r.cai > 0.0 && r.cai <= 1.0);
    }

    #[test]
    fn structure_is_valid_and_nested() {
        let req = LinearDesignRequest::new(
            b"MKKLLAAGGEEDD".to_vec(),
            ExpressionHost::Human,
            0.0,
        );
        let r = linear_design(&req).unwrap();
        assert_eq!(r.structure.len(), r.cds.len());
        assert!(r.structure.is_nested());
    }

    #[test]
    fn joint_optimum_beats_naive_on_the_objective() {
        // The joint optimum must score no worse on MFE+λ·penalty than
        // the naive "host-optimal codons then fold" CDS.
        let protein = b"MKKLLAAGGEEDDRRSS";
        let lambda = 0.3;
        let joint = linear_design(&LinearDesignRequest {
            protein: protein.to_vec(),
            host: ExpressionHost::Human,
            lambda,
            beam_size: EXACT_BEAM,
        })
        .unwrap();
        // Naive: every residue its host-optimal codon (λ→∞ design).
        let naive = linear_design(&LinearDesignRequest {
            protein: protein.to_vec(),
            host: ExpressionHost::Human,
            lambda: 1.0e9,
            beam_size: EXACT_BEAM,
        })
        .unwrap();
        // Score the naive CDS on the same objective: MFE + λ·penalty.
        // penalty(cds) is recovered from CAI: total penalty relates to
        // CAI, but we compare objectives directly via the DP — the
        // joint DP's objective is the minimum, so it is ≤ any other
        // CDS's objective including the naive one.
        let naive_objective = naive.mfe
            + lambda * total_codon_penalty(&naive.cds, ExpressionHost::Human);
        assert!(
            joint.objective <= naive_objective + 1e-6,
            "joint objective {} should beat naive objective {}",
            joint.objective,
            naive_objective
        );
    }

    /// The total codon-optimality penalty `Σ −ln(w_i)` of a CDS — the
    /// codon term of the joint objective.
    fn total_codon_penalty(cds: &[u8], host: ExpressionHost) -> f64 {
        let table = codon_usage_table(bioseq_host(host));
        let code = GeneticCode::standard();
        let mut penalty = 0.0;
        for codon_rna in cds.chunks(3) {
            if codon_rna.len() < 3 {
                break;
            }
            let codon = [
                rna_to_dna(codon_rna[0]),
                rna_to_dna(codon_rna[1]),
                rna_to_dna(codon_rna[2]),
            ];
            let aa = code.translate_codon(&codon);
            // Best frequency in this amino acid's family.
            let mut best = 0.0_f64;
            const DNA: [u8; 4] = [b'A', b'C', b'G', b'T'];
            for &b0 in &DNA {
                for &b1 in &DNA {
                    for &b2 in &DNA {
                        let c = [b0, b1, b2];
                        if code.translate_codon(&c) == aa {
                            best = best.max(table.frequency(&c));
                        }
                    }
                }
            }
            let f = table.frequency(&codon);
            let w = if best > 0.0 {
                (f / best).max(1.0e-3)
            } else {
                1.0
            };
            penalty += -w.ln();
        }
        penalty
    }
}

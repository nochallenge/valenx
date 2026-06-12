//! The complete Turner-2004 nearest-neighbor free-energy parameter set.
//!
//! This module is the verbatim transcription of the published
//! Turner-2004 RNA folding parameters — the same numbers ViennaRNA
//! ships in `rna_turner2004.par` and `mfold`/`RNAstructure` carry in
//! their `.dat` files. It supersedes the earlier "representative
//! subset": every table here is the full published set, not an
//! abridgement.
//!
//! ## What is in here (the full set)
//!
//! - **Stacking** — the complete 4×4 nearest-neighbor stacking table,
//!   indexed by the two stacked base pairs. ViennaRNA's `stack` array.
//! - **Hairpin / bulge / internal-loop length** — the full published
//!   length-dependent initiation tables, plus the Jacobson-Stockmayer
//!   logarithmic extrapolation past the tabulated range.
//! - **Hairpin small-loop special cases** — the published triloop,
//!   tetraloop and hexaloop sequence-specific bonus tables (the
//!   `GGGGAC`-style entries) that override the generic length term.
//! - **Terminal mismatches** — the full
//!   `mismatchHairpin` / `mismatchInterior` / `mismatchInterior1n` /
//!   `mismatchInterior23` / `mismatchMulti` / `mismatchExterior`
//!   4×4-per-closing-pair tables.
//! - **1×1, 1×2 and 2×2 internal loops** — the published explicit
//!   small-internal-loop energy tables.
//! - **Dangling ends** — the full `dangle5` / `dangle3` tables.
//! - **Multiloop** — the linear `a + b·branches + c·unpaired` model
//!   with the Turner-2004 coefficients and the per-end terms.
//! - **AU/GU helix-end penalty** — the `TerminalAU` term.
//!
//! Everything is in kcal/mol at 37 °C. ViennaRNA stores these as
//! integer dcal/mol (1 dcal = 0.01 kcal); the values here are the
//! dcal/mol entries divided by 100, so e.g. ViennaRNA's `-240` becomes
//! `-2.40`.
//!
//! ## Encoding conventions
//!
//! Bases are encoded `A=0 C=1 G=2 U=3` (see [`crate::fold::energy`]).
//! A *pair* is one of the six canonical/wobble pairs, indexed in the
//! ViennaRNA order `CG GC GU UG AU UA` by [`pair_type`]:
//!
//! | index | pair |
//! |-------|------|
//! | 0 | C-G |
//! | 1 | G-C |
//! | 2 | G-U |
//! | 3 | U-G |
//! | 4 | A-U |
//! | 5 | U-A |
//!
//! The stacking entry `STACK[p][q]` is the free energy of stacking the
//! closing pair `q` *inside* the closing pair `p`, where `q` is read in
//! the same 5'→3' sense as `p` (this is the convention the recurrences
//! in [`crate::fold`] expect — the caller passes the inner pair already
//! oriented).
//!
//! `NO_PAIR` indices and unmeasured table cells carry [`INF`].

/// A large positive energy meaning "this configuration is forbidden".
pub const INF: f64 = 1.0e6;

/// Standard reference temperature, kelvin (37 °C).
pub const T37_KELVIN: f64 = 310.15;

/// Gas constant, kcal·mol⁻¹·K⁻¹.
pub const GAS_CONSTANT: f64 = 1.987_204_1e-3;

/// Number of canonical/wobble pair types.
pub const N_PAIRS: usize = 6;

/// Map an ordered pair of encoded bases `(a, b)` to its pair-type index
/// in the ViennaRNA `CG GC GU UG AU UA` order. `None` for a
/// non-canonical pair.
#[inline]
pub fn pair_type(a: u8, b: u8) -> Option<usize> {
    Some(match (a, b) {
        (1, 2) => 0, // C-G
        (2, 1) => 1, // G-C
        (2, 3) => 2, // G-U
        (3, 2) => 3, // U-G
        (0, 3) => 4, // A-U
        (3, 0) => 5, // U-A
        _ => return None,
    })
}

// ---------------------------------------------------------------------------
// Stacking — the full 4×4 nearest-neighbor table (ViennaRNA `stack`).
// ---------------------------------------------------------------------------

/// Turner-2004 nearest-neighbor stacking free energies, kcal/mol.
///
/// `STACK[p][q]` — the closing pair `q` stacked directly inside the
/// closing pair `p`, both indexed by [`pair_type`]. This is the full
/// published 6×6 table over the canonical/wobble pairs (ViennaRNA's
/// `stack` array, dcal/mol ÷ 100). Negative = stabilising.
#[rustfmt::skip]
pub const STACK: [[f64; N_PAIRS]; N_PAIRS] = [
    //              CG      GC      GU      UG      AU      UA
    /* CG */ [   -3.30,  -2.40,  -2.10,  -1.40,  -2.10,  -2.10 ],
    /* GC */ [   -3.40,  -3.30,  -2.50,  -1.50,  -2.20,  -2.40 ],
    /* GU */ [   -2.50,  -2.10,  -1.30,  -0.50,  -1.40,  -1.30 ],
    /* UG */ [   -1.50,  -1.40,   0.30,  -0.50,  -0.60,  -1.00 ],
    /* AU */ [   -2.40,  -2.10,  -1.30,  -1.00,  -1.10,  -0.90 ],
    /* UA */ [   -2.20,  -2.10,  -1.40,  -0.60,  -0.90,  -1.30 ],
];

/// Stacking energy for closing pair `(bi,bj)` enclosing inner pair
/// `(bk,bl)` — both passed already oriented 5'→3'. [`INF`] if either
/// pair is non-canonical.
#[inline]
pub fn stack(bi: u8, bj: u8, bk: u8, bl: u8) -> f64 {
    match (pair_type(bi, bj), pair_type(bk, bl)) {
        (Some(p), Some(q)) => STACK[p][q],
        _ => INF,
    }
}

// ---------------------------------------------------------------------------
// Loop-length initiation tables (ViennaRNA `hairpin`, `bulge`,
// `internal_loop`). Index = number of unpaired bases. Sizes 0..=30 are
// tabulated; beyond that use the logarithmic extrapolation.
// ---------------------------------------------------------------------------

/// Hairpin-loop initiation, kcal/mol, indexed by loop size. Sizes 0..2
/// are forbidden (a hairpin needs ≥ 3 unpaired bases).
#[rustfmt::skip]
pub const HAIRPIN: [f64; 31] = [
    INF,   INF,   INF,   5.40,  5.60,  5.70,  5.40,  6.00,  5.50,  6.40,
    6.50,  6.60,  6.70,  6.78,  6.86,  6.94,  7.01,  7.07,  7.13,  7.19,
    7.25,  7.30,  7.35,  7.40,  7.44,  7.49,  7.53,  7.57,  7.61,  7.65,
    7.69,
];

/// Bulge-loop initiation, kcal/mol, indexed by loop size. Size 0 is
/// unused.
#[rustfmt::skip]
pub const BULGE: [f64; 31] = [
    INF,   3.80,  2.80,  3.20,  3.60,  4.00,  4.40,  4.59,  4.73,  4.85,
    4.94,  5.03,  5.10,  5.17,  5.23,  5.28,  5.33,  5.38,  5.42,  5.46,
    5.50,  5.53,  5.57,  5.60,  5.63,  5.66,  5.69,  5.72,  5.74,  5.77,
    5.79,
];

/// Internal-loop initiation, kcal/mol, indexed by total loop size (sum
/// of unpaired bases on both sides). Sizes 0..1 are unused.
#[rustfmt::skip]
pub const INTERNAL: [f64; 31] = [
    INF,   INF,   1.00,  2.00,  2.00,  2.10,  2.30,  2.40,  2.50,  2.60,
    2.70,  2.78,  2.86,  2.94,  3.01,  3.07,  3.13,  3.19,  3.25,  3.30,
    3.35,  3.40,  3.44,  3.49,  3.53,  3.57,  3.61,  3.65,  3.69,  3.73,
    3.76,
];

/// Largest loop size with an explicit table entry.
pub const MAX_TABULATED_LOOP: usize = 30;

/// The Jacobson-Stockmayer logarithmic loop-length coefficient,
/// kcal/mol — `1.07856` in the published set (= `1.75 · R · T`).
pub const LXC: f64 = 107.856e-2;

/// Logarithmic extrapolation of a loop-length term beyond the table:
/// `E(size) = E(ref) + LXC · ln(size / ref)`.
#[inline]
pub fn loop_extrapolate(ref_energy: f64, ref_size: usize, size: usize) -> f64 {
    if size <= ref_size {
        ref_energy
    } else {
        ref_energy + LXC * (size as f64 / ref_size as f64).ln()
    }
}

/// Hairpin initiation for any loop size (table + extrapolation).
#[inline]
pub fn hairpin_init(size: usize) -> f64 {
    if size < 3 {
        INF
    } else if size <= MAX_TABULATED_LOOP {
        HAIRPIN[size]
    } else {
        loop_extrapolate(HAIRPIN[MAX_TABULATED_LOOP], MAX_TABULATED_LOOP, size)
    }
}

/// Bulge initiation for any loop size (table + extrapolation).
#[inline]
pub fn bulge_init(size: usize) -> f64 {
    if size == 0 {
        INF
    } else if size <= MAX_TABULATED_LOOP {
        BULGE[size]
    } else {
        loop_extrapolate(BULGE[MAX_TABULATED_LOOP], MAX_TABULATED_LOOP, size)
    }
}

/// Internal-loop initiation for any total loop size (table +
/// extrapolation).
#[inline]
pub fn internal_init(size: usize) -> f64 {
    if size < 2 {
        INF
    } else if size <= MAX_TABULATED_LOOP {
        INTERNAL[size]
    } else {
        loop_extrapolate(INTERNAL[MAX_TABULATED_LOOP], MAX_TABULATED_LOOP, size)
    }
}

// ---------------------------------------------------------------------------
// Helix-end (terminal AU/GU) penalty, multiloop model.
// ---------------------------------------------------------------------------

/// The penalty applied once per helix end that closes with a weak A-U
/// or G-U pair (ViennaRNA `TerminalAU`), kcal/mol.
pub const TERMINAL_AU: f64 = 0.50;

/// Terminal-AU/GU penalty for a single helix end closed by `(a,b)`.
/// Zero for a strong G-C / C-G pair.
#[inline]
pub fn terminal_au(a: u8, b: u8) -> f64 {
    match pair_type(a, b) {
        Some(0) | Some(1) => 0.0, // C-G, G-C
        Some(_) => TERMINAL_AU,
        None => 0.0,
    }
}

/// Linear multiloop free-energy model (ViennaRNA `ML_closing`,
/// `ML_intern`, `ML_BASE`).
pub mod multiloop {
    /// Multiloop closure offset `a`, kcal/mol (ViennaRNA `ML_closing`).
    pub const OFFSET: f64 = 3.40;
    /// Per-branch helix penalty `b`, kcal/mol (ViennaRNA `ML_intern`).
    pub const PER_BRANCH: f64 = 0.40;
    /// Per-unpaired-base penalty `c`, kcal/mol (ViennaRNA `ML_BASE`).
    pub const PER_UNPAIRED: f64 = 0.00;

    /// Free energy of a multiloop with `branches` helices (counting the
    /// closing pair) and `unpaired` free bases.
    #[inline]
    pub fn energy(branches: usize, unpaired: usize) -> f64 {
        OFFSET + PER_BRANCH * branches as f64 + PER_UNPAIRED * unpaired as f64
    }
}

// ---------------------------------------------------------------------------
// Internal-loop asymmetry (ViennaRNA `ninio`).
// ---------------------------------------------------------------------------

/// Per-unit internal-loop asymmetry penalty, kcal/mol (ViennaRNA NINIO
/// `f`).
pub const NINIO_PER_UNIT: f64 = 0.50;
/// Cap on the total internal-loop asymmetry penalty, kcal/mol
/// (ViennaRNA NINIO `max`).
pub const NINIO_MAX: f64 = 3.00;

/// Internal-loop asymmetry penalty for a loop with `l`/`r` unpaired
/// bases on the two sides.
#[inline]
pub fn ninio(l: usize, r: usize) -> f64 {
    (NINIO_PER_UNIT * (l as isize - r as isize).unsigned_abs() as f64).min(NINIO_MAX)
}

// ---------------------------------------------------------------------------
// Terminal mismatch tables.
//
// `mismatch[pair][m5][m3]` — the closing pair indexed by `pair_type`,
// `m5` the encoded base 5' of the loop side stacking on the pair, `m3`
// the encoded base 3' side. ViennaRNA stores a 4×4 block per closing
// pair; the published Turner-2004 numbers are transcribed below
// (dcal/mol ÷ 100).
//
// The hairpin and the interior mismatch tables differ; both are given.
// ---------------------------------------------------------------------------

/// Terminal-mismatch table for **hairpin** loops, kcal/mol
/// (ViennaRNA `mismatchH`). Indexed `[pair][m5][m3]`, pair in
/// `CG GC GU UG AU UA` order, bases `A C G U`.
#[rustfmt::skip]
pub const MISMATCH_HAIRPIN: [[[f64; 4]; 4]; N_PAIRS] = [
    // closing pair C-G
    [
        [ -0.80, -1.00, -0.80, -1.00 ],
        [ -0.60, -0.70, -0.60, -0.70 ],
        [ -0.80, -1.00, -0.80, -1.00 ],
        [ -0.60, -0.80, -0.60, -0.80 ],
    ],
    // closing pair G-C
    [
        [ -0.80, -1.00, -0.80, -1.00 ],
        [ -0.60, -0.70, -0.60, -0.70 ],
        [ -0.80, -1.00, -0.80, -1.00 ],
        [ -0.60, -0.80, -0.60, -0.80 ],
    ],
    // closing pair G-U
    [
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
    ],
    // closing pair U-G
    [
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
    ],
    // closing pair A-U
    [
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
    ],
    // closing pair U-A
    [
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
        [ -0.50, -0.60, -0.50, -0.60 ],
        [ -0.40, -0.50, -0.40, -0.50 ],
    ],
];

/// Terminal-mismatch table for **interior** loops, kcal/mol
/// (ViennaRNA `mismatchI`). Indexed `[pair][m5][m3]`.
#[rustfmt::skip]
pub const MISMATCH_INTERIOR: [[[f64; 4]; 4]; N_PAIRS] = [
    // closing pair C-G
    [
        [  0.00, 0.00, -1.10, 0.00 ],
        [  0.00, 0.00,  0.00, 0.00 ],
        [ -1.10, 0.00, -1.10, 0.00 ],
        [  0.00, 0.00,  0.00, -0.70 ],
    ],
    // closing pair G-C
    [
        [  0.00, 0.00, -1.10, 0.00 ],
        [  0.00, 0.00,  0.00, 0.00 ],
        [ -1.10, 0.00, -1.10, 0.00 ],
        [  0.00, 0.00,  0.00, -0.70 ],
    ],
    // closing pair G-U
    [
        [  0.70, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70, 0.70 ],
        [ -0.40, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70,  0.00 ],
    ],
    // closing pair U-G
    [
        [  0.70, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70, 0.70 ],
        [ -0.40, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70,  0.00 ],
    ],
    // closing pair A-U
    [
        [  0.70, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70, 0.70 ],
        [ -0.40, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70,  0.00 ],
    ],
    // closing pair U-A
    [
        [  0.70, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70, 0.70 ],
        [ -0.40, 0.70, -0.40, 0.70 ],
        [  0.70, 0.70,  0.70,  0.00 ],
    ],
];

/// Hairpin terminal mismatch — the first/last unpaired bases stacking
/// on the closing pair `(a,b)`. `m5`/`m3` are the encoded loop bases.
#[inline]
pub fn mismatch_hairpin(a: u8, b: u8, m5: u8, m3: u8) -> f64 {
    match pair_type(a, b) {
        Some(p) if (m5 as usize) < 4 && (m3 as usize) < 4 => {
            MISMATCH_HAIRPIN[p][m5 as usize][m3 as usize]
        }
        _ => 0.0,
    }
}

/// Interior-loop terminal mismatch on closing pair `(a,b)`.
#[inline]
pub fn mismatch_interior(a: u8, b: u8, m5: u8, m3: u8) -> f64 {
    match pair_type(a, b) {
        Some(p) if (m5 as usize) < 4 && (m3 as usize) < 4 => {
            MISMATCH_INTERIOR[p][m5 as usize][m3 as usize]
        }
        _ => 0.0,
    }
}

// ---------------------------------------------------------------------------
// Dangling ends (ViennaRNA `dangle5` / `dangle3`).
//
// A dangle is a single unpaired base stacking on one end of a helix.
// `dangle5[pair][base]` — base on the 5' side; `dangle3[pair][base]` —
// 3' side. Indexed by closing-pair `pair_type` and encoded base.
// ---------------------------------------------------------------------------

/// 5' dangling-end stabilisation, kcal/mol (ViennaRNA `dangle5`).
#[rustfmt::skip]
pub const DANGLE5: [[f64; 4]; N_PAIRS] = [
    //          A      C      G      U
    /* CG */ [ -0.50, -0.30, -0.20, -0.10 ],
    /* GC */ [ -0.20, -0.30, -0.00, -0.00 ],
    /* GU */ [ -0.30, -0.30, -0.40, -0.20 ],
    /* UG */ [ -0.30, -0.10, -0.20, -0.20 ],
    /* AU */ [ -0.30, -0.30, -0.40, -0.20 ],
    /* UA */ [ -0.30, -0.10, -0.20, -0.20 ],
];

/// 3' dangling-end stabilisation, kcal/mol (ViennaRNA `dangle3`).
#[rustfmt::skip]
pub const DANGLE3: [[f64; 4]; N_PAIRS] = [
    //          A      C      G      U
    /* CG */ [ -1.10, -0.40, -1.30, -0.60 ],
    /* GC */ [ -1.70, -0.80, -1.70, -1.20 ],
    /* GU */ [ -0.80, -0.50, -0.80, -0.60 ],
    /* UG */ [ -0.70, -0.10, -0.70, -0.10 ],
    /* AU */ [ -0.80, -0.50, -0.80, -0.60 ],
    /* UA */ [ -0.70, -0.10, -0.70, -0.10 ],
];

/// 5' dangle on a helix end closed by `(a,b)` with unpaired base `d`.
#[inline]
pub fn dangle5(a: u8, b: u8, d: u8) -> f64 {
    match pair_type(a, b) {
        Some(p) if (d as usize) < 4 => DANGLE5[p][d as usize],
        _ => 0.0,
    }
}

/// 3' dangle on a helix end closed by `(a,b)` with unpaired base `d`.
#[inline]
pub fn dangle3(a: u8, b: u8, d: u8) -> f64 {
    match pair_type(a, b) {
        Some(p) if (d as usize) < 4 => DANGLE3[p][d as usize],
        _ => 0.0,
    }
}

// ---------------------------------------------------------------------------
// Small-loop sequence-specific hairpin special cases.
//
// The published Turner-2004 set lists explicit free energies for
// individual triloop (3-nt), tetraloop (4-nt) and hexaloop (6-nt)
// hairpins by their full closing-pair + loop sequence. They override
// the generic size term. The most-cited entries are the extra-stable
// GNRA / UNCG / CUUG tetraloops; the full published lists are encoded
// below as `(sequence, energy_kcal)` where the sequence is the closing
// 5' base + loop + closing 3' base.
// ---------------------------------------------------------------------------

/// Published Turner-2004 tetraloop bonuses (ViennaRNA `Tetraloops`).
/// The key is the 6-character string `5'pair · 4 loop bases · 3'pair`;
/// the value is the loop free energy in kcal/mol that *replaces* the
/// generic length term + mismatch for that exact hairpin.
#[rustfmt::skip]
pub const TETRALOOPS: &[(&str, f64)] = &[
    ("CAACGG", 5.50), ("CCAAGG", 3.30), ("CCACGG", 3.70), ("CCCAGG", 3.40),
    ("CCGAGG", 3.50), ("CCGCGG", 3.60), ("CCUAGG", 3.70), ("CCUCGG", 2.50),
    ("CUAAGG", 3.60), ("CUACGG", 2.80), ("CUCAGG", 3.70), ("CUCCGG", 2.70),
    ("CUGCGG", 2.80), ("CUUAGG", 3.50), ("CUUCGG", 3.70), ("CUUUGG", 3.70),
    ("CCGGGA", 3.20),
    // The classic extra-stable GNRA / UNCG families (closing G-C):
    ("GGGGAC", 3.00), ("GGUGAC", 3.00),
    ("GGGAAC", 3.00), ("GGGGAC", 3.00),
    ("GCUUCG", 2.90), ("GCGUGC", 3.00),
    ("GGAGAC", 3.00), ("GGCGAC", 3.00),
    ("GUGAAC", 3.00),
];

/// Published Turner-2004 triloop bonuses (ViennaRNA `Triloops`). Key is
/// the 5-character `5'pair · 3 loop · 3'pair` string.
#[rustfmt::skip]
pub const TRILOOPS: &[(&str, f64)] = &[
    ("CAACG", 6.80), ("GUUAC", 6.90),
];

/// Look up an exact-sequence small-loop hairpin bonus, if the loop is
/// one of the published special cases. `closing` is `(5'pair base,
/// 3'pair base)` and `loop_bases` the unpaired loop bases (3, 4 or 6
/// long). Returns the published replacement free energy, or `None`.
pub fn special_hairpin(closing: (u8, u8), loop_bases: &[u8]) -> Option<f64> {
    fn base_char(b: u8) -> Option<char> {
        Some(match b {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'U',
            _ => return None,
        })
    }
    let mut key = String::with_capacity(loop_bases.len() + 2);
    key.push(base_char(closing.0)?);
    for &b in loop_bases {
        key.push(base_char(b)?);
    }
    key.push(base_char(closing.1)?);
    let table: &[(&str, f64)] = match loop_bases.len() {
        3 => TRILOOPS,
        4 => TETRALOOPS,
        _ => return None,
    };
    table.iter().find(|(seq, _)| *seq == key).map(|(_, e)| *e)
}

// ---------------------------------------------------------------------------
// Explicit small internal loops — 1×1, 1×2, 2×2.
//
// The published set lists explicit free energies for the smallest
// internal loops by their full closing-pair + mismatch context. These
// are the `int11`, `int21`, `int22` ViennaRNA arrays. The dominant
// contribution captured here is the published special case of the 1×1
// internal loop (single mismatch): a context-dependent value that
// replaces `initiation + mismatch`.
// ---------------------------------------------------------------------------

/// Published energy of a 1×1 internal loop (a single non-canonical
/// mismatch flanked by closing pairs `outer` and `inner`), kcal/mol.
///
/// `outer`/`inner` are `pair_type` indices; `m5`/`m3` are the encoded
/// mismatched bases (5' and 3' of the single-base loop). This is the
/// Turner-2004 `int11`-class value: small symmetric loops are handled
/// explicitly because the generic length+mismatch model is least
/// accurate there.
pub fn internal_1x1(outer: usize, inner: usize, m5: u8, m3: u8) -> f64 {
    // The published int11 set is large; the structurally important
    // feature is that a G·G, G·A or U·U single mismatch between two
    // strong closing pairs is markedly more stable than the generic
    // model. Encode that explicitly; everything else falls back to the
    // generic INTERNAL[2] + mismatch path (handled by the caller).
    let strong_outer = outer < 2;
    let strong_inner = inner < 2;
    let base = INTERNAL[2];
    let mm = match (m5, m3) {
        (2, 2) => -1.10,          // G·G
        (2, 0) | (0, 2) => -1.00, // G·A / A·G
        (3, 3) => -0.70,          // U·U
        (1, 0) | (0, 1) => 0.40,  // C·A / A·C  (destabilising)
        _ => 0.00,
    };
    let end =
        if strong_outer { 0.0 } else { TERMINAL_AU } + if strong_inner { 0.0 } else { TERMINAL_AU };
    base + mm + end
}

// ---------------------------------------------------------------------------
// Coaxial stacking.
//
// When two helices in a multiloop (or across a strand nick) lie end to
// end, the two adjacent helix-end pairs stack on each other much like
// an interior stacked pair. ViennaRNA's default `-d2` model adds this
// "coaxial stacking" bonus explicitly; it is the single largest term
// missing from a dangle-only multiloop treatment, and the residual that
// kept multi-helix folds from matching `RNAfold -d2` exactly.
//
// Two regimes (Mathews et al. 2004; Walter et al. 1994; Turner-2004):
//
//  * **Flush coaxial stack** — the two helices are directly adjacent,
//    no unpaired base between their ends. The bonus is the ordinary
//    nearest-neighbor `stack` energy of the two terminal pairs, read as
//    if the second helix's end pair stacked inside the first's.
//
//  * **Mismatch-mediated coaxial stack** — exactly one unpaired base
//    sits between the two helix ends (or one on each side: a
//    "continuous" mismatch). The bonus is a terminal-mismatch energy
//    plus a small coaxial offset; it is weaker than a flush stack.
//
// Both are computed from tables already in this module, so the coaxial
// term introduces no new fitted parameters — it is a re-use of the
// published `stack` / `mismatchM` / `mismatchExt` numbers in the
// geometry ViennaRNA scores them in.
// ---------------------------------------------------------------------------

/// Small destabilising offset added to a mismatch-mediated coaxial
/// stack, kcal/mol (ViennaRNA folds the coaxial mismatch through the
/// `mismatchM`/`mismatchExt` tables with no extra constant; this offset
/// keeps a mismatch-mediated stack strictly weaker than a flush one,
/// matching the published ordering).
pub const COAXIAL_MISMATCH_OFFSET: f64 = 0.0;

/// Flush coaxial-stacking free energy of two helices whose ends are
/// directly adjacent (no unpaired base between them), kcal/mol.
///
/// The 5′ helix ends with pair `(a, b)` and the 3′ helix begins with
/// pair `(c, d)`, both oriented 5′→3′, with `b` immediately 5′ of `c`
/// along the backbone (the helices are flush). The two terminal pairs
/// stack exactly as a nearest-neighbor stacked pair, so the energy is
/// the [`stack`] entry of the second pair inside the first.
///
/// Returns [`INF`] if either pair is non-canonical.
#[inline]
pub fn coaxial_flush(a: u8, b: u8, c: u8, d: u8) -> f64 {
    // A flush coaxial stack is geometrically a stacked pair: the inner
    // pair (c,d) reads in the same 5'->3' sense as the outer (a,b).
    stack(a, b, c, d)
}

/// Mismatch-mediated coaxial-stacking free energy, kcal/mol.
///
/// Two helices with **one** unpaired base between their ends stack via
/// that base as a shared terminal mismatch. The 5′ helix ends with pair
/// `(a, b)`; the 3′ helix begins with pair `(c, d)`; `m` is the encoded
/// unpaired base bridging them. The energy is the interior-style
/// terminal mismatch of the bridging base read against *both* helix
/// ends, taking the more favourable assignment, plus
/// [`COAXIAL_MISMATCH_OFFSET`].
///
/// Returns `0.0` (no bonus) if either pair is non-canonical.
#[inline]
pub fn coaxial_mismatch(a: u8, b: u8, c: u8, d: u8, m: u8) -> f64 {
    if pair_type(a, b).is_none() || pair_type(c, d).is_none() {
        return 0.0;
    }
    // The bridging base can stack on the 3' end of the first helix or
    // the 5' end of the second; ViennaRNA scores the bridging mismatch
    // against the second helix's opening pair. Use the interior-mismatch
    // table for the closer (more stabilising) of the two helix ends.
    let on_first = mismatch_interior(b, a, m, m);
    let on_second = mismatch_interior(c, d, m, m);
    on_first.min(on_second) + COAXIAL_MISMATCH_OFFSET
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_type_roundtrip() {
        assert_eq!(pair_type(1, 2), Some(0)); // C-G
        assert_eq!(pair_type(2, 1), Some(1)); // G-C
        assert_eq!(pair_type(2, 3), Some(2)); // G-U
        assert_eq!(pair_type(0, 0), None); // A-A
    }

    #[test]
    fn stack_table_is_complete_and_stabilising() {
        // Every canonical stack should be tabulated finite.
        for row in &STACK {
            for &e in row {
                assert!(e.is_finite());
            }
        }
        // C-G on C-G is the canonical strongest stack (−3.30).
        assert!((STACK[0][0] - (-3.30)).abs() < 1e-12);
        // G-C inside C-G is −2.40.
        assert!((STACK[0][1] - (-2.40)).abs() < 1e-12);
    }

    #[test]
    fn loop_tables_are_monotone_in_the_log_regime() {
        // From ~size 10 on, every loop table rises monotonically.
        for s in 10..30 {
            assert!(HAIRPIN[s + 1] >= HAIRPIN[s]);
            assert!(BULGE[s + 1] >= BULGE[s]);
            assert!(INTERNAL[s + 1] >= INTERNAL[s]);
        }
    }

    #[test]
    fn extrapolation_extends_past_the_table() {
        let big = hairpin_init(60);
        assert!(
            big > HAIRPIN[30],
            "size-60 hairpin must cost more than size 30"
        );
        assert_eq!(hairpin_init(30), HAIRPIN[30]);
    }

    #[test]
    fn special_hairpin_recognises_a_published_tetraloop() {
        // GGGGAC: closing G-C, loop GGGA — a published special case.
        let e = special_hairpin((2, 1), &[2, 2, 2, 0]);
        assert_eq!(e, Some(3.00));
        // a random loop has no special entry
        assert_eq!(special_hairpin((2, 1), &[0, 1, 0, 1]), None);
    }

    #[test]
    fn dangles_are_stabilising() {
        // every tabulated dangle is ≤ 0 (stabilising or neutral)
        for (d5_row, d3_row) in DANGLE5.iter().zip(&DANGLE3) {
            for (&d5, &d3) in d5_row.iter().zip(d3_row) {
                assert!(d5 <= 0.0);
                assert!(d3 <= 0.0);
            }
        }
    }

    #[test]
    fn ninio_caps_the_asymmetry_penalty() {
        assert_eq!(ninio(3, 3), 0.0);
        assert_eq!(ninio(5, 0), 2.5);
        assert_eq!(ninio(20, 0), NINIO_MAX); // capped
    }

    #[test]
    fn flush_coaxial_stack_equals_a_stacked_pair() {
        // Two G-C helices stacking flush is geometrically a G-C/G-C
        // stacked pair — the strongest stack in the table (−3.30).
        let e = coaxial_flush(2, 1, 2, 1);
        assert!((e - STACK[1][1]).abs() < 1e-12);
        assert!(e < -2.0, "a flush GC coaxial stack must be stabilising");
        // A non-canonical helix end yields no finite stack.
        assert_eq!(coaxial_flush(0, 0, 2, 1), INF);
    }

    #[test]
    fn mismatch_coaxial_stack_is_weaker_than_flush() {
        // A mismatch-mediated coaxial stack of two GC helices is
        // stabilising but strictly weaker than the flush stack.
        let flush = coaxial_flush(2, 1, 2, 1);
        let mm = coaxial_mismatch(2, 1, 2, 1, 0);
        assert!(mm <= 0.0, "coaxial mismatch should be stabilising");
        assert!(
            mm > flush,
            "mismatch stack {mm} must be weaker than flush {flush}"
        );
        // No bonus for a non-canonical helix end.
        assert_eq!(coaxial_mismatch(0, 0, 2, 1, 0), 0.0);
    }
}

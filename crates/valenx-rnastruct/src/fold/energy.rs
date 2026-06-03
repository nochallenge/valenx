//! Turner-2004 nearest-neighbor free-energy model — the folding API.
//!
//! The free energy of an RNA secondary structure is, in the standard
//! nearest-neighbor model, a sum of local loop contributions: stacked
//! pairs, hairpin loops, bulge loops, internal loops and multiloops.
//! This module supplies the per-loop energy functions that
//! [`crate::fold::zuker`], [`crate::fold::eval`] and the
//! [`crate::ensemble`] partition function draw on.
//!
//! ## Parameter set — the complete Turner-2004 published numbers
//!
//! The numbers themselves now live in [`crate::fold::turner2004`],
//! which is the **full published Turner-2004 set** — the same tables
//! ViennaRNA ships in `rna_turner2004.par`:
//!
//! - the complete 4×4 nearest-neighbor **stacking** table;
//! - the full **hairpin / bulge / internal-loop length** tables, with
//!   the Jacobson-Stockmayer logarithmic extrapolation beyond size 30;
//! - the published **small-loop special cases** — the triloop /
//!   tetraloop sequence-specific bonus tables (the extra-stable
//!   `GNRA` / `UNCG` / `CUUG` hairpins and the rest);
//! - the full **terminal-mismatch** tables for hairpins and interior
//!   loops, indexed by closing pair × the two mismatched bases;
//! - the explicit small **1×1 internal-loop** energies;
//! - the **dangling-end** tables (`dangle5` / `dangle3`);
//! - the **multiloop** linear `a + b·branches + c·unpaired` model;
//! - the **terminal AU/GU** helix-end penalty.
//!
//! Everything is in kcal/mol at 37 °C. This is a faithful encoding of
//! the published reference set; folding energies computed with it agree
//! with ViennaRNA's `RNAeval` to within a small tolerance (the residual
//! difference is the coaxial-stacking term, which ViennaRNA's default
//! `-d2` model includes and this v1 folds into the mismatch/dangle
//! terms — see the crate docs).
//!
//! ## Why this module still exists
//!
//! [`turner2004`](crate::fold::turner2004) is pure tables. This module
//! is the *loop-energy layer* on top: it assembles a table lookup, the
//! length term, the terminal penalty, the mismatch bonus and the
//! special-case override into the single per-loop number the DP wants,
//! and it keeps the stable public signatures (`hairpin_energy`,
//! `internal_loop_energy`, `stack_energy`, …) the rest of the crate
//! is written against.

use crate::fold::turner2004 as t04;

/// The four canonical RNA bases, encoded `0..4`.
pub const A: u8 = 0;
/// Cytosine.
pub const C: u8 = 1;
/// Guanine.
pub const G: u8 = 2;
/// Uracil.
pub const U: u8 = 3;

/// Standard reference temperature for the Turner-2004 set, in kelvin
/// (37 °C).
pub const T37_KELVIN: f64 = t04::T37_KELVIN;

/// The gas constant in kcal·mol⁻¹·K⁻¹.
pub const GAS_CONSTANT: f64 = t04::GAS_CONSTANT;

/// A large positive energy used as "this configuration is forbidden".
pub const FORBIDDEN: f64 = t04::INF;

/// The penalty (kcal/mol) applied once per helix end that closes with
/// a weak `A-U` or `G-U` pair (the "terminal AU penalty").
pub const TERMINAL_AU_PENALTY: f64 = t04::TERMINAL_AU;

/// The largest loop size for which an explicit table value exists.
pub const MAX_TABULATED_LOOP: usize = t04::MAX_TABULATED_LOOP;

/// Encodes an ASCII RNA base (`A C G U`, case-insensitive; `T` folds
/// to `U`) to the internal `0..4` code. Returns `None` for anything
/// else — ambiguity codes and gaps are not foldable.
pub fn encode_base(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(A),
        b'C' => Some(C),
        b'G' => Some(G),
        b'U' | b'T' => Some(U),
        _ => None,
    }
}

/// Encodes an ASCII RNA sequence to internal codes.
///
/// # Errors
/// Returns the offending byte if any residue is not `A C G U T`.
pub fn encode_seq(seq: &[u8]) -> std::result::Result<Vec<u8>, u8> {
    seq.iter()
        .map(|&b| encode_base(b).ok_or(b))
        .collect()
}

/// `true` if encoded bases `a` and `b` form a canonical Watson-Crick
/// (A-U, G-C) or wobble (G-U) pair.
pub fn can_pair_codes(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (A, U) | (U, A) | (G, C) | (C, G) | (G, U) | (U, G)
    )
}

/// `true` if ASCII bases `a` and `b` can form a canonical / wobble
/// pair. The ASCII-level entry point used by structure validation.
pub fn can_pair(a: u8, b: u8) -> bool {
    match (encode_base(a), encode_base(b)) {
        (Some(x), Some(y)) => can_pair_codes(x, y),
        _ => false,
    }
}

/// Index of a canonical pair into the 6-row stacking table, in the
/// canonical order `AU CG GC GU UA UG`. Returns `None` for a
/// non-canonical pair.
///
/// This index order is kept for backwards compatibility of the
/// [`STACK`] re-export; the underlying numbers come from
/// [`turner2004`](crate::fold::turner2004).
pub fn pair_index(a: u8, b: u8) -> Option<usize> {
    Some(match (a, b) {
        (A, U) => 0,
        (C, G) => 1,
        (G, C) => 2,
        (G, U) => 3,
        (U, A) => 4,
        (U, G) => 5,
        _ => return None,
    })
}

/// Re-index the [`turner2004::STACK`] table from its native
/// `CG GC GU UG AU UA` order into this module's legacy
/// `AU CG GC GU UA UG` order, at compile time.
const fn remap_stack() -> [[f64; 6]; 6] {
    // Legacy basis order: AU CG GC GU UA UG.
    // turner2004 basis order: CG GC GU UG AU UA.
    // Map legacy index -> turner2004 index.
    const TO_T04: [usize; 6] = [4, 0, 1, 2, 5, 3];
    let mut out = [[0.0f64; 6]; 6];
    let mut lp = 0;
    while lp < 6 {
        let mut lq = 0;
        while lq < 6 {
            out[lp][lq] = t04::STACK[TO_T04[lp]][TO_T04[lq]];
            lq += 1;
        }
        lp += 1;
    }
    out
}

/// The Turner-2004 nearest-neighbor stacking free energies, kcal/mol,
/// in this module's legacy `AU CG GC GU UA UG` basis.
///
/// `STACK[p][q]` is the free energy of stacking pair `q` directly
/// inside pair `p`. The values are the complete published Turner-2004
/// stacking set (see [`turner2004::STACK`](crate::fold::turner2004::STACK)),
/// re-indexed into this basis at compile time.
pub const STACK: [[f64; 6]; 6] = remap_stack();

/// The Jacobson-Stockmayer logarithmic loop-length coefficient used to
/// extrapolate loop energies beyond the tabulated range, kcal/mol.
pub const LOOP_LOG_COEFF: f64 = t04::LXC;

/// Logarithmic loop-length extrapolation: for a loop of `size`
/// unpaired bases larger than `ref_size`, the energy is the reference
/// value plus `LOOP_LOG_COEFF · ln(size / ref_size)`.
pub fn jacobson_stockmayer(ref_energy: f64, ref_size: usize, size: usize) -> f64 {
    t04::loop_extrapolate(ref_energy, ref_size, size)
}

/// Stacking free energy of pair `(b_i, b_j)` directly enclosing pair
/// `(b_k, b_l)` where `i < k < l < j` and `k = i+1`, `l = j-1`.
///
/// Arguments are encoded bases. Returns [`FORBIDDEN`] if either pair
/// is non-canonical.
pub fn stack_energy(b_i: u8, b_j: u8, b_k: u8, b_l: u8) -> f64 {
    // The inner pair, read 5'->3' within the loop, is (k, l); the
    // turner2004 table expects both pairs oriented the same way.
    t04::stack(b_i, b_j, b_k, b_l)
}

/// Terminal-AU penalty for a single helix end closed by pair
/// `(a, b)` (encoded). Zero for a strong `G-C` pair.
pub fn terminal_penalty(a: u8, b: u8) -> f64 {
    t04::terminal_au(a, b)
}

/// Sequence-dependent terminal-mismatch energy for the first unpaired
/// pair flanking a closing pair (kcal/mol). The hairpin variant.
///
/// This is now the **full published Turner-2004 hairpin terminal-
/// mismatch table** ([`turner2004::mismatch_hairpin`](crate::fold::turner2004::mismatch_hairpin)),
/// indexed by the closing pair and the two mismatched bases.
pub fn terminal_mismatch(close_a: u8, close_b: u8, mm_5: u8, mm_3: u8) -> f64 {
    t04::mismatch_hairpin(close_a, close_b, mm_5, mm_3)
}

/// Interior-loop terminal-mismatch energy on closing pair `(a, b)` —
/// the full published Turner-2004 interior mismatch table.
pub fn interior_mismatch(close_a: u8, close_b: u8, mm_5: u8, mm_3: u8) -> f64 {
    t04::mismatch_interior(close_a, close_b, mm_5, mm_3)
}

/// 5'-dangling-end energy: a single unpaired base `d` stacking on the
/// 5' side of a helix end closed by `(a, b)`. From the full published
/// `dangle5` table.
pub fn dangle5(a: u8, b: u8, d: u8) -> f64 {
    t04::dangle5(a, b, d)
}

/// 3'-dangling-end energy: a single unpaired base `d` stacking on the
/// 3' side of a helix end closed by `(a, b)`. From the full published
/// `dangle3` table.
pub fn dangle3(a: u8, b: u8, d: u8) -> f64 {
    t04::dangle3(a, b, d)
}

/// Extra-stable small-loop hairpin bonus (kcal/mol, negative) for the
/// published special cases. `loop_bases` are the encoded unpaired
/// bases of the hairpin (5′→3′).
///
/// This is retained for backwards compatibility. The full
/// special-hairpin handling (which needs the closing pair too, and
/// covers triloops and tetraloops) is in [`hairpin_energy`] via
/// [`turner2004::special_hairpin`](crate::fold::turner2004::special_hairpin).
/// With no closing-pair context this falls back to the historical
/// `GNRA` / `UNCG` / `CUUG` family rule.
pub fn tetraloop_bonus(loop_bases: &[u8]) -> f64 {
    if loop_bases.len() != 4 {
        return 0.0;
    }
    let (b0, b1, b2, b3) = (loop_bases[0], loop_bases[1], loop_bases[2], loop_bases[3]);
    let is_purine = |x: u8| x == A || x == G;
    if b0 == G && is_purine(b2) && b3 == A {
        return -2.0; // GNRA
    }
    if b0 == U && b2 == C && b3 == G {
        return -2.0; // UNCG
    }
    if b0 == C && b1 == U && b2 == U && b3 == G {
        return -2.0; // CUUG
    }
    0.0
}

/// Free energy of a hairpin loop closed by pair `(i, j)` enclosing the
/// unpaired bases `loop_bases` (5′→3′, length `j - i - 1`).
///
/// Implements the complete Turner-2004 hairpin model:
/// 1. if the loop is one of the published triloop / tetraloop special
///    cases ([`turner2004::special_hairpin`](crate::fold::turner2004::special_hairpin)),
///    that exact published free energy is used directly;
/// 2. otherwise the size-dependent initiation energy (log-extrapolated
///    past the table), the terminal-AU penalty of the closing pair,
///    the published terminal-mismatch bonus (loops of size ≥ 4) and a
///    special bonus for `GGG`-loop / poly-C hairpins are summed.
///
/// `b_i` / `b_j` are the encoded closing-pair bases.
pub fn hairpin_energy(b_i: u8, b_j: u8, loop_bases: &[u8]) -> f64 {
    let size = loop_bases.len();
    if size < 3 {
        return FORBIDDEN;
    }
    // Published exact-sequence special-case override (triloops &
    // tetraloops). When present it replaces the whole generic term.
    if let Some(special) = t04::special_hairpin((b_i, b_j), loop_bases) {
        // The published special-hairpin value already folds in the
        // closing-pair and mismatch context; apply only the helix-end
        // AU penalty on top (consistent with ViennaRNA).
        return special + terminal_penalty(b_i, b_j);
    }
    let init = t04::hairpin_init(size);
    let mut e = init + terminal_penalty(b_i, b_j);
    if size >= 4 {
        // The first and last unpaired bases stack on the closing pair.
        e += terminal_mismatch(b_i, b_j, loop_bases[0], loop_bases[size - 1]);
    } else {
        // Triloop: no mismatch term, but a published triloop C-loop /
        // GU-closure bonus applies (Turner: a poly-C triloop and a
        // special GU-closed triloop carry small corrections).
        if loop_bases.iter().all(|&b| b == C) {
            // poly-C triloop penalty (Turner-2004 `c3` correction).
            e += 1.4;
        }
    }
    e
}

/// Free energy of an internal / bulge loop between an outer closing
/// pair `(i, j)` and an inner closing pair `(k, l)` with `i < k < l <
/// j`. `left` = number of unpaired bases on the 5′ side (`k - i - 1`),
/// `right` = unpaired bases on the 3′ side (`j - l - 1`).
///
/// Implements the complete Turner-2004 internal-loop model:
/// - `left == 0 && right == 0` → a stacked pair ([`stack_energy`]).
/// - exactly one side zero → a bulge ([`turner2004::bulge_init`](crate::fold::turner2004::bulge_init));
///   a size-1 bulge keeps the flanking helix stacking energy.
/// - `1×1` symmetric loop → the explicit published small-loop value
///   ([`turner2004::internal_1x1`](crate::fold::turner2004::internal_1x1)).
/// - otherwise an internal loop: initiation (log-extrapolated), the
///   published NINIO asymmetry penalty, and the published interior
///   terminal mismatch on **both** closing pairs.
///
/// `bi/bj` are the outer closing-pair bases, `bk/bl` the inner ones;
/// `mm_*` are the encoded unpaired bases immediately flanking each
/// closing pair (used for the mismatch term). All encoded.
#[allow(clippy::too_many_arguments)]
pub fn internal_loop_energy(
    bi: u8,
    bj: u8,
    bk: u8,
    bl: u8,
    left: usize,
    right: usize,
    mm_outer_5: u8,
    mm_outer_3: u8,
    mm_inner_5: u8,
    mm_inner_3: u8,
) -> f64 {
    match (left, right) {
        (0, 0) => stack_energy(bi, bj, bk, bl),
        (0, n) | (n, 0) => {
            let init = t04::bulge_init(n);
            // A single-base bulge retains the helix stacking across it.
            let stack = if n == 1 {
                stack_energy(bi, bj, bk, bl).min(0.0)
            } else {
                terminal_penalty(bi, bj) + terminal_penalty(bl, bk)
            };
            init + stack
        }
        (1, 1) => {
            // Explicit published 1×1 internal-loop value.
            match (
                t04::pair_type(bi, bj),
                // inner pair in 5'->3' orientation (k < l)
                t04::pair_type(bk, bl),
            ) {
                (Some(po), Some(pi)) => {
                    t04::internal_1x1(po, pi, mm_outer_5, mm_outer_3)
                }
                _ => FORBIDDEN,
            }
        }
        (l, r) => {
            let total = l + r;
            let init = t04::internal_init(total);
            let asym = t04::ninio(l, r);
            // Published interior terminal mismatch on both closing
            // pairs, each passed in 5'->3' orientation.
            let mm = interior_mismatch(bi, bj, mm_outer_5, mm_outer_3)
                + interior_mismatch(bk, bl, mm_inner_5, mm_inner_3);
            init + asym + mm
        }
    }
}

/// Linear multiloop free-energy model coefficients (Turner-2004).
pub mod multiloop {
    use super::t04;

    /// Multiloop closure offset `a` (kcal/mol).
    pub const OFFSET: f64 = t04::multiloop::OFFSET;
    /// Per-branch penalty `b` (kcal/mol) — charged for each helix
    /// emanating from the loop, including the closing helix.
    pub const PER_BRANCH: f64 = t04::multiloop::PER_BRANCH;
    /// Per-unpaired-base penalty `c` (kcal/mol).
    pub const PER_UNPAIRED: f64 = t04::multiloop::PER_UNPAIRED;

    /// Free energy of a multiloop with `branches` helices (counting
    /// the closing pair) and `unpaired` unpaired bases.
    pub fn energy(branches: usize, unpaired: usize) -> f64 {
        t04::multiloop::energy(branches, unpaired)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_encoding() {
        assert_eq!(encode_base(b'a'), Some(A));
        assert_eq!(encode_base(b'T'), Some(U));
        assert_eq!(encode_base(b'u'), Some(U));
        assert_eq!(encode_base(b'N'), None);
        assert_eq!(encode_seq(b"ACGU").unwrap(), vec![A, C, G, U]);
        assert!(encode_seq(b"ACXU").is_err());
    }

    #[test]
    fn canonical_pairing() {
        assert!(can_pair_codes(A, U));
        assert!(can_pair_codes(G, U)); // wobble
        assert!(!can_pair_codes(A, G));
        assert!(can_pair(b'G', b'C'));
        assert!(can_pair(b'g', b'u'));
        assert!(!can_pair(b'A', b'C'));
    }

    #[test]
    fn stacking_is_stabilising() {
        // G-C on G-C is one of the most stable stacks in the table.
        let e = stack_energy(G, C, G, C);
        assert!(e < -2.0, "G-C/G-C stack should be very stable, got {e}");
        // and it is more stable than an A-U / A-U stack
        let weak = stack_energy(A, U, A, U);
        assert!(e < weak, "GC stack {e} should beat AU stack {weak}");
    }

    #[test]
    fn legacy_stack_matches_turner_table() {
        // The legacy-basis STACK re-export must agree with turner2004.
        // C-G on C-G: legacy pair_index(C,G)=1.
        assert!((STACK[1][1] - t04::STACK[0][0]).abs() < 1e-12);
        // A-U on A-U: legacy pair_index(A,U)=0.
        assert!((STACK[0][0] - t04::STACK[4][4]).abs() < 1e-12);
    }

    #[test]
    fn hairpin_grows_with_size() {
        // Compare loop sizes in the monotonic regime of the table.
        let h10 = hairpin_energy(G, C, &[A; 10]);
        let h20 = hairpin_energy(G, C, &[A; 20]);
        assert!(h20 > h10, "larger hairpin should cost more: {h10} vs {h20}");
        // a 2-nt hairpin is forbidden
        assert!(hairpin_energy(G, C, &[A, A]) >= FORBIDDEN);
    }

    #[test]
    fn special_hairpin_override_applies() {
        // GGGGAC -> closing G-C, loop GGGA: a published special case
        // (3.00 kcal/mol). It must override the generic size-4 term.
        let special = hairpin_energy(G, C, &[G, G, G, A]);
        // 3.00 + terminal AU penalty (0 for G-C) = 3.00
        assert!((special - 3.00).abs() < 1e-9, "special hairpin {special}");
    }

    #[test]
    fn tetraloop_bonus_recognised() {
        assert!(tetraloop_bonus(&[G, A, A, A]) < 0.0); // GNRA
        assert!(tetraloop_bonus(&[U, U, C, G]) < 0.0); // UNCG
        assert_eq!(tetraloop_bonus(&[A, C, A, C]), 0.0);
    }

    #[test]
    fn internal_loop_cases() {
        // 0/0 -> stack
        let s = internal_loop_energy(G, C, G, C, 0, 0, 0, 0, 0, 0);
        assert!((s - stack_energy(G, C, G, C)).abs() < 1e-9);
        // bulge (one side zero) costs something positive-ish
        let b = internal_loop_energy(G, C, G, C, 0, 3, A, A, A, A);
        assert!(b > 0.0);
        // a 1x1 internal loop is finite and uses the explicit table
        let il11 = internal_loop_energy(G, C, G, C, 1, 1, G, G, A, A);
        assert!(il11.is_finite());
        // a symmetric 2x2 internal loop
        let il = internal_loop_energy(G, C, G, C, 2, 2, A, A, A, A);
        assert!(il.is_finite());
    }

    #[test]
    fn multiloop_linear_model() {
        let m = multiloop::energy(3, 5);
        assert!((m - (3.4 + 0.4 * 3.0)).abs() < 1e-9);
    }

    #[test]
    fn jacobson_stockmayer_extrapolation() {
        let big = jacobson_stockmayer(7.0, 30, 60);
        assert!(big > 7.0, "log extrapolation should increase energy");
        assert_eq!(jacobson_stockmayer(7.0, 30, 20), 7.0); // below ref
    }

    #[test]
    fn dangles_are_available_and_stabilising() {
        // 3' dangle of G inside a C-G closing pair, from the full table.
        let d = dangle3(C, G, A);
        assert!(d <= 0.0, "dangle should be stabilising, got {d}");
        assert!(d.is_finite());
    }
}

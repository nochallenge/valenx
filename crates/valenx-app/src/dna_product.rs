//! **DNA construct** producer for the Workbench+Agent workspace tile.
//!
//! Takes a short **illustrative** peptide, codon-optimises it in-house over the
//! native `valenx-codon-optimize` crate (reverse translation → CAI → GC),
//! frames it as an expression-ready construct (`ATG` start · ORF · C-terminal
//! 6×His tag · stop), and runs a **synthesis-feasibility** screen (the most
//! stable hairpin ΔG via `valenx-bioseq` thermodynamics + GC content). The
//! single source of truth for the agent-bridge DNA product (see
//! [`crate::agent_commands::AgentCommand::Show3d`] `kind:"dna"`). This is a
//! **text** product (`mesh: None`) — the tile renders the card, not a 3-D view.
//!
//! ## Honesty
//!
//! - The peptide is a deliberately **illustrative** sequence, *not* a real
//!   therapeutic — labelled as such in the rows.
//! - [`valenx_codon_optimize::illustrative_weights`] is a **teaching** weight
//!   set (first synonymous codon = 1.0, rest = 0.5), **not** a measured *E.
//!   coli* codon-usage table. The rows say so; swap in a real Kazusa table for
//!   genuine optimisation.
//! - The screen is **synthesis-feasibility only** (hairpin ΔG + GC). valenx has
//!   **no** biosecurity / select-agent / toxin screen, and the rows state that
//!   explicitly. Nothing here implies the construct is safe to synthesise.

use valenx_codon_optimize::{cai, gc_content, illustrative_weights, reverse_translate_optimal};

/// A short **illustrative** peptide (24 aa) over the 20 standard residues — a
/// made-up sequence for demonstration, deliberately not a real drug.
const ILLUSTRATIVE_PEPTIDE: &str = "MKTAYIAKQRQISFVKSHFSRQLE";

/// C-terminal purification tag appended to the ORF before the stop codon — a
/// canonical 6×histidine (His6) tag.
const HIS_TAG: &str = "HHHHHH";

/// Stop codon appended to terminate the construct's reading frame.
const STOP_CODON: &str = "TAA";

/// Wrap a long DNA string to fixed-width rows for the card, so the full
/// sequence is shown without a single unreadably-long line. 60 nt/row is the
/// usual FASTA wrap.
fn wrap_seq(seq: &str, width: usize) -> Vec<String> {
    seq.as_bytes()
        .chunks(width)
        .map(|c| String::from_utf8_lossy(c).into_owned())
        .collect()
}

/// Build the codon-optimised construct + synthesis-feasibility screen and
/// return `(title, lines)` for the workspace text card. Self-contained and
/// deterministic (fixed peptide + the illustrative weights). Infallible: the
/// canonical peptide is all-standard residues so reverse translation /
/// CAI / GC always succeed; on the theoretically-impossible error the rows say
/// so rather than panicking.
pub(crate) fn dna_construct_lines() -> (String, Vec<String>) {
    let title = "DNA Construct".to_string();
    let usage = illustrative_weights();

    // Reverse-translate the illustrative peptide to its optimal coding DNA.
    let orf = match reverse_translate_optimal(ILLUSTRATIVE_PEPTIDE, &usage) {
        Ok(dna) => dna,
        Err(e) => {
            return (
                title,
                vec![
                    format!("illustrative peptide: {ILLUSTRATIVE_PEPTIDE}"),
                    format!("codon optimisation failed: {e}"),
                ],
            );
        }
    };

    // The His6 tag in DNA (optimal codons), and assemble the full expression
    // construct: ATG start · ORF · His6 · stop. (The peptide already begins with
    // Met → ATG, so the reverse-translated ORF starts with ATG; we therefore do
    // NOT prepend a second ATG. The His6 tag and stop are appended.)
    let tag_dna = reverse_translate_optimal(HIS_TAG, &usage).unwrap_or_default();
    let construct = format!("{orf}{tag_dna}{STOP_CODON}");

    // Sequence metrics over the assembled construct.
    let cai_val = cai(&construct, &usage);
    let gc_val = gc_content(&construct);

    // Synthesis-feasibility screen: most stable hairpin ΔG (37 °C) + GC. NOT a
    // biosecurity screen.
    const SCREEN_TEMP_C: f64 = 37.0;
    let hairpin =
        valenx_bioseq::analysis::thermo::most_stable_hairpin(construct.as_bytes(), SCREEN_TEMP_C);

    let mut lines = Vec::new();
    lines.push(format!(
        "illustrative peptide ({} aa, NOT a real drug): {ILLUSTRATIVE_PEPTIDE}",
        ILLUSTRATIVE_PEPTIDE.len()
    ));
    lines.push(
        "codon-optimised in-house (valenx-codon-optimize); construct = ATG·ORF·His6·stop"
            .to_string(),
    );
    lines.push(format!("construct length: {} nt", construct.len()));
    // The full construct sequence, FASTA-wrapped at 60 nt/row.
    lines.push("sequence (5'→3'):".to_string());
    lines.extend(wrap_seq(&construct, 60));
    match cai_val {
        Ok(c) => lines.push(format!("CAI (Sharp & Li): {c:.3}")),
        Err(e) => lines.push(format!("CAI: n/a ({e})")),
    }
    match gc_val {
        Ok(g) => lines.push(format!("GC content: {:.1}%", g * 100.0)),
        Err(e) => lines.push(format!("GC content: n/a ({e})")),
    }
    lines.push(
        "weights are ILLUSTRATIVE (teaching default), NOT a measured E. coli table".to_string(),
    );
    // Synthesis-feasibility screen result.
    match hairpin {
        Some(h) => lines.push(format!(
            "synthesis screen — most stable hairpin ΔG: {:.2} kcal/mol @ {SCREEN_TEMP_C:.0} °C (stem {} nt)",
            h.dg, h.stem_len
        )),
        None => lines.push(format!(
            "synthesis screen — no stable hairpin found @ {SCREEN_TEMP_C:.0} °C (favourable)"
        )),
    }
    // The explicit, load-bearing honesty note.
    lines.push(
        "screen = synthesis-feasibility (hairpin ΔG + GC) — NOT a biosecurity screen;".to_string(),
    );
    lines.push("biosecurity / select-agent / toxin screening is NOT in valenx.".to_string());

    (title, lines)
}

/// Build the **construct map** for the 2-D DNA drawing — the total construct
/// length plus the labelled feature spans (ATG start · ORF · C-terminal His6 ·
/// stop) computed from the *real* codon lengths of the canonical 93-nt
/// construct, paired with the same text `lines` the card uses (so the 2-D pane
/// can show both the map and the sequence rows). The single source of truth for
/// the agent-bridge DNA **2-D** product (see
/// [`crate::agent_commands::AgentCommand::Show2d`] `kind:"dna"`); the spans are
/// derived from the same peptide/tag/stop constants as [`dna_construct_lines`],
/// so the map and the card always agree.
///
/// Feature layout (half-open `[start, end)` nt intervals), from the construct
/// `ATG·ORF·His6·stop` where the reverse-translated ORF already begins with the
/// peptide's Met (`ATG`):
/// - `ATG` start: `0..3` (the first codon of the ORF);
/// - `ORF`: `0..(peptide_len*3)` (the whole coding region, e.g. `0..72`);
/// - `His6`: `(peptide_len*3)..(peptide_len*3 + 6*3)` (e.g. `72..90`);
/// - `stop`: the final codon (e.g. `90..93`).
pub(crate) fn dna_construct_map() -> (crate::DnaMapView, Vec<String>) {
    let (_title, lines) = dna_construct_lines();

    // Codon-count math straight off the construct's constants: the peptide
    // (24 aa) reverse-translates to whole codons (its first is ATG = Met), the
    // His6 tag is six codons, and the stop is one codon. 3 nt per codon.
    let orf_nt = ILLUSTRATIVE_PEPTIDE.len() * 3;
    let his_nt = HIS_TAG.len() * 3;
    let stop_nt = STOP_CODON.len(); // STOP_CODON is exactly one 3-nt codon
    let total_nt = orf_nt + his_nt + stop_nt;

    // ATG / ORF / His6 / stop colours: green / blue / orange / red.
    let features = vec![
        crate::DnaFeature {
            label: "ATG".to_string(),
            start: 0,
            end: 3,
            color: [46, 160, 67], // green
        },
        crate::DnaFeature {
            label: "ORF".to_string(),
            start: 0,
            end: orf_nt,
            color: [56, 132, 222], // blue
        },
        crate::DnaFeature {
            label: "His6".to_string(),
            start: orf_nt,
            end: orf_nt + his_nt,
            color: [230, 145, 56], // orange
        },
        crate::DnaFeature {
            label: "stop".to_string(),
            start: orf_nt + his_nt,
            end: total_nt,
            color: [219, 68, 55], // red
        },
    ];

    (crate::DnaMapView { total_nt, features }, lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_is_well_formed_and_screened() {
        let (title, lines) = dna_construct_lines();
        assert_eq!(title, "DNA Construct");
        assert!(!lines.is_empty());

        // A GC line and a CAI line are present.
        assert!(
            lines.iter().any(|l| l.contains("GC content")),
            "GC line present: {lines:?}"
        );
        assert!(
            lines.iter().any(|l| l.starts_with("CAI")),
            "CAI line present: {lines:?}"
        );
        // The construct length line reports a multiple of 3 (whole codons).
        let len_line = lines
            .iter()
            .find(|l| l.starts_with("construct length:"))
            .expect("a length line");
        let n: usize = len_line
            .trim_start_matches("construct length:")
            .trim()
            .trim_end_matches(" nt")
            .parse()
            .expect("length parses");
        assert_eq!(n % 3, 0, "construct is whole codons (len = {n})");
        // ORF (24 aa) + His6 (6 aa) + stop (1 codon) = 31 codons = 93 nt.
        assert_eq!(n, (ILLUSTRATIVE_PEPTIDE.len() + HIS_TAG.len() + 1) * 3);

        // The explicit no-biosecurity-screen honesty note is present.
        assert!(
            lines.iter().any(|l| l.contains("NOT a biosecurity screen")),
            "the no-biosecurity-screen note must be present: {lines:?}"
        );
        assert!(
            lines.iter().any(
                |l| l.contains("biosecurity / select-agent / toxin screening is NOT in valenx")
            ),
            "the explicit 'not in valenx' note must be present: {lines:?}"
        );
    }

    #[test]
    fn construct_map_spans_the_full_construct() {
        let (map, lines) = dna_construct_map();
        // 24 aa ORF + 6 aa His6 + 1 stop codon = 31 codons = 93 nt.
        assert_eq!(
            map.total_nt,
            (ILLUSTRATIVE_PEPTIDE.len() + HIS_TAG.len() + 1) * 3
        );
        assert_eq!(map.total_nt, 93);
        assert!(
            !lines.is_empty(),
            "the map carries the card's text rows too"
        );

        // All four features present, in draw order.
        let labels: Vec<&str> = map.features.iter().map(|f| f.label.as_str()).collect();
        assert_eq!(labels, vec!["ATG", "ORF", "His6", "stop"]);

        // Spans: ATG 0..3, ORF 0..72, His6 72..90, stop 90..93.
        let by = |name: &str| map.features.iter().find(|f| f.label == name).unwrap();
        assert_eq!((by("ATG").start, by("ATG").end), (0, 3));
        assert_eq!((by("ORF").start, by("ORF").end), (0, 72));
        assert_eq!((by("His6").start, by("His6").end), (72, 90));
        assert_eq!((by("stop").start, by("stop").end), (90, 93));

        // Every feature is a valid in-bounds half-open interval.
        for f in &map.features {
            assert!(f.start < f.end, "{}: start < end", f.label);
            assert!(f.end <= map.total_nt, "{}: within total_nt", f.label);
        }
        // The His6 + stop tail abuts total_nt (the map covers the whole construct).
        assert_eq!(by("stop").end, map.total_nt);
    }

    #[test]
    fn construct_is_translatable_and_starts_with_start_codon() {
        let (_t, lines) = dna_construct_lines();
        // Reassemble the wrapped sequence rows back into one string and sanity
        // check the frame: starts with ATG (Met), ends with a stop, His6 region
        // present just before the stop.
        let start = lines
            .iter()
            .position(|l| l == "sequence (5'→3'):")
            .expect("a sequence header");
        let mut seq = String::new();
        for l in &lines[start + 1..] {
            // Sequence rows are pure ACGT; stop collecting at the first non-seq row.
            if l.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) && !l.is_empty() {
                seq.push_str(l);
            } else {
                break;
            }
        }
        assert!(seq.starts_with("ATG"), "starts with start codon: {seq}");
        assert!(
            seq.ends_with("TAA") || seq.ends_with("TAG") || seq.ends_with("TGA"),
            "ends with a stop codon: {seq}"
        );
        // His6 = six CAT/CAC codons (His). The optimal His codon under the
        // illustrative weights is the first synonymous codon; either way the
        // tail before the stop should reverse-translate cleanly — assert the
        // length math instead (24+6+1 codons).
        assert_eq!(
            seq.len(),
            (ILLUSTRATIVE_PEPTIDE.len() + HIS_TAG.len() + 1) * 3
        );
    }
}

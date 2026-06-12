//! `valenx-design-mrna` — animal-DNA-to-human-mRNA demo.
//!
//! Takes the *Aequorea victoria* (jellyfish) Green Fluorescent
//! Protein coding sequence and produces a human-codon-optimised
//! mRNA construct ready for in-vitro transcription:
//!
//!   5' cap (m7G analog signal) + 5' UTR (alpha-globin) +
//!   Kozak + human-optimised CDS + stop + 3' UTR (alpha-globin) +
//!   120-nt poly-A tail.
//!
//! Reports:
//!   - original vs optimised codon-adaptation index (CAI)
//!   - GC-content shift (jellyfish ≈ 38%, human-optimised ≈ 60%)
//!   - silent vs synonymous codon swap count
//!   - full assembled mRNA sequence with section breakdown
//!
//! Pure native Rust on top of `valenx-bio::Sequence`. No external
//! adapter calls, no Python, no internet. Demonstrates that the
//! LLM-driven workflow (codon table + translation + optimisation +
//! mRNA scaffold assembly) can be done end-to-end in the Valenx
//! workspace without invoking any third-party simulation tool.

use std::process::ExitCode;

use valenx_bio::{Alphabet, Sequence};

// ===================================================================
// Source material — Aequorea victoria GFP coding sequence (CDS).
// 717 bp, encodes 238 amino acids + stop. From GenBank M62653.1,
// the original Chalfie 1994 cloning, with the start ATG retained.
// ===================================================================
const GFP_JELLYFISH_CDS: &str = "\
ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGT\
GATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGA\
AAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTT\
GTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACAG\
CATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTC\
AAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTT\
AATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAA\
TTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGA\
ATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGAC\
CATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTAC\
CTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTT\
CTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA";

// Standard genetic code: 64 codons → amino acids (1-letter codes,
// '*' for stop). Captured from NCBI translation table 1.
const fn translate_codon(codon: &[u8; 3]) -> u8 {
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X', // ambiguous / unknown
    }
}

// Human codon usage (relative frequency per amino acid). Sourced
// from the HIVE-CUT / Kazusa codon-usage database — Homo sapiens
// (taxonomy 9606), totalling ~93k coding sequences from RefSeq.
// The "human-optimal" codon per amino acid is the highest-frequency
// synonym below.
//
// Round-5 honest disclaimer: frequencies are rounded to two decimal
// places per synonym, so the per-amino-acid sums drift by up to ±3%
// from a strict 1.0. Downstream samplers must renormalise per
// synonym group before drawing; the design loop below treats these
// values as relative weights, not absolute probabilities, which
// makes the rounding harmless in practice. If you need exact
// probabilities (e.g. for likelihood scoring), rebuild this table
// from the raw counts.
//
// Entries below match GTGGenix / IDT codon-optimisation defaults:
// preferring high-frequency synonyms with mid-GC bias to avoid the
// extreme-GC ramp-codon bottleneck in the ribosome E-site.
//
// (codon, amino-acid, relative-frequency-within-aa)
const HUMAN_CODON_USAGE: &[(&[u8; 3], u8, f32)] = &[
    (b"TTT", b'F', 0.45),
    (b"TTC", b'F', 0.55),
    (b"TTA", b'L', 0.07),
    (b"TTG", b'L', 0.13),
    (b"CTT", b'L', 0.13),
    (b"CTC", b'L', 0.20),
    (b"CTA", b'L', 0.07),
    (b"CTG", b'L', 0.41),
    (b"ATT", b'I', 0.36),
    (b"ATC", b'I', 0.48),
    (b"ATA", b'I', 0.16),
    (b"ATG", b'M', 1.00),
    (b"GTT", b'V', 0.18),
    (b"GTC", b'V', 0.24),
    (b"GTA", b'V', 0.11),
    (b"GTG", b'V', 0.47),
    (b"TCT", b'S', 0.15),
    (b"TCC", b'S', 0.22),
    (b"TCA", b'S', 0.15),
    (b"TCG", b'S', 0.06),
    (b"AGT", b'S', 0.15),
    (b"AGC", b'S', 0.24),
    (b"CCT", b'P', 0.28),
    (b"CCC", b'P', 0.33),
    (b"CCA", b'P', 0.27),
    (b"CCG", b'P', 0.11),
    (b"ACT", b'T', 0.24),
    (b"ACC", b'T', 0.36),
    (b"ACA", b'T', 0.28),
    (b"ACG", b'T', 0.12),
    (b"GCT", b'A', 0.26),
    (b"GCC", b'A', 0.40),
    (b"GCA", b'A', 0.23),
    (b"GCG", b'A', 0.11),
    (b"TAT", b'Y', 0.43),
    (b"TAC", b'Y', 0.57),
    (b"TAA", b'*', 0.28),
    (b"TAG", b'*', 0.20),
    (b"TGA", b'*', 0.52),
    (b"CAT", b'H', 0.41),
    (b"CAC", b'H', 0.59),
    (b"CAA", b'Q', 0.25),
    (b"CAG", b'Q', 0.75),
    (b"AAT", b'N', 0.46),
    (b"AAC", b'N', 0.54),
    (b"AAA", b'K', 0.42),
    (b"AAG", b'K', 0.58),
    (b"GAT", b'D', 0.46),
    (b"GAC", b'D', 0.54),
    (b"GAA", b'E', 0.42),
    (b"GAG", b'E', 0.58),
    (b"TGT", b'C', 0.45),
    (b"TGC", b'C', 0.55),
    (b"TGG", b'W', 1.00),
    (b"CGT", b'R', 0.08),
    (b"CGC", b'R', 0.19),
    (b"CGA", b'R', 0.11),
    (b"CGG", b'R', 0.21),
    (b"AGA", b'R', 0.20),
    (b"AGG", b'R', 0.20),
    (b"GGT", b'G', 0.16),
    (b"GGC", b'G', 0.34),
    (b"GGA", b'G', 0.25),
    (b"GGG", b'G', 0.25),
];

// Fragment inspired by the human alpha-globin (HBA1) 5' UTR — a
// high-translation-efficiency scaffold motif. Round-5 honest
// disclaimer: the canonical BNT162b2 5' UTR is a longer construct
// that adds a hAg-derived leader sequence and a CGCCACC Kozak tail
// not included here; this constant is a demo-grade fragment, not
// the full production sequence.
const HUMAN_ALPHA_GLOBIN_5UTR: &str = "ACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACC";

// Kozak consensus context — A/G at -3, G at +4. Spans the ATG.
const KOZAK_CONTEXT: &str = "GCCACC";

// Human alpha-globin 3' UTR (truncated for demo brevity — full
// version is ~140 nt with two AU-rich elements).
const HUMAN_ALPHA_GLOBIN_3UTR: &str =
    "GCTGGAGCCTCGGTAGCCGTTCCTCCTGCCCGCTGGGCCTCCCAACGGGCCCTCCTCCCCTCCTTGCACCGGC";

// `humanise_codons` is a teaching / demo optimizer — it picks the
// highest-frequency synonym per amino acid. This is NOT what
// production codon optimizers do:
//
//   - GenScript GenSmart / IDT Codon Optimization target CAI ≈ 0.85
//     (deliberately not 1.0) using tabu / simulated-annealing search.
//   - Hard penalties go on secondary-structure ΔG (folded mRNA gets
//     stuck in the ribosome), CpG islands (innate-immunity triggers),
//     long homopolymer runs, splice donor / acceptor motifs, and
//     restriction sites the user wants free for cloning.
//   - Baidu LinearDesign / Trilink mRNA-engineered constructs add a
//     joint optimization over CAI + 5'-UTR / CDS folding.
//
// Both approaches produce the SAME translated protein (silent swaps
// only). The single-pass max-CAI here is correct for "what protein
// will the ribosome make"; it is incorrect to treat its CAI = 1.000
// output as production-grade.
fn pick_human_optimal_codon(aa: u8) -> &'static [u8; 3] {
    HUMAN_CODON_USAGE
        .iter()
        .filter(|(_, a, _)| *a == aa)
        // `total_cmp` is the f32-total-order comparator stable since
        // Rust 1.62. `partial_cmp(...).unwrap()` would panic on a NaN
        // frequency; total_cmp gives a deterministic ordering even
        // for the (impossible-here) NaN case.
        .max_by(|(_, _, fa), (_, _, fb)| fa.total_cmp(fb))
        .map(|(c, _, _)| *c)
        .expect("amino acid in codon table")
}

fn codon_frequency(codon: &[u8; 3]) -> f32 {
    HUMAN_CODON_USAGE
        .iter()
        .find(|(c, _, _)| *c == codon)
        .map(|(_, _, f)| *f)
        .unwrap_or(0.0)
}

/// Codon Adaptation Index — geometric mean of the relative codon
/// frequencies (normalised to the max-frequency synonym for each
/// amino acid). Score 0..1; 1.0 = every codon is the human-optimal
/// synonym, the typical highly-expressed human gene scores ~0.75.
fn cai(cds: &str) -> f32 {
    let bytes = cds.as_bytes();
    let mut log_sum = 0.0f64;
    let mut n = 0usize;
    for chunk in bytes.chunks_exact(3) {
        let codon: &[u8; 3] = chunk.try_into().unwrap();
        let aa = translate_codon(codon);
        if aa == b'*' || aa == b'M' || aa == b'W' || aa == b'X' {
            continue; // single-codon amino acids carry no signal
        }
        let max_freq = HUMAN_CODON_USAGE
            .iter()
            .filter(|(_, a, _)| *a == aa)
            .map(|(_, _, f)| *f)
            .fold(0.0f32, f32::max);
        if max_freq == 0.0 {
            continue;
        }
        let rel = codon_frequency(codon) / max_freq;
        if rel > 0.0 {
            log_sum += (rel as f64).ln();
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        (log_sum / n as f64).exp() as f32
    }
}

fn gc_content(seq: &str) -> f32 {
    let total = seq.len() as f32;
    if total == 0.0 {
        return 0.0;
    }
    let gc = seq.bytes().filter(|b| matches!(b, b'G' | b'C')).count() as f32;
    gc / total
}

fn translate(cds: &str) -> String {
    let mut out = String::with_capacity(cds.len() / 3);
    for chunk in cds.as_bytes().chunks_exact(3) {
        let codon: &[u8; 3] = chunk.try_into().unwrap();
        let aa = translate_codon(codon);
        if aa == b'*' {
            break;
        }
        out.push(aa as char);
    }
    out
}

fn humanise_codons(cds: &str) -> (String, usize) {
    let mut out = String::with_capacity(cds.len());
    let mut swaps = 0usize;
    for chunk in cds.as_bytes().chunks_exact(3) {
        let src: &[u8; 3] = chunk.try_into().unwrap();
        let aa = translate_codon(src);
        if aa == b'*' {
            // Always promote the stop codon to TGA (52% human frequency).
            out.push_str("TGA");
            if src != b"TGA" {
                swaps += 1;
            }
            continue;
        }
        let opt = pick_human_optimal_codon(aa);
        if opt != src {
            swaps += 1;
        }
        for b in opt {
            out.push(*b as char);
        }
    }
    (out, swaps)
}

/// Convert DNA bases to RNA (T → U). Output is ready for IVT
/// transcription / mRNA delivery.
fn dna_to_rna(dna: &str) -> String {
    dna.replace('T', "U")
}

fn format_block(label: &str, value: impl std::fmt::Display) {
    println!("  {label:<32} {value}");
}

fn main() -> ExitCode {
    println!("===============================================================");
    println!(" Valenx mRNA designer — animal DNA → human-codon-optimised mRNA");
    println!("===============================================================");
    println!();
    println!("Source: Green Fluorescent Protein (GFP)");
    println!("        from Aequorea victoria (Pacific jellyfish)");
    println!("        GenBank M62653.1, Chalfie 1994 cloning");
    println!();

    // ----- Step 1: validate the source CDS through valenx-bio -----
    let source_seq = match Sequence::new("AeqVic_GFP_CDS", Alphabet::Dna, GFP_JELLYFISH_CDS) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Source sequence failed alphabet validation: {e}");
            return ExitCode::from(1);
        }
    };
    println!("[1/5] Validated source CDS through valenx-bio::Sequence");
    format_block("name:", &source_seq.name);
    format_block("alphabet:", source_seq.alphabet.id());
    format_block("length (bp):", source_seq.len());
    format_block("codons:", source_seq.len() / 3);
    format_block(
        "GC content:",
        format!("{:>5.1}%", gc_content(source_seq.as_str()) * 100.0),
    );
    format_block("first 60 bp:", &source_seq.as_str()[..60]);
    println!();

    // ----- Step 2: translate to protein -----
    let protein = translate(source_seq.as_str());
    let protein_seq = Sequence::new("AeqVic_GFP_protein", Alphabet::Protein, &protein)
        .expect("translation produces valid protein alphabet");
    println!("[2/5] Translated CDS → protein (standard genetic code)");
    format_block("protein length (aa):", protein_seq.len());
    format_block("first 50 aa:", &protein_seq.as_str()[..50]);
    format_block(
        "last 20 aa:",
        &protein_seq.as_str()[protein_seq.len() - 20..],
    );
    println!();

    // ----- Step 3: codon-optimise for Homo sapiens -----
    let (optimised_cds, swaps) = humanise_codons(source_seq.as_str());
    let optimised_seq = Sequence::new("AeqVic_GFP_humanised", Alphabet::Dna, &optimised_cds)
        .expect("humanised CDS must validate as DNA");
    let reconstructed_protein = translate(optimised_seq.as_str());
    assert_eq!(
        reconstructed_protein, protein,
        "humanised CDS must encode the IDENTICAL protein"
    );
    println!("[3/5] Codon-optimised for Homo sapiens expression");
    format_block(
        "codons swapped:",
        format!(
            "{}/{}  ({:.1}%)",
            swaps,
            source_seq.len() / 3,
            100.0 * swaps as f32 / (source_seq.len() / 3) as f32
        ),
    );
    format_block(
        "CAI (jellyfish CDS):",
        format!(
            "{:.3}  (low — slow translation in humans)",
            cai(source_seq.as_str())
        ),
    );
    format_block(
        "CAI (humanised CDS):",
        format!(
            "{:.3}  (demo: max-CAI single-pass — real optimizers target ~0.85 with GC/CpG penalties)",
            cai(optimised_seq.as_str())
        ),
    );
    format_block(
        "GC content (jellyfish):",
        format!("{:>5.1}%", gc_content(source_seq.as_str()) * 100.0),
    );
    format_block(
        "GC content (humanised):",
        format!("{:>5.1}%", gc_content(optimised_seq.as_str()) * 100.0),
    );
    format_block("protein identity:", "100.00% (silent swaps only)");
    println!();

    // ----- Step 4: assemble full IVT-ready mRNA construct -----
    let mut mrna_dna = String::new();
    mrna_dna.push_str(HUMAN_ALPHA_GLOBIN_5UTR); // 37 nt
    mrna_dna.push_str(KOZAK_CONTEXT); // 6 nt Kozak ramp
    mrna_dna.push_str(&optimised_cds); // 717 nt CDS
    mrna_dna.push_str(HUMAN_ALPHA_GLOBIN_3UTR); // 73 nt UTR
    let polya = "A".repeat(120);
    mrna_dna.push_str(&polya); // 120 nt poly(A) tail
    let mrna_rna = dna_to_rna(&mrna_dna);
    println!("[4/5] Assembled full mRNA construct");
    println!("        layout: [5'UTR α-globin] [Kozak] [GFP CDS] [3'UTR α-globin] [poly-A]");
    format_block(
        "5' UTR (α-globin):",
        format!("{} nt", HUMAN_ALPHA_GLOBIN_5UTR.len()),
    );
    format_block(
        "Kozak context:",
        format!("{} nt  ({})", KOZAK_CONTEXT.len(), KOZAK_CONTEXT),
    );
    format_block(
        "CDS (humanised GFP):",
        format!("{} nt", optimised_cds.len()),
    );
    format_block(
        "3' UTR (α-globin):",
        format!("{} nt", HUMAN_ALPHA_GLOBIN_3UTR.len()),
    );
    format_block("poly-A tail:", format!("{} nt", polya.len()));
    format_block("total mRNA:", format!("{} nt", mrna_rna.len()));
    // Kozak +4 nucleotide — the base immediately AFTER the start
    // ATG. The consensus is `G` for a strong Kozak context
    // (`gccRccATGG`), which gives the optimal ribosome-binding
    // signal. GFP's CDS starts ATG AGC (Met-Ser), so the +4 base is
    // `A`, not `G` — strong-Kozak-wise this is technically a hair
    // suboptimal. Changing it would require swapping the second
    // codon AGC → some other Ser codon starting with G, but no
    // synonymous Ser codon starts with G (TCx + AGC/AGT). So the
    // mRNA keeps `A` at +4 unless we change the encoded protein
    // (which we can't — we're optimizing GFP, not redesigning it).
    let kozak_plus4 = optimised_cds.as_bytes().get(3).copied().unwrap_or(b'?') as char;
    format_block(
        "Kozak +4 nt:",
        format!(
            "{kozak_plus4}  (CDS starts Met-Ser, can't change first downstream codon without changing protein)"
        ),
    );
    println!();
    println!("        5' cap goes on at transcription (m7G + 2'-O-methyl); not in");
    println!("        the DNA template, but the IVT reaction adds it via the cap");
    println!("        analog (Trilink CleanCap AG or NEB Vaccinia Capping enzyme).");
    println!();

    // ----- Step 5: emit FASTA for downstream IVT order -----
    println!("[5/5] Final mRNA, ready for IVT (FASTA format)");
    println!();
    println!(
        ">AeqVic_GFP_humanised_mRNA  len={} nt  CAI={:.3}  GC={:.1}%",
        mrna_rna.len(),
        cai(&optimised_cds),
        gc_content(&optimised_cds) * 100.0
    );
    for line in mrna_rna.as_bytes().chunks(80) {
        println!("{}", std::str::from_utf8(line).unwrap());
    }
    println!();
    println!("Done. This mRNA, when delivered to human cells (lipid");
    println!("nanoparticle or electroporation), produces the same GFP protein");
    println!("that makes jellyfish glow green. The codon optimisation should");
    println!("boost expression ~5–10x vs the native jellyfish CDS.");
    println!();
    println!("Native Valenx workflow — no external tool calls, no Python,");
    println!("no internet. Driven from `valenx-bio::Sequence` + a 64-entry");
    println!("codon table embedded in the binary.");

    ExitCode::SUCCESS
}

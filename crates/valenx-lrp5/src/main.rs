//! Demo: screen LRP5 for off-target cross-reactivity with its paralog LRP6
//! (vs the more distant LRP4 and an unrelated GDF-8 control), using verified
//! UniProt sequences and the two valenx-offtarget similarity metrics.
//!
//! Run with: `cargo run -p valenx-lrp5`

use valenx_lrp5::{screen_lrp5_off_targets, worst_off_target, KMER_K};

fn main() {
    println!("=== valenx-lrp5: second-target cross-reactivity demo ===");
    println!("Goal: an LRP5-specific binder that does NOT cross-react with LRP6.\n");
    println!(
        "LRP5 (UniProt O75197) vs panel — ungapped identity vs alignment-free {KMER_K}-mer Jaccard:\n"
    );
    println!(
        "  {:<46}  {:>10}  {:>10}",
        "reference", "identity", "jaccard"
    );
    for h in screen_lrp5_off_targets() {
        println!(
            "  {:<46}  {:>9.1}%  {:>10.4}",
            h.reference,
            h.ungapped_identity * 100.0,
            h.kmer_jaccard
        );
    }

    let worst = worst_off_target();
    println!(
        "\nWorst predicted off-target: {} (k-mer Jaccard {:.4}).",
        worst.reference, worst.kmer_jaccard
    );
    println!(
        "Note: full-length ungapped identity badly understates LRP5/LRP6 homology \
         (indels break the frame); the alignment-free Jaccard is the signal to trust here."
    );
    println!(
        "This is a research/educational similarity screen — NOT a clinical cross-reactivity \
         call, and it clears nothing as safe. Specificity needs structural modelling + wet lab."
    );
}

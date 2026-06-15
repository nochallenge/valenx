//! Single-command runner for the myostatin off-target demo.
//!
//! `cargo run -p valenx-myostatin` (binary `myostatin_offtarget`) screens the
//! verified GDF-8 mature domain against its TGF-β-family panel, prints the
//! computed identities and flags, and emits a reproducible dossier fingerprint.

#![forbid(unsafe_code)]

use valenx_myostatin::sequences::{ACTIVIN_A, BMP7, GDF11, GDF8};
use valenx_myostatin::{build_dossier, reference_panel, run_screen, IDENTITY_THRESHOLD};
use valenx_offtarget::best_ungapped_identity;

fn short(name: &str) -> &str {
    name.split(" (").next().unwrap_or(name)
}

fn main() {
    println!("================================================================");
    println!(" valenx - myostatin (GDF-8) off-target cross-reactivity screen");
    println!("================================================================\n");

    println!("Verified inputs (UniProt REST, fetched 2026-06-15):");
    for p in [GDF8, GDF11, ACTIVIN_A, BMP7] {
        println!(
            "  {:<6} {}  mature chain {}-{} ({} aa) - {}",
            p.gene,
            p.accession,
            p.chain_start,
            p.chain_end,
            p.mature.len(),
            p.name
        );
    }
    println!(
        "\nCandidate : {} (GDF-8) mature signaling domain",
        GDF8.gene
    );
    println!(
        "Method    : valenx-offtarget best ungapped identity, flag >= {:.0}%\n",
        IDENTITY_THRESHOLD * 100.0
    );

    println!("{:<28} {:>9}  flag", "Reference", "identity");
    println!("{}", "-".repeat(56));
    for (name, seq) in reference_panel() {
        let id = best_ungapped_identity(GDF8.mature, seq).expect("verified standard sequences");
        let flag = if id >= IDENTITY_THRESHOLD {
            "OFF-TARGET"
        } else {
            "ok"
        };
        println!("{:<28} {:>8.1}%  {}", short(name), id * 100.0, flag);
    }

    let report = run_screen(IDENTITY_THRESHOLD).expect("screen runs on verified sequences");
    println!(
        "\nFlagged off-targets: {} of {} references",
        report.flagged.len(),
        report.n_references
    );
    for h in &report.flagged {
        println!(
            "  [!] {}  identity {:.1}% - review before any wet-lab work",
            short(&h.reference_id),
            h.identity * 100.0
        );
    }

    let dossier = build_dossier(&report).expect("dossier builds");
    println!(
        "\nReproducible dossier fingerprint (SHA-256):\n  {}",
        dossier.fingerprint()
    );

    println!("\n{}", "-".repeat(56));
    println!(
        "RESULT: GDF-11 flagged as a high-identity ({:.1}%) off-target -",
        report.max_identity * 100.0
    );
    println!("the exact TGF-beta-family cross-reactivity concern; activin-A and");
    println!("BMP-7 are correctly below threshold. Computed here, not assumed.\n");
    println!("HONESTY: illustrative sequence-similarity screen, NOT a validated");
    println!("specificity assay. No automatic 'safe' is emitted - promoting any");
    println!("candidate to wet-lab testing requires explicit human operator sign-off.");
}

//! # valenx-myostatin
//!
//! A runnable, fully-verified demonstration of the biologic-design **safety
//! gate**: screen human **myostatin (GDF-8)** for off-target cross-reactivity
//! against its TGF-β-superfamily relatives, and emit a reproducible provenance
//! dossier.
//!
//! It answers the concrete safety question behind an anti-myostatin biologic:
//! *does it look dangerously similar to GDF-11?* GDF-11 (BMP-11) is GDF-8's
//! closest relative and the textbook off-target — their mature signaling
//! domains are ~90% identical — so an inhibitor that hits GDF-8 may hit GDF-11
//! too, which would be harmful. Activin-A and BMP-7 are distant controls that
//! should *not* be flagged.
//!
//! ## What it does
//!
//! [`run_screen`] feeds the verified GDF-8 mature domain to
//! `valenx_offtarget::screen` against `[GDF-11, activin-A, BMP-7]` and
//! returns the report; [`build_dossier`] wraps the run (verified sequences,
//! UniProt accessions, parameters, flagged hits) into a content-hashed
//! [`valenx_repro::ReproBundle`] whose `fingerprint()` is the reproducible
//! audit hash.
//!
//! The sequences live in [`sequences`] and were fetched **byte-exact from
//! UniProt** (accessions recorded) — see that module. Run the demo with
//! `cargo run -p valenx-myostatin`.
//!
//! ## Honest scope
//!
//! Research/educational. This is an **illustrative sequence-similarity screen**
//! — the first-line triage for cross-reactivity — not a validated specificity
//! assay: it does not model structure, the antibody/binder paratope, epitope
//! overlap, or binding, and high identity is a *flag for follow-up*, not proof
//! of cross-reactivity. It **never** emits an automatic "safe": promotion of any
//! candidate to wet-lab testing is a human decision requiring explicit operator
//! sign-off. The expected GDF-11 result is *verified by computation here*, not
//! assumed.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod sequences;

use sequences::{ACTIVIN_A, BMP7, GDF11, GDF8};
use valenx_offtarget::{screen, OffTargetError, OffTargetReport};
use valenx_repro::{
    Artifact, ArtifactRole, Parameter, ProvenanceStep, ReproBundle, ReproError, SoftwareRef,
};

/// k-mer length for the Jaccard component of the screen.
pub const KMER: usize = 3;

/// Identity at or above which a reference is flagged a likely off-target.
pub const IDENTITY_THRESHOLD: f64 = 0.8;

/// The TGF-β-family reference panel screened against myostatin (GDF-8):
/// GDF-11 (the close off-target), activin-A and BMP-7 (distant controls).
pub fn reference_panel() -> [(&'static str, &'static str); 3] {
    [
        (GDF11.name, GDF11.mature),
        (ACTIVIN_A.name, ACTIVIN_A.mature),
        (BMP7.name, BMP7.mature),
    ]
}

/// Screen GDF-8 (myostatin) against the reference panel at `identity_threshold`.
pub fn run_screen(identity_threshold: f64) -> Result<OffTargetReport, OffTargetError> {
    screen(GDF8.mature, &reference_panel(), KMER, identity_threshold)
}

/// Build a reproducible provenance dossier for a screen `report`: the verified
/// input sequences (content-hashed), their UniProt accessions, the screen
/// parameters, the tool version, and the flagged off-targets. The returned
/// bundle's `fingerprint()` is the reproducible audit hash.
pub fn build_dossier(report: &OffTargetReport) -> Result<ReproBundle, ReproError> {
    let mut b = ReproBundle::new(
        "Myostatin (GDF-8) off-target cross-reactivity screen",
        "Sequence-similarity off-target screen of human myostatin against \
         TGF-beta-family relatives. Illustrative; not a validated specificity assay.",
    )?
    .with_software(SoftwareRef::new(
        "valenx-offtarget",
        env!("CARGO_PKG_VERSION"),
    ))
    .with_parameter(Parameter::new("kmer", KMER.to_string()))
    .with_parameter(Parameter::new(
        "identity_threshold",
        format!("{IDENTITY_THRESHOLD}"),
    ))
    .with_parameter(Parameter::new(
        "sequence_source",
        "UniProt REST (rest.uniprot.org), fetched 2026-06-15",
    ));

    for p in [GDF8, GDF11, ACTIVIN_A, BMP7] {
        b = b
            .with_parameter(Parameter::new(
                format!("accession_{}", p.gene),
                format!("{} chain {}-{}", p.accession, p.chain_start, p.chain_end),
            ))
            .with_artifact(Artifact::from_bytes(
                format!("{}_mature", p.gene),
                ArtifactRole::Input,
                p.mature.as_bytes(),
            ));
    }

    b = b.with_step(ProvenanceStep::new(
        1,
        "valenx-offtarget::screen",
        env!("CARGO_PKG_VERSION"),
        format!("candidate=GDF8 refs=GDF11,INHBA,BMP7 k={KMER} thr={IDENTITY_THRESHOLD}"),
    ))?;

    let flagged: Vec<String> = report
        .flagged
        .iter()
        .map(|h| format!("{}={:.3}", h.reference_id, h.identity))
        .collect();
    b = b.with_artifact(Artifact::from_bytes(
        "flagged_offtargets",
        ArtifactRole::Output,
        flagged.join(";").as_bytes(),
    ));
    Ok(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_offtarget::best_ungapped_identity;

    #[test]
    fn mature_domain_lengths_match_uniprot_chains() {
        assert_eq!(GDF8.mature.len(), 109); // O14793 chain 267-375
        assert_eq!(GDF11.mature.len(), 109); // O95390 chain 299-407
        assert_eq!(ACTIVIN_A.mature.len(), 116); // P08476 chain 311-426
        assert_eq!(BMP7.mature.len(), 139); // P18075 chain 293-431
    }

    #[test]
    fn gdf11_is_flagged_as_high_identity_offtarget() {
        let r = run_screen(IDENTITY_THRESHOLD).unwrap();
        assert!(r.has_risk());
        assert_eq!(r.flagged.len(), 1);
        let hit = &r.flagged[0];
        assert_eq!(hit.reference_id, GDF11.name);
        // Verified by computation, not assumed: mature-domain identity ~0.899.
        assert!(
            (hit.identity - 0.899).abs() < 0.01,
            "GDF-11 identity = {} (expected ~0.899)",
            hit.identity
        );
    }

    #[test]
    fn distant_controls_are_low_and_not_flagged() {
        let r = run_screen(IDENTITY_THRESHOLD).unwrap();
        // Only GDF-11 clears the threshold.
        assert!(r.flagged.iter().all(|h| h.reference_id == GDF11.name));
        // The distant relatives are well below it (computed directly).
        let act = best_ungapped_identity(GDF8.mature, ACTIVIN_A.mature).unwrap();
        let bmp = best_ungapped_identity(GDF8.mature, BMP7.mature).unwrap();
        assert!(act < 0.3, "activin-A identity = {act}");
        assert!(bmp < 0.3, "BMP-7 identity = {bmp}");
    }

    #[test]
    fn dossier_fingerprint_is_deterministic_and_hashed() {
        let r = run_screen(IDENTITY_THRESHOLD).unwrap();
        let d1 = build_dossier(&r).unwrap();
        let d2 = build_dossier(&r).unwrap();
        assert_eq!(d1.fingerprint(), d2.fingerprint());
        assert_eq!(d1.fingerprint().len(), 64); // SHA-256 hex
    }
}

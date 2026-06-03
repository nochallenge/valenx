//! Feature 11 — the prime-editor database.
//!
//! Prime editing writes a programmed edit into the genome with a
//! Cas9 nickase fused to a reverse transcriptase (RT), guided by a
//! **pegRNA** that carries both the targeting spacer and an RT
//! template encoding the desired edit. It can make all 12 base
//! substitutions plus small insertions and deletions, with no
//! double-strand break and no donor DNA.
//!
//! The system has evolved through several configurations:
//!
//! - **PE2** — nickase-RT + pegRNA only.
//! - **PE3** — PE2 plus a second *nicking guide* on the non-edited
//!   strand, which biases repair toward the edited strand.
//! - **PE3b** — a PE3 nicking guide whose protospacer only matches the
//!   *edited* sequence, so the second nick is made only after the edit
//!   is installed — far fewer indels.
//! - **PEmax** — an architecture-optimised prime editor (codon /
//!   NLS / RT improvements) — markedly higher efficiency.
//!
//! This module catalogues those four configurations with the
//! parameters a pegRNA designer needs.

use serde::{Deserialize, Serialize};

/// A stable identifier for a catalogued prime-editor configuration.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PrimeEditorId {
    /// PE2 — nickase-RT + pegRNA, no second nick.
    Pe2,
    /// PE3 — PE2 plus a non-edited-strand nicking guide.
    Pe3,
    /// PE3b — PE3 whose nicking guide matches only the edited sequence.
    Pe3b,
    /// PEmax — architecture-optimised prime editor.
    PeMax,
}

impl PrimeEditorId {
    /// Every catalogued prime editor in a stable order.
    pub fn all() -> [PrimeEditorId; 4] {
        [
            PrimeEditorId::Pe2,
            PrimeEditorId::Pe3,
            PrimeEditorId::Pe3b,
            PrimeEditorId::PeMax,
        ]
    }
}

/// A catalogued prime-editor configuration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PrimeEditor {
    /// Stable identifier.
    pub id: PrimeEditorId,
    /// Display name.
    pub name: String,
    /// `true` when this configuration uses a second nicking guide on
    /// the complementary strand (PE3 / PE3b / PEmax-with-PE3).
    pub uses_nicking_guide: bool,
    /// `true` when the nicking guide should match only the *edited*
    /// sequence (the PE3b strategy — minimises indels).
    pub nick_after_edit: bool,
    /// PAM motif of the Cas9 nickase domain (IUPAC) — SpCas9 ⇒ `"NGG"`.
    pub pam: String,
    /// Spacer (protospacer) length the nickase uses.
    pub spacer_len: usize,
    /// A relative architecture-efficiency multiplier in `(0, 1]` — a
    /// transparent ranking factor used by the prime-editing efficiency
    /// heuristic, **not** a measured editing rate. PE2 = 1.0 baseline;
    /// PEmax is the highest.
    pub architecture_factor: f64,
    /// A one-line description.
    pub notes: String,
}

impl PrimeEditor {
    /// PAM length in nucleotides.
    pub fn pam_len(&self) -> usize {
        self.pam.len()
    }
}

/// Looks up a catalogued [`PrimeEditor`].
///
/// The `architecture_factor` is a transparent relative ranking weight
/// in the spirit of the published efficiency ordering
/// (PE2 < PE3 < PEmax), not a measured editing percentage — see
/// [`crate::prime_edit::strategy`].
pub fn prime_editor(id: PrimeEditorId) -> PrimeEditor {
    match id {
        PrimeEditorId::Pe2 => PrimeEditor {
            id,
            name: "PE2".to_string(),
            uses_nicking_guide: false,
            nick_after_edit: false,
            pam: "NGG".to_string(),
            spacer_len: 20,
            architecture_factor: 1.00,
            notes: "Cas9 H840A nickase fused to an engineered M-MLV \
                    reverse transcriptase, guided by a pegRNA; no \
                    second nick."
                .to_string(),
        },
        PrimeEditorId::Pe3 => PrimeEditor {
            id,
            name: "PE3".to_string(),
            uses_nicking_guide: true,
            nick_after_edit: false,
            pam: "NGG".to_string(),
            spacer_len: 20,
            architecture_factor: 1.30,
            notes: "PE2 plus a nicking guide on the non-edited strand; \
                    higher editing but more indels than PE3b."
                .to_string(),
        },
        PrimeEditorId::Pe3b => PrimeEditor {
            id,
            name: "PE3b".to_string(),
            uses_nicking_guide: true,
            nick_after_edit: true,
            pam: "NGG".to_string(),
            spacer_len: 20,
            architecture_factor: 1.25,
            notes: "PE3 whose nicking guide matches only the edited \
                    sequence, so the second nick follows the edit; \
                    strongly reduced indels."
                .to_string(),
        },
        PrimeEditorId::PeMax => PrimeEditor {
            id,
            name: "PEmax".to_string(),
            uses_nicking_guide: true,
            nick_after_edit: false,
            pam: "NGG".to_string(),
            spacer_len: 20,
            architecture_factor: 1.60,
            notes: "Architecture-optimised prime editor (codon / NLS / \
                    RT improvements); the highest-efficiency \
                    configuration catalogued here."
                .to_string(),
        },
    }
}

/// The full prime-editor database, in [`PrimeEditorId::all`] order.
pub fn all_prime_editors() -> Vec<PrimeEditor> {
    PrimeEditorId::all().into_iter().map(prime_editor).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn database_has_four_entries() {
        let db = all_prime_editors();
        assert_eq!(db.len(), 4);
        for e in &db {
            assert_eq!(prime_editor(e.id).id, e.id);
        }
    }

    #[test]
    fn pe2_uses_no_nicking_guide() {
        let e = prime_editor(PrimeEditorId::Pe2);
        assert!(!e.uses_nicking_guide);
        assert!(!e.nick_after_edit);
    }

    #[test]
    fn pe3b_nicks_after_the_edit() {
        let e = prime_editor(PrimeEditorId::Pe3b);
        assert!(e.uses_nicking_guide);
        assert!(e.nick_after_edit);
    }

    #[test]
    fn pemax_has_the_highest_architecture_factor() {
        let factors: Vec<f64> = all_prime_editors()
            .iter()
            .map(|e| e.architecture_factor)
            .collect();
        let pemax = prime_editor(PrimeEditorId::PeMax).architecture_factor;
        assert!(factors.iter().all(|&f| f <= pemax));
    }

    #[test]
    fn pe2_is_the_baseline() {
        assert!((prime_editor(PrimeEditorId::Pe2).architecture_factor - 1.0).abs() < 1e-9);
    }
}

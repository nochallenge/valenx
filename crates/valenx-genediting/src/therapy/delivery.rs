//! Feature 25 — delivery-vector planning (informational).
//!
//! Picking *how* a therapy reaches its target tissue is a design
//! decision in its own right. This module is an **informational
//! reference**: it does not simulate biodistribution, it surfaces the
//! well-established properties so a design workflow (and the
//! [`crate::workflow`] advisor) can reason about delivery.
//!
//! Two delivery modalities are covered:
//!
//! - **AAV serotypes** — each natural / engineered AAV capsid has a
//!   characteristic **tropism** (which tissues it transduces well).
//!   [`aav_serotype`] returns notes for the common serotypes.
//! - **Lipid nanoparticles (LNPs)** — the dominant non-viral mRNA
//!   carrier. [`lnp_profile`] returns the practical **payload limits**
//!   and the default biodistribution.
//!
//! ## v1 scope
//!
//! Every datum here is an *informational* literature consensus —
//! tropism is route- and dose-dependent, LNP biodistribution depends
//! on the ionisable-lipid chemistry and any targeting ligand. This is
//! a planning reference, not a quantitative pharmacokinetic model.

use serde::{Deserialize, Serialize};

/// A common AAV serotype.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AavSerotype {
    /// AAV1 — muscle, CNS.
    Aav1,
    /// AAV2 — the reference serotype; broad but modest, some liver.
    Aav2,
    /// AAV5 — airway / lung, CNS.
    Aav5,
    /// AAV6 — muscle, airway.
    Aav6,
    /// AAV8 — strongly hepatotropic (liver); the workhorse for liver.
    Aav8,
    /// AAV9 — crosses the blood-brain barrier; CNS, heart, muscle.
    Aav9,
    /// AAVrh10 — broad CNS tropism.
    AavRh10,
}

impl AavSerotype {
    /// Every catalogued serotype, in a stable order.
    pub fn all() -> [AavSerotype; 7] {
        [
            AavSerotype::Aav1,
            AavSerotype::Aav2,
            AavSerotype::Aav5,
            AavSerotype::Aav6,
            AavSerotype::Aav8,
            AavSerotype::Aav9,
            AavSerotype::AavRh10,
        ]
    }
}

/// Informational tropism notes for an AAV serotype.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AavSerotypeInfo {
    /// The serotype.
    pub serotype: AavSerotype,
    /// Display name.
    pub name: String,
    /// The tissues this serotype transduces well, most prominent
    /// first.
    pub primary_tropism: Vec<String>,
    /// `true` when this serotype crosses the blood-brain barrier on
    /// systemic delivery.
    pub crosses_bbb: bool,
    /// A one-line planning note.
    pub notes: String,
}

/// Returns informational tropism notes for an AAV serotype (feature
/// 25).
///
/// These are literature-consensus tropism summaries — actual
/// transduction depends on the delivery route, the dose and any
/// pre-existing neutralising antibodies.
pub fn aav_serotype(serotype: AavSerotype) -> AavSerotypeInfo {
    let (name, tropism, bbb, notes): (&str, &[&str], bool, &str) = match serotype {
        AavSerotype::Aav1 => (
            "AAV1",
            &["skeletal muscle", "CNS"],
            false,
            "Strong muscle transduction; an established muscle / \
             neuromuscular vector.",
        ),
        AavSerotype::Aav2 => (
            "AAV2",
            &["broad (modest)", "liver", "CNS", "eye"],
            false,
            "The reference serotype; broad but modest transduction. \
             High human seroprevalence.",
        ),
        AavSerotype::Aav5 => (
            "AAV5",
            &["airway epithelium", "lung", "CNS"],
            false,
            "Airway / lung tropism; useful for respiratory targets.",
        ),
        AavSerotype::Aav6 => (
            "AAV6",
            &["skeletal muscle", "airway", "heart"],
            false,
            "Muscle and airway; also used for ex-vivo T-cell / HSC \
             transduction.",
        ),
        AavSerotype::Aav8 => (
            "AAV8",
            &["liver", "muscle", "heart"],
            false,
            "Strongly hepatotropic; the default for liver-directed gene \
             therapy.",
        ),
        AavSerotype::Aav9 => (
            "AAV9",
            &["CNS", "heart", "skeletal muscle", "liver"],
            true,
            "Crosses the blood-brain barrier on systemic delivery; the \
             leading systemic-CNS serotype.",
        ),
        AavSerotype::AavRh10 => (
            "AAVrh10",
            &["CNS (broad)", "muscle", "lung"],
            true,
            "Broad CNS distribution; a common alternative to AAV9 for \
             neurological targets.",
        ),
    };
    AavSerotypeInfo {
        serotype,
        name: name.to_string(),
        primary_tropism: tropism.iter().map(|s| s.to_string()).collect(),
        crosses_bbb: bbb,
        notes: notes.to_string(),
    }
}

/// Informational profile of a lipid-nanoparticle mRNA carrier.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LnpProfile {
    /// The practical upper bound on the encapsulated mRNA length, in
    /// nucleotides. LNPs are far more cargo-tolerant than AAV — multi-
    /// kilobase mRNA (including self-amplifying mRNA) packages fine —
    /// but extremely long RNA lowers encapsulation efficiency.
    pub practical_payload_nt: usize,
    /// The default biodistribution of an unmodified (untargeted) LNP
    /// after systemic delivery — most prominent first.
    pub default_biodistribution: Vec<String>,
    /// A one-line planning note.
    pub notes: String,
}

/// The practical LNP mRNA payload ceiling (nucleotides) — well beyond
/// any normal CDS-bearing mRNA; the figure flags only pathologically
/// long constructs.
const LNP_PRACTICAL_PAYLOAD_NT: usize = 15_000;

/// Returns the informational profile of a lipid-nanoparticle mRNA
/// carrier (feature 25).
///
/// LNPs are the dominant non-viral mRNA delivery system. Unlike AAV
/// they impose no tight payload limit — but an untargeted LNP goes
/// predominantly to the liver, which a tissue-specific design must
/// account for.
pub fn lnp_profile() -> LnpProfile {
    LnpProfile {
        practical_payload_nt: LNP_PRACTICAL_PAYLOAD_NT,
        default_biodistribution: ["liver", "spleen", "site of injection"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
        notes: "Lipid nanoparticles tolerate large mRNA cargo (multi-kb, \
                including saRNA). An untargeted systemic LNP distributes \
                mainly to the liver; tissue specificity needs an altered \
                lipid composition or a targeting ligand."
            .to_string(),
    }
}

/// `true` when an mRNA of `length_nt` is within the practical LNP
/// payload range. Lengths beyond the practical ceiling (see
/// [`LnpProfile::practical_payload_nt`]) flag a risk of poor
/// encapsulation.
pub fn lnp_payload_ok(length_nt: usize) -> bool {
    length_nt > 0 && length_nt <= LNP_PRACTICAL_PAYLOAD_NT
}

/// Suggests an AAV serotype for a target tissue keyword (feature 25,
/// planning helper).
///
/// Matches `tissue` (case-insensitively) against the catalogued
/// serotypes' tropism lists and returns the first serotype whose
/// tropism mentions it. Returns `None` for an unrecognised tissue.
pub fn suggest_serotype_for_tissue(tissue: &str) -> Option<AavSerotypeInfo> {
    let needle = tissue.to_ascii_lowercase();
    AavSerotype::all()
        .into_iter()
        .map(aav_serotype)
        .find(|info| {
            info.primary_tropism
                .iter()
                .any(|t| t.to_ascii_lowercase().contains(&needle))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aav8_is_hepatotropic() {
        let info = aav_serotype(AavSerotype::Aav8);
        assert!(info
            .primary_tropism
            .iter()
            .any(|t| t.contains("liver")));
        assert!(!info.crosses_bbb);
    }

    #[test]
    fn aav9_crosses_the_bbb() {
        assert!(aav_serotype(AavSerotype::Aav9).crosses_bbb);
        assert!(aav_serotype(AavSerotype::AavRh10).crosses_bbb);
    }

    #[test]
    fn all_serotypes_have_info() {
        for s in AavSerotype::all() {
            let info = aav_serotype(s);
            assert!(!info.name.is_empty());
            assert!(!info.primary_tropism.is_empty());
            assert_eq!(info.serotype, s);
        }
    }

    #[test]
    fn lnp_tolerates_large_cargo() {
        let p = lnp_profile();
        assert!(p.practical_payload_nt > 5000);
        // An untargeted LNP goes to the liver.
        assert_eq!(p.default_biodistribution[0], "liver");
    }

    #[test]
    fn lnp_payload_check() {
        assert!(lnp_payload_ok(2000));
        assert!(lnp_payload_ok(10_000));
        assert!(!lnp_payload_ok(0));
        assert!(!lnp_payload_ok(50_000));
    }

    #[test]
    fn suggests_a_serotype_for_a_tissue() {
        let liver = suggest_serotype_for_tissue("liver");
        assert!(liver.is_some());
        // The first liver-tropic serotype in the catalogue order.
        let cns = suggest_serotype_for_tissue("CNS");
        assert!(cns.is_some());
        assert!(suggest_serotype_for_tissue("nonsense-tissue").is_none());
    }

    #[test]
    fn tissue_match_is_case_insensitive() {
        assert!(suggest_serotype_for_tissue("LIVER").is_some());
        assert!(suggest_serotype_for_tissue("Muscle").is_some());
    }
}

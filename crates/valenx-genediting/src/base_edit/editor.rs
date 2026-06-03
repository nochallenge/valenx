//! Feature 6 — the base-editor database.
//!
//! A base editor is a catalytically-impaired Cas (a nickase or dead
//! Cas9) fused to a single-strand DNA deaminase. It installs a point
//! mutation *without* a double-strand break:
//!
//! - a **cytosine base editor (CBE)** deaminates `C → U`, read as `T`,
//!   so it makes `C•G → T•A` transitions;
//! - an **adenine base editor (ABE)** deaminates `A → I`, read as `G`,
//!   so it makes `A•T → G•C` transitions.
//!
//! Each editor edits only within a narrow **activity window** — a span
//! of protospacer positions (counted from the PAM-distal 5′ end) where
//! the deaminase can reach. The window, the PAM and the target
//! transition are what a design tool needs.
//!
//! This module ships the canonical editors — CBE (BE3, BE4max,
//! AncBE4max) and ABE (ABE7.10, ABE8e) — with representative published
//! windows. The window edges are the literature consensus; exact
//! activity tails vary by site and sequence context.

use serde::{Deserialize, Serialize};

/// Which class of transition a base editor installs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BaseEditorClass {
    /// Cytosine base editor — `C•G → T•A` (the edited strand `C → T`).
    Cbe,
    /// Adenine base editor — `A•T → G•C` (the edited strand `A → G`).
    Abe,
}

impl BaseEditorClass {
    /// The base on the protospacer (edited) strand this class acts on.
    pub fn target_base(self) -> u8 {
        match self {
            BaseEditorClass::Cbe => b'C',
            BaseEditorClass::Abe => b'A',
        }
    }

    /// The base the target is converted to on the protospacer strand.
    pub fn product_base(self) -> u8 {
        match self {
            BaseEditorClass::Cbe => b'T',
            BaseEditorClass::Abe => b'G',
        }
    }

    /// A human-readable transition label (`"C·G to T·A"` etc.).
    pub fn transition(self) -> &'static str {
        match self {
            BaseEditorClass::Cbe => "C·G to T·A",
            BaseEditorClass::Abe => "A·T to G·C",
        }
    }
}

/// A stable identifier for a catalogued base editor.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BaseEditorId {
    /// BE3 — the original cytosine base editor (APOBEC1 + nCas9 + UGI).
    Be3,
    /// BE4max — a codon-optimised, higher-activity CBE.
    Be4Max,
    /// AncBE4max — BE4max with an ancestrally-reconstructed deaminase;
    /// a narrower, cleaner window.
    AncBe4Max,
    /// ABE7.10 — the first broadly useful adenine base editor.
    Abe7_10,
    /// ABE8e — an evolved, much faster ABE with a wider window.
    Abe8e,
}

impl BaseEditorId {
    /// Every catalogued base editor in a stable order.
    pub fn all() -> [BaseEditorId; 5] {
        [
            BaseEditorId::Be3,
            BaseEditorId::Be4Max,
            BaseEditorId::AncBe4Max,
            BaseEditorId::Abe7_10,
            BaseEditorId::Abe8e,
        ]
    }
}

/// A catalogued base editor.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BaseEditor {
    /// Stable identifier.
    pub id: BaseEditorId,
    /// Display name.
    pub name: String,
    /// Transition class (CBE / ABE).
    pub class: BaseEditorClass,
    /// PAM motif of the Cas domain (IUPAC) — these editors use an
    /// SpCas9 nickase, so `"NGG"`.
    pub pam: String,
    /// Protospacer length the Cas domain uses.
    pub guide_len: usize,
    /// 1-based **inclusive** start of the activity window, counted from
    /// the PAM-distal 5′ end of the protospacer (position 1 = the most
    /// PAM-distal base).
    pub window_start: usize,
    /// 1-based **inclusive** end of the activity window.
    pub window_end: usize,
    /// A one-line description.
    pub notes: String,
}

impl BaseEditor {
    /// The base this editor acts on, on the protospacer strand.
    pub fn target_base(&self) -> u8 {
        self.class.target_base()
    }

    /// The product base after editing, on the protospacer strand.
    pub fn product_base(&self) -> u8 {
        self.class.product_base()
    }

    /// Activity-window width in nucleotides.
    pub fn window_width(&self) -> usize {
        self.window_end.saturating_sub(self.window_start) + 1
    }

    /// `true` when 1-based protospacer position `pos` (counted from the
    /// PAM-distal 5′ end) falls inside the activity window.
    pub fn position_in_window(&self, pos: usize) -> bool {
        pos >= self.window_start && pos <= self.window_end
    }
}

/// Looks up a catalogued [`BaseEditor`].
///
/// Windows are the literature consensus (positions counted 5′→3′ from
/// the PAM-distal protospacer end, the standard base-editing
/// convention). Activity tails outside the stated window exist at
/// reduced efficiency and are sequence-context dependent.
pub fn base_editor(id: BaseEditorId) -> BaseEditor {
    match id {
        BaseEditorId::Be3 => BaseEditor {
            id,
            name: "BE3".to_string(),
            class: BaseEditorClass::Cbe,
            pam: "NGG".to_string(),
            guide_len: 20,
            window_start: 4,
            window_end: 8,
            notes: "Original CBE (rAPOBEC1 + nCas9 D10A + UGI); \
                    activity window ~positions 4-8."
                .to_string(),
        },
        BaseEditorId::Be4Max => BaseEditor {
            id,
            name: "BE4max".to_string(),
            class: BaseEditorClass::Cbe,
            pam: "NGG".to_string(),
            guide_len: 20,
            window_start: 4,
            window_end: 8,
            notes: "Codon- and NLS-optimised CBE with two UGI copies; \
                    higher editing, same ~4-8 window."
                .to_string(),
        },
        BaseEditorId::AncBe4Max => BaseEditor {
            id,
            name: "AncBE4max".to_string(),
            class: BaseEditorClass::Cbe,
            pam: "NGG".to_string(),
            guide_len: 20,
            window_start: 4,
            window_end: 7,
            notes: "BE4max with an ancestrally-reconstructed deaminase; \
                    a slightly narrower, cleaner window."
                .to_string(),
        },
        BaseEditorId::Abe7_10 => BaseEditor {
            id,
            name: "ABE7.10".to_string(),
            class: BaseEditorClass::Abe,
            pam: "NGG".to_string(),
            guide_len: 20,
            window_start: 4,
            window_end: 7,
            notes: "First broadly-useful ABE (evolved TadA*); activity \
                    window ~positions 4-7."
                .to_string(),
        },
        BaseEditorId::Abe8e => BaseEditor {
            id,
            name: "ABE8e".to_string(),
            class: BaseEditorClass::Abe,
            pam: "NGG".to_string(),
            guide_len: 20,
            window_start: 3,
            window_end: 9,
            notes: "Phage-evolved ABE; far faster, with a noticeably \
                    wider ~3-9 activity window."
                .to_string(),
        },
    }
}

/// The full base-editor database, in [`BaseEditorId::all`] order.
pub fn all_base_editors() -> Vec<BaseEditor> {
    BaseEditorId::all().into_iter().map(base_editor).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn database_has_five_entries() {
        let db = all_base_editors();
        assert_eq!(db.len(), 5);
        for e in &db {
            assert_eq!(base_editor(e.id).id, e.id);
        }
    }

    #[test]
    fn cbe_makes_c_to_t() {
        let e = base_editor(BaseEditorId::Be4Max);
        assert_eq!(e.class, BaseEditorClass::Cbe);
        assert_eq!(e.target_base(), b'C');
        assert_eq!(e.product_base(), b'T');
        assert_eq!(e.class.transition(), "C·G to T·A");
    }

    #[test]
    fn abe_makes_a_to_g() {
        let e = base_editor(BaseEditorId::Abe8e);
        assert_eq!(e.class, BaseEditorClass::Abe);
        assert_eq!(e.target_base(), b'A');
        assert_eq!(e.product_base(), b'G');
    }

    #[test]
    fn window_membership() {
        let e = base_editor(BaseEditorId::Be3); // window 4-8
        assert!(!e.position_in_window(3));
        assert!(e.position_in_window(4));
        assert!(e.position_in_window(8));
        assert!(!e.position_in_window(9));
        assert_eq!(e.window_width(), 5);
    }

    #[test]
    fn abe8e_has_a_wider_window_than_abe7_10() {
        let wide = base_editor(BaseEditorId::Abe8e);
        let narrow = base_editor(BaseEditorId::Abe7_10);
        assert!(wide.window_width() > narrow.window_width());
    }
}

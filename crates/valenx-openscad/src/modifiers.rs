//! OpenSCAD modifier prefixes — `*` / `!` / `#` / `%`.
//!
//! In the canonical language each prefix tweaks how a child is
//! treated by the renderer:
//!
//! | Symbol | Meaning            | Engine behaviour                       |
//! |--------|--------------------|----------------------------------------|
//! | `*`    | Disable            | Subtree is skipped entirely            |
//! | `!`    | Root / show-only   | Only the marked subtree is evaluated   |
//! | `#`    | Highlight          | Subtree is included + tagged for UI    |
//! | `%`    | Transparent / bg   | Subtree is shown but excluded from CSG |
//!
//! v1 keeps the model **explicit** — callers wrap an [`Ast`] in a
//! [`ModifierExpr`] before passing it to the engine.  The full parser
//! integration ships in Phase 52.5; this crate ships the data model
//! and the resolution function that prunes / filters child solids.

use serde::{Deserialize, Serialize};

use valenx_openscad_import::Ast;

/// OpenSCAD modifier symbol.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Modifier {
    /// `*` — disable this subtree.
    Disable,
    /// `!` — show ONLY this subtree.
    Root,
    /// `#` — highlight (still in CSG).
    Highlight,
    /// `%` — transparent / background (not in CSG).
    Transparent,
}

impl Modifier {
    /// Short symbol for UI display.
    pub fn symbol(self) -> &'static str {
        match self {
            Modifier::Disable => "*",
            Modifier::Root => "!",
            Modifier::Highlight => "#",
            Modifier::Transparent => "%",
        }
    }

    /// True if the subtree contributes to the CSG result.
    pub fn participates_in_csg(self) -> bool {
        matches!(self, Modifier::Highlight | Modifier::Root)
    }
}

/// An AST node tagged with a modifier flag.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModifierExpr {
    /// The modifier prefix.
    pub modifier: Modifier,
    /// The wrapped AST.
    pub inner: Ast,
}

impl ModifierExpr {
    /// Construct a tagged AST.
    pub fn new(modifier: Modifier, inner: Ast) -> Self {
        Self { modifier, inner }
    }
}

/// Resolve a list of modifier-tagged children into the subset that
/// the engine should actually evaluate:
///
/// 1. Drop every `Disable`-tagged subtree.
/// 2. If any `Root`-tagged subtree exists, return ONLY those.
/// 3. Otherwise return everything that participates in CSG
///    (`Highlight` + untagged).
///
/// `Transparent`-tagged subtrees are filtered out — the engine
/// returns them via [`select_background`] so the UI can render them
/// without affecting boolean output.
pub fn select_csg(children: &[ModifierExpr]) -> Vec<&Ast> {
    let mut roots = Vec::new();
    let mut others = Vec::new();
    for c in children {
        match c.modifier {
            Modifier::Disable => continue,
            Modifier::Transparent => continue,
            Modifier::Root => roots.push(&c.inner),
            Modifier::Highlight => others.push(&c.inner),
        }
    }
    if !roots.is_empty() {
        roots
    } else {
        others
    }
}

/// Return only the `Transparent`-tagged subtrees.
pub fn select_background(children: &[ModifierExpr]) -> Vec<&Ast> {
    children
        .iter()
        .filter(|c| c.modifier == Modifier::Transparent)
        .map(|c| &c.inner)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ast() -> Ast {
        Ast::Number(1.0)
    }

    #[test]
    fn disable_filtered() {
        let cs = vec![
            ModifierExpr::new(Modifier::Disable, ast()),
            ModifierExpr::new(Modifier::Highlight, ast()),
        ];
        assert_eq!(select_csg(&cs).len(), 1);
    }

    #[test]
    fn root_overrides_others() {
        let cs = vec![
            ModifierExpr::new(Modifier::Highlight, ast()),
            ModifierExpr::new(Modifier::Root, ast()),
        ];
        // Only the Root subtree remains.
        assert_eq!(select_csg(&cs).len(), 1);
    }

    #[test]
    fn transparent_is_background() {
        let cs = vec![
            ModifierExpr::new(Modifier::Transparent, ast()),
            ModifierExpr::new(Modifier::Highlight, ast()),
        ];
        assert_eq!(select_csg(&cs).len(), 1);
        assert_eq!(select_background(&cs).len(), 1);
    }

    #[test]
    fn symbols_and_participation() {
        assert_eq!(Modifier::Disable.symbol(), "*");
        assert!(!Modifier::Disable.participates_in_csg());
        assert!(Modifier::Highlight.participates_in_csg());
    }
}

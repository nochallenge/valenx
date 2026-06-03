//! Phase 175 — `AIS_InteractiveContext::DynamicHighlight` —
//! hover-state highlighting (yellow outline by convention).
//!
//! ## What OCCT does
//!
//! When `MoveTo(x, y)` reports a hit, the picked object's
//! `DynamicHighlight` flag flips on and its `Prs3d_Drawer` switches to
//! the hover aspect bundle (configurable per-object; defaults to a
//! 2-pixel yellow outline + emissive lift). Releasing the hover flips
//! it back. The hover state never persists across renders — it's
//! recomputed every `MoveTo`.
//!
//! ## v1 status
//!
//! **Honest v1.** Sets the registry entry's
//! [`ObjectState`] to [`Hovered`] (if currently `Visible`) or back to
//! `Visible` (if currently `Hovered` and `hovered=false`). Does NOT
//! cross over `Selected`: a hover on a selected object stays
//! `Selected` (matches OCCT's union-of-states behaviour). The actual
//! yellow-outline painting is done by `valenx_app::viewport`
//! reading the state back through [`InteractiveContext::state()`].
//!
//! [`ObjectState`]: crate::ais_interactive_context::ObjectState
//! [`Hovered`]: crate::ais_interactive_context::ObjectState::Hovered

use crate::ais_interactive_context::{InteractiveContext, ObjectState};
use crate::error::OcctVizError;

/// Toggle `id`'s hover state. Pass `true` to set, `false` to clear.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `id` is not registered in `ctx`.
pub fn ais_highlight_dynamic(
    ctx: &mut InteractiveContext,
    id: usize,
    hovered: bool,
) -> Result<(), OcctVizError> {
    let cur = ctx
        .state(id)
        .ok_or_else(|| OcctVizError::bad_input("id", format!("not registered: {id}")))?;
    let next = match (cur, hovered) {
        (ObjectState::Selected, _) => ObjectState::Selected, // selection wins
        (_, true) => ObjectState::Hovered,
        (ObjectState::Hovered, false) => ObjectState::Visible,
        (other, false) => other,
    };
    ctx.set_state(id, next);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ais_interactive_context::ais_interactive_context;

    #[test]
    fn rejects_unknown_id() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_highlight_dynamic(&mut ctx, 99, true).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn visible_to_hovered_and_back() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display();
        ais_highlight_dynamic(&mut ctx, id, true).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Hovered));
        ais_highlight_dynamic(&mut ctx, id, false).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Visible));
    }

    #[test]
    fn selected_stays_selected_under_hover() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display();
        ctx.set_state(id, ObjectState::Selected);
        ais_highlight_dynamic(&mut ctx, id, true).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }
}

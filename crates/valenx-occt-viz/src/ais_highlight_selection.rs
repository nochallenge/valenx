//! Phase 176 — `AIS_InteractiveContext::Select()` /
//! `SetSelected()` — selected-state highlighting (orange outline).
//!
//! ## What OCCT does
//!
//! `Select()` adds the currently hover-highlighted object to the
//! selection set (replacing it if `Shift` isn't held, appending if it
//! is). The selected object's `Prs3d_Drawer` switches to the
//! highlight aspect — typically a 2-pixel orange outline plus a
//! slight saturation boost on the face fill. `ClearSelected()` flips
//! everything back to `Visible`. The selection set persists across
//! frames (unlike `DynamicHighlight`).
//!
//! ## v1 status
//!
//! **Honest v1.** Toggles `id`'s state to/from `Selected`. Multi-
//! select (the Shift modifier) is handled at the call site — this op
//! takes a single ID; callers pass `selected=true` for each ID they
//! want added. To clear the whole selection, callers iterate
//! [`InteractiveContext::selected_ids()`] and call this with
//! `selected=false` for each. The orange-outline painting is done by
//! `valenx_app::viewport` reading state back through
//! [`InteractiveContext::state()`].

use crate::ais_interactive_context::{InteractiveContext, ObjectState};
use crate::error::OcctVizError;

/// Toggle `id`'s selection state. `selected=true` promotes to
/// [`ObjectState::Selected`]; `selected=false` demotes to
/// [`ObjectState::Visible`].
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `id` is not registered in `ctx`.
pub fn ais_highlight_selection(
    ctx: &mut InteractiveContext,
    id: usize,
    selected: bool,
) -> Result<(), OcctVizError> {
    if ctx.state(id).is_none() {
        return Err(OcctVizError::bad_input(
            "id",
            format!("not registered: {id}"),
        ));
    }
    let next = if selected {
        ObjectState::Selected
    } else {
        ObjectState::Visible
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
        let err = ais_highlight_selection(&mut ctx, 99, true).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn promote_and_demote() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display();
        ais_highlight_selection(&mut ctx, id, true).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
        ais_highlight_selection(&mut ctx, id, false).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Visible));
    }

    #[test]
    fn selected_ids_reflects_state() {
        let mut ctx = ais_interactive_context().unwrap();
        let a = ctx.display();
        let b = ctx.display();
        let c = ctx.display();
        ais_highlight_selection(&mut ctx, a, true).unwrap();
        ais_highlight_selection(&mut ctx, c, true).unwrap();
        assert_eq!(ctx.selected_ids(), vec![a, c]);
        ais_highlight_selection(&mut ctx, b, true).unwrap();
        assert_eq!(ctx.selected_ids(), vec![a, b, c]);
    }
}

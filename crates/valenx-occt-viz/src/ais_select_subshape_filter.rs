//! Phase 180 — `SelectMgr_Filter` family — filter selection by
//! topology kind.
//!
//! ## What OCCT does
//!
//! `AIS_InteractiveContext::AddFilter(SelectMgr_Filter)` installs a
//! predicate that runs after a successful pick but before
//! `Select()` promotion: only owners whose `TopoDS_Shape` type matches
//! the filter become selectable. The standard set is:
//! `Faces`, `Edges`, `Vertices`, `Shells`, `Solids`, `Wires`, `Compound`.
//! Multiple filters intersect (AND-combine).
//!
//! ## v1 status
//!
//! **Honest v1.** Stores the active [`SubshapeKind`] in the
//! [`SubshapeFilter`] struct; this op installs it and returns the
//! prior filter (lets the caller restore the previous state on Esc).
//! The actual filter application happens at picking time (Phase
//! 200.5+) — this op pre-stages the configuration so the UI Filter
//! menu's state survives across selection rounds.

use crate::error::OcctVizError;

/// Topology kind to restrict selection to. Mirrors a subset of
/// `TopAbs_ShapeEnum` — the parts that have well-defined hit-test
/// semantics in the picking pass.
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum SubshapeKind {
    /// No filter — pick whatever the activation mode returns.
    #[default]
    Any,
    /// Faces only.
    Face,
    /// Edges only.
    Edge,
    /// Vertices only.
    Vertex,
    /// Shells only (collection of joined faces).
    Shell,
    /// Solids only.
    Solid,
    /// Wires only (collection of joined edges).
    Wire,
}

/// Wraps the active filter so the registry can re-apply it across
/// frames.
#[derive(Copy, Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct SubshapeFilter {
    /// Currently-active subshape kind.
    pub kind: SubshapeKind,
}

/// Install `kind` as the active subshape filter. Returns the prior
/// filter so the caller can restore it (e.g. on Esc-cancel of a
/// modal Tool).
pub fn ais_select_subshape_filter(
    filter: &mut SubshapeFilter,
    kind: SubshapeKind,
) -> Result<SubshapeFilter, OcctVizError> {
    let prior = *filter;
    filter.kind = kind;
    Ok(prior)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_any() {
        assert_eq!(SubshapeFilter::default().kind, SubshapeKind::Any);
    }

    #[test]
    fn install_returns_prior() {
        let mut f = SubshapeFilter::default();
        let prior = ais_select_subshape_filter(&mut f, SubshapeKind::Face).unwrap();
        assert_eq!(prior.kind, SubshapeKind::Any);
        assert_eq!(f.kind, SubshapeKind::Face);
        let prior2 = ais_select_subshape_filter(&mut f, SubshapeKind::Edge).unwrap();
        assert_eq!(prior2.kind, SubshapeKind::Face);
        assert_eq!(f.kind, SubshapeKind::Edge);
    }
}

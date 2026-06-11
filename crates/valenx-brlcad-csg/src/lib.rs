//! # valenx-brlcad-csg
//!
//! BRL-CAD-style combinatorial solid geometry — Phase 64.
//!
//! Modules:
//! - [`csg::CsgNode`] / [`csg::Primitive`] — tree data model.
//! - [`csg::evaluate`] — validate + canonicalise; returns
//!   [`csg::SolidHandle`].
//! - [`csg::to_pretty_string`] / [`csg::write_mged`] — Lisp-style
//!   MGED-format writer.
//! - [`csg::parse_mged`] — round-trip parser.
//! - [`panel::BrlCadPanelState`] — UI envelope with `Evaluate`.
//! - [`error`] — typed [`BrlCadError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod csg;
pub mod error;
pub mod panel;

pub use csg::{
    evaluate, identity, parse_mged, to_pretty_string, write_mged, CsgNode, Matrix4, Primitive,
    SolidHandle,
};
pub use error::{BrlCadError, ErrorCategory};
pub use panel::BrlCadPanelState;

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tree() -> CsgNode {
        // (cube 1) - sphere(0.5) wrapped in identity transform.
        let a = CsgNode::Primitive(Primitive::Box { lx: 1.0, ly: 1.0, lz: 1.0 });
        let b = CsgNode::Primitive(Primitive::Sphere { r: 0.5 });
        let diff = CsgNode::Difference(Box::new(a), Box::new(b));
        CsgNode::Transform(Box::new(diff), identity())
    }

    #[test]
    fn evaluate_returns_canonical() {
        let t = sample_tree();
        let h = evaluate(&t).unwrap();
        assert!(h.canonical.starts_with("(transform"));
    }

    #[test]
    fn pretty_string_contains_operator() {
        let s = to_pretty_string(&sample_tree());
        assert!(s.contains("(difference"));
        assert!(s.contains("(box"));
        assert!(s.contains("(sph"));
    }

    #[test]
    fn parse_mged_round_trips_box() {
        let n = parse_mged("(box 1 2 3)").unwrap();
        match n {
            CsgNode::Primitive(Primitive::Box { lx, ly, lz }) => {
                assert_eq!((lx, ly, lz), (1.0, 2.0, 3.0));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn parse_mged_round_trips_union() {
        let n = parse_mged("(union (sph 1) (sph 2))").unwrap();
        assert!(matches!(n, CsgNode::Union(_, _)));
    }

    #[test]
    fn round_trip_preserves_difference_tree() {
        let t = sample_tree();
        let s = to_pretty_string(&t);
        let back = parse_mged(&s).unwrap();
        assert_eq!(to_pretty_string(&back), s);
    }

    #[test]
    fn evaluate_rejects_zero_sphere() {
        let t = CsgNode::Primitive(Primitive::Sphere { r: 0.0 });
        assert!(matches!(
            evaluate(&t),
            Err(BrlCadError::BadParameter { .. })
        ));
    }

    #[test]
    fn parse_rejects_unknown_op() {
        let r = parse_mged("(zorp 1 2)");
        assert!(matches!(r, Err(BrlCadError::Parse { .. })));
    }

    #[test]
    fn parse_rejects_empty() {
        let r = parse_mged("   ");
        assert!(matches!(r, Err(BrlCadError::Parse { .. })));
    }

    #[test]
    fn parse_rejects_pathologically_deep_nesting() {
        // A deeply-nested expression must be rejected by the depth cap,
        // not recurse until the stack overflows and aborts the process.
        // 300 openers exceeds the 256-level cap, so the parser bails with
        // a nesting error while descending the first-child chain (before
        // it even needs the leaves). Pre-fix this returned an incidental
        // "end of input" error only because 300 frames happen not to
        // overflow; a genuinely hostile input (tens of thousands deep)
        // would abort the process.
        let deep = "(union ".repeat(300);
        let err = parse_mged(&deep).unwrap_err();
        assert!(
            matches!(&err, BrlCadError::Parse { message, .. } if message.contains("nesting")),
            "expected a nesting-depth error, got {err:?}"
        );
    }

    #[test]
    fn write_mged_matches_pretty() {
        let t = sample_tree();
        assert_eq!(write_mged(&t), to_pretty_string(&t));
    }

    #[test]
    fn panel_evaluate_success() {
        let mut p = BrlCadPanelState::new();
        p.text = "(sph 1.5)".into();
        p.evaluate().unwrap();
        assert!(p.result.is_some());
        assert!(p.last_status.is_some());
    }

    #[test]
    fn panel_evaluate_parse_error() {
        let mut p = BrlCadPanelState::new();
        p.text = "(zorp 1)".into();
        assert!(p.evaluate().is_err());
    }
}

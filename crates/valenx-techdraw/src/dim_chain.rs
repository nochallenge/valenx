//! Auto-dimension chains — sequences of related linear dimensions
//! generated from a set of points / edges with a common reference.
//!
//! Three chain kinds are supported per Phase 18B:
//!
//! - **Ordinate** — every dimension measures from a single shared
//!   origin (the first point). Common for hole patterns where each
//!   hole is called out relative to a known corner.
//! - **Baseline** — every dimension shares the same baseline edge but
//!   measures the next point's distance from it. Equivalent to
//!   ordinate when the baseline is axis-aligned, but the API keeps
//!   them distinct because the visual layout (multiple parallel dim
//!   lines vs. stacked dim lines) differs.
//! - **Chain** — end-to-end "running" dims: each dim measures the gap
//!   between two consecutive points.
//!
//! A [`DimChain`] is **pure data** — it carries the entries and the
//! offset to put the dim line at, and `expand` produces the concrete
//! [`crate::Dimension`] list at render time. This means a chain
//! survives serialization and re-renders identically across SVG /
//! PDF / DXF.

use serde::{Deserialize, Serialize};

use crate::dimension::Dimension;

/// Which arrangement to use when expanding [`DimChain::expand`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DimChainKind {
    /// Every dim measured from `entries[0]` (the origin).
    Ordinate,
    /// Every dim measured from the shared baseline point (also
    /// `entries[0]`), with stacked offsets so the dim lines don't
    /// overlap. Differs from [`Self::Ordinate`] only in visual layout.
    Baseline,
    /// Consecutive pair-wise dims: `(entries[0]→entries[1])`,
    /// `(entries[1]→entries[2])`, etc.
    Chain,
}

impl DimChainKind {
    /// Short label for UI dropdowns.
    pub fn label(self) -> &'static str {
        match self {
            DimChainKind::Ordinate => "Ordinate",
            DimChainKind::Baseline => "Baseline",
            DimChainKind::Chain => "Chain",
        }
    }
}

/// A run of related linear dimensions.
///
/// `entries` is the ordered list of measurement points in sheet
/// millimeters. `chain_kind` controls the expansion strategy; `offset`
/// is the perpendicular distance from the witness lines to the dim
/// line (positive = right of the witness, matching the convention used
/// by [`Dimension::Linear`]).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DimChain {
    /// Points to dimension between, in order.
    pub entries: Vec<[f64; 2]>,
    /// Layout strategy.
    pub chain_kind: DimChainKind,
    /// Perpendicular offset of the dim line(s) from the witness lines.
    pub offset: f64,
}

impl DimChain {
    /// Empty chain stub.
    pub fn new(chain_kind: DimChainKind, offset: f64) -> Self {
        Self {
            entries: Vec::new(),
            chain_kind,
            offset,
        }
    }

    /// Expand to a flat list of [`Dimension::Linear`] entries the
    /// renderer can consume directly.
    pub fn expand(&self) -> Vec<Dimension> {
        if self.entries.len() < 2 {
            return Vec::new();
        }
        match self.chain_kind {
            DimChainKind::Ordinate => self
                .entries
                .iter()
                .skip(1)
                .map(|p| {
                    let from = self.entries[0];
                    let value = distance(from, *p);
                    Dimension::Linear {
                        from,
                        to: *p,
                        offset: self.offset,
                        value,
                    }
                })
                .collect(),
            DimChainKind::Baseline => self
                .entries
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, p)| {
                    let from = self.entries[0];
                    let value = distance(from, *p);
                    // Stack each dim line by `offset * i` so visually
                    // they don't overlap. Same numeric value as
                    // Ordinate; layout only.
                    Dimension::Linear {
                        from,
                        to: *p,
                        offset: self.offset * i as f64,
                        value,
                    }
                })
                .collect(),
            DimChainKind::Chain => self
                .entries
                .windows(2)
                .map(|w| {
                    let value = distance(w[0], w[1]);
                    Dimension::Linear {
                        from: w[0],
                        to: w[1],
                        offset: self.offset,
                        value,
                    }
                })
                .collect(),
        }
    }
}

fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    (dx * dx + dy * dy).sqrt()
}

/// Build a chain from an existing list of edges (each edge a 2-point
/// segment) by taking the *first* endpoint of every edge as a chain
/// entry plus the last endpoint of the final edge.
///
/// This is the workhorse used by the UI's "Add chain" button when the
/// user has selected a contiguous strip of edges — turn the selection
/// into a [`DimChain`] without making the user click each endpoint.
pub fn auto_chain(edges: &[[(f64, f64); 2]], kind: DimChainKind, offset: f64) -> DimChain {
    let mut entries: Vec<[f64; 2]> = Vec::new();
    if let Some(first) = edges.first() {
        entries.push([first[0].0, first[0].1]);
    }
    for e in edges {
        entries.push([e[1].0, e[1].1]);
    }
    DimChain {
        entries,
        chain_kind: kind,
        offset,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pts(v: &[[f64; 2]]) -> Vec<[f64; 2]> {
        v.to_vec()
    }

    #[test]
    fn ordinate_chain_dims_share_origin() {
        let mut c = DimChain::new(DimChainKind::Ordinate, 5.0);
        c.entries = pts(&[[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]]);
        let dims = c.expand();
        assert_eq!(dims.len(), 2);
        match dims[0] {
            Dimension::Linear {
                from, to, value, ..
            } => {
                assert_eq!(from, [0.0, 0.0]);
                assert_eq!(to, [10.0, 0.0]);
                assert!((value - 10.0).abs() < 1e-9);
            }
            _ => panic!("expected Linear"),
        }
        match dims[1] {
            Dimension::Linear {
                from, to, value, ..
            } => {
                assert_eq!(from, [0.0, 0.0]);
                assert_eq!(to, [20.0, 0.0]);
                assert!((value - 20.0).abs() < 1e-9);
            }
            _ => panic!("expected Linear"),
        }
    }

    #[test]
    fn baseline_chain_stacks_offset_per_entry() {
        let mut c = DimChain::new(DimChainKind::Baseline, 3.0);
        c.entries = pts(&[[0.0, 0.0], [5.0, 0.0], [10.0, 0.0]]);
        let dims = c.expand();
        match dims[0] {
            Dimension::Linear { offset, .. } => assert!((offset - 3.0).abs() < 1e-9),
            _ => panic!("expected Linear"),
        }
        match dims[1] {
            Dimension::Linear { offset, .. } => assert!((offset - 6.0).abs() < 1e-9),
            _ => panic!("expected Linear"),
        }
    }

    #[test]
    fn running_chain_pairs_consecutive_points() {
        let mut c = DimChain::new(DimChainKind::Chain, 5.0);
        c.entries = pts(&[[0.0, 0.0], [3.0, 0.0], [7.0, 0.0]]);
        let dims = c.expand();
        assert_eq!(dims.len(), 2);
        match dims[0] {
            Dimension::Linear { from, to, .. } => {
                assert_eq!(from, [0.0, 0.0]);
                assert_eq!(to, [3.0, 0.0]);
            }
            _ => panic!("expected Linear"),
        }
        match dims[1] {
            Dimension::Linear { from, to, .. } => {
                assert_eq!(from, [3.0, 0.0]);
                assert_eq!(to, [7.0, 0.0]);
            }
            _ => panic!("expected Linear"),
        }
    }

    #[test]
    fn expand_empty_chain_returns_empty() {
        let c = DimChain::new(DimChainKind::Chain, 1.0);
        assert!(c.expand().is_empty());
    }

    #[test]
    fn auto_chain_uses_edge_endpoints() {
        let edges = vec![
            [(0.0, 0.0), (5.0, 0.0)],
            [(5.0, 0.0), (10.0, 0.0)],
            [(10.0, 0.0), (15.0, 0.0)],
        ];
        let c = auto_chain(&edges, DimChainKind::Chain, 2.0);
        assert_eq!(c.entries.len(), 4);
        assert_eq!(c.entries[0], [0.0, 0.0]);
        assert_eq!(c.entries[3], [15.0, 0.0]);
        let dims = c.expand();
        assert_eq!(dims.len(), 3);
    }

    #[test]
    fn chain_kind_labels_unique() {
        let labels = [
            DimChainKind::Ordinate.label(),
            DimChainKind::Baseline.label(),
            DimChainKind::Chain.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }
}

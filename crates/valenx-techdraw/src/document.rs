//! Top-level [`Drawing`] document.
//!
//! A drawing owns a [`Sheet`] (paper + title block), a list of
//! [`View`]s, and a list of [`Dimension`]s. The drawing itself does
//! not own the source 3D solids — each `View` carries a numeric
//! `source_solid_index`. The UI layer (or any caller) maintains a
//! parallel `Vec<Solid>` and resolves the index at generation /
//! render time.
//!
//! Why a numeric index instead of an `Arc<Solid>` reference? It
//! survives serialization round-trips (RON / JSON) without surprise,
//! and matches the way `Feature` already references sketches by
//! `SketchRef(usize)` in the feature tree.

use serde::{Deserialize, Serialize};

use crate::balloon::Balloon;
use crate::bom::Bom;
use crate::detail_view::DetailView;
use crate::dim_chain::DimChain;
use crate::dimension::Dimension;
use crate::error::TechDrawError;
use crate::gdt::{Datum, GdtSymbol};
use crate::leader::Leader;
use crate::parametric_view::{self, ParametricView};
use crate::projection_group::ProjectionGroup;
use crate::revision_block::RevisionBlock;
use crate::sheet::Sheet;
use crate::surface_finish::SurfaceFinish;
use crate::view::View;
use crate::weld::WeldSymbol;

/// A complete engineering drawing — sheet, views, dimensions, plus
/// Phase 18 overlays (balloons, leaders, weld / surface-finish / GD&T
/// symbols, dim chains, parametric-view links).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Drawing {
    /// The paper sheet (size + title-block metadata).
    pub sheet: Sheet,
    /// All views on the sheet, in insertion order. Index is the
    /// view's stable identifier within this drawing.
    pub views: Vec<View>,
    /// All dimensions overlaid on the views. Drawing-global rather
    /// than per-view so a single dimension can span multiple views
    /// (chained linear dimensions, for example).
    pub dimensions: Vec<Dimension>,
    /// Auto-dim chains (Phase 18B). Each chain expands to a sequence
    /// of [`Dimension`]s at render time so existing renderers keep
    /// working untouched.
    #[serde(default)]
    pub dim_chains: Vec<DimChain>,
    /// Balloons (Phase 18C) — item callouts referencing a BOM line.
    #[serde(default)]
    pub balloons: Vec<Balloon>,
    /// Leader lines (Phase 18C, refined in 18G) with arrowheads + text.
    #[serde(default)]
    pub leaders: Vec<Leader>,
    /// Weld symbols (Phase 18D) per ISO 2553.
    #[serde(default)]
    pub welds: Vec<WeldSymbol>,
    /// Surface-finish callouts (Phase 18E) per ISO 1302.
    #[serde(default)]
    pub surface_finishes: Vec<SurfaceFinish>,
    /// GD&T feature control frames (Phase 18F) per ASME Y14.5.
    #[serde(default)]
    pub gdt: Vec<GdtSymbol>,
    /// Datum-feature symbols (Phase 18F).
    #[serde(default)]
    pub datums: Vec<Datum>,
    /// Parametric-view links (Phase 18A). Each entry ties a [`View`]
    /// back to a [`valenx_feature_tree::feature::FeatureId`] in a
    /// feature tree so the view can re-extract edges when the tree
    /// replays.
    #[serde(default)]
    pub parametric_views: Vec<ParametricView>,
    /// Orthographic projection groups (Phase 19). Each group bundles
    /// Front / Top / Right (and optional Iso) views with auto-layout
    /// and optional parametric backing.
    #[serde(default)]
    pub projection_groups: Vec<ProjectionGroup>,
    /// Detail views (Phase 19) — magnified close-ups of a parent
    /// view's region.
    #[serde(default)]
    pub detail_views: Vec<DetailView>,
    /// BOM tables placed on the drawing (Phase 19). Each entry pairs
    /// a [`Bom`] payload with a sheet-space placement.
    #[serde(default)]
    pub bom_placements: Vec<BomPlacement>,
    /// Revision-history blocks (Phase 19). Usually exactly one, but
    /// the vec form makes it composable for multi-sheet drawings later.
    #[serde(default)]
    pub revision_blocks: Vec<RevisionBlock>,
}

/// A [`Bom`] placed at a specific position on the drawing sheet
/// (Phase 19).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BomPlacement {
    /// The BOM payload.
    pub bom: Bom,
    /// Lower-left corner of the rendered table in sheet mm.
    pub origin: [f64; 2],
}

impl BomPlacement {
    /// Construct a placement.
    pub fn new(bom: Bom, origin: [f64; 2]) -> Self {
        Self { bom, origin }
    }
}

impl Drawing {
    /// Empty drawing for the given sheet.
    pub fn new(sheet: Sheet) -> Self {
        Self {
            sheet,
            views: Vec::new(),
            dimensions: Vec::new(),
            dim_chains: Vec::new(),
            balloons: Vec::new(),
            leaders: Vec::new(),
            welds: Vec::new(),
            surface_finishes: Vec::new(),
            gdt: Vec::new(),
            datums: Vec::new(),
            parametric_views: Vec::new(),
            projection_groups: Vec::new(),
            detail_views: Vec::new(),
            bom_placements: Vec::new(),
            revision_blocks: Vec::new(),
        }
    }

    /// Regenerate every parametric view whose `auto_update` flag is
    /// `true`, plus every [`ProjectionGroup`] whose `feature_id` is
    /// set. The caller supplies a closure `solve_solid` that resolves
    /// a [`valenx_feature_tree::feature::FeatureId`] to the solid it
    /// should produce — this is decoupled from the feature tree to
    /// avoid pulling the entire tree as a dependency at every call
    /// site (the UI layer already keeps the resolved solids cached).
    ///
    /// Errors are collected per-view: a single bad view never blocks
    /// the others. The returned `Vec` carries `(view_idx, error)`
    /// entries for every view that failed (whether from a parametric
    /// view link or from a projection-group view).
    pub fn regenerate_all<F>(&mut self, mut solve_solid: F) -> Vec<(usize, TechDrawError)>
    where
        F: FnMut(valenx_feature_tree::feature::FeatureId) -> Option<valenx_cad::Solid>,
    {
        let mut errs = Vec::new();
        for pv in &self.parametric_views {
            if !pv.auto_update {
                continue;
            }
            let Some(solid) = solve_solid(pv.feature_id) else {
                continue;
            };
            if let Some(view) = self.views.get_mut(pv.view_idx) {
                if let Err(e) = parametric_view::regenerate(view, &solid) {
                    errs.push((pv.view_idx, e));
                }
            }
        }
        // Phase 19 — regenerate every projection group whose feature_id is set.
        // Clone the groups vector so we can mutably borrow `self.views`
        // inside `rebuild` while iterating.
        let groups_snapshot: Vec<ProjectionGroup> = self.projection_groups.clone();
        for g in &groups_snapshot {
            let Some(fid) = g.feature_id else {
                continue;
            };
            let Some(solid) = solve_solid(fid) else {
                continue;
            };
            errs.extend(g.rebuild(self, &solid));
        }
        errs
    }

    /// Append `view`, assign it the next stable id, return that id.
    pub fn add_view(&mut self, mut view: View) -> usize {
        let id = self.views.len();
        view.id = id;
        self.views.push(view);
        id
    }

    /// Remove the view at index `idx`. Later views shift down by one.
    /// Returns [`TechDrawError::UnknownView`] when the index is out
    /// of range.
    pub fn remove_view(&mut self, idx: usize) -> Result<View, TechDrawError> {
        if idx >= self.views.len() {
            return Err(TechDrawError::UnknownView(idx));
        }
        Ok(self.views.remove(idx))
    }

    /// Borrow a view by index.
    pub fn get_view(&self, idx: usize) -> Result<&View, TechDrawError> {
        self.views.get(idx).ok_or(TechDrawError::UnknownView(idx))
    }

    /// Mutable borrow of a view by index.
    pub fn get_view_mut(&mut self, idx: usize) -> Result<&mut View, TechDrawError> {
        let n = self.views.len();
        self.views.get_mut(idx).ok_or(TechDrawError::UnknownView(n))
    }

    /// Append a dimension. Returns its new index.
    pub fn add_dimension(&mut self, d: Dimension) -> usize {
        let id = self.dimensions.len();
        self.dimensions.push(d);
        id
    }

    /// Append a dimension chain (Phase 18B). Returns the new index.
    pub fn add_dim_chain(&mut self, c: DimChain) -> usize {
        let id = self.dim_chains.len();
        self.dim_chains.push(c);
        id
    }

    /// Append a balloon (Phase 18C).
    pub fn add_balloon(&mut self, b: Balloon) -> usize {
        let id = self.balloons.len();
        self.balloons.push(b);
        id
    }

    /// Append a leader (Phase 18C / 18G).
    pub fn add_leader(&mut self, l: Leader) -> usize {
        let id = self.leaders.len();
        self.leaders.push(l);
        id
    }

    /// Append a weld symbol (Phase 18D).
    pub fn add_weld(&mut self, w: WeldSymbol) -> usize {
        let id = self.welds.len();
        self.welds.push(w);
        id
    }

    /// Append a surface-finish symbol (Phase 18E).
    pub fn add_surface_finish(&mut self, s: SurfaceFinish) -> usize {
        let id = self.surface_finishes.len();
        self.surface_finishes.push(s);
        id
    }

    /// Append a GD&T feature control frame (Phase 18F).
    pub fn add_gdt(&mut self, g: GdtSymbol) -> usize {
        let id = self.gdt.len();
        self.gdt.push(g);
        id
    }

    /// Append a datum-feature symbol (Phase 18F).
    pub fn add_datum(&mut self, d: Datum) -> usize {
        let id = self.datums.len();
        self.datums.push(d);
        id
    }

    /// Register a parametric-view link (Phase 18A). The
    /// [`ParametricView::view_idx`] must already point at an existing
    /// view in `self.views`; the caller is responsible for that
    /// ordering (the helper does not validate to keep the link
    /// pure-data).
    pub fn add_parametric_view(&mut self, pv: ParametricView) -> usize {
        let id = self.parametric_views.len();
        self.parametric_views.push(pv);
        id
    }

    /// Build a [`ProjectionGroup`] from `group`, attach it to the
    /// drawing, and return the index. Mutates `group.view_indices` to
    /// record the freshly-created views — re-use this `ProjectionGroup`
    /// later by index to refer to its constituent views.
    ///
    /// Errors forward from [`ProjectionGroup::build_into`] (bad scale,
    /// empty solid). On failure no group is added.
    pub fn add_projection_group(
        &mut self,
        mut group: ProjectionGroup,
        solid: &valenx_cad::Solid,
    ) -> Result<usize, TechDrawError> {
        group.build_into(self, solid)?;
        let id = self.projection_groups.len();
        self.projection_groups.push(group);
        Ok(id)
    }

    /// Append a detail view. Auto-assigns its `id` to the next stable
    /// index in this drawing. If the caller's `label` is empty,
    /// auto-pick the next letter (A → B → …).
    pub fn add_detail_view(&mut self, mut dv: DetailView) -> usize {
        let id = self.detail_views.len();
        dv.id = id;
        if dv.label.is_empty() {
            dv.label = next_detail_label(self.detail_views.len());
        }
        self.detail_views.push(dv);
        id
    }

    /// Place a [`Bom`] on the drawing at `origin`.
    pub fn add_bom_placement(&mut self, bom: Bom, origin: [f64; 2]) -> usize {
        let id = self.bom_placements.len();
        self.bom_placements.push(BomPlacement::new(bom, origin));
        id
    }

    /// Add a revision-history block. The `position` of the block is
    /// taken as-is; callers commonly use [`RevisionBlock::standard_position`]
    /// to place it above the title block.
    pub fn add_revision_block(&mut self, blk: RevisionBlock) -> usize {
        let id = self.revision_blocks.len();
        self.revision_blocks.push(blk);
        id
    }

    /// Total number of annotation objects on this drawing — the sum of dimensions,
    /// dim-chains, balloons, leaders, weld symbols, surface-finish callouts, GD&T frames,
    /// and datum features. A drawing-complexity / annotation-density diagnostic, distinct
    /// from the view, BOM, and revision-block counts.
    pub fn annotation_count(&self) -> usize {
        self.dimensions.len()
            + self.dim_chains.len()
            + self.balloons.len()
            + self.leaders.len()
            + self.welds.len()
            + self.surface_finishes.len()
            + self.gdt.len()
            + self.datums.len()
    }
}

/// Map an integer to the next detail-view letter (`0 → "A"`, `25 →
/// "Z"`, `26 → "AA"`).
fn next_detail_label(n: usize) -> String {
    let mut k = n;
    let mut out: Vec<char> = Vec::new();
    loop {
        let digit = (k % 26) as u8;
        out.push((b'A' + digit) as char);
        k /= 26;
        if k == 0 {
            break;
        }
        k -= 1;
    }
    out.reverse();
    out.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sheet::Sheet;
    use crate::view::{View, ViewKind};

    #[test]
    fn empty_drawing_has_no_views_or_dimensions() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        assert!(d.views.is_empty());
        assert!(d.dimensions.is_empty());
    }

    #[test]
    fn add_view_assigns_sequential_ids() {
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let id_a = d.add_view(View::new(ViewKind::Front, 1.0, [0.0, 0.0]));
        let id_b = d.add_view(View::new(ViewKind::Top, 1.0, [0.0, 100.0]));
        assert_eq!(id_a, 0);
        assert_eq!(id_b, 1);
        assert_eq!(d.views[0].id, 0);
        assert_eq!(d.views[1].id, 1);
    }

    #[test]
    fn remove_view_unknown_returns_error() {
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let e = d.remove_view(7).unwrap_err();
        match e {
            TechDrawError::UnknownView(7) => {}
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn add_dimension_appends() {
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let id = d.add_dimension(Dimension::Linear {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 5.0,
            value: 10.0,
        });
        assert_eq!(id, 0);
        assert_eq!(d.dimensions.len(), 1);
    }

    /// Phase 18A Task 5 — editing the source solid then calling
    /// [`Drawing::regenerate_all`] updates a parametric view's
    /// visible_edges in place.
    #[test]
    fn parametric_view_auto_updates_when_solid_changes() {
        use crate::parametric_view::ParametricView;
        use valenx_cad::primitives::box_solid;
        use valenx_feature_tree::feature::FeatureId;

        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let idx = d.add_view(View::new(ViewKind::Front, 1.0, [10.0, 10.0]));
        d.add_parametric_view(ParametricView::new(idx, FeatureId(0)));

        // First "solid" is a 1×1×1 cube.
        let small = box_solid(1.0, 1.0, 1.0).unwrap();
        d.regenerate_all(|fid| {
            if fid == FeatureId(0) {
                Some(small.clone())
            } else {
                None
            }
        });
        // Edge bounds of a 1×1 front view → bbox is [0,1] in both axes.
        let bbox_small = d.views[0].bbox().unwrap();
        let w_small = bbox_small.1[0] - bbox_small.0[0];
        let h_small = bbox_small.1[1] - bbox_small.0[1];

        // Second "solid" is a 4×4×4 cube — same FeatureId, just a
        // bigger output. After regenerate, the view's edges should
        // reflect the new size.
        let big = box_solid(4.0, 4.0, 4.0).unwrap();
        d.regenerate_all(|fid| {
            if fid == FeatureId(0) {
                Some(big.clone())
            } else {
                None
            }
        });
        let bbox_big = d.views[0].bbox().unwrap();
        let w_big = bbox_big.1[0] - bbox_big.0[0];
        let h_big = bbox_big.1[1] - bbox_big.0[1];

        assert!(
            w_big > w_small * 2.0,
            "front view should grow with the source solid (got {w_small} → {w_big})"
        );
        assert!(
            h_big > h_small * 2.0,
            "front view height should grow too (got {h_small} → {h_big})"
        );
    }

    /// Parametric views with `auto_update = false` are skipped.
    #[test]
    fn parametric_view_with_auto_update_false_is_frozen() {
        use crate::parametric_view::ParametricView;
        use valenx_cad::primitives::box_solid;
        use valenx_feature_tree::feature::FeatureId;

        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let idx = d.add_view(View::new(ViewKind::Front, 1.0, [10.0, 10.0]));
        let mut pv = ParametricView::new(idx, FeatureId(0));
        pv.auto_update = false;
        d.add_parametric_view(pv);
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        d.regenerate_all(|_| Some(cube.clone()));
        // Edge lists stay empty because auto_update is off.
        assert!(d.views[0].visible_edges.is_empty());
    }

    /// Phase 19 — `add_projection_group` registers a group + its
    /// generated views.
    #[test]
    fn add_projection_group_attaches_group_and_views() {
        use crate::projection_group::{Projection, ProjectionGroup};
        use valenx_cad::primitives::box_solid;
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let cube = box_solid(20.0, 20.0, 20.0).unwrap();
        let group = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        let gid = d.add_projection_group(group, &cube).unwrap();
        assert_eq!(gid, 0);
        assert_eq!(d.projection_groups.len(), 1);
        assert_eq!(d.views.len(), 3);
        assert_eq!(d.projection_groups[0].view_indices.len(), 3);
    }

    /// Phase 19 — `regenerate_all` rebuilds projection-group views
    /// when the backing FeatureId is set.
    #[test]
    fn regenerate_all_rebuilds_projection_group_when_feature_changes() {
        use crate::projection_group::{Projection, ProjectionGroup};
        use valenx_cad::primitives::box_solid;
        use valenx_feature_tree::feature::FeatureId;

        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let small = box_solid(10.0, 10.0, 10.0).unwrap();
        let big = box_solid(50.0, 50.0, 50.0).unwrap();
        let mut g = ProjectionGroup::new([150.0, 150.0], 1.0, Projection::ThirdAngle);
        g.feature_id = Some(FeatureId(0));
        d.add_projection_group(g, &small).unwrap();
        let w_before = d.views[0].bbox().unwrap().1[0] - d.views[0].bbox().unwrap().0[0];
        // Replay with the bigger solid.
        d.regenerate_all(|fid| if fid == FeatureId(0) { Some(big.clone()) } else { None });
        let w_after = d.views[0].bbox().unwrap().1[0] - d.views[0].bbox().unwrap().0[0];
        assert!(w_after > w_before * 2.0, "front view should grow ({w_before} → {w_after})");
    }

    /// Phase 19 — auto-assigned detail labels: first → "A", second → "B".
    #[test]
    fn add_detail_view_auto_assigns_labels() {
        use crate::detail_view::DetailView;
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let id_a = d.add_detail_view(DetailView::new(0, [10.0, 10.0], 5.0, [100.0, 100.0], 2.0, ""));
        let id_b = d.add_detail_view(DetailView::new(0, [20.0, 20.0], 5.0, [150.0, 100.0], 4.0, ""));
        assert_eq!(d.detail_views[id_a].label, "A");
        assert_eq!(d.detail_views[id_b].label, "B");
    }

    /// Phase 19 — manual label is preserved.
    #[test]
    fn add_detail_view_keeps_caller_label() {
        use crate::detail_view::DetailView;
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let id = d.add_detail_view(DetailView::new(0, [10.0, 10.0], 5.0, [100.0, 100.0], 2.0, "Z"));
        assert_eq!(d.detail_views[id].label, "Z");
    }

    /// Phase 19 — `add_bom_placement` stores BOM + origin.
    #[test]
    fn add_bom_placement_stores_bom_and_origin() {
        use crate::bom::{Bom, BomItem};
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let mut bom = Bom::new();
        bom.add(BomItem::new("Foo", 1));
        bom.renumber_items();
        let id = d.add_bom_placement(bom, [10.0, 100.0]);
        assert_eq!(id, 0);
        assert_eq!(d.bom_placements[0].origin, [10.0, 100.0]);
        assert_eq!(d.bom_placements[0].bom.items[0].item_number, 1);
    }

    /// Phase 19 — `add_revision_block` stores a revision history.
    #[test]
    fn add_revision_block_stores_entries() {
        use crate::revision_block::{RevisionBlock, RevisionEntry};
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let mut blk = RevisionBlock::new([10.0, 70.0]);
        blk.add_entry(RevisionEntry::new("A", "2026-05-23", "initial", "GH", ""));
        d.add_revision_block(blk);
        assert_eq!(d.revision_blocks.len(), 1);
        assert_eq!(d.revision_blocks[0].entries[0].rev, "A");
    }

    #[test]
    fn annotation_count_sums_annotations() {
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        assert_eq!(d.annotation_count(), 0);
        d.add_dimension(Dimension::Linear {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 5.0,
            value: 10.0,
        });
        d.add_dimension(Dimension::Linear {
            from: [0.0, 0.0],
            to: [0.0, 20.0],
            offset: 5.0,
            value: 20.0,
        });
        // Two dimensions and no other annotation type → 2.
        assert_eq!(d.annotation_count(), 2);
        assert_eq!(d.annotation_count(), d.dimensions.len());
    }
}

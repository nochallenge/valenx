//! RON-based persistence for [`Drawing`].
//!
//! Identical envelope shape to `valenx-draft::DraftFile`: an outer
//! `version: u32` plus the payload, written / read via
//! `ron::ser::to_string_pretty` and `ron::from_str`.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::document::Drawing;
use crate::error::TechDrawError;

/// On-disk envelope wrapping a drawing with a format-version tag.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TechDrawFile {
    /// Format version — bumped when on-disk schema changes.
    pub version: u32,
    /// The drawing payload.
    pub drawing: Drawing,
}

impl TechDrawFile {
    /// Current on-disk format version.
    ///
    /// **Version 2** (Phase 18): adds `parametric_views`, `balloons`,
    /// `leaders`, `welds`, `surface_finishes`, `gdt`, `datums`,
    /// `dim_chains` fields. All new fields are gated by
    /// `#[serde(default)]` so v1 RON files round-trip without
    /// migration.
    ///
    /// **Version 3** (Phase 19): adds `projection_groups`,
    /// `detail_views`, `bom_placements`, `revision_blocks`, plus the
    /// `BomItem` extended columns (`item_number`, `part_number`,
    /// `description`). All new fields use `#[serde(default)]` so v1 +
    /// v2 RON files round-trip cleanly.
    pub const VERSION: u32 = 3;

    /// Wrap a drawing with the current version tag.
    pub fn from_drawing(drawing: &Drawing) -> Self {
        Self {
            version: Self::VERSION,
            drawing: drawing.clone(),
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, TechDrawError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| TechDrawError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), TechDrawError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, TechDrawError> {
        ron::from_str(s).map_err(|e| TechDrawError::Ron(e.to_string()))
    }

    /// Read from a file.
    pub fn read_from(path: &Path) -> Result<Self, TechDrawError> {
        // R29 D: canonical valenx_core::io_caps::read_capped_to_string at
        // MAX_DOC_FILE_BYTES (16 MiB), replacing the private dupe.
        let s = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?;
        Self::from_ron(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dimension::Dimension;
    use crate::sheet::Sheet;
    use crate::view::{View, ViewKind};

    #[test]
    fn round_trips_empty_drawing() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let f = TechDrawFile::from_drawing(&d);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 3"));
        let back = TechDrawFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 3);
        assert_eq!(back.drawing.views.len(), 0);
    }

    /// Phase 18 — Phase 5 v1 RON files (no parametric_views /
    /// balloons / etc) must still parse cleanly under the v2 schema.
    #[test]
    fn legacy_v1_ron_round_trips() {
        // v1 RON shape: only `sheet`, `views`, `dimensions` exist.
        let legacy = r#"(
    version: 1,
    drawing: (
        sheet: (
            size: A4,
            title: "Legacy",
            author: "Me",
            revision: "A",
        ),
        views: [],
        dimensions: [],
    ),
)"#;
        let back = TechDrawFile::from_ron(legacy).unwrap();
        assert_eq!(back.version, 1);
        assert!(back.drawing.parametric_views.is_empty());
        assert!(back.drawing.balloons.is_empty());
        assert!(back.drawing.welds.is_empty());
    }

    /// Phase 18 — Round-trip every new annotation type to make sure
    /// the v2 envelope keeps them.
    #[test]
    fn v2_round_trips_all_new_annotations() {
        use crate::balloon::Balloon;
        use crate::dim_chain::{DimChain, DimChainKind};
        use crate::gdt::{Datum, DatumRef, GdtSymbol, GeometricCharacteristic};
        use crate::leader::{ArrowKind, Leader};
        use crate::parametric_view::ParametricView;
        use crate::surface_finish::SurfaceFinish;
        use crate::weld::WeldSymbol;
        use valenx_feature_tree::feature::FeatureId;

        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        d.add_view(View::new(ViewKind::Front, 1.0, [10.0, 10.0]));
        d.parametric_views
            .push(ParametricView::new(0, FeatureId(5)));
        d.balloons.push(Balloon::new([5.0, 5.0], "1", [10.0, 10.0]));
        let mut leader = Leader::new([0.0, 0.0], [10.0, 0.0], "note");
        leader.arrow_kind = ArrowKind::Open;
        d.leaders.push(leader);
        d.welds
            .push(WeldSymbol::new_fillet([10.0, 10.0], [30.0, 30.0], "8"));
        d.surface_finishes
            .push(SurfaceFinish::new([20.0, 20.0], 1.6));
        let mut g = GdtSymbol::new([40.0, 40.0], GeometricCharacteristic::Position, "0.1");
        g.datums.push(DatumRef::new("A"));
        d.gdt.push(g);
        d.datums.push(Datum::new([50.0, 50.0], "A", [60.0, 60.0]));
        let mut c = DimChain::new(DimChainKind::Chain, 5.0);
        c.entries = vec![[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]];
        d.dim_chains.push(c);

        let tmp = std::env::temp_dir().join("valenx_techdraw_v2.ron");
        TechDrawFile::from_drawing(&d).write_to(&tmp).unwrap();
        let back = TechDrawFile::read_from(&tmp).unwrap();
        assert_eq!(back.drawing, d);
        let _ = std::fs::remove_file(&tmp);
    }

    /// Task 38 — full round-trip with views (all 8 ViewKinds) and
    /// dimensions of every variant.
    #[test]
    fn round_trips_full_drawing_via_disk() {
        let mut d = Drawing::new(Sheet::a4_landscape("Bracket", "A. Engineer", "A"));
        d.add_view(View::new(ViewKind::Front, 1.0, [50.0, 50.0]));
        d.add_view(View::new(ViewKind::Top, 0.5, [50.0, 150.0]));
        d.add_view(View::new(ViewKind::Right, 1.0, [200.0, 50.0]));
        d.add_view(View::new(ViewKind::Back, 1.0, [200.0, 150.0]));
        d.add_view(View::new(ViewKind::Bottom, 1.0, [50.0, 100.0]));
        d.add_view(View::new(ViewKind::Left, 1.0, [200.0, 100.0]));
        d.add_view(View::new(ViewKind::Isometric, 1.0, [300.0, 100.0]));
        d.add_view(View::new(
            ViewKind::Custom {
                eye: nalgebra::Vector3::new(50.0, 50.0, 50.0),
                target: nalgebra::Vector3::zeros(),
                up: nalgebra::Vector3::new(0.0, 0.0, 1.0),
            },
            1.0,
            [350.0, 100.0],
        ));
        d.add_dimension(Dimension::Linear {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 5.0,
            value: 10.0,
        });
        d.add_dimension(Dimension::Angular {
            vertex: [0.0, 0.0],
            a: [1.0, 0.0],
            b: [0.0, 1.0],
            offset: 5.0,
            value: 90.0,
        });
        d.add_dimension(Dimension::Radial {
            center: [0.0, 0.0],
            radius: 3.0,
            label_pos: [10.0, 0.0],
            value: 3.0,
        });
        d.add_dimension(Dimension::Diameter {
            center: [0.0, 0.0],
            radius: 4.0,
            label_pos: [10.0, 0.0],
            value: 8.0,
        });

        let tmp = std::env::temp_dir().join("valenx_techdraw_full.ron");
        TechDrawFile::from_drawing(&d).write_to(&tmp).unwrap();
        let back = TechDrawFile::read_from(&tmp).unwrap();
        assert_eq!(back.drawing, d);
        let _ = std::fs::remove_file(&tmp);
    }

    /// Phase 19 — v3 round-trips projection groups + detail views +
    /// BOM placements + revision blocks.
    #[test]
    fn v3_round_trips_phase19_artifacts() {
        use crate::bom::{Bom, BomItem};
        use crate::detail_view::DetailView;
        use crate::projection_group::{Projection, ProjectionGroup};
        use crate::revision_block::{RevisionBlock, RevisionEntry};

        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        // Projection group (just the data — no live views here so we
        // don't depend on solid edge extraction in a persist test).
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        g.include_isometric = true;
        g.view_indices = vec![0, 1, 2, 3];
        d.projection_groups.push(g);
        // Detail view.
        d.detail_views.push(DetailView::new(0, [10.0, 10.0], 5.0, [200.0, 100.0], 4.0, "A"));
        // BOM with extended columns.
        let mut bom = Bom::new();
        bom.add(BomItem::full("Bracket", 2, "P-001", "desc", "Steel"));
        bom.renumber_items();
        d.add_bom_placement(bom, [10.0, 70.0]);
        // Revision block.
        let mut blk = RevisionBlock::new([100.0, 70.0]);
        blk.add_entry(RevisionEntry::new("A", "2026-05-23", "init", "GH", ""));
        d.add_revision_block(blk);

        let ron = TechDrawFile::from_drawing(&d).to_ron().unwrap();
        assert!(ron.contains("version: 3"));
        let back = TechDrawFile::from_ron(&ron).unwrap();
        assert_eq!(back.drawing, d);
    }

    #[test]
    fn bad_ron_returns_ron_error() {
        let bad = "not a ron document";
        match TechDrawFile::from_ron(bad).unwrap_err() {
            TechDrawError::Ron(_) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }
}

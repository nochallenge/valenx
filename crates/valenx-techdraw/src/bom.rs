//! Bill of materials (BOM).
//!
//! A [`Bom`] is an ordered list of [`BomItem`]s — one row per
//! distinct part. Counts are aggregated when building from a feature
//! tree (each feature kind becomes one row with `quantity` equal to
//! the number of instances of that kind).
//!
//! [`Bom::to_text`] renders a fixed-width text table suitable for
//! pasting into a title-block annotation, the log, or a CI artifact.

use serde::{Deserialize, Serialize};

use valenx_feature_tree::feature::Feature;

/// One line of the bill of materials.
///
/// Extended (Phase 19) with the standard drawing-table columns:
/// `part_number` and `description`. Both default to empty so existing
/// callers + serialized data round-trip unchanged. The `item_number`
/// is auto-assigned by [`Bom::renumber_items`] / [`Bom::from_parts`]
/// — callers should not set it directly.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BomItem {
    /// Human-readable part name (e.g. `"Pad"`, `"Bracket"`).
    pub part_name: String,
    /// Number of instances of this part.
    pub quantity: u32,
    /// Material specification (e.g. `"Steel 1018"`, `"Al 6061-T6"`).
    /// Blank by default — the UI fills it in.
    pub material: String,
    /// Drawing item number (rendered as the first column on the
    /// table). Auto-assigned starting at 1 by
    /// [`Bom::renumber_items`].
    #[serde(default)]
    pub item_number: u32,
    /// Part number (drawing column, e.g. `"P-1001"`). Blank by
    /// default.
    #[serde(default)]
    pub part_number: String,
    /// Description (drawing column, e.g. `"Mounting bracket, 6061-T6"`).
    /// Blank by default.
    #[serde(default)]
    pub description: String,
}

impl BomItem {
    /// Convenience constructor with an empty material / part-number /
    /// description.
    pub fn new(part_name: &str, quantity: u32) -> Self {
        Self {
            part_name: part_name.into(),
            quantity,
            material: String::new(),
            item_number: 0,
            part_number: String::new(),
            description: String::new(),
        }
    }

    /// Full-row constructor matching the standard drawing-table
    /// columns. `item_number` is set to 0 — let [`Bom::renumber_items`]
    /// assign it when the row joins a [`Bom`].
    pub fn full(
        part_name: &str,
        quantity: u32,
        part_number: &str,
        description: &str,
        material: &str,
    ) -> Self {
        Self {
            part_name: part_name.into(),
            quantity,
            material: material.into(),
            item_number: 0,
            part_number: part_number.into(),
            description: description.into(),
        }
    }
}

/// A complete bill of materials.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Default)]
pub struct Bom {
    /// Rows, in insertion order.
    pub items: Vec<BomItem>,
}

impl Bom {
    /// Empty BOM.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a BOM row.
    pub fn add(&mut self, item: BomItem) {
        self.items.push(item);
    }

    /// Build a BOM by counting feature kinds. Each distinct
    /// [`Feature::kind_label`] becomes one row; quantity = count of
    /// occurrences in the input list.
    ///
    /// Rows are returned in the order each kind first appears so the
    /// BOM is deterministic regardless of `HashMap` iteration order.
    pub fn from_features(features: &[Feature]) -> Self {
        let mut order: Vec<&'static str> = Vec::new();
        let mut counts: std::collections::HashMap<&'static str, u32> =
            std::collections::HashMap::new();
        for f in features {
            let k = f.kind_label();
            let entry = counts.entry(k).or_insert_with(|| {
                order.push(k);
                0
            });
            *entry += 1;
        }
        let items = order
            .into_iter()
            .map(|k| BomItem::new(k, counts[k]))
            .collect();
        let mut b = Self { items };
        b.renumber_items();
        b
    }

    /// Renumber every row's `item_number` starting at 1 in current
    /// order. Called automatically by [`Bom::from_parts`] and
    /// [`Bom::from_assembly_parts`]; available manually so callers
    /// who mutate `items` directly can keep the column consistent.
    pub fn renumber_items(&mut self) {
        for (i, it) in self.items.iter_mut().enumerate() {
            it.item_number = (i + 1) as u32;
        }
    }

    /// Build a BOM from a caller-supplied list of parts. Each tuple is
    /// `(part_name, quantity, part_number, description, material)`.
    /// Rows are renumbered starting at 1 in input order.
    pub fn from_parts(parts: &[(&str, u32, &str, &str, &str)]) -> Self {
        let mut b = Self::new();
        for (name, qty, pn, desc, mat) in parts {
            b.items.push(BomItem::full(name, *qty, pn, desc, mat));
        }
        b.renumber_items();
        b
    }

    /// Build a BOM from a list of [`valenx_assembly::Part`]s. Parts
    /// sharing the same `name` are aggregated into a single row whose
    /// quantity is the count. Materials, part numbers, and
    /// descriptions are left blank — the UI / caller fills them in
    /// after import (assembly parts don't carry those fields).
    pub fn from_assembly_parts(parts: &[valenx_assembly::Part]) -> Self {
        let mut order: Vec<String> = Vec::new();
        let mut counts: std::collections::HashMap<String, u32> = std::collections::HashMap::new();
        for p in parts {
            let k = p.name.clone();
            counts.entry(k.clone()).or_insert_with(|| {
                order.push(k.clone());
                0
            });
            *counts.get_mut(&k).unwrap() += 1;
        }
        let mut b = Self::new();
        for name in order {
            let qty = counts[&name];
            b.items.push(BomItem::new(&name, qty));
        }
        b.renumber_items();
        b
    }

    /// Render the BOM as line + label segments suitable for inclusion
    /// in a drawing. `origin` is the **lower-left corner** of the
    /// table in sheet mm; the table grows up + right.
    ///
    /// Standard column widths (mm):
    /// - Item       — 12
    /// - Qty        — 12
    /// - PartNumber — 32
    /// - Description — 60
    /// - Material   — 40
    ///
    /// Total table width: 156 mm. Row height is 7 mm. The header row
    /// is one row tall and labelled.
    ///
    /// Returns `(grid_segments, labels)` where `labels` is `(x, y,
    /// text)` tuples for the cell content. The renderer is responsible
    /// for choosing font size / clipping long strings.
    pub fn render_table(&self, origin: [f64; 2]) -> crate::revision_block::RenderedTable {
        let col_widths = [12.0, 12.0, 32.0, 60.0, 40.0];
        let row_h = 7.0;
        let header_h = 7.0;
        let total_width: f64 = col_widths.iter().sum();
        let n_rows = self.items.len();
        // +1 for the header.
        let total_height = header_h + row_h * n_rows as f64;
        let [ox, oy] = origin;
        let mut grid: Vec<[(f64, f64); 2]> = vec![
            // Outer rectangle.
            [(ox, oy), (ox + total_width, oy)],
            [
                (ox + total_width, oy),
                (ox + total_width, oy + total_height),
            ],
            [
                (ox + total_width, oy + total_height),
                (ox, oy + total_height),
            ],
            [(ox, oy + total_height), (ox, oy)],
        ];
        // Column dividers.
        let mut acc = 0.0;
        for w in &col_widths[..col_widths.len() - 1] {
            acc += w;
            grid.push([(ox + acc, oy), (ox + acc, oy + total_height)]);
        }
        // Row dividers (one above each data row + one above the header).
        // Header sits at the *top* of the table; data rows below.
        let header_y = oy + total_height - header_h;
        grid.push([(ox, header_y), (ox + total_width, header_y)]);
        for i in 1..n_rows {
            // Row i (counting from header below it) sits between
            // (top - header - i*row_h) and (top - header - (i+1)*row_h).
            let y = header_y - row_h * i as f64;
            grid.push([(ox, y), (ox + total_width, y)]);
        }

        // Labels.
        let mut labels: Vec<(f64, f64, String)> = Vec::new();
        let pad_x = 1.0;
        let baseline_y_offset = row_h * 0.3;
        let headers = ["Item", "Qty", "Part No.", "Description", "Material"];
        let mut cx = ox;
        for (i, h) in headers.iter().enumerate() {
            labels.push((cx + pad_x, header_y + baseline_y_offset, (*h).into()));
            cx += col_widths[i];
        }
        // Data rows: row 0 is the *topmost* data row (just below the header).
        for (row_idx, it) in self.items.iter().enumerate() {
            let y = header_y - row_h * (row_idx + 1) as f64 + baseline_y_offset;
            let cells = [
                it.item_number.to_string(),
                it.quantity.to_string(),
                it.part_number.clone(),
                if it.description.is_empty() {
                    it.part_name.clone()
                } else {
                    it.description.clone()
                },
                it.material.clone(),
            ];
            let mut cx = ox;
            for (i, txt) in cells.iter().enumerate() {
                labels.push((cx + pad_x, y, txt.clone()));
                cx += col_widths[i];
            }
        }

        (grid, labels)
    }

    /// Render as a fixed-width text table. Three columns — part /
    /// qty / material — with column widths computed from the widest
    /// row plus the header.
    pub fn to_text(&self) -> String {
        let header = ("Part", "Qty", "Material");
        let mut w_part = header.0.len();
        let mut w_qty = header.1.len();
        let mut w_mat = header.2.len();
        for it in &self.items {
            w_part = w_part.max(it.part_name.len());
            w_qty = w_qty.max(it.quantity.to_string().len());
            w_mat = w_mat.max(it.material.len());
        }
        let mut out = String::new();
        // Header row.
        out.push_str(&format!(
            "{:<wp$}  {:>wq$}  {:<wm$}\n",
            header.0,
            header.1,
            header.2,
            wp = w_part,
            wq = w_qty,
            wm = w_mat
        ));
        // Separator.
        out.push_str(&format!(
            "{}  {}  {}\n",
            "-".repeat(w_part),
            "-".repeat(w_qty),
            "-".repeat(w_mat),
        ));
        // Rows.
        for it in &self.items {
            out.push_str(&format!(
                "{:<wp$}  {:>wq$}  {:<wm$}\n",
                it.part_name,
                it.quantity,
                it.material,
                wp = w_part,
                wq = w_qty,
                wm = w_mat
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_feature_tree::feature::{
        Feature, FeatureId, MirrorParams, PadParams, PocketParams, SketchRef,
    };

    #[test]
    fn from_features_counts_each_kind() {
        let pad = Feature::Pad(PadParams {
            sketch: SketchRef(0),
            depth: 10.0.into(),
            direction_positive: true,
        });
        let pocket = Feature::Pocket(PocketParams {
            sketch: SketchRef(1),
            depth: 5.0.into(),
            direction_positive: false,
        });
        let mirror = Feature::Mirror(MirrorParams {
            target: FeatureId(0),
            plane_origin: Vector3::zeros(),
            plane_normal: Vector3::new(1.0, 0.0, 0.0),
            keep_original: true,
        });
        let bom = Bom::from_features(&[pad.clone(), pad, pocket, mirror]);
        // Three distinct kinds: Pad (×2), Pocket (×1), Mirror (×1).
        assert_eq!(bom.items.len(), 3);
        let pad_row = bom.items.iter().find(|i| i.part_name == "Pad").unwrap();
        assert_eq!(pad_row.quantity, 2);
        let mirror_row = bom.items.iter().find(|i| i.part_name == "Mirror").unwrap();
        assert_eq!(mirror_row.quantity, 1);
    }

    #[test]
    fn from_features_preserves_first_seen_order() {
        let pocket = Feature::Pocket(PocketParams {
            sketch: SketchRef(1),
            depth: 5.0.into(),
            direction_positive: true,
        });
        let pad = Feature::Pad(PadParams {
            sketch: SketchRef(0),
            depth: 10.0.into(),
            direction_positive: true,
        });
        let bom = Bom::from_features(&[pocket, pad]);
        assert_eq!(bom.items[0].part_name, "Pocket");
        assert_eq!(bom.items[1].part_name, "Pad");
    }

    #[test]
    fn to_text_has_header_separator_and_rows() {
        let mut bom = Bom::new();
        bom.add(BomItem {
            part_name: "Bracket".into(),
            quantity: 2,
            material: "Steel".into(),
            item_number: 0,
            part_number: String::new(),
            description: String::new(),
        });
        bom.add(BomItem::new("Pin", 4));
        let txt = bom.to_text();
        let lines: Vec<_> = txt.lines().collect();
        assert!(lines[0].contains("Part"));
        assert!(lines[0].contains("Qty"));
        assert!(lines[0].contains("Material"));
        assert!(lines[1].starts_with('-'));
        assert!(lines[2].contains("Bracket"));
        assert!(lines[2].contains("Steel"));
        assert!(lines[3].contains("Pin"));
    }

    #[test]
    fn empty_bom_emits_header_only() {
        let bom = Bom::new();
        let txt = bom.to_text();
        assert!(txt.contains("Part"));
        assert_eq!(txt.lines().count(), 2);
    }

    /// Phase 19 — `renumber_items` assigns 1..N in order.
    #[test]
    fn renumber_items_assigns_sequential_ids() {
        let mut b = Bom::new();
        b.add(BomItem::new("Foo", 1));
        b.add(BomItem::new("Bar", 2));
        b.add(BomItem::new("Baz", 3));
        b.renumber_items();
        assert_eq!(b.items[0].item_number, 1);
        assert_eq!(b.items[1].item_number, 2);
        assert_eq!(b.items[2].item_number, 3);
    }

    /// Phase 19 — `from_parts` builds a BOM with renumbered rows and
    /// every column populated.
    #[test]
    fn from_parts_builds_renumbered_table() {
        let b = Bom::from_parts(&[
            ("Bracket", 2, "P-1001", "Mounting bracket", "Al 6061-T6"),
            ("Pin", 4, "P-2034", "Dowel pin Ø6 × 20", "Steel 1018"),
        ]);
        assert_eq!(b.items.len(), 2);
        assert_eq!(b.items[0].item_number, 1);
        assert_eq!(b.items[0].part_number, "P-1001");
        assert_eq!(b.items[0].description, "Mounting bracket");
        assert_eq!(b.items[0].material, "Al 6061-T6");
        assert_eq!(b.items[1].item_number, 2);
        assert_eq!(b.items[1].quantity, 4);
    }

    /// Phase 19 — `from_assembly_parts` aggregates duplicates by name.
    #[test]
    fn from_assembly_parts_aggregates_by_name() {
        use valenx_cad::primitives::box_solid;
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let parts = vec![
            valenx_assembly::Part::new(0, "Bracket", cube.clone()),
            valenx_assembly::Part::new(1, "Pin", cube.clone()),
            valenx_assembly::Part::new(2, "Bracket", cube.clone()),
            valenx_assembly::Part::new(3, "Pin", cube.clone()),
            valenx_assembly::Part::new(4, "Bracket", cube),
        ];
        let b = Bom::from_assembly_parts(&parts);
        // Two unique kinds.
        assert_eq!(b.items.len(), 2);
        let bracket = b.items.iter().find(|i| i.part_name == "Bracket").unwrap();
        assert_eq!(bracket.quantity, 3);
        let pin = b.items.iter().find(|i| i.part_name == "Pin").unwrap();
        assert_eq!(pin.quantity, 2);
        // Renumbered.
        assert!(b.items.iter().all(|i| i.item_number > 0));
    }

    /// Phase 19 — `render_table` emits a 5-column outline + per-row +
    /// per-column dividers.
    #[test]
    fn render_table_emits_grid_and_labels() {
        let b = Bom::from_parts(&[
            ("Bracket", 2, "P-1001", "Mounting bracket", "Al 6061-T6"),
            ("Pin", 4, "P-2034", "Dowel pin", "Steel"),
        ]);
        let (grid, labels) = b.render_table([10.0, 10.0]);
        // 4 outer + 4 column dividers + 2 row dividers (1 header + 1 between rows) = 10.
        // Outer (4) + column dividers (4) + header divider + (n_rows - 1) row dividers.
        // With n_rows = 2 → 4 + 4 + 1 + 1 = 10.
        assert_eq!(grid.len(), 10, "expected 10 grid lines, got {}", grid.len());
        // Labels: 5 headers + 5 per row × 2 rows = 15.
        assert_eq!(labels.len(), 15);
        // First label is the "Item" header.
        assert!(labels.iter().any(|(_, _, t)| t == "Item"));
        assert!(labels
            .iter()
            .any(|(_, _, t)| t == "Bracket" || t == "Mounting bracket"));
    }

    /// Phase 19 — `BomItem::full` populates every column.
    #[test]
    fn bom_item_full_populates_every_column() {
        let it = BomItem::full("Foo", 3, "PN-1", "desc", "Steel");
        assert_eq!(it.part_name, "Foo");
        assert_eq!(it.quantity, 3);
        assert_eq!(it.part_number, "PN-1");
        assert_eq!(it.description, "desc");
        assert_eq!(it.material, "Steel");
        assert_eq!(it.item_number, 0); // unset until renumber
    }
}

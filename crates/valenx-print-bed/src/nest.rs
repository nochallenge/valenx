//! 2D bin packing for the printer bed.

use nalgebra::Vector3;

use crate::error::PrintBedError;
use crate::printer::{Part, Printer};

/// Run a coarse first-fit decreasing 2D bin pack on `parts`. Returns
/// the same parts with `bed_position` updated to non-overlapping
/// (x, y) coordinates inside the printer bed. Parts that don't fit
/// surface as [`PrintBedError::PartTooLarge`].
///
/// The algorithm:
/// 1. Compute each part's XY-projected bounding box.
/// 2. Sort by descending bbox area.
/// 3. Greedy row-fill: place each part at the next free slot in the
///    current row; when the next slot would exceed bed width, advance
///    to a new row whose height = the previous row's tallest part.
pub fn auto_pack(parts: Vec<Part>, printer: &Printer) -> Result<Vec<Part>, PrintBedError> {
    let (bed_w, bed_d, _) = printer.bed_size;
    if bed_w <= 0.0 || bed_d <= 0.0 {
        return Err(PrintBedError::BadParameter {
            name: "printer.bed_size",
            reason: "x and y dimensions must be > 0".into(),
        });
    }

    // Compute footprints.
    let mut indexed: Vec<(usize, (f64, f64))> = parts
        .iter()
        .enumerate()
        .map(|(i, p)| (i, footprint(&p.mesh.nodes)))
        .collect();
    // Sort indices by descending footprint area.
    indexed.sort_by(|a, b| {
        let aa = a.1 .0 * a.1 .1;
        let bb = b.1 .0 * b.1 .1;
        bb.partial_cmp(&aa).unwrap_or(std::cmp::Ordering::Equal)
    });

    // First pass: validate each part fits.
    for (i, (w, h)) in &indexed {
        if *w > bed_w || *h > bed_d {
            return Err(PrintBedError::PartTooLarge {
                name: parts[*i].name.clone(),
                w: *w,
                h: *h,
                bw: bed_w,
                bh: bed_d,
            });
        }
    }

    // Place.
    let mut placed = parts.clone();
    let mut cursor_x = 0.0;
    let mut cursor_y = 0.0;
    let mut row_height = 0.0;
    for (i, (w, h)) in indexed {
        if cursor_x + w > bed_w {
            cursor_x = 0.0;
            cursor_y += row_height;
            row_height = 0.0;
        }
        if cursor_y + h > bed_d {
            // Out of vertical room — flag the offending part.
            return Err(PrintBedError::PartTooLarge {
                name: placed[i].name.clone(),
                w,
                h,
                bw: bed_w,
                bh: bed_d,
            });
        }
        placed[i].bed_position = [cursor_x + 0.5 * w, cursor_y + 0.5 * h];
        cursor_x += w;
        if h > row_height {
            row_height = h;
        }
    }
    Ok(placed)
}

/// XY-projected bounding box of a point cloud, returned as
/// `(width, depth)` in mm.
fn footprint(nodes: &[Vector3<f64>]) -> (f64, f64) {
    if nodes.is_empty() {
        return (0.0, 0.0);
    }
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for n in nodes {
        if n.x < min_x {
            min_x = n.x;
        }
        if n.y < min_y {
            min_y = n.y;
        }
        if n.x > max_x {
            max_x = n.x;
        }
        if n.y > max_y {
            max_y = n.y;
        }
    }
    (max_x - min_x, max_y - min_y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::printer::{BedMaterial, BedType};
    use valenx_mesh::Mesh;

    fn cube_part(name: &str, side: f64) -> Part {
        let mut m = Mesh::new(name);
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(side, 0.0, 0.0));
        m.nodes.push(Vector3::new(side, side, 0.0));
        m.nodes.push(Vector3::new(0.0, side, 0.0));
        Part::new(name, m)
    }

    #[test]
    fn auto_pack_places_two_parts_side_by_side() {
        let printer = Printer::new(
            (220.0, 220.0, 250.0),
            BedType::Heated,
            BedMaterial::Pei,
        );
        let parts = vec![cube_part("a", 50.0), cube_part("b", 50.0)];
        let placed = auto_pack(parts, &printer).unwrap();
        // Two parts shouldn't share an x position.
        assert_ne!(placed[0].bed_position[0], placed[1].bed_position[0]);
    }

    #[test]
    fn auto_pack_rejects_oversized_part() {
        let printer = Printer::new(
            (100.0, 100.0, 100.0),
            BedType::Unheated,
            BedMaterial::Glass,
        );
        let parts = vec![cube_part("big", 200.0)];
        let err = auto_pack(parts, &printer).unwrap_err();
        assert!(matches!(err, PrintBedError::PartTooLarge { .. }));
    }
}

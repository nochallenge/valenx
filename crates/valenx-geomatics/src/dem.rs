//! Digital Elevation Model (regular grid).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::GeomaticsError;

/// A regular-grid Digital Elevation Model.
///
/// Layout:
/// - `(origin_x, origin_y)` is the lower-left grid cell.
/// - `grid_size_m` is the spacing of one cell along both X and Y
///   (square cells in v1).
/// - `nx` × `ny` cells (so `(nx + 1) * (ny + 1)` corner samples or
///   `nx * ny` cell-centre samples — the convention here is corner-
///   sampled, so `data.len() == nx * ny` describes `nx-1`/`ny-1` cells).
/// - `data` is row-major, `data[y * nx + x]` reads the sample at
///   `(origin_x + x * grid_size_m, origin_y + y * grid_size_m)`.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Dem {
    /// Origin x (any planar coord — typically UTM easting in m).
    pub origin_x: f64,
    /// Origin y (UTM northing in m, etc).
    pub origin_y: f64,
    /// Cell side length.
    pub grid_size_m: f64,
    /// Number of samples along x.
    pub nx: usize,
    /// Number of samples along y.
    pub ny: usize,
    /// Row-major samples, length `nx * ny`.
    pub data: Vec<f32>,
}

impl Dem {
    /// Construct from a row-major `data` Vec. Sanity-checks the
    /// length against `nx * ny`.
    pub fn from_grid(
        origin_x: f64,
        origin_y: f64,
        grid_size_m: f64,
        nx: usize,
        ny: usize,
        data: Vec<f32>,
    ) -> Result<Self, GeomaticsError> {
        if grid_size_m <= 0.0 {
            return Err(GeomaticsError::BadParameter {
                name: "grid_size_m",
                reason: format!("must be > 0, got {grid_size_m}"),
            });
        }
        if nx == 0 || ny == 0 {
            return Err(GeomaticsError::BadParameter {
                name: "nx|ny",
                reason: format!("must be > 0, got {nx}x{ny}"),
            });
        }
        if data.len() != nx * ny {
            return Err(GeomaticsError::BadParameter {
                name: "data",
                reason: format!("expected {}, got {}", nx * ny, data.len()),
            });
        }
        Ok(Self {
            origin_x,
            origin_y,
            grid_size_m,
            nx,
            ny,
            data,
        })
    }

    /// Read the sample at `(ix, iy)`. Returns `None` for out-of-
    /// range indices.
    pub fn at(&self, ix: usize, iy: usize) -> Option<f32> {
        if ix >= self.nx || iy >= self.ny {
            return None;
        }
        Some(self.data[iy * self.nx + ix])
    }
}

/// Parse `x y z` whitespace-separated ASCII DEM data into a [`Dem`].
///
/// The function infers grid spacing from the first row's `dx` and
/// errors if subsequent rows don't match (within 1 % tolerance).
pub fn from_xyz_ascii(text: &str) -> Result<Dem, GeomaticsError> {
    let mut rows: Vec<(f64, f64, f64)> = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(GeomaticsError::Parse {
                line: i + 1,
                msg: format!("expected 3 columns, got {}", parts.len()),
            });
        }
        let x: f64 =
            parts[0]
                .parse()
                .map_err(|e: std::num::ParseFloatError| GeomaticsError::Parse {
                    line: i + 1,
                    msg: format!("bad x: {e}"),
                })?;
        let y: f64 =
            parts[1]
                .parse()
                .map_err(|e: std::num::ParseFloatError| GeomaticsError::Parse {
                    line: i + 1,
                    msg: format!("bad y: {e}"),
                })?;
        let z: f64 =
            parts[2]
                .parse()
                .map_err(|e: std::num::ParseFloatError| GeomaticsError::Parse {
                    line: i + 1,
                    msg: format!("bad z: {e}"),
                })?;
        rows.push((x, y, z));
    }
    if rows.is_empty() {
        return Err(GeomaticsError::IrregularGrid("no samples".into()));
    }
    // Sort by (y, x) so the row-major reshape is well-defined.
    rows.sort_by(|a, b| {
        a.1.partial_cmp(&b.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
    });
    // Find unique y's → ny.
    let mut ys: Vec<f64> = rows.iter().map(|r| r.1).collect();
    ys.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let ny = ys.len();
    let nx = rows.len() / ny;
    if nx * ny != rows.len() {
        return Err(GeomaticsError::IrregularGrid(format!(
            "{} samples don't fit a {}x{} grid",
            rows.len(),
            nx,
            ny
        )));
    }
    let origin_x = rows[0].0;
    let origin_y = rows[0].1;
    let grid_size_m = if nx >= 2 {
        rows[1].0 - rows[0].0
    } else if ny >= 2 {
        rows[nx].1 - rows[0].1
    } else {
        1.0
    };
    if grid_size_m <= 0.0 {
        return Err(GeomaticsError::IrregularGrid(format!(
            "non-positive grid spacing inferred: {grid_size_m}"
        )));
    }
    // `Dem` carries a single spacing for both axes (square cells, v1). A
    // regular but ANISOTROPIC grid (dy != dx) would be stored with the
    // X-derived spacing and then silently misplace every row in Y, so reject
    // it rather than return wrong coordinates.
    if nx >= 2 && ny >= 2 {
        let dy = rows[nx].1 - rows[0].1;
        if (dy - grid_size_m).abs() > 1e-6 * grid_size_m {
            return Err(GeomaticsError::IrregularGrid(format!(
                "non-square grid: x spacing {grid_size_m} != y spacing {dy} \
                 (only square cells are supported)"
            )));
        }
    }
    let mut data = Vec::with_capacity(rows.len());
    for r in &rows {
        data.push(r.2 as f32);
    }
    Dem::from_grid(origin_x, origin_y, grid_size_m, nx, ny, data)
}

/// Build a triangulated surface mesh from a [`Dem`]. Each grid cell
/// is split into two triangles; vertex Z carries the elevation.
pub fn to_mesh(dem: &Dem) -> Mesh {
    let mut out = Mesh::new("dem_surface");
    for iy in 0..dem.ny {
        for ix in 0..dem.nx {
            let x = dem.origin_x + (ix as f64) * dem.grid_size_m;
            let y = dem.origin_y + (iy as f64) * dem.grid_size_m;
            let z = dem.data[iy * dem.nx + ix] as f64;
            out.nodes.push(Vector3::new(x, y, z));
        }
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    for iy in 0..dem.ny.saturating_sub(1) {
        for ix in 0..dem.nx.saturating_sub(1) {
            let nx = dem.nx as u32;
            let a = (iy as u32) * nx + ix as u32;
            let b = a + 1;
            let c = a + nx;
            let d = c + 1;
            block.connectivity.extend_from_slice(&[a, b, d]);
            block.connectivity.extend_from_slice(&[a, d, c]);
        }
    }
    out.element_blocks.push(block);
    out.recompute_stats();
    out
}

/// Bilinear interpolation at planar `(x, y)`. Returns `0.0` when the
/// query is outside the grid.
pub fn sample(dem: &Dem, x: f64, y: f64) -> f32 {
    if dem.grid_size_m <= 0.0 || dem.nx == 0 || dem.ny == 0 {
        return 0.0;
    }
    let fx = (x - dem.origin_x) / dem.grid_size_m;
    let fy = (y - dem.origin_y) / dem.grid_size_m;
    if fx < 0.0 || fy < 0.0 {
        return 0.0;
    }
    let ix = fx.floor() as usize;
    let iy = fy.floor() as usize;
    if ix >= dem.nx - 1 || iy >= dem.ny - 1 {
        return 0.0;
    }
    let dx = fx - ix as f64;
    let dy = fy - iy as f64;
    let v00 = dem.data[iy * dem.nx + ix] as f64;
    let v10 = dem.data[iy * dem.nx + ix + 1] as f64;
    let v01 = dem.data[(iy + 1) * dem.nx + ix] as f64;
    let v11 = dem.data[(iy + 1) * dem.nx + ix + 1] as f64;
    let v0 = v00 * (1.0 - dx) + v10 * dx;
    let v1 = v01 * (1.0 - dx) + v11 * dx;
    (v0 * (1.0 - dy) + v1 * dy) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_grid_validates() {
        let r = Dem::from_grid(0.0, 0.0, 1.0, 2, 2, vec![1.0, 2.0, 3.0]);
        assert!(matches!(r, Err(GeomaticsError::BadParameter { .. })));
        let ok = Dem::from_grid(0.0, 0.0, 1.0, 2, 2, vec![1.0, 2.0, 3.0, 4.0]).unwrap();
        assert_eq!(ok.at(1, 1), Some(4.0));
    }

    #[test]
    fn parse_xyz_basic() {
        let text = "0 0 10\n1 0 20\n0 1 30\n1 1 40\n";
        let dem = from_xyz_ascii(text).unwrap();
        assert_eq!(dem.nx, 2);
        assert_eq!(dem.ny, 2);
        assert_eq!(dem.data, vec![10.0, 20.0, 30.0, 40.0]);
    }

    #[test]
    fn parse_xyz_bad_row_errors() {
        let text = "0 0\n";
        assert!(matches!(
            from_xyz_ascii(text),
            Err(GeomaticsError::Parse { .. })
        ));
    }

    #[test]
    fn parse_xyz_irregular_errors() {
        let text = "0 0 10\n1 0 20\n0 1 30\n"; // 3 samples → can't reshape
        assert!(matches!(
            from_xyz_ascii(text),
            Err(GeomaticsError::IrregularGrid(_))
        ));
    }

    #[test]
    fn parse_xyz_rejects_non_square_grid() {
        // A valid regular grid but with dx = 10, dy = 5. `Dem` is square-cells
        // only; accepting it would store grid_size_m = 10 and misplace every
        // row in Y. Must fail loud, not silently return wrong coordinates.
        let text = "0 0 10\n10 0 20\n0 5 30\n10 5 40\n";
        assert!(matches!(
            from_xyz_ascii(text),
            Err(GeomaticsError::IrregularGrid(_))
        ));
    }

    #[test]
    fn to_mesh_basic() {
        let dem = Dem::from_grid(0.0, 0.0, 1.0, 3, 3, vec![0.0; 9]).unwrap();
        let m = to_mesh(&dem);
        assert_eq!(m.nodes.len(), 9);
        // 2 × 2 cells * 2 tris = 8 triangles.
        assert_eq!(m.total_elements(), 8);
    }

    #[test]
    fn bilinear_sampling_interior() {
        let dem = Dem::from_grid(0.0, 0.0, 1.0, 2, 2, vec![0.0, 10.0, 20.0, 30.0]).unwrap();
        // Midpoint should be average of all 4 corners = 15.0.
        let v = sample(&dem, 0.5, 0.5);
        assert!((v - 15.0).abs() < 1e-4);
    }

    #[test]
    fn bilinear_sampling_outside_returns_zero() {
        let dem = Dem::from_grid(0.0, 0.0, 1.0, 2, 2, vec![0.0, 10.0, 20.0, 30.0]).unwrap();
        assert_eq!(sample(&dem, -1.0, 0.0), 0.0);
        assert_eq!(sample(&dem, 10.0, 0.0), 0.0);
    }
}

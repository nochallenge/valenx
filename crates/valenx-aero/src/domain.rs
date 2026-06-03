//! The virtual-wind-tunnel domain builder and the boundary-condition
//! set.
//!
//! [`WindTunnel`] auto-sizes a rectangular box domain around a body —
//! enough clearance upstream, downstream and to the sides that the
//! tunnel walls do not contaminate the flow over the body — places the
//! body inside it, voxelizes it, and builds the Cartesian grid. It is
//! the one call that turns *"this triangle mesh + this wind"* into a
//! ready-to-solve case.
//!
//! [`BoundaryConditions`] names the role of each of the six tunnel
//! faces: a uniform-velocity inlet on the upstream face, a convective
//! / zero-gradient outlet downstream, slip (symmetry) far-field on the
//! sides and top, and either a slip or a moving no-slip wall on the
//! floor (the latter models a road under a car).

use nalgebra::Vector3;

use crate::error::AeroError;
use crate::geometry::TriMesh;
use crate::grid::Grid3;
use crate::immersed::{voxelize, ImmersedBody};
use crate::wind::Wind;

/// The condition on one tunnel face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FaceBc {
    /// A uniform-velocity inlet — the free-stream enters here.
    Inlet,
    /// A convective / zero-gradient outlet — the flow leaves here, the
    /// pressure is pinned to the reference.
    Outlet,
    /// A symmetry / slip plane — the normal velocity and the normal
    /// gradient of the tangential velocity vanish. Models an
    /// undisturbed far-field with no spurious wall friction.
    Slip,
    /// A no-slip wall moving with the given tangential velocity. A
    /// stationary road is `MovingWall(0)`; a *moving-ground* road
    /// under a car carries the road at the free-stream speed so the
    /// relative wind at the floor is correct.
    MovingWall(f64),
}

/// The six-faced boundary specification of the rectangular tunnel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryConditions {
    /// Upstream face (`-x`).
    pub x_min: FaceBc,
    /// Downstream face (`+x`).
    pub x_max: FaceBc,
    /// Side face (`-y`).
    pub y_min: FaceBc,
    /// Side face (`+y`).
    pub y_max: FaceBc,
    /// Floor (`-z`).
    pub z_min: FaceBc,
    /// Ceiling (`+z`).
    pub z_max: FaceBc,
}

impl BoundaryConditions {
    /// The standard external-aero tunnel: inlet upstream, outlet
    /// downstream, slip on the two sides and the ceiling, slip on the
    /// floor (a free-air run with no ground).
    pub fn external_aero() -> BoundaryConditions {
        BoundaryConditions {
            x_min: FaceBc::Inlet,
            x_max: FaceBc::Outlet,
            y_min: FaceBc::Slip,
            y_max: FaceBc::Slip,
            z_min: FaceBc::Slip,
            z_max: FaceBc::Slip,
        }
    }

    /// The automotive tunnel: as [`external_aero`](Self::external_aero)
    /// but the floor is a no-slip wall moving at `road_speed` — the
    /// moving-ground condition that models a car on a road correctly
    /// (the road and the car move together, the air is the reference).
    pub fn automotive(road_speed: f64) -> BoundaryConditions {
        let mut bc = BoundaryConditions::external_aero();
        bc.z_min = FaceBc::MovingWall(road_speed);
        bc
    }
}

/// A fully-built virtual wind tunnel — a grid, a placed body
/// voxelization, the wind, the boundary set, and the reference
/// quantities derived from the body.
#[derive(Clone, Debug)]
pub struct WindTunnel {
    /// The Cartesian solver grid.
    pub grid: Grid3,
    /// The voxelized immersed body.
    pub body: ImmersedBody,
    /// The on-coming wind.
    pub wind: Wind,
    /// The six-face boundary specification.
    pub bc: BoundaryConditions,
    /// The reference (frontal) area for coefficient normalisation
    /// (m²).
    pub reference_area: f64,
    /// The characteristic body length along the wind (m) — used for
    /// the Reynolds number and the moment-coefficient normalisation.
    pub reference_length: f64,
}

/// How big to make the auto-sized tunnel box, in multiples of the
/// body's bounding-box size.
#[derive(Clone, Copy, Debug)]
pub struct TunnelSizing {
    /// Clearance upstream of the body (multiples of body length).
    pub upstream: f64,
    /// Clearance downstream of the body — the wake needs room
    /// (multiples of body length).
    pub downstream: f64,
    /// Lateral / vertical clearance on each side (multiples of the
    /// body's lateral size).
    pub lateral: f64,
    /// Target number of cells across the body's smallest dimension —
    /// the grid is sized to resolve the body at roughly this
    /// resolution.
    pub cells_across_body: usize,
    /// Hard cap on the total cell count so an over-ambitious request
    /// cannot allocate an unbounded grid.
    pub max_cells: usize,
}

impl Default for TunnelSizing {
    /// Conservative, broadly-valid defaults: 3 body-lengths upstream,
    /// 8 downstream (a generous wake), 3 lateral, ~16 cells across the
    /// body, capped at 2 million cells.
    fn default() -> Self {
        TunnelSizing {
            upstream: 3.0,
            downstream: 8.0,
            lateral: 3.0,
            cells_across_body: 16,
            max_cells: 2_000_000,
        }
    }
}

impl WindTunnel {
    /// Build a virtual wind tunnel around a body with the default
    /// sizing.
    ///
    /// `body` is the triangle-mesh geometry in any coordinate frame;
    /// it is translated so its bounding box sits in the placed tunnel.
    /// `wind` defines the on-coming flow. The reference area is the
    /// body's projected frontal area onto the plane normal to the
    /// wind.
    pub fn build(body: &TriMesh, wind: Wind) -> Result<WindTunnel, AeroError> {
        WindTunnel::build_with(body, wind, BoundaryConditions::external_aero(), TunnelSizing::default())
    }

    /// Build a virtual wind tunnel with an explicit boundary set and
    /// sizing policy.
    pub fn build_with(
        body: &TriMesh,
        wind: Wind,
        bc: BoundaryConditions,
        sizing: TunnelSizing,
    ) -> Result<WindTunnel, AeroError> {
        if body.is_empty() {
            return Err(AeroError::BadGeometry("body has no triangles".into()));
        }
        let bb = body
            .aabb()
            .ok_or_else(|| AeroError::BadGeometry("body bounding box is empty".into()))?;
        let extent = bb.extent();
        if !(extent.x.is_finite() && extent.y.is_finite() && extent.z.is_finite())
            || extent.x <= 0.0
            || extent.y <= 0.0
            || extent.z <= 0.0
        {
            return Err(AeroError::BadGeometry(format!(
                "body bounding box has zero / non-finite extent: {extent:?}"
            )));
        }

        // The wind blows nominally along +x; the tunnel's long axis is
        // x. Size the box from the body extent.
        let body_len = extent.x.max(extent.y).max(extent.z);
        let lx = extent.x + (sizing.upstream + sizing.downstream) * body_len;
        let ly = extent.y + 2.0 * sizing.lateral * extent.y.max(0.25 * body_len);
        let lz = extent.z + 2.0 * sizing.lateral * extent.z.max(0.25 * body_len);

        // Cell size: resolve the body's smallest dimension with the
        // requested cell count.
        let smallest = extent.x.min(extent.y).min(extent.z);
        let mut h = smallest / sizing.cells_across_body.max(2) as f64;
        // Enforce the total-cell cap by growing the cell size if the
        // naive grid would be too large.
        loop {
            let nx = (lx / h).ceil() as usize;
            let ny = (ly / h).ceil() as usize;
            let nz = (lz / h).ceil() as usize;
            if nx.saturating_mul(ny).saturating_mul(nz) <= sizing.max_cells.max(8) {
                break;
            }
            h *= 1.26; // grow ~2× the cell volume each step
            if !h.is_finite() || h > 1e12 {
                return Err(AeroError::DomainTooSmall(
                    "could not fit the tunnel under the cell cap".into(),
                ));
            }
        }
        // Round the cell counts up to even numbers so the geometric-
        // multigrid solver can coarsen the grid.
        let round_even = |n: usize| -> usize {
            let n = n.max(4);
            if n % 2 == 0 {
                n
            } else {
                n + 1
            }
        };
        let nx = round_even((lx / h).ceil() as usize);
        let ny = round_even((ly / h).ceil() as usize);
        let nz = round_even((lz / h).ceil() as usize);
        // Final domain lengths follow from the rounded cell counts.
        let lx = nx as f64 * h;
        let ly = ny as f64 * h;
        let lz = nz as f64 * h;

        // Place the body: upstream clearance puts the body min-x at
        // `upstream·body_len` from the inlet; centre it laterally.
        let x0 = bb.min.x - sizing.upstream * body_len;
        let y0 = bb.centre().y - 0.5 * ly;
        let z0 = bb.centre().z - 0.5 * lz;
        let grid = Grid3::new(nx, ny, nz, lx, ly, lz, x0, y0, z0);

        // Voxelize the body in place (no translation — the grid was
        // anchored around the body's existing coordinates).
        let voxel = voxelize(&grid, body);
        let (_, cut, solid) = voxel.tag_counts();
        if cut + solid == 0 {
            return Err(AeroError::DomainTooSmall(
                "body voxelized to zero cells — grid too coarse for this body".into(),
            ));
        }

        let dir = wind.direction();
        let reference_area = body.frontal_area(dir).max(1e-9);
        let reference_length = extent.x.max(1e-9);

        Ok(WindTunnel {
            grid,
            body: voxel,
            wind,
            bc,
            reference_area,
            reference_length,
        })
    }

    /// The Reynolds number of this case, on the reference length.
    pub fn reynolds_number(&self) -> f64 {
        self.wind.reynolds_number(self.reference_length)
    }

    /// The free-stream dynamic pressure of this case (Pa).
    pub fn dynamic_pressure(&self) -> f64 {
        self.wind.dynamic_pressure()
    }

    /// The world-space free-stream velocity vector.
    pub fn free_stream(&self) -> Vector3<f64> {
        self.wind.velocity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{box_body, sphere_body};

    #[test]
    fn external_aero_bc_has_inlet_outlet_and_slip_walls() {
        let bc = BoundaryConditions::external_aero();
        assert_eq!(bc.x_min, FaceBc::Inlet);
        assert_eq!(bc.x_max, FaceBc::Outlet);
        assert_eq!(bc.y_min, FaceBc::Slip);
        assert_eq!(bc.z_min, FaceBc::Slip);
    }

    #[test]
    fn automotive_bc_has_a_moving_floor() {
        let bc = BoundaryConditions::automotive(30.0);
        assert_eq!(bc.z_min, FaceBc::MovingWall(30.0));
        // The other faces are unchanged.
        assert_eq!(bc.x_min, FaceBc::Inlet);
        assert_eq!(bc.z_max, FaceBc::Slip);
    }

    #[test]
    fn build_auto_sizes_a_box_around_the_body() {
        // A 4×2×1.5 box body — a car-ish bluff body.
        let body = box_body(Vector3::zeros(), Vector3::new(4.0, 2.0, 1.5));
        let wind = Wind::straight(30.0).unwrap();
        let tunnel = WindTunnel::build(&body, wind).unwrap();
        // The tunnel must be much longer than the body (upstream +
        // downstream clearance).
        assert!(tunnel.grid.lx > 4.0 * 4.0, "tunnel too short");
        // The body must have voxelized to a non-trivial solid region.
        let (_, cut, solid) = tunnel.body.tag_counts();
        assert!(solid > 0 && cut > 0);
        // The cell count is under the cap.
        assert!(tunnel.grid.cell_count() <= TunnelSizing::default().max_cells);
        // Grid cell counts are even (multigrid-coarsenable).
        assert_eq!(tunnel.grid.nx % 2, 0);
        assert_eq!(tunnel.grid.ny % 2, 0);
        assert_eq!(tunnel.grid.nz % 2, 0);
    }

    #[test]
    fn build_rejects_an_empty_body() {
        let wind = Wind::straight(20.0).unwrap();
        assert!(WindTunnel::build(&TriMesh::new(), wind).is_err());
    }

    #[test]
    fn reference_area_is_the_frontal_silhouette() {
        // A 2×3×4 box, wind +x: frontal area = 3×4 = 12.
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 3.0, 4.0));
        let tunnel = WindTunnel::build(&body, Wind::straight(25.0).unwrap()).unwrap();
        assert!(
            (tunnel.reference_area - 12.0).abs() < 0.5,
            "reference area {} should be ~12",
            tunnel.reference_area
        );
    }

    #[test]
    fn reynolds_number_is_in_the_expected_regime() {
        let body = sphere_body(Vector3::zeros(), 0.5, 24, 48);
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        // Sphere diameter 1 m, U = 20 m/s → Re ~ 1.3e6.
        let re = tunnel.reynolds_number();
        assert!(re > 1.0e5 && re < 1.0e7, "sphere Re {re} out of range");
    }

    #[test]
    fn cell_cap_grows_the_cell_size_for_a_demanding_request() {
        // A demanding sizing — many cells across a small body — must
        // still respect the cap by coarsening.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let sizing = TunnelSizing {
            cells_across_body: 200,
            max_cells: 100_000,
            ..TunnelSizing::default()
        };
        let tunnel = WindTunnel::build_with(
            &body,
            Wind::straight(10.0).unwrap(),
            BoundaryConditions::external_aero(),
            sizing,
        )
        .unwrap();
        assert!(tunnel.grid.cell_count() <= 100_000);
    }
}

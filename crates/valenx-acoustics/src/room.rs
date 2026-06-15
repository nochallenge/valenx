//! Rectangular-room standing-wave (eigenmode) frequencies.
//!
//! ## Model
//!
//! A rigid-walled rectangular room of interior dimensions `Lx * Ly * Lz`
//! supports a discrete set of acoustic standing waves. The Rayleigh modal
//! frequencies are
//!
//! ```text
//! f(nx, ny, nz) = (c / 2) * sqrt( (nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2 )
//! ```
//!
//! where `c` is the speed of sound and `nx, ny, nz` are non-negative
//! integer mode orders (not all zero). Modes are conventionally classed
//! by how many indices are non-zero:
//!
//! - **Axial** — one non-zero index (between one opposite wall pair); the
//!   strongest and lowest. `f(1,0,0) = c / (2 Lx)`.
//! - **Tangential** — two non-zero indices.
//! - **Oblique** — all three non-zero.
//!
//! ## Honest scope
//!
//! This is the textbook rigid-wall eigenfrequency formula. It models only
//! the modal *frequencies* of an idealised empty rectangular box with
//! perfectly reflecting walls. It says nothing about modal amplitude,
//! damping / Q, wall admittance, absorption, the Schroeder transition, or
//! non-rectangular geometry. Research/educational grade — a first-order
//! room-mode estimator, not a boundary-element room solver.

use serde::{Deserialize, Serialize};

use crate::error::{AcousticsError, Result};

/// Interior dimensions of a rectangular room, in metres.
///
/// Construct with [`RoomDimensions::new`], which validates that every edge
/// is finite and strictly positive.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RoomDimensions {
    /// Length along x, in metres (`Lx`).
    pub length_x: f64,
    /// Length along y, in metres (`Ly`).
    pub length_y: f64,
    /// Length along z, in metres (`Lz`).
    pub length_z: f64,
}

/// Classification of a room mode by how many of its indices are non-zero.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum ModeKind {
    /// Exactly one non-zero index — between one pair of opposite walls.
    Axial,
    /// Exactly two non-zero indices.
    Tangential,
    /// All three indices non-zero.
    Oblique,
}

/// A single resolved room mode: its integer orders, class and frequency.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RoomMode {
    /// Mode order along x.
    pub nx: u32,
    /// Mode order along y.
    pub ny: u32,
    /// Mode order along z.
    pub nz: u32,
    /// Axial / tangential / oblique classification.
    pub kind: ModeKind,
    /// Modal frequency, in hertz.
    pub frequency_hz: f64,
}

impl RoomDimensions {
    /// Build a validated set of room dimensions (metres).
    ///
    /// # Errors
    ///
    /// Returns [`AcousticsError::InvalidDimension`] for the first edge
    /// that is not finite and strictly positive.
    pub fn new(length_x: f64, length_y: f64, length_z: f64) -> Result<Self> {
        check_dimension("length_x", length_x)?;
        check_dimension("length_y", length_y)?;
        check_dimension("length_z", length_z)?;
        Ok(Self {
            length_x,
            length_y,
            length_z,
        })
    }

    /// Frequency (hertz) of the mode with integer orders `(nx, ny, nz)` at
    /// a speed of sound `speed_of_sound` (m/s):
    ///
    /// `f = (c/2) * sqrt((nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2)`.
    ///
    /// # Errors
    ///
    /// - [`AcousticsError::InvalidSpeedOfSound`] if `speed_of_sound` is not
    ///   finite and strictly positive.
    /// - [`AcousticsError::TrivialMode`] for the `(0, 0, 0)` mode, which is
    ///   the zero-frequency DC mode.
    pub fn mode_frequency(&self, nx: u32, ny: u32, nz: u32, speed_of_sound: f64) -> Result<f64> {
        if !speed_of_sound.is_finite() || speed_of_sound <= 0.0 {
            return Err(AcousticsError::InvalidSpeedOfSound {
                name: "speed_of_sound",
                value: speed_of_sound,
            });
        }
        if nx == 0 && ny == 0 && nz == 0 {
            return Err(AcousticsError::TrivialMode);
        }
        let tx = nx as f64 / self.length_x;
        let ty = ny as f64 / self.length_y;
        let tz = nz as f64 / self.length_z;
        Ok(0.5 * speed_of_sound * (tx * tx + ty * ty + tz * tz).sqrt())
    }

    /// Resolve a single [`RoomMode`] (orders + class + frequency) for
    /// `(nx, ny, nz)` at the given speed of sound.
    ///
    /// # Errors
    ///
    /// Same as [`mode_frequency`](Self::mode_frequency).
    pub fn mode(&self, nx: u32, ny: u32, nz: u32, speed_of_sound: f64) -> Result<RoomMode> {
        let frequency_hz = self.mode_frequency(nx, ny, nz, speed_of_sound)?;
        Ok(RoomMode {
            nx,
            ny,
            nz,
            kind: classify(nx, ny, nz),
            frequency_hz,
        })
    }

    /// Enumerate every mode with each index in `0..=max_order` (excluding
    /// the trivial `(0,0,0)` mode), sorted by ascending frequency.
    ///
    /// Useful for listing the low-frequency modal "stack" of a room.
    ///
    /// # Errors
    ///
    /// Returns [`AcousticsError::InvalidSpeedOfSound`] if `speed_of_sound`
    /// is not finite and strictly positive.
    pub fn modes_up_to(&self, max_order: u32, speed_of_sound: f64) -> Result<Vec<RoomMode>> {
        if !speed_of_sound.is_finite() || speed_of_sound <= 0.0 {
            return Err(AcousticsError::InvalidSpeedOfSound {
                name: "speed_of_sound",
                value: speed_of_sound,
            });
        }
        let mut modes = Vec::new();
        for nx in 0..=max_order {
            for ny in 0..=max_order {
                for nz in 0..=max_order {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }
                    // Speed already validated; this cannot fail.
                    let m = self.mode(nx, ny, nz, speed_of_sound)?;
                    modes.push(m);
                }
            }
        }
        modes.sort_by(|a, b| {
            a.frequency_hz
                .partial_cmp(&b.frequency_hz)
                .unwrap_or(core::cmp::Ordering::Equal)
        });
        Ok(modes)
    }
}

/// Classify a mode by its number of non-zero indices. Panics-free; the
/// all-zero case is reported as [`ModeKind::Axial`] but is never produced
/// by the public API, which rejects `(0,0,0)`.
fn classify(nx: u32, ny: u32, nz: u32) -> ModeKind {
    let nonzero = (nx > 0) as u8 + (ny > 0) as u8 + (nz > 0) as u8;
    match nonzero {
        3 => ModeKind::Oblique,
        2 => ModeKind::Tangential,
        _ => ModeKind::Axial,
    }
}

fn check_dimension(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidDimension { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;
    const C: f64 = 343.0;

    /// The (1,0,0) axial mode is exactly c/(2*Lx).
    #[test]
    fn first_axial_mode_is_c_over_two_lx() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let f = room.mode_frequency(1, 0, 0, C).unwrap();
        let expected = C / (2.0 * 5.0);
        assert!((f - expected).abs() < EPS, "got {f}, expected {expected}");
    }

    /// (0,1,0) and (0,0,1) axial modes pick out Ly and Lz respectively.
    #[test]
    fn axial_modes_pick_out_each_dimension() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let fy = room.mode_frequency(0, 1, 0, C).unwrap();
        let fz = room.mode_frequency(0, 0, 1, C).unwrap();
        assert!((fy - C / (2.0 * 4.0)).abs() < EPS, "fy = {fy}");
        assert!((fz - C / (2.0 * 3.0)).abs() < EPS, "fz = {fz}");
    }

    /// The second-order axial mode (2,0,0) is exactly twice the first.
    #[test]
    fn axial_overtone_is_harmonic() {
        let room = RoomDimensions::new(6.0, 4.0, 2.5).unwrap();
        let f1 = room.mode_frequency(1, 0, 0, C).unwrap();
        let f2 = room.mode_frequency(2, 0, 0, C).unwrap();
        assert!((f2 - 2.0 * f1).abs() < EPS, "{f2} vs {}", 2.0 * f1);
    }

    /// Tangential (1,1,0) matches the Pythagorean sum of two axials.
    #[test]
    fn tangential_mode_matches_formula() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let f = room.mode_frequency(1, 1, 0, C).unwrap();
        let expected = 0.5 * C * ((1.0f64 / 5.0).powi(2) + (1.0f64 / 4.0).powi(2)).sqrt();
        assert!((f - expected).abs() < EPS, "got {f}, expected {expected}");
    }

    /// A cube's (1,1,1) oblique mode is sqrt(3) times its (1,0,0) axial.
    #[test]
    fn cube_oblique_is_sqrt3_times_axial() {
        let room = RoomDimensions::new(3.0, 3.0, 3.0).unwrap();
        let axial = room.mode_frequency(1, 0, 0, C).unwrap();
        let oblique = room.mode_frequency(1, 1, 1, C).unwrap();
        assert!(
            (oblique - axial * 3.0_f64.sqrt()).abs() < EPS,
            "{oblique} vs {}",
            axial * 3.0_f64.sqrt()
        );
    }

    /// Classification by non-zero index count.
    #[test]
    fn mode_classification() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        assert_eq!(room.mode(1, 0, 0, C).unwrap().kind, ModeKind::Axial);
        assert_eq!(room.mode(1, 1, 0, C).unwrap().kind, ModeKind::Tangential);
        assert_eq!(room.mode(1, 1, 1, C).unwrap().kind, ModeKind::Oblique);
    }

    /// Enumerated modes are sorted ascending and exclude (0,0,0); the
    /// lowest mode of this room is the (1,0,0) axial along the longest
    /// wall.
    #[test]
    fn enumeration_is_sorted_and_excludes_dc() {
        let room = RoomDimensions::new(7.0, 5.0, 3.0).unwrap();
        let modes = room.modes_up_to(2, C).unwrap();
        // 3*3*3 - 1 = 26 modes for max_order = 2.
        assert_eq!(modes.len(), 26);
        for pair in modes.windows(2) {
            assert!(
                pair[0].frequency_hz <= pair[1].frequency_hz,
                "not sorted: {} then {}",
                pair[0].frequency_hz,
                pair[1].frequency_hz
            );
        }
        // Longest wall (7 m) gives the lowest mode.
        let first = modes[0];
        assert_eq!((first.nx, first.ny, first.nz), (1, 0, 0));
        assert!((first.frequency_hz - C / (2.0 * 7.0)).abs() < EPS);
    }

    /// The (0,0,0) DC mode is rejected.
    #[test]
    fn trivial_mode_rejected() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let err = room.mode_frequency(0, 0, 0, C).unwrap_err();
        assert_eq!(err.code(), "acoustics.trivial_mode");
    }

    /// Bad dimensions and bad speeds are rejected.
    #[test]
    fn invalid_inputs_rejected() {
        assert!(RoomDimensions::new(0.0, 4.0, 3.0).is_err());
        assert!(RoomDimensions::new(5.0, -1.0, 3.0).is_err());
        assert!(RoomDimensions::new(5.0, 4.0, f64::NAN).is_err());

        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        assert!(room.mode_frequency(1, 0, 0, 0.0).is_err());
        assert!(room.modes_up_to(2, -1.0).is_err());
    }

    /// `RoomDimensions` and `RoomMode` round-trip through JSON.
    #[test]
    fn serde_round_trip() {
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let json = serde_json::to_string(&room).unwrap();
        let back: RoomDimensions = serde_json::from_str(&json).unwrap();
        assert_eq!(room, back);

        let mode = room.mode(1, 1, 0, C).unwrap();
        let json = serde_json::to_string(&mode).unwrap();
        let back: RoomMode = serde_json::from_str(&json).unwrap();
        assert_eq!(mode, back);
    }
}

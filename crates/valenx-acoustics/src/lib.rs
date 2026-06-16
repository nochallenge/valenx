//! # valenx-acoustics
//!
//! Closed-form acoustics calculators: sound-pressure level and decibel
//! arithmetic, the temperature-dependent speed of sound, the classical
//! Doppler frequency shift, and rectangular-room standing-wave modes.
//!
//! ## What
//!
//! A small, dependency-light toolbox of the everyday textbook acoustics
//! formulas, each behind a validated function that rejects non-physical
//! inputs rather than silently returning `NaN`:
//!
//! - **Sound-pressure level** ([`mod@spl`]) — [`spl::spl`] gives
//!   `L = 20*log10(p/p_ref)` against the standard 20 micropascal air
//!   reference, with the exact inverse [`spl::pressure_from_spl`] and the
//!   incoherent decibel sum [`spl::combine_incoherent_levels`]
//!   (`10*log10(sum 10^(Li/10))`).
//! - **Speed of sound & Doppler** ([`doppler`]) —
//!   [`doppler::speed_of_sound`] gives `c = 331.3*sqrt(1+T/273.15)` and
//!   [`doppler::doppler_shift`] gives the classical
//!   `f_obs = f*(c+vo)/(c-vs)`.
//! - **Room modes** ([`room`]) — [`room::RoomDimensions`] computes the
//!   rigid-wall eigenfrequencies
//!   `f = (c/2)*sqrt((nx/Lx)^2+(ny/Ly)^2+(nz/Lz)^2)`, classifies each
//!   mode as axial / tangential / oblique, and enumerates the low-end
//!   modal stack.
//! - **Reverberation time** ([`reverberation`]) — the statistical
//!   diffuse-field `RT60`: [`reverberation::sabine_reverberation_time`]
//!   (`24 ln10 * V / (c*A)`, the textbook `~0.161 V/A`),
//!   [`reverberation::eyring_reverberation_time`] for absorptive rooms,
//!   and the [`reverberation::total_absorption`] helper `A = sum S_i a_i`.
//!
//! ## Model
//!
//! Each value is the standard closed-form definition:
//!
//! - SPL is the base-10 logarithm of an RMS pressure ratio scaled by 20;
//!   levels of incoherent sources combine on a mean-square (power) basis.
//! - The speed of sound is the first-order dry-air relation in degrees
//!   Celsius, with the `331.3` m/s value at 0 degrees Celsius.
//! - The Doppler shift is the 1-D still-medium form, with
//!   observer-velocity positive toward the source and source-velocity
//!   positive toward the observer.
//! - The room modes are the Rayleigh eigenfrequencies of an idealised
//!   rigid-walled rectangular box.
//!
//! These reproduce the standard rules of thumb exactly: the reference
//! pressure is `0` dB, doubling pressure adds `≈ 6` dB, two equal
//! incoherent sources add `≈ 3` dB, an approaching source raises the
//! pitch above the emitted frequency and a receding one lowers it,
//! `c(20 degrees C) ≈ 343` m/s, and the first axial room mode is exactly
//! `c/(2*Lx)`.
//!
//! ## Honest scope
//!
//! Research/educational grade. Every routine is a textbook closed-form,
//! well-established model — not a clinical / medical instrument and not a
//! production engineering acoustics tool. There is no frequency weighting
//! (A / C), no octave-band or fractional-octave machinery, no coherent
//! interference term, no atmospheric absorption or humidity correction in
//! the speed of sound, no transverse (angle-of-arrival) Doppler term, and
//! no modal damping / absorption / wall-admittance in the room solver.
//! Each module documents its own simplifications inline. Use it to learn,
//! prototype and sanity-check — not to certify a measurement or a design.
//!
//! ```
//! use valenx_acoustics::{combine_incoherent_levels, doppler_shift, speed_of_sound, spl};
//!
//! let c = speed_of_sound(20.0).unwrap();           // ~343 m/s
//! let level = spl(1.0).unwrap();                   // ~94 dB at 1 Pa
//! let combined = combine_incoherent_levels(&[60.0, 60.0]); // ~63 dB
//! let pitch = doppler_shift(440.0, 0.0, 30.0, c).unwrap(); // > 440 Hz
//! assert!(level > 93.0 && level < 95.0);
//! assert!(combined > 62.9 && combined < 63.1);
//! assert!(pitch > 440.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod doppler;
pub mod error;
pub mod reverberation;
pub mod room;
pub mod spl;

// --- Convenience re-exports of the most-used items --------------------

pub use error::{AcousticsError, ErrorCategory, Result};

pub use spl::{
    combine_incoherent_levels, incoherent_sum_excess, pressure_from_spl, pressure_from_spl_ref,
    spl, spl_ref, DOUBLE_POWER_DB, DOUBLE_PRESSURE_DB, P_REF,
};

pub use doppler::{doppler_shift, speed_of_sound, C0};

pub use room::{ModeKind, RoomDimensions, RoomMode};

pub use reverberation::{eyring_reverberation_time, sabine_reverberation_time, total_absorption};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// End-to-end: derive the speed of sound at room temperature, build a
    /// room, and confirm its first axial mode uses exactly that speed.
    #[test]
    fn speed_of_sound_drives_room_modes() {
        let c = speed_of_sound(20.0).unwrap();
        let room = RoomDimensions::new(5.0, 4.0, 3.0).unwrap();
        let f = room.mode_frequency(1, 0, 0, c).unwrap();
        assert!((f - c / (2.0 * 5.0)).abs() < EPS, "got {f}");
    }

    /// End-to-end: a pressure SPL round-trips and combining a level with
    /// itself adds the expected ~3 dB.
    #[test]
    fn spl_pipeline_round_trip() {
        let p = pressure_from_spl(85.0);
        let level = spl(p).unwrap();
        assert!((level - 85.0).abs() < EPS, "got {level}");

        let doubled = combine_incoherent_levels(&[level, level]);
        assert!(
            (doubled - (85.0 + DOUBLE_POWER_DB)).abs() < EPS,
            "got {doubled}"
        );
    }

    /// The crate-level re-exports point at the same constants the modules
    /// expose.
    #[test]
    fn reexports_are_consistent() {
        assert!((P_REF - spl::P_REF).abs() < f64::EPSILON);
        assert!((C0 - doppler::C0).abs() < f64::EPSILON);
        assert!((DOUBLE_PRESSURE_DB - 20.0 * 2.0_f64.log10()).abs() < 1e-12);
    }
}

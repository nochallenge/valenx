//! OS mouse injection from a validated pointer event (compiled only with the
//! `live-capture` feature).
//!
//! Maps a phone [`InputEvent`] — already validated to carry finite, in-range
//! `[0, 1]` coordinates by [`crate::parse_input`] — onto an absolute mouse
//! position and a left-button action via the cross-platform `enigo` crate:
//!
//! | kind   | action                                  |
//! |--------|-----------------------------------------|
//! | `down` | move to point, press left              |
//! | `move` | move to point (drag if already pressed) |
//! | `up`   | release left                            |
//! | `tap`  | move to point, click left              |
//!
//! Normalized coordinates are scaled against the OS main-display size. When the
//! captured target is a sub-window rather than the whole screen this is an
//! approximation (the coordinate space differs); it is good enough for the
//! companion-screen use case and is documented as such. Mapping the touch to
//! the exact captured window rectangle is a possible refinement.

use crate::{Config, InputEvent, InputKind};
use enigo::{Button, Coordinate, Direction, Enigo, Mouse, Settings};

/// Relay one validated [`InputEvent`] as OS mouse input.
///
/// Best-effort: any `enigo` failure is returned as an error string for the
/// caller to swallow, so a transient injection problem never crashes the
/// server or leaks details to the phone.
pub fn dispatch(ev: InputEvent, cfg: &Config) -> Result<(), String> {
    let _ = cfg; // reserved for future window-rectangle-aware mapping
    let mut enigo = Enigo::new(&Settings::default()).map_err(|e| format!("enigo init: {e}"))?;

    let (w, h) = enigo
        .main_display()
        .map_err(|e| format!("main display size: {e}"))?;
    // `ev.nx`/`ev.ny` are guaranteed finite and within [0, 1] by validation.
    let x = (ev.nx * f64::from(w - 1).max(0.0)).round() as i32;
    let y = (ev.ny * f64::from(h - 1).max(0.0)).round() as i32;

    match ev.kind {
        InputKind::Down => {
            move_to(&mut enigo, x, y)?;
            button(&mut enigo, Direction::Press)?;
        }
        InputKind::Move => {
            move_to(&mut enigo, x, y)?;
        }
        InputKind::Up => {
            button(&mut enigo, Direction::Release)?;
        }
        InputKind::Tap => {
            move_to(&mut enigo, x, y)?;
            button(&mut enigo, Direction::Click)?;
        }
    }
    Ok(())
}

fn move_to(enigo: &mut Enigo, x: i32, y: i32) -> Result<(), String> {
    enigo
        .move_mouse(x, y, Coordinate::Abs)
        .map_err(|e| format!("move_mouse: {e}"))
}

fn button(enigo: &mut Enigo, dir: Direction) -> Result<(), String> {
    enigo
        .button(Button::Left, dir)
        .map_err(|e| format!("button: {e}"))
}

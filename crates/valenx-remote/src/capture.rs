//! Screen capture → JPEG (compiled only with the `live-capture` feature).
//!
//! Picks the on-screen window whose title contains [`Config::title`] (case
//! insensitive) and encodes its current contents as a JPEG at
//! [`Config::quality`]. If no matching, non-minimized window is found it falls
//! back to capturing the primary monitor, so the phone still shows *something*
//! useful rather than an error.
//!
//! Uses the cross-platform `xcap` crate (X11 / Wayland / Windows / macOS), so
//! there is nothing Windows-specific here; the whole module simply does not
//! exist in the default (std-only) build.

use crate::Config;
use image::codecs::jpeg::JpegEncoder;
use image::{ExtendedColorType, RgbaImage};
use xcap::{Monitor, Window};

/// Capture the target window (or the primary monitor) as JPEG bytes.
///
/// Returns a human-readable error string on any failure; callers turn that
/// into a `501`/`500` response rather than panicking.
pub fn capture_jpeg(cfg: &Config) -> Result<Vec<u8>, String> {
    let rgba = capture_target(cfg)?;
    encode_jpeg(&rgba, cfg.quality)
}

/// Grab the best-matching window's pixels, falling back to the primary monitor.
fn capture_target(cfg: &Config) -> Result<RgbaImage, String> {
    let needle = cfg.title.to_lowercase();

    // Prefer a visible window whose title matches the configured substring.
    if let Ok(windows) = Window::all() {
        for w in windows {
            let minimized = w.is_minimized().unwrap_or(true);
            if minimized {
                continue;
            }
            let title = w.title().unwrap_or_default().to_lowercase();
            if !needle.is_empty() && title.contains(&needle) {
                if let Ok(img) = w.capture_image() {
                    return Ok(img);
                }
            }
        }
    }

    // Fallback: the primary monitor (or the first one we can find).
    let monitors = Monitor::all().map_err(|e| format!("enumerate monitors: {e}"))?;
    let primary = monitors
        .iter()
        .find(|m| m.is_primary().unwrap_or(false))
        .or_else(|| monitors.first())
        .ok_or_else(|| "no monitors found".to_string())?;
    primary
        .capture_image()
        .map_err(|e| format!("capture monitor: {e}"))
}

/// Encode an RGBA frame as JPEG at the given quality (1..=100).
fn encode_jpeg(rgba: &RgbaImage, quality: u8) -> Result<Vec<u8>, String> {
    let (w, h) = rgba.dimensions();
    // Clamp the dimensions to the JPEG-valid range (the original binary did the
    // same 1..=65535 clamp); a zero- or over-sized frame is a hard error.
    if w == 0 || h == 0 || w > 65535 || h > 65535 {
        return Err(format!("frame size out of range: {w}x{h}"));
    }

    // JPEG has no alpha channel, so drop it.
    let rgb = image::DynamicImage::ImageRgba8(rgba.clone()).to_rgb8();
    let mut out = Vec::new();
    let mut encoder = JpegEncoder::new_with_quality(&mut out, quality.clamp(1, 100));
    encoder
        .encode(rgb.as_raw(), w, h, ExtendedColorType::Rgb8)
        .map_err(|e| format!("jpeg encode: {e}"))?;
    Ok(out)
}

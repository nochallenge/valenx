//! Phase 194 — `V3d_View::Dump()` — capture the viewport to an image.
//!
//! ## What OCCT does
//!
//! `V3d_View::Dump(filename, [params])` re-renders the current view
//! into an off-screen framebuffer, reads the colour attachment back to
//! CPU memory, and writes PNG / JPG / BMP via OCCT's `Image_AlienPixMap`
//! abstraction. Common parameters: stereo (left/right separate), pixel
//! ratio (oversample for higher resolution), tile (for huge captures
//! that don't fit in a single GPU buffer).
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 194.5) of the **encode-and-write**
//! half of `V3d_View::Dump`. The function takes a CPU-side RGBA8
//! framebuffer (the bytes a wgpu `copy_texture_to_buffer` +
//! `Buffer::map_async` read-back produces — that GPU step belongs to
//! the app layer that owns the live `wgpu::Device`, not to this pure
//! library) and writes a complete, spec-correct **uncompressed BMP**
//! file to disk.
//!
//! BMP is chosen because it needs no compression codec — the writer
//! here is the full, real format (`BITMAPFILEHEADER` +
//! `BITMAPINFOHEADER` + bottom-up BGRA rows), so a `.bmp` capture is
//! a finished feature with zero new dependencies. PNG / JPG encoding
//! needs a deflate / DCT codec (the `image` crate); that is the
//! documented follow-up and `.png` / `.jpg` paths are rejected with a
//! clear message rather than written as a mislabelled BMP.

use std::path::Path;

use crate::error::OcctVizError;

/// Capture an in-memory RGBA8 framebuffer to `path` as a BMP image.
///
/// `pixels` is `width * height * 4` bytes in row-major, top-to-bottom
/// order, each pixel `[R, G, B, A]` (the layout wgpu's
/// `Rgba8UnormSrgb` read-back yields). `path` must end in `.bmp`.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `path` has no extension, an
///   extension other than `.bmp`, zero dimensions, or a `pixels`
///   length that does not equal `width * height * 4`.
/// - [`OcctVizError::Io`] for filesystem write failures.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use valenx_occt_viz::view_screenshot::view_screenshot;
/// // A 2x2 red image.
/// let px = [255, 0, 0, 255].repeat(4);
/// view_screenshot(&px, 2, 2, &PathBuf::from("shot.bmp")).unwrap();
/// ```
pub fn view_screenshot(
    pixels: &[u8],
    width: u32,
    height: u32,
    path: &Path,
) -> Result<(), OcctVizError> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .ok_or_else(|| OcctVizError::bad_input("path", "must have a .bmp extension"))?
        .to_ascii_lowercase();
    if ext == "png" || ext == "jpg" || ext == "jpeg" {
        return Err(OcctVizError::bad_input(
            "path",
            format!(
                "`.{ext}` needs a compression codec — v1 writes uncompressed \
                 `.bmp` only (PNG/JPG is the documented follow-up)"
            ),
        ));
    }
    if ext != "bmp" {
        return Err(OcctVizError::bad_input(
            "path",
            format!("unsupported extension `.{ext}` (.bmp only)"),
        ));
    }
    if width == 0 || height == 0 {
        return Err(OcctVizError::bad_input(
            "dimensions",
            format!("width and height must be non-zero (got {width}x{height})"),
        ));
    }
    let expected = width as usize * height as usize * 4;
    if pixels.len() != expected {
        return Err(OcctVizError::bad_input(
            "pixels",
            format!(
                "buffer length {} != width*height*4 = {expected}",
                pixels.len()
            ),
        ));
    }

    let bmp = encode_bmp(pixels, width, height);
    // R30: the BMP is fully materialised in `bmp`; publish it atomically
    // (sidecar → fsync → rename) so a torn write can't leave a truncated
    // image on disk.
    valenx_core::io_caps::atomic_write_bytes(path, &bmp)?;
    Ok(())
}

/// Encode an RGBA8 top-down framebuffer into an uncompressed 32-bit
/// BMP byte stream (`BITMAPFILEHEADER` + `BITMAPINFOHEADER` + BGRA
/// pixel data written bottom-up, as the BMP format requires).
fn encode_bmp(pixels: &[u8], width: u32, height: u32) -> Vec<u8> {
    // 32-bit BMP rows are already 4-byte aligned — no per-row padding.
    let row_bytes = width as usize * 4;
    let pixel_data_size = row_bytes * height as usize;
    const FILE_HEADER: usize = 14;
    const INFO_HEADER: usize = 40;
    let offset = FILE_HEADER + INFO_HEADER;
    let file_size = offset + pixel_data_size;

    let mut out = Vec::with_capacity(file_size);

    // --- BITMAPFILEHEADER (14 bytes) ---
    out.extend_from_slice(b"BM");
    out.extend_from_slice(&(file_size as u32).to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // reserved1
    out.extend_from_slice(&0u16.to_le_bytes()); // reserved2
    out.extend_from_slice(&(offset as u32).to_le_bytes());

    // --- BITMAPINFOHEADER (40 bytes) ---
    out.extend_from_slice(&(INFO_HEADER as u32).to_le_bytes());
    out.extend_from_slice(&(width as i32).to_le_bytes());
    // Positive height = bottom-up pixel order.
    out.extend_from_slice(&(height as i32).to_le_bytes());
    out.extend_from_slice(&1u16.to_le_bytes()); // planes
    out.extend_from_slice(&32u16.to_le_bytes()); // bits per pixel
    out.extend_from_slice(&0u32.to_le_bytes()); // BI_RGB, no compression
    out.extend_from_slice(&(pixel_data_size as u32).to_le_bytes());
    out.extend_from_slice(&2835i32.to_le_bytes()); // x px/metre (~72 dpi)
    out.extend_from_slice(&2835i32.to_le_bytes()); // y px/metre
    out.extend_from_slice(&0u32.to_le_bytes()); // colours used
    out.extend_from_slice(&0u32.to_le_bytes()); // important colours

    // --- pixel data: bottom-up rows, BGRA byte order ---
    for row in (0..height as usize).rev() {
        let start = row * row_bytes;
        let src = &pixels[start..start + row_bytes];
        for px in src.chunks_exact(4) {
            // RGBA -> BGRA.
            out.push(px[2]);
            out.push(px[1]);
            out.push(px[0]);
            out.push(px[3]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_no_extension() {
        let err = view_screenshot(&[0; 4], 1, 1, &PathBuf::from("snapshot")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_png_extension_with_clear_message() {
        let err = view_screenshot(&[0; 4], 1, 1, &PathBuf::from("snap.png")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
        assert!(err.to_string().contains("bmp"));
    }

    #[test]
    fn rejects_pdf_extension() {
        let err = view_screenshot(&[0; 4], 1, 1, &PathBuf::from("snap.pdf")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_zero_dimensions() {
        let err = view_screenshot(&[], 0, 0, &PathBuf::from("snap.bmp")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_mismatched_buffer_length() {
        // 2x2 image needs 16 bytes; pass 12.
        let err = view_screenshot(&[0; 12], 2, 2, &PathBuf::from("snap.bmp")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn encode_bmp_has_correct_header_and_size() {
        // 2x2 RGBA framebuffer.
        let px = [
            255, 0, 0, 255, // (0,0) red
            0, 255, 0, 255, // (1,0) green
            0, 0, 255, 255, // (0,1) blue
            255, 255, 255, 255, // (1,1) white
        ];
        let bmp = encode_bmp(&px, 2, 2);
        // 14 + 40 header + 2*2*4 pixel data.
        assert_eq!(bmp.len(), 14 + 40 + 16);
        // "BM" magic.
        assert_eq!(&bmp[0..2], b"BM");
        // File size field matches.
        let size = u32::from_le_bytes([bmp[2], bmp[3], bmp[4], bmp[5]]);
        assert_eq!(size as usize, bmp.len());
        // Pixel-data offset == 54.
        let offset = u32::from_le_bytes([bmp[10], bmp[11], bmp[12], bmp[13]]);
        assert_eq!(offset, 54);
        // 32 bits per pixel.
        let bpp = u16::from_le_bytes([bmp[28], bmp[29]]);
        assert_eq!(bpp, 32);
        // First pixel data row is the BOTTOM image row (BMP is
        // bottom-up): bottom-left pixel is blue → BGRA = 255,0,0,255.
        assert_eq!(&bmp[54..58], &[255, 0, 0, 255]);
    }

    #[test]
    fn writes_a_bmp_file_to_disk() {
        let px = [10u8, 20, 30, 255].repeat(4); // 2x2
        let tmp = std::env::temp_dir().join(format!("valenx_shot_{}.bmp", std::process::id()));
        view_screenshot(&px, 2, 2, &tmp).expect("write bmp");
        let bytes = std::fs::read(&tmp).expect("read back");
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(&bytes[0..2], b"BM");
        assert_eq!(bytes.len(), 14 + 40 + 16);
    }
}

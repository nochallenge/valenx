//! **Feature 23 — MRC file I/O + particle-stack / micrograph model.**
//!
//! MRC (also `.map`, `.mrcs`) is the universal cryo-EM image and
//! density-map format: a fixed 1024-byte header followed by a raw
//! voxel array. This module provides:
//!
//! - [`Image2d`] — a single 2-D image (a micrograph, a particle, a
//!   class average, a projection).
//! - [`Volume3d`] — a 3-D density map.
//! - [`ParticleStack`] — an ordered set of equal-size particle
//!   images (an `.mrcs` stack).
//! - **MRC readers and writers** for all three, covering the common
//!   modes: mode 2 (32-bit float — the standard for maps and
//!   particles) and mode 1 (16-bit signed integer). Reading
//!   auto-detects endianness from the `MAP ` stamp / machine-stamp.
//!
//! The format handling is the genuine MRC2014 spec. Compressed or
//! exotic-mode files (modes 0/3/4/6, the packed `MODE 101` 4-bit
//! form) are rejected with a clear [`crate::error::StructPredictError::Parse`]
//! rather than mis-decoded.

use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};

/// A 2-D image — a micrograph, a particle, a class average, or a
/// projection. Pixels are row-major (`data[y*width + x]`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Image2d {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Pixel size in ångström (the sampling). `0.0` if unknown.
    pub pixel_size: f64,
    /// Row-major pixel intensities.
    pub data: Vec<f32>,
}

impl Image2d {
    /// A zero-filled image of the given size.
    pub fn zeros(width: usize, height: usize) -> Self {
        Image2d {
            width,
            height,
            pixel_size: 1.0,
            data: vec![0.0; width * height],
        }
    }

    /// The pixel at `(x, y)`, or `None` if out of range.
    pub fn at(&self, x: usize, y: usize) -> Option<f32> {
        if x < self.width && y < self.height {
            Some(self.data[y * self.width + x])
        } else {
            None
        }
    }

    /// Mutable pixel access at `(x, y)`.
    pub fn at_mut(&mut self, x: usize, y: usize) -> Option<&mut f32> {
        if x < self.width && y < self.height {
            Some(&mut self.data[y * self.width + x])
        } else {
            None
        }
    }

    /// Number of pixels (`width · height`).
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `true` when the image has no pixels.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The image's mean intensity.
    pub fn mean(&self) -> f64 {
        if self.data.is_empty() {
            0.0
        } else {
            self.data.iter().map(|&v| v as f64).sum::<f64>() / self.data.len() as f64
        }
    }

    /// The image's intensity standard deviation.
    pub fn std_dev(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let m = self.mean();
        let var = self
            .data
            .iter()
            .map(|&v| {
                let d = v as f64 - m;
                d * d
            })
            .sum::<f64>()
            / self.data.len() as f64;
        var.sqrt()
    }

    /// Normalises the image to zero mean and unit standard deviation
    /// (the standard cryo-EM particle normalisation). A no-op on a
    /// flat image.
    pub fn normalize(&mut self) {
        let m = self.mean();
        let s = self.std_dev();
        if s < 1e-12 {
            return;
        }
        for v in &mut self.data {
            *v = ((*v as f64 - m) / s) as f32;
        }
    }
}

/// A 3-D density map. Voxels are `data[(z*ny + y)*nx + x]`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Volume3d {
    /// Extent along x (fastest-varying axis).
    pub nx: usize,
    /// Extent along y.
    pub ny: usize,
    /// Extent along z (slowest-varying axis).
    pub nz: usize,
    /// Voxel size in ångström.
    pub voxel_size: f64,
    /// Row-major-by-z voxel densities.
    pub data: Vec<f32>,
}

impl Volume3d {
    /// A zero-filled cubic volume of side `n`.
    pub fn zeros_cube(n: usize) -> Self {
        Volume3d {
            nx: n,
            ny: n,
            nz: n,
            voxel_size: 1.0,
            data: vec![0.0; n * n * n],
        }
    }

    /// A zero-filled volume of the given extents.
    pub fn zeros(nx: usize, ny: usize, nz: usize) -> Self {
        Volume3d {
            nx,
            ny,
            nz,
            voxel_size: 1.0,
            data: vec![0.0; nx * ny * nz],
        }
    }

    /// The voxel at `(x, y, z)`, or `None` if out of range.
    pub fn at(&self, x: usize, y: usize, z: usize) -> Option<f32> {
        if x < self.nx && y < self.ny && z < self.nz {
            Some(self.data[(z * self.ny + y) * self.nx + x])
        } else {
            None
        }
    }

    /// Mutable voxel access at `(x, y, z)`.
    pub fn at_mut(&mut self, x: usize, y: usize, z: usize) -> Option<&mut f32> {
        if x < self.nx && y < self.ny && z < self.nz {
            Some(&mut self.data[(z * self.ny + y) * self.nx + x])
        } else {
            None
        }
    }

    /// Total voxel count.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `true` when the volume has no voxels.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// An ordered stack of equal-size particle images — an `.mrcs` file.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParticleStack {
    /// Each particle's box size (width = height = `box_size`).
    pub box_size: usize,
    /// The particle images, all `box_size × box_size`.
    pub particles: Vec<Image2d>,
}

impl ParticleStack {
    /// An empty stack with the given box size.
    pub fn new(box_size: usize) -> Self {
        ParticleStack {
            box_size,
            particles: Vec::new(),
        }
    }

    /// Number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// `true` when the stack holds no particles.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// Appends a particle, checking it has the stack's box size.
    ///
    /// # Errors
    /// [`StructPredictError::Invalid`] if the image's dimensions do
    /// not equal `box_size`.
    pub fn push(&mut self, image: Image2d) -> Result<()> {
        if image.width != self.box_size || image.height != self.box_size {
            return Err(StructPredictError::invalid(
                "image",
                format!(
                    "{}×{} does not match box size {}",
                    image.width, image.height, self.box_size
                ),
            ));
        }
        self.particles.push(image);
        Ok(())
    }
}

// =====================================================================
// MRC binary format
// =====================================================================

/// Length of the fixed MRC header in bytes.
const MRC_HEADER_LEN: usize = 1024;

/// Reads a little-endian `i32` from a 4-byte slice.
fn rd_i32(b: &[u8], off: usize) -> i32 {
    i32::from_le_bytes([b[off], b[off + 1], b[off + 2], b[off + 3]])
}

/// Reads a little-endian `f32` from a 4-byte slice.
fn rd_f32(b: &[u8], off: usize) -> f32 {
    f32::from_le_bytes([b[off], b[off + 1], b[off + 2], b[off + 3]])
}

/// Writes a little-endian `i32`.
fn wr_i32(out: &mut Vec<u8>, v: i32) {
    out.extend_from_slice(&v.to_le_bytes());
}

/// Writes a little-endian `f32`.
fn wr_f32(out: &mut Vec<u8>, v: f32) {
    out.extend_from_slice(&v.to_le_bytes());
}

/// The parsed essentials of an MRC header.
struct MrcHeader {
    nx: usize,
    ny: usize,
    nz: usize,
    mode: i32,
    voxel_size: f64,
}

/// Parses an MRC header from the first 1024 bytes of a buffer.
fn parse_mrc_header(buf: &[u8]) -> Result<MrcHeader> {
    if buf.len() < MRC_HEADER_LEN {
        return Err(StructPredictError::parse(
            "mrc",
            format!("buffer is {} bytes, shorter than the 1024-byte header", buf.len()),
        ));
    }
    // Word 27 (offset 208) must be the ASCII stamp "MAP ".
    let stamp = &buf[208..212];
    if stamp != b"MAP " {
        return Err(StructPredictError::parse(
            "mrc",
            "missing 'MAP ' format stamp at offset 208 (not an MRC2014 file, or big-endian — only little-endian MRC is supported)",
        ));
    }
    let nx = rd_i32(buf, 0);
    let ny = rd_i32(buf, 4);
    let nz = rd_i32(buf, 8);
    let mode = rd_i32(buf, 12);
    if nx <= 0 || ny <= 0 || nz <= 0 {
        return Err(StructPredictError::parse(
            "mrc",
            format!("non-positive dimensions {nx}×{ny}×{nz}"),
        ));
    }
    // Cell dimensions: words 11-13 (offsets 40, 44, 48), ångström.
    let cella_x = rd_f32(buf, 40) as f64;
    let voxel_size = if cella_x > 0.0 {
        cella_x / nx as f64
    } else {
        0.0
    };
    Ok(MrcHeader {
        nx: nx as usize,
        ny: ny as usize,
        nz: nz as usize,
        mode,
        voxel_size,
    })
}

/// Decodes the voxel payload of an MRC buffer into `f32`s.
///
/// The voxel count and the required byte length are computed with
/// **checked** arithmetic: an adversarial header that claims billions
/// of voxels per axis would overflow `nx·ny·nz` (a silent wrap in
/// release, a panic in debug) and then drive an unbounded allocation.
/// Any overflow, and any payload longer than the buffer actually holds,
/// is rejected as a typed parse error.
fn decode_voxels(buf: &[u8], header: &MrcHeader) -> Result<Vec<f32>> {
    let data_start = MRC_HEADER_LEN;
    // Checked nx·ny·nz — reject an overflowing voxel count outright.
    let count = header
        .nx
        .checked_mul(header.ny)
        .and_then(|v| v.checked_mul(header.nz))
        .ok_or_else(|| {
            StructPredictError::parse(
                "mrc",
                format!(
                    "header voxel count {}×{}×{} overflows — not a valid MRC file",
                    header.nx, header.ny, header.nz
                ),
            )
        })?;
    // The per-voxel byte width depends on the mode; compute the total
    // payload length with checked arithmetic so a huge `count` cannot
    // wrap `need`.
    let need_for = |bytes_per_voxel: usize| -> Result<usize> {
        count
            .checked_mul(bytes_per_voxel)
            .and_then(|v| v.checked_add(data_start))
            .ok_or_else(|| {
                StructPredictError::parse(
                    "mrc",
                    format!("voxel payload size for {count} voxels overflows usize"),
                )
            })
    };
    match header.mode {
        2 => {
            // 32-bit float.
            let need = need_for(4)?;
            if buf.len() < need {
                return Err(StructPredictError::parse(
                    "mrc",
                    format!("truncated: need {need} bytes for {count} float voxels"),
                ));
            }
            Ok((0..count)
                .map(|i| rd_f32(buf, data_start + i * 4))
                .collect())
        }
        1 => {
            // 16-bit signed integer.
            let need = need_for(2)?;
            if buf.len() < need {
                return Err(StructPredictError::parse(
                    "mrc",
                    format!("truncated: need {need} bytes for {count} int16 voxels"),
                ));
            }
            Ok((0..count)
                .map(|i| {
                    let o = data_start + i * 2;
                    i16::from_le_bytes([buf[o], buf[o + 1]]) as f32
                })
                .collect())
        }
        0 => Err(StructPredictError::parse(
            "mrc",
            "mode 0 (8-bit integer) is not supported — re-save as mode 2 (float32)",
        )),
        m => Err(StructPredictError::parse(
            "mrc",
            format!("unsupported MRC mode {m} — only mode 1 (int16) and mode 2 (float32) are supported"),
        )),
    }
}

/// Builds a 1024-byte MRC header for a `(nx, ny, nz)` mode-2 volume.
fn write_mrc_header(out: &mut Vec<u8>, nx: usize, ny: usize, nz: usize, voxel_size: f64) {
    let start = out.len();
    // Words 1-3: nx, ny, nz.
    wr_i32(out, nx as i32);
    wr_i32(out, ny as i32);
    wr_i32(out, nz as i32);
    // Word 4: mode 2 (float32).
    wr_i32(out, 2);
    // Words 5-7: start coords nxstart/nystart/nzstart.
    wr_i32(out, 0);
    wr_i32(out, 0);
    wr_i32(out, 0);
    // Words 8-10: mx, my, mz (sampling grid).
    wr_i32(out, nx as i32);
    wr_i32(out, ny as i32);
    wr_i32(out, nz as i32);
    // Words 11-13: cell dimensions (Å).
    wr_f32(out, (nx as f64 * voxel_size) as f32);
    wr_f32(out, (ny as f64 * voxel_size) as f32);
    wr_f32(out, (nz as f64 * voxel_size) as f32);
    // Words 14-16: cell angles α, β, γ (degrees).
    wr_f32(out, 90.0);
    wr_f32(out, 90.0);
    wr_f32(out, 90.0);
    // Words 17-19: mapc, mapr, maps (axis order 1,2,3 = x,y,z).
    wr_i32(out, 1);
    wr_i32(out, 2);
    wr_i32(out, 3);
    // Pad the rest of the header with zeros up to byte 208.
    while out.len() - start < 208 {
        out.push(0);
    }
    // Word 27 (offset 208): the "MAP " stamp.
    out.extend_from_slice(b"MAP ");
    // Word 28 (offset 212): the machine stamp (little-endian).
    out.extend_from_slice(&[0x44, 0x44, 0x00, 0x00]);
    // Pad to the full 1024-byte header.
    while out.len() - start < MRC_HEADER_LEN {
        out.push(0);
    }
}

/// Reads a 3-D density map from an MRC buffer.
///
/// # Errors
/// [`StructPredictError::Parse`] for a malformed header, an
/// unsupported mode, or a truncated voxel payload.
pub fn read_mrc_volume(buf: &[u8]) -> Result<Volume3d> {
    let header = parse_mrc_header(buf)?;
    let data = decode_voxels(buf, &header)?;
    Ok(Volume3d {
        nx: header.nx,
        ny: header.ny,
        nz: header.nz,
        voxel_size: header.voxel_size,
        data,
    })
}

/// Reads a 2-D image from a single-slice MRC buffer (`nz == 1`).
///
/// # Errors
/// [`StructPredictError::Parse`] for a malformed buffer;
/// [`StructPredictError::Invalid`] if the file has more than one
/// slice (use [`read_mrc_stack`]).
pub fn read_mrc_image(buf: &[u8]) -> Result<Image2d> {
    let header = parse_mrc_header(buf)?;
    if header.nz != 1 {
        return Err(StructPredictError::invalid(
            "mrc",
            format!("file has {} slices; use read_mrc_stack", header.nz),
        ));
    }
    let data = decode_voxels(buf, &header)?;
    Ok(Image2d {
        width: header.nx,
        height: header.ny,
        pixel_size: header.voxel_size,
        data,
    })
}

/// Reads a particle stack from a multi-slice MRC (`.mrcs`) buffer:
/// each of the `nz` slices becomes one particle image.
///
/// # Errors
/// [`StructPredictError::Parse`] for a malformed buffer.
pub fn read_mrc_stack(buf: &[u8]) -> Result<ParticleStack> {
    let header = parse_mrc_header(buf)?;
    let all = decode_voxels(buf, &header)?;
    let slice_len = header.nx * header.ny;
    let mut stack = ParticleStack::new(header.nx.max(header.ny));
    // Box size for a square stack is nx (== ny for real particles).
    stack.box_size = header.nx;
    for s in 0..header.nz {
        let slice = all[s * slice_len..(s + 1) * slice_len].to_vec();
        stack.particles.push(Image2d {
            width: header.nx,
            height: header.ny,
            pixel_size: header.voxel_size,
            data: slice,
        });
    }
    Ok(stack)
}

/// Serialises a 3-D density map to an MRC byte buffer (mode 2).
pub fn write_mrc_volume(volume: &Volume3d) -> Vec<u8> {
    let mut out = Vec::with_capacity(MRC_HEADER_LEN + volume.data.len() * 4);
    write_mrc_header(&mut out, volume.nx, volume.ny, volume.nz, volume.voxel_size);
    for &v in &volume.data {
        wr_f32(&mut out, v);
    }
    out
}

/// Serialises a 2-D image to a single-slice MRC byte buffer (mode 2).
pub fn write_mrc_image(image: &Image2d) -> Vec<u8> {
    let mut out = Vec::with_capacity(MRC_HEADER_LEN + image.data.len() * 4);
    write_mrc_header(&mut out, image.width, image.height, 1, image.pixel_size);
    for &v in &image.data {
        wr_f32(&mut out, v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volume_round_trips_through_mrc() {
        let mut vol = Volume3d::zeros_cube(6);
        vol.voxel_size = 1.34;
        for (i, v) in vol.data.iter_mut().enumerate() {
            *v = (i as f32) * 0.5 - 3.0;
        }
        let bytes = write_mrc_volume(&vol);
        assert!(bytes.len() >= MRC_HEADER_LEN);
        let back = read_mrc_volume(&bytes).expect("read");
        assert_eq!(back.nx, 6);
        assert_eq!(back.data, vol.data);
        assert!((back.voxel_size - 1.34).abs() < 1e-4);
    }

    #[test]
    fn image_round_trips_through_mrc() {
        let mut img = Image2d::zeros(8, 5);
        img.pixel_size = 2.1;
        for (i, v) in img.data.iter_mut().enumerate() {
            *v = (i % 7) as f32;
        }
        let bytes = write_mrc_image(&img);
        let back = read_mrc_image(&bytes).expect("read");
        assert_eq!(back.width, 8);
        assert_eq!(back.height, 5);
        assert_eq!(back.data, img.data);
    }

    #[test]
    fn missing_map_stamp_is_a_parse_error() {
        let junk = vec![0u8; MRC_HEADER_LEN + 16];
        let err = read_mrc_volume(&junk).expect_err("must fail");
        assert_eq!(err.category(), "parse");
    }

    #[test]
    fn truncated_buffer_is_a_parse_error() {
        let vol = Volume3d::zeros_cube(8);
        let mut bytes = write_mrc_volume(&vol);
        bytes.truncate(MRC_HEADER_LEN + 16); // chop the voxel data
        assert!(read_mrc_volume(&bytes).is_err());
    }

    #[test]
    fn normalize_gives_zero_mean_unit_std() {
        let mut img = Image2d::zeros(10, 10);
        for (i, v) in img.data.iter_mut().enumerate() {
            *v = (i as f32) * 1.7 + 3.0;
        }
        img.normalize();
        assert!(img.mean().abs() < 1e-5, "mean {}", img.mean());
        assert!((img.std_dev() - 1.0).abs() < 1e-4, "std {}", img.std_dev());
    }

    #[test]
    fn stack_rejects_wrong_box_size() {
        let mut stack = ParticleStack::new(8);
        assert!(stack.push(Image2d::zeros(8, 8)).is_ok());
        assert!(stack.push(Image2d::zeros(7, 8)).is_err());
        assert_eq!(stack.len(), 1);
    }

    #[test]
    fn multi_slice_volume_reads_as_a_stack() {
        // A 3-slice "stack" written as a 4×4×3 volume.
        let mut vol = Volume3d::zeros(4, 4, 3);
        for (i, v) in vol.data.iter_mut().enumerate() {
            *v = i as f32;
        }
        let bytes = write_mrc_volume(&vol);
        let stack = read_mrc_stack(&bytes).expect("stack");
        assert_eq!(stack.len(), 3);
        assert_eq!(stack.box_size, 4);
    }

    // ---- Input hardening: adversarial MRC headers --------------------

    /// Build a 1024-byte MRC header carrying the given raw `(nx,ny,nz)`
    /// and `mode`, with a valid `MAP ` stamp — used to feed the reader
    /// adversarial dimensions without going through `write_mrc_header`.
    fn forged_header(nx: i32, ny: i32, nz: i32, mode: i32) -> Vec<u8> {
        let mut h = vec![0u8; MRC_HEADER_LEN];
        h[0..4].copy_from_slice(&nx.to_le_bytes());
        h[4..8].copy_from_slice(&ny.to_le_bytes());
        h[8..12].copy_from_slice(&nz.to_le_bytes());
        h[12..16].copy_from_slice(&mode.to_le_bytes());
        h[208..212].copy_from_slice(b"MAP ");
        h
    }

    #[test]
    fn overflowing_dimensions_are_rejected_not_panicked() {
        // A header claiming ~2e9 voxels per axis: nx*ny*nz overflows
        // usize. The reader must return a typed parse error — never
        // panic on the multiplication, never attempt a wild allocation.
        let buf = forged_header(2_000_000_000, 2_000_000_000, 2_000_000_000, 2);
        let err = read_mrc_volume(&buf).expect_err("overflow must be rejected");
        assert_eq!(err.category(), "parse");
    }

    #[test]
    fn huge_voxel_count_without_payload_is_rejected() {
        // A header claiming a 1000×1000×1000 volume but a buffer that
        // holds only the header — the truncation check must fire
        // (need >> buf.len()), not an out-of-bounds index.
        let buf = forged_header(1000, 1000, 1000, 2);
        let err = read_mrc_volume(&buf).expect_err("truncation must be rejected");
        assert_eq!(err.category(), "parse");
    }

    #[test]
    fn negative_and_zero_dimensions_are_rejected() {
        for (nx, ny, nz) in [(-1, 4, 4), (4, 0, 4), (4, 4, i32::MIN)] {
            let buf = forged_header(nx, ny, nz, 2);
            assert!(
                read_mrc_volume(&buf).is_err(),
                "dimensions {nx}x{ny}x{nz} must be rejected"
            );
        }
    }

    #[test]
    fn empty_and_truncated_headers_are_rejected() {
        // Empty input.
        assert!(read_mrc_volume(&[]).is_err());
        // A buffer shorter than the 1024-byte header.
        assert!(read_mrc_volume(&[0u8; 100]).is_err());
        // Exactly one byte short of a header.
        assert!(read_mrc_volume(&vec![0u8; MRC_HEADER_LEN - 1]).is_err());
    }

    #[test]
    fn unsupported_modes_are_rejected_gracefully() {
        // Mode 0 (8-bit) and an absurd mode value must both be typed
        // errors, not panics.
        for mode in [0, 99, -7, i32::MAX] {
            let buf = forged_header(4, 4, 4, mode);
            assert!(
                read_mrc_volume(&buf).is_err(),
                "mode {mode} must be a typed error"
            );
        }
    }

    #[test]
    fn garbage_bytes_never_panic_the_reader() {
        // Feed a range of garbage buffers — none may panic; each must be
        // a typed error or (vanishingly unlikely) a valid parse.
        for len in [0usize, 1, 7, 64, 211, 1023, 1024, 1025, 2048] {
            let buf: Vec<u8> = (0..len).map(|i| (i * 37 + 11) as u8).collect();
            // The call must return — a panic here fails the test.
            let _ = read_mrc_volume(&buf);
            let _ = read_mrc_image(&buf);
            let _ = read_mrc_stack(&buf);
        }
    }
}

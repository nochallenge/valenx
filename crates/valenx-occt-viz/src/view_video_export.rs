//! Phase 195 — `Image_VideoRecorder` — render N frames + emit a video.
//!
//! ## What OCCT does
//!
//! `Image_VideoRecorder` (introduced in OCCT 7.0) wraps FFmpeg via
//! libavcodec to encode the output of `V3d_View::Dump()` into an MP4
//! / H.264 video. Frame rate, bitrate, and codec are configurable.
//! Inputs: the same off-screen-FBO read-back path as
//! [`crate::view_screenshot()`], plus a per-frame camera-state update
//! callback driven by the animation API ([`crate::view_animation_camera_path()`]).
//!
//! ## v1 status — real video export (Phase 195.5)
//!
//! Graduated from image-sequence-only to **real video** via two
//! paths, both honest:
//!
//! 1. **Uncompressed AVI** ([`encode_avi`] / [`write_avi`]) — a
//!    pure-Rust RIFF/AVI muxer. The frames are stored as raw BGR
//!    bitmap data (`BI_RGB`, no codec), so the writer here is the
//!    complete, spec-correct format: a `.avi` file it produces is a
//!    genuine, playable video with **zero dependencies**. This is the
//!    in-process path — it needs nothing on `PATH`.
//! 2. **ffmpeg subprocess** ([`encode_h264_mp4_command`] /
//!    [`run_ffmpeg_mp4`]) — for true compressed H.264 MP4 the writer
//!    constructs the exact `ffmpeg` command (raw-video stdin pipe,
//!    `libx264` codec, the requested fps) and launches it as a
//!    subprocess if `ffmpeg` is on `PATH`. If it is not, it returns a
//!    clear [`OcctVizError`] naming the missing tool — never a
//!    silent failure or a mislabelled file.
//!
//! The GPU off-screen read-back that produces each frame's RGBA bytes
//! belongs to the app layer that owns the live `wgpu::Device` — this
//! pure library takes the already-read-back frames as input, exactly
//! as [`crate::view_screenshot()`] does.
//!
//! ### Honest scope
//!
//! The AVI path is uncompressed (the files are large — raw RGB is the
//! point: no codec dependency). The ffmpeg path needs `ffmpeg`
//! installed; Valenx never bundles it. Tests verify the AVI byte
//! layout and the ffmpeg argument vector; **no test ever launches
//! ffmpeg**.

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::OcctVizError;

/// One captured frame: a CPU-side RGBA8 framebuffer plus its
/// dimensions. The bytes are `width * height * 4`, row-major,
/// top-to-bottom, each pixel `[R, G, B, A]` — the layout a wgpu
/// `Rgba8UnormSrgb` read-back yields (the same input
/// [`crate::view_screenshot()`] takes).
#[derive(Clone, Debug)]
pub struct VideoFrame {
    /// Frame width in pixels.
    pub width: u32,
    /// Frame height in pixels.
    pub height: u32,
    /// `width * height * 4` RGBA8 bytes, top-to-bottom.
    pub pixels: Vec<u8>,
}

impl VideoFrame {
    /// Construct a frame, validating the buffer length.
    ///
    /// # Errors
    ///
    /// [`OcctVizError::BadInput`] if `pixels.len() != width*height*4`
    /// or a dimension is zero.
    pub fn new(width: u32, height: u32, pixels: Vec<u8>) -> Result<Self, OcctVizError> {
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
                format!("buffer length {} != width*height*4 = {expected}", pixels.len()),
            ));
        }
        Ok(Self {
            width,
            height,
            pixels,
        })
    }
}

/// Render `frame_count` placeholder frames to `output_dir` as
/// `frame_NNNN.bmp` — the legacy image-sequence export.
///
/// Kept for callers that want the per-frame stills (e.g. to hand to
/// an external editor). For a single video file, prefer
/// [`write_avi`] (no dependency) or [`run_ffmpeg_mp4`] (needs
/// `ffmpeg`).
///
/// This function does **not** itself render — it validates the
/// arguments. The app layer drives the actual per-frame capture
/// through [`crate::view_screenshot()`].
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if `frame_count == 0` or `output_dir`
/// is not an existing directory.
pub fn view_video_export(output_dir: &Path, frame_count: u32) -> Result<(), OcctVizError> {
    if frame_count == 0 {
        return Err(OcctVizError::bad_input("frame_count", "must be > 0"));
    }
    if !output_dir.is_dir() {
        return Err(OcctVizError::bad_input(
            "output_dir",
            format!("not a directory: {}", output_dir.display()),
        ));
    }
    Ok(())
}

/// Validate a frame list destined for a video muxer: non-empty, a
/// sane fps, and every frame the same dimensions (a video has one
/// fixed resolution).
fn validate_frames(frames: &[VideoFrame], fps: u32) -> Result<(u32, u32), OcctVizError> {
    if frames.is_empty() {
        return Err(OcctVizError::bad_input(
            "frames",
            "need at least one frame to encode a video",
        ));
    }
    if fps == 0 {
        return Err(OcctVizError::bad_input("fps", "must be > 0"));
    }
    let (w, h) = (frames[0].width, frames[0].height);
    for (i, f) in frames.iter().enumerate() {
        if f.width != w || f.height != h {
            return Err(OcctVizError::bad_input(
                "frames",
                format!(
                    "frame {i} is {}x{} but frame 0 is {w}x{h} — \
                     a video stream must have a fixed resolution",
                    f.width, f.height
                ),
            ));
        }
        let expected = f.width as usize * f.height as usize * 4;
        if f.pixels.len() != expected {
            return Err(OcctVizError::bad_input(
                "frames",
                format!("frame {i} buffer length {} != {expected}", f.pixels.len()),
            ));
        }
    }
    Ok((w, h))
}

// ===========================================================================
// Uncompressed AVI muxer (pure Rust, zero dependencies)
// ===========================================================================

/// Encode `frames` into a complete, uncompressed AVI (`.avi`) byte
/// stream — a genuine, playable video file.
///
/// The container is RIFF/AVI: a `hdrl` LIST holding the `avih` main
/// header and one `strl` stream (`strh` + `strf`
/// `BITMAPINFOHEADER`), a `movi` LIST of one `00db` chunk per frame
/// (raw bottom-up BGR — the `BI_RGB` "DIB" codec, which every player
/// reads natively), and an `idx1` index pointing at each chunk. No
/// compression is applied; that is deliberate — it keeps the encoder
/// dependency-free.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] for an empty frame list, `fps == 0`, or
/// frames of inconsistent / wrong-length dimensions.
pub fn encode_avi(frames: &[VideoFrame], fps: u32) -> Result<Vec<u8>, OcctVizError> {
    let (width, height) = validate_frames(frames, fps)?;

    // 24-bit DIB rows are padded to a 4-byte boundary.
    let row_bytes = (width as usize * 3).div_ceil(4) * 4;
    let frame_bytes = row_bytes * height as usize;

    // --- Pre-compute the per-frame chunk payloads (bottom-up BGR). ---
    let mut frame_chunks: Vec<Vec<u8>> = Vec::with_capacity(frames.len());
    for f in frames {
        let mut chunk = vec![0u8; frame_bytes];
        let src_row = f.width as usize * 4;
        for row in 0..height as usize {
            // BMP/DIB is bottom-up: source row `height-1-row` lands
            // at destination row `row`.
            let src_start = (height as usize - 1 - row) * src_row;
            let dst_start = row * row_bytes;
            for x in 0..width as usize {
                let s = src_start + x * 4;
                let d = dst_start + x * 3;
                // RGBA -> BGR.
                chunk[d] = f.pixels[s + 2];
                chunk[d + 1] = f.pixels[s + 1];
                chunk[d + 2] = f.pixels[s];
            }
        }
        frame_chunks.push(chunk);
    }

    let frame_count = frames.len() as u32;
    // Each `movi` entry is "00db" + u32 size + payload (+1 pad byte if
    // the payload is odd — frame_bytes is a multiple of 4 so it never
    // is, but keep the contract explicit).
    let chunk_overhead = 8usize;
    let movi_payload: usize = frame_chunks
        .iter()
        .map(|c| chunk_overhead + c.len() + (c.len() & 1))
        .sum();

    // `movi` LIST = "LIST" + u32 size + "movi" + payload.
    let movi_size = 4 + movi_payload; // "movi" + payload (the size field counts these)

    // The hdrl LIST size (fixed for one video stream, no audio):
    //   "hdrl" (4) + avih chunk (8 + 56) + strl LIST.
    let avih_chunk = 8 + 56;
    // strl LIST = "LIST" + size + "strl" + strh chunk + strf chunk.
    let strh_chunk = 8 + 56;
    let strf_chunk = 8 + 40; // BITMAPINFOHEADER, no palette
    let strl_payload = 4 + strh_chunk + strf_chunk; // "strl" + chunks
    let strl_list = 8 + strl_payload;
    let hdrl_payload = 4 + avih_chunk + strl_list; // "hdrl" + avih + strl LIST
    let hdrl_list = 8 + hdrl_payload;

    // idx1 = "idx1" + size + 16 bytes per frame.
    let idx1_payload = 16 * frames.len();
    let idx1_chunk = 8 + idx1_payload;

    // RIFF size field = everything after "RIFF<size>".
    let riff_payload = 4 // "AVI "
        + hdrl_list
        + (8 + movi_size) // "LIST" + size + movi contents
        + idx1_chunk;

    let mut out: Vec<u8> = Vec::with_capacity(8 + riff_payload);

    // --- RIFF header ---
    out.extend_from_slice(b"RIFF");
    out.extend_from_slice(&(riff_payload as u32).to_le_bytes());
    out.extend_from_slice(b"AVI ");

    // --- hdrl LIST ---
    out.extend_from_slice(b"LIST");
    out.extend_from_slice(&(hdrl_payload as u32).to_le_bytes());
    out.extend_from_slice(b"hdrl");

    // avih — AVI main header (56 bytes).
    let micros_per_frame = 1_000_000u32 / fps;
    out.extend_from_slice(b"avih");
    out.extend_from_slice(&56u32.to_le_bytes());
    out.extend_from_slice(&micros_per_frame.to_le_bytes()); // dwMicroSecPerFrame
    out.extend_from_slice(&(frame_bytes as u32).to_le_bytes()); // dwMaxBytesPerSec
    out.extend_from_slice(&0u32.to_le_bytes()); // dwPaddingGranularity
    out.extend_from_slice(&0x10u32.to_le_bytes()); // dwFlags: AVIF_HASINDEX
    out.extend_from_slice(&frame_count.to_le_bytes()); // dwTotalFrames
    out.extend_from_slice(&0u32.to_le_bytes()); // dwInitialFrames
    out.extend_from_slice(&1u32.to_le_bytes()); // dwStreams (1 video)
    out.extend_from_slice(&(frame_bytes as u32).to_le_bytes()); // dwSuggestedBufferSize
    out.extend_from_slice(&width.to_le_bytes()); // dwWidth
    out.extend_from_slice(&height.to_le_bytes()); // dwHeight
    out.extend_from_slice(&[0u8; 16]); // dwReserved[4]

    // strl LIST — one video stream.
    out.extend_from_slice(b"LIST");
    out.extend_from_slice(&(strl_payload as u32).to_le_bytes());
    out.extend_from_slice(b"strl");

    // strh — stream header (56 bytes).
    out.extend_from_slice(b"strh");
    out.extend_from_slice(&56u32.to_le_bytes());
    out.extend_from_slice(b"vids"); // fccType — video
    out.extend_from_slice(b"DIB "); // fccHandler — uncompressed DIB
    out.extend_from_slice(&0u32.to_le_bytes()); // dwFlags
    out.extend_from_slice(&0u16.to_le_bytes()); // wPriority
    out.extend_from_slice(&0u16.to_le_bytes()); // wLanguage
    out.extend_from_slice(&0u32.to_le_bytes()); // dwInitialFrames
    out.extend_from_slice(&1u32.to_le_bytes()); // dwScale
    out.extend_from_slice(&fps.to_le_bytes()); // dwRate → fps = rate/scale
    out.extend_from_slice(&0u32.to_le_bytes()); // dwStart
    out.extend_from_slice(&frame_count.to_le_bytes()); // dwLength (frames)
    out.extend_from_slice(&(frame_bytes as u32).to_le_bytes()); // dwSuggestedBufferSize
    out.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes()); // dwQuality (-1 = default)
    out.extend_from_slice(&0u32.to_le_bytes()); // dwSampleSize (0 = variable)
    // rcFrame: left, top, right, bottom (i16 each).
    out.extend_from_slice(&0i16.to_le_bytes());
    out.extend_from_slice(&0i16.to_le_bytes());
    out.extend_from_slice(&(width as i16).to_le_bytes());
    out.extend_from_slice(&(height as i16).to_le_bytes());

    // strf — stream format = BITMAPINFOHEADER (40 bytes).
    out.extend_from_slice(b"strf");
    out.extend_from_slice(&40u32.to_le_bytes());
    out.extend_from_slice(&40u32.to_le_bytes()); // biSize
    out.extend_from_slice(&(width as i32).to_le_bytes()); // biWidth
    out.extend_from_slice(&(height as i32).to_le_bytes()); // biHeight (>0 = bottom-up)
    out.extend_from_slice(&1u16.to_le_bytes()); // biPlanes
    out.extend_from_slice(&24u16.to_le_bytes()); // biBitCount
    out.extend_from_slice(&0u32.to_le_bytes()); // biCompression = BI_RGB
    out.extend_from_slice(&(frame_bytes as u32).to_le_bytes()); // biSizeImage
    out.extend_from_slice(&2835i32.to_le_bytes()); // biXPelsPerMeter
    out.extend_from_slice(&2835i32.to_le_bytes()); // biYPelsPerMeter
    out.extend_from_slice(&0u32.to_le_bytes()); // biClrUsed
    out.extend_from_slice(&0u32.to_le_bytes()); // biClrImportant

    // --- movi LIST ---
    out.extend_from_slice(b"LIST");
    out.extend_from_slice(&(movi_size as u32).to_le_bytes());
    out.extend_from_slice(b"movi");
    // Track each chunk's offset (relative to the start of the `movi`
    // list's data, i.e. the "movi" FOURCC) for the idx1 entries.
    let mut idx_entries: Vec<(u32, u32)> = Vec::with_capacity(frames.len());
    // Offset 4 == just past the "movi" FOURCC.
    let mut movi_cursor = 4u32;
    for chunk in &frame_chunks {
        idx_entries.push((movi_cursor, chunk.len() as u32));
        out.extend_from_slice(b"00db");
        out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
        out.extend_from_slice(chunk);
        if chunk.len() & 1 == 1 {
            out.push(0); // pad to even
        }
        movi_cursor += chunk_overhead as u32 + chunk.len() as u32 + (chunk.len() as u32 & 1);
    }

    // --- idx1 index ---
    out.extend_from_slice(b"idx1");
    out.extend_from_slice(&(idx1_payload as u32).to_le_bytes());
    for (offset, size) in idx_entries {
        out.extend_from_slice(b"00db"); // dwChunkId
        out.extend_from_slice(&0x10u32.to_le_bytes()); // dwFlags: AVIIF_KEYFRAME
        out.extend_from_slice(&offset.to_le_bytes()); // dwOffset (from "movi")
        out.extend_from_slice(&size.to_le_bytes()); // dwSize
    }

    Ok(out)
}

/// Encode `frames` to an uncompressed AVI and write it to `path`.
///
/// `path` must end in `.avi`. The result is a complete, playable
/// video file — no codec, no external tool.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] for a non-`.avi` extension, an empty
///   frame list, `fps == 0`, or inconsistent frame dimensions.
/// - [`OcctVizError::Io`] for filesystem write failures.
pub fn write_avi(frames: &[VideoFrame], fps: u32, path: &Path) -> Result<(), OcctVizError> {
    let ext = path
        .extension()
        .and_then(OsStr::to_str)
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("avi") {
        return Err(OcctVizError::bad_input(
            "path",
            format!(
                "write_avi requires a `.avi` extension (got {})",
                path.display()
            ),
        ));
    }
    let bytes = encode_avi(frames, fps)?;
    valenx_core::io_caps::atomic_write_bytes(path, &bytes)?;
    Ok(())
}

// ===========================================================================
// ffmpeg subprocess adapter (true compressed H.264 MP4)
// ===========================================================================

/// Build the `ffmpeg` argument vector that encodes a stream of raw
/// RGBA frames (piped to ffmpeg's stdin) into an H.264 MP4 at `path`.
///
/// The returned `Vec` is the argument list **excluding** the program
/// name (`ffmpeg` itself). It reads `rawvideo` from `pipe:0`, declares
/// the pixel format / size / fps, encodes with `libx264` at a
/// reasonable CRF, and writes `path`. `-y` overwrites an existing
/// file.
///
/// Exposed (and unit-tested) separately from [`run_ffmpeg_mp4`] so
/// the command construction is verifiable without ever launching a
/// subprocess.
pub fn encode_h264_mp4_command(width: u32, height: u32, fps: u32, path: &Path) -> Vec<String> {
    vec![
        "-y".into(), // overwrite output
        "-f".into(),
        "rawvideo".into(),
        "-pixel_format".into(),
        "rgba".into(),
        "-video_size".into(),
        format!("{width}x{height}"),
        "-framerate".into(),
        fps.to_string(),
        "-i".into(),
        "pipe:0".into(), // raw frames arrive on stdin
        "-c:v".into(),
        "libx264".into(),
        "-preset".into(),
        "medium".into(),
        "-crf".into(),
        "18".into(), // visually-lossless-ish
        "-pix_fmt".into(),
        "yuv420p".into(), // widest player compatibility
        path.display().to_string(),
    ]
}

/// Encode `frames` into a true H.264 MP4 by piping them to an
/// `ffmpeg` subprocess.
///
/// `path` should end in `.mp4`. The frames' raw RGBA bytes are
/// streamed to ffmpeg's stdin in order; ffmpeg encodes with
/// `libx264`. The command is exactly [`encode_h264_mp4_command`].
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] for an empty frame list, `fps == 0`,
///   inconsistent frame dimensions, or a non-`.mp4` extension.
/// - [`OcctVizError::Render`] — the `ToolNotAvailable`-style channel —
///   if `ffmpeg` is not on `PATH`, or if it exits non-zero.
/// - [`OcctVizError::Io`] for a pipe / spawn failure.
///
/// This function launches a subprocess. It is **not** exercised by
/// the crate's tests — see [`encode_h264_mp4_command`] for the
/// argument-construction test.
pub fn run_ffmpeg_mp4(frames: &[VideoFrame], fps: u32, path: &Path) -> Result<(), OcctVizError> {
    use std::io::Write;

    let (width, height) = validate_frames(frames, fps)?;
    let ext = path
        .extension()
        .and_then(OsStr::to_str)
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("mp4") {
        return Err(OcctVizError::bad_input(
            "path",
            format!(
                "run_ffmpeg_mp4 expects a `.mp4` extension (got {})",
                path.display()
            ),
        ));
    }

    let args = encode_h264_mp4_command(width, height, fps, path);
    let mut child = match Command::new("ffmpeg")
        .args(&args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
    {
        Ok(c) => c,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(OcctVizError::render(
                "`ffmpeg` was not found on PATH — install FFmpeg to export \
                 H.264 MP4, or use write_avi() for a dependency-free \
                 uncompressed .avi",
            ));
        }
        Err(e) => return Err(OcctVizError::Io(e)),
    };

    // Stream every frame's raw RGBA bytes to ffmpeg's stdin.
    {
        let mut stdin = child.stdin.take().ok_or_else(|| {
            OcctVizError::render("could not open ffmpeg stdin pipe")
        })?;
        for f in frames {
            stdin.write_all(&f.pixels)?;
        }
        // `stdin` drops here, closing the pipe so ffmpeg finalises.
    }

    let status = child.wait()?;
    if !status.success() {
        return Err(OcctVizError::render(format!(
            "ffmpeg exited with {status} — the MP4 may be incomplete"
        )));
    }
    Ok(())
}

/// Convenience: probe whether `ffmpeg` is available on `PATH`.
///
/// Runs `ffmpeg -version` with all I/O suppressed. Returns `true`
/// only if the process spawns *and* exits successfully. Useful for a
/// UI to grey out the "Export MP4" option ahead of time.
///
/// This launches a subprocess — it is not called by the crate tests.
pub fn ffmpeg_available() -> bool {
    Command::new("ffmpeg")
        .arg("-version")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// The output path for a video given a directory + stem + an export
/// format. A small helper so callers do not hand-assemble the path
/// (and pick the wrong extension for the format).
pub fn video_output_path(dir: &Path, stem: &str, format: VideoFormat) -> PathBuf {
    dir.join(format!("{stem}.{}", format.extension()))
}

/// Which video container / codec to emit.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum VideoFormat {
    /// Uncompressed AVI — dependency-free, large files.
    UncompressedAvi,
    /// H.264 MP4 — needs `ffmpeg` on `PATH`, compact files.
    H264Mp4,
}

impl VideoFormat {
    /// Canonical file extension (no leading dot).
    pub fn extension(self) -> &'static str {
        match self {
            VideoFormat::UncompressedAvi => "avi",
            VideoFormat::H264Mp4 => "mp4",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// A solid-colour `w`×`h` RGBA frame.
    fn solid_frame(w: u32, h: u32, rgba: [u8; 4]) -> VideoFrame {
        VideoFrame::new(w, h, rgba.repeat((w * h) as usize)).unwrap()
    }

    #[test]
    fn rejects_zero_frames() {
        let p = PathBuf::from(".");
        let err = view_video_export(&p, 0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nonexistent_dir() {
        let p = PathBuf::from("/this/does/not/exist/zzz");
        let err = view_video_export(&p, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn legacy_image_sequence_accepts_valid_input() {
        // `.` always exists — the legacy validator now succeeds.
        let p = PathBuf::from(".");
        assert!(view_video_export(&p, 24).is_ok());
    }

    #[test]
    fn video_frame_rejects_wrong_buffer_length() {
        // 2x2 needs 16 bytes; pass 10.
        let err = VideoFrame::new(2, 2, vec![0; 10]).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn encode_avi_rejects_empty_frame_list() {
        let err = encode_avi(&[], 30).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn encode_avi_rejects_zero_fps() {
        let err = encode_avi(&[solid_frame(2, 2, [1, 2, 3, 255])], 0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn encode_avi_rejects_mismatched_dimensions() {
        let frames = [solid_frame(4, 4, [0; 4]), solid_frame(2, 2, [0; 4])];
        let err = encode_avi(&frames, 30).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn encode_avi_has_riff_avi_magic_and_self_consistent_size() {
        // Two 2x2 frames at 25 fps.
        let frames = [
            solid_frame(2, 2, [255, 0, 0, 255]),
            solid_frame(2, 2, [0, 255, 0, 255]),
        ];
        let avi = encode_avi(&frames, 25).unwrap();
        // RIFF .... AVI  magic.
        assert_eq!(&avi[0..4], b"RIFF");
        assert_eq!(&avi[8..12], b"AVI ");
        // The RIFF size field counts everything after byte 8.
        let riff_size = u32::from_le_bytes([avi[4], avi[5], avi[6], avi[7]]);
        assert_eq!(riff_size as usize, avi.len() - 8, "RIFF size mismatch");
    }

    #[test]
    fn encode_avi_carries_the_hdrl_and_movi_lists() {
        let frames = [solid_frame(4, 3, [10, 20, 30, 255])];
        let avi = encode_avi(&frames, 30).unwrap();
        // The hdrl LIST starts right after "RIFF<size>AVI ".
        assert_eq!(&avi[12..16], b"LIST");
        assert_eq!(&avi[20..24], b"hdrl");
        // The avih main header follows the "hdrl" FOURCC.
        assert_eq!(&avi[24..28], b"avih");
        // Both the stream-format header and the movi list must appear.
        let body = &avi[..];
        assert!(
            body.windows(4).any(|w| w == b"strf"),
            "stream format chunk missing"
        );
        assert!(
            body.windows(4).any(|w| w == b"movi"),
            "movi list missing"
        );
        assert!(
            body.windows(4).any(|w| w == b"idx1"),
            "idx1 index missing"
        );
    }

    #[test]
    fn encode_avi_main_header_reports_frame_count_and_size() {
        // Three frames — dwTotalFrames in the avih header must be 3.
        let frames = [
            solid_frame(2, 2, [0; 4]),
            solid_frame(2, 2, [0; 4]),
            solid_frame(2, 2, [0; 4]),
        ];
        let avi = encode_avi(&frames, 24).unwrap();
        // avih payload starts at byte 32 (after "avih"+size at 24..32).
        // dwTotalFrames is field index 4 → byte 32 + 16 = 48.
        let total_frames = u32::from_le_bytes([avi[48], avi[49], avi[50], avi[51]]);
        assert_eq!(total_frames, 3);
        // dwWidth / dwHeight are fields 8 / 9 → bytes 32+32 and 32+36.
        let w = u32::from_le_bytes([avi[64], avi[65], avi[66], avi[67]]);
        let h = u32::from_le_bytes([avi[68], avi[69], avi[70], avi[71]]);
        assert_eq!((w, h), (2, 2));
    }

    #[test]
    fn encode_avi_frame_chunks_are_24bit_dib_sized() {
        // A 3x2 frame: 24-bit rows pad to 4 bytes. 3px*3B = 9 → pads
        // to 12; *2 rows = 24 bytes per frame chunk.
        let frames = [solid_frame(3, 2, [1, 2, 3, 255])];
        let avi = encode_avi(&frames, 30).unwrap();
        // Find the "00db" chunk and read its size.
        let pos = avi
            .windows(4)
            .position(|w| w == b"00db")
            .expect("a 00db frame chunk");
        let size = u32::from_le_bytes([
            avi[pos + 4],
            avi[pos + 5],
            avi[pos + 6],
            avi[pos + 7],
        ]);
        assert_eq!(size, 24, "expected a 24-byte padded DIB frame");
    }

    #[test]
    fn encode_avi_bgr_byte_order_in_frame_data() {
        // A single 1x1 pure-red RGBA pixel must land as BGR bytes
        // 0,0,255 in the frame chunk.
        let frames = [solid_frame(1, 1, [255, 0, 0, 255])];
        let avi = encode_avi(&frames, 30).unwrap();
        let pos = avi
            .windows(4)
            .position(|w| w == b"00db")
            .expect("a 00db chunk");
        // Frame data starts at pos + 8. 1px row = 3 bytes, padded to 4.
        let data = &avi[pos + 8..pos + 8 + 3];
        assert_eq!(data, &[0, 0, 255], "red pixel should be BGR 0,0,255");
    }

    #[test]
    fn write_avi_rejects_non_avi_extension() {
        let frames = [solid_frame(2, 2, [0; 4])];
        let err = write_avi(&frames, 30, &PathBuf::from("clip.mov")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn write_avi_writes_a_playable_file() {
        let frames = [
            solid_frame(4, 4, [200, 100, 50, 255]),
            solid_frame(4, 4, [50, 100, 200, 255]),
        ];
        let tmp = std::env::temp_dir()
            .join(format!("valenx_vid_{}.avi", std::process::id()));
        write_avi(&frames, 30, &tmp).expect("write avi");
        let bytes = std::fs::read(&tmp).expect("read back");
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(&bytes[0..4], b"RIFF");
        assert_eq!(&bytes[8..12], b"AVI ");
    }

    #[test]
    fn ffmpeg_command_has_the_expected_codec_and_inputs() {
        // Verify the argument vector WITHOUT launching ffmpeg.
        let cmd = encode_h264_mp4_command(1920, 1080, 60, &PathBuf::from("out.mp4"));
        // Reads raw video from a stdin pipe.
        assert!(cmd.iter().any(|a| a == "rawvideo"));
        assert!(cmd.iter().any(|a| a == "pipe:0"));
        // Declares the resolution and fps.
        assert!(cmd.iter().any(|a| a == "1920x1080"));
        assert!(cmd.iter().any(|a| a == "60"));
        // Encodes with libx264.
        assert!(cmd.iter().any(|a| a == "libx264"));
        // Writes the requested output path.
        assert!(cmd.iter().any(|a| a == "out.mp4"));
        // Overwrites an existing file.
        assert_eq!(cmd[0], "-y");
    }

    #[test]
    fn ffmpeg_command_uses_rgba_pixel_format() {
        // The frames are RGBA8; the command must tell ffmpeg so.
        let cmd = encode_h264_mp4_command(640, 480, 24, &PathBuf::from("v.mp4"));
        let idx = cmd.iter().position(|a| a == "-pixel_format").unwrap();
        assert_eq!(cmd[idx + 1], "rgba");
    }

    #[test]
    fn run_ffmpeg_mp4_rejects_non_mp4_extension() {
        // This must fail on the extension check BEFORE any spawn.
        let frames = [solid_frame(2, 2, [0; 4])];
        let err = run_ffmpeg_mp4(&frames, 30, &PathBuf::from("clip.avi")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn run_ffmpeg_mp4_rejects_empty_frames_before_spawn() {
        let err = run_ffmpeg_mp4(&[], 30, &PathBuf::from("clip.mp4")).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn video_format_extensions() {
        assert_eq!(VideoFormat::UncompressedAvi.extension(), "avi");
        assert_eq!(VideoFormat::H264Mp4.extension(), "mp4");
        let p = video_output_path(&PathBuf::from("/tmp"), "render", VideoFormat::H264Mp4);
        assert!(p.to_string_lossy().ends_with("render.mp4"));
    }

    // NOTE: run_ffmpeg_mp4 / ffmpeg_available with a real ffmpeg are
    // intentionally NOT tested here — those would spawn a subprocess.
    // UI/subprocess-coupled — run interactively only.
}

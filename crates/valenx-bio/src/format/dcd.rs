//! DCD trajectory reader (NAMD / CHARMM / VMD interchange).
//!
//! DCD is a Fortran-unformatted binary record format: every logical
//! record is wrapped by a 4-byte little-endian length prefix and a
//! matching 4-byte trailer, with the record bytes in between. The
//! Valenx scope target is "produced by recent OpenMM / NAMD" — we
//! handle the common case (little-endian 32-bit float coordinates,
//! no per-frame unit-cell metadata) and stop short of the exotic
//! variants (CHARMM extra blocks, big-endian VAX-era files, 64-bit
//! atom counts, fixed-atom subsets, etc.).
//!
//! Layout we accept (each numbered item is one Fortran record):
//!
//! 1. **Header** (84 bytes content): 4-byte magic `"CORD"`, then
//!    9 × `i32` (nframes, istart, nsavc, nstep, then 4 zero pads,
//!    then ndegf), then `f32` DELTA, then `i32` cell-info flag, then
//!    8 × `i32` zero pads, then `i32` version. We read only the frame
//!    count + cell-info flag and discard the rest.
//! 2. **Titles** (variable): `i32` ntitle followed by ntitle × 80
//!    bytes of ASCII title text. Discarded.
//! 3. **Atom count** (4 bytes): `i32` natom.
//! 4. **Per frame**: optional 48-byte cell record (skipped if the
//!    cell-info flag is non-zero in the header), then three records
//!    of `natom × f32` for X, Y, Z respectively. We materialise each
//!    frame as `Vec<Vector3<f64>>` for the canonical trajectory type.
//!
//! Spec reference:
//! <https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html>

use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use nalgebra::Vector3;
use thiserror::Error;

use crate::trajectory::Trajectory;

/// Errors surfaced by [`read`].
#[derive(Debug, Error)]
pub enum DcdError {
    #[error("dcd I/O error: {0}")]
    Io(#[from] io::Error),
    /// The 4-byte magic at offset 4 of the header record was not the
    /// expected ASCII `"CORD"`.
    #[error("dcd header magic mismatch: found {found:?}, expected {expected:?}")]
    Magic {
        found: [u8; 4],
        expected: &'static str,
    },
    /// The reader hit EOF before consuming a full record. `offset` is
    /// the byte offset where the read attempt started.
    #[error("dcd input truncated at byte offset {offset}")]
    Truncated { offset: u64 },
    /// A per-frame coordinate record carried a different atom count
    /// than the file's declared `natom`.
    #[error(
        "dcd atom-count drift: file declares {first} atoms; frame {frame} \
         carries {found}"
    )]
    InconsistentAtomCount {
        first: usize,
        frame: usize,
        found: usize,
    },
    /// A length read from the DCD header / record prefix exceeded the
    /// adapter's hard cap. The cap exists to defend against malformed
    /// or hostile inputs that declare e.g. `nframes = i32::MAX` and
    /// trick the reader into a multi-gigabyte `Vec::with_capacity`
    /// allocation before any real data has been read.
    #[error("dcd input too large: {what} = {got}, max = {max}")]
    TooLarge {
        what: &'static str,
        got: u64,
        max: u64,
    },
}

/// Maximum frame count we'll honour from a DCD header. Real-world
/// trajectories sit in the 1e3–1e6 frame range; 1e7 is a generous
/// ceiling that's still small enough that an attacker can't force a
/// huge `Vec::with_capacity` allocation before we've read a single
/// coordinate.
const MAX_DCD_FRAMES: usize = 10_000_000;

/// Maximum single Fortran record length we'll honour. 100 MiB is large
/// enough for any plausible atom-count × f32 frame record (25 million
/// atoms × 4 bytes) and small enough to bound the [`read_record`]
/// allocator's worst case.
const MAX_DCD_RECORD_LEN: u32 = 100 * 1024 * 1024;

/// Round-4 DoS hardening: cap natom at 25 million (about 100 MiB of
/// f64 coordinate data per axis) to match `MAX_DCD_RECORD_LEN`. The
/// existing record-length cap protects per-frame `read_f32_vec` from
/// allocating a huge slice for one frame, but the `natom` value also
/// feeds `Vec::with_capacity(natom)` directly — bound it here so a
/// malicious DCD setting `natom = i32::MAX` doesn't OOM the host
/// before the per-record cap fires.
pub const MAX_NATOM: usize = 25_000_000;

/// Parse a DCD byte stream into a [`Trajectory`].
///
/// `id` becomes the trajectory's identifier — typically the source
/// filename's stem. The returned trajectory is pre-validated: every
/// frame is guaranteed to carry the same atom count as the file's
/// declared `natom` header.
pub fn read(bytes: &[u8], id: impl Into<String>) -> Result<Trajectory, DcdError> {
    let mut cursor = io::Cursor::new(bytes);

    // Record 1: header. Fortran wraps with 4-byte len prefix +
    // record bytes + 4-byte len trailer; the content for a v1 DCD
    // header is exactly 84 bytes.
    let header = read_record(&mut cursor)?;
    if header.len() < 84 {
        return Err(DcdError::Truncated {
            offset: cursor.position(),
        });
    }
    if &header[0..4] != b"CORD" {
        let mut found = [0u8; 4];
        found.copy_from_slice(&header[0..4]);
        return Err(DcdError::Magic {
            found,
            expected: "CORD",
        });
    }
    // Header layout after magic: 9 × i32, 1 × f32 (DELTA), 1 × i32
    // (cell-info flag), 8 × i32 zero pads, 1 × i32 version. We pull
    // nframes from offset 4 and the cell-info flag from offset 48 (=
    // 4 magic + 9*4 + 4 delta), and otherwise treat the body as
    // opaque.
    let nframes_signed = i32_le(&header[4..8])?;
    if nframes_signed < 0 {
        return Err(DcdError::TooLarge {
            what: "nframes",
            got: nframes_signed as i64 as u64,
            max: MAX_DCD_FRAMES as u64,
        });
    }
    let nframes = nframes_signed as usize;
    if nframes > MAX_DCD_FRAMES {
        return Err(DcdError::TooLarge {
            what: "nframes",
            got: nframes as u64,
            max: MAX_DCD_FRAMES as u64,
        });
    }
    let has_cell = i32_le(&header[48..52])? != 0;

    // Record 2: titles. Discarded — we don't surface them in the
    // canonical Trajectory type yet.
    let _titles = read_record(&mut cursor)?;

    // Record 3: atom count. Exactly 4 bytes of i32.
    let natom_rec = read_record(&mut cursor)?;
    if natom_rec.len() < 4 {
        return Err(DcdError::Truncated {
            offset: cursor.position(),
        });
    }
    // Round-4: reject negative natom before the i32 -> usize cast
    // silently wraps to a huge positive value. Without the check, a
    // malicious DCD setting natom = -1 would emerge as ~18 quintillion
    // and downstream `Vec::with_capacity(natom)` would either OOM the
    // host or panic with allocation-too-large.
    let natom_signed = i32_le(&natom_rec[0..4])?;
    if natom_signed < 0 {
        return Err(DcdError::TooLarge {
            what: "natom",
            got: natom_signed as i64 as u64,
            max: MAX_NATOM as u64,
        });
    }
    let natom = natom_signed as usize;
    if natom > MAX_NATOM {
        return Err(DcdError::TooLarge {
            what: "natom",
            got: natom as u64,
            max: MAX_NATOM as u64,
        });
    }

    // Frame records. Each frame: optional 48-byte cell record, then
    // three N×f32 coordinate records (X, Y, Z).
    let mut frames: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(nframes);
    for frame_idx in 0..nframes {
        if has_cell {
            // Cell record: 6 × f64 = 48 bytes for a-, b-, c-axes +
            // alpha-, beta-, gamma-cos. We discard it.
            let _cell = read_record(&mut cursor)?;
        }
        let xs = read_f32_vec(&mut cursor)?;
        let ys = read_f32_vec(&mut cursor)?;
        let zs = read_f32_vec(&mut cursor)?;
        if xs.len() != natom || ys.len() != natom || zs.len() != natom {
            return Err(DcdError::InconsistentAtomCount {
                first: natom,
                frame: frame_idx,
                found: xs.len().min(ys.len()).min(zs.len()),
            });
        }
        let mut frame = Vec::with_capacity(natom);
        for i in 0..natom {
            frame.push(Vector3::new(xs[i] as f64, ys[i] as f64, zs[i] as f64));
        }
        frames.push(frame);
    }

    Ok(Trajectory::new(id, frames))
}

/// Read one Fortran-unformatted record from `cursor`, returning the
/// inner bytes (the leading + trailing length wrappers are consumed
/// but discarded).
fn read_record<R: io::Read + io::Seek>(cursor: &mut R) -> Result<Vec<u8>, DcdError> {
    let start = cursor.stream_position()?;
    let leading = cursor.read_u32::<LittleEndian>().map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            DcdError::Truncated { offset: start }
        } else {
            DcdError::Io(e)
        }
    })?;
    // Cap the leading-length prefix before allocating. Otherwise a
    // hostile DCD could declare a 4 GB record and we'd issue a single
    // `vec![0u8; 4_294_967_295]` allocation before reading any data.
    if leading > MAX_DCD_RECORD_LEN {
        return Err(DcdError::TooLarge {
            what: "record length",
            got: leading as u64,
            max: MAX_DCD_RECORD_LEN as u64,
        });
    }
    let mut buf = vec![0u8; leading as usize];
    cursor.read_exact(&mut buf).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            DcdError::Truncated { offset: start }
        } else {
            DcdError::Io(e)
        }
    })?;
    let trailing = cursor.read_u32::<LittleEndian>().map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            DcdError::Truncated { offset: start }
        } else {
            DcdError::Io(e)
        }
    })?;
    if trailing != leading {
        // Length prefix and trailer didn't match — file is corrupt
        // or the byte order isn't little-endian. Treat as truncation
        // (the alternative is a dedicated `DcdError::CorruptRecord`,
        // but every well-formed DCD wrapper matches by construction).
        return Err(DcdError::Truncated { offset: start });
    }
    Ok(buf)
}

/// Read one Fortran-unformatted record and reinterpret its body as a
/// `Vec<f32>`. Errors out on non-multiple-of-4 record bodies.
fn read_f32_vec<R: io::Read + io::Seek>(cursor: &mut R) -> Result<Vec<f32>, DcdError> {
    let bytes = read_record(cursor)?;
    if !bytes.len().is_multiple_of(4) {
        return Err(DcdError::Truncated {
            offset: cursor.stream_position()?,
        });
    }
    let n = bytes.len() / 4;
    let mut out = Vec::with_capacity(n);
    for chunk in bytes.chunks_exact(4) {
        let mut buf = [0u8; 4];
        buf.copy_from_slice(chunk);
        out.push(f32::from_le_bytes(buf));
    }
    Ok(out)
}

/// Pull a little-endian `i32` out of a 4-byte slice. Returns
/// [`DcdError::Truncated`] if the slice is short — callers always
/// pass a window they've already length-checked, so the path is a
/// defensive belt-and-suspenders rather than a hot error path.
fn i32_le(bytes: &[u8]) -> Result<i32, DcdError> {
    if bytes.len() < 4 {
        return Err(DcdError::Truncated { offset: 0 });
    }
    let mut buf = [0u8; 4];
    buf.copy_from_slice(&bytes[0..4]);
    Ok(i32::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthesise a Fortran-wrapped record from raw bytes — i.e.
    /// prepend / append the LE u32 length wrappers.
    fn wrap(body: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(body.len() + 8);
        out.extend_from_slice(&(body.len() as u32).to_le_bytes());
        out.extend_from_slice(body);
        out.extend_from_slice(&(body.len() as u32).to_le_bytes());
        out
    }

    /// Build the 84-byte DCD header body for a given frame count
    /// with no cell metadata. Layout: 4 bytes magic + 9×i32 + 1×f32
    /// (DELTA) + 1×i32 (has_cell flag) + 8×i32 zero pad + 1×i32
    /// version = 4 + 36 + 4 + 4 + 32 + 4 = 84 bytes.
    fn build_header(nframes: u32) -> Vec<u8> {
        let mut h = Vec::with_capacity(84);
        h.extend_from_slice(b"CORD");
        h.extend_from_slice(&(nframes as i32).to_le_bytes()); // nframes
        h.extend_from_slice(&0i32.to_le_bytes()); // istart
        h.extend_from_slice(&1i32.to_le_bytes()); // nsavc
        h.extend_from_slice(&(nframes as i32).to_le_bytes()); // nstep
        for _ in 0..4 {
            h.extend_from_slice(&0i32.to_le_bytes()); // 4 zero pads
        }
        h.extend_from_slice(&0i32.to_le_bytes()); // ndegf
        h.extend_from_slice(&0.001f32.to_le_bytes()); // DELTA
        h.extend_from_slice(&0i32.to_le_bytes()); // cell-info flag (no cell)
        for _ in 0..8 {
            h.extend_from_slice(&0i32.to_le_bytes()); // 8 zero pads
        }
        h.extend_from_slice(&24i32.to_le_bytes()); // version
        debug_assert_eq!(h.len(), 84);
        h
    }

    /// Build the title block: 1 title with placeholder padding.
    fn build_titles() -> Vec<u8> {
        let mut t = Vec::with_capacity(4 + 80);
        t.extend_from_slice(&1i32.to_le_bytes()); // ntitle
        let title = b"valenx synthetic test trajectory";
        let mut padded = [b' '; 80];
        padded[..title.len()].copy_from_slice(title);
        t.extend_from_slice(&padded);
        t
    }

    /// Build a contiguous f32 byte array from a slice.
    fn pack_f32s(values: &[f32]) -> Vec<u8> {
        let mut out = Vec::with_capacity(values.len() * 4);
        for v in values {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }

    #[test]
    fn read_synthesised_two_frame_dcd() {
        // 3-atom, 2-frame trajectory.
        let frame0_x = [0.0f32, 1.0, 2.0];
        let frame0_y = [0.0f32, 0.0, 0.0];
        let frame0_z = [0.0f32, 0.0, 0.0];
        let frame1_x = [0.1f32, 1.1, 2.1];
        let frame1_y = [0.2f32, 0.2, 0.2];
        let frame1_z = [0.3f32, 0.3, 0.3];

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&build_header(2)));
        bytes.extend_from_slice(&wrap(&build_titles()));
        bytes.extend_from_slice(&wrap(&3i32.to_le_bytes()));
        // Frame 0
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame0_x)));
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame0_y)));
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame0_z)));
        // Frame 1
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame1_x)));
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame1_y)));
        bytes.extend_from_slice(&wrap(&pack_f32s(&frame1_z)));

        let traj = read(&bytes, "synth").expect("synthetic DCD parses");
        assert_eq!(traj.id, "synth");
        assert_eq!(traj.frame_count(), 2);
        assert_eq!(traj.atom_count(), Some(3));

        let f0 = traj.frame(0).expect("frame 0");
        assert_eq!(f0[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(f0[1], Vector3::new(1.0, 0.0, 0.0));
        assert_eq!(f0[2], Vector3::new(2.0, 0.0, 0.0));

        let f1 = traj.frame(1).expect("frame 1");
        // Compare with epsilon — single-precision round-trip leaks
        // a few ulps of imprecision into the f64 promotion.
        let eps = 1e-5;
        assert!((f1[0].x - 0.1).abs() < eps);
        assert!((f1[0].y - 0.2).abs() < eps);
        assert!((f1[0].z - 0.3).abs() < eps);
        assert!((f1[2].x - 2.1).abs() < eps);
    }

    #[test]
    fn read_rejects_non_cord_magic() {
        // Build a header but corrupt the magic.
        let mut header = build_header(0);
        header[0..4].copy_from_slice(b"BAAD");

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&header));
        bytes.extend_from_slice(&wrap(&build_titles()));
        bytes.extend_from_slice(&wrap(&0i32.to_le_bytes()));

        let err = read(&bytes, "bad").expect_err("magic mismatch");
        match err {
            DcdError::Magic { found, expected } => {
                assert_eq!(&found, b"BAAD");
                assert_eq!(expected, "CORD");
            }
            other => panic!("expected Magic, got {other:?}"),
        }
    }

    #[test]
    fn read_handles_zero_frames() {
        // Header declares 0 frames; titles + atom-count records still
        // present. Reader must yield an empty `frames` vec without
        // attempting any per-frame reads.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&build_header(0)));
        bytes.extend_from_slice(&wrap(&build_titles()));
        bytes.extend_from_slice(&wrap(&5i32.to_le_bytes()));

        let traj = read(&bytes, "empty").expect("zero-frame DCD parses");
        assert_eq!(traj.frame_count(), 0);
        // atom_count() reads from frame 0 — None when no frames.
        assert_eq!(traj.atom_count(), None);
        // But the validate contract still holds (vacuously).
        assert!(traj.validate().is_ok());
    }

    /// Round-3 fix: hostile DCD with `nframes = i32::MAX` would
    /// previously trigger a multi-gigabyte `Vec::with_capacity`
    /// allocation before any payload data was read.
    #[test]
    fn read_rejects_unbounded_frame_count() {
        let header = build_header(i32::MAX as u32);
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&header));
        bytes.extend_from_slice(&wrap(&build_titles()));
        bytes.extend_from_slice(&wrap(&0i32.to_le_bytes()));

        let err = read(&bytes, "hostile").expect_err("frame cap rejects huge nframes");
        match err {
            DcdError::TooLarge { what, got, max } => {
                assert_eq!(what, "nframes");
                assert!(got > max, "expected got > max; got={got} max={max}");
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    /// Round-4 fix: hostile DCD with `natom = -1` would silently wrap
    /// to `~18 quintillion` after the `as usize` cast, then
    /// `Vec::with_capacity(natom)` would OOM the host. The signed
    /// check must reject before the cast.
    #[test]
    fn read_rejects_negative_natom() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&build_header(1)));
        bytes.extend_from_slice(&wrap(&build_titles()));
        // natom = -1 — the canonical attack value. The reader's i32
        // cast would otherwise re-emerge as 0xFFFFFFFF = 4294967295
        // which sails past MAX_NATOM unless we check the sign first.
        bytes.extend_from_slice(&wrap(&(-1i32).to_le_bytes()));

        let err = read(&bytes, "negative-natom").expect_err("natom cap rejects negative");
        match err {
            DcdError::TooLarge { what, max, .. } => {
                assert_eq!(what, "natom");
                assert_eq!(max, MAX_NATOM as u64);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    /// Round-4 fix: hostile DCD with `natom = i32::MAX` must also be
    /// rejected — the positive-overflow path against `MAX_NATOM`.
    #[test]
    fn read_rejects_natom_above_cap() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&build_header(1)));
        bytes.extend_from_slice(&wrap(&build_titles()));
        bytes.extend_from_slice(&wrap(&i32::MAX.to_le_bytes()));

        let err = read(&bytes, "huge-natom").expect_err("natom cap rejects huge value");
        match err {
            DcdError::TooLarge { what, got, max } => {
                assert_eq!(what, "natom");
                assert_eq!(max, MAX_NATOM as u64);
                assert!(got > max, "expected got > max; got={got} max={max}");
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    /// Round-3 fix: a single Fortran record claiming a 4 GB length
    /// must be rejected before `vec![0u8; len]` allocates anything.
    /// We craft a stream with a header whose leading-length prefix
    /// exceeds the cap.
    #[test]
    fn read_rejects_oversized_record() {
        let mut bytes = Vec::new();
        // Fake a leading-length prefix of MAX_DCD_RECORD_LEN + 1, then
        // garbage. The reader should reject as soon as it sees the
        // prefix.
        let too_big = MAX_DCD_RECORD_LEN as u64 + 1;
        bytes.extend_from_slice(&(too_big as u32).to_le_bytes());
        // Stuff a single zero byte so the cursor doesn't immediately
        // hit EOF on the prefix read; the cap check fires first.
        bytes.push(0);
        let err = read(&bytes, "oversized").expect_err("record cap rejects huge prefix");
        match err {
            DcdError::TooLarge { what, got, .. } => {
                assert_eq!(what, "record length");
                assert_eq!(got, too_big);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }
}

//! GROMACS **TRR** trajectory reader (uncompressed XDR binary).
//!
//! TRR is GROMACS's *lossless* trajectory format: every frame may carry
//! the simulation box, coordinates, velocities and forces, plus the
//! step, time and λ. Unlike [`xtc`](super) (the *compressed* sibling,
//! deferred — see below), TRR stores its reals uncompressed, so it is
//! the clean first binary trajectory reader to implement from scratch.
//!
//! ## Wire format
//!
//! TRR is **XDR**-encoded: every scalar is **big-endian**, integers are
//! 4 bytes, and a "real" is 4 bytes (`f32`) or 8 bytes (`f64`) depending
//! on the precision the file was written in (see below). Each frame is a
//! header followed by the present data blocks. The header, as emitted by
//! GROMACS `do_trnheader` (`src/gromacs/fileio/trrio.cpp`), is:
//!
//! | field        | type   | notes                                      |
//! |--------------|--------|--------------------------------------------|
//! | `magic`      | i32    | **1993** for TRR                           |
//! | `version`    | string | XDR string: i32 length, then padded bytes  |
//! | `ir_size`    | i32    | input-record block size (bytes)            |
//! | `e_size`     | i32    | energy block size                          |
//! | `box_size`   | i32    | box block size (0 ⇒ no box this frame)     |
//! | `vir_size`   | i32    | virial block size                          |
//! | `pres_size`  | i32    | pressure block size                        |
//! | `top_size`   | i32    | topology block size                        |
//! | `sym_size`   | i32    | symmetry block size                        |
//! | `x_size`     | i32    | coordinate block size (0 ⇒ no x)           |
//! | `v_size`     | i32    | velocity block size (0 ⇒ no v)             |
//! | `f_size`     | i32    | force block size (0 ⇒ no f)                |
//! | `natoms`     | i32    | atom count                                 |
//! | `step`       | i32    | MD step index                              |
//! | `nre`        | i32    | number of energy terms                     |
//! | `t`          | real   | simulation time (ps)                       |
//! | `lambda`     | real   | free-energy λ                              |
//!
//! After the header the present blocks follow **in this fixed order**,
//! each present iff its size field is non-zero: `box` (`3×3` reals),
//! `vir` (`3×3`), `pres` (`3×3`), `x` (`natoms×3`), `v` (`natoms×3`),
//! `f` (`natoms×3`). A file is a back-to-back sequence of such frames.
//!
//! ## Precision detection
//!
//! TRR does not store a precision flag; it is *inferred* from a block
//! size. GROMACS uses the box: if `box_size != 0`, the size of one real
//! is `box_size / 9` (the box is `DIM·DIM = 9` reals). If there is no
//! box this frame we fall back to `x_size / (natoms·3)`, then `v`, then
//! `f`. The result must be exactly `4` (single) or `8` (double); any
//! other value is a parse error. This is exactly GROMACS's own rule.
//!
//! ## Scope, assumptions and what is parsed
//!
//! - **Positions only.** [`read_trr`] / [`read_trr_frames`] build a
//!   [`Trajectory`] of coordinate frames, capturing each frame's
//!   **time** and **box** as per-frame metadata
//!   ([`Trajectory::frame_time`] / [`Trajectory::frame_box`]). Velocity
//!   and force blocks are *correctly skipped* (their byte spans are
//!   consumed so framing stays aligned) but not returned; the trajectory
//!   container is coordinate-centric. A frame with `x_size == 0` (a
//!   velocity/force-only frame) is rejected as having no coordinates.
//! - **Box → [`SimBox`].** The `3×3` box is read in the GROMACS row
//!   convention (`box[i]` is lattice vector `i`); the rows become the
//!   columns of the lattice matrix via [`SimBox::triclinic`]. A box
//!   whose vectors are degenerate (e.g. all-zero) is treated as *no
//!   box* for that frame rather than a hard error, since GROMACS will
//!   write a zero box for a non-periodic system.
//! - **Endianness.** Big-endian only, per the XDR spec. Little-endian
//!   TRR files do not occur in practice (GROMACS always writes XDR); a
//!   little-endian file would fail the magic check loudly.
//! - **Bounds-checked.** Every read is length-checked against the
//!   buffer; a truncated or garbage file yields [`MdError::Parse`],
//!   never an out-of-bounds index or panic. Header-declared counts are
//!   capped before any allocation so a hostile `natoms`/`nframes`
//!   cannot drive a giant `Vec` (mirrors the binary DCD-class reader in
//!   [`super::trajectory`]).
//!
//! ### Deferred: XTC (the compressed format)
//!
//! XTC packs coordinates with GROMACS's bespoke `xdr3dfcoord`
//! integer-compression scheme (a per-frame precision factor, smallest
//! bounding-box bit-packing of integerised coordinates). That codec is
//! substantially more involved than the uncompressed XDR reads here and
//! is intentionally **left as a follow-up**; TRR is the clean
//! uncompressed first reader.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::io::trajectory::Trajectory;
use crate::pbc::SimBox;

/// The TRR header magic integer.
const TRR_MAGIC: i32 = 1993;

/// Spatial dimension. TRR boxes are `DIM×DIM`; coordinates `natoms×DIM`.
const DIM: usize = 3;

/// Upper bound on the XDR version string length (bytes). The real
/// string is `"GMX_trn_file"` (13 with the trailing NUL); this generous
/// cap rejects a corrupt length field before it is used to skip bytes.
const MAX_VERSION_LEN: usize = 1024;

/// Hard cap on the per-frame atom count, applied before any allocation.
/// Matches the cap the sibling DCD-class reader uses.
const MAX_TRR_ATOMS: usize = 25_000_000;

/// Hard cap on the number of frames read from one TRR buffer.
const MAX_TRR_FRAMES: usize = 10_000_000;

/// A forward-only, bounds-checked cursor over a TRR byte buffer.
///
/// Every accessor advances the cursor and returns [`MdError::Parse`]
/// rather than panicking when the buffer is too short, so a truncated
/// or hostile file can never trigger an out-of-bounds read.
struct Xdr<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Xdr<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Xdr { bytes, pos: 0 }
    }

    /// Whether the cursor is exactly at end-of-buffer.
    fn at_end(&self) -> bool {
        self.pos >= self.bytes.len()
    }

    /// Returns the next `n` bytes and advances, or a parse error if
    /// fewer than `n` bytes remain.
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or_else(|| MdError::parse("trr", "byte offset overflow"))?;
        if end > self.bytes.len() {
            return Err(MdError::parse(
                "trr",
                format!(
                    "truncated: need {n} bytes at offset {}, only {} remain",
                    self.pos,
                    self.bytes.len().saturating_sub(self.pos)
                ),
            ));
        }
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    /// Reads a big-endian `i32`.
    fn read_i32(&mut self) -> Result<i32> {
        let b = self.take(4)?;
        Ok(i32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// Reads a big-endian `f32`.
    fn read_f32(&mut self) -> Result<f32> {
        let b = self.take(4)?;
        Ok(f32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// Reads a big-endian `f64`.
    fn read_f64(&mut self) -> Result<f64> {
        let b = self.take(8)?;
        Ok(f64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    /// Reads one real (`f32` or `f64`) as `f64`, per the frame's
    /// precision.
    fn read_real(&mut self, prec: Precision) -> Result<f64> {
        match prec {
            Precision::Single => Ok(self.read_f32()? as f64),
            Precision::Double => self.read_f64(),
        }
    }

    /// Skips an XDR string: a 4-byte length followed by that many bytes
    /// padded up to a multiple of 4. The length is bounds-checked so a
    /// corrupt field cannot skip past the buffer.
    fn skip_xdr_string(&mut self) -> Result<()> {
        let len = self.read_i32()?;
        if len < 0 || len as usize > MAX_VERSION_LEN {
            return Err(MdError::parse(
                "trr",
                format!("implausible version-string length {len}"),
            ));
        }
        // XDR pads opaque/string data to a 4-byte boundary.
        let padded = (len as usize).div_ceil(4) * 4;
        self.take(padded)?;
        Ok(())
    }

    /// Skips `n` bytes (used to consume blocks we do not retain, e.g.
    /// virial / pressure / velocity / force).
    fn skip(&mut self, n: usize) -> Result<()> {
        self.take(n)?;
        Ok(())
    }
}

/// Real-number precision of a TRR frame.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Precision {
    /// 4-byte `f32` reals.
    Single,
    /// 8-byte `f64` reals.
    Double,
}

impl Precision {
    /// Bytes per real.
    fn size(self) -> usize {
        match self {
            Precision::Single => 4,
            Precision::Double => 8,
        }
    }

    /// Maps a measured byte-size-of-one-real to a precision.
    fn from_real_size(sz: usize) -> Result<Self> {
        match sz {
            4 => Ok(Precision::Single),
            8 => Ok(Precision::Double),
            other => Err(MdError::parse(
                "trr",
                format!("real size {other} is neither single (4) nor double (8)"),
            )),
        }
    }
}

/// The parsed fixed part of a TRR frame header (block byte-sizes plus
/// the scalar fields), with precision already resolved.
#[derive(Debug)]
struct TrrHeader {
    box_size: usize,
    vir_size: usize,
    pres_size: usize,
    x_size: usize,
    v_size: usize,
    f_size: usize,
    natoms: usize,
    t: f64,
    prec: Precision,
}

/// Parses one TRR frame header from the cursor (magic already implied
/// by the caller having peeked a frame). Returns the resolved header.
fn read_header(xdr: &mut Xdr) -> Result<TrrHeader> {
    let magic = xdr.read_i32()?;
    if magic != TRR_MAGIC {
        return Err(MdError::parse(
            "trr",
            format!("bad magic {magic} (expected {TRR_MAGIC})"),
        ));
    }
    // version string — informational, skipped.
    xdr.skip_xdr_string()?;

    // The 10 block-size fields. Sizes are byte counts and must be >= 0.
    let read_size = |xdr: &mut Xdr, what: &'static str| -> Result<usize> {
        let v = xdr.read_i32()?;
        if v < 0 {
            return Err(MdError::parse(
                "trr",
                format!("negative {what} block size {v}"),
            ));
        }
        Ok(v as usize)
    };
    let _ir_size = read_size(xdr, "ir")?;
    let _e_size = read_size(xdr, "e")?;
    let box_size = read_size(xdr, "box")?;
    let vir_size = read_size(xdr, "vir")?;
    let pres_size = read_size(xdr, "pres")?;
    let _top_size = read_size(xdr, "top")?;
    let _sym_size = read_size(xdr, "sym")?;
    let x_size = read_size(xdr, "x")?;
    let v_size = read_size(xdr, "v")?;
    let f_size = read_size(xdr, "f")?;

    let natoms_i = xdr.read_i32()?;
    if natoms_i <= 0 {
        return Err(MdError::parse(
            "trr",
            format!("non-positive natoms {natoms_i}"),
        ));
    }
    let natoms = natoms_i as usize;
    if natoms > MAX_TRR_ATOMS {
        return Err(MdError::parse(
            "trr",
            format!("natoms {natoms} exceeds the {MAX_TRR_ATOMS} cap"),
        ));
    }
    let _step = xdr.read_i32()?;
    let _nre = xdr.read_i32()?;

    // Resolve precision from a present block BEFORE reading t / lambda,
    // which are themselves reals.
    let prec = resolve_precision(box_size, x_size, v_size, f_size, natoms)?;

    let t = xdr.read_real(prec)?;
    let _lambda = xdr.read_real(prec)?;

    Ok(TrrHeader {
        box_size,
        vir_size,
        pres_size,
        x_size,
        v_size,
        f_size,
        natoms,
        t,
        prec,
    })
}

/// Infers the real-size (and thus precision) from whichever block is
/// present, in GROMACS's priority order: box, then x, v, f.
fn resolve_precision(
    box_size: usize,
    x_size: usize,
    v_size: usize,
    f_size: usize,
    natoms: usize,
) -> Result<Precision> {
    if box_size != 0 {
        // box is DIM*DIM reals.
        if box_size % (DIM * DIM) != 0 {
            return Err(MdError::parse(
                "trr",
                format!("box_size {box_size} is not a multiple of {}", DIM * DIM),
            ));
        }
        return Precision::from_real_size(box_size / (DIM * DIM));
    }
    let per_atom = natoms
        .checked_mul(DIM)
        .ok_or_else(|| MdError::parse("trr", "natoms*DIM overflow"))?;
    if per_atom == 0 {
        return Err(MdError::parse("trr", "cannot infer precision: zero atoms"));
    }
    for sz in [x_size, v_size, f_size] {
        if sz != 0 {
            if sz % per_atom != 0 {
                return Err(MdError::parse(
                    "trr",
                    format!("block size {sz} is not a multiple of natoms*{DIM}"),
                ));
            }
            return Precision::from_real_size(sz / per_atom);
        }
    }
    Err(MdError::parse(
        "trr",
        "frame has no box / x / v / f block — cannot infer precision",
    ))
}

/// Reads the `DIM×DIM` box block into a [`SimBox`], or `None` if the box
/// is degenerate (a zero box for a non-periodic system).
fn read_box(xdr: &mut Xdr, prec: Precision) -> Result<Option<SimBox>> {
    let mut rows = [Vector3::<f64>::zeros(); DIM];
    for row in rows.iter_mut() {
        let x = xdr.read_real(prec)?;
        let y = xdr.read_real(prec)?;
        let z = xdr.read_real(prec)?;
        *row = Vector3::new(x, y, z);
    }
    // GROMACS box[i] is lattice vector i; SimBox columns are the lattice
    // vectors. A degenerate (e.g. all-zero) box means "no periodicity"
    // for this frame, which is valid — surface it as None rather than an
    // error.
    match SimBox::triclinic(rows[0], rows[1], rows[2]) {
        Ok(b) => Ok(Some(b)),
        Err(_) => Ok(None),
    }
}

/// Reads the coordinate block (`natoms` × `DIM` reals) into positions.
fn read_coords(xdr: &mut Xdr, natoms: usize, prec: Precision) -> Result<Vec<Vector3<f64>>> {
    // `natoms` is already capped in the header, so this allocation is
    // bounded; the per-real reads below are individually bounds-checked.
    let mut out = Vec::with_capacity(natoms);
    for _ in 0..natoms {
        let x = xdr.read_real(prec)?;
        let y = xdr.read_real(prec)?;
        let z = xdr.read_real(prec)?;
        out.push(Vector3::new(x, y, z));
    }
    Ok(out)
}

/// Parses **all** frames of a TRR buffer into a [`Trajectory`].
///
/// Each frame contributes its coordinate block as a frame, with the
/// frame's time and box recorded as per-frame metadata. Velocity and
/// force blocks are skipped (their bytes are consumed to keep framing
/// aligned) but not retained.
///
/// The trajectory's nominal frame spacing ([`Trajectory::dt`]) is set
/// from the first two frames' times when at least two frames are
/// present and their times differ; otherwise it defaults to `1.0` ps.
/// Per-frame absolute times are always available via
/// [`Trajectory::frame_time`].
///
/// # Errors
/// [`MdError::Parse`] on a bad magic, a truncated / garbage buffer, an
/// unresolvable precision, a frame without coordinates, an atom-count
/// that drifts between frames, or counts exceeding the safety caps —
/// never a panic or out-of-bounds read.
pub fn read_trr(bytes: &[u8]) -> Result<Trajectory> {
    if bytes.is_empty() {
        return Err(MdError::parse("trr", "empty input"));
    }
    let mut xdr = Xdr::new(bytes);

    // Frame-0 header fixes natoms and seeds the trajectory.
    let h0 = read_header(&mut xdr)?;
    let mut traj = Trajectory::new(h0.natoms, 1.0)?;

    let mut natoms_first = h0.natoms;
    let mut first_two_times: Vec<f64> = Vec::with_capacity(2);
    let mut header = Some(h0);

    let mut nframes = 0usize;
    loop {
        // Read the next header unless we already have frame 0's.
        let h = match header.take() {
            Some(h) => h,
            None => {
                if xdr.at_end() {
                    break;
                }
                let h = read_header(&mut xdr)?;
                if h.natoms != natoms_first {
                    return Err(MdError::parse(
                        "trr",
                        format!(
                            "atom-count drift: frame 0 has {natoms_first} atoms, \
                             a later frame has {}",
                            h.natoms
                        ),
                    ));
                }
                h
            }
        };

        if nframes >= MAX_TRR_FRAMES {
            return Err(MdError::parse(
                "trr",
                format!("exceeded the {MAX_TRR_FRAMES}-frame cap"),
            ));
        }

        // Blocks in fixed order; box / vir / pres precede coordinates.
        let cell = if h.box_size != 0 {
            read_box(&mut xdr, h.prec)?
        } else {
            None
        };
        // vir and pres are full DIM*DIM real blocks we do not retain;
        // skip exactly their declared byte spans.
        xdr.skip(h.vir_size)?;
        xdr.skip(h.pres_size)?;

        if h.x_size == 0 {
            return Err(MdError::parse(
                "trr",
                "frame carries no coordinate (x) block",
            ));
        }
        // Coordinate block size must match natoms * DIM * real-size.
        let expect_x = h
            .natoms
            .checked_mul(DIM)
            .and_then(|n| n.checked_mul(h.prec.size()))
            .ok_or_else(|| MdError::parse("trr", "coordinate block size overflow"))?;
        if h.x_size != expect_x {
            return Err(MdError::parse(
                "trr",
                format!(
                    "x_size {} disagrees with natoms*{DIM}*{} = {expect_x}",
                    h.x_size,
                    h.prec.size()
                ),
            ));
        }
        let coords = read_coords(&mut xdr, h.natoms, h.prec)?;

        // Velocity / force blocks follow; consumed but not retained.
        xdr.skip(h.v_size)?;
        xdr.skip(h.f_size)?;

        if first_two_times.len() < 2 {
            first_two_times.push(h.t);
        }
        traj.push_frame_with(coords, Some(h.t), cell)?;
        nframes += 1;
        natoms_first = h.natoms;
    }

    // Derive a nominal dt from the first two frame times if we can.
    if let [t0, t1] = first_two_times.as_slice() {
        let span = t1 - t0;
        if span.is_finite() && span > 0.0 {
            // Rebuild with the inferred dt; cheap (metadata move).
            traj = with_dt(traj, span)?;
        }
    }

    Ok(traj)
}

/// Returns every TRR frame's coordinates (nm), discarding the
/// [`Trajectory`] wrapper — a convenience mirroring
/// [`super::xyz::read_xyz_frames`].
///
/// # Errors
/// Same conditions as [`read_trr`].
pub fn read_trr_frames(bytes: &[u8]) -> Result<Vec<Vec<Vector3<f64>>>> {
    let traj = read_trr(bytes)?;
    Ok(traj.frames().to_vec())
}

/// Rebuilds a trajectory with a new nominal `dt`, preserving every
/// frame and its per-frame time / box metadata.
fn with_dt(src: Trajectory, dt: f64) -> Result<Trajectory> {
    let mut out = Trajectory::new(src.n_atoms(), dt)?;
    for (i, frame) in src.frames().iter().enumerate() {
        out.push_frame_with(frame.clone(), src.frame_time(i), src.frame_box(i).cloned())?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Appends a big-endian `i32` to `buf`.
    fn put_i32(buf: &mut Vec<u8>, v: i32) {
        buf.extend_from_slice(&v.to_be_bytes());
    }

    /// Appends a big-endian `f32` to `buf`.
    fn put_f32(buf: &mut Vec<u8>, v: f32) {
        buf.extend_from_slice(&v.to_be_bytes());
    }

    /// Appends a big-endian `f64` to `buf`.
    fn put_f64(buf: &mut Vec<u8>, v: f64) {
        buf.extend_from_slice(&v.to_be_bytes());
    }

    /// Appends an XDR string (i32 length incl. NUL, then padded bytes),
    /// as GROMACS writes the version string.
    fn put_xdr_string(buf: &mut Vec<u8>, s: &str) {
        let bytes = s.as_bytes();
        let len = bytes.len() + 1; // include trailing NUL, GROMACS convention
        put_i32(buf, len as i32);
        buf.extend_from_slice(bytes);
        buf.push(0);
        // pad to a 4-byte boundary.
        while buf.len() % 4 != 0 {
            buf.push(0);
        }
    }

    /// Builds one synthetic single-precision TRR frame with a box,
    /// coordinates, and (optionally) velocities, matching the GROMACS
    /// `do_trnheader` layout exactly.
    fn build_trr_frame_single(
        natoms: usize,
        step: i32,
        t: f32,
        box_diag: Option<f32>,
        coords: &[[f32; 3]],
        vels: Option<&[[f32; 3]]>,
    ) -> Vec<u8> {
        assert_eq!(coords.len(), natoms);
        let real = 4usize;
        let box_size = if box_diag.is_some() {
            DIM * DIM * real
        } else {
            0
        };
        let x_size = natoms * DIM * real;
        let v_size = if vels.is_some() {
            natoms * DIM * real
        } else {
            0
        };

        let mut buf = Vec::new();
        put_i32(&mut buf, TRR_MAGIC);
        put_xdr_string(&mut buf, "GMX_trn_file");
        put_i32(&mut buf, 0); // ir_size
        put_i32(&mut buf, 0); // e_size
        put_i32(&mut buf, box_size as i32);
        put_i32(&mut buf, 0); // vir_size
        put_i32(&mut buf, 0); // pres_size
        put_i32(&mut buf, 0); // top_size
        put_i32(&mut buf, 0); // sym_size
        put_i32(&mut buf, x_size as i32);
        put_i32(&mut buf, v_size as i32);
        put_i32(&mut buf, 0); // f_size
        put_i32(&mut buf, natoms as i32);
        put_i32(&mut buf, step);
        put_i32(&mut buf, 0); // nre
        put_f32(&mut buf, t); // time
        put_f32(&mut buf, 0.0); // lambda

        if let Some(d) = box_diag {
            // Diagonal box: rows (d,0,0)(0,d,0)(0,0,d).
            for i in 0..DIM {
                for j in 0..DIM {
                    put_f32(&mut buf, if i == j { d } else { 0.0 });
                }
            }
        }
        for c in coords {
            put_f32(&mut buf, c[0]);
            put_f32(&mut buf, c[1]);
            put_f32(&mut buf, c[2]);
        }
        if let Some(vs) = vels {
            for v in vs {
                put_f32(&mut buf, v[0]);
                put_f32(&mut buf, v[1]);
                put_f32(&mut buf, v[2]);
            }
        }
        buf
    }

    /// Builds one synthetic double-precision TRR frame (box + coords).
    fn build_trr_frame_double(
        natoms: usize,
        step: i32,
        t: f64,
        box_diag: f64,
        coords: &[[f64; 3]],
    ) -> Vec<u8> {
        assert_eq!(coords.len(), natoms);
        let real = 8usize;
        let box_size = DIM * DIM * real;
        let x_size = natoms * DIM * real;

        let mut buf = Vec::new();
        put_i32(&mut buf, TRR_MAGIC);
        put_xdr_string(&mut buf, "GMX_trn_file");
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, box_size as i32);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, x_size as i32);
        put_i32(&mut buf, 0); // v_size
        put_i32(&mut buf, 0); // f_size
        put_i32(&mut buf, natoms as i32);
        put_i32(&mut buf, step);
        put_i32(&mut buf, 0);
        put_f64(&mut buf, t);
        put_f64(&mut buf, 0.0);

        for i in 0..DIM {
            for j in 0..DIM {
                put_f64(&mut buf, if i == j { box_diag } else { 0.0 });
            }
        }
        for c in coords {
            put_f64(&mut buf, c[0]);
            put_f64(&mut buf, c[1]);
            put_f64(&mut buf, c[2]);
        }
        buf
    }

    #[test]
    fn reads_single_precision_frame_exactly() {
        let coords = [[0.0f32, 1.0, 2.0], [3.0, 4.0, 5.0]];
        let bytes = build_trr_frame_single(2, 7, 1.5, Some(4.0), &coords, None);
        let traj = read_trr(&bytes).unwrap();

        assert_eq!(traj.len(), 1);
        assert_eq!(traj.n_atoms(), 2);
        let f = traj.frame(0).unwrap();
        assert!((f[0] - Vector3::new(0.0, 1.0, 2.0)).norm() < 1e-6);
        assert!((f[1] - Vector3::new(3.0, 4.0, 5.0)).norm() < 1e-6);
        // Per-frame time captured.
        assert!((traj.frame_time(0).unwrap() - 1.5).abs() < 1e-6);
        // Box captured (cubic edge 4.0 nm -> volume 64).
        let b = traj.frame_box(0).expect("box present");
        assert!((b.volume() - 64.0).abs() < 1e-4);
    }

    #[test]
    fn reads_double_precision_frame_exactly() {
        let coords = [[0.25f64, -1.5, 2.75], [10.0, 20.0, 30.0]];
        let bytes = build_trr_frame_double(2, 3, 0.002, 5.0, &coords);
        let traj = read_trr(&bytes).unwrap();

        assert_eq!(traj.len(), 1);
        let f = traj.frame(0).unwrap();
        // Double precision is exact for these values.
        assert!((f[0] - Vector3::new(0.25, -1.5, 2.75)).norm() < 1e-12);
        assert!((f[1] - Vector3::new(10.0, 20.0, 30.0)).norm() < 1e-12);
        assert!((traj.frame_time(0).unwrap() - 0.002).abs() < 1e-12);
    }

    #[test]
    fn reads_multi_frame_with_velocities_and_skips_them() {
        // Three frames, each with a box, coordinates AND velocities.
        // The velocities must be skipped without disturbing framing.
        let mut bytes = Vec::new();
        for f in 0i32..3 {
            let c = [[f as f32, 0.0, 0.0], [f as f32 + 0.5, 1.0, 1.0]];
            let v = [[9.0f32, 9.0, 9.0], [8.0, 8.0, 8.0]];
            bytes.extend(build_trr_frame_single(
                2,
                f,
                0.1 * f as f32,
                Some(3.0),
                &c,
                Some(&v),
            ));
        }
        let traj = read_trr(&bytes).unwrap();
        assert_eq!(traj.len(), 3);
        // Frame 2's first atom x == 2.0 confirms framing stayed aligned
        // across the skipped velocity blocks.
        assert!((traj.frame(2).unwrap()[0].x - 2.0).abs() < 1e-6);
        assert!((traj.frame(1).unwrap()[1].x - 1.5).abs() < 1e-6);
        // dt inferred from the first two times (0.1 - 0.0 = 0.1).
        assert!((traj.dt() - 0.1).abs() < 1e-6);
        assert!((traj.frame_time(2).unwrap() - 0.2).abs() < 1e-6);
    }

    #[test]
    fn read_trr_frames_returns_raw_coordinates() {
        let coords = [[1.0f32, 2.0, 3.0]];
        let bytes = build_trr_frame_single(1, 0, 0.0, Some(2.0), &coords, None);
        let frames = read_trr_frames(&bytes).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].len(), 1);
        assert!((frames[0][0] - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-6);
    }

    #[test]
    fn truncated_buffer_fails_loud_no_panic() {
        let coords = [[0.0f32, 1.0, 2.0], [3.0, 4.0, 5.0]];
        let full = build_trr_frame_single(2, 0, 1.0, Some(4.0), &coords, None);
        // Truncate at every prefix length; every one must Err, never panic.
        for cut in 0..full.len() {
            let res = read_trr(&full[..cut]);
            assert!(
                res.is_err(),
                "prefix of len {cut} should fail to parse, not succeed"
            );
        }
        // The full buffer, by contrast, parses.
        assert!(read_trr(&full).is_ok());
    }

    #[test]
    fn bad_magic_is_rejected() {
        let coords = [[0.0f32, 0.0, 0.0]];
        let mut bytes = build_trr_frame_single(1, 0, 0.0, Some(1.0), &coords, None);
        // Corrupt the magic (first 4 bytes).
        bytes[0] = 0xFF;
        bytes[1] = 0xFF;
        let err = read_trr(&bytes).unwrap_err();
        assert_eq!(err.code(), "md.parse");
    }

    #[test]
    fn garbage_input_never_panics() {
        for bad in [
            b"".as_slice(),
            b"\0\0\0\0",
            b"xxxx",
            &[0xFF; 16],
            &[0u8; 1],
            &[0u8; 64],
        ] {
            // Must return (Ok or Err) without panicking / OOB.
            let _ = read_trr(bad);
            let _ = read_trr_frames(bad);
        }
    }

    #[test]
    fn frame_without_coordinates_is_rejected() {
        // A velocity-only frame (x_size = 0) has no coordinates to load.
        let natoms = 2usize;
        let real = 4usize;
        let v_size = natoms * DIM * real;
        let mut buf = Vec::new();
        put_i32(&mut buf, TRR_MAGIC);
        put_xdr_string(&mut buf, "GMX_trn_file");
        for _ in 0..2 {
            put_i32(&mut buf, 0);
        }
        put_i32(&mut buf, 0); // box_size
        put_i32(&mut buf, 0); // vir
        put_i32(&mut buf, 0); // pres
        put_i32(&mut buf, 0); // top
        put_i32(&mut buf, 0); // sym
        put_i32(&mut buf, 0); // x_size = 0
        put_i32(&mut buf, v_size as i32); // v_size present
        put_i32(&mut buf, 0); // f_size
        put_i32(&mut buf, natoms as i32);
        put_i32(&mut buf, 0); // step
        put_i32(&mut buf, 0); // nre
        put_f32(&mut buf, 0.0); // t (precision inferred from v block)
        put_f32(&mut buf, 0.0); // lambda
        for _ in 0..(natoms * DIM) {
            put_f32(&mut buf, 1.0);
        }
        assert!(read_trr(&buf).is_err());
    }

    #[test]
    fn bad_precision_size_is_rejected() {
        // A box_size that implies a 6-byte real (neither 4 nor 8) must
        // be rejected, not silently mis-decoded.
        let mut buf = Vec::new();
        put_i32(&mut buf, TRR_MAGIC);
        put_xdr_string(&mut buf, "GMX_trn_file");
        for _ in 0..2 {
            put_i32(&mut buf, 0);
        }
        put_i32(&mut buf, (DIM * DIM * 6) as i32); // box_size => real size 6
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, (DIM * 4) as i32); // x_size for 1 atom
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 1); // natoms
        put_i32(&mut buf, 0);
        put_i32(&mut buf, 0);
        // (no point appending more — it must fail at precision resolution)
        assert!(read_trr(&buf).is_err());
    }
}

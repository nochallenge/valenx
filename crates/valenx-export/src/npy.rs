//! `.npy` writer (and the companion `.npz` writer that wraps multiple
//! `.npy` arrays into a single zip archive).
//!
//! Both formats are standard NumPy on-disk layouts that pandas /
//! NumPy / PyTorch / TensorFlow can load directly.

use std::path::Path;

use crate::ExportError;

/// Write a 1D `f64` array as a NumPy `.npy` file (version 1.0).
///
/// Format reference:
/// <https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html>
///
/// On-disk layout:
///
/// 1. Magic `\x93NUMPY` (6 bytes).
/// 2. Major / minor version: `0x01 0x00` (2 bytes).
/// 3. Header length as little-endian `u16` (2 bytes).
/// 4. Header: an ASCII Python-dict literal padded with spaces and
///    a trailing `\n`, total length aligned so that the first data
///    byte sits at a 64-byte multiple from the file start.
/// 5. Raw little-endian f64 bytes (`shape.0 × 8` bytes).
///
/// Loaded from Python with one line:
/// ```python
/// import numpy as np
/// arr = np.load("p.npy")     # → np.ndarray, dtype=float64
/// ```
pub fn write_npy_f64(path: &Path, data: &[f64]) -> Result<(), ExportError> {
    write_npy_f64_nd(path, data, &[data.len()])
}

/// N-dimensional generalisation of [`write_npy_f64`]. Same on-disk
/// layout, but the shape tuple in the header reflects whatever
/// `shape` the caller passes in. `data.len()` must equal the product
/// of `shape`; mismatches panic on debug builds and are clamped on
/// release (the writer trusts the caller — this is a low-level
/// primitive).
///
/// Used by the dataset-batch exporter to emit `(n_samples, n_inputs)`
/// arrays that ML loaders expect.
pub fn write_npy_f64_nd(path: &Path, data: &[f64], shape: &[usize]) -> Result<(), ExportError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| ExportError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    // R30: the `.npy` payload is fully materialised in `bytes`; publish
    // it atomically (sidecar → fsync → rename), matching the sibling
    // `write_npz_f64` writer — a torn write can't leave a truncated array
    // an ML loader would mis-parse.
    let bytes = build_npy_bytes_f64_nd(data, shape);
    valenx_core::io_caps::atomic_write_bytes(path, &bytes).map_err(|e| ExportError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    Ok(())
}

/// Build the on-disk byte sequence for a 1D `.npy` v1 file. Pulled
/// out of the writer so unit tests can verify the byte layout without
/// touching the filesystem.
pub fn build_npy_bytes_f64(data: &[f64]) -> Vec<u8> {
    build_npy_bytes_f64_nd(data, &[data.len()])
}

/// N-dimensional version of [`build_npy_bytes_f64`]. Renders the
/// shape tuple as `(d0, d1, ..., dN,)` per NumPy's format spec and
/// expects `data` to be C-order ravelled (last axis varies fastest).
pub fn build_npy_bytes_f64_nd(data: &[f64], shape: &[usize]) -> Vec<u8> {
    let expected: usize = shape.iter().product();
    debug_assert_eq!(
        data.len(),
        expected,
        "shape {shape:?} expects {expected} elements; got {}",
        data.len()
    );
    let mut header = String::new();
    use std::fmt::Write;
    let _ = write!(
        &mut header,
        "{{'descr': '<f8', 'fortran_order': False, 'shape': ("
    );
    // Python tuple form: (N,) for 1D so it parses as a tuple not a
    // bare expression; (d0, d1, …) for ND with no trailing comma.
    if shape.len() == 1 {
        let _ = write!(&mut header, "{},", shape[0]);
    } else {
        for (i, d) in shape.iter().enumerate() {
            if i > 0 {
                let _ = write!(&mut header, ", ");
            }
            let _ = write!(&mut header, "{d}");
        }
    }
    let _ = write!(&mut header, "), }}");
    // Pad with spaces + trailing newline so the total preamble
    // (magic 6 + version 2 + header_len 2 + header bytes) is a
    // multiple of 64. Numpy's official format reference asks for
    // this for alignment-friendly mmap.
    let preamble_fixed = 6 + 2 + 2;
    let mut total = preamble_fixed + header.len() + 1; // +1 for trailing \n
    let pad = (64 - (total % 64)) % 64;
    for _ in 0..pad {
        header.push(' ');
    }
    header.push('\n');
    total = preamble_fixed + header.len();
    debug_assert_eq!(total % 64, 0, "preamble alignment broke");
    let header_len = header.len() as u16;

    let mut out: Vec<u8> = Vec::with_capacity(total + data.len() * 8);
    out.extend_from_slice(b"\x93NUMPY");
    out.push(1); // major version
    out.push(0); // minor version
    out.extend_from_slice(&header_len.to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    for v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

// ---------------------------------------------------------------------------
// .npz writer (zip-store wrapping multiple .npy arrays)
// ---------------------------------------------------------------------------

/// One named array bound for an `.npz` archive. Owns the array name
/// (will become the entry's filename inside the zip) and a reference
/// to the f64 data; ND shape is required so the wrapped `.npy` header
/// reflects the original tensor layout.
pub struct NpzEntry<'a> {
    pub name: String,
    pub data: &'a [f64],
    pub shape: Vec<usize>,
}

/// Write multiple f64 arrays into a single `.npz` archive (the
/// zip-of-`.npy` format that NumPy's `np.savez` produces). Loaded
/// from Python with one line:
///
/// ```python
/// import numpy as np
/// arrs = np.load("dataset.npz")
/// inputs = arrs["inputs"]
/// outputs = arrs["outputs"]
/// ```
///
/// Uses STORE compression (no DEFLATE) so we don't pull in `flate2`
/// for the export crate's hot path. Files are larger but loading is
/// faster + the archive structure is simpler to validate. NumPy's
/// `savez_compressed` exists separately for the DEFLATE variant.
pub fn write_npz_f64(path: &Path, entries: &[NpzEntry<'_>]) -> Result<(), ExportError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| ExportError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    let bytes = build_npz_bytes(entries);
    valenx_core::io_caps::atomic_write_bytes(path, &bytes).map_err(|e| ExportError::Io {
        path: path.to_path_buf(),
        source: e,
    })
}

/// Build the on-disk byte sequence for a `.npz` archive. Pulled out
/// of the writer so unit tests can verify the zip layout without
/// touching the filesystem.
pub fn build_npz_bytes(entries: &[NpzEntry<'_>]) -> Vec<u8> {
    // Pre-build each entry's `.npy` payload + CRC, then assemble
    // local file headers + central directory in two passes.
    struct PreparedEntry {
        name_bytes: Vec<u8>,
        npy_bytes: Vec<u8>,
        crc: u32,
    }

    let prepared: Vec<PreparedEntry> = entries
        .iter()
        .map(|e| {
            // Each archive entry's filename is `<name>.npy` per
            // NumPy convention (np.load uses the bare name as the
            // dict key after stripping `.npy`).
            let name = format!("{}.npy", e.name);
            let npy_bytes = build_npy_bytes_f64_nd(e.data, &e.shape);
            let crc = crc32_ieee(&npy_bytes);
            PreparedEntry {
                name_bytes: name.into_bytes(),
                npy_bytes,
                crc,
            }
        })
        .collect();

    let mut out: Vec<u8> = Vec::new();
    let mut local_header_offsets: Vec<u32> = Vec::with_capacity(prepared.len());

    // Pass 1: local file headers + raw data.
    for entry in &prepared {
        local_header_offsets.push(out.len() as u32);
        write_local_file_header(&mut out, entry.crc, &entry.npy_bytes, &entry.name_bytes);
        out.extend_from_slice(&entry.npy_bytes);
    }

    // Pass 2: central directory.
    let cd_start = out.len() as u32;
    for (entry, &offset) in prepared.iter().zip(local_header_offsets.iter()) {
        write_central_directory_header(
            &mut out,
            entry.crc,
            &entry.npy_bytes,
            &entry.name_bytes,
            offset,
        );
    }
    let cd_size = out.len() as u32 - cd_start;

    // EOCD record.
    write_end_of_central_directory(&mut out, prepared.len() as u16, cd_size, cd_start);
    out
}

fn write_local_file_header(out: &mut Vec<u8>, crc: u32, payload: &[u8], name_bytes: &[u8]) {
    let size = payload.len() as u32;
    out.extend_from_slice(&0x04034b50u32.to_le_bytes()); // signature
    out.extend_from_slice(&20u16.to_le_bytes()); // version needed
    out.extend_from_slice(&0u16.to_le_bytes()); // flags
    out.extend_from_slice(&0u16.to_le_bytes()); // compression: STORE
    out.extend_from_slice(&0u16.to_le_bytes()); // mod time
    out.extend_from_slice(&0x21u16.to_le_bytes()); // mod date (1980-01-01)
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&size.to_le_bytes()); // compressed size = uncompressed
    out.extend_from_slice(&size.to_le_bytes()); // uncompressed size
    out.extend_from_slice(&(name_bytes.len() as u16).to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // extra field length
    out.extend_from_slice(name_bytes);
}

fn write_central_directory_header(
    out: &mut Vec<u8>,
    crc: u32,
    payload: &[u8],
    name_bytes: &[u8],
    local_header_offset: u32,
) {
    let size = payload.len() as u32;
    out.extend_from_slice(&0x02014b50u32.to_le_bytes()); // signature
    out.extend_from_slice(&20u16.to_le_bytes()); // version made by
    out.extend_from_slice(&20u16.to_le_bytes()); // version needed
    out.extend_from_slice(&0u16.to_le_bytes()); // flags
    out.extend_from_slice(&0u16.to_le_bytes()); // compression: STORE
    out.extend_from_slice(&0u16.to_le_bytes()); // mod time
    out.extend_from_slice(&0x21u16.to_le_bytes()); // mod date
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&size.to_le_bytes()); // compressed size
    out.extend_from_slice(&size.to_le_bytes()); // uncompressed size
    out.extend_from_slice(&(name_bytes.len() as u16).to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // extra field length
    out.extend_from_slice(&0u16.to_le_bytes()); // file comment length
    out.extend_from_slice(&0u16.to_le_bytes()); // disk number start
    out.extend_from_slice(&0u16.to_le_bytes()); // internal attrs
    out.extend_from_slice(&0u32.to_le_bytes()); // external attrs
    out.extend_from_slice(&local_header_offset.to_le_bytes());
    out.extend_from_slice(name_bytes);
}

fn write_end_of_central_directory(out: &mut Vec<u8>, n_entries: u16, cd_size: u32, cd_offset: u32) {
    out.extend_from_slice(&0x06054b50u32.to_le_bytes()); // signature
    out.extend_from_slice(&0u16.to_le_bytes()); // disk number
    out.extend_from_slice(&0u16.to_le_bytes()); // disk with CD start
    out.extend_from_slice(&n_entries.to_le_bytes()); // entries on this disk
    out.extend_from_slice(&n_entries.to_le_bytes()); // total entries
    out.extend_from_slice(&cd_size.to_le_bytes());
    out.extend_from_slice(&cd_offset.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes()); // comment length
}

/// CRC-32 (IEEE polynomial 0xEDB88320) — required by the ZIP spec.
/// Inlined so we don't take a dep on `crc32fast` for one user.
pub fn crc32_ieee(bytes: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &b in bytes {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ 0xEDB8_8320
            } else {
                crc >> 1
            };
        }
    }
    crc ^ 0xFFFF_FFFF
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn npy_bytes_start_with_numpy_magic_and_v1_version() {
        let bytes = build_npy_bytes_f64(&[1.0, 2.0, 3.0]);
        // Magic + version.
        assert_eq!(&bytes[0..6], b"\x93NUMPY");
        assert_eq!(bytes[6], 1, "major version should be 1");
        assert_eq!(bytes[7], 0, "minor version should be 0");
    }

    #[test]
    fn npy_header_length_field_is_little_endian_u16() {
        let bytes = build_npy_bytes_f64(&[0.0]);
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]);
        // Header text must declare the f8 little-endian dtype +
        // the (1,) shape.
        let header =
            std::str::from_utf8(&bytes[10..10 + header_len as usize]).expect("header is ascii");
        assert!(header.contains("'descr': '<f8'"), "header was: {header:?}");
        assert!(header.contains("'shape': (1,)"), "header was: {header:?}");
        assert!(header.contains("'fortran_order': False"));
    }

    #[test]
    fn npy_preamble_aligns_to_64_bytes_for_mmap_friendliness() {
        // Numpy's format reference asks the data section to start
        // at a 64-byte boundary. Verify the writer pads the header
        // appropriately for several sizes.
        for n in [0, 1, 7, 1000, 1_000_000] {
            let dummy: Vec<f64> = vec![0.0; n];
            let bytes = build_npy_bytes_f64(&dummy);
            let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
            let preamble = 6 + 2 + 2 + header_len;
            assert_eq!(
                preamble % 64,
                0,
                "preamble for n={n} not 64-aligned (got {preamble})"
            );
        }
    }

    #[test]
    fn npy_data_section_is_raw_le_f64_bytes() {
        let data = [1.5_f64, 2.5, -3.75];
        let bytes = build_npy_bytes_f64(&data);
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let data_start = 6 + 2 + 2 + header_len;
        // 3 × 8 = 24 bytes of data follow.
        assert_eq!(bytes.len() - data_start, 24);
        let v0 = f64::from_le_bytes(bytes[data_start..data_start + 8].try_into().unwrap());
        let v2 = f64::from_le_bytes(bytes[data_start + 16..data_start + 24].try_into().unwrap());
        assert_eq!(v0, 1.5);
        assert_eq!(v2, -3.75);
    }

    #[test]
    fn write_npy_f64_creates_a_loadable_file() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-export-npy-{}.npy",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_npy_f64(&tmp, &[42.0, 43.0, 44.0]).expect("write");
        let bytes = std::fs::read(&tmp).expect("read");
        // Round-trip the magic + first value.
        assert_eq!(&bytes[0..6], b"\x93NUMPY");
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let data_start = 6 + 2 + 2 + header_len;
        let first = f64::from_le_bytes(bytes[data_start..data_start + 8].try_into().unwrap());
        assert_eq!(first, 42.0);
        let _ = std::fs::remove_file(&tmp);
    }

    // -----------------------------------------------------------------
    // .npy nd writer
    // -----------------------------------------------------------------

    #[test]
    fn npy_nd_header_renders_the_full_shape_tuple() {
        let bytes = build_npy_bytes_f64_nd(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[2, 3]);
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let header = std::str::from_utf8(&bytes[10..10 + header_len]).expect("header is ascii");
        assert!(header.contains("'shape': (2, 3)"), "got header: {header:?}");
    }

    #[test]
    fn npy_nd_data_section_is_row_major_packed_f64() {
        // 2x3 matrix in row-major order: row 0 then row 1.
        let data = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
        let bytes = build_npy_bytes_f64_nd(&data, &[2, 3]);
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let data_start = 6 + 2 + 2 + header_len;
        // 6 elements x 8 bytes = 48 bytes.
        assert_eq!(bytes.len() - data_start, 48);
        let v0 = f64::from_le_bytes(bytes[data_start..data_start + 8].try_into().unwrap());
        let v5 = f64::from_le_bytes(bytes[data_start + 40..data_start + 48].try_into().unwrap());
        assert_eq!(v0, 10.0);
        assert_eq!(v5, 60.0);
    }

    // -----------------------------------------------------------------
    // .npz writer (zip-of-npy)
    // -----------------------------------------------------------------

    #[test]
    fn crc32_ieee_matches_known_test_vector() {
        // Reference vector: CRC-32 of "123456789" = 0xCBF43926
        // (see RFC 1952, ZIP appnote, etc.).
        assert_eq!(crc32_ieee(b"123456789"), 0xCBF4_3926);
    }

    #[test]
    fn crc32_ieee_of_empty_input_is_zero() {
        // Standard convention: CRC-32 of an empty buffer is 0 because
        // the seed and final XOR cancel out.
        assert_eq!(crc32_ieee(b""), 0);
    }

    #[test]
    fn npz_starts_with_local_file_header_signature() {
        let data: Vec<f64> = vec![1.0, 2.0, 3.0];
        let entries = vec![NpzEntry {
            name: "alpha".into(),
            data: &data,
            shape: vec![3],
        }];
        let bytes = build_npz_bytes(&entries);
        // First 4 bytes = LFH signature 0x04034b50 (little-endian).
        assert_eq!(&bytes[0..4], &0x04034b50u32.to_le_bytes());
    }

    #[test]
    fn npz_ends_with_end_of_central_directory_signature() {
        let data: Vec<f64> = vec![1.0, 2.0, 3.0];
        let entries = vec![NpzEntry {
            name: "alpha".into(),
            data: &data,
            shape: vec![3],
        }];
        let bytes = build_npz_bytes(&entries);
        // Last 22 bytes are the EOCD record (no comment); check
        // signature at the start of that block.
        let eocd_start = bytes.len() - 22;
        assert_eq!(
            &bytes[eocd_start..eocd_start + 4],
            &0x06054b50u32.to_le_bytes()
        );
    }

    #[test]
    fn npz_central_directory_lists_all_entries() {
        // 3 entries -> EOCD reports 3 in both "entries on this disk"
        // and "total entries" fields.
        let d1 = vec![1.0_f64; 4];
        let d2 = vec![2.0_f64; 8];
        let d3 = vec![3.0_f64; 16];
        let entries = vec![
            NpzEntry {
                name: "a".into(),
                data: &d1,
                shape: vec![4],
            },
            NpzEntry {
                name: "b".into(),
                data: &d2,
                shape: vec![2, 4],
            },
            NpzEntry {
                name: "c".into(),
                data: &d3,
                shape: vec![16],
            },
        ];
        let bytes = build_npz_bytes(&entries);
        let eocd_start = bytes.len() - 22;
        let entries_disk =
            u16::from_le_bytes(bytes[eocd_start + 8..eocd_start + 10].try_into().unwrap());
        let entries_total =
            u16::from_le_bytes(bytes[eocd_start + 10..eocd_start + 12].try_into().unwrap());
        assert_eq!(entries_disk, 3);
        assert_eq!(entries_total, 3);
    }

    #[test]
    fn npz_entry_filenames_get_npy_suffix() {
        // np.load uses the bare name (sans .npy) as the key, so each
        // entry inside the zip must end with `.npy`.
        let d = vec![1.0_f64];
        let entries = vec![NpzEntry {
            name: "drag".into(),
            data: &d,
            shape: vec![1],
        }];
        let bytes = build_npz_bytes(&entries);
        // Local file header filename starts at offset 30; its length
        // is at offset 26 (u16 LE).
        let name_len = u16::from_le_bytes(bytes[26..28].try_into().unwrap()) as usize;
        let name = std::str::from_utf8(&bytes[30..30 + name_len]).unwrap();
        assert_eq!(name, "drag.npy");
    }

    #[test]
    fn write_npz_creates_a_loadable_file() {
        let path = std::env::temp_dir().join(format!(
            "valenx-export-npz-{}.npz",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let inputs = vec![10.0_f64, 20.0, 30.0];
        let outputs = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0];
        let entries = vec![
            NpzEntry {
                name: "inputs".into(),
                data: &inputs,
                shape: vec![3, 1],
            },
            NpzEntry {
                name: "outputs".into(),
                data: &outputs,
                shape: vec![3, 2],
            },
        ];
        write_npz_f64(&path, &entries).expect("write");
        let bytes = std::fs::read(&path).expect("read");
        // Smoke check: starts with LFH, ends with EOCD.
        assert_eq!(&bytes[0..4], &0x04034b50u32.to_le_bytes());
        let eocd_start = bytes.len() - 22;
        assert_eq!(
            &bytes[eocd_start..eocd_start + 4],
            &0x06054b50u32.to_le_bytes()
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn npz_local_header_crc_field_matches_payload_crc() {
        // The local file header records the CRC of the npy-payload
        // bytes; computing the CRC of those bytes a second time
        // should match. Catches off-by-one in the LFH layout.
        let data = vec![1.5_f64, 2.5, -3.75];
        let entries = vec![NpzEntry {
            name: "x".into(),
            data: &data,
            shape: vec![3],
        }];
        let npy_bytes = build_npy_bytes_f64_nd(&data, &[3]);
        let expected_crc = crc32_ieee(&npy_bytes);
        let bytes = build_npz_bytes(&entries);
        // CRC field starts at offset 14 in the LFH.
        let crc_in_header = u32::from_le_bytes(bytes[14..18].try_into().unwrap());
        assert_eq!(crc_in_header, expected_crc);
    }
}

//! Sparse volumetric grid — `SparseVolume<T>`.
//!
//! Dense 3-D scalar/value fields (CFD gas clouds, voxel CT/MRI scans,
//! level sets, occupancy volumes) are usually *mostly empty*: a wisp of
//! smoke fills a tiny fraction of its bounding box, a scanned object is
//! surrounded by air. A flat dense `Vec<T>` over the bounding box pays
//! RAM for every cell whether occupied or not — `O(bbox)` — which is
//! ruinous for large, sparse volumes (a 1024³ grid is ~1e9 cells = GBs
//! even at one byte each).
//!
//! `SparseVolume<T>` stores **only the occupied voxels**, OpenVDB-style:
//! a hash map from *block coordinates* to dense fixed-size **leaf
//! blocks**, each carrying a per-voxel occupancy bitmask. Unset voxels
//! read back as the *background* (`None`). Memory is `O(occupied)` at
//! block granularity — and because real volumetric data clusters
//! spatially, the touched blocks track the data, not the bounding box.
//!
//! This is an in-house structure (no external dep). The well-known
//! `vdb-rs` crate is *read/parse-only* — it loads `.vdb` files into
//! memory but, per its own docs, "supports reading the data and nothing
//! more" (no `set`, no in-memory mutation, VDB writing unimplemented),
//! and it drags in a heavy C/codec dependency stack (`blosc-src`,
//! `bytemuck`, `flate2`, `half`, `glam`, …). The valuable, reusable
//! piece for Valenx is exactly the manipulable in-memory sparse grid
//! that `vdb-rs` lacks, so we build it directly.
//!
//! # Example
//! ```
//! use valenx_fields::sparse_volume::SparseVolume;
//!
//! // 1 mm isotropic voxels.
//! let mut vol: SparseVolume<f32> = SparseVolume::new([1.0, 1.0, 1.0]);
//! vol.set(10, 20, 30, 1.5);
//! assert_eq!(vol.get(10, 20, 30), Some(&1.5));
//! assert_eq!(vol.get(0, 0, 0), None); // background
//! assert_eq!(vol.active_count(), 1);
//! ```

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// Edge length (in voxels) of a leaf block: `8` → `8×8×8 = 512` voxels
/// per block. A power of two so block/in-block coordinates are pure
/// bit-twiddling. Matches OpenVDB's default leaf dimension.
const LEAF_DIM: i32 = 8;
/// `log2(LEAF_DIM)` — shift to map a global coordinate to its block.
const LEAF_LOG2: u32 = 3;
/// Mask to extract the in-block offset of a global coordinate.
const LEAF_MASK: i32 = LEAF_DIM - 1;
/// Voxels per leaf block.
const LEAF_SIZE: usize = (LEAF_DIM * LEAF_DIM * LEAF_DIM) as usize;

/// A dense `LEAF_DIM³` block of values plus a per-voxel occupancy
/// bitmask. Only allocated for blocks that contain at least one set
/// voxel, so empty regions cost nothing.
#[derive(Clone, Debug, Serialize, Deserialize)]
struct LeafBlock<T> {
    /// Dense storage, indexed by [`local_index`]. Slots whose `active`
    /// bit is clear are background and their value is meaningless.
    values: Vec<T>,
    /// One bit per voxel: set ⇒ that voxel is occupied (foreground).
    active: [u64; LEAF_SIZE / 64],
    /// Count of set bits in `active`, cached so `active_count` stays
    /// `O(blocks)` rather than `O(voxels)`.
    count: u32,
}

impl<T: Clone + Default> LeafBlock<T> {
    fn new() -> Self {
        Self {
            values: vec![T::default(); LEAF_SIZE],
            active: [0; LEAF_SIZE / 64],
            count: 0,
        }
    }

    #[inline]
    fn is_active(&self, idx: usize) -> bool {
        (self.active[idx >> 6] >> (idx & 63)) & 1 == 1
    }

    /// Set a voxel; returns `true` if this newly activated it (so the
    /// owner can keep the global active count in sync).
    #[inline]
    fn set(&mut self, idx: usize, val: T) -> bool {
        let was = self.is_active(idx);
        self.values[idx] = val;
        if !was {
            self.active[idx >> 6] |= 1u64 << (idx & 63);
            self.count += 1;
        }
        !was
    }

    /// Clear a voxel; returns `true` if it had been active.
    #[inline]
    fn unset(&mut self, idx: usize) -> bool {
        if self.is_active(idx) {
            self.active[idx >> 6] &= !(1u64 << (idx & 63));
            self.count -= 1;
            self.values[idx] = T::default();
            true
        } else {
            false
        }
    }
}

/// Map a global voxel coordinate to its `(block_coord, local_index)`.
///
/// `block_coord` keys the hash map; `local_index` indexes within the
/// leaf's dense array. Uses arithmetic shift / mask so it is correct
/// for negative coordinates too (`>>` on `i32` is arithmetic in Rust).
#[inline]
fn split(i: i32, j: i32, k: i32) -> ([i32; 3], usize) {
    let block = [i >> LEAF_LOG2, j >> LEAF_LOG2, k >> LEAF_LOG2];
    let li = (i & LEAF_MASK) as usize;
    let lj = (j & LEAF_MASK) as usize;
    let lk = (k & LEAF_MASK) as usize;
    let local = (lk * LEAF_DIM as usize + lj) * LEAF_DIM as usize + li;
    (block, local)
}

/// Inverse of the `local` part of [`split`]: the block origin plus a
/// local index gives back the global coordinate.
#[inline]
fn unsplit(block: [i32; 3], local: usize) -> (i32, i32, i32) {
    let li = (local % LEAF_DIM as usize) as i32;
    let lj = ((local / LEAF_DIM as usize) % LEAF_DIM as usize) as i32;
    let lk = (local / (LEAF_DIM as usize * LEAF_DIM as usize)) as i32;
    (
        block[0] * LEAF_DIM + li,
        block[1] * LEAF_DIM + lj,
        block[2] * LEAF_DIM + lk,
    )
}

/// (De)serialize the block map as a sequence of `(coord, block)` pairs.
///
/// `serde_json` requires object keys to be strings; an `[i32; 3]` key
/// is not, so a bare `HashMap` cannot round-trip. Representing it as a
/// list of pairs sidesteps the restriction and works in every format.
mod block_map_serde {
    use super::LeafBlock;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::collections::HashMap;

    pub fn serialize<S, T>(
        blocks: &HashMap<[i32; 3], LeafBlock<T>>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
        T: Serialize,
    {
        let pairs: Vec<(&[i32; 3], &LeafBlock<T>)> = blocks.iter().collect();
        pairs.serialize(serializer)
    }

    pub fn deserialize<'de, D, T>(
        deserializer: D,
    ) -> Result<HashMap<[i32; 3], LeafBlock<T>>, D::Error>
    where
        D: Deserializer<'de>,
        T: Deserialize<'de>,
    {
        let pairs: Vec<([i32; 3], LeafBlock<T>)> = Vec::deserialize(deserializer)?;
        Ok(pairs.into_iter().collect())
    }
}

/// A hierarchical, hashed **sparse volumetric grid** over integer voxel
/// coordinates `(i, j, k) ∈ ℤ³` (any sign). Stores only occupied
/// voxels; unset cells read as the background `None`.
///
/// Memory is `O(occupied)` at leaf-block granularity, *not* `O(bbox)`.
/// `spacing` records the physical voxel size so world-space extents can
/// be recovered, but indexing is purely by integer voxel coordinate.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))]
pub struct SparseVolume<T> {
    /// Occupied leaf blocks, keyed by block coordinate. Absent key ⇒
    /// that whole `LEAF_DIM³` region is background.
    ///
    /// Serialized as a flat list of `(coord, block)` pairs rather than a
    /// JSON object: a `[i32; 3]` is not a string, and JSON map keys must
    /// be strings, so a plain `HashMap` would fail to round-trip through
    /// `serde_json`.
    #[serde(with = "block_map_serde")]
    blocks: HashMap<[i32; 3], LeafBlock<T>>,
    /// Physical size of one voxel, `[dx, dy, dz]` (world units, e.g. m).
    spacing: [f64; 3],
    /// Total occupied voxels across all blocks (cached for `O(1)`).
    active: usize,
    /// Inclusive min corner of the occupied voxel extent, if non-empty.
    min: Option<[i32; 3]>,
    /// Inclusive max corner of the occupied voxel extent, if non-empty.
    max: Option<[i32; 3]>,
}

impl<T: Clone + Default> SparseVolume<T> {
    /// Create an empty volume with the given physical voxel `spacing`
    /// `[dx, dy, dz]` (world units per voxel along each axis).
    pub fn new(spacing: [f64; 3]) -> Self {
        Self {
            blocks: HashMap::new(),
            spacing,
            active: 0,
            min: None,
            max: None,
        }
    }

    /// Physical voxel size `[dx, dy, dz]` this volume was created with.
    #[inline]
    pub fn spacing(&self) -> [f64; 3] {
        self.spacing
    }

    /// Set voxel `(i, j, k)` to `val`, allocating its leaf block on
    /// first touch. Overwriting an already-set voxel keeps the active
    /// count unchanged.
    pub fn set(&mut self, i: i32, j: i32, k: i32, val: T) {
        let (block, local) = split(i, j, k);
        let leaf = self.blocks.entry(block).or_insert_with(LeafBlock::new);
        if leaf.set(local, val) {
            self.active += 1;
            self.extend_bounds(i, j, k);
        }
    }

    /// Read voxel `(i, j, k)`. Returns `None` for background (never-set
    /// or cleared) cells, `Some(&val)` for occupied ones.
    pub fn get(&self, i: i32, j: i32, k: i32) -> Option<&T> {
        let (block, local) = split(i, j, k);
        let leaf = self.blocks.get(&block)?;
        if leaf.is_active(local) {
            Some(&leaf.values[local])
        } else {
            None
        }
    }

    /// Clear voxel `(i, j, k)` back to background. Returns the previous
    /// value if it had been occupied, else `None`.
    ///
    /// Note: the recorded [`bounds`](Self::bounds) are *not* shrunk on
    /// unset (it would cost a full rescan); they stay a conservative
    /// superset until the next [`recompute_bounds`](Self::recompute_bounds)
    /// or [`clear`](Self::clear).
    pub fn unset(&mut self, i: i32, j: i32, k: i32) -> Option<T> {
        let (block, local) = split(i, j, k);
        let leaf = self.blocks.get_mut(&block)?;
        if leaf.is_active(local) {
            let prev = leaf.values[local].clone();
            leaf.unset(local);
            self.active -= 1;
            // Drop the whole block once it empties, reclaiming its RAM.
            if leaf.count == 0 {
                self.blocks.remove(&block);
            }
            Some(prev)
        } else {
            None
        }
    }

    /// Number of occupied voxels. `O(1)`. For a sparse fill this is the
    /// count of cells actually set, *not* the bounding-box volume.
    #[inline]
    pub fn active_count(&self) -> usize {
        self.active
    }

    /// `true` if no voxel is occupied.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.active == 0
    }

    /// Number of allocated leaf blocks — a proxy for the actual memory
    /// footprint (each block is `LEAF_DIM³` slots). Useful for asserting
    /// the structure stays sparse.
    #[inline]
    pub fn block_count(&self) -> usize {
        self.blocks.len()
    }

    /// Inclusive `(min, max)` corners of the occupied voxel extent, or
    /// `None` while empty. The bounds *grow* to cover every `set` voxel;
    /// they are not shrunk by `unset` (see [`unset`](Self::unset)) until
    /// [`recompute_bounds`](Self::recompute_bounds) is called.
    #[inline]
    pub fn bounds(&self) -> Option<([i32; 3], [i32; 3])> {
        match (self.min, self.max) {
            (Some(lo), Some(hi)) => Some((lo, hi)),
            _ => None,
        }
    }

    /// World-space axis-aligned bounding box of the occupied extent,
    /// `(min_corner, max_corner)` in physical units, accounting for
    /// `spacing` and the full extent of the max voxel. `None` if empty.
    pub fn world_bounds(&self) -> Option<([f64; 3], [f64; 3])> {
        let (lo, hi) = self.bounds()?;
        let min = [
            lo[0] as f64 * self.spacing[0],
            lo[1] as f64 * self.spacing[1],
            lo[2] as f64 * self.spacing[2],
        ];
        let max = [
            (hi[0] + 1) as f64 * self.spacing[0],
            (hi[1] + 1) as f64 * self.spacing[1],
            (hi[2] + 1) as f64 * self.spacing[2],
        ];
        Some((min, max))
    }

    /// Remove every voxel and reset bounds to empty, freeing all block
    /// storage.
    pub fn clear(&mut self) {
        self.blocks.clear();
        self.active = 0;
        self.min = None;
        self.max = None;
    }

    /// Iterate over occupied voxels as `((i, j, k), &value)`, in no
    /// particular order. Yields exactly [`active_count`](Self::active_count)
    /// items.
    pub fn iter(&self) -> impl Iterator<Item = ((i32, i32, i32), &T)> {
        self.blocks.iter().flat_map(|(&block, leaf)| {
            (0..LEAF_SIZE).filter_map(move |local| {
                if leaf.is_active(local) {
                    Some((unsplit(block, local), &leaf.values[local]))
                } else {
                    None
                }
            })
        })
    }

    /// Recompute [`bounds`](Self::bounds) tightly from the currently
    /// occupied voxels. `O(occupied)`. Call after a batch of `unset`s if
    /// a tight extent is needed.
    pub fn recompute_bounds(&mut self) {
        let mut min: Option<[i32; 3]> = None;
        let mut max: Option<[i32; 3]> = None;
        for ((i, j, k), _) in self.iter() {
            match &mut min {
                Some(m) => {
                    m[0] = m[0].min(i);
                    m[1] = m[1].min(j);
                    m[2] = m[2].min(k);
                }
                None => min = Some([i, j, k]),
            }
            match &mut max {
                Some(m) => {
                    m[0] = m[0].max(i);
                    m[1] = m[1].max(j);
                    m[2] = m[2].max(k);
                }
                None => max = Some([i, j, k]),
            }
        }
        self.min = min;
        self.max = max;
    }

    /// Grow the cached bounds to include `(i, j, k)`.
    #[inline]
    fn extend_bounds(&mut self, i: i32, j: i32, k: i32) {
        match &mut self.min {
            Some(m) => {
                m[0] = m[0].min(i);
                m[1] = m[1].min(j);
                m[2] = m[2].min(k);
            }
            None => self.min = Some([i, j, k]),
        }
        match &mut self.max {
            Some(m) => {
                m[0] = m[0].max(i);
                m[1] = m[1].max(j);
                m[2] = m[2].max(k);
            }
            None => self.max = Some([i, j, k]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_get_round_trip() {
        let mut v: SparseVolume<f32> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(3, 4, 5, 2.5);
        v.set(-7, 100, -2, -9.0);
        v.set(0, 0, 0, 1.0);
        assert_eq!(v.get(3, 4, 5), Some(&2.5));
        assert_eq!(v.get(-7, 100, -2), Some(&-9.0));
        assert_eq!(v.get(0, 0, 0), Some(&1.0));
        assert_eq!(v.active_count(), 3);
    }

    #[test]
    fn unset_cells_read_as_background_none() {
        let mut v: SparseVolume<i32> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(10, 10, 10, 42);
        // Never-set neighbours are background.
        assert_eq!(v.get(10, 10, 11), None);
        assert_eq!(v.get(9, 10, 10), None);
        assert_eq!(v.get(1000, 0, 0), None);
        // Negative, far-away cell is background.
        assert_eq!(v.get(-1, -1, -1), None);
    }

    #[test]
    fn overwrite_keeps_count() {
        let mut v: SparseVolume<f64> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(2, 2, 2, 1.0);
        v.set(2, 2, 2, 9.0);
        assert_eq!(v.active_count(), 1);
        assert_eq!(v.get(2, 2, 2), Some(&9.0));
    }

    #[test]
    fn unset_clears_to_background_and_returns_prev() {
        let mut v: SparseVolume<i32> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(5, 6, 7, 11);
        assert_eq!(v.unset(5, 6, 7), Some(11));
        assert_eq!(v.get(5, 6, 7), None);
        assert_eq!(v.active_count(), 0);
        // Unsetting a background cell is a no-op returning None.
        assert_eq!(v.unset(5, 6, 7), None);
        assert_eq!(v.unset(0, 0, 0), None);
    }

    #[test]
    fn sparse_fill_stores_about_one_percent() {
        // A large bbox: 200×200×200 = 8e6 cells. A dense array would
        // hold all 8e6; a 1% scattered fill should store ~1% of cells.
        let mut v: SparseVolume<f32> = SparseVolume::new([0.5, 0.5, 0.5]);
        let side: i32 = 200;
        let total = (side as usize).pow(3);
        let target = total / 100; // 1%
                                  // Scatter `target` distinct voxels deterministically across the
                                  // box using a coprime stride so they spread (not all clustered).
        let mut placed = 0usize;
        let mut idx: u64 = 1;
        let modulus = total as u64;
        let stride: u64 = 7_919; // prime, coprime to modulus
        while placed < target {
            let lin = (idx.wrapping_mul(stride)) % modulus;
            let i = (lin % side as u64) as i32;
            let j = ((lin / side as u64) % side as u64) as i32;
            let k = (lin / (side as u64 * side as u64)) as i32;
            if v.get(i, j, k).is_none() {
                v.set(i, j, k, lin as f32);
                placed += 1;
            }
            idx += 1;
        }
        // active_count is EXACTLY the number of distinct cells set.
        assert_eq!(v.active_count(), target);
        // And that is ~1% of the bounding-box volume — far less than
        // dense. (Sanity: it is well under 2% and well over 0.5%.)
        let frac = v.active_count() as f64 / total as f64;
        assert!((0.009..=0.011).contains(&frac), "fraction {frac} not ~1%");
        // Memory proxy: blocks allocated must be far fewer than the
        // dense block count for the full box.
        let dense_blocks = ((side / LEAF_DIM + 1) as usize).pow(3);
        assert!(
            v.block_count() < dense_blocks,
            "blocks {} should be < dense {}",
            v.block_count(),
            dense_blocks
        );
    }

    #[test]
    fn bounds_track_occupied_extent() {
        let mut v: SparseVolume<u8> = SparseVolume::new([1.0, 1.0, 1.0]);
        assert_eq!(v.bounds(), None);
        v.set(-3, 5, 0, 1);
        v.set(10, -2, 7, 1);
        v.set(4, 4, 4, 1);
        let (lo, hi) = v.bounds().unwrap();
        assert_eq!(lo, [-3, -2, 0]);
        assert_eq!(hi, [10, 5, 7]);
    }

    #[test]
    fn world_bounds_use_spacing() {
        let mut v: SparseVolume<f32> = SparseVolume::new([2.0, 0.5, 4.0]);
        v.set(0, 0, 0, 1.0);
        v.set(1, 2, 1, 1.0);
        let (min, max) = v.world_bounds().unwrap();
        assert_eq!(min, [0.0, 0.0, 0.0]);
        // Max voxel extent = (hi+1) * spacing.
        assert_eq!(max, [4.0, 1.5, 8.0]); // (1+1)*2, (2+1)*0.5, (1+1)*4
    }

    #[test]
    fn iter_yields_all_occupied() {
        let mut v: SparseVolume<i32> = SparseVolume::new([1.0, 1.0, 1.0]);
        let cells = [
            ((0, 0, 0), 1),
            ((1, 0, 0), 2),
            ((0, 9, 0), 3), // different block in j
            ((-5, -5, -5), 4),
        ];
        for ((i, j, k), val) in cells {
            v.set(i, j, k, val);
        }
        let mut seen: Vec<((i32, i32, i32), i32)> = v.iter().map(|(c, &val)| (c, val)).collect();
        seen.sort();
        let mut expect: Vec<((i32, i32, i32), i32)> = cells.to_vec();
        expect.sort();
        assert_eq!(seen, expect);
        assert_eq!(v.iter().count(), 4);
    }

    #[test]
    fn clear_empties_everything() {
        let mut v: SparseVolume<f32> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(1, 2, 3, 1.0);
        v.set(4, 5, 6, 2.0);
        assert_eq!(v.active_count(), 2);
        v.clear();
        assert_eq!(v.active_count(), 0);
        assert!(v.is_empty());
        assert_eq!(v.block_count(), 0);
        assert_eq!(v.bounds(), None);
        assert_eq!(v.get(1, 2, 3), None);
    }

    #[test]
    fn emptied_block_is_reclaimed() {
        let mut v: SparseVolume<i32> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(0, 0, 0, 1);
        v.set(1, 1, 1, 2); // same 8³ block
        assert_eq!(v.block_count(), 1);
        v.unset(0, 0, 0);
        assert_eq!(v.block_count(), 1); // block still has one voxel
        v.unset(1, 1, 1);
        assert_eq!(v.block_count(), 0); // now empty → reclaimed
        assert!(v.is_empty());
    }

    #[test]
    fn recompute_bounds_tightens_after_unset() {
        let mut v: SparseVolume<u8> = SparseVolume::new([1.0, 1.0, 1.0]);
        v.set(0, 0, 0, 1);
        v.set(50, 50, 50, 1);
        v.unset(50, 50, 50);
        // Conservative bounds still cover the removed voxel.
        assert_eq!(v.bounds(), Some(([0, 0, 0], [50, 50, 50])));
        v.recompute_bounds();
        assert_eq!(v.bounds(), Some(([0, 0, 0], [0, 0, 0])));
    }

    #[test]
    fn serde_round_trip() {
        let mut v: SparseVolume<f32> = SparseVolume::new([1.0, 2.0, 3.0]);
        v.set(1, 2, 3, 4.5);
        v.set(-1, -2, -3, -6.5);
        let json = serde_json::to_string(&v).unwrap();
        let back: SparseVolume<f32> = serde_json::from_str(&json).unwrap();
        assert_eq!(back.active_count(), 2);
        assert_eq!(back.get(1, 2, 3), Some(&4.5));
        assert_eq!(back.get(-1, -2, -3), Some(&-6.5));
        assert_eq!(back.spacing(), [1.0, 2.0, 3.0]);
        assert_eq!(back.bounds(), Some(([-1, -2, -3], [1, 2, 3])));
    }
}

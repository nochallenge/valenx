//! Incremental mapper: the capstone that ties Stages 1–5 into one
//! Structure-from-Motion pipeline (the **orchestration** stage).
//!
//! Stages 1–5 each solve one sub-problem of SfM in isolation: feature
//! detection/description ([`crate::fast`] / [`crate::descriptor`]), descriptor
//! matching ([`crate::matching`]), two-view verification + the fundamental
//! matrix ([`crate::geometry`]), the calibrated two-view seed pose plus
//! triangulation ([`crate::twoview`]), per-view registration by resectioning
//! ([`crate::pnp`]), and the joint nonlinear refinement ([`crate::bundle`]).
//! This module is the **conductor** that drives them, in the order COLMAP
//! popularized, to turn a set of images and their pairwise verified matches
//! into a single consistent reconstruction — a set of registered camera poses
//! and a sparse 3-D point cloud:
//!
//! ```text
//!   per-image features  +  per-pair verified inlier matches
//!            │
//!            ├─ build multi-view tracks (union-find over the pairwise matches)
//!            │
//!            ├─ pick the best initial pair → recover_pose → seed 2 cameras +
//!            │     triangulated points + their tracks                (Stage 3)
//!            │
//!            ├─ loop:  next-best unregistered view (most 2D–3D corr.)
//!            │           → solve_pnp_ransac → register it            (Stage 4)
//!            │           → triangulate the new tracks it shares with
//!            │             already-registered views                  (Stage 3)
//!            │           → periodically bundle_adjust                (Stage 5)
//!            │
//!            └─ final global bundle_adjust                           (Stage 5)
//!            ▼
//!         Reconstruction { cameras: [Option<CameraPose>], points, tracks }
//! ```
//!
//! ## Track building (union-find)
//!
//! A *track* is one 3-D scene point together with every 2-D observation of it
//! across the views. The pairwise matcher only ever produces **two-view**
//! correspondences: a [`Match`] inside a [`PairMatches`] links keypoint `q` of
//! image `i` to keypoint `t` of image `j`. To get a *multi-view* track we must
//! *chain* these: if feature `(i, q)` matches `(j, t)` in pair `(i, j)` and
//! `(j, t)` matches `(k, s)` in pair `(j, k)`, then `(i, q)`, `(j, t)` and
//! `(k, s)` are all observations of the **same** 3-D point and belong in one
//! track. We compute these connected components with a classic
//! **union-find** (disjoint-set) over the global feature identifiers
//! `(image_idx, keypoint_idx)`: each match unions its two endpoints, and every
//! resulting set becomes a candidate track.
//!
//! **Inconsistent tracks are dropped.** A correct track observes each image at
//! most once. Chaining can nonetheless merge a component that contains *two
//! different keypoints of the same image* (e.g. via a wrong match somewhere in
//! the chain, or two image features that both — incorrectly — got linked to the
//! same world point). Such a component cannot be a single rigid 3-D point, so
//! [`build_tracks`] discards any component with a repeated `image_idx`. This is
//! a deliberately conservative filter: it removes a genuinely contradictory
//! track outright rather than trying to split it, trading a little recall for
//! correctness (see the honesty notes).
//!
//! ## Greedy incremental registration (stated honestly)
//!
//! The mapper is **greedy**, exactly like COLMAP's: at each step it registers
//! the single unregistered view with the most 2D–3D correspondences into the
//! *current* model, and never backtracks on that choice. Greedy next-view
//! selection is the standard, robust heuristic, but it is a heuristic — a
//! different order can in principle yield a slightly different reconstruction,
//! and a poor initial pair can hurt everything downstream. We mitigate the
//! initial-pair risk by scoring all candidate pairs and taking the strongest
//! (most verified inliers, subject to a minimum), but we do **not** implement
//! COLMAP's more elaborate next-view scoring (which also weighs the spatial
//! distribution of the correspondences), nor re-triangulation, track
//! completion, or the redundant-view removal of a full production mapper.
//!
//! ## Bundle-adjustment cadence (stated honestly)
//!
//! Incremental SfM accumulates drift, so we re-run [`bundle_adjust`]
//! periodically — after every [`MapperParams::ba_every`] newly registered views
//! (a *local*, here actually *global-over-the-current-model* BA) — and once more
//! at the very end (the final global BA). This is the usual local-then-global
//! pattern; the honest caveat is that our intermediate BA already optimizes the
//! whole current model (it is not a windowed/local BA over only the most recent
//! cameras), because this crate's [`bundle_adjust`] is a *dense* solver intended
//! for the modest problem sizes here. For large scenes the windowed-local +
//! sparse-global split would be required (a deliberate future extension, noted
//! in the [`crate::bundle`] docs).
//!
//! ## Gauge + scale freedom (stated honestly)
//!
//! The whole reconstruction is expressed in the **first registered camera's**
//! frame (that camera is the bundle-adjustment gauge — see [`crate::bundle`]).
//! Because the pipeline is seeded from a two-view essential-matrix
//! decomposition, the absolute **scale** is unrecoverable from the pixels
//! alone: the entire scene (every camera centre and every point) may be scaled
//! by one common positive factor without changing any image. The reconstruction
//! is therefore correct only **up to a global similarity** (the rigid gauge is
//! pinned by camera 0; the residual scale is free). Fixing absolute scale needs
//! external information (a known baseline, a scale bar, GPS). See the
//! [`crate::twoview`] and [`crate::bundle`] module docs.
//!
//! ## Scope / failure modes (stated honestly)
//!
//! - Returns [`None`] when no adequate initial pair exists (fewer than
//!   [`MapperParams::min_initial_inliers`] inliers in every pair), when the
//!   seed pose cannot be recovered, or when fewer than two images / no matches
//!   are supplied. A single image with no pairs yields [`None`] (there is
//!   nothing to triangulate against).
//! - All-collinear (or coplanar) scene points are degenerate for the linear
//!   DLT-PnP used to register later views (see [`crate::pnp`]); such views
//!   simply fail to register, and the mapper stops when no view can be added,
//!   returning whatever consistent partial model it built (possibly just the
//!   initial pair).
//! - It never panics on degenerate input.

use crate::bundle::{bundle_adjust, BundleParams, BundleProblem, Observation};
use crate::matching::Match;
use crate::pnp::{solve_pnp_ransac, CameraPose, PnpRansacParams};
use crate::twoview::{
    essential_from_fundamental, recover_pose, triangulate_point, CameraIntrinsics,
};
use crate::Keypoint;
use nalgebra::{Matrix3, Matrix3x4, Vector3};

/// The detected features (keypoints) of a single image.
///
/// This is the per-image input to [`incremental_sfm`]: index `i` of the
/// `images` slice is image `i`, and `keypoints[k]` is its `k`-th feature. The
/// pixel coordinates of these keypoints are what every later stage consumes;
/// the descriptors used to *produce* the matches are not needed here (matching
/// already happened — see [`PairMatches`]).
#[derive(Debug, Clone)]
pub struct ImageFeatures {
    /// The image's detected keypoints, in a fixed order. A [`Match`]'s
    /// `query_idx` / `train_idx` index into this slice for the relevant image.
    pub keypoints: Vec<Keypoint>,
}

impl ImageFeatures {
    /// Wrap a keypoint list as one image's features.
    #[inline]
    #[must_use]
    pub fn new(keypoints: Vec<Keypoint>) -> Self {
        Self { keypoints }
    }
}

/// The geometrically **verified** inlier matches between one ordered image
/// pair `(i, j)`.
///
/// `i` and `j` are indices into the `images` slice passed to
/// [`incremental_sfm`]. Each [`Match`] in `matches` pairs
/// `images[i].keypoints[m.query_idx]` with `images[j].keypoints[m.train_idx]`
/// — i.e. `query_idx` belongs to image `i` and `train_idx` to image `j`. These
/// are expected to be the *inliers* of the two-view verification
/// ([`crate::TwoViewResult::inliers`]), not raw putative matches.
#[derive(Debug, Clone)]
pub struct PairMatches {
    /// Index of the first image of the pair (the `query` side of each match).
    pub i: usize,
    /// Index of the second image of the pair (the `train` side of each match).
    pub j: usize,
    /// The verified inlier correspondences between image `i` and image `j`.
    pub matches: Vec<Match>,
}

/// One multi-view **track**: a single 3-D scene point and the 2-D observations
/// of it across the views.
///
/// A track is the unit the bundle adjuster optimizes a point over: each
/// `(image_idx, keypoint_idx)` in [`Self::observations`] says "feature
/// `keypoint_idx` of image `image_idx` is a projection of this point." A valid
/// track observes any given image **at most once** (the invariant
/// [`build_tracks`] enforces by dropping inconsistent components).
#[derive(Debug, Clone)]
pub struct Track {
    /// The 2-D observations forming the track, as `(image_idx, keypoint_idx)`
    /// pairs. Every `image_idx` is distinct within a single track.
    pub observations: Vec<(usize, usize)>,
    /// The triangulated 3-D point for this track once it has been triangulated
    /// from at least two registered views, in the first registered camera's
    /// frame; [`None`] until then. When `Some(p)`, `p == points[point_idx]` for
    /// the [`Self::point_idx`] this track was assigned.
    pub point3d: Option<Vector3<f64>>,
    /// Index into [`Reconstruction::points`] of this track's 3-D point once
    /// triangulated; [`None`] until then.
    pub point_idx: Option<usize>,
}

impl Track {
    /// Look up this track's observed keypoint index in `image_idx`, if the
    /// track observes that image.
    #[must_use]
    fn keypoint_in(&self, image_idx: usize) -> Option<usize> {
        self.observations
            .iter()
            .find(|(im, _)| *im == image_idx)
            .map(|(_, kp)| *kp)
    }
}

/// A sparse Structure-from-Motion reconstruction: the per-image camera poses
/// (registered views are `Some`), the triangulated 3-D points, and the tracks
/// tying points to their 2-D observations.
///
/// The geometry is expressed in the **first registered camera's** frame and is
/// determined only up to a global scale (see the [module docs](self)).
#[derive(Debug, Clone)]
pub struct Reconstruction {
    /// Per-image camera pose, indexed identically to the input `images` slice.
    /// `cameras[i]` is `Some(pose)` once image `i` has been registered, and
    /// `None` for an image that was never registered (e.g. it shared too few
    /// 2D–3D correspondences with the model, or its points were degenerate).
    pub cameras: Vec<Option<CameraPose>>,
    /// The triangulated 3-D scene points, in the first registered camera's
    /// frame. Index `p` is referenced by the [`Track`] whose `point_idx` is
    /// `Some(p)`.
    pub points: Vec<Vector3<f64>>,
    /// The multi-view tracks. A track with `point_idx == Some(p)` has been
    /// triangulated into `points[p]`; an un-triangulated track (too few
    /// registered views observe it, or its rays were degenerate) keeps
    /// `point_idx == None` and is not part of the point cloud.
    pub tracks: Vec<Track>,
}

impl Reconstruction {
    /// Number of registered (posed) cameras.
    #[must_use]
    pub fn num_registered(&self) -> usize {
        self.cameras.iter().filter(|c| c.is_some()).count()
    }
}

/// Tuning parameters for [`incremental_sfm`].
#[derive(Debug, Clone, Copy)]
pub struct MapperParams {
    /// Minimum number of verified inliers a pair must have to be eligible as
    /// the **initial** pair. Pairs below this are never used to seed the
    /// reconstruction; if no pair clears it, [`incremental_sfm`] returns
    /// [`None`]. The initial pair is the strongest pair (most inliers) among
    /// those that clear this bar.
    pub min_initial_inliers: usize,
    /// Minimum number of 2D–3D correspondences an unregistered view must share
    /// with the current model to be eligible for PnP registration. Must be at
    /// least 6 (the linear DLT-PnP minimum); values below 6 are raised to 6
    /// internally.
    pub min_pnp_correspondences: usize,
    /// Run a (global-over-the-current-model) [`bundle_adjust`] after every this
    /// many newly registered views. `0` disables intermediate BA (only the
    /// final global BA runs). The seed pair counts as the first registration
    /// event for this cadence.
    pub ba_every: usize,
    /// RANSAC parameters for the per-view PnP registration ([`solve_pnp_ransac`]).
    pub pnp: PnpRansacParams,
    /// Parameters for every [`bundle_adjust`] call (intermediate and final).
    pub bundle: BundleParams,
}

impl Default for MapperParams {
    /// Sensible defaults: a 30-inlier initial-pair floor, a 6-correspondence
    /// PnP floor, bundle adjustment every 3 registered views, and the default
    /// PnP-RANSAC / bundle parameters.
    fn default() -> Self {
        Self {
            min_initial_inliers: 30,
            min_pnp_correspondences: 6,
            ba_every: 3,
            pnp: PnpRansacParams::default(),
            bundle: BundleParams::default(),
        }
    }
}

/// Minimum correspondences for the linear DLT-PnP solver (mirrors the private
/// constant in [`crate::pnp`]); a view needs at least this many 2D–3D matches
/// to be registered.
const MIN_PNP: usize = 6;

/// The 2D–3D correspondences an unregistered view shares with the current
/// model, returned aligned element-wise: the world points, their observed
/// pixels, and the track index each came from.
type Correspondences = (Vec<Vector3<f64>>, Vec<(f64, f64)>, Vec<usize>);

// ---------------------------------------------------------------------------
// Union-find (disjoint-set) over global feature ids, for track building.
// ---------------------------------------------------------------------------

/// A minimal union-find with path compression and union-by-rank, over a fixed
/// number of elements (the flattened `(image, keypoint)` feature ids).
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        // Iterative find with full path compression (no recursion, so a long
        // chain cannot overflow the stack).
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        // Union by rank.
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

/// Build multi-view [`Track`]s from the pairwise matches by **union-find**
/// (see the [module docs](self)).
///
/// Every feature is given a global id; each [`Match`] of each [`PairMatches`]
/// unions its two endpoints' ids; and every resulting connected component
/// becomes a candidate track. A component that observes **the same image
/// twice** (two different keypoints of one image) is **inconsistent** — it
/// cannot be a single rigid 3-D point — and is **dropped**. Singleton
/// components (a feature that was never matched) are not tracks and are
/// omitted.
///
/// `keypoint_counts[i]` must be the number of keypoints in image `i` (i.e.
/// `images[i].keypoints.len()`); it bounds the per-image id ranges. Out-of-
/// range indices in a [`Match`] (defensive — a well-formed pipeline never
/// produces them) are skipped rather than panicking.
///
/// The returned tracks have `point3d == None` and `point_idx == None`
/// (triangulation happens later, during [`incremental_sfm`]); their
/// observations are sorted by `image_idx` for determinism.
#[must_use]
pub fn build_tracks(keypoint_counts: &[usize], pairwise: &[PairMatches]) -> Vec<Track> {
    // Prefix offsets so feature (image, kp) maps to a single global id.
    let num_images = keypoint_counts.len();
    let mut offset = vec![0usize; num_images + 1];
    for i in 0..num_images {
        offset[i + 1] = offset[i] + keypoint_counts[i];
    }
    let total = offset[num_images];
    if total == 0 {
        return Vec::new();
    }

    // Map a (image, keypoint) to its global id, guarding ranges defensively.
    let gid = |im: usize, kp: usize| -> Option<usize> {
        if im >= num_images || kp >= keypoint_counts[im] {
            return None;
        }
        Some(offset[im] + kp)
    };

    let mut uf = UnionFind::new(total);
    // Track which feature ids actually participate in at least one match, so we
    // can skip the never-matched singletons when forming tracks.
    let mut touched = vec![false; total];

    for pm in pairwise {
        for m in &pm.matches {
            let (Some(a), Some(b)) = (gid(pm.i, m.query_idx), gid(pm.j, m.train_idx)) else {
                continue;
            };
            uf.union(a, b);
            touched[a] = true;
            touched[b] = true;
        }
    }

    // Group touched feature ids by their union-find root.
    use std::collections::HashMap;
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for (id, &is_touched) in touched.iter().enumerate() {
        if is_touched {
            let root = uf.find(id);
            groups.entry(root).or_default().push(id);
        }
    }

    // Decode each group back to (image, keypoint) observations and keep only
    // the consistent (one-keypoint-per-image) components of size >= 2.
    //
    // Inverse map global id -> (image, keypoint): find the image whose offset
    // range contains the id (offsets are sorted ascending).
    let to_obs = |id: usize| -> (usize, usize) {
        // Binary search the image owning this id.
        let im = match offset.binary_search(&id) {
            // Exact hit on offset[im] means id is the first feature of image im.
            Ok(k) => {
                // offset can repeat for zero-keypoint images; advance to the
                // image that actually owns the slot (the last offset == id that
                // is followed by a larger one).
                let mut k = k;
                while k + 1 < offset.len() && offset[k + 1] == id {
                    k += 1;
                }
                k
            }
            // Err(k): id falls between offset[k-1] and offset[k], so image k-1.
            Err(k) => k - 1,
        };
        (im, id - offset[im])
    };

    let mut tracks: Vec<Track> = Vec::new();
    for (_root, members) in groups {
        if members.len() < 2 {
            continue; // a singleton is not a track
        }
        let mut observations: Vec<(usize, usize)> = members.iter().map(|&id| to_obs(id)).collect();
        observations.sort_unstable();

        // Reject an inconsistent track: the same image appearing twice.
        let mut consistent = true;
        for w in observations.windows(2) {
            if w[0].0 == w[1].0 {
                consistent = false;
                break;
            }
        }
        if !consistent {
            continue;
        }

        tracks.push(Track {
            observations,
            point3d: None,
            point_idx: None,
        });
    }

    // Deterministic order (HashMap iteration is unordered): sort tracks by
    // their first observation.
    tracks.sort_by(|a, b| a.observations.cmp(&b.observations));
    tracks
}

/// Build the `3×4` projection matrix `P = K [R | t]`.
fn projection_matrix(k: &Matrix3<f64>, r: &Matrix3<f64>, t: &Vector3<f64>) -> Matrix3x4<f64> {
    let mut rt = Matrix3x4::<f64>::zeros();
    rt.fixed_view_mut::<3, 3>(0, 0).copy_from(r);
    rt.fixed_view_mut::<3, 1>(0, 3).copy_from(t);
    k * rt
}

/// Pixel `(u, v)` of `images[image_idx].keypoints[kp]`, if in range.
fn pixel_of(images: &[ImageFeatures], image_idx: usize, kp: usize) -> Option<(f64, f64)> {
    let p = images.get(image_idx)?.keypoints.get(kp)?;
    Some((f64::from(p.x), f64::from(p.y)))
}

/// Pick the initial image pair: the pair with the most verified inliers among
/// those clearing [`MapperParams::min_initial_inliers`].
///
/// Returns the index into `pairwise` of the chosen pair, or [`None`] if no pair
/// clears the bar. Ties are broken by the lowest `(i, j)` for determinism.
fn pick_initial_pair(pairwise: &[PairMatches], min_inliers: usize) -> Option<usize> {
    let mut best: Option<(usize, usize, (usize, usize))> = None; // (inliers, idx, (i,j))
    for (idx, pm) in pairwise.iter().enumerate() {
        let n = pm.matches.len();
        if n < min_inliers {
            continue;
        }
        let take = match best {
            None => true,
            Some((bn, _, (bi, bj))) => n > bn || (n == bn && (pm.i, pm.j) < (bi, bj)),
        };
        if take {
            best = Some((n, idx, (pm.i, pm.j)));
        }
    }
    best.map(|(_, idx, _)| idx)
}

/// Internal mutable mapper state during incremental reconstruction.
struct Mapper<'a> {
    images: &'a [ImageFeatures],
    k: CameraIntrinsics,
    /// Per-image pose (registered = `Some`).
    cameras: Vec<Option<CameraPose>>,
    /// Triangulated 3-D points.
    points: Vec<Vector3<f64>>,
    /// Tracks (with their point assignments updated as we triangulate).
    tracks: Vec<Track>,
    /// For each image, the indices of tracks that observe it (built once).
    tracks_of_image: Vec<Vec<usize>>,
}

impl<'a> Mapper<'a> {
    fn new(images: &'a [ImageFeatures], k: CameraIntrinsics, tracks: Vec<Track>) -> Self {
        let num_images = images.len();
        let mut tracks_of_image = vec![Vec::new(); num_images];
        for (ti, tr) in tracks.iter().enumerate() {
            for &(im, _) in &tr.observations {
                if im < num_images {
                    tracks_of_image[im].push(ti);
                }
            }
        }
        Self {
            images,
            k,
            cameras: vec![None; num_images],
            points: Vec::new(),
            tracks,
            tracks_of_image,
        }
    }

    /// Is `image_idx` registered?
    fn is_registered(&self, image_idx: usize) -> bool {
        self.cameras
            .get(image_idx)
            .map(Option::is_some)
            .unwrap_or(false)
    }

    /// Triangulate one track from two specific registered views, returning the
    /// 3-D point (in the reference frame) if the rays were well-conditioned and
    /// the point lies in front of both cameras.
    fn triangulate_track_from(
        &self,
        track_idx: usize,
        im_a: usize,
        im_b: usize,
    ) -> Option<Vector3<f64>> {
        let tr = &self.tracks[track_idx];
        let kp_a = tr.keypoint_in(im_a)?;
        let kp_b = tr.keypoint_in(im_b)?;
        let pa = self.cameras[im_a].as_ref()?;
        let pb = self.cameras[im_b].as_ref()?;
        let x_a = pixel_of(self.images, im_a, kp_a)?;
        let x_b = pixel_of(self.images, im_b, kp_b)?;

        let km = self.k.matrix();
        let proj_a = projection_matrix(&km, &pa.rotation, &pa.translation);
        let proj_b = projection_matrix(&km, &pb.rotation, &pb.translation);
        let x = triangulate_point(&proj_a, &proj_b, x_a, x_b);

        // Cheirality: positive depth in both cameras.
        let depth_a = (pa.rotation * x + pa.translation).z;
        let depth_b = (pb.rotation * x + pb.translation).z;
        if depth_a > 0.0 && depth_b > 0.0 && x.iter().all(|c| c.is_finite()) {
            Some(x)
        } else {
            None
        }
    }

    /// Assign a freshly triangulated point to a track (appending to the point
    /// cloud and linking the track).
    fn set_track_point(&mut self, track_idx: usize, x: Vector3<f64>) {
        let pi = self.points.len();
        self.points.push(x);
        self.tracks[track_idx].point3d = Some(x);
        self.tracks[track_idx].point_idx = Some(pi);
    }

    /// Triangulate every not-yet-triangulated track that is observed by at
    /// least two **registered** views (using the first two such views). Returns
    /// the number of new points created.
    fn triangulate_all_ready(&mut self) -> usize {
        let mut created = 0usize;
        for ti in 0..self.tracks.len() {
            if self.tracks[ti].point_idx.is_some() {
                continue;
            }
            // Collect the registered views observing this track.
            let regs: Vec<usize> = self.tracks[ti]
                .observations
                .iter()
                .map(|(im, _)| *im)
                .filter(|&im| self.is_registered(im))
                .collect();
            if regs.len() < 2 {
                continue;
            }
            if let Some(x) = self.triangulate_track_from(ti, regs[0], regs[1]) {
                self.set_track_point(ti, x);
                created += 1;
            }
        }
        created
    }

    /// Gather the 2D–3D correspondences an unregistered `image_idx` shares with
    /// the current model: for each track that already has a 3-D point and
    /// observes this image, pair the image's keypoint pixel with the 3-D point.
    /// Returns `(points3d, points2d, track_indices)` aligned element-wise.
    fn correspondences_for(&self, image_idx: usize) -> Correspondences {
        let mut p3d = Vec::new();
        let mut p2d = Vec::new();
        let mut tix = Vec::new();
        for &ti in &self.tracks_of_image[image_idx] {
            let tr = &self.tracks[ti];
            let Some(x) = tr.point3d else { continue };
            let Some(kp) = tr.keypoint_in(image_idx) else {
                continue;
            };
            let Some(px) = pixel_of(self.images, image_idx, kp) else {
                continue;
            };
            p3d.push(x);
            p2d.push(px);
            tix.push(ti);
        }
        (p3d, p2d, tix)
    }

    /// Choose the next unregistered view to register: the one with the most
    /// 2D–3D correspondences into the current model (ties → lowest image
    /// index). Returns `(image_idx, correspondence_count)` or [`None`] if no
    /// unregistered view has at least `min_corr` correspondences.
    fn pick_next_view(&self, min_corr: usize) -> Option<(usize, usize)> {
        let mut best: Option<(usize, usize)> = None; // (count, image_idx)
        for im in 0..self.images.len() {
            if self.is_registered(im) {
                continue;
            }
            let (p3d, _, _) = self.correspondences_for(im);
            let count = p3d.len();
            if count < min_corr {
                continue;
            }
            let take = match best {
                None => true,
                Some((bc, bi)) => count > bc || (count == bc && im < bi),
            };
            if take {
                best = Some((count, im));
            }
        }
        best.map(|(count, im)| (im, count))
    }

    /// Run the (global-over-current-model) bundle adjustment on the registered
    /// cameras + triangulated points, writing the refined geometry back.
    fn run_bundle(&mut self, params: &BundleParams) {
        // Compact the registered cameras into a contiguous list, remembering
        // the mapping back to image indices. Camera 0 of the bundle is the
        // FIRST registered view (the gauge).
        let mut cam_image: Vec<usize> = Vec::new();
        let mut cam_index_of_image: Vec<Option<usize>> = vec![None; self.images.len()];
        let mut bundle_cameras: Vec<CameraPose> = Vec::new();
        for (im, c) in self.cameras.iter().enumerate() {
            if let Some(pose) = c {
                cam_index_of_image[im] = Some(bundle_cameras.len());
                cam_image.push(im);
                bundle_cameras.push(*pose);
            }
        }
        if bundle_cameras.len() < 2 || self.points.is_empty() {
            return; // nothing meaningful to refine
        }

        // Observations: for every triangulated track, every registered view
        // that observes it contributes one observation of that point.
        let mut observations: Vec<Observation> = Vec::new();
        for tr in &self.tracks {
            let Some(pi) = tr.point_idx else { continue };
            for &(im, kp) in &tr.observations {
                let Some(ci) = cam_index_of_image[im] else {
                    continue;
                };
                let Some(px) = pixel_of(self.images, im, kp) else {
                    continue;
                };
                observations.push(Observation {
                    camera_idx: ci,
                    point_idx: pi,
                    pixel: px,
                });
            }
        }
        if observations.is_empty() {
            return;
        }

        let problem = BundleProblem {
            cameras: bundle_cameras,
            points: self.points.clone(),
            intrinsics: self.k,
            observations,
        };
        let result = bundle_adjust(&problem, params);

        // Write refined cameras back to their image slots.
        for (ci, &im) in cam_image.iter().enumerate() {
            self.cameras[im] = Some(result.cameras[ci]);
        }
        // Write refined points back (and keep the tracks' cached point3d in
        // sync for any subsequent triangulation comparisons).
        self.points = result.points;
        for tr in &mut self.tracks {
            if let Some(pi) = tr.point_idx {
                tr.point3d = self.points.get(pi).copied();
            }
        }
    }

    /// Consume the mapper into a [`Reconstruction`].
    fn into_reconstruction(self) -> Reconstruction {
        Reconstruction {
            cameras: self.cameras,
            points: self.points,
            tracks: self.tracks,
        }
    }
}

/// Run the full **incremental Structure-from-Motion** pipeline: turn per-image
/// features and pairwise verified matches into a single reconstruction (the
/// capstone tying Stages 1–5 together — see the [module docs](self)).
///
/// # Arguments
///
/// - `images`: the per-image [`ImageFeatures`]; index `i` is image `i`.
/// - `pairwise`: the verified inlier [`PairMatches`] for image pairs. Each
///   carries the pair `(i, j)` and the inlier [`Match`]es (with `query_idx`
///   indexing image `i` and `train_idx` indexing image `j`).
/// - `k`: the shared camera [`CameraIntrinsics`] (this crate uses one calibrated
///   camera moving through the scene; per-camera intrinsics are a future
///   extension — see [`crate::bundle`]).
/// - `params`: the [`MapperParams`] tuning the initial-pair / next-view
///   thresholds and the PnP-RANSAC / bundle-adjustment behaviour.
///
/// # Pipeline
///
/// 1. **Tracks.** Chain the pairwise matches into multi-view tracks by
///    union-find ([`build_tracks`]); drop inconsistent (one-image-twice)
///    tracks.
/// 2. **Seed.** Pick the strongest eligible initial pair, recover its relative
///    pose with [`recover_pose`] (placing the first camera at the identity gauge
///    and the second at `(R, t)`), and triangulate the tracks the two views
///    share.
/// 3. **Grow.** Repeatedly register the unregistered view with the most 2D–3D
///    correspondences via [`solve_pnp_ransac`], then triangulate the new tracks
///    it shares with already-registered views, running [`bundle_adjust`] every
///    [`MapperParams::ba_every`] registrations.
/// 4. **Polish.** A final global [`bundle_adjust`].
///
/// # Returns
///
/// `Some(reconstruction)` with at least the two seed cameras registered, or
/// [`None`] when no adequate initial pair exists, the seed pose cannot be
/// recovered, or fewer than two images / no usable matches are supplied. Views
/// that cannot be registered (too few correspondences, or a degenerate —
/// collinear/coplanar — local configuration for the linear PnP) are simply left
/// `None` in [`Reconstruction::cameras`]; the mapper stops when no further view
/// can be added. Never panics.
///
/// The reconstruction is correct only **up to a global similarity** (rigid
/// gauge pinned by the first registered camera; the scale is free) — see the
/// [module docs](self).
#[must_use]
pub fn incremental_sfm(
    images: &[ImageFeatures],
    pairwise: &[PairMatches],
    k: &CameraIntrinsics,
    params: &MapperParams,
) -> Option<Reconstruction> {
    // Need at least two images and at least one pair to do anything.
    if images.len() < 2 || pairwise.is_empty() {
        return None;
    }

    // Stage: build the multi-view tracks.
    let keypoint_counts: Vec<usize> = images.iter().map(|im| im.keypoints.len()).collect();
    let tracks = build_tracks(&keypoint_counts, pairwise);
    if tracks.is_empty() {
        return None;
    }

    // Stage: choose the initial pair (strongest geometry / inlier support).
    let init_idx = pick_initial_pair(pairwise, params.min_initial_inliers)?;
    let init = &pairwise[init_idx];
    let (im0, im1) = (init.i, init.j);
    if im0 >= images.len() || im1 >= images.len() || im0 == im1 {
        return None;
    }

    // Seed pose via the two-view essential-matrix path. We re-derive E from a
    // fundamental matrix built off the inlier correspondences. Rather than
    // re-run Stage-2 RANSAC here, we recover the relative pose directly from the
    // essential matrix estimated from the verified inliers (the mapper assumes
    // the pairwise matches are already verified inliers).
    let recon_seed = seed_two_view(images, init, k)?;

    let mut mapper = Mapper::new(images, *k, tracks);
    // Place camera im0 at the identity gauge, im1 at the recovered (R, t).
    mapper.cameras[im0] = Some(CameraPose {
        rotation: Matrix3::identity(),
        translation: Vector3::zeros(),
    });
    mapper.cameras[im1] = Some(CameraPose {
        rotation: recon_seed.rotation,
        translation: recon_seed.translation,
    });

    // Triangulate the tracks shared by the seed pair.
    mapper.triangulate_all_ready();
    if mapper.points.is_empty() {
        // The seed produced no usable structure: report the bare two-camera
        // reconstruction (honest minimal result) rather than None — both
        // cameras ARE registered.
        return Some(mapper.into_reconstruction());
    }

    let min_corr = params.min_pnp_correspondences.max(MIN_PNP);
    let mut registrations_since_ba = 0usize;

    // Stage: incremental growth loop.
    loop {
        let Some((next_im, _count)) = mapper.pick_next_view(min_corr) else {
            break; // no registerable view remains
        };

        let (p3d, p2d, _tix) = mapper.correspondences_for(next_im);
        // Defensive: pick_next_view guaranteed >= min_corr >= 6.
        let Some(result) = solve_pnp_ransac(&p3d, &p2d, k, &params.pnp) else {
            // Could not register this view (degenerate local geometry, e.g.
            // collinear/coplanar points). Mark it so we do not retry forever:
            // we cannot register it now, and its correspondence count will not
            // grow until more points appear — but to avoid an infinite loop we
            // must make progress. Triangulate any newly-ready tracks (none new
            // here) and, since pick_next_view will keep returning this same
            // best view, break out.
            //
            // To be robust we attempt the *next* candidate by temporarily
            // skipping: re-run selection excluding failures.
            if !try_register_any_other(&mut mapper, next_im, min_corr, params) {
                break;
            }
            continue;
        };

        mapper.cameras[next_im] = Some(result.pose);
        registrations_since_ba += 1;

        // Triangulate the new tracks this view shares with registered views.
        mapper.triangulate_all_ready();

        // Periodic (global-over-current-model) bundle adjustment.
        if params.ba_every > 0 && registrations_since_ba >= params.ba_every {
            mapper.run_bundle(&params.bundle);
            registrations_since_ba = 0;
        }
    }

    // Stage: final global bundle adjustment.
    mapper.run_bundle(&params.bundle);

    Some(mapper.into_reconstruction())
}

/// Recover the seed relative pose for the initial pair from its verified inlier
/// matches, by the essential-matrix path.
///
/// We estimate the essential matrix directly from the calibrated correspondences
/// (the inliers are assumed already geometrically verified), then decompose +
/// cheirality-select via [`recover_pose`]. To do that we first need a
/// fundamental matrix; we obtain it from the normalized 8-point algorithm
/// applied to the inliers. Rather than depend on Stage-2's RANSAC wrapper
/// (which re-filters), we reuse the public [`crate::verify_two_view`] on the
/// inlier set with a permissive threshold so essentially all of them survive,
/// yielding `F`, and proceed.
fn seed_two_view(
    images: &[ImageFeatures],
    pair: &PairMatches,
    k: &CameraIntrinsics,
) -> Option<crate::twoview::TwoViewReconstruction> {
    let kp_i = &images.get(pair.i)?.keypoints;
    let kp_j = &images.get(pair.j)?.keypoints;

    // Estimate F from the inliers via the existing Stage-2 verifier. A generous
    // RANSAC threshold + enough iterations make the verifier keep the (already
    // clean) inliers and return their fundamental matrix.
    let ransac = crate::geometry::RansacParams {
        iterations: 1000,
        threshold: 3.0,
        seed: 0x5EED_5EED_5EED_5EED,
    };
    let tv = crate::geometry::verify_two_view(kp_i, kp_j, &pair.matches, &ransac)?;
    let e = essential_from_fundamental(&tv.fundamental, k, k);
    recover_pose(&e, k, k, kp_i, kp_j, &tv.inliers)
}

/// Fallback used when the *best* next view fails PnP: try to register any other
/// eligible unregistered view (excluding `failed`). Returns `true` if some view
/// was registered (so the main loop should continue), `false` if none could be
/// — in which case the main loop breaks.
///
/// This keeps a single degenerate view (e.g. collinear local points) from
/// stalling the whole mapper: we skip it and let the others in. A view skipped
/// here may still become registerable on a later pass once more points exist,
/// because selection is recomputed from scratch each outer iteration.
fn try_register_any_other(
    mapper: &mut Mapper,
    failed: usize,
    min_corr: usize,
    params: &MapperParams,
) -> bool {
    // Collect candidate views (unregistered, enough correspondences), best
    // first, excluding the one that just failed.
    let mut candidates: Vec<(usize, usize)> = Vec::new(); // (count, image)
    for im in 0..mapper.images.len() {
        if im == failed || mapper.is_registered(im) {
            continue;
        }
        let (p3d, _, _) = mapper.correspondences_for(im);
        if p3d.len() >= min_corr {
            candidates.push((p3d.len(), im));
        }
    }
    candidates.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));

    for (_, im) in candidates {
        let (p3d, p2d, _) = mapper.correspondences_for(im);
        if let Some(result) = solve_pnp_ransac(&p3d, &p2d, &mapper.k, &params.pnp) {
            mapper.cameras[im] = Some(result.pose);
            mapper.triangulate_all_ready();
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::twoview::CameraIntrinsics;

    fn test_intrinsics() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// Rotation from yaw (about y) then pitch (about x), radians.
    fn rot_yaw_pitch(yaw: f64, pitch: f64) -> Matrix3<f64> {
        let ry = Matrix3::new(
            yaw.cos(),
            0.0,
            yaw.sin(),
            0.0,
            1.0,
            0.0,
            -yaw.sin(),
            0.0,
            yaw.cos(),
        );
        let rx = Matrix3::new(
            1.0,
            0.0,
            0.0,
            0.0,
            pitch.cos(),
            -pitch.sin(),
            0.0,
            pitch.sin(),
            pitch.cos(),
        );
        ry * rx
    }

    /// Twelve non-coplanar 3-D world points spanning a range of depths.
    fn scene_points() -> Vec<Vector3<f64>> {
        let raw = [
            (-0.4, -0.3, 5.0),
            (0.3, -0.2, 6.0),
            (0.1, 0.4, 4.5),
            (-0.2, 0.25, 7.0),
            (0.45, 0.35, 5.5),
            (-0.35, -0.15, 6.5),
            (0.2, -0.35, 4.8),
            (-0.1, 0.1, 8.0),
            (0.05, -0.05, 5.2),
            (0.4, 0.05, 6.8),
            (-0.45, 0.3, 4.2),
            (0.15, 0.2, 7.5),
        ];
        raw.iter().map(|&(x, y, z)| Vector3::new(x, y, z)).collect()
    }

    /// Project a world point through `K, R, t` to a pixel.
    fn project(
        k: &Matrix3<f64>,
        r: &Matrix3<f64>,
        t: &Vector3<f64>,
        x: &Vector3<f64>,
    ) -> (f64, f64) {
        let cam = r * x + t;
        let px = k * cam;
        (px.x / px.z, px.y / px.z)
    }

    /// A camera rig of `n` poses. Camera 0 is the reference (identity, origin);
    /// the rest have small yaw/pitch + a sideways baseline.
    fn camera_rig(n: usize) -> Vec<CameraPose> {
        let mut cams = vec![CameraPose {
            rotation: Matrix3::identity(),
            translation: Vector3::zeros(),
        }];
        // A few generic poses; cycle if n is large.
        let presets = [
            (0.15, -0.08, Vector3::new(-0.8, 0.05, 0.10)),
            (-0.10, 0.12, Vector3::new(0.7, -0.06, 0.20)),
            (0.22, 0.05, Vector3::new(-0.3, 0.20, -0.15)),
            (-0.18, -0.10, Vector3::new(0.5, 0.10, 0.25)),
        ];
        for idx in 1..n {
            let (yaw, pitch, t) = presets[(idx - 1) % presets.len()];
            cams.push(CameraPose {
                rotation: rot_yaw_pitch(yaw, pitch),
                translation: t,
            });
        }
        cams.truncate(n);
        cams
    }

    /// Build `images` (one `ImageFeatures` per camera) and the *exact* pairwise
    /// matches: every camera observes every point, keypoint index == point
    /// index, so the matches between any two views are the identity 1:1 set.
    fn synth_scene(
        k: &CameraIntrinsics,
        cameras: &[CameraPose],
        points: &[Vector3<f64>],
    ) -> (Vec<ImageFeatures>, Vec<PairMatches>) {
        let km = k.matrix();
        // Per-image keypoints (one per point, in point order).
        let images: Vec<ImageFeatures> = cameras
            .iter()
            .map(|cam| {
                let kps = points
                    .iter()
                    .map(|x| {
                        let (u, v) = project(&km, &cam.rotation, &cam.translation, x);
                        Keypoint::new(u as f32, v as f32, 1.0)
                    })
                    .collect();
                ImageFeatures::new(kps)
            })
            .collect();

        // All-pairs exact matches.
        let mut pairwise = Vec::new();
        for i in 0..cameras.len() {
            for j in (i + 1)..cameras.len() {
                let matches = (0..points.len())
                    .map(|p| Match {
                        query_idx: p,
                        train_idx: p,
                        distance: 0,
                    })
                    .collect();
                pairwise.push(PairMatches { i, j, matches });
            }
        }
        (images, pairwise)
    }

    /// RMSE (px) of a reconstruction over all observations of triangulated
    /// tracks, using the registered cameras.
    fn reconstruction_rmse(
        recon: &Reconstruction,
        images: &[ImageFeatures],
        k: &CameraIntrinsics,
    ) -> f64 {
        use crate::pnp::project_point;
        let mut sum = 0.0;
        let mut n = 0usize;
        for tr in &recon.tracks {
            let Some(pi) = tr.point_idx else { continue };
            let x = recon.points[pi];
            for &(im, kp) in &tr.observations {
                let Some(Some(pose)) = recon.cameras.get(im) else {
                    continue;
                };
                let Some(px) = pixel_of(images, im, kp) else {
                    continue;
                };
                if let Some((pu, pv)) = project_point(k, &pose.rotation, &pose.translation, &x) {
                    sum += (pu - px.0).powi(2) + (pv - px.1).powi(2);
                    n += 1;
                }
            }
        }
        if n == 0 {
            0.0
        } else {
            (sum / n as f64).sqrt()
        }
    }

    // ===== Test 1: synthetic 3-view scene -> all cameras + points recovered. =
    #[test]
    fn three_view_full_reconstruction() {
        let k = test_intrinsics();
        let cameras = camera_rig(3);
        let points = scene_points();
        let (images, pairwise) = synth_scene(&k, &cameras, &points);

        let params = MapperParams {
            // 12 points → 12 inliers per pair; lower the floor so the pair is
            // eligible.
            min_initial_inliers: 8,
            ..MapperParams::default()
        };
        let recon = incremental_sfm(&images, &pairwise, &k, &params)
            .expect("3-view synthetic scene must reconstruct");

        // All three cameras registered.
        assert_eq!(
            recon.num_registered(),
            3,
            "all 3 cameras should be registered, got {}",
            recon.num_registered()
        );
        for (i, c) in recon.cameras.iter().enumerate() {
            assert!(c.is_some(), "camera {i} not registered");
        }

        // All 12 points triangulated.
        let triangulated = recon
            .tracks
            .iter()
            .filter(|t| t.point_idx.is_some())
            .count();
        assert_eq!(
            triangulated,
            points.len(),
            "all 12 tracks should triangulate"
        );
        assert_eq!(recon.points.len(), points.len());

        // Reprojection RMSE ~ 0 after BA.
        let rmse = reconstruction_rmse(&recon, &images, &k);
        assert!(
            rmse < 1e-3,
            "post-BA reprojection RMSE {rmse} px should be ~0"
        );

        // Structure recovered UP TO A GLOBAL SIMILARITY: rotations match truth
        // directly; translations + points match under one common scale, with
        // the first registered camera as the gauge.
        //
        // The first registered camera is the initial pair's `i`. With our
        // all-pairs scene the strongest pair is (0,1) (ties broken low), so the
        // gauge is camera 0 = identity, matching `cameras[0]`.
        let cam0 = recon.cameras[0].expect("cam0 registered");
        assert!(
            (cam0.rotation - Matrix3::identity()).norm() < 1e-9,
            "gauge camera 0 should be identity rotation"
        );
        assert!(cam0.translation.norm() < 1e-9, "gauge camera 0 at origin");

        // Recover the global scale from camera 1's translation vs truth.
        let cam1 = recon.cameras[1].expect("cam1 registered");
        let s = cam1.translation.norm() / cameras[1].translation.norm();
        assert!(s.is_finite() && s > 0.0, "scale must be positive, got {s}");

        // Rotations exact; translations consistent under that one scale.
        for (j, c) in recon.cameras.iter().enumerate() {
            let pose = c.expect("registered");
            let r_err = (pose.rotation - cameras[j].rotation).norm();
            assert!(r_err < 1e-5, "camera {j} rotation off by {r_err}");
            let t_err = (pose.translation - cameras[j].translation * s).norm();
            assert!(
                t_err < 1e-4,
                "camera {j} translation inconsistent with scale {s}: err {t_err}"
            );
        }

        // Every triangulated point matches its truth under the SAME scale.
        for tr in &recon.tracks {
            let Some(pi) = tr.point_idx else { continue };
            // The track's single observation in image 0 identifies the point id
            // (keypoint index == point index in the synthetic scene).
            let kp0 = tr.keypoint_in(0).expect("track observes image 0");
            let truth = points[kp0];
            let err = (recon.points[pi] - truth * s).norm();
            assert!(
                err < 1e-4,
                "point {kp0} off truth under scale {s}: err {err}"
            );
        }
    }

    // ===== Test 2: track building — correct merge + inconsistent rejection. ==
    #[test]
    fn track_building_chains_and_rejects() {
        // Three images, keypoint counts large enough for the indices used.
        let counts = vec![5usize, 5, 5];

        // A consistent chain: (img0,kp0)-(img1,kp1) and (img1,kp1)-(img2,kp2)
        // must merge into ONE track over {(0,0),(1,1),(2,2)}.
        let pair01 = PairMatches {
            i: 0,
            j: 1,
            matches: vec![Match {
                query_idx: 0,
                train_idx: 1,
                distance: 0,
            }],
        };
        let pair12 = PairMatches {
            i: 1,
            j: 2,
            matches: vec![Match {
                query_idx: 1,
                train_idx: 2,
                distance: 0,
            }],
        };

        let tracks = build_tracks(&counts, &[pair01.clone(), pair12.clone()]);
        // Exactly one track, with the three chained observations.
        assert_eq!(tracks.len(), 1, "chained matches should form ONE track");
        assert_eq!(
            tracks[0].observations,
            vec![(0, 0), (1, 1), (2, 2)],
            "track should chain across all three views"
        );

        // Now add an INCONSISTENT link: (img2,kp2)-(img0,kp3). This drags a
        // SECOND keypoint of image 0 (kp3) into the same component as kp0,
        // so the component observes image 0 twice -> must be REJECTED.
        let pair20_bad = PairMatches {
            i: 2,
            j: 0,
            matches: vec![Match {
                query_idx: 2,
                train_idx: 3,
                distance: 0,
            }],
        };
        let tracks_bad = build_tracks(&counts, &[pair01, pair12, pair20_bad]);
        assert!(
            tracks_bad.is_empty(),
            "a component observing one image twice must be dropped, got {tracks_bad:?}"
        );

        // A separate independent two-view match forms its own clean track and
        // is unaffected by the rejected one.
        let counts2 = vec![3usize, 3];
        let indep = PairMatches {
            i: 0,
            j: 1,
            matches: vec![
                Match {
                    query_idx: 0,
                    train_idx: 0,
                    distance: 0,
                },
                Match {
                    query_idx: 2,
                    train_idx: 1,
                    distance: 0,
                },
            ],
        };
        let tracks2 = build_tracks(&counts2, &[indep]);
        assert_eq!(tracks2.len(), 2, "two independent matches -> two tracks");
        // Each is a clean two-view track.
        for t in &tracks2 {
            assert_eq!(t.observations.len(), 2);
            assert_ne!(t.observations[0].0, t.observations[1].0);
        }
    }

    // ===== Test 3: best-initial-pair selection picks strongest support. ======
    #[test]
    fn initial_pair_picks_strongest() {
        // Three pairs with different inlier counts; the selector must pick the
        // richest one (pair (1,2) with 40), and respect the min floor.
        let mk = |i: usize, j: usize, n: usize| PairMatches {
            i,
            j,
            matches: (0..n)
                .map(|p| Match {
                    query_idx: p,
                    train_idx: p,
                    distance: 0,
                })
                .collect(),
        };
        let pairwise = vec![mk(0, 1, 15), mk(1, 2, 40), mk(0, 2, 25)];

        // With a floor of 10, the strongest (40) wins: that is index 1.
        let idx = pick_initial_pair(&pairwise, 10).expect("a pair clears the floor");
        assert_eq!(idx, 1, "should pick the 40-inlier pair (1,2)");
        assert_eq!((pairwise[idx].i, pairwise[idx].j), (1, 2));

        // With a floor of 30, only the 40-inlier pair is eligible.
        let idx2 = pick_initial_pair(&pairwise, 30).expect("the 40-pair clears 30");
        assert_eq!(idx2, 1);

        // With a floor above all of them, no pair is eligible.
        assert!(
            pick_initial_pair(&pairwise, 100).is_none(),
            "no pair clears a floor of 100"
        );
    }

    // ===== Test 4: degenerate inputs -> None or minimal, never panic. ========
    #[test]
    fn degenerate_inputs_are_graceful() {
        let k = test_intrinsics();
        let params = MapperParams {
            min_initial_inliers: 4,
            ..MapperParams::default()
        };

        // (a) A single image, no pairs -> None.
        let one = vec![ImageFeatures::new(vec![Keypoint::new(1.0, 2.0, 1.0)])];
        assert!(
            incremental_sfm(&one, &[], &k, &params).is_none(),
            "single image / no pairs -> None"
        );

        // (b) Two images but NO matches -> None.
        let two = vec![
            ImageFeatures::new(vec![Keypoint::new(1.0, 2.0, 1.0)]),
            ImageFeatures::new(vec![Keypoint::new(3.0, 4.0, 1.0)]),
        ];
        assert!(
            incremental_sfm(&two, &[], &k, &params).is_none(),
            "two images, no pairwise -> None"
        );

        // (c) Pairs present but all below the initial-inlier floor -> None.
        let weak = vec![PairMatches {
            i: 0,
            j: 1,
            matches: vec![Match {
                query_idx: 0,
                train_idx: 0,
                distance: 0,
            }],
        }];
        assert!(
            incremental_sfm(&two, &weak, &k, &params).is_none(),
            "all pairs below the inlier floor -> None"
        );

        // (d) All-COLLINEAR scene points: the two-view seed may still place two
        // cameras, but later PnP registration is degenerate. Must not panic;
        // returns Some (>=2 cams) or None, never a crash.
        let cameras = camera_rig(3);
        // Points on a single 3-D line (collinear), all in front of the cameras.
        let line: Vec<Vector3<f64>> = (0..12)
            .map(|i| {
                let s = i as f64;
                Vector3::new(-0.5 + 0.1 * s, 0.2, 5.0 + 0.3 * s)
            })
            .collect();
        let (images, pairwise) = synth_scene(&k, &cameras, &line);
        let out = incremental_sfm(&images, &pairwise, &k, &params);
        if let Some(recon) = out {
            // Whatever it managed, it must be well-formed and finite.
            assert!(
                recon.points.iter().all(|p| p.iter().all(|c| c.is_finite())),
                "collinear-scene points must stay finite"
            );
            for c in recon.cameras.iter().flatten() {
                assert!(
                    (c.rotation.determinant() - 1.0).abs() < 1e-6,
                    "registered rotation must be proper"
                );
                assert!(c.translation.iter().all(|x| x.is_finite()));
            }
        }
    }

    // ===== Test 5: a 4th view registers incrementally via PnP. ===============
    #[test]
    fn fourth_view_registers_incrementally() {
        let k = test_intrinsics();
        let cameras = camera_rig(4);
        let points = scene_points();
        let (images, pairwise) = synth_scene(&k, &cameras, &points);

        let params = MapperParams {
            min_initial_inliers: 8,
            ..MapperParams::default()
        };
        let recon = incremental_sfm(&images, &pairwise, &k, &params)
            .expect("4-view synthetic scene must reconstruct");

        // All four cameras registered (the 4th necessarily via incremental PnP,
        // since only two are seeded).
        assert_eq!(
            recon.num_registered(),
            4,
            "all 4 cameras should register, got {}",
            recon.num_registered()
        );

        // Reprojection RMSE ~ 0 after BA.
        let rmse = reconstruction_rmse(&recon, &images, &k);
        assert!(
            rmse < 1e-3,
            "post-BA RMSE {rmse} px should be ~0 for 4 views"
        );

        // Up-to-similarity check on the 4th camera specifically: recover scale
        // from camera 1, verify camera 3's rotation (exact) and translation
        // (under that scale).
        let cam1 = recon.cameras[1].expect("cam1");
        let s = cam1.translation.norm() / cameras[1].translation.norm();
        let cam3 = recon.cameras[3].expect("cam3 registered via PnP");
        let r_err = (cam3.rotation - cameras[3].rotation).norm();
        assert!(r_err < 1e-5, "4th camera rotation off by {r_err}");
        let t_err = (cam3.translation - cameras[3].translation * s).norm();
        assert!(
            t_err < 1e-4,
            "4th camera translation inconsistent with scale {s}: err {t_err}"
        );
    }

    // ===== Bonus: reconstruction is expressed in the first camera's frame. ===
    #[test]
    fn gauge_is_first_registered_camera() {
        let k = test_intrinsics();
        let cameras = camera_rig(3);
        let points = scene_points();
        let (images, pairwise) = synth_scene(&k, &cameras, &points);
        let params = MapperParams {
            min_initial_inliers: 8,
            ..MapperParams::default()
        };
        let recon = incremental_sfm(&images, &pairwise, &k, &params).expect("reconstructs");

        // The gauge camera (image 0 here) is exactly identity/origin.
        let cam0 = recon.cameras[0].expect("cam0");
        assert!((cam0.rotation - Matrix3::identity()).norm() < 1e-9);
        assert!(cam0.translation.norm() < 1e-9);
    }
}

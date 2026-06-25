//! The right-side **Photogrammetry Workbench** panel — a native front-end over
//! the in-house [`valenx_photogrammetry`] crate (a COLMAP-style
//! Structure-from-Motion / scan pipeline).
//!
//! Real photogrammetry starts from a folder of overlapping photographs, which
//! a desktop app loads through a file dialog. A file dialog cannot run in the
//! headless test / CI environment, so this workbench instead drives the **real**
//! SfM solver on a fully-synthetic, fully-transparent **demo scene**: it
//! generates `N` camera views of a *known* 3-D point cloud (a cube grid),
//! projects every point into every view through a shared pinhole camera, adds
//! optional Gaussian pixel **noise**, and feeds those synthetic 2-D
//! observations through the genuine
//! [`valenx_photogrammetry::incremental_sfm`] mapper (two-view seed → PnP
//! registration → **bundle adjustment**). The recovered sparse points + camera
//! poses + reprojection error are exactly what the solver returns — the
//! workbench invents no number.
//!
//! The user sets three parameters — the number of cameras, the number of grid
//! points (rounded to a cube), and the per-pixel observation noise — clicks
//! **Run**, and valenx-photogrammetry then:
//!
//! * builds the synthetic camera rig + cube point cloud (the **ground truth**),
//!   projects + perturbs the observations, and runs the incremental mapper;
//! * reports the solver's own **quality metrics**
//!   ([`valenx_photogrammetry::reconstruction_quality`]): mean / median / max
//!   reprojection error (px), the registered-camera count, and the
//!   triangulated-point count;
//! * additionally measures the **recovery error** — the RMS distance between
//!   each recovered 3-D point and its known truth — *after* a least-squares
//!   **similarity alignment** (scale + rotation + translation) of the
//!   reconstruction onto the truth. A pure two-view-seeded SfM reconstruction
//!   is only determined **up to a global similarity** (the rigid gauge is
//!   pinned by the first camera, the absolute scale is free — see the
//!   [`valenx_photogrammetry::mapper`] / [`valenx_photogrammetry::bundle`]
//!   module docs), so comparing raw coordinates to truth would be meaningless;
//!   we factor that gauge out with an Umeyama-style closed-form similarity fit
//!   and report the residual.
//!
//! Two painter views are drawn (2-D orthographic projections of the recovered
//! 3-D cloud, *after* the alignment so they sit in the truth frame): an **XY**
//! (top-down) and an **XZ** (side) scatter of the recovered points plus the
//! recovered camera centres, with the known truth points underlaid for visual
//! comparison. A readout grid below gives the reprojection error, the #cameras
//! / #points, and the recovery RMS.
//!
//! Mirrors the other workbenches (`uq_workbench`, `rom_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_photogrammetry_workbench`], toggled from the View
//! menu and openable by the agent bridge under the workbench id
//! `"photogrammetry"` (aliases `"sfm"` / `"scan"`; see
//! [`crate::project_tabs::TabKind`]). Every numeric control is `.labelled_by`
//! an accessible caption so the panel is AI-drivable by name.
//!
//! Honesty: valenx-photogrammetry is **research / educational-grade** textbook
//! SfM (FAST/ORB features, the normalized 8-point + RANSAC fundamental matrix,
//! DLT triangulation / PnP, and dense Levenberg–Marquardt bundle adjustment) —
//! NOT a calibrated, survey-grade metrology product, and the reconstruction is
//! fixed only up to a global similarity. This workbench drives that solver on a
//! **synthetic** scene with **known** intrinsics (no auto-calibration); a real
//! photo set would additionally need feature detection + matching from pixels
//! and EXIF / calibrated intrinsics. On a **noise-free** scene the recovery is
//! essentially exact (reprojection error ≈ 0 and a near-perfect similarity fit
//! to truth), which the tests pin. Degenerate inputs (a single camera, or zero
//! points) surface an in-panel error, **not** a panic.

use eframe::egui;
use nalgebra::{Matrix3, Vector3};
use valenx_photogrammetry::{
    incremental_sfm, reconstruction_quality, CameraIntrinsics, ImageFeatures, Keypoint,
    MapperParams, Match, PairMatches, QualityParams,
};

use crate::agent_commands::AgentValue;
use crate::ValenxApp;

/// Fixed seed for the deterministic synthetic-noise PRNG so a run is
/// reproducible across machines (valenx-photogrammetry itself takes no `rand`
/// dependency; we add the observation noise here with a seeded SplitMix64-style
/// generator).
const PHOTOGRAMMETRY_SEED: u64 = 0x5CA1_AB1E_F00D_2025;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable inputs shown in the workbench: the synthetic scene size (number of
/// cameras + number of cube-grid points) and the observation noise.
#[derive(Clone, Copy, Debug)]
pub struct PhotogrammetryParams {
    /// Number of synthetic camera views to generate. Structure-from-motion
    /// needs at least **2** views (a single view cannot triangulate); the
    /// solver seeds two and registers the rest by PnP.
    pub num_cameras: usize,
    /// Requested number of 3-D scene points. The synthetic cloud is a cube
    /// grid, so the realised count is rounded **up** to the next perfect cube
    /// (`⌈∛n⌉³`); a clearly-non-coplanar cube is what the linear PnP needs.
    pub num_points: usize,
    /// Standard deviation of the zero-mean Gaussian noise (in **pixels**) added
    /// to every synthetic 2-D observation before the solver sees it. `0`
    /// reproduces the exact noise-free scene (analytic pin).
    pub noise_px: f64,
}

impl Default for PhotogrammetryParams {
    fn default() -> Self {
        Self {
            // A modest, well-conditioned default: 4 views of a 27-point cube,
            // a little sub-pixel noise so the default run exercises bundle
            // adjustment rather than hitting the exact-input fast path.
            num_cameras: 4,
            num_points: 27,
            noise_px: 0.3,
        }
    }
}

impl PhotogrammetryParams {
    /// Side length of the cube grid: `⌈∛(num_points)⌉`, at least 2 so the cloud
    /// is genuinely 3-D (not a single point / line / plane).
    fn grid_side(&self) -> usize {
        let n = self.num_points.max(1) as f64;
        (n.cbrt().ceil() as usize).max(2)
    }

    /// The realised number of cube-grid points (`grid_side³`).
    fn realised_points(&self) -> usize {
        let s = self.grid_side();
        s * s * s
    }
}

// ---------------------------------------------------------------------------
// Deterministic seeded noise (no `rand` dependency)
// ---------------------------------------------------------------------------

/// A tiny SplitMix64 generator + Box–Muller normal draw, so the synthetic
/// pixel noise is reproducible and dependency-free (mirrors the in-crate PRNG
/// style used elsewhere in valenx). Only used to *perturb* the synthetic
/// observations; the SfM solver itself is fully deterministic.
struct Splitmix {
    state: u64,
}

impl Splitmix {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next `u64`.
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform `f64` in the half-open `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        // Top 53 bits → a double in [0,1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// One standard-normal `N(0, 1)` draw by the Box–Muller transform.
    fn next_normal(&mut self) -> f64 {
        // Guard u1 away from 0 so ln() is finite.
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// A recovered camera centre + recovered 3-D point cloud, expressed in the
/// **truth frame** (after the similarity alignment), for the painter.
#[derive(Default, Clone)]
pub struct PhotogrammetryResult {
    /// Recovered 3-D points, aligned into the ground-truth frame by the
    /// similarity fit. Used for the scatter painters and the recovery RMS.
    pub recovered_points: Vec<Vector3<f64>>,
    /// The known ground-truth points (the synthetic cube grid). Underlaid in
    /// the painters for visual comparison.
    pub truth_points: Vec<Vector3<f64>>,
    /// Recovered camera centres `C = −Rᵀ t`, aligned into the truth frame.
    pub camera_centres: Vec<Vector3<f64>>,
    /// Number of registered (posed) cameras the mapper recovered.
    pub num_registered_cameras: usize,
    /// Number of triangulated 3-D points in the recovered cloud.
    pub num_points: usize,
    /// Mean reprojection error over all valid observations, in pixels
    /// (the solver's own [`valenx_photogrammetry::QualityMetrics`]).
    pub mean_reprojection_error: f64,
    /// Median reprojection error (px).
    pub median_reprojection_error: f64,
    /// Maximum reprojection error (px).
    pub max_reprojection_error: f64,
    /// RMS distance (in the truth's units) between each recovered point and its
    /// known truth point, **after** the closed-form similarity alignment that
    /// factors out the reconstruction's gauge freedom (scale + rigid).
    pub recovery_rms: f64,
    /// The scale factor the similarity alignment applied to the reconstruction
    /// to best match truth (the residual gauge freedom; ≈ 1 only by luck — the
    /// absolute scale of a two-view-seeded SfM is genuinely unrecoverable).
    pub alignment_scale: f64,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the photogrammetry workbench.
#[derive(Default)]
pub struct PhotogrammetryWorkbenchState {
    /// User-editable parameters.
    pub params: PhotogrammetryParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<PhotogrammetryResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl PhotogrammetryWorkbenchState {
    /// Run the full synthetic SfM pipeline, fail-loud.
    ///
    /// Builds a synthetic camera rig + cube point cloud (the ground truth),
    /// projects every point into every view through a shared pinhole camera,
    /// adds the requested Gaussian pixel noise, runs the **real**
    /// [`incremental_sfm`] mapper, reads the solver's quality metrics, and
    /// measures the recovery error after a similarity alignment onto truth.
    ///
    /// Every failure path returns an `Err(String)` — never a panic, never an
    /// invented number. Degenerate inputs (`num_cameras < 2`, or `num_points
    /// == 0`) are rejected up front; a mapper that cannot reconstruct (returns
    /// [`None`]) or quality metrics that fail are surfaced verbatim.
    pub fn run(&self) -> Result<PhotogrammetryResult, String> {
        let p = &self.params;

        // --- Degenerate-input guards (fail loud, no panic) ------------------
        if p.num_cameras < 2 {
            return Err(format!(
                "structure-from-motion needs >= 2 camera views (a single view cannot \
                 triangulate); got {}",
                p.num_cameras
            ));
        }
        if p.num_points == 0 {
            return Err("the synthetic scene needs >= 1 point (got 0)".to_string());
        }

        // --- Synthetic ground truth: intrinsics, cameras, cube point cloud --
        let k = synthetic_intrinsics();
        let cameras = synthetic_camera_rig(p.num_cameras);
        let truth_points = cube_point_cloud(p.grid_side());

        // --- Project + perturb into per-image features + all-pairs matches --
        let (images, pairwise) = synthesize_observations(&k, &cameras, &truth_points, p.noise_px)?;

        // --- Run the REAL incremental SfM mapper ---------------------------
        // Lower the initial-inlier floor below the point count so the
        // all-pairs synthetic scene is eligible to seed.
        let min_initial = truth_points.len().saturating_sub(2).clamp(6, 30);
        let params = MapperParams {
            min_initial_inliers: min_initial,
            ..MapperParams::default()
        };
        let recon = incremental_sfm(&images, &pairwise, &k, &params).ok_or_else(|| {
            "the SfM mapper could not reconstruct this scene (no adequate initial pair, \
             or the seed pose could not be recovered). Try more cameras / points or less \
             noise."
                .to_string()
        })?;

        // --- Solver's own quality metrics (fail-loud) ----------------------
        let quality = reconstruction_quality(&recon, &images, &k, &QualityParams::default())
            .map_err(|e| format!("quality metrics failed: {e}"))?;

        // --- Recovery error vs truth, after a similarity alignment ---------
        // Collect (recovered, truth) point pairs. Each triangulated track's
        // observation in image 0 identifies which truth point it is (keypoint
        // index == point index in the synthetic scene).
        let mut rec_pts: Vec<Vector3<f64>> = Vec::new();
        let mut tru_pts: Vec<Vector3<f64>> = Vec::new();
        for tr in &recon.tracks {
            let Some(pi) = tr.point_idx else { continue };
            let Some(&rec) = recon.points.get(pi) else {
                continue;
            };
            // The point id is the keypoint index in any observed image (all
            // images share point order in the synthetic scene).
            let Some(&(_im, kp)) = tr.observations.first() else {
                continue;
            };
            let Some(&tru) = truth_points.get(kp) else {
                continue;
            };
            rec_pts.push(rec);
            tru_pts.push(tru);
        }

        // Fit a similarity (scale s, rotation R, translation t) mapping the
        // recovered cloud onto truth, factoring out the gauge freedom, then
        // measure the residual RMS. With < 3 correspondences the fit is
        // ill-posed; report scale 1 / RMS 0 (the painters still draw).
        let (s, rot, trans, recovery_rms) = if rec_pts.len() >= 3 {
            similarity_align(&rec_pts, &tru_pts)
                .ok_or_else(|| "similarity alignment was ill-conditioned".to_string())?
        } else {
            (1.0, Matrix3::identity(), Vector3::zeros(), 0.0)
        };

        // Map the recovered points + camera centres into the truth frame for
        // display (so the scatter overlays line up with the truth underlay).
        let recovered_points: Vec<Vector3<f64>> =
            rec_pts.iter().map(|x| s * rot * x + trans).collect();
        let camera_centres: Vec<Vector3<f64>> = recon
            .cameras
            .iter()
            .flatten()
            .map(|pose| {
                // Camera centre in the reconstruction frame: C = −Rᵀ t.
                let c = -pose.rotation.transpose() * pose.translation;
                s * rot * c + trans
            })
            .collect();

        Ok(PhotogrammetryResult {
            recovered_points,
            truth_points,
            camera_centres,
            num_registered_cameras: quality.num_registered_cameras,
            num_points: quality.num_points,
            mean_reprojection_error: quality.mean_reprojection_error,
            median_reprojection_error: quality.median_reprojection_error,
            max_reprojection_error: quality.max_reprojection_error,
            recovery_rms,
            alignment_scale: s,
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. The captions match exactly what
    /// the workbench form draws (and what each control is `labelled_by`). The
    /// `realised points` readout is a read-only echo, not a control, so it is
    /// intentionally absent.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["camera views", "scene points", "image noise (px)"]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the wrong
    /// type returns `Err(String)` (the bridge turns it into a `warn` feed note) —
    /// never a panic, and no field is written on error. The two count captions
    /// read [`AgentValue::as_i64`]; the noise caption reads [`AgentValue::as_f64`].
    /// Range validation (e.g. `< 2` cameras) stays in [`run`](Self::run).
    pub fn agent_set(&mut self, name: &str, value: &AgentValue) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            "camera views" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("camera views must be >= 0, got {n}"));
                }
                p.num_cameras = n as usize;
            }
            "scene points" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("scene points must be >= 0, got {n}"));
                }
                p.num_points = n as usize;
            }
            "image noise (px)" => p.noise_px = value.as_f64()?,
            other => return Err(format!("unknown photogrammetry control: {other:?}")),
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Synthetic-scene construction
// ---------------------------------------------------------------------------

/// A fixed pinhole camera for the synthetic scene (640×480-ish, focal ≈ width).
fn synthetic_intrinsics() -> CameraIntrinsics {
    CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
}

/// A `n`-pose camera rig looking at the scene. Camera 0 is the reference
/// (identity rotation, origin), the rest have small yaw/pitch and a sideways
/// baseline (mirrors the in-crate test rig so the geometry is well-conditioned).
fn synthetic_camera_rig(n: usize) -> Vec<CameraPoseLite> {
    let mut cams = vec![CameraPoseLite {
        rotation: Matrix3::identity(),
        translation: Vector3::zeros(),
    }];
    let presets = [
        (0.15, -0.08, Vector3::new(-0.8, 0.05, 0.10)),
        (-0.10, 0.12, Vector3::new(0.7, -0.06, 0.20)),
        (0.22, 0.05, Vector3::new(-0.3, 0.20, -0.15)),
        (-0.18, -0.10, Vector3::new(0.5, 0.10, 0.25)),
        (0.08, 0.18, Vector3::new(-0.6, -0.12, 0.05)),
        (-0.22, -0.04, Vector3::new(0.4, 0.18, -0.20)),
    ];
    for idx in 1..n {
        let (yaw, pitch, t) = presets[(idx - 1) % presets.len()];
        cams.push(CameraPoseLite {
            rotation: rot_yaw_pitch(yaw, pitch),
            translation: t,
        });
    }
    cams.truncate(n);
    cams
}

/// A local camera-pose holder for synthesis (the crate's `CameraPose` is the
/// same shape; this avoids depending on the exact constructor and keeps the
/// synthesis self-contained).
#[derive(Clone, Copy)]
struct CameraPoseLite {
    rotation: Matrix3<f64>,
    translation: Vector3<f64>,
}

/// Rotation from yaw (about y) then pitch (about x), in radians.
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

/// A `side × side × side` cube grid of 3-D points centred on the scene and set
/// a few units in front of the cameras (positive z, so cheirality holds).
fn cube_point_cloud(side: usize) -> Vec<Vector3<f64>> {
    let mut pts = Vec::with_capacity(side * side * side);
    // Span roughly [-0.5, 0.5] in x/y and [4.5, 6.0] in z (in front of the
    // cameras), so the whole cube is comfortably imaged by every view.
    let denom = (side.max(2) - 1) as f64;
    for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                let fx = ix as f64 / denom - 0.5;
                let fy = iy as f64 / denom - 0.5;
                let fz = iz as f64 / denom;
                pts.push(Vector3::new(fx, fy, 4.5 + 1.5 * fz));
            }
        }
    }
    pts
}

/// Project a world point through `K, R, t` to a pixel `(u, v)`.
fn project(
    k: &Matrix3<f64>,
    r: &Matrix3<f64>,
    t: &Vector3<f64>,
    x: &Vector3<f64>,
) -> Option<(f64, f64)> {
    let cam = r * x + t;
    if cam.z.abs() < 1e-9 {
        return None;
    }
    let px = k * cam;
    Some((px.x / px.z, px.y / px.z))
}

/// Build the per-image features and the all-pairs verified matches from the
/// synthetic scene: every camera observes every point (keypoint index == point
/// index), with optional Gaussian pixel noise added to each observation.
///
/// Fails loud if a ground-truth point fails to project in front of some camera
/// (it should not for the well-conditioned synthetic rig, but we never silently
/// drop geometry into a `NaN`).
#[allow(clippy::type_complexity)]
fn synthesize_observations(
    k: &CameraIntrinsics,
    cameras: &[CameraPoseLite],
    points: &[Vector3<f64>],
    noise_px: f64,
) -> Result<(Vec<ImageFeatures>, Vec<PairMatches>), String> {
    let km = k.matrix();
    let mut rng = Splitmix::new(PHOTOGRAMMETRY_SEED);
    let sigma = noise_px.max(0.0);

    // Per-image keypoints, one per point in point order.
    let mut images: Vec<ImageFeatures> = Vec::with_capacity(cameras.len());
    for (ci, cam) in cameras.iter().enumerate() {
        let mut kps = Vec::with_capacity(points.len());
        for (pi, x) in points.iter().enumerate() {
            let (u, v) = project(&km, &cam.rotation, &cam.translation, x).ok_or_else(|| {
                format!(
                    "synthetic point {pi} does not project in front of camera {ci} \
                     (degenerate scene)"
                )
            })?;
            let (nu, nv) = if sigma > 0.0 {
                (u + sigma * rng.next_normal(), v + sigma * rng.next_normal())
            } else {
                (u, v)
            };
            kps.push(Keypoint::new(nu as f32, nv as f32, 1.0));
        }
        images.push(ImageFeatures::new(kps));
    }

    // All-pairs exact (identity) matches: keypoint p of image i ↔ keypoint p of
    // image j, for every i < j. These are the "verified inliers" the mapper
    // expects.
    let mut pairwise = Vec::new();
    for i in 0..cameras.len() {
        for j in (i + 1)..cameras.len() {
            let matches: Vec<Match> = (0..points.len())
                .map(|pt| Match {
                    query_idx: pt,
                    train_idx: pt,
                    distance: 0,
                })
                .collect();
            pairwise.push(PairMatches { i, j, matches });
        }
    }

    Ok((images, pairwise))
}

// ---------------------------------------------------------------------------
// Similarity alignment (Umeyama-style closed form) — factors the gauge freedom
// ---------------------------------------------------------------------------

/// Fit the similarity transform `(s, R, t)` minimizing `Σ‖ s·R·src + t − dst ‖²`
/// (Umeyama 1991, with scale), and return `(s, R, t, rms)` where `rms` is the
/// residual root-mean-square distance after the fit.
///
/// This factors out the **gauge freedom** of a two-view-seeded SfM
/// reconstruction (free scale + rigid frame) so the recovery error against the
/// known truth cloud is meaningful. Returns [`None`] if the source cloud is
/// degenerate (zero variance → the scale is undefined) or the SVD produced a
/// non-finite result.
fn similarity_align(
    src: &[Vector3<f64>],
    dst: &[Vector3<f64>],
) -> Option<(f64, Matrix3<f64>, Vector3<f64>, f64)> {
    let n = src.len();
    if n < 3 || n != dst.len() {
        return None;
    }
    let inv_n = 1.0 / n as f64;

    // Centroids.
    let mu_s: Vector3<f64> = src.iter().sum::<Vector3<f64>>() * inv_n;
    let mu_d: Vector3<f64> = dst.iter().sum::<Vector3<f64>>() * inv_n;

    // Cross-covariance Σ = (1/n) Σ (dst−μ_d)(src−μ_s)ᵀ, and the source variance.
    let mut cov = Matrix3::<f64>::zeros();
    let mut var_s = 0.0;
    for (a, b) in src.iter().zip(dst.iter()) {
        let ds = a - mu_s;
        let dd = b - mu_d;
        cov += dd * ds.transpose();
        var_s += ds.dot(&ds);
    }
    cov *= inv_n;
    var_s *= inv_n;

    if var_s < 1e-18 {
        // The recovered cloud is a single point (no scale information).
        return None;
    }

    // SVD of the cross-covariance: Σ = U S Vᵀ → R = U D Vᵀ with the reflection
    // fix (D = diag(1,1,det(UVᵀ))) so R stays a proper rotation.
    let svd = cov.svd(true, true);
    let u = svd.u?;
    let v_t = svd.v_t?;
    let mut d = Matrix3::<f64>::identity();
    let det = (u * v_t).determinant();
    if det < 0.0 {
        d[(2, 2)] = -1.0;
    }
    let rot = u * d * v_t;

    // Scale: s = trace(S·D) / var_s (Umeyama eq. 41).
    let sv = svd.singular_values;
    let trace_sd = sv[0] * d[(0, 0)] + sv[1] * d[(1, 1)] + sv[2] * d[(2, 2)];
    let s = trace_sd / var_s;

    // Translation: t = μ_d − s·R·μ_s.
    let trans = mu_d - s * rot * mu_s;

    if !(s.is_finite() && rot.iter().all(|x| x.is_finite()) && trans.iter().all(|x| x.is_finite()))
    {
        return None;
    }

    // Residual RMS after the fit.
    let mut sum = 0.0;
    for (a, b) in src.iter().zip(dst.iter()) {
        let mapped = s * rot * a + trans;
        sum += (mapped - b).norm_squared();
    }
    let rms = (sum * inv_n).sqrt();

    Some((s, rot, trans, rms))
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the photogrammetry workbench. A no-op unless toggled on via
/// View → Photogrammetry.
///
/// Mirrors [`crate::uq_workbench::draw_uq_workbench`].
pub fn draw_photogrammetry_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_photogrammetry_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_photogrammetry_workbench",
        "Photogrammetry / SfM scan",
        photogrammetry_workbench_body,
    );
    if close {
        app.show_photogrammetry_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn photogrammetry_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Photogrammetry / structure-from-motion \u{2014} a synthetic demo scene (known cube \
             point cloud + camera rig) is projected, noised, and fed through the REAL \
             valenx-photogrammetry mapper (two-view seed \u{2192} PnP registration \u{2192} \
             bundle adjustment). Recovers a sparse point cloud + camera poses. \
             [research / educational \u{2014} textbook SfM; up to a global similarity]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.photogrammetry;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Synthetic scene").strong());
        egui::Grid::new("photogrammetry_scene_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("camera views");
                ui.add(
                    egui::DragValue::new(&mut p.num_cameras)
                        .speed(1)
                        .range(0..=64),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Number of synthetic camera views generated of the scene. SfM needs \
                     >= 2 (a single view cannot triangulate); the mapper seeds two and \
                     registers the rest by PnP.",
                );
                ui.end_row();

                let lbl = ui.label("scene points");
                ui.add(
                    egui::DragValue::new(&mut p.num_points)
                        .speed(1)
                        .range(0..=4096),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Requested number of 3-D scene points. The cloud is a cube grid, so the \
                     realised count is rounded up to the next perfect cube (a non-coplanar \
                     cube is what the linear PnP needs).",
                );
                ui.end_row();

                let lbl = ui.label("image noise (px)");
                ui.add(
                    egui::DragValue::new(&mut p.noise_px)
                        .speed(0.05)
                        .range(0.0..=20.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Standard deviation of the zero-mean Gaussian noise added to every \
                     synthetic 2-D observation, in pixels. 0 reproduces the exact noise-free \
                     scene (near-perfect recovery).",
                );
                ui.end_row();

                // Read-only echo of the realised cube size (helps the user see
                // the rounding-to-a-cube). Not a control — just a caption.
                ui.label("realised points");
                ui.label(format!(
                    "{} ({}\u{00B3} cube)",
                    p.realised_points(),
                    p.grid_side()
                ));
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Generate the synthetic views, run the SfM mapper + bundle adjustment, \
                     and report the recovered cloud, camera poses, reprojection error, and \
                     the recovery error vs the known truth (after a similarity alignment).",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.photogrammetry;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_photogrammetry_viz(s, ui);
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.photogrammetry;
    match s.run() {
        Ok(res) => {
            s.status = format!(
                "\u{2714} {} cams \u{00B7} {} pts \u{00B7} reproj mean {:.4} px (max {:.4}) \
                 \u{00B7} recovery RMS {:.4} (scale {:.3})",
                res.num_registered_cameras,
                res.num_points,
                res.mean_reprojection_error,
                res.max_reprojection_error,
                res.recovery_rms,
                res.alignment_scale,
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (orthographic XY + XZ scatter of the recovered cloud)
// ---------------------------------------------------------------------------

fn draw_photogrammetry_viz(s: &PhotogrammetryWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to reconstruct the synthetic scene and visualise the recovered \
                 point cloud + camera positions",
            )
            .weak(),
        );
        return;
    };

    ui.label(egui::RichText::new("Recovered sparse cloud (orthographic)").strong());
    ui.label(
        egui::RichText::new(
            "grey dots = known truth points \u{00B7} cyan dots = recovered points \u{00B7} \
             amber squares = recovered camera centres \u{00B7} aligned to truth by a similarity fit",
        )
        .weak()
        .small(),
    );

    // Two side-by-side orthographic projections: XY (top-down) and XZ (side).
    ui.horizontal_wrapped(|ui| {
        draw_ortho_scatter(res, ui, Axis::X, Axis::Y, "XY (top-down)");
        draw_ortho_scatter(res, ui, Axis::X, Axis::Z, "XZ (side)");
    });

    // Readouts grid below the painters.
    ui.add_space(6.0);
    egui::Grid::new("photogrammetry_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(
                ui,
                "registered cameras",
                format!("{}", res.num_registered_cameras),
            );
            row(ui, "triangulated points", format!("{}", res.num_points));
            row(
                ui,
                "reproj error mean / median / max (px)",
                format!(
                    "{:.4} / {:.4} / {:.4}",
                    res.mean_reprojection_error,
                    res.median_reprojection_error,
                    res.max_reprojection_error
                ),
            );
            row(
                ui,
                "recovery RMS vs truth (post-align)",
                format!("{:.5}", res.recovery_rms),
            );
            row(
                ui,
                "similarity-alignment scale",
                format!("{:.5}", res.alignment_scale),
            );
        });
}

/// Which world axis maps to a painter axis in an orthographic projection.
#[derive(Clone, Copy)]
enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    fn of(self, p: &Vector3<f64>) -> f64 {
        match self {
            Axis::X => p.x,
            Axis::Y => p.y,
            Axis::Z => p.z,
        }
    }
}

/// Draw one orthographic scatter (truth + recovered points + camera centres),
/// projecting world `(horiz, vert)` axes onto the painter, auto-scaled to fit.
fn draw_ortho_scatter(
    res: &PhotogrammetryResult,
    ui: &mut egui::Ui,
    horiz: Axis,
    vert: Axis,
    title: &str,
) {
    ui.vertical(|ui| {
        ui.label(egui::RichText::new(title).small().weak());

        let (rect, _) = ui.allocate_exact_size(egui::vec2(220.0, 200.0), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

        // Compute the data bounds over everything we'll draw, so the scatter
        // auto-fits. Guard against an empty / degenerate (zero-extent) range.
        let mut lo_h = f64::INFINITY;
        let mut hi_h = f64::NEG_INFINITY;
        let mut lo_v = f64::INFINITY;
        let mut hi_v = f64::NEG_INFINITY;
        let mut acc = |p: &Vector3<f64>| {
            let h = horiz.of(p);
            let v = vert.of(p);
            if h.is_finite() && v.is_finite() {
                lo_h = lo_h.min(h);
                hi_h = hi_h.max(h);
                lo_v = lo_v.min(v);
                hi_v = hi_v.max(v);
            }
        };
        for p in &res.truth_points {
            acc(p);
        }
        for p in &res.recovered_points {
            acc(p);
        }
        for p in &res.camera_centres {
            acc(p);
        }

        if !(lo_h.is_finite() && hi_h.is_finite() && lo_v.is_finite() && hi_v.is_finite()) {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "no finite geometry",
                egui::FontId::monospace(12.0),
                egui::Color32::from_gray(120),
            );
            return;
        }
        // Pad zero-extent ranges so a single column of points still draws.
        if (hi_h - lo_h).abs() < 1e-9 {
            lo_h -= 0.5;
            hi_h += 0.5;
        }
        if (hi_v - lo_v).abs() < 1e-9 {
            lo_v -= 0.5;
            hi_v += 0.5;
        }

        let margin = 12.0_f32;
        let inner = rect.shrink(margin);
        let span_h = (hi_h - lo_h) as f32;
        let span_v = (hi_v - lo_v) as f32;

        // Map a world point to a painter position (vert axis flipped so +y/+z
        // points up on screen).
        let to_screen = |p: &Vector3<f64>| -> egui::Pos2 {
            let fh = ((horiz.of(p) - lo_h) as f32 / span_h).clamp(0.0, 1.0);
            let fv = ((vert.of(p) - lo_v) as f32 / span_v).clamp(0.0, 1.0);
            egui::pos2(
                inner.left() + fh * inner.width(),
                inner.bottom() - fv * inner.height(),
            )
        };

        // Truth points first (grey, underlaid).
        for p in &res.truth_points {
            painter.circle_filled(to_screen(p), 1.6, egui::Color32::from_gray(110));
        }
        // Recovered points (cyan).
        for p in &res.recovered_points {
            painter.circle_filled(to_screen(p), 1.8, egui::Color32::from_rgb(70, 200, 210));
        }
        // Recovered camera centres (amber squares).
        for c in &res.camera_centres {
            let cs = to_screen(c);
            let r = 3.0;
            painter.rect_filled(
                egui::Rect::from_center_size(cs, egui::vec2(r * 2.0, r * 2.0)),
                0.0,
                egui::Color32::from_rgb(230, 180, 70),
            );
        }
    });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring uq_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = PhotogrammetryWorkbenchState::default();
        let res = s.run().expect("default photogrammetry run should succeed");
        // All requested cameras register on the well-conditioned synthetic rig.
        assert_eq!(
            res.num_registered_cameras, s.params.num_cameras,
            "all cameras should register"
        );
        // The realised cube point count is triangulated.
        assert_eq!(
            res.num_points,
            s.params.realised_points(),
            "all cube points should triangulate"
        );
        assert_eq!(res.recovered_points.len(), res.num_points);
        assert_eq!(res.camera_centres.len(), res.num_registered_cameras);
        // Reprojection error is small (sub-pixel-noise scene → BA cleans it).
        assert!(
            res.mean_reprojection_error.is_finite() && res.mean_reprojection_error < 2.0,
            "mean reproj error {} px should be small",
            res.mean_reprojection_error
        );
        assert!(res.recovery_rms.is_finite(), "recovery RMS must be finite");
        assert!(
            res.alignment_scale.is_finite() && res.alignment_scale > 0.0,
            "alignment scale must be positive, got {}",
            res.alignment_scale
        );
    }

    #[test]
    fn noise_free_recovery_is_near_exact_pin() {
        // PIN (analytic): a noise-free synthetic scene recovers the cloud to
        // ~zero reprojection error AND matches the ground-truth geometry up to
        // the similarity gauge within a tight tolerance.
        let mut s = PhotogrammetryWorkbenchState::default();
        s.params.noise_px = 0.0;
        s.params.num_cameras = 4;
        s.params.num_points = 27; // 3x3x3 cube

        let res = s.run().expect("noise-free run should succeed");

        // Every camera + every point recovered.
        assert_eq!(res.num_registered_cameras, 4, "all 4 cameras register");
        assert_eq!(res.num_points, 27, "all 27 cube points triangulate");

        // Reprojection error ~ 0 (clean data, bundle adjustment converges).
        assert!(
            res.mean_reprojection_error < 1e-4,
            "noise-free mean reproj error {} px should be ~0",
            res.mean_reprojection_error
        );
        assert!(
            res.max_reprojection_error < 1e-3,
            "noise-free max reproj error {} px should be ~0",
            res.max_reprojection_error
        );

        // Geometry matches truth up to the similarity gauge: the residual RMS
        // after the similarity alignment is ~0 relative to the scene extent
        // (the cube spans ~1.5 units, so an RMS < 1e-3 is a near-perfect fit).
        assert!(
            res.recovery_rms < 1e-3,
            "noise-free recovery RMS {} (after similarity align) should be ~0",
            res.recovery_rms
        );
    }

    #[test]
    fn similarity_alignment_recovers_a_known_transform() {
        // Sanity-pin the gauge-handling math itself: apply a known similarity
        // (scale 2, a rotation, a translation) to a cloud and confirm the fit
        // inverts it with ~zero residual.
        let src = vec![
            Vector3::new(-0.4, -0.3, 5.0),
            Vector3::new(0.3, -0.2, 6.0),
            Vector3::new(0.1, 0.4, 4.5),
            Vector3::new(-0.2, 0.25, 7.0),
            Vector3::new(0.45, 0.35, 5.5),
        ];
        let s_true = 2.0;
        let r_true = rot_yaw_pitch(0.3, -0.2);
        let t_true = Vector3::new(1.0, -2.0, 0.5);
        let dst: Vec<Vector3<f64>> = src.iter().map(|x| s_true * r_true * x + t_true).collect();

        let (s, rot, trans, rms) = similarity_align(&src, &dst).expect("alignment should succeed");
        assert!((s - s_true).abs() < 1e-9, "scale off: {s} vs {s_true}");
        assert!((rot - r_true).norm() < 1e-9, "rotation off");
        assert!((trans - t_true).norm() < 1e-9, "translation off");
        assert!(
            rms < 1e-9,
            "residual RMS {rms} should be ~0 for an exact fit"
        );
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn single_camera_returns_err() {
        let mut s = PhotogrammetryWorkbenchState::default();
        s.params.num_cameras = 1;
        assert!(
            s.run().is_err(),
            "1 camera must return Err (SfM needs >= 2 views), not panic"
        );
    }

    #[test]
    fn zero_points_returns_err() {
        let mut s = PhotogrammetryWorkbenchState::default();
        s.params.num_points = 0;
        assert!(s.run().is_err(), "0 points must return Err, not panic");
    }

    #[test]
    fn noisy_scene_still_reconstructs() {
        // A moderately noisy scene should still reconstruct (more cameras /
        // points give the bundle adjuster the redundancy it needs).
        let mut s = PhotogrammetryWorkbenchState::default();
        s.params.noise_px = 1.0;
        s.params.num_cameras = 4;
        s.params.num_points = 27;
        let res = s.run().expect("noisy scene should still reconstruct");
        assert!(res.num_registered_cameras >= 2, "at least the seed pair");
        assert!(res.num_points > 0, "some points triangulated");
        assert!(
            res.mean_reprojection_error.is_finite(),
            "reproj error must stay finite under noise"
        );
    }

    // ---- agent_set / agent_control_names (the SetControl bridge) ----

    #[test]
    fn agent_set_sets_params_and_rejects_unknown_and_typemismatch() {
        let mut s = PhotogrammetryWorkbenchState::default();

        // Representative count param, verified via state.
        s.agent_set("camera views", &AgentValue::Int(6))
            .expect("set camera views");
        assert_eq!(s.params.num_cameras, 6);
        // Float noise param.
        s.agent_set("image noise (px)", &AgentValue::Float(0.75))
            .expect("set image noise");
        assert!((s.params.noise_px - 0.75).abs() < 1e-12);

        // Unknown caption -> Err (not a panic).
        assert!(s.agent_set("nope", &AgentValue::Int(1)).is_err());
        // Type mismatch: a count caption fed a string -> Err.
        assert!(s
            .agent_set("camera views", &AgentValue::Str("six".into()))
            .is_err());
        // Negative count -> Err.
        assert!(s.agent_set("scene points", &AgentValue::Int(-3)).is_err());

        // Every advertised control name is settable with a value of its type.
        for name in PhotogrammetryWorkbenchState::agent_control_names() {
            let v = if *name == "image noise (px)" {
                AgentValue::Float(0.1)
            } else {
                AgentValue::Int(4)
            };
            assert!(
                s.agent_set(name, &v).is_ok(),
                "advertised control '{name}' must be settable"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_photogrammetry_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_photogrammetry_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_photogrammetry_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_photogrammetry_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_photogrammetry_workbench = true;
        let res = app.photogrammetry.run().expect("run should succeed");
        app.photogrammetry.result = Some(res);
        app.photogrammetry.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_photogrammetry_workbench = true;
        // Trigger an error state (1 camera is fail-loud in run()).
        app.photogrammetry.params.num_cameras = 1;
        let result = app.photogrammetry.run();
        app.photogrammetry.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.photogrammetry.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_photogrammetry_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // Three numeric DragValues (cameras, points, noise) — all MUST carry an
        // accessible name (be labelled_by a caption) so the panel is AI-drivable.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 3,
            "expected at least 3 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check the specific captions are present as named accessibility nodes.
        for caption in ["camera views", "scene points", "image noise (px)"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name); each `labelled_by`
        // target must RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_photogrammetry_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 3,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        for caption in ["camera views", "scene points"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn noise_free_recovery_pin_from_ui_state() {
        // Mirror of the unit pin, exercised from the UI-state struct: a
        // noise-free scene recovers to ~zero reprojection error and a
        // near-perfect similarity fit to truth.
        let mut s = PhotogrammetryWorkbenchState::default();
        s.params.noise_px = 0.0;
        let res = s.run().expect("noise-free run");
        assert!(
            res.mean_reprojection_error < 1e-4,
            "noise-free mean reproj error {} should be ~0",
            res.mean_reprojection_error
        );
        assert!(
            res.recovery_rms < 1e-3,
            "noise-free recovery RMS {} should be ~0",
            res.recovery_rms
        );
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // A single camera (or zero points) must surface the error in-panel, not
        // panic.
        let mut state = PhotogrammetryWorkbenchState::default();
        state.params.num_cameras = 1;
        assert!(state.run().is_err(), "1 camera must produce Err, not panic");
        state.params.num_cameras = 4;
        state.params.num_points = 0;
        assert!(state.run().is_err(), "0 points must produce Err, not panic");
    }

    #[test]
    fn agent_bridge_photogrammetry_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "photogrammetry" }`:
        //   1. TabKind::from_id("photogrammetry") -> Some(TabKind::Photogrammetry)
        //      (plus the aliases "sfm" / "scan")
        //   2. set_workbench_flag(app, "photogrammetry", true)
        //      -> show_photogrammetry_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup (canonical + aliases).
        assert_eq!(
            TabKind::from_id("photogrammetry"),
            Some(TabKind::Photogrammetry),
            "\"photogrammetry\" must resolve to TabKind::Photogrammetry"
        );
        assert_eq!(TabKind::from_id("sfm"), Some(TabKind::Photogrammetry));
        assert_eq!(TabKind::from_id("scan"), Some(TabKind::Photogrammetry));
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(
            TabKind::from_id("  Photogrammetry  "),
            Some(TabKind::Photogrammetry)
        );

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_photogrammetry_workbench);
        set_workbench_flag(&mut app, "photogrammetry", true);
        assert!(
            app.show_photogrammetry_workbench,
            "set_workbench_flag(\"photogrammetry\", true) must set the flag"
        );
        set_workbench_flag(&mut app, "photogrammetry", false);
        assert!(!app.show_photogrammetry_workbench);
    }
}

# Third-Party Notices — valenx-photogrammetry

`valenx-photogrammetry` is licensed under `MIT OR Apache-2.0`, like the
rest of valenx. It contains **no third-party source code**. This file
records the external works that served as *method references* for the
clean-room reimplementation in this crate, and the licenses under which
those reference works are distributed.

## Method reference: COLMAP

The overall Structure-from-Motion (SfM) pipeline that this crate's
photogrammetry track targets — feature detection/description → matching →
two-view geometry → incremental mapping → bundle adjustment → sparse
reconstruction — follows the well-established design popularized by
**COLMAP**.

- Project: COLMAP — https://github.com/colmap/colmap
- Author: Johannes L. Schönberger and contributors
- License: **BSD 3-Clause** (permissive)
- Key publications:
  - J. L. Schönberger and J.-M. Frahm, "Structure-from-Motion Revisited,"
    CVPR 2016.
  - J. L. Schönberger, E. Zheng, M. Pollefeys, J.-M. Frahm, "Pixelwise
    View Selection for Unstructured Multi-View Stereo," ECCV 2016.

COLMAP is used here **only as an architectural / algorithmic reference**.
No COLMAP code is copied, linked, or vendored. The BSD-3-Clause license is
permissive and compatible with valenx's `MIT OR Apache-2.0` licensing;
this notice is provided as attribution for the method reference.

## Algorithm references (textbook computer vision)

The feature front end and the matching / two-view-verification stage are
implemented directly from the original published algorithms, all of which
are standard, widely re-implemented computer-vision methods.

### Stage 1 — feature detection & description

- **FAST** corner detector (FAST-9 variant):
  E. Rosten and T. Drummond, "Machine Learning for High-Speed Corner
  Detection," ECCV 2006.

- **ORB** oriented binary descriptor (oriented-FAST + steered/rotation-
  aware BRIEF):
  E. Rublee, V. Rabaud, K. Konolige, G. Bradski, "ORB: An efficient
  alternative to SIFT or SURF," ICCV 2011.
  - Built on BRIEF: M. Calonder, V. Lepetit, C. Strecha, P. Fua,
    "BRIEF: Binary Robust Independent Elementary Features," ECCV 2010.

  Implementation note: this crate ships oriented-FAST + **steered BRIEF**
  (a deterministic Gaussian test pattern, steered by the intensity-centroid
  orientation). It does **not** include ORB's *learned* rBRIEF test-
  selection / decorrelation step, and is not bit-compatible with OpenCV's
  ORB. See the `descriptor` module documentation for details.

### Stage 2 — descriptor matching & two-view geometric verification

- **Lowe ratio test** for distinctive nearest-neighbour matching
  (the `matching` module's keep-if-`best < ratio·second` rule):
  D. G. Lowe, "Distinctive Image Features from Scale-Invariant
  Keypoints," IJCV 60(2), 2004. (The mutual cross-check is a standard
  symmetric-consistency filter built on top.)

- **Normalized 8-point algorithm** for the fundamental matrix, with
  isotropic coordinate conditioning (the `geometry` module's estimator):
  R. I. Hartley, "In Defense of the Eight-Point Algorithm," IEEE TPAMI
  19(6), 1997. The underlying linear formulation is the eight-point method
  of H. C. Longuet-Higgins, "A computer algorithm for reconstructing a
  scene from two projections," Nature 293, 1981. General epipolar-geometry
  background: R. Hartley and A. Zisserman, *Multiple View Geometry in
  Computer Vision*, 2nd ed., 2004 (including the Sampson-distance inlier
  cost used for scoring).

- **RANSAC** robust estimation wrapping the 8-point fit:
  M. A. Fischler and R. C. Bolles, "Random Sample Consensus: A Paradigm
  for Model Fitting with Applications to Image Analysis and Automated
  Cartography," Communications of the ACM 24(6), 1981.

  Implementation note: this stage estimates the **fundamental matrix and
  the inlier correspondence set only**. The 8-point estimator is linear
  (algebraic residual), not the maximum-likelihood geometric estimate. See
  the `geometry` module documentation for details.

### Stage 3 — two-view geometry: essential matrix, relative pose, triangulation

The `twoview` module upgrades the fundamental matrix to the essential
matrix using known camera intrinsics, decomposes it into the relative
camera pose, and triangulates correspondences. All three are standard
results from the same epipolar-geometry literature.

- **Essential matrix from the fundamental matrix and intrinsics**
  (`E = K₂ᵀ F K₁`) and the **decomposition of `E` into `(R, t)`** via the
  SVD with the orthogonal `W` matrix, yielding the four candidate poses
  (Result 9.19) and the cheirality (positive-depth) test that selects the
  physically correct one (§9.6.3):
  R. Hartley and A. Zisserman, *Multiple View Geometry in Computer Vision*,
  2nd ed., Cambridge University Press, 2004, Chapter 9.

- **Linear triangulation (DLT — Direct Linear Transform)**, assembling the
  homogeneous `A X = 0` system from the two cameras' projection matrices and
  solving by SVD:
  R. Hartley and A. Zisserman, *Multiple View Geometry in Computer Vision*,
  2nd ed., §12.2. See also R. I. Hartley and P. Sturm, "Triangulation,"
  Computer Vision and Image Understanding 68(2), 1997.

  Implementation note: Stage 3 is the **calibrated** path — the camera
  intrinsics must be supplied (from EXIF, a prior calibration, or a sensible
  default). Full uncalibrated **auto-calibration** is out of scope. The
  recovered translation, and hence the triangulated points, is determined
  only **up to a global scale** (the essential matrix is homogeneous); the
  translation is returned as a unit direction. The four-fold rotation/
  translation ambiguity is resolved by the cheirality test, but the scale
  ambiguity is fundamental to two views. See the `twoview` module
  documentation for details.

### Stage 4 — incremental mapper: camera resectioning / PnP

The `pnp` module registers a new view into an existing reconstruction by
recovering its camera pose from 2D–3D correspondences (already-triangulated
scene points matched to the new image's features) — the core step of the
incremental mapper.

- **DLT camera resectioning / Perspective-n-Point (PnP)**, the linear
  Direct-Linear-Transform pose recovery: normalize image points to calibrated
  rays `x̂ = K⁻¹ x`, assemble the homogeneous `2n × 12` system from the
  cross-product constraint `x̂ᵢ × ([R | t] Xᵢ) = 0`, solve by SVD, reshape to
  the `3×4` `[R | t]`, then orthonormalize the rotation block (nearest proper
  rotation by SVD) and rescale the translation for metric consistency:
  R. Hartley and A. Zisserman, *Multiple View Geometry in Computer Vision*,
  2nd ed., Cambridge University Press, 2004, Chapter 7 (camera resectioning,
  the DLT algorithm). General PnP background: the calibrated minimal case is
  P3P; the linear DLT used here needs `n ≥ 6` points.

- **RANSAC** robust estimation wrapping the linear PnP fit, scored by
  **reprojection error** (the same Fischler & Bolles 1981 paradigm cited for
  Stage 2):
  M. A. Fischler and R. C. Bolles, "Random Sample Consensus," Communications
  of the ACM 24(6), 1981.

  Implementation note: this stage solves PnP by the **linear DLT** — an
  algebraic, not maximum-likelihood, estimate. It is a fast, accurate
  *initializer* (essentially exact on clean correspondences) for the later
  bundle-adjustment stage, but is below EPnP / iterative (Levenberg–Marquardt)
  accuracy on noisy data. It requires at least **six** correspondences, and a
  **coplanar** (or collinear) point configuration is degenerate for the DLT (a
  planar scene needs a homography-based pose, not implemented here). Because
  the 3-D points are supplied in a fixed world frame, the recovered pose is at
  that frame's **metric scale** — unlike Stage 3, there is no free scale
  factor. See the `pnp` module documentation for details.

### Stage 5 — bundle adjustment (joint nonlinear refinement)

The `bundle` module is the final solver stage: it jointly refines all camera
poses and all 3-D points to minimize the total reprojection error, by a dense
Levenberg–Marquardt optimization. These are standard results from the same
multiple-view-geometry / optimization literature.

- **Bundle adjustment** as the maximum-likelihood refinement of structure and
  motion under Gaussian reprojection noise, and its solution by the
  **Levenberg–Marquardt** damped-Gauss–Newton method:
  R. Hartley and A. Zisserman, *Multiple View Geometry in Computer Vision*,
  2nd ed., Cambridge University Press, 2004, Appendix 6 (iterative estimation
  methods), and
  B. Triggs, P. F. McLauchlan, R. I. Hartley, A. W. Fitzgibbon, "Bundle
  Adjustment — A Modern Synthesis," in *Vision Algorithms: Theory and
  Practice* (LNCS 1883), Springer, 2000.
  The Levenberg–Marquardt algorithm itself: K. Levenberg, "A Method for the
  Solution of Certain Non-Linear Problems in Least Squares," *Quarterly of
  Applied Mathematics* 2(2), 1944; D. W. Marquardt, "An Algorithm for
  Least-Squares Estimation of Nonlinear Parameters," *Journal of SIAM* 11(2),
  1963.

- **Angle-axis (Rodrigues) rotation parametrization** used for the 6-DOF
  camera blocks (the `exp` / `log` maps between a rotation matrix and its
  axis-angle 3-vector): O. Rodrigues, "Des lois géométriques qui régissent les
  déplacements d'un système solide…," *Journal de Mathématiques Pures et
  Appliquées* 5, 1840; standard modern treatment in Hartley & Zisserman
  (above) and in the bundle-adjustment literature.

  Implementation note: this is a **dense** Levenberg–Marquardt — it forms and
  factorizes the full normal matrix `JᵀJ`, costing roughly
  `O((6·#cameras + 3·#points)³)` per iteration. It is appropriate for the
  small problems this crate targets (tens of cameras, hundreds of points) and
  does **not** implement the sparse **Schur-complement** reduction that
  production bundle adjustment uses to scale to thousands of cameras (a
  deliberate future extension). The Jacobian is built by **numerical finite
  differences** (an analytic Jacobian is a future optimization); the camera
  intrinsics `K` are **shared and held fixed** (per-camera / refined intrinsics
  are a future extension); the **gauge** is fixed by holding camera 0's pose,
  which removes the 6-DOF rigid freedom but not the global **scale** of a
  two-view-seeded reconstruction. Levenberg–Marquardt is a **local** optimizer
  and requires the kind of initialization that Stages 3–4 provide. See the
  `bundle` module documentation for details.

These publications describe the methods; the Rust code in this crate is an
independent implementation. The matching, verification, two-view-geometry,
camera-resectioning (PnP), and bundle-adjustment code is written from the
mathematics in the references above — **no COLMAP or OpenCV source is used or
copied.**

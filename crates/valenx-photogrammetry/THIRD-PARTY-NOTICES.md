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
  the inlier correspondence set only**. The essential matrix, relative
  camera pose `(R, t)`, and triangulation (which need camera intrinsics)
  are a later stage. The 8-point estimator is linear (algebraic residual),
  not the maximum-likelihood geometric estimate. See the `geometry` module
  documentation for details.

These publications describe the methods; the Rust code in this crate is an
independent implementation. The matching and verification code is written
from the mathematics in the references above — **no COLMAP or OpenCV source
is used or copied.**

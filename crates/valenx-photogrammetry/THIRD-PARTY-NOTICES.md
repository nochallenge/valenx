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

The Stage-1 feature front end is implemented directly from the original
published algorithms, both of which are standard, widely re-implemented
computer-vision methods:

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

These publications describe the methods; the Rust code in this crate is an
independent implementation.

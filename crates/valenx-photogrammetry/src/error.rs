//! Error type for the photogrammetry crate.

use thiserror::Error;

/// Shorthand for `Result<T, PhotogrammetryError>`.
pub type Result<T> = core::result::Result<T, PhotogrammetryError>;

/// Anything that can go wrong constructing or operating on photogrammetry
/// inputs.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PhotogrammetryError {
    /// An image dimension was zero. A `GrayImage` must have a strictly
    /// positive width and height so that pixel indexing is well defined.
    #[error("invalid image dimensions: width = {width}, height = {height} (both must be > 0)")]
    ZeroDimension {
        /// The supplied width.
        width: usize,
        /// The supplied height.
        height: usize,
    },

    /// The length of the pixel buffer did not equal `width * height`.
    /// A `GrayImage` is row-major with exactly one byte per pixel.
    #[error("pixel buffer length {actual} does not match width * height = {expected} ({width}x{height})")]
    PixelCountMismatch {
        /// The expected length, `width * height`.
        expected: usize,
        /// The actual buffer length that was supplied.
        actual: usize,
        /// The supplied width.
        width: usize,
        /// The supplied height.
        height: usize,
    },

    /// `width * height` overflowed `usize`. Guards the validated
    /// constructor against pathological dimensions on 32-bit targets.
    #[error("image dimensions {width}x{height} overflow usize when multiplied")]
    DimensionOverflow {
        /// The supplied width.
        width: usize,
        /// The supplied height.
        height: usize,
    },

    /// A reconstruction-quality operation was asked to work on a reconstruction
    /// with nothing to measure or export — no registered cameras *and* no 3-D
    /// points (or, for the exporters, an empty point cloud / no registered
    /// cameras). Reported loudly rather than returning vacuous zero metrics or
    /// an empty file. See the `quality` module.
    #[error("empty reconstruction: nothing to measure, filter, or export")]
    EmptyReconstruction,

    /// A non-finite (`NaN` or `±∞`) value was encountered where a finite number
    /// is required — a stored 3-D point coordinate, or a non-finite outlier
    /// rule parameter. Guards the quality metrics / filtering against silently
    /// propagating `NaN`.
    #[error("non-finite value (NaN or infinity) where a finite number is required")]
    NonFiniteValue,

    /// A negative threshold was supplied to the reprojection-outlier filter
    /// (the pixel threshold or the robust `k·MAD` multiplier must be `>= 0`).
    #[error("negative threshold supplied to the outlier filter (must be >= 0)")]
    NegativeThreshold,

    /// A per-point colour array was supplied to the PLY exporter whose length
    /// did not match the number of points.
    #[error("colour count {actual} does not match the {expected} points in the cloud")]
    ColorCountMismatch {
        /// The expected length (number of points in the cloud).
        expected: usize,
        /// The actual colour-array length supplied.
        actual: usize,
    },

    /// An I/O error occurred writing an export file (PLY / camera poses). The
    /// payload is the underlying error's display string (kept as `String` so the
    /// error type stays `Clone`/`PartialEq`, which `std::io::Error` is not).
    #[error("I/O error writing export file: {0}")]
    Io(String),
}

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
}

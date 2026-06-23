//! 8-bit grayscale image buffer.
//!
//! [`GrayImage`] is the single input type for Stage 1 of the
//! structure-from-motion pipeline. It deliberately avoids the `image`
//! crate: it is just a row-major `Vec<u8>` of intensities plus its
//! dimensions, validated on construction so the rest of the crate can
//! index it without bounds anxiety.

use crate::error::{PhotogrammetryError, Result};

/// A row-major 8-bit grayscale image: one `u8` intensity per pixel.
///
/// Pixels are stored row by row, so the intensity at integer coordinate
/// `(x, y)` lives at index `y * width + x`. The origin `(0, 0)` is the
/// top-left corner, `x` increases to the right and `y` increases
/// downward — the usual image-processing convention.
///
/// Construct one with [`GrayImage::new`], which validates that the
/// dimensions are non-zero and that the buffer length is exactly
/// `width * height`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GrayImage {
    /// Image width in pixels (number of columns). Always `> 0`.
    pub width: usize,
    /// Image height in pixels (number of rows). Always `> 0`.
    pub height: usize,
    /// Row-major intensity buffer, length `width * height`.
    pub pixels: Vec<u8>,
}

impl GrayImage {
    /// Build a validated grayscale image from its dimensions and a
    /// row-major intensity buffer.
    ///
    /// # Errors
    ///
    /// - [`PhotogrammetryError::ZeroDimension`] if `width` or `height` is
    ///   zero.
    /// - [`PhotogrammetryError::DimensionOverflow`] if `width * height`
    ///   overflows `usize`.
    /// - [`PhotogrammetryError::PixelCountMismatch`] if
    ///   `pixels.len() != width * height`.
    pub fn new(width: usize, height: usize, pixels: Vec<u8>) -> Result<Self> {
        if width == 0 || height == 0 {
            return Err(PhotogrammetryError::ZeroDimension { width, height });
        }
        let expected = width
            .checked_mul(height)
            .ok_or(PhotogrammetryError::DimensionOverflow { width, height })?;
        if pixels.len() != expected {
            return Err(PhotogrammetryError::PixelCountMismatch {
                expected,
                actual: pixels.len(),
                width,
                height,
            });
        }
        Ok(Self {
            width,
            height,
            pixels,
        })
    }

    /// Total number of pixels (`width * height`). Never overflows because
    /// [`GrayImage::new`] already rejected dimensions whose product does.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.width * self.height
    }

    /// Always `false`: a validated [`GrayImage`] has non-zero dimensions
    /// and therefore at least one pixel. Provided for API completeness
    /// alongside [`GrayImage::len`].
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Intensity at `(x, y)` without bounds checking on the caller's part.
    ///
    /// # Panics
    ///
    /// Panics if `x >= width` or `y >= height` (the underlying slice
    /// index is out of range). Internal hot loops that have already
    /// confirmed they stay inside a border use this; external callers
    /// should prefer [`GrayImage::get`].
    #[inline]
    #[must_use]
    pub fn at(&self, x: usize, y: usize) -> u8 {
        debug_assert!(
            x < self.width && y < self.height,
            "GrayImage::at out of bounds"
        );
        self.pixels[y * self.width + x]
    }

    /// Intensity at `(x, y)`, or `None` if the coordinate lies outside the
    /// image. The bounds-checked counterpart of [`GrayImage::at`].
    #[inline]
    #[must_use]
    pub fn get(&self, x: usize, y: usize) -> Option<u8> {
        if x < self.width && y < self.height {
            Some(self.pixels[y * self.width + x])
        } else {
            None
        }
    }
}

//! Print-bed-layout error taxonomy.

use thiserror::Error;

/// Errors raised by the print-bed crate.
#[derive(Debug, Error)]
pub enum PrintBedError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// A part is bigger than the printer bed and can't be placed.
    #[error("part `{name}` ({w:.1}x{h:.1} mm) does not fit in bed ({bw:.1}x{bh:.1} mm)")]
    PartTooLarge {
        /// Part name.
        name: String,
        /// Part footprint width (mm).
        w: f64,
        /// Part footprint depth (mm).
        h: f64,
        /// Bed width (mm).
        bw: f64,
        /// Bed depth (mm).
        bh: f64,
    },

    /// I/O failure during bundle export.
    #[error("io ({path}): {reason}")]
    Io {
        /// Path.
        path: String,
        /// Reason.
        reason: String,
    },

    /// Mesh write failed during bundle export.
    #[error("stl-write ({path}): {reason}")]
    StlWrite {
        /// Path.
        path: String,
        /// Reason.
        reason: String,
    },
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Algorithm / packing.
    Algorithm,
    /// I/O.
    Io,
}

impl PrintBedError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            PrintBedError::BadParameter { .. } => "print_bed.bad_parameter",
            PrintBedError::PartTooLarge { .. } => "print_bed.part_too_large",
            PrintBedError::Io { .. } => "print_bed.io",
            PrintBedError::StlWrite { .. } => "print_bed.stl_write",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PrintBedError::BadParameter { .. } | PrintBedError::PartTooLarge { .. } => {
                ErrorCategory::Input
            }
            PrintBedError::Io { .. } | PrintBedError::StlWrite { .. } => ErrorCategory::Io,
        }
    }
}

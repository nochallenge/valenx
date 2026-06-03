//! Typed error taxonomy for the CAMotics-style simulator.

use thiserror::Error;

/// Errors raised by `valenx-camotics-sim`.
#[derive(Debug, Error)]
pub enum CamoticsError {
    /// Bad parameter (invalid frame count, NaN tool radius, etc.).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Frame index out of range.
    #[error("frame index {0} out of range (n_frames = {1})")]
    FrameOutOfRange(usize, usize),

    /// Underlying voxel / cam error forwarded from `valenx-cam`.
    #[error("cam: {0}")]
    Cam(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input (parameters, frame indices).
    Input,
    /// Underlying CAM kernel failure.
    Algorithm,
}

impl CamoticsError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            CamoticsError::BadParameter { .. } => "camotics.bad_parameter",
            CamoticsError::FrameOutOfRange(_, _) => "camotics.frame_out_of_range",
            CamoticsError::Cam(_) => "camotics.cam",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CamoticsError::BadParameter { .. } | CamoticsError::FrameOutOfRange(_, _) => {
                ErrorCategory::Input
            }
            CamoticsError::Cam(_) => ErrorCategory::Algorithm,
        }
    }
}

//! The crate-wide error type.

use thiserror::Error;

/// Errors produced by the autonomy V&V framework.
///
/// Like `valenx-sensors`, every constructor and runner validates its inputs and
/// returns a [`Result`] rather than panicking or silently producing a `NaN`, so
/// a bad configuration (an empty suite, a non-finite scenario parameter, a
/// requirement evaluated against a trace it cannot apply to) is a recoverable
/// error caught early — *fail loud, fail early*.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum VnvError {
    /// A configuration parameter was out of its valid range or otherwise
    /// unusable — e.g. an empty [`crate::ScenarioSuite`], a non-positive time
    /// step, an empty command sequence, or a sweep grid with a zero-length axis.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    /// A supplied value was non-finite (`NaN` / `±∞`) where the framework
    /// requires finite numbers — e.g. a scenario parameter, a requirement
    /// threshold, or an axis-of-variation value.
    #[error("non-finite value: {0}")]
    NonFinite(String),

    /// A [`crate::Requirement`] could not be applied to the [`crate::Trace`] it
    /// was given — e.g. a requirement that inspects LiDAR returns evaluated
    /// against a trace whose frames carry no LiDAR, or an empty trace. This is a
    /// *mismatch* between what the requirement needs and what the run produced,
    /// distinct from an ordinary pass/fail.
    #[error("requirement/trace mismatch: {0}")]
    RequirementMismatch(String),

    /// An error bubbled up from the underlying `valenx-sensors` harness while
    /// driving a scenario (a bad `dt`, a non-finite command, …).
    #[error("sensor harness error: {0}")]
    Harness(#[from] valenx_sensors::SensorError),
}

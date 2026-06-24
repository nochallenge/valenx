//! [`VnvReport`] — the per-requirement verdict for one run, and [`evaluate`].
//!
//! [`evaluate`] scores a [`RequirementSet`] against a single [`Trace`],
//! producing a [`VnvReport`]: one [`RequirementOutcome`] per requirement plus an
//! overall pass (the logical AND of all of them — every requirement must hold).

use crate::error::VnvError;
use crate::requirement::{RequirementOutcome, RequirementSet};
use crate::trace::Trace;

/// The verdict for one [`Trace`] against one [`RequirementSet`].
#[derive(Debug, Clone, PartialEq)]
pub struct VnvReport {
    /// The scenario name this report is for.
    pub scenario: String,
    /// One outcome per requirement, in the requirement set's order.
    pub outcomes: Vec<RequirementOutcome>,
    /// Overall pass: every requirement passed. An empty requirement set passes
    /// vacuously (`true`) — recorded honestly so the reader knows nothing was
    /// actually checked.
    pub overall_pass: bool,
}

impl VnvReport {
    /// The number of requirements that passed.
    #[must_use]
    pub fn num_passed(&self) -> usize {
        self.outcomes.iter().filter(|o| o.pass).count()
    }

    /// The number of requirements that failed.
    #[must_use]
    pub fn num_failed(&self) -> usize {
        self.outcomes.iter().filter(|o| !o.pass).count()
    }

    /// The worst (minimum) margin across all requirements, or `None` if the
    /// requirement set was empty. With the uniform sign convention this is the
    /// closest-to-failing (or deepest-failing) requirement.
    #[must_use]
    pub fn worst_margin(&self) -> Option<f64> {
        self.outcomes
            .iter()
            .map(|o| o.margin)
            .fold(None, |acc, m| Some(acc.map_or(m, |a: f64| a.min(m))))
    }
}

/// Evaluate a [`RequirementSet`] against a [`Trace`], producing a [`VnvReport`].
///
/// Each requirement is evaluated independently against the trace (which carries
/// its own scene, so clearance/collision checks are self-contained). The overall
/// verdict is the AND of all per-requirement passes.
///
/// # Errors
/// Propagates the first requirement evaluation error — a parameter validation
/// failure ([`VnvError::InvalidConfig`] / [`VnvError::NonFinite`]) or a
/// trace/requirement mismatch ([`VnvError::RequirementMismatch`], e.g. a
/// `DetectByTime` over a LiDAR-less trace, or any requirement over an empty
/// trace). A mismatch is a *setup* error worth failing loud on, distinct from an
/// ordinary requirement *failure* (which is reported as `pass = false`).
pub fn evaluate(requirements: &RequirementSet, trace: &Trace) -> Result<VnvReport, VnvError> {
    let mut outcomes = Vec::with_capacity(requirements.len());
    for r in &requirements.requirements {
        outcomes.push(r.evaluate(trace)?);
    }
    let overall_pass = outcomes.iter().all(|o| o.pass);
    Ok(VnvReport {
        scenario: trace.scenario.clone(),
        outcomes,
        overall_pass,
    })
}

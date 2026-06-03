//! Time axis for fields.
//!
//! A result may be **steady** (one snapshot), **iteration-indexed**
//! (SIMPLE / coupled loops count without physical time), or
//! **transient** (physical time carrying its own unit). The enum
//! below covers all three uniformly so `FieldCatalog` can index
//! series the same way regardless of physics.

use serde::{Deserialize, Serialize};

use crate::units::{Units, SECOND};

/// A time key identifying a snapshot within a field's time series.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum TimeKey {
    /// A single, time-independent snapshot.
    Steady,
    /// A dimensionless iteration count (SIMPLE outer loop, coupled
    /// iteration, etc.).
    Iteration(u64),
    /// A physical time value with units.
    Time { value: f64, units: Units },
}

impl TimeKey {
    /// A physical time in seconds — shorthand for
    /// `TimeKey::Time { value, units: SECOND }`.
    pub fn seconds(value: f64) -> Self {
        Self::Time {
            value,
            units: SECOND,
        }
    }
}

impl PartialEq for TimeKey {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (TimeKey::Steady, TimeKey::Steady) => true,
            (TimeKey::Iteration(a), TimeKey::Iteration(b)) => a == b,
            (
                TimeKey::Time {
                    value: a,
                    units: ua,
                },
                TimeKey::Time {
                    value: b,
                    units: ub,
                },
            ) => a == b && ua.is_compatible(ub),
            _ => false,
        }
    }
}

impl Eq for TimeKey {}

impl PartialOrd for TimeKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        use std::cmp::Ordering::*;
        match (self, other) {
            (TimeKey::Steady, TimeKey::Steady) => Some(Equal),
            (TimeKey::Iteration(a), TimeKey::Iteration(b)) => a.partial_cmp(b),
            (
                TimeKey::Time {
                    value: a,
                    units: ua,
                },
                TimeKey::Time {
                    value: b,
                    units: ub,
                },
            ) => {
                // Only comparable if dimensions match.
                if !ua.is_compatible(ub) {
                    return None;
                }
                // Normalise `b` into `a`'s scale so we compare the
                // same physical quantity, not raw display numbers.
                let factor = ub.convert_factor(ua)?;
                a.partial_cmp(&(b * factor))
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{MILLISECOND, SECOND};

    #[test]
    fn seconds_shorthand() {
        let t = TimeKey::seconds(1.5);
        assert_eq!(
            t,
            TimeKey::Time {
                value: 1.5,
                units: SECOND
            }
        );
    }

    #[test]
    fn iteration_equality() {
        assert_eq!(TimeKey::Iteration(42), TimeKey::Iteration(42));
        assert_ne!(TimeKey::Iteration(1), TimeKey::Iteration(2));
    }

    #[test]
    fn time_compares_across_scales() {
        let a = TimeKey::Time {
            value: 1.0,
            units: SECOND,
        };
        let b = TimeKey::Time {
            value: 1_000.0,
            units: MILLISECOND,
        };
        // 1 s == 1000 ms physically, so `partial_cmp` returns Equal.
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Equal));
    }

    #[test]
    fn incompatible_times_not_comparable() {
        let a = TimeKey::Time {
            value: 1.0,
            units: SECOND,
        };
        let b = TimeKey::Iteration(1);
        assert!(a.partial_cmp(&b).is_none());
    }
}

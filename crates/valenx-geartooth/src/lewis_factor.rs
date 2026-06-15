//! Lewis form factor `Y` for involute spur-gear teeth.
//!
//! ## Model
//!
//! The Lewis form factor is a dimensionless number that captures the
//! tooth profile geometry in the Lewis bending equation. It is derived
//! by treating the tooth as a cantilever beam of parabolic (uniform
//! bending stress) outline inscribed in the tooth, and tabulated as a
//! function of the number of teeth `N`. As `N` grows the tooth becomes
//! wider and stiffer at the root, so `Y` increases monotonically toward
//! the rack limit (`N -> infinity`).
//!
//! The table below is the standard one for **20-degree full-depth
//! involute** teeth measured at the **highest point of single-tooth
//! contact** as reproduced in Shigley's *Mechanical Engineering
//! Design*. Values for `N` between table rows are obtained by linear
//! interpolation; values above the largest tabulated `N` saturate at
//! the rack value.
//!
//! Note: this `Y` is the *load-applied-at-the-tip*-equivalent profile
//! factor as Shigley tabulates it for use with `sigma = Wt / (F m Y)`
//! (SI, module `m`). It is dimensionless. Do not confuse it with the
//! older Lewis tip-load factor `y = Y / pi`.

use crate::error::GearToothError;

/// One row of the Lewis form-factor table: a tooth count paired with
/// the corresponding form factor `Y`.
#[derive(Copy, Clone, Debug, PartialEq)]
struct LewisRow {
    /// Number of teeth.
    teeth: u32,
    /// Lewis form factor `Y` (dimensionless).
    y: f64,
}

/// Lewis form factors for 20-degree full-depth involute teeth, load at
/// the highest point of single-tooth contact (Shigley, Table 14-2).
///
/// Strictly increasing in both columns, which the unit tests assert.
const LEWIS_TABLE: &[LewisRow] = &[
    LewisRow {
        teeth: 12,
        y: 0.245,
    },
    LewisRow {
        teeth: 13,
        y: 0.261,
    },
    LewisRow {
        teeth: 14,
        y: 0.277,
    },
    LewisRow {
        teeth: 15,
        y: 0.290,
    },
    LewisRow {
        teeth: 16,
        y: 0.296,
    },
    LewisRow {
        teeth: 17,
        y: 0.303,
    },
    LewisRow {
        teeth: 18,
        y: 0.309,
    },
    LewisRow {
        teeth: 19,
        y: 0.314,
    },
    LewisRow {
        teeth: 20,
        y: 0.322,
    },
    LewisRow {
        teeth: 21,
        y: 0.328,
    },
    LewisRow {
        teeth: 22,
        y: 0.331,
    },
    LewisRow {
        teeth: 24,
        y: 0.337,
    },
    LewisRow {
        teeth: 26,
        y: 0.346,
    },
    LewisRow {
        teeth: 28,
        y: 0.353,
    },
    LewisRow {
        teeth: 30,
        y: 0.359,
    },
    LewisRow {
        teeth: 34,
        y: 0.371,
    },
    LewisRow {
        teeth: 38,
        y: 0.384,
    },
    LewisRow {
        teeth: 43,
        y: 0.397,
    },
    LewisRow {
        teeth: 50,
        y: 0.409,
    },
    LewisRow {
        teeth: 60,
        y: 0.422,
    },
    LewisRow {
        teeth: 75,
        y: 0.435,
    },
    LewisRow {
        teeth: 100,
        y: 0.447,
    },
    LewisRow {
        teeth: 150,
        y: 0.460,
    },
    LewisRow {
        teeth: 300,
        y: 0.472,
    },
    LewisRow {
        teeth: 400,
        y: 0.480,
    },
];

/// Smallest tooth count present in the embedded Lewis table. Querying
/// below this is rejected because such pinions undercut without a
/// profile shift and the tabulated factor no longer applies.
pub const MIN_TABULATED_TEETH: u32 = 12;

/// Largest tooth count present in the embedded Lewis table. Above this
/// the factor is taken as the rack value (the final table entry).
pub const MAX_TABULATED_TEETH: u32 = 400;

/// Look up the Lewis form factor `Y` for a 20-degree full-depth
/// involute tooth with `teeth` teeth.
///
/// Behaviour:
///
/// - An exact table row returns its value directly.
/// - A count between two rows returns the linear interpolation.
/// - A count above [`MAX_TABULATED_TEETH`] saturates at the rack value.
/// - A count below [`MIN_TABULATED_TEETH`] returns
///   [`GearToothError::OutOfDomain`].
///
/// # Errors
///
/// Returns [`GearToothError::OutOfDomain`] when `teeth` is below the
/// minimum tabulated count.
pub fn lewis_form_factor(teeth: u32) -> Result<f64, GearToothError> {
    if teeth < MIN_TABULATED_TEETH {
        return Err(GearToothError::OutOfDomain(format!(
            "tooth count {teeth} below minimum tabulated {MIN_TABULATED_TEETH} \
             (such pinions undercut; apply a profile shift)"
        )));
    }

    // Saturate above the table.
    let last = LEWIS_TABLE[LEWIS_TABLE.len() - 1];
    if teeth >= last.teeth {
        return Ok(last.y);
    }

    // Find the bracketing rows. The table is sorted ascending by teeth.
    for window in LEWIS_TABLE.windows(2) {
        let lo = window[0];
        let hi = window[1];
        if teeth == lo.teeth {
            return Ok(lo.y);
        }
        if teeth > lo.teeth && teeth < hi.teeth {
            let span = (hi.teeth - lo.teeth) as f64;
            let frac = (teeth - lo.teeth) as f64 / span;
            return Ok(lo.y + frac * (hi.y - lo.y));
        }
    }

    // Unreachable: teeth is in [first.teeth, last.teeth) so the loop
    // above always returns. Fall back to the first row defensively.
    Ok(LEWIS_TABLE[0].y)
}

/// Number of rows in the embedded table. Exposed mainly for tests and
/// introspection.
pub fn table_len() -> usize {
    LEWIS_TABLE.len()
}

/// The `(teeth, Y)` pair at table row `index`, or `None` if out of
/// range. Rows are ordered by ascending tooth count.
pub fn table_row(index: usize) -> Option<(u32, f64)> {
    LEWIS_TABLE.get(index).map(|r| (r.teeth, r.y))
}

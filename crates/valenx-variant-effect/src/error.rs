//! Error taxonomy for `valenx-variant-effect`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, VariantError>`](crate::Result). The variants separate the
//! three failure modes a variant-effect caller cares about:
//!
//! 1. The variant string did not parse ([`VariantError::Parse`]).
//! 2. The variant does not match the reference sequence — a wrong
//!    wild-type residue ([`VariantError::WildTypeMismatch`]), a position
//!    past the end of the sequence ([`VariantError::PositionOutOfRange`]),
//!    or an illegal residue token ([`VariantError::InvalidResidue`]).
//! 3. The underlying [`valenx_bioseq`] layer rejected an operation
//!    ([`VariantError::Bioseq`]) — e.g. building or translating a `Seq`.

use thiserror::Error;

/// Errors produced by `valenx-variant-effect`.
///
/// Derives [`thiserror::Error`]; each variant carries a human-readable
/// `Display` message via its `#[error(...)]` attribute.
///
/// Marked `#[non_exhaustive]`: more failure modes may be added as this
/// orchestration layer grows, so downstream matches must include a
/// wildcard arm.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum VariantError {
    /// A variant string failed to parse. `input` is the offending text
    /// (verbatim); `reason` is a human-readable explanation.
    #[error("could not parse variant `{input}`: {reason}")]
    Parse {
        /// The variant string that failed to parse.
        input: String,
        /// Why parsing failed.
        reason: String,
    },

    /// The reference residue at `position` did not match the wild-type
    /// residue named by the variant. Positions are 1-based, matching the
    /// HGVS convention used by the variant syntax.
    #[error(
        "wild-type mismatch at position {position}: variant expects `{expected}` \
         but the reference has `{found}`"
    )]
    WildTypeMismatch {
        /// 1-based position of the mismatch.
        position: usize,
        /// Residue the variant declared as wild-type.
        expected: char,
        /// Residue actually present in the reference at that position.
        found: char,
    },

    /// The variant references a 1-based `position` that lies past the end
    /// of a reference sequence of length `len`.
    #[error("position {position} is out of range for a sequence of length {len}")]
    PositionOutOfRange {
        /// 1-based position requested by the variant.
        position: usize,
        /// Length of the reference sequence.
        len: usize,
    },

    /// A residue token was not a recognised single-letter amino acid or
    /// nucleotide (or a 3-letter amino-acid code outside the standard 20).
    #[error("invalid residue `{residue}`")]
    InvalidResidue {
        /// The unrecognised residue character.
        residue: char,
    },

    /// The underlying [`valenx_bioseq`] layer rejected an operation. The
    /// string is the wrapped error's `Display` text.
    #[error("bioseq error: {0}")]
    Bioseq(String),
}

impl VariantError {
    /// Convenience constructor for [`VariantError::Parse`].
    pub fn parse(input: impl Into<String>, reason: impl Into<String>) -> Self {
        VariantError::Parse {
            input: input.into(),
            reason: reason.into(),
        }
    }
}

/// Wraps any [`valenx_bioseq`] error into [`VariantError::Bioseq`] by
/// capturing its `Display` text. `valenx-bioseq`'s `BioseqError` does not
/// derive `thiserror::Error`, so we bridge through its `Display` impl
/// rather than `#[from]`.
impl From<valenx_bioseq::BioseqError> for VariantError {
    fn from(e: valenx_bioseq::BioseqError) -> Self {
        VariantError::Bioseq(e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_display_contains_input_and_reason() {
        let e = VariantError::parse("p.X1Y", "unknown amino acid");
        let msg = e.to_string();
        assert!(msg.contains("p.X1Y"), "got: {msg}");
        assert!(msg.contains("unknown amino acid"), "got: {msg}");
    }

    #[test]
    fn wild_type_mismatch_display_names_residues_and_position() {
        let e = VariantError::WildTypeMismatch {
            position: 273,
            expected: 'R',
            found: 'C',
        };
        let msg = e.to_string();
        assert!(msg.contains("273"), "got: {msg}");
        assert!(msg.contains('R'), "got: {msg}");
        assert!(msg.contains('C'), "got: {msg}");
    }

    #[test]
    fn position_out_of_range_display() {
        let e = VariantError::PositionOutOfRange {
            position: 999,
            len: 100,
        };
        let msg = e.to_string();
        assert!(msg.contains("999"), "got: {msg}");
        assert!(msg.contains("100"), "got: {msg}");
    }

    #[test]
    fn invalid_residue_display() {
        let e = VariantError::InvalidResidue { residue: 'Z' };
        assert!(e.to_string().contains('Z'));
    }

    #[test]
    fn bioseq_error_is_wrapped_via_display() {
        let be = valenx_bioseq::BioseqError::alphabet('U', "DNA");
        let ve: VariantError = be.clone().into();
        match &ve {
            VariantError::Bioseq(s) => {
                // The wrapped Display text is preserved.
                assert_eq!(s, &be.to_string());
                assert!(s.contains('U'), "got: {s}");
            }
            other => panic!("expected Bioseq variant, got {other:?}"),
        }
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(VariantError::InvalidResidue { residue: 'J' });
        assert!(err.to_string().contains('J'));
    }
}

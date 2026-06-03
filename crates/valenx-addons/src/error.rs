//! Error taxonomy for `valenx-addons`.

use thiserror::Error;

/// Errors surfaced by the add-on manager.
#[derive(Debug, Error)]
pub enum AddonError {
    /// File-system failure (read, write, copy).
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// Manifest parse error.
    #[error("manifest parse: {0}")]
    Manifest(String),

    /// Manifest references missing fields or invalid values.
    #[error("manifest invalid: {0}")]
    InvalidManifest(String),

    /// Source directory or repository didn't exist.
    #[error("source not found: {0}")]
    SourceNotFound(String),

    /// Install destination already has an add-on with that name.
    #[error("add-on already installed: {0}")]
    AlreadyInstalled(String),

    /// Operation failed because GitHub network search is not wired in
    /// v1. Surface as a soft error so the UI can route to the manual
    /// install path instead of crashing.
    #[error("github search is not available (use manual install path)")]
    GithubSearchUnavailable,
}

impl AddonError {
    /// Stable kebab-cased identifier — mirrors the other crate error
    /// codes.
    pub fn code(&self) -> &'static str {
        match self {
            AddonError::Io(_) => "addons.io",
            AddonError::Manifest(_) => "addons.manifest",
            AddonError::InvalidManifest(_) => "addons.manifest_invalid",
            AddonError::SourceNotFound(_) => "addons.source_not_found",
            AddonError::AlreadyInstalled(_) => "addons.already_installed",
            AddonError::GithubSearchUnavailable => "addons.github_unavailable",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_a_code() {
        let cases = [
            (AddonError::Manifest("x".into()), "addons.manifest"),
            (
                AddonError::InvalidManifest("y".into()),
                "addons.manifest_invalid",
            ),
            (
                AddonError::SourceNotFound("z".into()),
                "addons.source_not_found",
            ),
            (
                AddonError::AlreadyInstalled("w".into()),
                "addons.already_installed",
            ),
            (
                AddonError::GithubSearchUnavailable,
                "addons.github_unavailable",
            ),
        ];
        for (e, expected) in cases {
            assert_eq!(e.code(), expected, "wrong code for {e:?}");
        }
    }

    #[test]
    fn io_error_has_io_code() {
        let inner = std::io::Error::new(std::io::ErrorKind::NotFound, "x");
        let e = AddonError::from(inner);
        assert_eq!(e.code(), "addons.io");
    }
}

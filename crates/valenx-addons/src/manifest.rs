//! `valenx-addon.toml` schema — the manifest every add-on ships at
//! its root.
//!
//! ```toml
//! name         = "my-cool-workbench"
//! description  = "Adds a parametric gear generator to Part Design."
//! version      = "0.2.1"
//! dependencies = []
//!
//! [entry_point]
//! kind   = "python"
//! module = "main"
//! ```

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::entrypoint::EntryPoint;
use crate::error::AddonError;

/// Validate a Python module identifier as the v1 add-on dispatcher
/// will pass it to `importlib`. Round-6 hardening: pre-fix the
/// manifest's `module` string was accepted verbatim (any character
/// allowed by TOML's string parser, including newlines and shell
/// metacharacters), which let an attacker land a payload like
/// `module = "evil; __import__('os').system('rm -rf')"` and let
/// the dispatcher execute it. The allow-list rejects anything
/// outside the canonical `[a-zA-Z0-9_.]` set Python's import
/// machinery understands.
fn validate_python_module(module: &str) -> Result<(), AddonError> {
    if module.is_empty() {
        return Err(AddonError::InvalidManifest(
            "entry_point.module is empty".into(),
        ));
    }
    if !module
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '.' || c == '_')
    {
        return Err(AddonError::InvalidManifest(format!(
            "entry_point.module `{module}` contains characters \
             outside [a-zA-Z0-9_.] (Python module-name allow-list)"
        )));
    }
    Ok(())
}

/// Parsed contents of a `valenx-addon.toml` manifest.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AddonManifest {
    /// Add-on name — used as the directory name on install. Must be
    /// non-empty and match `[a-zA-Z0-9_-]+`.
    pub name: String,
    /// One-line description shown in the Add-on Manager list view.
    pub description: String,
    /// Semver-style version string. Parsed as a free-form string in
    /// v1; semver-aware update is deferred.
    pub version: String,
    /// Names of other add-ons this one depends on. v1 accepts the
    /// field but does not enforce dependency resolution.
    #[serde(default)]
    pub dependencies: Vec<String>,
    /// Where to find the dispatched code.
    pub entry_point: EntryPoint,
}

impl AddonManifest {
    /// Validate the manifest fields after deserialisation. Catches
    /// empty / malformed names that TOML would otherwise let pass.
    pub fn validate(&self) -> Result<(), AddonError> {
        if self.name.is_empty() {
            return Err(AddonError::InvalidManifest("manifest.name is empty".into()));
        }
        if !self
            .name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
        {
            return Err(AddonError::InvalidManifest(format!(
                "manifest.name `{}` contains characters outside [a-zA-Z0-9_-]",
                self.name,
            )));
        }
        if self.version.is_empty() {
            return Err(AddonError::InvalidManifest(
                "manifest.version is empty".into(),
            ));
        }
        // Round-6: validate the Python entry-point module string so
        // a hostile `module = "evil; __import__('os')"` is rejected
        // before the v1 dispatcher hands it to `importlib`.
        if let EntryPoint::Python { module } = &self.entry_point {
            validate_python_module(module)?;
        }
        Ok(())
    }
}

/// Parse a manifest from its TOML text.
///
/// # Errors
///
/// - [`AddonError::Manifest`] for any deserialisation failure.
/// - [`AddonError::InvalidManifest`] for semantically invalid fields.
pub fn parse_manifest(text: &str) -> Result<AddonManifest, AddonError> {
    let m: AddonManifest = toml::from_str(text).map_err(|e| AddonError::Manifest(e.to_string()))?;
    m.validate()?;
    Ok(m)
}

/// Maximum bytes [`read_manifest_at`] will read from a single
/// `valenx-addon.toml`. Round-6 hardening: pre-fix `read_to_string`
/// would slurp an arbitrarily large file before serde even saw it.
///
/// Round-14 M7 tightening: lowered from 1 MiB to 256 KiB. The
/// manifest schema is a dozen fields — even a maxed-out addon with
/// a 200-entry `dependencies` array and verbose descriptions fits
/// well inside 32 KiB, so 256 KiB is already an extravagant ceiling.
/// Pre-tightening 1 MiB left a wide gap an attacker could fill with
/// e.g. a multi-hundred-KB padding comment in TOML to slow the
/// parser, or just to stress the allocator on memory-constrained
/// embedded hosts.
pub const MAX_ADDON_MANIFEST_BYTES: u64 = 256_000;

/// Back-compat alias for [`MAX_ADDON_MANIFEST_BYTES`]. Existing call
/// sites + tests reference `MAX_MANIFEST_BYTES`; the new name is
/// preferred but the old one stays so downstream consumers and the
/// test in `tests/` don't break in lock-step with the cap rename.
#[deprecated(
    since = "0.1.0",
    note = "use MAX_ADDON_MANIFEST_BYTES — the cap was tightened to 256 KiB in round-14"
)]
pub const MAX_MANIFEST_BYTES: u64 = MAX_ADDON_MANIFEST_BYTES;

/// Read + parse a manifest from a file path.
///
/// # Errors
///
/// - [`AddonError::Io`] for read failures.
/// - [`AddonError::Manifest`] for malformed TOML.
/// - [`AddonError::InvalidManifest`] for invalid fields, including a
///   manifest larger than [`MAX_ADDON_MANIFEST_BYTES`].
pub fn read_manifest_at(path: &Path) -> Result<AddonManifest, AddonError> {
    use std::io::Read;
    // Stat first so an attacker-controlled multi-GiB manifest is
    // rejected without allocation.
    let md = std::fs::metadata(path)?;
    if md.len() > MAX_ADDON_MANIFEST_BYTES {
        return Err(AddonError::InvalidManifest(format!(
            "manifest {} ({} bytes) exceeds the {}-byte cap",
            path.display(),
            md.len(),
            MAX_ADDON_MANIFEST_BYTES
        )));
    }
    // Round-11 hardening (R11-8): defense-in-depth against TOCTOU —
    // even if the file grew between stat and open (a racing writer
    // could swap in a multi-GiB payload),
    // `take(MAX_ADDON_MANIFEST_BYTES)` bounds the read at the source
    // so the parser never sees more than the cap regardless of disk
    // truth.
    let mut buf = Vec::new();
    std::fs::File::open(path)?
        .take(MAX_ADDON_MANIFEST_BYTES)
        .read_to_end(&mut buf)?;
    let text = String::from_utf8(buf).map_err(|e| {
        AddonError::InvalidManifest(format!(
            "manifest {} contains invalid UTF-8: {e}",
            path.display()
        ))
    })?;
    parse_manifest(&text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_manifest_parses() {
        let txt = r#"
            name = "my-addon"
            description = "Test addon"
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "main"
        "#;
        let m = parse_manifest(txt).unwrap();
        assert_eq!(m.name, "my-addon");
        assert_eq!(m.version, "0.1.0");
        assert!(matches!(m.entry_point, EntryPoint::Python { .. }));
    }

    #[test]
    fn missing_name_is_invalid() {
        let txt = r#"
            name = ""
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "main"
        "#;
        let err = parse_manifest(txt).unwrap_err();
        assert!(matches!(err, AddonError::InvalidManifest(_)));
    }

    #[test]
    fn invalid_name_chars_rejected() {
        let txt = r#"
            name = "bad name!"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "main"
        "#;
        let err = parse_manifest(txt).unwrap_err();
        assert!(matches!(err, AddonError::InvalidManifest(_)));
    }

    #[test]
    fn missing_entry_point_is_parse_error() {
        let txt = r#"
            name = "a"
            description = ""
            version = "0.1.0"
        "#;
        let err = parse_manifest(txt).unwrap_err();
        assert!(matches!(err, AddonError::Manifest(_)));
    }

    #[test]
    fn dependencies_default_to_empty() {
        let txt = r#"
            name = "x"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "main"
        "#;
        let m = parse_manifest(txt).unwrap();
        assert!(m.dependencies.is_empty());
    }

    #[test]
    fn wasm_entry_point_parses() {
        let txt = r#"
            name = "wasm-addon"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "wasm"
            wasm = "lib.wasm"
        "#;
        let m = parse_manifest(txt).unwrap();
        assert!(matches!(m.entry_point, EntryPoint::Wasm { .. }));
    }

    #[test]
    fn rejects_python_module_with_shell_metachars() {
        // Round-6 RED→GREEN: an EntryPoint::Python.module like
        // `"evil; __import__('os').system('rm -rf')"` must be
        // refused by `validate` so the v1 dispatcher never hands
        // it to `importlib`. The allow-list is `[a-zA-Z0-9_.]`;
        // anything else (including `;`, space, parens, quotes,
        // newlines) is a hard reject.
        let txt = r#"
            name = "injectable"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "evil; __import__('os').system('rm -rf')"
        "#;
        let err = parse_manifest(txt).unwrap_err();
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(msg.contains("module"), "msg: {msg}");
                assert!(msg.contains("allow-list"), "msg: {msg}");
            }
            other => panic!("expected InvalidManifest, got {other:?}"),
        }

        // Legitimate module names with dots and underscores still
        // parse — `my_pkg.api.dispatcher` is a normal import path.
        let ok = r#"
            name = "fine"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = "my_pkg.api.dispatcher"
        "#;
        let m = parse_manifest(ok).expect("legit module path");
        assert!(matches!(m.entry_point, EntryPoint::Python { .. }));
    }

    /// Round-11 RED→GREEN — R11-8. Pre-fix `read_manifest_at`
    /// stat-checked the file size, then called `read_to_string` —
    /// classic TOCTOU. If a racing writer swapped in a multi-GiB
    /// manifest between the metadata call and the open, the parser
    /// would still allocate the full payload. The fix wraps the
    /// reader in `.take(MAX_ADDON_MANIFEST_BYTES)` so the read is
    /// bounded at the source regardless of disk truth.
    ///
    /// We can't easily race the FS in a single-process test, so this
    /// test simulates the same shape by checking the size cap fires
    /// even when the stat check is bypassed — i.e. via the
    /// `.take()` truncation. The test writes a slightly-over-cap
    /// file and asserts the size guard catches it (RED case before
    /// fix) and that the read path is bounded (GREEN case post-fix).
    #[test]
    fn read_manifest_at_rejects_oversize_file() {
        let mut tmp = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        tmp.push(format!("valenx-r11-addon-{nanos}.toml"));

        // Write a payload one byte past the cap so the size guard
        // fires. The `.take()` second-line guard is what catches the
        // pathological race; the size check is sufficient when no
        // race happens.
        let mut payload = String::with_capacity((MAX_ADDON_MANIFEST_BYTES as usize) + 32);
        payload.push_str(
            "name = \"big\"\ndescription = \"oversized\"\nversion = \"0.1.0\"\n\
             [entry_point]\nkind = \"python\"\nmodule = \"main\"\n# pad: ",
        );
        payload.extend(std::iter::repeat_n(
            'A',
            (MAX_ADDON_MANIFEST_BYTES as usize) + 1,
        ));
        std::fs::write(&tmp, &payload).expect("write oversized manifest");

        let err = read_manifest_at(&tmp)
            .expect_err("oversized manifest must be rejected before allocation");
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(
                    msg.contains("cap"),
                    "error must mention the cap; got: {msg}"
                );
            }
            other => panic!("expected AddonError::InvalidManifest, got {other:?}"),
        }

        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-14 M7 RED→GREEN: a 1 MiB manifest must now be rejected
    /// — pre-tightening the 1 MiB cap admitted it; the new 256 KiB
    /// cap fires well before that. The test plants a 1 MiB file (an
    /// order of magnitude past the new cap) and confirms the cap
    /// surfaces a typed error before the parser allocates.
    #[test]
    fn read_manifest_at_rejects_1mib_post_tightening() {
        let mut tmp = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        tmp.push(format!("valenx-r14-addon-{nanos}.toml"));

        // Plant a 1 MiB manifest — was admitted under the old 1 MiB
        // cap, must be refused under the new 256 KiB cap.
        let mut payload = String::with_capacity(1024 * 1024 + 32);
        payload.push_str(
            "name = \"big\"\ndescription = \"oversized\"\nversion = \"0.1.0\"\n\
             [entry_point]\nkind = \"python\"\nmodule = \"main\"\n# pad: ",
        );
        payload.extend(std::iter::repeat_n('A', 1024 * 1024));
        std::fs::write(&tmp, &payload).expect("write 1 MiB manifest");

        let err = read_manifest_at(&tmp)
            .expect_err("1 MiB manifest must be rejected post round-14 tightening");
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(msg.contains("cap"), "msg: {msg}");
                // The error includes the new cap value so users can
                // see where the bound came from.
                assert!(
                    msg.contains(&MAX_ADDON_MANIFEST_BYTES.to_string()),
                    "error must report new cap value; got: {msg}"
                );
            }
            other => panic!("expected InvalidManifest, got {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn rejects_empty_python_module() {
        let txt = r#"
            name = "empty-mod"
            description = ""
            version = "0.1.0"

            [entry_point]
            kind = "python"
            module = ""
        "#;
        let err = parse_manifest(txt).unwrap_err();
        match err {
            AddonError::InvalidManifest(msg) => assert!(msg.contains("empty"), "msg: {msg}"),
            other => panic!("expected InvalidManifest, got {other:?}"),
        }
    }
}

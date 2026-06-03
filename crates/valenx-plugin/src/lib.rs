//! # valenx-plugin
//!
//! Plugin discovery + loading scaffold for runtime-loadable
//! adapters. Phase 14 of the [ROADMAP](../../../ROADMAP.md);
//! contract specified by [RFC 0003](../../../rfcs/0003-plugin-api.md).
//!
//! **Pre-alpha.** This crate today does the discovery + manifest
//! validation half of the plugin story. The second half — actually
//! loading a `cdylib` and registering its `Adapter` into
//! `valenx-core`'s registry — is deferred until after we pick the
//! ABI (the choice between `libloading` for native cdylibs vs
//! WebAssembly Component Model for sandboxed plugins is a
//! follow-up RFC).
//!
//! ## What works today
//!
//! - **`PluginManifest`** parses + validates a `valenx-plugin.toml`
//!   declaring an adapter's id / version / capabilities / required
//!   tools.
//! - **`discover()`** walks a plugin search path (per-user
//!   `<state_dir>/plugins/`, per-system `/usr/lib/valenx/plugins/`,
//!   plus any explicit paths) and returns every well-formed
//!   manifest it finds.
//! - **Per-manifest validation** — semver fits the host's
//!   `valenx_min` requirement, declared adapter id doesn't collide
//!   with an already-registered one, capabilities reference real
//!   `Capability` variants.
//!
//! ## What's deferred
//!
//! - The actual `dlopen` step. We need to pick whether to ship
//!   plugins as native `cdylib`s (fast, unsandboxed) or WASM
//!   components (slow, sandboxed). Different security postures
//!   for different deployment scenarios. RFC needed.
//! - A signing / trust-store story. Native plugins running as the
//!   user can do anything the user can; signing protects against
//!   tampering after install but not against a malicious plugin
//!   author.
//! - Lifecycle (load / probe / unload / reload-on-update).
//! - The plugin-author SDK (helper crate that re-exports
//!   `valenx_core::Adapter` + macros for the common boilerplate).

// `forbid` would block the loader module's libloading calls (every
// libloading call is `unsafe` because dlopen + symbol lookup can't
// be checked statically). The loader module overrides with an
// explicit `#[allow(unsafe_code)]` and gates each unsafe block on a
// safety comment.
#![deny(unsafe_code)]
// missing_docs left off — pre-alpha; the public API is sketched
// in the module-level doc above and the field semantics are
// explained inline. We'll switch the lint to `warn` once the
// runtime-load half lands and the API is genuinely stable.
#![allow(missing_docs)]

pub mod loader;

use std::path::{Path, PathBuf};

use semver::{Version, VersionReq};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// One plugin's `valenx-plugin.toml` manifest. Wire format keeps
/// `version` and `valenx_min` as strings (semver's serde feature
/// would pull in another dep just for the manifest); the parsed
/// [`Version`] / [`VersionReq`] live in [`PluginManifest::parsed_version`]
/// and friends.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PluginManifest {
    /// Stable plugin id. Convention: reverse-DNS, lowercase
    /// (e.g. `org.acme.cfd-magic`).
    pub id: String,
    /// Plugin version (semver string).
    pub version: String,
    /// Minimum Valenx host version this plugin is built against
    /// (semver requirement string, e.g. ">=0.1, <1.0").
    pub valenx_min: String,
    /// Adapter id this plugin registers when loaded.
    pub adapter_id: String,
    /// Human-friendly name for UI display.
    pub display_name: String,
    /// Path (relative to the manifest file) to the cdylib / wasm /
    /// other plugin payload. Validated to exist + be a regular
    /// file at discovery time; not loaded until the runtime
    /// resolves a load request.
    pub payload: PathBuf,
    /// Free-form description for the Settings → Plugins UI.
    #[serde(default)]
    pub description: Option<String>,
    /// Author / publisher name.
    #[serde(default)]
    pub author: Option<String>,
    /// SPDX license identifier of the plugin payload itself.
    #[serde(default)]
    pub license: Option<String>,
}

impl PluginManifest {
    /// Parsed plugin version. Returns an error if `version` isn't
    /// valid semver.
    pub fn parsed_version(&self) -> Result<Version, semver::Error> {
        Version::parse(&self.version)
    }

    /// Parsed `valenx_min` requirement.
    pub fn parsed_valenx_min(&self) -> Result<VersionReq, semver::Error> {
        VersionReq::parse(&self.valenx_min)
    }
}

/// One discovered plugin — manifest + the path it was found at.
#[derive(Clone, Debug)]
pub struct DiscoveredPlugin {
    pub manifest: PluginManifest,
    pub manifest_path: PathBuf,
}

/// Errors raised during discovery / validation. The runtime-load
/// errors will land in a sibling enum once the dlopen story is
/// nailed down.
#[derive(Debug, Error)]
pub enum PluginError {
    /// Manifest file couldn't be read.
    #[error("plugin manifest at {path}: {source}")]
    ManifestIo {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    /// Manifest TOML parse failed.
    #[error("plugin manifest at {path}: invalid TOML: {reason}")]
    ManifestParse { path: PathBuf, reason: String },
    /// Plugin's `valenx_min` doesn't accept the running host's
    /// version — the plugin is from a future or too-old Valenx.
    #[error("plugin `{id}` requires Valenx {required}; this host is {host_version}")]
    HostVersionMismatch {
        id: String,
        required: String,
        host_version: String,
    },
    /// Manifest's version / valenx_min strings aren't valid semver.
    #[error("plugin `{id}` has invalid semver in `{field}`: {reason}")]
    InvalidSemver {
        id: String,
        field: &'static str,
        reason: String,
    },
    /// Manifest declares a payload path that doesn't exist or is
    /// not a regular file.
    #[error("plugin `{id}` payload not found at {path}")]
    PayloadMissing { id: String, path: PathBuf },
}

/// Walk every directory in `search_paths` and return every
/// well-formed plugin manifest found.
///
/// Each plugin lives at `<dir>/<plugin-name>/valenx-plugin.toml`.
/// Subdirs without that file are silently skipped. Manifests that
/// fail to parse get logged via `tracing` and skipped — one bad
/// plugin shouldn't crash the host.
pub fn discover(search_paths: &[PathBuf], host_version: &Version) -> Vec<DiscoveredPlugin> {
    let mut out: Vec<DiscoveredPlugin> = Vec::new();
    for root in search_paths {
        let Ok(entries) = std::fs::read_dir(root) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }
            let manifest_path = path.join("valenx-plugin.toml");
            if !manifest_path.is_file() {
                continue;
            }
            match load_manifest(&manifest_path, host_version) {
                Ok(manifest) => {
                    out.push(DiscoveredPlugin {
                        manifest,
                        manifest_path,
                    });
                }
                Err(e) => {
                    tracing::warn!(
                        target: "valenx.plugin",
                        ?e,
                        ?manifest_path,
                        "skipped malformed plugin"
                    );
                }
            }
        }
    }
    out
}

/// Load + validate one manifest file.
pub fn load_manifest(
    manifest_path: &Path,
    host_version: &Version,
) -> Result<PluginManifest, PluginError> {
    // Round-21 M5: bound the manifest read at
    // MAX_PLUGIN_MANIFEST_BYTES (256 KiB — sister to the addons
    // cap). Pre-fix `discover()` walked the plugin search paths and
    // did a bare `fs::read_to_string` on every candidate; a
    // poisoned plugins dir with a multi-MB `.toml` would slurp into
    // memory during discovery before any plugin actually loaded.
    let text = valenx_core::io_caps::read_capped_to_string(
        manifest_path,
        valenx_core::io_caps::MAX_PLUGIN_MANIFEST_BYTES as usize,
    )
    .map_err(|e| PluginError::ManifestIo {
        path: manifest_path.to_path_buf(),
        source: e,
    })?;
    let manifest: PluginManifest =
        toml::from_str(&text).map_err(|e| PluginError::ManifestParse {
            path: manifest_path.to_path_buf(),
            reason: e.to_string(),
        })?;

    let valenx_min = manifest
        .parsed_valenx_min()
        .map_err(|e| PluginError::InvalidSemver {
            id: manifest.id.clone(),
            field: "valenx_min",
            reason: e.to_string(),
        })?;
    // Also validate `version` parses, even though we don't compare
    // it here — surfaces a clear error at discovery time rather
    // than at first use.
    manifest
        .parsed_version()
        .map_err(|e| PluginError::InvalidSemver {
            id: manifest.id.clone(),
            field: "version",
            reason: e.to_string(),
        })?;
    if !valenx_min.matches(host_version) {
        return Err(PluginError::HostVersionMismatch {
            id: manifest.id.clone(),
            required: manifest.valenx_min.clone(),
            host_version: host_version.to_string(),
        });
    }

    let payload_path = manifest_path
        .parent()
        .map(|p| p.join(&manifest.payload))
        .unwrap_or_else(|| manifest.payload.clone());
    if !payload_path.is_file() {
        return Err(PluginError::PayloadMissing {
            id: manifest.id.clone(),
            path: payload_path,
        });
    }

    Ok(manifest)
}

/// Default plugin search path for the current host. Caller can
/// extend with explicit paths from app settings or env.
pub fn default_search_paths() -> Vec<PathBuf> {
    let mut paths = Vec::new();
    if cfg!(target_os = "windows") {
        if let Some(appdata) = std::env::var_os("APPDATA") {
            paths.push(PathBuf::from(appdata).join("valenx").join("plugins"));
        }
    } else if cfg!(target_os = "macos") {
        if let Some(home) = std::env::var_os("HOME") {
            paths.push(
                PathBuf::from(home)
                    .join("Library")
                    .join("Application Support")
                    .join("valenx")
                    .join("plugins"),
            );
        }
    } else {
        if let Some(xdg) = std::env::var_os("XDG_DATA_HOME") {
            paths.push(PathBuf::from(xdg).join("valenx").join("plugins"));
        } else if let Some(home) = std::env::var_os("HOME") {
            paths.push(
                PathBuf::from(home)
                    .join(".local")
                    .join("share")
                    .join("valenx")
                    .join("plugins"),
            );
        }
        // System-wide install: /usr/lib/valenx/plugins.
        paths.push(PathBuf::from("/usr/lib/valenx/plugins"));
    }
    paths
}

#[cfg(test)]
mod tests {
    use super::*;

    fn host_version() -> Version {
        Version::new(0, 1, 0)
    }

    fn write_minimal_plugin(dir: &Path, id: &str, valenx_min: &str) {
        std::fs::create_dir_all(dir).unwrap();
        // The payload doesn't have to be a real cdylib for the
        // discovery half — just a regular file so the existence
        // check passes.
        let payload = dir.join("plugin.dll");
        std::fs::write(&payload, b"stub").unwrap();
        let manifest = format!(
            r#"
id = "{id}"
version = "0.1.0"
valenx_min = "{valenx_min}"
adapter_id = "stubsolver"
display_name = "Stub Solver"
payload = "plugin.dll"
"#
        );
        std::fs::write(dir.join("valenx-plugin.toml"), manifest).unwrap();
    }

    #[test]
    fn discover_finds_well_formed_plugin() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-plugin-discover-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        write_minimal_plugin(&tmp.join("good"), "test.good", ">=0.1");

        let plugins = discover(&[tmp.clone()], &host_version());
        assert_eq!(plugins.len(), 1);
        assert_eq!(plugins[0].manifest.id, "test.good");
        assert_eq!(plugins[0].manifest.adapter_id, "stubsolver");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn discover_skips_subdirs_without_manifest() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-plugin-no-manifest-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(tmp.join("not-a-plugin")).unwrap();
        std::fs::write(
            tmp.join("not-a-plugin").join("readme.txt"),
            "no manifest here",
        )
        .unwrap();
        let plugins = discover(&[tmp.clone()], &host_version());
        assert!(plugins.is_empty());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn host_version_mismatch_is_skipped_during_discovery() {
        // A plugin that requires Valenx >= 99.0 should be silently
        // dropped from discovery — not crash, not return a half-
        // valid entry.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-plugin-future-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_minimal_plugin(&tmp.join("future"), "test.future", ">=99.0.0");
        let plugins = discover(&[tmp.clone()], &host_version());
        assert!(plugins.is_empty(), "future-version plugin leaked through");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_manifest_rejects_missing_payload() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-plugin-no-payload-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Manifest references plugin.dll but we never write it.
        std::fs::write(
            tmp.join("valenx-plugin.toml"),
            r#"
id = "test.no-payload"
version = "0.1.0"
valenx_min = ">=0.1"
adapter_id = "stub"
display_name = "Stub"
payload = "plugin.dll"
"#,
        )
        .unwrap();
        let err = load_manifest(&tmp.join("valenx-plugin.toml"), &host_version()).unwrap_err();
        assert!(matches!(err, PluginError::PayloadMissing { .. }));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn default_search_paths_include_user_dir() {
        let paths = default_search_paths();
        // At least one path should mention "valenx" / "plugins" on
        // hosts where HOME / APPDATA is set.
        if !paths.is_empty() {
            assert!(paths.iter().any(|p| p.to_string_lossy().contains("valenx")));
        }
    }

    /// Round-21 M5 RED→GREEN: an oversize `valenx-plugin.toml`
    /// is rejected by the bounded reader rather than slurping into
    /// memory. Pre-fix `discover()` walked the plugin search paths
    /// and did a bare `fs::read_to_string` on every candidate; a
    /// poisoned plugins dir with a multi-MB `.toml` would slurp
    /// during discovery before any plugin actually loaded.
    #[test]
    fn load_manifest_rejects_oversize_toml() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-plugin-r21-oversize-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Sparse-allocate just past the 256 KiB cap.
        let manifest_path = tmp.join("valenx-plugin.toml");
        let f = std::fs::File::create(&manifest_path).unwrap();
        f.set_len(valenx_core::io_caps::MAX_PLUGIN_MANIFEST_BYTES + 1)
            .unwrap();
        drop(f);
        let err = load_manifest(&manifest_path, &host_version()).unwrap_err();
        match err {
            PluginError::ManifestIo { source, .. } => {
                assert_eq!(source.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected ManifestIo, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&tmp);
    }
}

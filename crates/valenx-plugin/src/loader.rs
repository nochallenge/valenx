//! `libloading`-based dlopen path for plugin payloads.
//!
//! ## Scope
//!
//! This is the v0 dlopen scaffold from RFC 0003. The loader can:
//!
//! - dlopen a plugin's payload (a `.so` / `.dylib` / `.dll`)
//! - look up the well-known `valenx_plugin_abi_version` C symbol
//! - validate the ABI version against the host's
//!   [`HOST_PLUGIN_ABI_VERSION`]
//! - hand the caller back a [`LoadedPlugin`] that owns the
//!   `libloading::Library` handle so subsequent `Symbol` lookups
//!   stay valid for the lifetime of the loaded plugin
//!
//! What's deliberately NOT here yet:
//!
//! - The full `Adapter` FFI table — `Adapter` is a Rust trait with
//!   `Box<dyn Adapter>` + `Arc<dyn Adapter>` boxed everywhere; calling
//!   it across an FFI boundary needs a careful C-vtable design that
//!   is its own follow-up RFC. The v0 ABI just exposes the version
//!   number so we can prove the dlopen path works end-to-end.
//! - Plugin signing / trust store. A loaded plugin runs as the host
//!   user; signing protects against tampering after install but not
//!   against a malicious plugin author. Out of scope for this commit.
//!
//! ## Safety story
//!
//! Every `libloading` call is `unsafe` because dlopen + symbol lookup
//! can't be checked statically. We tightly bound the unsafe scope:
//!
//! 1. `Library::new(path)` — unsafe because the OS dynamic loader
//!    runs the library's init/static-constructor code as a side
//!    effect. We mitigate by doing this only on payload paths the
//!    `PluginManifest` declared and we already verified exist on
//!    disk.
//! 2. `lib.get(symbol)` — unsafe because the returned `Symbol` is
//!    type-asserted by the caller without runtime checks. We assert
//!    a *single* well-known symbol shape (`extern "C" fn() -> u32`)
//!    so the surface is small and reviewable.
//!
//! Both calls happen inside an `unsafe { }` block immediately
//! preceded by a `// SAFETY: ...` comment that names the
//! preconditions.

#![allow(unsafe_code)]

use std::path::{Path, PathBuf};

use libloading::{Library, Symbol};
use thiserror::Error;

use crate::DiscoveredPlugin;

/// Plugin ABI version the host expects. Plugins must export a C
/// function `valenx_plugin_abi_version` returning this exact value
/// (or an error before any further calls happen).
///
/// Bumped manually whenever the FFI surface changes. Plugins built
/// against an older host see a fast `AbiMismatch` error and refuse
/// to load.
pub const HOST_PLUGIN_ABI_VERSION: u32 = 1;

/// Symbol name plugins must export. Convention: lowercase,
/// underscore-separated.
pub const ABI_VERSION_SYMBOL: &[u8] = b"valenx_plugin_abi_version";

/// One successfully-loaded plugin. Owns the `libloading::Library`
/// handle so subsequent symbol lookups stay valid as long as the
/// `LoadedPlugin` is alive.
pub struct LoadedPlugin {
    /// The underlying dlopened library. Public so callers can do
    /// further `lib.get::<...>()` lookups (when the future Adapter
    /// FFI surface lands).
    pub library: Library,
    /// ABI version the plugin reports. Always equal to
    /// [`HOST_PLUGIN_ABI_VERSION`] when load succeeds.
    pub abi_version: u32,
    /// Path the library was loaded from. Useful for diagnostics +
    /// the plugin-management UI.
    pub path: PathBuf,
}

/// Errors surfaced by the plugin loader.
#[derive(Debug, Error)]
pub enum LoaderError {
    /// The plugin's payload file (the `.dll`/`.so`/`.dylib`) is
    /// missing or unreadable.
    #[error("plugin payload at {path} not found or unreadable")]
    PayloadMissing {
        /// Resolved payload path.
        path: PathBuf,
    },
    /// `libloading::Library::new` failed to load the shared object.
    #[error("dlopen of {path} failed: {reason}")]
    DlopenFailed {
        /// Resolved payload path.
        path: PathBuf,
        /// Underlying OS reason.
        reason: String,
    },
    /// The shared object loaded but doesn't export the required
    /// `valenx_plugin_abi_version` symbol.
    #[error(
        "plugin {path} doesn't export the `valenx_plugin_abi_version` symbol — \
         either the build flags miss `crate-type = [\"cdylib\"]` or the symbol \
         is hidden by visibility rules"
    )]
    SymbolMissing {
        /// Resolved payload path.
        path: PathBuf,
    },
    /// Plugin's ABI version doesn't match what this host expects.
    #[error("plugin {path} reports ABI version {got}; this host expects {expected}")]
    AbiMismatch {
        /// Resolved payload path.
        path: PathBuf,
        /// ABI version the host was built against.
        expected: u32,
        /// ABI version the plugin advertised.
        got: u32,
    },
}

/// Resolve a [`DiscoveredPlugin`]'s payload path against the
/// manifest file's directory and return the absolute path.
pub fn resolve_payload(disc: &DiscoveredPlugin) -> PathBuf {
    if disc.manifest.payload.is_absolute() {
        disc.manifest.payload.clone()
    } else {
        disc.manifest_path
            .parent()
            .map(|p| p.join(&disc.manifest.payload))
            .unwrap_or_else(|| disc.manifest.payload.clone())
    }
}

/// Load a discovered plugin via dlopen + ABI-version check.
pub fn load(disc: &DiscoveredPlugin) -> Result<LoadedPlugin, LoaderError> {
    let path = resolve_payload(disc);
    if !path.is_file() {
        return Err(LoaderError::PayloadMissing { path });
    }
    load_from_path(&path)
}

/// Lower-level: dlopen a payload directly without going through a
/// `DiscoveredPlugin`. Used by tests + by callers that already have
/// an absolute path.
pub fn load_from_path(path: &Path) -> Result<LoadedPlugin, LoaderError> {
    if !path.is_file() {
        return Err(LoaderError::PayloadMissing {
            path: path.to_path_buf(),
        });
    }
    // SAFETY: dlopen runs the library's init / static-constructor
    // code as a side effect. We mitigate by only opening paths the
    // caller has already verified exist on disk; bad init code in a
    // user-installed plugin is a trust failure, not a memory-safety
    // bug we can guard against from here.
    let library = unsafe { Library::new(path) }.map_err(|e| LoaderError::DlopenFailed {
        path: path.to_path_buf(),
        reason: e.to_string(),
    })?;

    // SAFETY: symbol lookup is unsafe because libloading can't check
    // the function signature matches what the plugin actually exports.
    // We narrow this to a single well-known symbol with an
    // unambiguous C ABI: `extern "C" fn() -> u32`. A plugin that
    // exports the symbol with the wrong signature will dereference
    // through a mismatched fn pointer + UB; that's the trust failure
    // the future signing story closes.
    let abi_version: u32 = unsafe {
        let symbol: Symbol<unsafe extern "C" fn() -> u32> = library
            .get(ABI_VERSION_SYMBOL)
            .map_err(|_| LoaderError::SymbolMissing {
                path: path.to_path_buf(),
            })?;
        symbol()
    };

    if abi_version != HOST_PLUGIN_ABI_VERSION {
        return Err(LoaderError::AbiMismatch {
            path: path.to_path_buf(),
            expected: HOST_PLUGIN_ABI_VERSION,
            got: abi_version,
        });
    }

    Ok(LoadedPlugin {
        library,
        abi_version,
        path: path.to_path_buf(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fresh_path(label: &str, ext: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "valenx-loader-{label}-{}.{ext}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn load_from_missing_path_returns_payload_missing() {
        let p = fresh_path("missing", "so");
        match load_from_path(&p) {
            Err(LoaderError::PayloadMissing { path }) => {
                assert_eq!(path, p);
            }
            Err(other) => panic!("wrong error variant: {other:?}"),
            Ok(_) => panic!("expected Err, got Ok"),
        }
    }

    #[test]
    fn load_from_garbage_payload_returns_dlopen_failed() {
        // Write some bytes that aren't a valid shared library on
        // any platform. dlopen will refuse it; the error must
        // surface as DlopenFailed, not a panic or success.
        let p = fresh_path("garbage", "bin");
        std::fs::write(&p, b"this is definitely not a shared library").unwrap();
        match load_from_path(&p) {
            Err(LoaderError::DlopenFailed { path, .. }) => {
                assert_eq!(path, p);
            }
            Err(other) => panic!("wrong error variant: {other:?}"),
            Ok(_) => panic!("expected Err, got Ok"),
        }
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn host_abi_version_constant_is_stable() {
        // The v0 ABI is hard-coded to 1. Bumping this is a
        // breaking change for the plugin world — the constant lives
        // here so a casual edit stands out in code review.
        assert_eq!(HOST_PLUGIN_ABI_VERSION, 1);
    }

    #[test]
    fn abi_version_symbol_name_is_lowercase_underscore() {
        // Symbol-name convention is asserted via test so accidental
        // renames don't silently break every plugin in the wild.
        assert_eq!(ABI_VERSION_SYMBOL, b"valenx_plugin_abi_version");
    }

    #[test]
    fn resolve_payload_handles_relative_path() {
        let disc = DiscoveredPlugin {
            manifest: crate::PluginManifest {
                id: "test".into(),
                version: "0.1.0".into(),
                valenx_min: ">=0.1, <1.0".into(),
                adapter_id: "test".into(),
                display_name: "test".into(),
                payload: PathBuf::from("plugin.so"),
                description: None,
                author: None,
                license: None,
            },
            manifest_path: PathBuf::from("/plugins/test/valenx-plugin.toml"),
        };
        let resolved = resolve_payload(&disc);
        assert_eq!(resolved, PathBuf::from("/plugins/test/plugin.so"));
    }

    #[test]
    fn resolve_payload_passes_absolute_path_through() {
        let disc = DiscoveredPlugin {
            manifest: crate::PluginManifest {
                id: "test".into(),
                version: "0.1.0".into(),
                valenx_min: ">=0.1, <1.0".into(),
                adapter_id: "test".into(),
                display_name: "test".into(),
                payload: if cfg!(windows) {
                    PathBuf::from(r"C:\opt\plugins\test\plugin.dll")
                } else {
                    PathBuf::from("/opt/plugins/test/plugin.so")
                },
                description: None,
                author: None,
                license: None,
            },
            manifest_path: PathBuf::from("/plugins/test/valenx-plugin.toml"),
        };
        let resolved = resolve_payload(&disc);
        assert!(resolved.is_absolute());
        assert!(
            resolved.ends_with("plugin.so") || resolved.ends_with("plugin.dll"),
            "got: {resolved:?}"
        );
    }

    #[test]
    fn host_and_sdk_abi_versions_must_be_synced() {
        // Cross-crate sync: bumping one side without the other
        // would silently produce plugins the host can't load (or
        // vice versa). This test fails until both constants move
        // together.
        assert_eq!(
            HOST_PLUGIN_ABI_VERSION,
            valenx_plugin_sdk::PLUGIN_ABI_VERSION,
            "HOST_PLUGIN_ABI_VERSION ({HOST_PLUGIN_ABI_VERSION}) must equal \
             valenx_plugin_sdk::PLUGIN_ABI_VERSION ({}); bump both atomically",
            valenx_plugin_sdk::PLUGIN_ABI_VERSION,
        );
    }

    #[test]
    fn host_and_sdk_abi_symbol_names_must_be_synced() {
        // The host's loader probes for `b"valenx_plugin_abi_version"`;
        // the SDK macro emits a `#[no_mangle] pub extern "C" fn`
        // with the matching name. Drift here means the dlopen
        // succeeds but `lib.get(...)` returns SymbolMissing.
        assert_eq!(
            std::str::from_utf8(ABI_VERSION_SYMBOL).unwrap(),
            valenx_plugin_sdk::PLUGIN_ABI_VERSION_SYMBOL,
            "host loader symbol must match the SDK's emitted name"
        );
    }

    #[test]
    fn load_from_path_with_nonregular_path_errors_cleanly() {
        // Pass a directory instead of a file — same outcome as a
        // missing file: PayloadMissing.
        let dir = std::env::temp_dir();
        match load_from_path(&dir) {
            Err(LoaderError::PayloadMissing { .. }) => {}
            Err(other) => panic!("wrong error: {other:?}"),
            Ok(_) => panic!("expected Err, got Ok"),
        }
    }
}

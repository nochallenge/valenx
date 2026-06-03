//! # valenx-plugin-sdk
//!
//! Helper crate for **plugin authors**. Re-exports the
//! [`valenx_core::Adapter`] trait + supporting types so a plugin's
//! `Cargo.toml` only needs to depend on this one crate, and provides
//! the [`valenx_plugin_abi_version!`] macro that emits the
//! C-FFI symbol the host's loader looks up to validate the plugin
//! was built against a compatible version of Valenx.
//!
//! ## Why a separate crate
//!
//! The host's [`valenx-plugin`] crate carries `libloading` +
//! filesystem-discovery code that plugin authors don't want to pull
//! in. Splitting the small "vocabulary you need on the plugin side"
//! into its own crate keeps the plugin-side dep graph tight.
//!
//! Both this crate AND `valenx-plugin::loader` carry their own copy
//! of [`PLUGIN_ABI_VERSION`] (set to **`1`** today). A small
//! sync-test in each crate asserts they're equal — bumping one side
//! without the other breaks the test, forcing the change to be
//! atomic.
//!
//! ## Minimal plugin example
//!
//! `Cargo.toml`:
//!
//! ```toml
//! [package]
//! name = "my-valenx-plugin"
//!
//! [lib]
//! crate-type = ["cdylib"]
//!
//! [dependencies]
//! valenx-plugin-sdk = "0.1"
//! ```
//!
//! `src/lib.rs`:
//!
//! ```ignore
//! use valenx_plugin_sdk::valenx_plugin_abi_version;
//!
//! // Emits an `extern "C" fn valenx_plugin_abi_version() -> u32`
//! // that the host's loader probes to validate compatibility.
//! valenx_plugin_abi_version!();
//!
//! // (Future: `valenx_plugin_register!(MyAdapter::new())` to wire
//! // the adapter into the host registry.)
//! ```

#![deny(unsafe_code)]
#![allow(missing_docs)]

/// Plugin ABI version this SDK was built against. The host loader
/// (`valenx_plugin::loader::HOST_PLUGIN_ABI_VERSION`) MUST match;
/// the macro in this crate emits a function that returns this
/// value, and the host probes it on `dlopen`.
///
/// Bumped manually whenever the FFI surface changes. Cross-crate
/// sync test in `valenx-plugin::loader` keeps the two constants
/// from drifting.
pub const PLUGIN_ABI_VERSION: u32 = 1;

/// Symbol name the host's loader looks up. Convention: lowercase,
/// underscore-separated.
pub const PLUGIN_ABI_VERSION_SYMBOL: &str = "valenx_plugin_abi_version";

/// Emit the `extern "C" fn valenx_plugin_abi_version() -> u32`
/// symbol the host's loader probes. Call exactly once at the crate
/// root of the plugin's `lib.rs`.
///
/// The macro is needed (instead of plugin authors writing the
/// function by hand) to keep the symbol name + signature
/// authoritative on the SDK side — a typo in either would compile
/// but produce a plugin the host can't load. Codegen via macro
/// makes the wire shape impossible to get wrong.
#[macro_export]
macro_rules! valenx_plugin_abi_version {
    () => {
        /// Host loader entry point. Returns the ABI version the
        /// plugin was built against; a mismatch with the host's
        /// `HOST_PLUGIN_ABI_VERSION` produces a structured
        /// `LoaderError::AbiMismatch` and prevents the plugin from
        /// being registered.
        ///
        /// `#[no_mangle]` is required so the host loader can find
        /// the symbol by name. Recent Rust versions classify
        /// `#[no_mangle]` as an unsafe attribute (multiple crates
        /// could collide on the same exported symbol); we allow it
        /// here on the macro's behalf so plugin-author crates can
        /// opt into `#![forbid(unsafe_code)]` everywhere else.
        #[allow(unsafe_code)]
        #[no_mangle]
        pub extern "C" fn valenx_plugin_abi_version() -> u32 {
            $crate::PLUGIN_ABI_VERSION
        }
    };
}

// Re-exports for plugin authors. `use valenx_plugin_sdk::Adapter`
// etc.; the trait + supporting types live in valenx-core because
// the host registry uses them too.
pub use valenx_core::{
    Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plugin_abi_version_constant_is_one() {
        // Locked at 1 for the v0 ABI. Bumping this needs a
        // matching bump in valenx-plugin::loader::HOST_PLUGIN_ABI_VERSION
        // — the cross-crate sync test there will surface the drift.
        assert_eq!(PLUGIN_ABI_VERSION, 1);
    }

    #[test]
    fn plugin_abi_version_symbol_matches_host_convention() {
        // Same string the host's loader probes for. Don't tweak
        // without updating valenx-plugin::loader::ABI_VERSION_SYMBOL.
        assert_eq!(PLUGIN_ABI_VERSION_SYMBOL, "valenx_plugin_abi_version");
    }

    // Exercise the macro by emitting the symbol in this crate's test
    // binary. We can't call extern "C" functions through a normal
    // Rust import, but the compile + link of the test binary
    // confirms the macro produces a well-formed function definition
    // (no `pub` collision, no signature error).
    valenx_plugin_abi_version!();

    #[test]
    fn macro_emitted_function_returns_the_constant() {
        // The macro emits a `#[no_mangle] pub extern "C" fn
        // valenx_plugin_abi_version() -> u32` at the call site.
        // Coerce it to a function pointer + call to confirm the
        // macro's output shape compiles + behaves; calling an
        // `extern "C" fn` from Rust doesn't require `unsafe` (only
        // `unsafe extern "C" fn` would).
        let f: extern "C" fn() -> u32 = valenx_plugin_abi_version;
        assert_eq!(f(), PLUGIN_ABI_VERSION);
    }

    #[test]
    fn re_exports_compile_for_plugin_authors() {
        // Sanity check that the public Adapter-trait ergonomics are
        // available behind a single `use valenx_plugin_sdk::*;`
        // path. If any of these fail to import the crate is broken.
        let _: Option<Box<dyn Adapter>> = None;
        let _: Option<AdapterError> = None;
        let _: Option<AdapterInfo> = None;
        let _: Option<Capabilities> = None;
        let _: Option<Capability> = None;
        let _: Option<Case> = None;
        let _: Option<LicenseMode> = None;
        let _: Option<Physics> = None;
        let _: Option<PreparedJob> = None;
        let _: Option<ProbeReport> = None;
        let _: Option<RunContext> = None;
        let _: Option<RunReport> = None;
        let _: Option<VersionRange> = None;
    }
}

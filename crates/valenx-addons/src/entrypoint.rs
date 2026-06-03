//! Add-on entry-point variants.
//!
//! v1 ships Python entry-points — the desktop shell's embedded
//! valenx-py interpreter imports the add-on's `main.py` (or whichever
//! module the manifest names) at session start. The WASM variant is
//! parked for a future phase that will land an in-process
//! `wasmtime`-backed sandbox.

use serde::{Deserialize, Serialize};

/// What kind of dispatched code the add-on contributes.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum EntryPoint {
    /// Python module path inside the add-on's directory, e.g.
    /// `"main"` or `"my_pkg.api"`. The desktop shell imports this
    /// module via the Phase 11 `valenx` interpreter.
    Python {
        /// Module name (no `.py`).
        module: String,
    },
    /// WASM bytecode path (deferred — v1 doesn't dispatch).
    Wasm {
        /// Path to the `.wasm` file, relative to the add-on dir.
        wasm: String,
    },
}

impl EntryPoint {
    /// True if this entry-point can be dispatched by the v1 runtime
    /// (Python only).
    pub fn is_supported(&self) -> bool {
        matches!(self, EntryPoint::Python { .. })
    }

    /// User-facing label for the Add-on Manager UI.
    pub fn label(&self) -> &'static str {
        match self {
            EntryPoint::Python { .. } => "Python",
            EntryPoint::Wasm { .. } => "WASM (deferred)",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_is_supported() {
        let ep = EntryPoint::Python {
            module: "main".into(),
        };
        assert!(ep.is_supported());
        assert_eq!(ep.label(), "Python");
    }

    #[test]
    fn wasm_is_not_supported_in_v1() {
        let ep = EntryPoint::Wasm {
            wasm: "lib.wasm".into(),
        };
        assert!(!ep.is_supported());
        assert_eq!(ep.label(), "WASM (deferred)");
    }

    #[test]
    fn round_trips_via_toml() {
        let ep = EntryPoint::Python {
            module: "my_pkg.api".into(),
        };
        let t = toml::to_string(&ep).unwrap();
        let ep2: EntryPoint = toml::from_str(&t).unwrap();
        assert_eq!(ep, ep2);
    }
}

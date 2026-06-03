//! # valenx-fonts
//!
//! Bundled fonts: **Inter** (SIL OFL, UI body) and **JetBrains Mono**
//! (Apache 2.0, code and numerics). Shipped inside the binary so
//! rendering is identical on every platform regardless of installed
//! system fonts.
//!
//! See [DESIGN.md § 5 Typography specifics](../DESIGN.md#5-the-design-system--three-layers).
//!
//! Font files live under `assets/` and are embedded via `include_bytes!`.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

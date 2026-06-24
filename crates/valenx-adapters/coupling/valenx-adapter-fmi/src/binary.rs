//! Binary FMU loading — **off by default**, behind the `binary-fmu`
//! feature.
//!
//! The DEFAULT and tested path of this crate is the native in-process
//! [`crate::cosim::Subsystem`]. This module is the optional bridge that
//! `dlopen`s a real compiled co-simulation FMU shared library
//! (`.so` / `.dll` / `.dylib`) so it too can act as a `Subsystem`.
//!
//! **Honest status:** this path is compiled but **not exercised by this
//! crate's test suite** — running it needs an actual FMU binary on disk,
//! which is not available in CI. It therefore loads the library fail-loud
//! and exposes the FMU's interface via the parsed
//! [`crate::fmi::ModelDescription`]; wiring the full
//! `fmi2DoStep` / `fmi3DoStep` C entry-point call sequence into
//! [`crate::cosim::Subsystem::step`] is intentionally left as the next
//! increment rather than faked. Anything not yet implemented returns
//! [`crate::error::FmiError::BinaryFmu`], never a fabricated output.

use std::path::Path;

use libloading::Library;

use crate::error::{FmiError, Result};
use crate::fmi::ModelDescription;

/// A handle to a dynamically-loaded co-simulation FMU shared library plus
/// its parsed model description.
pub struct BinaryFmu {
    /// Loaded shared library (kept alive for the FMU's lifetime).
    _library: Library,
    /// Parsed `modelDescription.xml` interface.
    pub model_description: ModelDescription,
}

impl BinaryFmu {
    /// Load an FMU's shared library from `library_path` and pair it with an
    /// already-parsed [`ModelDescription`].
    ///
    /// Fail-loud: a `dlopen` / `LoadLibrary` failure becomes
    /// [`FmiError::BinaryFmu`] rather than a panic.
    ///
    /// # Safety
    ///
    /// Loading an arbitrary native library and later calling its C entry
    /// points is inherently unsafe: the caller must trust the FMU binary.
    /// `unsafe` is required here because `libloading::Library::new` runs
    /// the library's initializers.
    #[allow(unsafe_code)] // the one audited unsafe site: dlopen an FMU binary
    pub unsafe fn load(library_path: &Path, model_description: ModelDescription) -> Result<Self> {
        // SAFETY: loading a native library runs its initializers; the caller
        // has promised (via this fn being `unsafe`) that the FMU binary at
        // `library_path` is trusted.
        let library = unsafe {
            Library::new(library_path)
                .map_err(|e| FmiError::BinaryFmu(format!("dlopen {library_path:?}: {e}")))?
        };
        Ok(Self {
            _library: library,
            model_description,
        })
    }

    /// The number of input ports declared by the FMU's interface.
    pub fn n_inputs(&self) -> usize {
        self.model_description.inputs().len()
    }

    /// The number of output ports declared by the FMU's interface.
    pub fn n_outputs(&self) -> usize {
        self.model_description.outputs().len()
    }
}

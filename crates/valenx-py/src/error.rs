//! Error conversion: domain errors → Python exceptions.
//!
//! Each Valenx domain crate exposes its own typed error
//! ([`valenx_cad::CadError`], [`valenx_sketch::SketchError`],
//! [`valenx_feature_tree::FeatureError`]). PyO3 needs every fallible
//! Rust function to return [`pyo3::PyResult`], so we map each typed
//! error to a [`pyo3::exceptions::PyValueError`] preserving the
//! original message.
//!
//! Why `PyValueError` everywhere
//! ------------------------------
//!
//! All three error families are caused by bad inputs from the Python
//! caller (negative dimensions, malformed constraints, unknown feature
//! IDs). `PyValueError` is the closest standard Python exception. If
//! we ever surface I/O errors from persistence wrappers we'll add a
//! second mapping to `PyIOError`.

use pyo3::exceptions::PyValueError;
use pyo3::PyErr;

/// Convert a [`valenx_cad::CadError`] into a Python exception.
pub fn cad_err(e: valenx_cad::CadError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Convert a [`valenx_sketch::SketchError`] into a Python exception.
pub fn sketch_err(e: valenx_sketch::SketchError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Convert a [`valenx_feature_tree::FeatureError`] into a Python
/// exception.
pub fn feature_err(e: valenx_feature_tree::FeatureError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Convenience: any error implementing [`std::fmt::Display`] becomes
/// a `PyValueError`. Used by the constraint-string parser and other
/// ad-hoc validations.
pub fn value_err(msg: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(msg.to_string())
}

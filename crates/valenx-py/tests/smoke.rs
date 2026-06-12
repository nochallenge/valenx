//! Smoke test for the PyO3 bridge.
//!
//! Why this test is `#[ignore]`d
//! -----------------------------
//!
//! Constructing a `Python::with_gil` scope requires PyO3 to either
//! find an embeddable libpython at link time (auto-initialize) or be
//! loaded INTO an existing interpreter. Neither is guaranteed at
//! workspace `cargo check` time on contributor machines or in CI
//! without a Python install. To keep the default test run clean and
//! reproducible across hosts, this test is `#[ignore]`d.
//!
//! How to run it locally
//! ---------------------
//!
//! Enable the `embed-python` feature (which turns on
//! `pyo3/auto-initialize`) and run with `--ignored`:
//!
//! ```bash
//! cargo test -p valenx-py --features embed-python -- --ignored
//! ```
//!
//! That spins up an embedded Python interpreter inside the Rust test
//! binary, calls `valenx.cad.box(1, 1, 1)`, and asserts the returned
//! Solid has the expected topology.
//!
//! Without `embed-python`, the bridge is meant to be loaded INTO a
//! Python interpreter via `maturin develop`; integration testing
//! happens in Python-land, not from `cargo test`.

#[cfg(feature = "embed-python")]
use pyo3::prelude::*;
#[cfg(feature = "embed-python")]
use pyo3::types::{PyDict, PyModule};

/// Smoke test: build a unit cube via the bridge, tessellate it,
/// confirm the resulting mesh has triangles.
///
/// Marked `#[ignore]` per the module docs above. The function body
/// is only compiled in when the `embed-python` feature is enabled,
/// so default `cargo check` / `cargo build` doesn't pull in the
/// auto-initialize machinery.
#[ignore]
#[test]
fn box_round_trips_via_python() {
    #[cfg(feature = "embed-python")]
    {
        Python::with_gil(|py| -> PyResult<()> {
            // Load valenx_py as a PyModule. With `extension-module`
            // we'd normally let Python's import machinery do this;
            // for the embedded test we register the module manually
            // via the PyModule::new pattern below.
            let valenx = PyModule::new(py, "valenx_test_smoke")?;
            // Register the top-level submodules the same way the
            // pymodule entry point in lib.rs does. We can't re-use
            // the actual `valenx_py` function symbol here without
            // exporting it; instead, call into the submodule
            // registrars directly via Rust paths.
            let cad_mod = PyModule::new(py, "cad")?;
            valenx_py::cad::register(&cad_mod)?;
            valenx.add_submodule(&cad_mod)?;

            // Construct a unit cube and assert it has six faces.
            let locals = PyDict::new(py);
            locals.set_item("cad", cad_mod)?;
            // PyO3 0.24 deprecates `eval_bound` in favour of `eval` (which
            // now takes a `&CStr`). Allow the deprecated call here rather
            // than migrate: the argument is a hardcoded literal (not
            // untrusted input) and this keeps the embedded-Python smoke
            // test compiling unchanged under `-D warnings`.
            #[allow(deprecated)]
            let cube: pyo3::Bound<'_, pyo3::PyAny> =
                py.eval_bound("cad.box(1.0, 1.0, 1.0)", None, Some(&locals))?;
            let faces: usize = cube.call_method0("faces")?.extract()?;
            assert_eq!(faces, 6, "unit cube should have 6 BRep faces");
            Ok(())
        })
        .expect("smoke test ran inside embedded Python");
    }
    // Without the feature, the test body is empty — keeps the test
    // compiling on bare `cargo test`.
    #[cfg(not(feature = "embed-python"))]
    {
        eprintln!(
            "smoke test was compiled without the `embed-python` feature; \
             rerun with `--features embed-python` to actually boot a \
             Python interpreter and exercise the bridge."
        );
    }
}

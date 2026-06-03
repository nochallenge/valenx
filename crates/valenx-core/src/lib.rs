//! # valenx-core
//!
//! Shared runtime for Valenx: the canonical case / project types,
//! the adapter registry, the error taxonomy, the license-mode guard,
//! and the helpers every other workspace crate builds on.
//!
//! See:
//! - [ARCHITECTURE.md](../ARCHITECTURE.md) — layer cake, canonical
//!   types, workflow DAG
//! - [RFC 0002](../rfcs/0002-adapter-contract.md) — adapter contract
//! - [RFC 0004](../rfcs/0004-results-and-fields.md) — `Results` type

#![forbid(unsafe_code)]
// See valenx-fields/src/lib.rs for why `missing_docs` is relaxed
// during pre-alpha.
#![allow(missing_docs)]

pub mod adapter;
pub mod adapter_helpers;
pub mod error;
pub mod executor;
pub mod init_templates;
pub mod io_caps;
pub mod license;
pub mod physics;
pub mod project;
pub mod registry;
pub mod subprocess;
pub mod workflow;

pub use adapter::{
    Adapter, AdapterInfo, CancellationToken, Capabilities, Case, LogLevel, LogSink, PreparedJob,
    ProbeReport, ProgressSink, ResidualSample, RunContext, RunReport, VersionRange,
};
pub use error::{AdapterError, RunPhase, TranslateError};
pub use executor::{Executor, ExecutorError, ExecutorHandle, LocalExecutor, RunStatus};
pub use license::{assert_spawn_allowed, LicenseMode, LicenseModeViolation};
pub use physics::{Capability, Physics};
pub use project::{
    CaseDef, CaseHeader, LoadedProject, Project, ProjectHeader, ProjectLoadError, ProjectSaveError,
    ToolEntry, ToolsLock,
};
pub use registry::{AdapterEntry, AdapterRegistry, AdapterStatus, StatusCounts};
pub use subprocess::{FailurePhase, Hint, SubprocessReport};
pub use workflow::{PortType, Workflow, WorkflowEdge, WorkflowError, WorkflowNode};

/// Process-wide flag: when true, the Vina adapter always shells out
/// to the upstream `vina` binary, even if the case picked
/// `engine = "native"`. The app's Settings dialog flips this on
/// startup; the adapter reads it in `run()`. Using an atomic flag
/// rather than an env var avoids the pitfalls of `std::env::set_var`
/// (which is `unsafe` to call once threads exist on Linux).
mod force_external_vina {
    use std::sync::atomic::{AtomicBool, Ordering};
    static FLAG: AtomicBool = AtomicBool::new(false);
    /// Set the process-wide flag. Called once during app startup
    /// after settings are loaded.
    pub fn set(force: bool) {
        FLAG.store(force, Ordering::Relaxed);
    }
    /// Read the process-wide flag. Called by the Vina adapter on
    /// every run to decide whether to skip the native engine.
    pub fn get() -> bool {
        FLAG.load(Ordering::Relaxed)
    }
}
pub use force_external_vina::get as force_external_vina;
pub use force_external_vina::set as set_force_external_vina;

/// JSON descriptor of one registered adapter — what an LLM needs to
/// know to call it intelligently. Returned (in bulk) by
/// `valenx_app::ValenxApp::list_capabilities()`.
#[derive(Debug, Clone, serde::Serialize)]
pub struct AdapterDescriptor {
    /// Stable id, e.g. `"vina"`.
    pub id: &'static str,
    /// Human-friendly name.
    pub display_name: &'static str,
    /// "Subprocess" | "Native" | "Library" — see [`LicenseMode`].
    pub license_mode: String,
    /// Upstream tool license (e.g. "Apache-2.0", "GPL-3.0-only").
    pub tool_license: &'static str,
    /// Physics tags this adapter covers.
    pub physics: Vec<String>,
    /// Documentation URL.
    pub docs_url: &'static str,
    /// Homepage URL.
    pub homepage_url: &'static str,
}

impl AdapterDescriptor {
    /// Build from the existing [`AdapterInfo`].
    pub fn from_info(info: &AdapterInfo) -> Self {
        Self {
            id: info.id,
            display_name: info.display_name,
            license_mode: format!("{:?}", info.license_mode),
            tool_license: info.tool_license,
            physics: info.physics.iter().map(|p| format!("{p:?}")).collect(),
            docs_url: info.docs_url,
            homepage_url: info.homepage_url,
        }
    }
}

#[cfg(test)]
mod descriptor_tests {
    use super::*;

    #[test]
    fn descriptor_round_trips_to_json() {
        let info = AdapterInfo {
            id: "demo",
            display_name: "Demo Adapter",
            version_range: VersionRange {
                min_inclusive: semver::Version::new(1, 0, 0),
                max_exclusive: semver::Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Native,
            tool_license: "Apache-2.0",
            docs_url: "https://example.com/docs",
            homepage_url: "https://example.com",
        };
        let desc = AdapterDescriptor::from_info(&info);
        let json = serde_json::to_value(&desc).unwrap();
        assert_eq!(json["id"], "demo");
        assert_eq!(json["display_name"], "Demo Adapter");
        assert_eq!(json["license_mode"], "Native");
        assert!(json["physics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v == "Bio"));
    }
}

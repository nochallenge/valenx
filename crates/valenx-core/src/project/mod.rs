//! `.valenx` project model: manifest, tools lock, case definitions,
//! loader, and writer.
//!
//! Spec: [RFC 0001](../../rfcs/0001-project-file-format.md).

pub mod case_def;
pub mod loader;
pub mod manifest;
pub mod tools_lock;
pub mod validator;
pub mod writer;

pub use case_def::{CaseDef, CaseHeader};
pub use loader::{LoadedProject, ProjectLoadError, SUPPORTED_MAJOR};
pub use manifest::{
    CasesSection, GeometryEntry, GeometrySection, MeshEntry, Project, ProjectHeader, UiSection,
    UnitsConfig,
};
pub use tools_lock::{LockedIntegrationMode, ToolEntry, ToolsLock};
pub use validator::{validate as validate_project, ProjectWarning};
pub use writer::ProjectSaveError;

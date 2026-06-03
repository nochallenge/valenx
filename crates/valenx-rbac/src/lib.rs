//! # valenx-rbac
//!
//! Role-based access control. Second concrete chunk of
//! [RFC 0013](../../../rfcs/0013-enterprise-audit-rbac.md) (the
//! audit half landed last commit).
//!
//! Three roles, ordered weakest → strongest:
//!
//! | role     | read project | edit case.toml | run / prepare | manage adapters |
//! |----------|:------------:|:--------------:|:-------------:|:---------------:|
//! | `viewer` |     ✓        |       ✗        |       ✗       |        ✗        |
//! | `runner` |     ✓        |       ✓        |       ✓       |        ✗        |
//! | `admin`  |     ✓        |       ✓        |       ✓       |        ✓        |
//!
//! Per-user assignment lives in `<state_dir>/rbac.json`:
//!
//! ```json
//! {
//!   "users": {
//!     "alice@example.com": "admin",
//!     "bob@example.com": "runner",
//!     "guest@example.com": "viewer"
//!   },
//!   "default_role": "runner"
//! }
//! ```
//!
//! Unmatched users get `default_role`. Default-default is `runner`
//! (sensible single-user-machine default; admins override globally
//! or per-project).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// One of the three roles, ordered weakest first.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Role {
    /// Read-only — can browse projects and inspect results.
    Viewer,
    /// Default — can run / prepare / edit cases. Most users.
    #[default]
    Runner,
    /// Full access including registry / settings / plugin management.
    Admin,
}

/// One Action the user can attempt. Closed enum — every workflow
/// method gates on one of these and any new action goes through PR
/// review (which is the chance to discuss whether it should be
/// gated or not).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Action {
    /// Open / read a project.
    ProjectOpen,
    /// Save edits to project.toml.
    ProjectSave,
    /// Edit a case.toml.
    CaseModify,
    /// Run a case.
    RunCase,
    /// Prepare-only (no execute).
    PrepareCase,
    /// Cancel a running case.
    CancelRun,
    /// Re-probe adapters / install / remove tools.
    ManageAdapters,
    /// Settings / plugins / state-dir editing.
    ManageSettings,
}

impl Action {
    /// Minimum role required to perform this action.
    pub fn required_role(self) -> Role {
        match self {
            Self::ProjectOpen => Role::Viewer,
            Self::ProjectSave | Self::CaseModify => Role::Runner,
            Self::RunCase | Self::PrepareCase | Self::CancelRun => Role::Runner,
            Self::ManageAdapters | Self::ManageSettings => Role::Admin,
        }
    }

    /// Stable string id for audit-log emission. Matches the action
    /// vocabulary in RFC 0013 §"Format" so the same string flows
    /// through the audit log.
    pub fn audit_id(self) -> &'static str {
        match self {
            Self::ProjectOpen => "project.open",
            Self::ProjectSave => "project.save",
            Self::CaseModify => "case.modify",
            Self::RunCase => "run.start",
            Self::PrepareCase => "prepare.start",
            Self::CancelRun => "run.cancel",
            Self::ManageAdapters => "adapter.manage",
            Self::ManageSettings => "settings.modify",
        }
    }
}

/// Per-user role map persisted at `<state_dir>/rbac.json`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RbacConfig {
    /// User-id → Role. User ids match the audit log's `actor.id`
    /// (which today is `$USER` / `$USERNAME`; real auth replaces
    /// these with OIDC subs / SSO group claims).
    #[serde(default)]
    pub users: BTreeMap<String, Role>,
    /// Role assigned to any user-id not present in `users`.
    ///
    /// On the wire this field is `Option<Role>` so we can distinguish
    /// "operator explicitly set default_role" from "operator omitted
    /// the field". The distinction matters for project-level
    /// overrides: a project that only wants to add a per-user mapping
    /// shouldn't silently overwrite a hardened global `Viewer` with
    /// the `Option::None`-serializes-as-absent default. See
    /// [`Self::merge_with_project_override`].
    ///
    /// The Default for the struct still resolves this to `Runner` so
    /// single-user-machine installs just work without a config file.
    #[serde(default)]
    pub default_role: Option<Role>,
}

/// The "default-default" for an `RbacConfig` with `default_role:
/// None` — used by `role_for` to ensure the unconfigured-user path
/// always returns *some* concrete role. Public so test code can pin
/// the value.
pub const DEFAULT_DEFAULT_ROLE: Role = Role::Runner;

impl Default for RbacConfig {
    fn default() -> Self {
        Self {
            users: BTreeMap::new(),
            default_role: Some(DEFAULT_DEFAULT_ROLE),
        }
    }
}

impl RbacConfig {
    /// Resolve the effective default role for unconfigured users.
    pub fn effective_default_role(&self) -> Role {
        self.default_role.unwrap_or(DEFAULT_DEFAULT_ROLE)
    }

    /// Resolve the role for a given user id.
    pub fn role_for(&self, user_id: &str) -> Role {
        self.users
            .get(user_id)
            .copied()
            .unwrap_or_else(|| self.effective_default_role())
    }

    /// Quick yes/no for "can this user do this action?"
    pub fn can(&self, user_id: &str, action: Action) -> bool {
        self.role_for(user_id) >= action.required_role()
    }

    /// Merge a project-level override on top of this config. The
    /// override's `users` map is unioned over the base (project-level
    /// per-user roles win on conflict); the override's `default_role`
    /// replaces the base's ONLY when the override explicitly set the
    /// field (`Some(_)`) — an omitted field (`None`) preserves the
    /// global. Without this distinction a hardened global `Viewer`
    /// would be silently demoted to `Runner` whenever a project added
    /// even one per-user mapping (round-3 fix).
    ///
    /// Use case: a shared cluster install has a global rbac.json that
    /// lists `alice = Runner`; a sensitive project bumps `alice` to
    /// `Admin` for that project alone via its `project.toml` — and a
    /// separate locked-down project that ONLY adds a per-user
    /// mapping doesn't silently drop the global from Viewer to
    /// Runner.
    pub fn merge_with_project_override(&self, override_cfg: &Self) -> Self {
        let mut merged = self.clone();
        for (user, role) in &override_cfg.users {
            merged.users.insert(user.clone(), *role);
        }
        // Only override the default role when the project explicitly
        // wrote one. `Some(role)` means "operator picked this"; `None`
        // means "operator was silent, keep global".
        if let Some(role) = override_cfg.default_role {
            merged.default_role = Some(role);
        }
        merged
    }
}

/// Errors surfaced by the RBAC loader and authorisation checks.
#[derive(Debug, Error)]
pub enum RbacError {
    /// Filesystem IO error reading or writing the `rbac.json` config.
    #[error("rbac config at {path}: {source}")]
    Io {
        /// The offending path.
        path: std::path::PathBuf,
        /// Underlying [`std::io::Error`].
        #[source]
        source: std::io::Error,
    },
    /// `rbac.json` was present but malformed.
    #[error("rbac config at {path}: invalid JSON: {reason}")]
    Parse {
        /// The offending path.
        path: std::path::PathBuf,
        /// Human-readable explanation.
        reason: String,
    },
    /// `rbac.json` exceeds the `MAX_RBAC_FILE_BYTES` cap. Round-9
    /// hardening: a hostile actor with state-dir write access could
    /// otherwise hand us a multi-GB file and force the loader to
    /// allocate before parsing rejects it.
    #[error("rbac config at {path}: file is {size} bytes (cap: {cap})")]
    FileTooLarge {
        /// The offending path.
        path: std::path::PathBuf,
        /// Observed size in bytes.
        size: u64,
        /// Configured cap in bytes.
        cap: u64,
    },
    /// User attempted an action requiring a role they don't have.
    /// The workflow method translates this into a UI error and an
    /// audit-log "permission denied" entry.
    #[error(
        "user `{user}` (role `{user_role:?}`) cannot perform `{action:?}` (requires at least `{required:?}`)"
    )]
    PermissionDenied {
        /// Username that was rejected.
        user: String,
        /// Role the user actually holds.
        user_role: Role,
        /// Action they tried to perform.
        action: Action,
        /// Minimum role required to perform that action.
        required: Role,
    },
}

/// Maximum acceptable size of `rbac.json`. Round-9 hardening: real
/// rbac configs are tiny (a few dozen entries); a hostile or
/// accidental multi-GB file would otherwise force the loader to
/// allocate before parsing rejects it. 1 MiB is more than 4 orders
/// of magnitude above any realistic config.
pub const MAX_RBAC_FILE_BYTES: u64 = 1024 * 1024;

/// Load + parse `<state_dir>/rbac.json`. Returns `Ok(default)` when
/// the file doesn't exist (first-run / single-user-machine case);
/// returns `Err(Parse)` when the file is present but malformed
/// (operator misconfiguration — fail loud); returns
/// `Err(FileTooLarge)` when the file exceeds [`MAX_RBAC_FILE_BYTES`].
pub fn load(path: &std::path::Path) -> Result<RbacConfig, RbacError> {
    use std::io::Read;
    if !path.exists() {
        return Ok(RbacConfig::default());
    }
    // Round-9 hardening: file-size cap before allocating to read.
    let md = std::fs::metadata(path).map_err(|e| RbacError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    if md.len() > MAX_RBAC_FILE_BYTES {
        return Err(RbacError::FileTooLarge {
            path: path.to_path_buf(),
            size: md.len(),
            cap: MAX_RBAC_FILE_BYTES,
        });
    }
    let mut buf = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| RbacError::Io {
            path: path.to_path_buf(),
            source: e,
        })?
        .take(MAX_RBAC_FILE_BYTES)
        .read_to_end(&mut buf)
        .map_err(|e| RbacError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
    let text = std::str::from_utf8(&buf).map_err(|e| RbacError::Parse {
        path: path.to_path_buf(),
        reason: format!("invalid UTF-8: {e}"),
    })?;
    serde_json::from_str(text).map_err(|e| RbacError::Parse {
        path: path.to_path_buf(),
        reason: e.to_string(),
    })
}

/// Convenience guard: returns `Ok(())` if the user can perform the
/// action, `Err(PermissionDenied)` otherwise. Workflow methods
/// `?`-propagate this and the caller maps the error onto a UI
/// notification + an audit "permission denied" entry.
pub fn require(config: &RbacConfig, user: &str, action: Action) -> Result<(), RbacError> {
    if config.can(user, action) {
        Ok(())
    } else {
        let user_role = config.role_for(user);
        Err(RbacError::PermissionDenied {
            user: user.to_string(),
            user_role,
            action,
            required: action.required_role(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn role_ordering_is_weakest_first() {
        assert!(Role::Viewer < Role::Runner);
        assert!(Role::Runner < Role::Admin);
        // Default is Runner — middle ground for fresh installs.
        assert_eq!(Role::default(), Role::Runner);
    }

    #[test]
    fn action_required_role_matches_rfc_table() {
        assert_eq!(Action::ProjectOpen.required_role(), Role::Viewer);
        assert_eq!(Action::CaseModify.required_role(), Role::Runner);
        assert_eq!(Action::RunCase.required_role(), Role::Runner);
        assert_eq!(Action::ManageAdapters.required_role(), Role::Admin);
        assert_eq!(Action::ManageSettings.required_role(), Role::Admin);
    }

    #[test]
    fn action_audit_id_matches_rfc_vocabulary() {
        // Lock down the strings that flow into the audit log so
        // the schema stays stable across releases.
        assert_eq!(Action::RunCase.audit_id(), "run.start");
        assert_eq!(Action::PrepareCase.audit_id(), "prepare.start");
        assert_eq!(Action::CancelRun.audit_id(), "run.cancel");
        assert_eq!(Action::ProjectOpen.audit_id(), "project.open");
    }

    #[test]
    fn unconfigured_user_gets_default_role() {
        let config = RbacConfig::default();
        // Empty config + no env override = Runner default.
        assert_eq!(config.role_for("anyone"), Role::Runner);
        assert!(config.can("anyone", Action::RunCase));
        assert!(!config.can("anyone", Action::ManageAdapters));
    }

    #[test]
    fn explicit_user_overrides_default() {
        let mut config = RbacConfig::default();
        config.users.insert("alice@example.com".into(), Role::Admin);
        config.users.insert("bob@example.com".into(), Role::Viewer);
        config.default_role = Some(Role::Runner);

        assert!(config.can("alice@example.com", Action::ManageAdapters));
        assert!(!config.can("bob@example.com", Action::RunCase));
        assert!(config.can("anyone-else", Action::RunCase));
    }

    #[test]
    fn require_returns_permission_denied_with_full_detail() {
        let config = RbacConfig::default();
        let err = require(&config, "alice", Action::ManageAdapters).unwrap_err();
        match err {
            RbacError::PermissionDenied {
                user,
                user_role,
                action,
                required,
            } => {
                assert_eq!(user, "alice");
                assert_eq!(user_role, Role::Runner);
                assert_eq!(action, Action::ManageAdapters);
                assert_eq!(required, Role::Admin);
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn load_missing_file_returns_default_config() {
        let path = std::env::temp_dir().join(format!(
            "valenx-rbac-missing-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // File does not exist.
        let config = load(&path).expect("load");
        assert_eq!(config.effective_default_role(), Role::Runner);
        assert!(config.users.is_empty());
    }

    #[test]
    fn load_parses_well_formed_json() {
        let path = std::env::temp_dir().join(format!(
            "valenx-rbac-parse-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(
            &path,
            r#"{
              "users": {
                "alice": "admin",
                "bob": "viewer"
              },
              "default_role": "runner"
            }"#,
        )
        .unwrap();
        let config = load(&path).expect("parse");
        assert_eq!(config.users.get("alice"), Some(&Role::Admin));
        assert_eq!(config.users.get("bob"), Some(&Role::Viewer));
        assert_eq!(config.default_role, Some(Role::Runner));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn load_rejects_malformed_json() {
        let path = std::env::temp_dir().join(format!(
            "valenx-rbac-bad-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, "not json at all").unwrap();
        let err = load(&path).unwrap_err();
        assert!(matches!(err, RbacError::Parse { .. }));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn merge_with_project_override_promotes_a_user_to_admin_for_one_project() {
        let mut global = RbacConfig::default();
        global.users.insert("alice".into(), Role::Runner);
        global.users.insert("bob".into(), Role::Runner);

        let mut project_override = RbacConfig::default();
        project_override.users.insert("alice".into(), Role::Admin);

        let merged = global.merge_with_project_override(&project_override);
        // Alice is bumped to Admin within this project's context.
        assert_eq!(merged.role_for("alice"), Role::Admin);
        // Bob keeps the global Runner.
        assert_eq!(merged.role_for("bob"), Role::Runner);
        // Unknown user falls through to the merged default.
        assert_eq!(merged.role_for("eve"), Role::Runner);
    }

    #[test]
    fn merge_with_project_override_can_demote_a_user() {
        let mut global = RbacConfig::default();
        global.users.insert("alice".into(), Role::Admin);

        let mut project_override = RbacConfig::default();
        project_override.users.insert("alice".into(), Role::Viewer);

        let merged = global.merge_with_project_override(&project_override);
        assert_eq!(merged.role_for("alice"), Role::Viewer);
        // Verify the demote actually blocks RunCase (Viewer < Runner).
        assert!(!merged.can("alice", Action::RunCase));
    }

    #[test]
    fn merge_with_project_override_changes_the_default_role() {
        let global = RbacConfig::default(); // default = Runner
        let project_override = RbacConfig {
            default_role: Some(Role::Viewer), // lock down
            ..RbacConfig::default()
        };

        let merged = global.merge_with_project_override(&project_override);
        // No explicit user mapping for "anyone" => merged default = Viewer.
        assert_eq!(merged.role_for("anyone"), Role::Viewer);
        assert!(!merged.can("anyone", Action::RunCase));
    }

    /// Round-3 fix: a project-level override that omits `default_role`
    /// must preserve the global default. Previously the `default_role:
    /// Role` shape (no Option) coupled with serde-default = Runner
    /// silently demoted a hardened global `Viewer` to `Runner` the
    /// moment any project added a per-user mapping.
    #[test]
    fn project_override_without_default_role_preserves_global() {
        let global = RbacConfig {
            default_role: Some(Role::Viewer),
            ..Default::default()
        };
        // Project override that only adds a per-user mapping, leaves
        // default_role omitted on the wire — None after parse.
        let override_json = r#"{"users": {}}"#;
        let override_cfg: RbacConfig = serde_json::from_str(override_json).expect("parse");
        assert_eq!(override_cfg.default_role, None);

        let merged = global.merge_with_project_override(&override_cfg);
        assert_eq!(
            merged.effective_default_role(),
            Role::Viewer,
            "global Viewer must NOT be silently demoted to Runner"
        );
        assert!(!merged.can("anyone", Action::RunCase));
    }

    #[test]
    fn merge_with_project_override_preserves_global_users_not_mentioned() {
        let mut global = RbacConfig::default();
        global.users.insert("admin".into(), Role::Admin);
        global.users.insert("alice".into(), Role::Runner);

        let mut project_override = RbacConfig::default();
        project_override.users.insert("alice".into(), Role::Viewer);

        let merged = global.merge_with_project_override(&project_override);
        // The global admin user is still Admin (project didn't touch them).
        assert_eq!(merged.role_for("admin"), Role::Admin);
        // Alice flipped to Viewer per the project override.
        assert_eq!(merged.role_for("alice"), Role::Viewer);
    }

    /// Round-9 RED→GREEN: `load` used to call `read_to_string` with
    /// no upper bound on file size. A hostile actor with state-dir
    /// write access could otherwise hand us a multi-GB file and
    /// force the loader to allocate before parsing rejects it. The
    /// fix gates the read with [`MAX_RBAC_FILE_BYTES`] (1 MiB).
    #[test]
    fn load_rejects_files_above_max_rbac_file_bytes_cap() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rbac-big-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // 5 MiB of valid JSON-shaped padding inside a string field.
        // Even though it's syntactically valid JSON, the cap fires
        // before the parser ever sees it.
        let padding = "A".repeat(5 * 1024 * 1024);
        let body = format!(r#"{{"users": {{"alice": "viewer"}}, "padding": "{padding}"}}"#);
        std::fs::write(&tmp, body).unwrap();
        let err = load(&tmp).expect_err("must reject 5 MiB rbac.json");
        assert!(matches!(err, RbacError::FileTooLarge { .. }), "{err:?}");
        let _ = std::fs::remove_file(&tmp);
    }
}

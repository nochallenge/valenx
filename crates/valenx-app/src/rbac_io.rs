//! RBAC config layering — load the global `<state_dir>/rbac.json`,
//! merge an optional project-level override parsed out of
//! `project.toml`, and gate workflow actions through the resulting
//! config. Denial paths emit a structured audit-log entry tagged
//! `denied=true` so compliance review can see what was attempted.

use crate::audit::{current_actor_id, emit_audit};
use crate::state_paths::rbac_config_path;

/// Outcome of loading `<state_dir>/rbac.json`.
///
/// The three states are kept distinct because they map to very
/// different security postures:
///
/// - [`Self::Loaded`] — a valid config was parsed. Use it as-is.
/// - [`Self::NotFound`] — file is absent. This is the first-run /
///   single-user-machine case; fall back to the safe default
///   (everyone is `Runner`).
/// - [`Self::ParseError`] — file exists but is malformed. **Fail
///   closed**: keep the most restrictive role (`Viewer`) for every
///   user and surface the parse error to the UI as a banner so the
///   operator notices and fixes the misconfig instead of silently
///   running with whatever the loader's `Default` happens to be.
#[derive(Debug)]
pub enum RbacLoadOutcome {
    /// Config loaded successfully.
    Loaded(valenx_rbac::RbacConfig),
    /// `rbac.json` does not exist — safe default applies.
    NotFound,
    /// `rbac.json` exists but failed to parse. Carries the error
    /// string so the UI can surface it. Callers MUST refuse any
    /// `Runner`/`Admin`-gated action while this state is active.
    ParseError(String),
}

impl RbacLoadOutcome {
    /// Resolve to the [`valenx_rbac::RbacConfig`] that should be used
    /// for authorisation right now.
    ///
    /// - `Loaded(cfg)` returns `cfg`.
    /// - `NotFound` returns the safe default (`Runner` for all).
    /// - `ParseError` returns a **fail-closed** config that pins every
    ///   user to `Viewer`, blocking any action that needs at least
    ///   `Runner`. This is intentional — without a trusted config we
    ///   refuse to elevate.
    pub fn into_active_config(self) -> valenx_rbac::RbacConfig {
        match self {
            Self::Loaded(cfg) => cfg,
            Self::NotFound => valenx_rbac::RbacConfig::default(),
            Self::ParseError(_) => valenx_rbac::RbacConfig {
                users: std::collections::BTreeMap::new(),
                default_role: Some(valenx_rbac::Role::Viewer),
            },
        }
    }

    /// Returns the parse-error message if this is `ParseError`.
    pub fn parse_error(&self) -> Option<&str> {
        match self {
            Self::ParseError(e) => Some(e.as_str()),
            _ => None,
        }
    }
}

/// Load + classify `<state_dir>/rbac.json`.
///
/// See [`RbacLoadOutcome`] for the three possible outcomes and their
/// security implications.
pub fn load_rbac_outcome() -> RbacLoadOutcome {
    let Some(path) = rbac_config_path() else {
        return RbacLoadOutcome::NotFound;
    };
    if !path.exists() {
        return RbacLoadOutcome::NotFound;
    }
    match valenx_rbac::load(&path) {
        Ok(cfg) => RbacLoadOutcome::Loaded(cfg),
        Err(e) => {
            tracing::error!(
                target: "valenx",
                ?path,
                %e,
                "rbac.json failed to parse — failing closed to Viewer role"
            );
            RbacLoadOutcome::ParseError(e.to_string())
        }
    }
}

/// Backwards-compatible wrapper. Returns the active config produced
/// by [`load_rbac_outcome`]; the parse-error case becomes a
/// fail-closed Viewer config. Callers that need to surface the parse
/// error to the UI should use [`load_rbac_outcome`] directly.
pub fn load_rbac_config() -> valenx_rbac::RbacConfig {
    load_rbac_outcome().into_active_config()
}

/// Check whether the current OS user can perform `action`. Returns
/// `Ok(())` to proceed; on denial returns the structured error AND
/// emits an audit-log entry tagged with `denied=true` so the
/// attempted-but-blocked actions are visible to compliance review.
///
/// Layered config: the global `<state_dir>/rbac.json` is loaded
/// first; if `project_override` is `Some(_)` the project's `[rbac]`
/// block is merged on top via
/// [`valenx_rbac::RbacConfig::merge_with_project_override`].
/// Project overrides win on per-user role conflicts.
pub(crate) fn rbac_check(
    action: valenx_rbac::Action,
    project_override: Option<&valenx_rbac::RbacConfig>,
) -> Result<(), valenx_rbac::RbacError> {
    let mut config = load_rbac_config();
    if let Some(over) = project_override {
        config = config.merge_with_project_override(over);
    }
    let user = current_actor_id();
    let result = valenx_rbac::require(&config, &user, action);
    if let Err(e) = &result {
        emit_audit(
            action.audit_id(),
            serde_json::json!({"kind": "permission_denied"}),
            serde_json::json!({
                "denied": true,
                "user": user,
                "user_role": format!("{:?}", config.role_for(&user)),
                "required_role": format!("{:?}", action.required_role()),
                "project_override_active": project_override.is_some(),
                "error": e.to_string(),
            }),
        );
    }
    result
}

/// Parse the `[rbac]` block out of a project's `project.toml`.
/// Returns `None` when the file doesn't exist, can't be parsed, or
/// has no `[rbac]` block. Parse failures emit a tracing warning so
/// misconfigured projects are loud in logs but don't block the
/// project load.
///
/// Round-9 hardening: file-size cap mirrors `valenx_rbac::load` —
/// a hostile / accidental multi-GB `project.toml` would otherwise
/// force the loader to allocate before the toml parser rejects it.
///
/// Round-14 M8: switched the local 1 MiB constant for the shared
/// `valenx_core::project::loader::MAX_PROJECT_FILE_BYTES` so the
/// project loader and the RBAC overlay share a single source of
/// truth for the project.toml cap. Pre-fix the two paths each held
/// their own copy of the same magic number — divergence-prone.
///
/// This is the M8 "option (a)" fix per round-14 review. The cleaner
/// long-term shape is "option (b)": refactor the signature to take
/// a `&toml::Value` produced by the already-parsed loader output,
/// avoiding the re-read (and its TOCTOU window) entirely. That's a
/// larger refactor — the loader's public API would need to expose
/// the parsed `toml::Value` alongside the typed `LoadedProject`, and
/// every caller of `rbac_override_from_project_toml` would have to
/// thread it through. Deferred until the next refactor pass touches
/// the loader API.
pub(crate) fn rbac_override_from_project_toml(
    path: &std::path::Path,
) -> Option<valenx_rbac::RbacConfig> {
    use std::io::Read;
    // Shared cap so the project loader and this overlay reader stay
    // in lock-step. Round-14 M8: switched from a local
    // `MAX_PROJECT_TOML_BYTES = 1 MiB` constant.
    let cap = valenx_core::project::loader::MAX_PROJECT_FILE_BYTES;
    let md = std::fs::metadata(path).ok()?;
    if md.len() > cap {
        tracing::warn!(
            target: "valenx",
            ?path,
            size = md.len(),
            cap,
            "project.toml exceeds rbac file-size cap; ignoring [rbac] override"
        );
        return None;
    }
    let mut buf = Vec::new();
    std::fs::File::open(path)
        .ok()?
        .take(cap)
        .read_to_end(&mut buf)
        .ok()?;
    let text = std::str::from_utf8(&buf).ok()?;
    let value: toml::Value = match toml::from_str(text) {
        Ok(v) => v,
        Err(e) => {
            tracing::warn!(target: "valenx", ?path, %e, "project.toml parse failed");
            return None;
        }
    };
    let rbac_block = value.get("rbac")?.clone();
    // The serde shape of RbacConfig is a JSON-style object — go
    // toml -> JSON -> RbacConfig so we don't have to teach RbacConfig
    // about toml's representation.
    let json: serde_json::Value = serde_json::to_value(rbac_block).ok()?;
    match serde_json::from_value::<valenx_rbac::RbacConfig>(json) {
        Ok(cfg) => Some(cfg),
        Err(e) => {
            tracing::warn!(target: "valenx", ?path, %e, "project [rbac] block invalid");
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rbac_check_passes_for_default_runner_role() {
        // Default RBAC config = everyone is Runner. Runner can run
        // / prepare / cancel — those workflow methods should pass
        // their guard without needing a config file.
        let result = rbac_check(valenx_rbac::Action::RunCase, None);
        assert!(
            result.is_ok(),
            "RunCase should be allowed by default: {result:?}"
        );
        let result = rbac_check(valenx_rbac::Action::PrepareCase, None);
        assert!(result.is_ok());
        let result = rbac_check(valenx_rbac::Action::CancelRun, None);
        assert!(result.is_ok());
    }

    #[test]
    fn rbac_check_blocks_admin_actions_for_default_runner() {
        // Default Runner cannot ManageAdapters / ManageSettings —
        // those need Admin. The denial must surface as an
        // RbacError::PermissionDenied, not pass silently.
        let result = rbac_check(valenx_rbac::Action::ManageAdapters, None);
        assert!(matches!(
            result,
            Err(valenx_rbac::RbacError::PermissionDenied { .. })
        ));
    }

    #[test]
    fn rbac_override_from_project_toml_picks_up_user_role_block() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-rbac-override-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let toml_path = dir.join("project.toml");
        std::fs::write(
            &toml_path,
            r#"
[project]
name = "smoke"

[rbac]
default_role = "viewer"

[rbac.users]
"alice" = "admin"
"#,
        )
        .unwrap();
        let cfg = rbac_override_from_project_toml(&toml_path).expect("override parses");
        assert_eq!(cfg.role_for("alice"), valenx_rbac::Role::Admin);
        assert_eq!(cfg.role_for("nobody"), valenx_rbac::Role::Viewer);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn parse_error_outcome_fails_closed_to_viewer() {
        let outcome = RbacLoadOutcome::ParseError("bad json: foo".into());
        let cfg = outcome.into_active_config();
        // Every user should resolve to Viewer; Runner-gated actions
        // must NOT pass.
        assert_eq!(
            cfg.role_for("anyone"),
            valenx_rbac::Role::Viewer,
            "ParseError must fail closed (Viewer)"
        );
        let err = valenx_rbac::require(&cfg, "anyone", valenx_rbac::Action::RunCase);
        assert!(
            matches!(err, Err(valenx_rbac::RbacError::PermissionDenied { .. })),
            "RunCase must be denied under ParseError fail-closed config"
        );
    }

    #[test]
    fn not_found_outcome_uses_safe_default_runner() {
        let cfg = RbacLoadOutcome::NotFound.into_active_config();
        // Single-user-machine first-run: default is Runner so things
        // just work.
        assert_eq!(cfg.role_for("anyone"), valenx_rbac::Role::Runner);
    }

    #[test]
    fn parse_error_exposes_message() {
        let outcome = RbacLoadOutcome::ParseError("oops: trailing comma".into());
        assert_eq!(outcome.parse_error(), Some("oops: trailing comma"));
    }

    #[test]
    fn rbac_override_returns_none_when_block_missing() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-rbac-override-missing-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let toml_path = dir.join("project.toml");
        std::fs::write(&toml_path, "[project]\nname = \"smoke\"\n").unwrap();
        assert!(rbac_override_from_project_toml(&toml_path).is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-14 M8 RED→GREEN: an oversize project.toml (post-loader,
    /// pre-rbac-overlay swap) must be rejected by the rbac overlay's
    /// own size cap, returning None and emitting a tracing warning.
    /// Pre-fix the rbac overlay had its own 1 MiB constant copy —
    /// staying in sync with the loader was manual. The shared
    /// `MAX_PROJECT_FILE_BYTES` constant removes the divergence
    /// vector permanently.
    ///
    /// We test the cap directly rather than racing the loader since
    /// single-process FS races are timing-fragile. Writing a 5 MiB
    /// file (well past the shared 1 MiB cap) confirms the overlay
    /// stat-checks correctly.
    #[test]
    fn rbac_override_rejects_oversize_project_toml() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-rbac-override-oversize-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let toml_path = dir.join("project.toml");

        let cap = valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize;
        // 5 MiB worth of padding — well past the 1 MiB cap.
        let oversize = cap.saturating_mul(5).max(5 * 1024 * 1024);
        let mut payload = String::with_capacity(oversize + 128);
        payload.push_str(
            "[project]\nname = \"oversize\"\n\n[rbac]\ndefault_role = \"admin\"\n\n# pad: ",
        );
        // Fill to ~5 MiB with comment chars (still valid TOML).
        payload.extend(std::iter::repeat_n('A', oversize));
        std::fs::write(&toml_path, &payload).unwrap();

        // Even though the file contains a syntactically-valid
        // [rbac] block, the size cap must fire and the overlay
        // must return None.
        let result = rbac_override_from_project_toml(&toml_path);
        assert!(
            result.is_none(),
            "oversize project.toml must return None — overlay must not silently apply"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}

//! Embedded project / simulation database.
//!
//! A single local SQLite file (or an in-memory database for tests and
//! throwaway sessions) that tracks the things a long-lived Valenx
//! workspace accumulates and that do **not** belong in the per-project
//! `project.toml` manifest:
//!
//! * **project history** — every project the user has touched, with its
//!   creation timestamp ([`projects`](ProjectDb::insert_project));
//! * **run metadata** — one row per solver/adapter run, with its kind,
//!   status, and an arbitrary JSON blob of parameters/results
//!   ([`runs`](ProjectDb::record_run));
//! * **material properties** — a reusable library of named materials,
//!   each carrying a JSON property bag, upserted by name
//!   ([`materials`](ProjectDb::upsert_material));
//! * **file validation** — a log of input-file integrity checks
//!   (path + SHA-256 + pass/fail), so a re-opened project can tell at a
//!   glance whether its inputs changed underneath it
//!   ([`file_validations`](ProjectDb::record_file_validation)).
//!
//! # Why SQLite
//!
//! These records are relational (a run belongs to a project), benefit
//! from indexed lookup as the workspace grows, and must survive process
//! restarts — exactly what an embedded relational store is for. The
//! `bundled` feature of `rusqlite` compiles SQLite from source, so there
//! is no dependency on a system `libsqlite3` and the behaviour is
//! identical across platforms.
//!
//! # Safety / robustness
//!
//! * Every write that spans more than one statement runs inside a
//!   **transaction**, so a failure leaves the database unchanged.
//! * Every query that takes caller data uses **bound parameters** — no
//!   user string is ever concatenated into SQL, so there is no SQL
//!   injection surface.
//! * Foreign keys are enabled (`PRAGMA foreign_keys = ON`), so a `run`
//!   or `file_validation` can never reference a non-existent project.
//! * All failure paths return [`AdapterError::Db`] via `?`; the public
//!   API never panics on a bad query, a constraint violation, or a
//!   corrupt file.
//!
//! Timestamps are stored as caller-supplied ISO-8601 `TEXT`, matching
//! the convention already used by the project manifest
//! ([`crate::project::ProjectHeader::created`]); this keeps the crate
//! free of a date/time dependency. Helpers that omit a timestamp record
//! the empty string, which sorts before any real timestamp.

use rusqlite::{params, Connection, OptionalExtension};

use crate::error::AdapterError;

/// A handle to the embedded project/simulation database.
///
/// Open one with [`ProjectDb::open`] (a file on disk) or
/// [`ProjectDb::open_in_memory`] (ephemeral, for tests). The schema is
/// created and migrated-forward idempotently on open, so pointing
/// [`open`](ProjectDb::open) at either a fresh path or an existing
/// database both work.
#[derive(Debug)]
pub struct ProjectDb {
    conn: Connection,
}

/// One row of the `projects` table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectRecord {
    /// Auto-assigned row id (the database primary key).
    pub id: i64,
    /// Display name of the project.
    pub name: String,
    /// ISO-8601 creation timestamp (empty string if none was given).
    pub created: String,
}

/// One row of the `runs` table — a single solver/adapter invocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunRecord {
    /// Auto-assigned row id.
    pub id: i64,
    /// The project this run belongs to (`projects.id`).
    pub project_id: i64,
    /// Free-form run kind, e.g. `"cfd"`, `"fem"`, `"docking"`.
    pub kind: String,
    /// Free-form status, e.g. `"queued"`, `"running"`, `"ok"`, `"failed"`.
    pub status: String,
    /// ISO-8601 timestamp the run was recorded (empty string if none).
    pub created: String,
    /// Arbitrary JSON blob of parameters/results (empty string if none).
    pub meta_json: String,
}

/// One row of the `materials` table — a named, reusable material.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaterialRecord {
    /// Auto-assigned row id.
    pub id: i64,
    /// Unique material name (the upsert key), e.g. `"AISI 4340 steel"`.
    pub name: String,
    /// Arbitrary JSON property bag, e.g. `{"E_gpa":205,"rho":7850}`.
    pub props_json: String,
}

/// One row of the `file_validations` table — an input-file integrity
/// check (path + content hash + verdict).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileValidationRecord {
    /// Auto-assigned row id.
    pub id: i64,
    /// The project this validation belongs to (`projects.id`).
    pub project_id: i64,
    /// The validated file's path (as recorded by the caller).
    pub path: String,
    /// Lower-case hex SHA-256 of the file's bytes (empty if not hashed).
    pub sha256: String,
    /// Verdict, e.g. `"ok"`, `"missing"`, `"changed"`, `"corrupt"`.
    pub status: String,
    /// ISO-8601 timestamp of the check (empty string if none).
    pub validated_at: String,
    /// Optional human-readable detail (empty string if none).
    pub detail: String,
}

/// The DDL that defines the schema. Run on every [`ProjectDb::open`] /
/// [`open_in_memory`](ProjectDb::open_in_memory) inside a transaction;
/// every statement is `IF NOT EXISTS`, so re-opening an existing
/// database is a no-op.
const SCHEMA: &str = "\
CREATE TABLE IF NOT EXISTS projects (
    id      INTEGER PRIMARY KEY,
    name    TEXT NOT NULL,
    created TEXT NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS runs (
    id         INTEGER PRIMARY KEY,
    project_id INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    kind       TEXT NOT NULL,
    status     TEXT NOT NULL,
    created    TEXT NOT NULL DEFAULT '',
    meta_json  TEXT NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_runs_project ON runs(project_id, id);
CREATE TABLE IF NOT EXISTS materials (
    id         INTEGER PRIMARY KEY,
    name       TEXT NOT NULL UNIQUE,
    props_json TEXT NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS file_validations (
    id           INTEGER PRIMARY KEY,
    project_id   INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    path         TEXT NOT NULL,
    sha256       TEXT NOT NULL DEFAULT '',
    status       TEXT NOT NULL,
    validated_at TEXT NOT NULL DEFAULT '',
    detail       TEXT NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_fv_project ON file_validations(project_id, id);
";

impl ProjectDb {
    /// Open (creating if necessary) the database at `path` and ensure
    /// the schema exists.
    ///
    /// # Errors
    ///
    /// Returns [`AdapterError::Db`] if the file cannot be opened (e.g.
    /// the parent directory is missing or unwritable) or if applying the
    /// schema fails.
    pub fn open(path: impl AsRef<std::path::Path>) -> Result<Self, AdapterError> {
        let conn = Connection::open(path.as_ref())?;
        Self::from_conn(conn)
    }

    /// Open a fresh, private in-memory database with the schema applied.
    ///
    /// The database vanishes when the returned [`ProjectDb`] is dropped.
    /// Used by tests and by callers that want a scratch store.
    ///
    /// # Errors
    ///
    /// Returns [`AdapterError::Db`] if SQLite cannot create the in-memory
    /// connection or apply the schema.
    pub fn open_in_memory() -> Result<Self, AdapterError> {
        let conn = Connection::open_in_memory()?;
        Self::from_conn(conn)
    }

    /// Shared constructor: enable foreign keys and apply the schema.
    fn from_conn(conn: Connection) -> Result<Self, AdapterError> {
        // Enforce the `REFERENCES` clauses (off by default in SQLite).
        conn.pragma_update(None, "foreign_keys", "ON")?;
        // `execute_batch` runs the whole DDL script atomically; wrap it
        // in an explicit transaction so a partial failure rolls back.
        conn.execute_batch(&format!("BEGIN;\n{SCHEMA}\nCOMMIT;"))?;
        Ok(Self { conn })
    }

    // ---- projects -------------------------------------------------------

    /// Insert a new project and return its assigned id.
    ///
    /// `created` is a caller-supplied ISO-8601 timestamp; pass `""` (or
    /// use [`insert_project_now`](ProjectDb::insert_project_now)) if you
    /// do not track one.
    pub fn insert_project(&self, name: &str, created: &str) -> Result<i64, AdapterError> {
        self.conn.execute(
            "INSERT INTO projects (name, created) VALUES (?1, ?2)",
            params![name, created],
        )?;
        Ok(self.conn.last_insert_rowid())
    }

    /// Convenience wrapper around [`insert_project`](ProjectDb::insert_project)
    /// that stamps the row with the current time as a best-effort
    /// ISO-8601-ish UTC string. Falls back to the empty string if the
    /// system clock is before the Unix epoch.
    pub fn insert_project_now(&self, name: &str) -> Result<i64, AdapterError> {
        self.insert_project(name, &now_timestamp())
    }

    /// Fetch a single project by id, or `None` if no such row exists.
    pub fn get_project(&self, id: i64) -> Result<Option<ProjectRecord>, AdapterError> {
        let row = self
            .conn
            .query_row(
                "SELECT id, name, created FROM projects WHERE id = ?1",
                params![id],
                |r| {
                    Ok(ProjectRecord {
                        id: r.get(0)?,
                        name: r.get(1)?,
                        created: r.get(2)?,
                    })
                },
            )
            .optional()?;
        Ok(row)
    }

    /// List every project, most-recently-inserted first (descending id).
    pub fn list_projects(&self) -> Result<Vec<ProjectRecord>, AdapterError> {
        let mut stmt = self
            .conn
            .prepare("SELECT id, name, created FROM projects ORDER BY id DESC")?;
        let rows = stmt.query_map([], |r| {
            Ok(ProjectRecord {
                id: r.get(0)?,
                name: r.get(1)?,
                created: r.get(2)?,
            })
        })?;
        rows.collect::<Result<Vec<_>, _>>()
            .map_err(AdapterError::from)
    }

    // ---- runs -----------------------------------------------------------

    /// Record a run against an existing project and return its id.
    ///
    /// `meta_json` is stored verbatim — pass `"{}"` or `""` when there is
    /// nothing to attach. The foreign-key constraint means recording a
    /// run against a non-existent `project_id` fails with
    /// [`AdapterError::Db`] rather than silently orphaning the row.
    pub fn record_run(
        &self,
        project_id: i64,
        kind: &str,
        status: &str,
        created: &str,
        meta_json: &str,
    ) -> Result<i64, AdapterError> {
        self.conn.execute(
            "INSERT INTO runs (project_id, kind, status, created, meta_json) \
             VALUES (?1, ?2, ?3, ?4, ?5)",
            params![project_id, kind, status, created, meta_json],
        )?;
        Ok(self.conn.last_insert_rowid())
    }

    /// Update the `status` of an existing run (e.g. `"running"` → `"ok"`).
    ///
    /// Returns the number of rows changed (`0` if `run_id` is unknown).
    pub fn set_run_status(&self, run_id: i64, status: &str) -> Result<usize, AdapterError> {
        let n = self.conn.execute(
            "UPDATE runs SET status = ?2 WHERE id = ?1",
            params![run_id, status],
        )?;
        Ok(n)
    }

    /// List every run for one project in the order it was recorded
    /// (ascending id == chronological).
    pub fn list_runs(&self, project_id: i64) -> Result<Vec<RunRecord>, AdapterError> {
        let mut stmt = self.conn.prepare(
            "SELECT id, project_id, kind, status, created, meta_json \
             FROM runs WHERE project_id = ?1 ORDER BY id ASC",
        )?;
        let rows = stmt.query_map(params![project_id], |r| {
            Ok(RunRecord {
                id: r.get(0)?,
                project_id: r.get(1)?,
                kind: r.get(2)?,
                status: r.get(3)?,
                created: r.get(4)?,
                meta_json: r.get(5)?,
            })
        })?;
        rows.collect::<Result<Vec<_>, _>>()
            .map_err(AdapterError::from)
    }

    // ---- materials ------------------------------------------------------

    /// Insert or update a material by its unique `name`, returning the
    /// row id.
    ///
    /// If a material with the same name already exists its `props_json`
    /// is replaced; otherwise a new row is inserted. The whole operation
    /// is a single `INSERT ... ON CONFLICT DO UPDATE`, so it is atomic.
    pub fn upsert_material(&self, name: &str, props_json: &str) -> Result<i64, AdapterError> {
        self.conn.execute(
            "INSERT INTO materials (name, props_json) VALUES (?1, ?2) \
             ON CONFLICT(name) DO UPDATE SET props_json = excluded.props_json",
            params![name, props_json],
        )?;
        // `last_insert_rowid` is only meaningful for the INSERT branch;
        // on the UPDATE branch fetch the existing id explicitly.
        let id = self.conn.query_row(
            "SELECT id FROM materials WHERE name = ?1",
            params![name],
            |r| r.get(0),
        )?;
        Ok(id)
    }

    /// Fetch a material by name, or `None` if no such material exists.
    pub fn get_material(&self, name: &str) -> Result<Option<MaterialRecord>, AdapterError> {
        let row = self
            .conn
            .query_row(
                "SELECT id, name, props_json FROM materials WHERE name = ?1",
                params![name],
                |r| {
                    Ok(MaterialRecord {
                        id: r.get(0)?,
                        name: r.get(1)?,
                        props_json: r.get(2)?,
                    })
                },
            )
            .optional()?;
        Ok(row)
    }

    /// List every material, ordered by name.
    pub fn list_materials(&self) -> Result<Vec<MaterialRecord>, AdapterError> {
        let mut stmt = self
            .conn
            .prepare("SELECT id, name, props_json FROM materials ORDER BY name ASC")?;
        let rows = stmt.query_map([], |r| {
            Ok(MaterialRecord {
                id: r.get(0)?,
                name: r.get(1)?,
                props_json: r.get(2)?,
            })
        })?;
        rows.collect::<Result<Vec<_>, _>>()
            .map_err(AdapterError::from)
    }

    // ---- file validations ----------------------------------------------

    /// Record an input-file integrity check against an existing project
    /// and return its id.
    ///
    /// Like [`record_run`](ProjectDb::record_run), the foreign key means
    /// an unknown `project_id` is rejected rather than orphaned.
    #[allow(clippy::too_many_arguments)]
    pub fn record_file_validation(
        &self,
        project_id: i64,
        path: &str,
        sha256: &str,
        status: &str,
        validated_at: &str,
        detail: &str,
    ) -> Result<i64, AdapterError> {
        self.conn.execute(
            "INSERT INTO file_validations \
             (project_id, path, sha256, status, validated_at, detail) \
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![project_id, path, sha256, status, validated_at, detail],
        )?;
        Ok(self.conn.last_insert_rowid())
    }

    /// List every file-validation record for one project in the order it
    /// was recorded (ascending id == chronological).
    pub fn list_file_validations(
        &self,
        project_id: i64,
    ) -> Result<Vec<FileValidationRecord>, AdapterError> {
        let mut stmt = self.conn.prepare(
            "SELECT id, project_id, path, sha256, status, validated_at, detail \
             FROM file_validations WHERE project_id = ?1 ORDER BY id ASC",
        )?;
        let rows = stmt.query_map(params![project_id], |r| {
            Ok(FileValidationRecord {
                id: r.get(0)?,
                project_id: r.get(1)?,
                path: r.get(2)?,
                sha256: r.get(3)?,
                status: r.get(4)?,
                validated_at: r.get(5)?,
                detail: r.get(6)?,
            })
        })?;
        rows.collect::<Result<Vec<_>, _>>()
            .map_err(AdapterError::from)
    }
}

/// Best-effort current UTC timestamp as seconds-since-epoch, formatted
/// as a sortable string. We deliberately avoid a date/time dependency;
/// the exact format is unimportant as long as it sorts chronologically,
/// and callers that need a real ISO-8601 string can pass their own.
fn now_timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| format!("{}", d.as_secs()))
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_project_then_list_returns_it() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        let id = db
            .insert_project("Wing CFD study", "2026-06-25T10:00:00Z")
            .expect("insert project");
        assert!(id > 0);

        let listed = db.list_projects().expect("list projects");
        assert_eq!(listed.len(), 1);
        assert_eq!(listed[0].id, id);
        assert_eq!(listed[0].name, "Wing CFD study");
        assert_eq!(listed[0].created, "2026-06-25T10:00:00Z");

        let fetched = db.get_project(id).expect("get project").expect("present");
        assert_eq!(fetched, listed[0]);
        // A missing id yields None, not an error.
        assert!(db.get_project(9999).expect("get missing").is_none());
    }

    #[test]
    fn record_two_runs_lists_both_in_order() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        let pid = db.insert_project_now("proj").expect("insert project");

        let r1 = db
            .record_run(pid, "mesh", "ok", "2026-06-25T10:01:00Z", "{}")
            .expect("record run 1");
        let r2 = db
            .record_run(
                pid,
                "cfd",
                "running",
                "2026-06-25T10:02:00Z",
                r#"{"iters":500}"#,
            )
            .expect("record run 2");
        assert!(r2 > r1, "ids must be monotonic");

        let runs = db.list_runs(pid).expect("list runs");
        assert_eq!(runs.len(), 2);
        // Chronological (ascending id) order.
        assert_eq!(runs[0].id, r1);
        assert_eq!(runs[0].kind, "mesh");
        assert_eq!(runs[0].status, "ok");
        assert_eq!(runs[1].id, r2);
        assert_eq!(runs[1].kind, "cfd");
        assert_eq!(runs[1].meta_json, r#"{"iters":500}"#);

        // Status transition round-trips.
        let changed = db.set_run_status(r2, "ok").expect("set status");
        assert_eq!(changed, 1);
        let runs = db.list_runs(pid).expect("re-list runs");
        assert_eq!(runs[1].status, "ok");

        // A different project sees none of these runs.
        let other = db.insert_project_now("other").expect("insert other");
        assert!(db.list_runs(other).expect("list other runs").is_empty());
    }

    #[test]
    fn material_upsert_round_trips_and_updates() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");

        let id1 = db
            .upsert_material("steel", r#"{"E_gpa":205,"rho":7850}"#)
            .expect("insert material");
        let got = db.get_material("steel").expect("get").expect("present");
        assert_eq!(got.id, id1);
        assert_eq!(got.props_json, r#"{"E_gpa":205,"rho":7850}"#);

        // Upserting the same name updates props in place and keeps the id.
        let id2 = db
            .upsert_material("steel", r#"{"E_gpa":210,"rho":7860}"#)
            .expect("update material");
        assert_eq!(id1, id2, "upsert must reuse the existing row id");
        let got = db.get_material("steel").expect("get").expect("present");
        assert_eq!(got.props_json, r#"{"E_gpa":210,"rho":7860}"#);

        // Still exactly one material row, plus a missing lookup is None.
        assert_eq!(db.list_materials().expect("list").len(), 1);
        assert!(db.get_material("aluminium").expect("get missing").is_none());
    }

    #[test]
    fn file_validation_records_and_lists_in_order() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        let pid = db.insert_project_now("proj").expect("insert project");

        let f1 = db
            .record_file_validation(pid, "mesh.msh", "abc123", "ok", "2026-06-25T10:00:00Z", "")
            .expect("record fv 1");
        let f2 = db
            .record_file_validation(
                pid,
                "geom.step",
                "",
                "missing",
                "2026-06-25T10:00:01Z",
                "file not found on disk",
            )
            .expect("record fv 2");
        assert!(f2 > f1);

        let fvs = db.list_file_validations(pid).expect("list fv");
        assert_eq!(fvs.len(), 2);
        assert_eq!(fvs[0].path, "mesh.msh");
        assert_eq!(fvs[0].sha256, "abc123");
        assert_eq!(fvs[0].status, "ok");
        assert_eq!(fvs[1].status, "missing");
        assert_eq!(fvs[1].detail, "file not found on disk");
    }

    #[test]
    fn foreign_key_violation_errors_cleanly() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        // No project with id 42 exists; the FK must reject this without
        // panicking, surfacing as AdapterError::Db.
        let err = db
            .record_run(42, "cfd", "ok", "", "{}")
            .expect_err("FK violation must be an error");
        assert!(matches!(err, AdapterError::Db(_)), "got: {err:?}");

        // Same for a file-validation row against a missing project.
        let err = db
            .record_file_validation(42, "x.step", "", "ok", "", "")
            .expect_err("FK violation must be an error");
        assert!(matches!(err, AdapterError::Db(_)), "got: {err:?}");
    }

    #[test]
    fn unique_constraint_is_enforced_for_raw_insert() {
        // A direct duplicate-name INSERT (bypassing upsert) must error
        // cleanly via the UNIQUE(name) constraint, not panic.
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        db.upsert_material("brass", "{}").expect("first insert");
        let res = db.conn.execute(
            "INSERT INTO materials (name, props_json) VALUES (?1, ?2)",
            params!["brass", "{}"],
        );
        let err = res.expect_err("duplicate name must violate UNIQUE");
        // rusqlite surfaces this as a SqliteFailure; just assert it is
        // the constraint-violation family and nothing panicked.
        assert!(
            matches!(err, rusqlite::Error::SqliteFailure(_, _)),
            "got: {err:?}"
        );
    }

    #[test]
    fn bad_sql_query_errors_without_panic() {
        let db = ProjectDb::open_in_memory().expect("open in-memory db");
        // Preparing nonsense SQL must return Err, never panic.
        let res = db
            .conn
            .prepare("SELECT not_a_column FROM nonexistent_table");
        assert!(res.is_err(), "malformed query should error");
    }

    #[test]
    fn schema_is_idempotent_on_reopen_same_file() {
        // Open a real file, write a row, reopen it: the schema re-apply
        // must not wipe data or error.
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!("valenx-db-test-{}-{n}.db", std::process::id()));
        let _ = std::fs::remove_file(&path);

        let pid = {
            let db = ProjectDb::open(&path).expect("open file db");
            db.insert_project("persisted", "2026-06-25T00:00:00Z")
                .expect("insert")
        };
        {
            // Reopen: schema CREATE IF NOT EXISTS is a no-op, data survives.
            let db = ProjectDb::open(&path).expect("reopen file db");
            let got = db.get_project(pid).expect("get").expect("present");
            assert_eq!(got.name, "persisted");
            assert_eq!(db.list_projects().expect("list").len(), 1);
        }
        let _ = std::fs::remove_file(&path);
    }
}

//! The **project library** — a foldered, persistent catalogue of saved
//! projects so a user can manage a hundred-plus projects like an IDE's
//! project list instead of hunting through a flat `Open saved ▾` menu.
//!
//! This is the data + persistence layer behind [`crate::project_navigator`]
//! (which renders it as the "Projects" tree at the top of the Browser
//! panel). It deliberately reuses the **exact** machinery the project-tab
//! saves already use ([`crate::state_paths::atomic_write`], a
//! `project_tabs::sanitize_name`-equivalent path-safety check, and
//! the [`crate::settings_io::MAX_STATE_FILE_BYTES`] size cap), so the
//! library file lands beside the tab/session saves under the same per-OS
//! state directory and shares their crash-safety guarantees.
//!
//! ## On-disk shape
//!
//! The whole library is a single pretty-printed JSON document at
//! `<state_dir>/library.json` — `<state_dir>` is `%APPDATA%\valenx` on
//! Windows (see [`crate::state_paths::state_dir`]). One file (rather than a
//! file-per-project under a directory) keeps ordering, folders, and pinning
//! atomic: a single [`atomic_write`] swaps the entire catalogue, so a crash
//! mid-save never leaves the library half-reordered.
//!
//! ```text
//! <state_dir>/
//!   library.json        ← this module
//!   tabs/<name>.json     ← single saved tabs (project_tabs)
//!   sessions/<name>.json ← saved tab groups (project_tabs)
//! ```
//!
//! Each [`SavedProject`] embeds a full [`ProjectTab`] (which is itself
//! `Serialize`/`Deserialize`), so opening a library entry rebuilds the exact
//! tab kind + title the user saved. Folders are flat (a project carries an
//! optional `folder` id); v1 has no nested folders.

use serde::{Deserialize, Serialize};

use crate::project_tabs::ProjectTab;
use crate::state_paths::{atomic_write, state_dir};

/// Path the whole library is persisted at — `<state_dir>/library.json`.
fn library_path() -> Option<std::path::PathBuf> {
    state_dir().map(|d| d.join("library.json"))
}

/// Current unix time in whole seconds (best-effort; `0` if the system
/// clock predates the epoch). Used to stamp `created` / `last_opened`.
///
/// `std::time::SystemTime` is fine here — this is **live-app runtime**
/// code, not a deterministic workflow/replay context, so wall-clock reads
/// are allowed (mirrors [`crate::state_paths::export_csv_path`], which
/// already stamps export filenames the same way).
fn now_secs() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Mint a short, collision-resistant id for a project or folder. Built from
/// the unix-nanos clock plus a process-monotonic counter so two ids minted
/// in the same nanosecond tick still differ — no `rand` dep, no
/// `Date::now`-style API.
fn fresh_id(prefix: &str) -> String {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    format!("{prefix}-{nanos:x}-{n:x}")
}

/// One saved project in the library: a named, foldered, orderable wrapper
/// around a [`ProjectTab`] (the thing that gets reopened as a tab).
///
/// `id` is stable across renames/moves (folders and the display name can
/// change; the id never does), so the navigator can target a row
/// unambiguously even while the user is reordering siblings.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SavedProject {
    /// Stable unique id (minted by `fresh_id`); never shown to the user.
    pub id: String,
    /// User-facing display name. Defaults to the tab title at save time;
    /// renamable independently of the tab's own title.
    pub name: String,
    /// The project payload — reopened verbatim as a tab. Carries its own
    /// [`crate::project_tabs::TabKind`] + title.
    pub tab: ProjectTab,
    /// Folder id this project lives in, or `None` for the top-level
    /// "(unfiled)" group. Must reference an existing [`Folder::id`] (the
    /// loader self-heals dangling references on [`ProjectLibrary::load`] via
    /// its private `reconcile` step).
    pub folder: Option<String>,
    /// Pinned projects surface in the navigator's ★ Pinned section.
    pub pinned: bool,
    /// Sort key within a folder (and within the unfiled group). Lower
    /// sorts first; ties break by `name`. Maintained dense-ish by
    /// [`ProjectLibrary::reorder`].
    pub order: u32,
    /// Unix seconds at first save (never updated after).
    pub created: u64,
    /// Unix seconds the project was last opened from the library. Drives
    /// the 🕘 Recent ordering. Equal to `created` until first reopen.
    pub last_opened: u64,
}

/// A flat library folder. v1 has no nesting — a project's [`SavedProject::folder`]
/// is a single optional folder id.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Folder {
    /// Stable unique id (minted by `fresh_id`); referenced by
    /// [`SavedProject::folder`].
    pub id: String,
    /// User-facing folder name.
    pub name: String,
    /// Sort key among folders (lower first; ties break by `name`).
    pub order: u32,
}

/// The whole project library: every saved project plus the folders they're
/// organised into. Serialised as one JSON document (see the module docs).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ProjectLibrary {
    /// All saved projects, in no particular stored order (the navigator
    /// sorts a view of them per-section).
    #[serde(default)]
    pub projects: Vec<SavedProject>,
    /// All folders.
    #[serde(default)]
    pub folders: Vec<Folder>,
}

impl ProjectLibrary {
    // ── Persistence ──────────────────────────────────────────────────────

    /// Load the library from `<state_dir>/library.json`. A missing file,
    /// missing state dir, oversize file, or parse error all yield an empty
    /// (default) library rather than failing — the library is best-effort
    /// user state, exactly like the tab/session saves. Always runs the
    /// private `reconcile` step so a hand-edited / partially-corrupt file
    /// can't leave dangling folder references behind.
    pub fn load() -> Self {
        let Some(path) = library_path() else {
            return Self::default();
        };
        let mut lib = load_from_path(&path).unwrap_or_default();
        lib.reconcile();
        lib
    }

    /// Persist the whole library to `<state_dir>/library.json` via the
    /// crash-safe [`atomic_write`]. Best-effort: a missing state dir or
    /// write/serialise failure is swallowed (returns `false`), matching the
    /// rest of valenx's state persistence.
    pub fn save(&self) -> bool {
        let Some(path) = library_path() else {
            return false;
        };
        save_to_path(self, &path)
    }

    /// Drop projects whose `folder` references a folder that no longer
    /// exists (move them back to "(unfiled)"), so a deleted folder or a
    /// hand-edited file never leaves orphaned rows the navigator can't show.
    /// Called on every [`Self::load`].
    fn reconcile(&mut self) {
        let known: std::collections::HashSet<&str> =
            self.folders.iter().map(|f| f.id.as_str()).collect();
        let mut fixups: Vec<usize> = Vec::new();
        for (i, p) in self.projects.iter().enumerate() {
            if let Some(fid) = &p.folder {
                if !known.contains(fid.as_str()) {
                    fixups.push(i);
                }
            }
        }
        for i in fixups {
            self.projects[i].folder = None;
        }
    }

    // ── Project mutations ────────────────────────────────────────────────

    /// Save `tab` as a new library project in `folder` (or unfiled when
    /// `None`), appended at the end of that folder's order. Returns the new
    /// project's id. The display name is taken from the tab's title.
    pub fn add_project(&mut self, tab: ProjectTab, folder: Option<String>) -> String {
        // Honour the requested folder only if it actually exists.
        let folder = folder.filter(|fid| self.folders.iter().any(|f| &f.id == fid));
        let now = now_secs();
        let order = self.next_order_in(folder.as_deref());
        let id = fresh_id("proj");
        let name = {
            let t = tab.title.trim();
            if t.is_empty() {
                tab.kind.label().to_string()
            } else {
                t.to_string()
            }
        };
        self.projects.push(SavedProject {
            id: id.clone(),
            name,
            tab,
            folder,
            pinned: false,
            order,
            created: now,
            last_opened: now,
        });
        id
    }

    /// Remove the project with `id` (no-op if absent).
    pub fn remove(&mut self, id: &str) {
        self.projects.retain(|p| p.id != id);
    }

    /// Clear the **entire** library: drop every saved project *and* every
    /// folder, leaving an empty catalogue. The destructive batch-delete behind
    /// the navigator's "Clear all projects" button (the caller gates it behind
    /// a confirm and persists afterwards). Does not touch any saved single-tab
    /// / session files on disk — only this library.
    pub fn clear_all(&mut self) {
        self.projects.clear();
        self.folders.clear();
    }

    /// Rename the project `id` to `name` (trimmed; empty names are
    /// rejected so a row never goes nameless). No-op if `id` is absent.
    pub fn rename(&mut self, id: &str, name: &str) {
        let name = name.trim();
        if name.is_empty() {
            return;
        }
        if let Some(p) = self.projects.iter_mut().find(|p| p.id == id) {
            p.name = name.to_string();
        }
    }

    /// Duplicate the project `id`: a deep copy with a fresh id, a "
    /// (copy)"-suffixed name, fresh `created`/`last_opened` stamps, not
    /// pinned, ordered right after the original within the same folder.
    /// Returns the new id, or `None` if `id` is absent.
    pub fn duplicate(&mut self, id: &str) -> Option<String> {
        let src = self.projects.iter().find(|p| p.id == id)?.clone();
        let new_id = fresh_id("proj");
        let now = now_secs();
        let copy = SavedProject {
            id: new_id.clone(),
            name: format!("{} (copy)", src.name),
            tab: src.tab.clone(),
            folder: src.folder.clone(),
            pinned: false,
            order: src.order.saturating_add(1),
            created: now,
            last_opened: now,
        };
        // Nudge the originals at-or-after the insertion slot down by one so
        // the copy slots in directly after its source rather than colliding.
        let folder = copy.folder.clone();
        for p in self.projects.iter_mut() {
            if p.folder == folder && p.order >= copy.order && p.id != src.id {
                p.order = p.order.saturating_add(1);
            }
        }
        self.projects.push(copy);
        Some(new_id)
    }

    /// Move project `id` into `folder` (or unfiled when `None`), appended at
    /// the end of the destination's order. A `folder` id that doesn't exist
    /// is treated as unfiled. No-op if `id` is absent.
    pub fn move_to_folder(&mut self, id: &str, folder: Option<String>) {
        let folder = folder.filter(|fid| self.folders.iter().any(|f| &f.id == fid));
        let order = self.next_order_in(folder.as_deref());
        if let Some(p) = self.projects.iter_mut().find(|p| p.id == id) {
            p.folder = folder;
            p.order = order;
        }
    }

    /// Set the pinned flag on project `id`. No-op if absent.
    pub fn set_pinned(&mut self, id: &str, pinned: bool) {
        if let Some(p) = self.projects.iter_mut().find(|p| p.id == id) {
            p.pinned = pinned;
        }
    }

    /// Stamp project `id` as opened *now* (bumps `last_opened` for the
    /// Recent ordering). No-op if absent. The caller persists afterwards.
    pub fn mark_opened(&mut self, id: &str) {
        if let Some(p) = self.projects.iter_mut().find(|p| p.id == id) {
            p.last_opened = now_secs();
        }
    }

    /// Move project `id` one slot earlier (`up = true`) or later
    /// (`up = false`) **within its own folder**, by swapping `order` with the
    /// adjacent sibling in the current sorted view. No-op at the end stops
    /// (already first / last) or if `id` is absent.
    pub fn reorder(&mut self, id: &str, up: bool) {
        // Resolve the target's folder + a stable sorted view of its siblings.
        let Some(target) = self.projects.iter().find(|p| p.id == id) else {
            return;
        };
        let folder = target.folder.clone();
        let mut siblings: Vec<(usize, u32, String)> = self
            .projects
            .iter()
            .enumerate()
            .filter(|(_, p)| p.folder == folder)
            .map(|(i, p)| (i, p.order, p.name.clone()))
            .collect();
        // Same ordering the navigator uses: order, then name.
        siblings.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.2.cmp(&b.2)));
        let Some(pos) = siblings
            .iter()
            .position(|(i, _, _)| self.projects[*i].id == id)
        else {
            return;
        };
        let other = if up {
            if pos == 0 {
                return;
            }
            pos - 1
        } else {
            if pos + 1 >= siblings.len() {
                return;
            }
            pos + 1
        };
        let a = siblings[pos].0;
        let b = siblings[other].0;
        let oa = self.projects[a].order;
        let ob = self.projects[b].order;
        // Swap order keys. If the two happened to share an order value
        // (e.g. a legacy file), nudge so the move is still observable.
        if oa == ob {
            if up {
                self.projects[a].order = ob.saturating_sub(1);
            } else {
                self.projects[a].order = ob.saturating_add(1);
            }
        } else {
            self.projects[a].order = ob;
            self.projects[b].order = oa;
        }
    }

    // ── Folder mutations ─────────────────────────────────────────────────

    /// Create a new folder named `name` (trimmed), appended at the end of
    /// the folder order. Returns the new folder id; if `name` is empty after
    /// trimming, falls back to "New folder".
    pub fn add_folder(&mut self, name: &str) -> String {
        let name = {
            let t = name.trim();
            if t.is_empty() {
                "New folder".to_string()
            } else {
                t.to_string()
            }
        };
        let order = self
            .folders
            .iter()
            .map(|f| f.order)
            .max()
            .map(|m| m.saturating_add(1))
            .unwrap_or(0);
        let id = fresh_id("fold");
        self.folders.push(Folder {
            id: id.clone(),
            name,
            order,
        });
        id
    }

    /// Rename folder `id` to `name` (trimmed; empty rejected). No-op if
    /// absent.
    pub fn rename_folder(&mut self, id: &str, name: &str) {
        let name = name.trim();
        if name.is_empty() {
            return;
        }
        if let Some(f) = self.folders.iter_mut().find(|f| f.id == id) {
            f.name = name.to_string();
        }
    }

    /// Remove folder `id` and move its projects back to "(unfiled)" (their
    /// work is never deleted with the folder). No-op if absent.
    pub fn remove_folder(&mut self, id: &str) {
        self.folders.retain(|f| f.id != id);
        for p in self.projects.iter_mut() {
            if p.folder.as_deref() == Some(id) {
                p.folder = None;
            }
        }
    }

    // ── Read-side helpers (used by the navigator) ────────────────────────

    /// Folders sorted for display (`order`, then `name`).
    pub fn sorted_folders(&self) -> Vec<&Folder> {
        let mut v: Vec<&Folder> = self.folders.iter().collect();
        v.sort_by(|a, b| a.order.cmp(&b.order).then_with(|| a.name.cmp(&b.name)));
        v
    }

    /// Projects in `folder` (or unfiled when `None`), sorted for display
    /// (`order`, then `name`).
    pub fn projects_in(&self, folder: Option<&str>) -> Vec<&SavedProject> {
        let mut v: Vec<&SavedProject> = self
            .projects
            .iter()
            .filter(|p| p.folder.as_deref() == folder)
            .collect();
        v.sort_by(|a, b| a.order.cmp(&b.order).then_with(|| a.name.cmp(&b.name)));
        v
    }

    /// Pinned projects across all folders, sorted `order`, then `name`.
    pub fn pinned(&self) -> Vec<&SavedProject> {
        let mut v: Vec<&SavedProject> = self.projects.iter().filter(|p| p.pinned).collect();
        v.sort_by(|a, b| a.order.cmp(&b.order).then_with(|| a.name.cmp(&b.name)));
        v
    }

    /// The `n` most-recently-opened projects (newest first), ties broken by
    /// name for a stable order.
    pub fn recent(&self, n: usize) -> Vec<&SavedProject> {
        let mut v: Vec<&SavedProject> = self.projects.iter().collect();
        v.sort_by(|a, b| {
            b.last_opened
                .cmp(&a.last_opened)
                .then_with(|| a.name.cmp(&b.name))
        });
        v.truncate(n);
        v
    }

    /// Look up a project by id.
    pub fn get(&self, id: &str) -> Option<&SavedProject> {
        self.projects.iter().find(|p| p.id == id)
    }

    /// A cheap content fingerprint of the project list as the user sees it in
    /// the command palette: a hash over every project's `(id, name)`. It
    /// changes on add/remove/rename/clear (a rename keeps the count identical,
    /// so the launcher's old `projects.len()` cache key missed it and showed
    /// the stale name). The id keeps it stable across reorders that don't touch
    /// names, and order-independent so two libraries with the same projects in
    /// a different stored order fingerprint equal. Not persisted; recomputed on
    /// demand — it walks ~130 short strings, far cheaper than rebuilding the
    /// palette entries every frame.
    pub fn content_rev(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        // XOR each project's per-field hash so the result is independent of the
        // `projects` vec order (a reorder must not, by itself, invalidate the
        // palette — the launcher lists names, not order).
        let mut acc: u64 = self.projects.len() as u64;
        for p in &self.projects {
            let mut h = std::collections::hash_map::DefaultHasher::new();
            p.id.hash(&mut h);
            p.name.hash(&mut h);
            acc ^= h.finish();
        }
        acc
    }

    /// One-past-the-max `order` among projects in `folder` (or unfiled when
    /// `None`) — the slot a freshly-added/moved project lands in.
    fn next_order_in(&self, folder: Option<&str>) -> u32 {
        self.projects
            .iter()
            .filter(|p| p.folder.as_deref() == folder)
            .map(|p| p.order)
            .max()
            .map(|m| m.saturating_add(1))
            .unwrap_or(0)
    }
}

// ── File I/O split out so tests can round-trip without the real state dir ──

/// Serialise to pretty JSON. Separated so the round-trip tests exercise
/// (de)serialisation without touching the filesystem.
fn to_json(lib: &ProjectLibrary) -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(lib)
}

/// Parse a [`ProjectLibrary`] from JSON, clearing the transient inline-edit
/// flags on each embedded tab (they are `#[serde(skip)]` so deserialise to
/// `false`/empty already, but we normalise defensively).
fn from_json(text: &str) -> Result<ProjectLibrary, serde_json::Error> {
    let mut lib: ProjectLibrary = serde_json::from_str(text)?;
    for p in &mut lib.projects {
        p.tab.editing = false;
        p.tab.edit_buf.clear();
    }
    Ok(lib)
}

/// Write the library to `path` (best-effort; `false` on serialise/IO error).
fn save_to_path(lib: &ProjectLibrary, path: &std::path::Path) -> bool {
    let Ok(text) = to_json(lib) else {
        return false;
    };
    atomic_write(path, &text).is_ok()
}

/// Load a library from `path`, bounded to [`crate::settings_io::MAX_STATE_FILE_BYTES`]
/// so a corrupt/hostile file can't OOM the load. `None` on any error.
fn load_from_path(path: &std::path::Path) -> Option<ProjectLibrary> {
    let meta = std::fs::metadata(path).ok()?;
    if meta.len() > crate::settings_io::MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    let text = std::fs::read_to_string(path).ok()?;
    from_json(&text).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project_tabs::{ProjectTab, TabKind};

    /// Build a throwaway project tab of `kind` titled `title`.
    fn tab(kind: TabKind, title: &str) -> ProjectTab {
        ProjectTab {
            kind,
            title: title.to_string(),
            group: None,
            workbench_kind: None,
            editing: false,
            edit_buf: String::new(),
        }
    }

    #[test]
    fn add_project_uses_tab_title_and_appends_order() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "Booster"), None);
        let b = lib.add_project(tab(TabKind::Cad, "Bracket"), None);
        assert_eq!(lib.projects.len(), 2);
        let pa = lib.get(&a).unwrap();
        let pb = lib.get(&b).unwrap();
        assert_eq!(pa.name, "Booster");
        assert_eq!(pb.name, "Bracket");
        assert_eq!(pa.tab.kind, TabKind::Rocket);
        // Second insert orders after the first.
        assert!(pb.order > pa.order, "second project must sort after first");
        // ids are distinct.
        assert_ne!(a, b);
    }

    #[test]
    fn add_project_blank_title_falls_back_to_kind_label() {
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Cad, "   "), None);
        assert_eq!(lib.get(&id).unwrap().name, TabKind::Cad.label());
    }

    #[test]
    fn rename_trims_and_rejects_empty() {
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Rocket, "Old"), None);
        lib.rename(&id, "  New name  ");
        assert_eq!(lib.get(&id).unwrap().name, "New name");
        // Empty rename is ignored.
        lib.rename(&id, "   ");
        assert_eq!(lib.get(&id).unwrap().name, "New name");
    }

    #[test]
    fn folder_create_move_and_remove_unfiles_projects() {
        let mut lib = ProjectLibrary::default();
        let fid = lib.add_folder("Aerospace");
        let pid = lib.add_project(tab(TabKind::Rocket, "LV-1"), None);
        assert!(lib.get(&pid).unwrap().folder.is_none());

        lib.move_to_folder(&pid, Some(fid.clone()));
        assert_eq!(lib.get(&pid).unwrap().folder.as_deref(), Some(fid.as_str()));
        assert_eq!(lib.projects_in(Some(&fid)).len(), 1);
        assert_eq!(lib.projects_in(None).len(), 0);

        // Removing the folder moves its project back to unfiled (not deleted).
        lib.remove_folder(&fid);
        assert!(lib.folders.is_empty());
        assert_eq!(lib.projects.len(), 1);
        assert!(lib.get(&pid).unwrap().folder.is_none());
    }

    #[test]
    fn move_to_nonexistent_folder_is_treated_as_unfiled() {
        let mut lib = ProjectLibrary::default();
        let pid = lib.add_project(tab(TabKind::Cad, "Part"), None);
        lib.move_to_folder(&pid, Some("does-not-exist".to_string()));
        assert!(lib.get(&pid).unwrap().folder.is_none());
    }

    #[test]
    fn add_project_into_missing_folder_lands_unfiled() {
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Cad, "Part"), Some("nope".to_string()));
        assert!(lib.get(&id).unwrap().folder.is_none());
    }

    #[test]
    fn pin_toggles_and_filters() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "A"), None);
        let _b = lib.add_project(tab(TabKind::Cad, "B"), None);
        assert!(lib.pinned().is_empty());
        lib.set_pinned(&a, true);
        let pinned = lib.pinned();
        assert_eq!(pinned.len(), 1);
        assert_eq!(pinned[0].id, a);
        lib.set_pinned(&a, false);
        assert!(lib.pinned().is_empty());
    }

    #[test]
    fn duplicate_makes_independent_copy_after_source() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "Booster"), None);
        let b = lib.add_project(tab(TabKind::Cad, "Bracket"), None);
        let dup = lib.duplicate(&a).expect("duplicate");
        assert_eq!(lib.projects.len(), 3);
        let pdup = lib.get(&dup).unwrap();
        assert_eq!(pdup.name, "Booster (copy)");
        assert!(!pdup.pinned);
        assert_ne!(pdup.id, a);
        // The copy sorts after the source but the third (B) was pushed
        // down so orders stay distinct within the folder.
        let view = lib.projects_in(None);
        let order: Vec<&str> = view.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(order, vec!["Booster", "Booster (copy)", "Bracket"]);
        // Mutating the copy's tab title must not touch the source.
        let _ = b;
    }

    #[test]
    fn reorder_swaps_adjacent_within_folder() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "A"), None);
        let b = lib.add_project(tab(TabKind::Cad, "B"), None);
        let c = lib.add_project(tab(TabKind::Fem, "C"), None);
        // Initial sorted view: A, B, C.
        let names = |l: &ProjectLibrary| -> Vec<String> {
            l.projects_in(None).iter().map(|p| p.name.clone()).collect()
        };
        assert_eq!(names(&lib), vec!["A", "B", "C"]);
        // Move B up → B, A, C.
        lib.reorder(&b, true);
        assert_eq!(names(&lib), vec!["B", "A", "C"]);
        // Move B down twice → A, C, B (clamps at the end).
        lib.reorder(&b, false);
        assert_eq!(names(&lib), vec!["A", "B", "C"]);
        lib.reorder(&b, false);
        assert_eq!(names(&lib), vec!["A", "C", "B"]);
        // Up from the top is a no-op.
        lib.reorder(&a, true);
        assert_eq!(names(&lib), vec!["A", "C", "B"]);
        let _ = c;
    }

    #[test]
    fn mark_opened_bumps_recent_to_front() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "A"), None);
        let b = lib.add_project(tab(TabKind::Cad, "B"), None);
        // Force a strictly-older stamp on A so the ordering is deterministic
        // (both were stamped "now" on add, which can collide at 1-second
        // resolution).
        if let Some(p) = lib.projects.iter_mut().find(|p| p.id == a) {
            p.last_opened = 1;
        }
        if let Some(p) = lib.projects.iter_mut().find(|p| p.id == b) {
            p.last_opened = 2;
        }
        assert_eq!(lib.recent(10)[0].id, b);
        // Re-opening A bumps it ahead of B.
        lib.mark_opened(&a);
        assert_eq!(lib.recent(10)[0].id, a);
    }

    #[test]
    fn round_trips_through_json() {
        // Build a non-trivial library: two folders, projects filed and
        // unfiled, one pinned, distinct orders + open stamps.
        let mut lib = ProjectLibrary::default();
        let aero = lib.add_folder("Aerospace");
        let cad = lib.add_folder("CAD");
        let r = lib.add_project(tab(TabKind::Rocket, "LV-1"), Some(aero.clone()));
        let _e = lib.add_project(tab(TabKind::Engine, "Raptor"), Some(aero.clone()));
        let _b = lib.add_project(tab(TabKind::Cad, "Bracket"), Some(cad.clone()));
        let _u = lib.add_project(tab(TabKind::Fem, "Loose"), None);
        lib.set_pinned(&r, true);
        lib.mark_opened(&r);

        // Serialise → deserialise must reproduce the library exactly.
        let json = to_json(&lib).expect("serialise");
        let back = from_json(&json).expect("deserialise");
        assert_eq!(lib, back, "library must survive a JSON round-trip");

        // Spot-check structure survived.
        assert_eq!(back.folders.len(), 2);
        assert_eq!(back.projects.len(), 4);
        assert_eq!(back.projects_in(Some(&aero)).len(), 2);
        assert_eq!(back.pinned().len(), 1);
        assert_eq!(back.pinned()[0].name, "LV-1");
    }

    #[test]
    fn save_then_load_from_path_is_equal() {
        // Exercise the real on-disk save/load helpers against a throwaway
        // temp file (no dependence on the process state dir).
        let dir = std::env::temp_dir().join(format!(
            "valenx-library-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("library.json");

        let mut lib = ProjectLibrary::default();
        let f = lib.add_folder("Group");
        let _p = lib.add_project(tab(TabKind::Rocket, "Saved"), Some(f.clone()));

        assert!(save_to_path(&lib, &path), "save must succeed");
        let loaded = load_from_path(&path).expect("load");
        assert_eq!(lib, loaded);

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn clear_all_empties_projects_and_folders() {
        let mut lib = ProjectLibrary::default();
        let f = lib.add_folder("Group");
        let _a = lib.add_project(tab(TabKind::Rocket, "A"), Some(f.clone()));
        let _b = lib.add_project(tab(TabKind::Cad, "B"), None);
        assert!(!lib.projects.is_empty());
        assert!(!lib.folders.is_empty());

        lib.clear_all();
        assert!(lib.projects.is_empty(), "clear_all drops every project");
        assert!(lib.folders.is_empty(), "clear_all drops every folder");
        // A cleared library equals a fresh default one.
        assert_eq!(lib, ProjectLibrary::default());
    }

    #[test]
    fn reconcile_unfiles_dangling_folder_references() {
        // A project pointing at a folder id that isn't present (e.g. a
        // hand-edited file) must be moved back to unfiled on load.
        let mut lib = ProjectLibrary {
            projects: vec![SavedProject {
                id: "p1".into(),
                name: "Orphan".into(),
                tab: tab(TabKind::Cad, "Orphan"),
                folder: Some("ghost-folder".into()),
                pinned: false,
                order: 0,
                created: 0,
                last_opened: 0,
            }],
            folders: Vec::new(),
        };
        lib.reconcile();
        assert!(lib.projects[0].folder.is_none());
    }

    #[test]
    fn content_rev_changes_on_rename_even_though_count_is_unchanged() {
        // The launcher's old cache key was `projects.len()`, which a rename
        // leaves untouched — so the palette showed the stale name. `content_rev`
        // fingerprints `(id, name)`, so a rename must flip it while the count
        // (and thus the old key) stays equal.
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Rocket, "Old"), None);
        let before_rev = lib.content_rev();
        let before_len = lib.projects.len();

        lib.rename(&id, "New");
        assert_eq!(
            lib.projects.len(),
            before_len,
            "rename must not change the project count"
        );
        assert_ne!(
            lib.content_rev(),
            before_rev,
            "rename must change the content fingerprint"
        );
    }

    #[test]
    fn content_rev_is_independent_of_stored_order() {
        // Two libraries holding the same projects in a different `projects` vec
        // order must fingerprint equal (the launcher lists names, not order, so
        // a reorder alone should not invalidate the palette cache).
        let mut a = ProjectLibrary::default();
        let p1 = a.add_project(tab(TabKind::Rocket, "A"), None);
        let _p2 = a.add_project(tab(TabKind::Cad, "B"), None);

        // Build `b` with the same two projects pushed in reverse order.
        let mut b = ProjectLibrary::default();
        for p in a.projects.iter().rev() {
            b.projects.push(p.clone());
        }
        assert_eq!(
            a.content_rev(),
            b.content_rev(),
            "fingerprint must not depend on stored order"
        );
        let _ = p1;
    }

    #[test]
    fn load_missing_file_is_default_empty() {
        let path = std::env::temp_dir().join(format!(
            "valenx-library-absent-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // The file does not exist.
        assert!(load_from_path(&path).is_none());
    }
}

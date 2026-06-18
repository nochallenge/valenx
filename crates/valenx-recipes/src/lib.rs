//! # valenx-recipes — a reusable design-recipe library
//!
//! A small, dependency-light store for **design recipes**: a workbench
//! produces a design (a rocket, an engine, a CAD part, a gear train, …),
//! serializes it into a [`Recipe`], and saves it to disk; a later
//! session — or a different workbench — lists and loads recipes back to
//! reuse or compose them. This is what lets an AI-generated design
//! *persist* and be reused **modularly** rather than evaporating when the
//! app closes.
//!
//! ## What a recipe is
//!
//! A [`Recipe`] is a portable envelope:
//!
//! - `name` — a human label (`"LV-1 baseline"`);
//! - `domain` — which kind of design it is (`"rocket"`, `"engine"`, …),
//!   used to group recipes and to route a recipe back to the workbench
//!   that understands it;
//! - `created_utc` — an ISO-8601 timestamp string (supplied by the
//!   caller, so this crate stays free of any clock dependency and is
//!   trivially testable);
//! - `spec` — a free-form [`serde_json::Value`] holding the actual
//!   design parameters. Keeping the payload generic JSON is deliberate:
//!   the store never needs to know a workbench's schema, so any domain
//!   can save and load recipes without touching this crate.
//!
//! ## The store
//!
//! [`RecipeStore`] is a thin wrapper over a directory. [`save`] writes one
//! pretty-printed JSON file per recipe (filename derived from the domain
//! and a slug of the name); [`list`] reads every `*.json` back, sorted by
//! name; [`load`] reads a single file. A missing directory simply lists
//! as empty.
//!
//! [`save`]: RecipeStore::save
//! [`list`]: RecipeStore::list
//! [`load`]: RecipeStore::load
//!
//! ```
//! use valenx_recipes::{Recipe, RecipeStore};
//!
//! let dir = std::env::temp_dir().join("valenx_recipes_doctest");
//! let _ = std::fs::remove_dir_all(&dir);
//! let store = RecipeStore::new(&dir);
//!
//! let spec = serde_json::json!({ "chamber_pressure_bar": 97.0, "expansion_ratio": 16.0 });
//! let recipe = Recipe::new("LV-1 baseline", "rocket", "2026-06-18T00:00:00Z", spec).unwrap();
//! let path = store.save(&recipe).unwrap();
//! assert!(path.exists());
//!
//! let loaded = RecipeStore::load(&path).unwrap();
//! assert_eq!(loaded, recipe);
//! assert_eq!(store.list().unwrap().len(), 1);
//! # let _ = std::fs::remove_dir_all(&dir);
//! ```
//!
//! ## Honest scope
//!
//! This crate is a **storage and serialization layer only**. It does not
//! generate, validate, simulate, or interpret designs — interpreting a
//! recipe's `spec` is the job of the workbench that created it. It makes
//! no claim that a saved recipe is physically valid; that is the
//! originating tool's responsibility.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// Errors returned by the recipe store. Each carries enough context
/// (the offending path) to surface a useful message in the UI.
#[derive(Debug, thiserror::Error)]
pub enum RecipeError {
    /// A recipe was created with an empty `name`.
    #[error("recipe name must not be empty")]
    EmptyName,
    /// A recipe was created with an empty `domain`.
    #[error("recipe domain must not be empty")]
    EmptyDomain,
    /// A filesystem operation failed.
    #[error("i/o error at {path:?}: {source}")]
    Io {
        /// The path the operation was attempted on.
        path: PathBuf,
        /// The underlying I/O error.
        source: std::io::Error,
    },
    /// A recipe file could not be (de)serialized as JSON.
    #[error("JSON error at {path:?}: {source}")]
    Json {
        /// The recipe file involved.
        path: PathBuf,
        /// The underlying serde_json error.
        source: serde_json::Error,
    },
}

/// Convenience result alias for the recipe API.
pub type Result<T> = std::result::Result<T, RecipeError>;

/// A saved design recipe: a portable, JSON-serializable envelope around a
/// workbench's design output. See the [crate-level docs](crate) for the
/// field semantics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Recipe {
    /// Human-facing label, e.g. `"LV-1 baseline"`.
    pub name: String,
    /// Design domain, e.g. `"rocket"` or `"engine"` — groups recipes and
    /// routes one back to the workbench that understands it.
    pub domain: String,
    /// ISO-8601 creation timestamp, supplied by the caller.
    pub created_utc: String,
    /// Free-form design payload; opaque to this crate.
    pub spec: serde_json::Value,
}

impl Recipe {
    /// Build a recipe, rejecting an empty `name` or `domain` (both are
    /// trimmed before the check, so whitespace-only is also rejected).
    pub fn new(
        name: impl Into<String>,
        domain: impl Into<String>,
        created_utc: impl Into<String>,
        spec: serde_json::Value,
    ) -> Result<Self> {
        let name = name.into();
        let domain = domain.into();
        if name.trim().is_empty() {
            return Err(RecipeError::EmptyName);
        }
        if domain.trim().is_empty() {
            return Err(RecipeError::EmptyDomain);
        }
        Ok(Self {
            name,
            domain,
            created_utc: created_utc.into(),
            spec,
        })
    }

    /// The filename stem this recipe saves under: `"<domain>__<name>"`,
    /// each part slugged to lowercase ASCII alphanumerics joined by `-`.
    /// Two recipes with the same domain + name share a stem and so
    /// overwrite each other (a deliberate "named slot" semantics).
    pub fn file_stem(&self) -> String {
        format!("{}__{}", slug(&self.domain), slug(&self.name))
    }
}

/// Slug a label to a filesystem-safe token: lowercase ASCII
/// alphanumerics kept as-is, every other run collapsed to a single `-`,
/// leading/trailing `-` trimmed. An empty result falls back to `"item"`
/// so a stem is never blank.
fn slug(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut prev_dash = false;
    for ch in s.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            prev_dash = false;
        } else if !prev_dash {
            out.push('-');
            prev_dash = true;
        }
    }
    let trimmed = out.trim_matches('-');
    if trimmed.is_empty() {
        "item".to_string()
    } else {
        trimmed.to_string()
    }
}

/// A filesystem-backed recipe store rooted at a single directory.
pub struct RecipeStore {
    root: PathBuf,
}

impl RecipeStore {
    /// Create a store rooted at `root`. The directory is created lazily
    /// on the first [`save`](Self::save); it need not exist yet.
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self { root: root.into() }
    }

    /// The directory this store reads and writes.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// Save `recipe` as a pretty-printed JSON file under the store's
    /// directory (creating the directory if needed). Returns the path
    /// written. A recipe with the same domain + name overwrites the
    /// previous file (see [`Recipe::file_stem`]).
    pub fn save(&self, recipe: &Recipe) -> Result<PathBuf> {
        fs::create_dir_all(&self.root).map_err(|source| RecipeError::Io {
            path: self.root.clone(),
            source,
        })?;
        let path = self.root.join(format!("{}.json", recipe.file_stem()));
        let json = serde_json::to_string_pretty(recipe).map_err(|source| RecipeError::Json {
            path: path.clone(),
            source,
        })?;
        fs::write(&path, json).map_err(|source| RecipeError::Io {
            path: path.clone(),
            source,
        })?;
        Ok(path)
    }

    /// Load every `*.json` recipe in the store, sorted by name. A missing
    /// store directory lists as empty rather than erroring; a malformed
    /// file is a hard error (it names the bad path).
    pub fn list(&self) -> Result<Vec<Recipe>> {
        if !self.root.exists() {
            return Ok(Vec::new());
        }
        let entries = fs::read_dir(&self.root).map_err(|source| RecipeError::Io {
            path: self.root.clone(),
            source,
        })?;
        let mut recipes = Vec::new();
        for entry in entries {
            let entry = entry.map_err(|source| RecipeError::Io {
                path: self.root.clone(),
                source,
            })?;
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) == Some("json") {
                recipes.push(Self::load(&path)?);
            }
        }
        recipes.sort_by(|a, b| a.name.cmp(&b.name));
        Ok(recipes)
    }

    /// Load a single recipe from a JSON file.
    pub fn load(path: impl AsRef<Path>) -> Result<Recipe> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|source| RecipeError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        serde_json::from_str(&text).map_err(|source| RecipeError::Json {
            path: path.to_path_buf(),
            source,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_root(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("valenx_recipes_test_{tag}"));
        let _ = fs::remove_dir_all(&dir);
        dir
    }

    fn sample(name: &str) -> Recipe {
        Recipe::new(
            name,
            "rocket",
            "2026-06-18T00:00:00Z",
            serde_json::json!({ "chamber_pressure_bar": 97.0, "expansion_ratio": 16.0 }),
        )
        .unwrap()
    }

    #[test]
    fn new_rejects_empty_name_and_domain() {
        let spec = serde_json::json!({});
        assert!(matches!(
            Recipe::new("  ", "rocket", "t", spec.clone()),
            Err(RecipeError::EmptyName)
        ));
        assert!(matches!(
            Recipe::new("ok", "", "t", spec),
            Err(RecipeError::EmptyDomain)
        ));
    }

    #[test]
    fn slug_is_filesystem_safe() {
        assert_eq!(slug("LV-1 baseline"), "lv-1-baseline");
        assert_eq!(slug("  Heavy/Lift  "), "heavy-lift");
        assert_eq!(slug("***"), "item");
        assert_eq!(slug("Already_ok123"), "already-ok123");
    }

    #[test]
    fn file_stem_combines_domain_and_name() {
        let r = sample("LV-1 baseline");
        assert_eq!(r.file_stem(), "rocket__lv-1-baseline");
    }

    #[test]
    fn save_then_load_round_trips() {
        let store = RecipeStore::new(temp_root("roundtrip"));
        let r = sample("LV-1 baseline");
        let path = store.save(&r).unwrap();
        assert!(path.exists());
        let loaded = RecipeStore::load(&path).unwrap();
        assert_eq!(loaded, r);
        let _ = fs::remove_dir_all(store.root());
    }

    #[test]
    fn save_then_list_returns_sorted_recipes() {
        let store = RecipeStore::new(temp_root("list"));
        store.save(&sample("Zeta")).unwrap();
        store.save(&sample("Alpha")).unwrap();
        let listed = store.list().unwrap();
        let names: Vec<&str> = listed.iter().map(|r| r.name.as_str()).collect();
        assert_eq!(names, ["Alpha", "Zeta"], "listed sorted by name");
        let _ = fs::remove_dir_all(store.root());
    }

    #[test]
    fn same_domain_and_name_overwrites_the_slot() {
        let store = RecipeStore::new(temp_root("overwrite"));
        store.save(&sample("Falcon")).unwrap();
        // Re-save same domain+name with a different spec.
        let mut r = sample("Falcon");
        r.spec = serde_json::json!({ "chamber_pressure_bar": 250.0 });
        store.save(&r).unwrap();
        let listed = store.list().unwrap();
        assert_eq!(listed.len(), 1, "same slot overwrites, not duplicates");
        assert_eq!(
            listed[0].spec,
            serde_json::json!({ "chamber_pressure_bar": 250.0 })
        );
        let _ = fs::remove_dir_all(store.root());
    }

    #[test]
    fn list_on_missing_directory_is_empty() {
        let store = RecipeStore::new(temp_root("missing_never_created"));
        assert!(store.list().unwrap().is_empty());
    }
}

//! Add-on registry — keeps the in-memory list of installed and remote
//! add-ons and dispatches install/update/uninstall.
//!
//! The desktop shell holds a single `AddonRegistry` and renders it
//! into the Add-on Manager panel.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::AddonError;
use crate::install;
use crate::manifest::{read_manifest_at, AddonManifest};

/// Round-19 L1 cap on the number of subdirs we'll walk under
/// `install_dir`. A typical install carries 0–50 add-ons; 10k is
/// well past anything sane, while still catching a poisoned dir
/// that holds millions of placeholder subdirs (the read_dir walk
/// would otherwise allocate a Vec per entry before
/// `read_manifest_at` could refuse any of them).
pub const MAX_ADDONS: usize = 10_000;

/// A remote add-on candidate — what a network search would return if
/// the GitHub path were wired. v1 builds these manually from user
/// input ("paste a repo URL").
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RemoteAddon {
    /// User-friendly name (typically the repo name).
    pub name: String,
    /// One-line description.
    pub description: String,
    /// Source URL — typically `https://github.com/owner/repo`.
    pub source_url: String,
    /// Latest tag / version, if known.
    pub version: String,
}

impl RemoteAddon {
    /// Build a minimal candidate from a URL alone, for the "paste a
    /// URL" install flow when the user hasn't filled in name /
    /// version.
    pub fn from_url(url: impl Into<String>) -> Self {
        let url = url.into();
        let name = url
            .rsplit('/')
            .next()
            .unwrap_or("addon")
            .trim_end_matches(".git")
            .to_string();
        RemoteAddon {
            name,
            description: String::new(),
            source_url: url,
            version: "unknown".to_string(),
        }
    }
}

/// Installed add-on record — pairs the parsed manifest with the
/// on-disk path.
#[derive(Clone, Debug)]
pub struct LocalAddon {
    /// Manifest contents.
    pub manifest: AddonManifest,
    /// Filesystem path to the install directory.
    pub path: PathBuf,
}

/// Global registry — owned by the desktop shell, ticked when the user
/// opens the Add-on Manager panel.
#[derive(Clone, Debug)]
pub struct AddonRegistry {
    /// Top-level install directory (e.g. `~/.valenx/addons/`).
    pub install_dir: PathBuf,
    installed: Vec<LocalAddon>,
}

impl Default for AddonRegistry {
    fn default() -> Self {
        // Default to `~/.valenx/addons/` when the home directory is
        // resolvable; otherwise a sentinel path the desktop shell can
        // overwrite via [`AddonRegistry::set_install_dir`].
        let dir = std::env::var_os("USERPROFILE")
            .or_else(|| std::env::var_os("HOME"))
            .map(|h| PathBuf::from(h).join(".valenx").join("addons"))
            .unwrap_or_else(|| PathBuf::from("addons"));
        Self::new(dir)
    }
}

impl AddonRegistry {
    /// Build a fresh registry pointing at `install_dir`.
    pub fn new(install_dir: PathBuf) -> Self {
        Self {
            install_dir,
            installed: Vec::new(),
        }
    }

    /// Override the install directory at runtime.
    pub fn set_install_dir(&mut self, dir: PathBuf) {
        self.install_dir = dir;
    }

    /// Re-scan the install directory and rebuild the in-memory list.
    ///
    /// # Errors
    ///
    /// - [`AddonError::Io`] on read failure other than `NotFound`.
    /// - [`AddonError::Manifest`] / [`AddonError::InvalidManifest`] for
    ///   any add-on whose manifest is unreadable / malformed (logged
    ///   via `tracing` and skipped in the returned list; the function
    ///   itself succeeds).
    pub fn refresh(&mut self) -> Result<(), AddonError> {
        self.installed.clear();
        if !self.install_dir.exists() {
            return Ok(());
        }
        // Round-19 L1: cap the read_dir iteration at MAX_ADDONS. A
        // poisoned install dir with millions of subdirs would
        // otherwise let the refresh loop run unbounded — every
        // `read_manifest_at` would do a stat + manifest read before
        // the loop returned, freezing the Add-on Manager panel.
        for entry in std::fs::read_dir(&self.install_dir)?.take(MAX_ADDONS) {
            let entry = entry?;
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }
            let manifest_path = path.join("valenx-addon.toml");
            if !manifest_path.exists() {
                continue;
            }
            match read_manifest_at(&manifest_path) {
                Ok(manifest) => {
                    self.installed.push(LocalAddon { manifest, path });
                }
                Err(e) => {
                    tracing::warn!(
                        target: "valenx-addons",
                        "skipping addon at {}: {}",
                        path.display(),
                        e,
                    );
                }
            }
        }
        self.installed
            .sort_by(|a, b| a.manifest.name.cmp(&b.manifest.name));
        Ok(())
    }

    /// Currently-installed add-ons.
    pub fn installed(&self) -> &[LocalAddon] {
        &self.installed
    }

    /// Find an installed add-on by name.
    pub fn find(&self, name: &str) -> Option<&LocalAddon> {
        self.installed.iter().find(|a| a.manifest.name == name)
    }

    /// Install from a local source directory + refresh.
    ///
    /// # Errors
    ///
    /// Whatever [`install::install_from_dir`] returns.
    pub fn install_from_dir(&mut self, source_dir: &Path) -> Result<LocalAddon, AddonError> {
        let local = install::install_from_dir(source_dir, &self.install_dir)?;
        self.refresh()?;
        Ok(local)
    }

    /// Update a named add-on from a local source directory + refresh.
    ///
    /// # Errors
    ///
    /// Whatever [`install::update_from_dir`] returns.
    pub fn update_from_dir(&mut self, source_dir: &Path) -> Result<LocalAddon, AddonError> {
        let local = install::update_from_dir(source_dir, &self.install_dir)?;
        self.refresh()?;
        Ok(local)
    }

    /// Uninstall an add-on by name + refresh. Returns `Ok(false)` if
    /// the add-on wasn't installed.
    ///
    /// # Errors
    ///
    /// Whatever [`install::uninstall`] returns.
    pub fn uninstall(&mut self, name: &str) -> Result<bool, AddonError> {
        let Some(local) = self.find(name).cloned() else {
            return Ok(false);
        };
        let removed = install::uninstall(&local.path)?;
        self.refresh()?;
        Ok(removed)
    }

    /// Placeholder for the future GitHub-search code path. Returns
    /// [`AddonError::GithubSearchUnavailable`] in v1 so callers can
    /// route to the manual install flow.
    ///
    /// # Errors
    ///
    /// Always returns [`AddonError::GithubSearchUnavailable`] in v1.
    pub fn github_search(&self, _query: &str) -> Result<Vec<RemoteAddon>, AddonError> {
        Err(AddonError::GithubSearchUnavailable)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(name: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("valenx_addons_reg_{name}"));
        let _ = std::fs::remove_dir_all(&p);
        p
    }

    fn write_manifest(path: &Path, name: &str) {
        std::fs::create_dir_all(path).unwrap();
        std::fs::write(
            path.join("valenx-addon.toml"),
            format!(
                r#"
                    name = "{name}"
                    description = "{name} desc"
                    version = "1.2.3"

                    [entry_point]
                    kind = "python"
                    module = "main"
                "#
            ),
        )
        .unwrap();
    }

    #[test]
    fn empty_registry_lists_nothing() {
        let dir = tmpdir("empty");
        let mut r = AddonRegistry::new(dir.clone());
        r.refresh().unwrap();
        assert!(r.installed().is_empty());
    }

    #[test]
    fn refresh_picks_up_installed_addons() {
        let dir = tmpdir("refresh");
        std::fs::create_dir_all(&dir).unwrap();
        write_manifest(&dir.join("foo"), "foo");
        write_manifest(&dir.join("bar"), "bar");
        let mut r = AddonRegistry::new(dir.clone());
        r.refresh().unwrap();
        // Lex sort -> bar then foo.
        assert_eq!(r.installed().len(), 2);
        assert_eq!(r.installed()[0].manifest.name, "bar");
        assert_eq!(r.installed()[1].manifest.name, "foo");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn find_returns_matching_addon() {
        let dir = tmpdir("find");
        std::fs::create_dir_all(&dir).unwrap();
        write_manifest(&dir.join("alpha"), "alpha");
        let mut r = AddonRegistry::new(dir.clone());
        r.refresh().unwrap();
        assert!(r.find("alpha").is_some());
        assert!(r.find("missing").is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn github_search_returns_unavailable_in_v1() {
        let r = AddonRegistry::new(std::env::temp_dir());
        let err = r.github_search("anything").unwrap_err();
        assert!(matches!(err, AddonError::GithubSearchUnavailable));
    }

    #[test]
    fn remote_addon_from_url_extracts_name() {
        let a = RemoteAddon::from_url("https://github.com/owner/cool-addon");
        assert_eq!(a.name, "cool-addon");
        let b = RemoteAddon::from_url("https://github.com/owner/repo.git");
        assert_eq!(b.name, "repo");
    }

    #[test]
    fn install_from_dir_via_registry_picks_up_addon() {
        let install_dir = tmpdir("via_reg_install");
        let src = tmpdir("via_reg_src");
        write_manifest(&src, "via-reg");
        let mut r = AddonRegistry::new(install_dir.clone());
        r.install_from_dir(&src).unwrap();
        assert_eq!(r.installed().len(), 1);
        let _ = std::fs::remove_dir_all(&install_dir);
        let _ = std::fs::remove_dir_all(&src);
    }

    #[test]
    fn uninstall_missing_returns_false() {
        let install_dir = tmpdir("uninstall_missing");
        let mut r = AddonRegistry::new(install_dir.clone());
        assert!(!r.uninstall("nothing").unwrap());
        let _ = std::fs::remove_dir_all(&install_dir);
    }

    #[test]
    fn uninstall_removes_addon() {
        let install_dir = tmpdir("uninstall_real");
        let src = tmpdir("uninstall_src");
        write_manifest(&src, "rmme");
        let mut r = AddonRegistry::new(install_dir.clone());
        r.install_from_dir(&src).unwrap();
        assert!(r.uninstall("rmme").unwrap());
        assert!(r.installed().is_empty());
        let _ = std::fs::remove_dir_all(&install_dir);
        let _ = std::fs::remove_dir_all(&src);
    }
}

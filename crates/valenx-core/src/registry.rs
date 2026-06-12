//! The adapter registry — the runtime catalog `valenx-app` queries
//! to decide what physics is available on this install.
//!
//! Spec: [ARCHITECTURE.md § 5](../../ARCHITECTURE.md) and
//! [RFC 0002](../../rfcs/0002-adapter-contract.md).

use std::collections::HashMap;
use std::sync::Arc;

use crate::adapter::{Adapter, ProbeReport};
use crate::error::AdapterError;
use crate::physics::{Capability, Physics};

/// Status of a registered adapter after `probe()`.
#[derive(Clone, Debug)]
pub enum AdapterStatus {
    /// Present and at a supported version.
    Ready { report: ProbeReport },
    /// Tool not found on this machine.
    Missing { hint: String },
    /// Found but outside the adapter's version range.
    Outdated { expected: String, found: String },
    /// Present but probing failed for another reason.
    Broken { error: String },
    /// User explicitly disabled this adapter in Settings.
    Disabled,
}

impl AdapterStatus {
    /// `true` iff the adapter probed cleanly.
    pub fn is_ready(&self) -> bool {
        matches!(self, AdapterStatus::Ready { .. })
    }
}

/// One entry in the registry.
pub struct AdapterEntry {
    pub adapter: Arc<dyn Adapter>,
    pub status: AdapterStatus,
}

/// Runtime catalog of every integrated adapter.
#[derive(Default)]
pub struct AdapterRegistry {
    by_id: HashMap<&'static str, AdapterEntry>,
}

impl AdapterRegistry {
    /// New, empty registry.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an adapter. `status` should reflect a completed probe —
    /// call `probe_all()` afterwards if you want the registry to
    /// classify everything.
    pub fn register(&mut self, adapter: Arc<dyn Adapter>) {
        let info = adapter.info();
        self.by_id.insert(
            info.id,
            AdapterEntry {
                adapter,
                status: AdapterStatus::Broken {
                    error: "not yet probed".into(),
                },
            },
        );
    }

    /// Probe every registered adapter in sequence. Returns how many
    /// ended up `Ready`. For parallel probe under `rayon`, a caller
    /// can compose their own driver; this sequential version is
    /// always available.
    pub fn probe_all(&mut self) -> usize {
        let mut ready = 0usize;
        for entry in self.by_id.values_mut() {
            entry.status = classify(&*entry.adapter);
            if entry.status.is_ready() {
                ready += 1;
            }
        }
        ready
    }

    /// Probe every registered adapter on a BACKGROUND thread, streaming
    /// each adapter's classified status back through the returned channel
    /// as it completes. Non-blocking: the registry is not mutated here —
    /// the caller applies results with [`AdapterRegistry::apply_probe_result`]
    /// as they arrive (e.g. once per UI frame). Adapters are `Send + Sync`,
    /// so moving clones onto the thread is safe.
    ///
    /// Use this on hot startup paths instead of [`AdapterRegistry::probe_all`]:
    /// probing ~150 external tools (a PATH search + a version-spawn each),
    /// sequentially on the main thread, blocked the first frame for seconds
    /// — tens of seconds on a cold filesystem cache.
    pub fn spawn_probe_all(&self) -> std::sync::mpsc::Receiver<(&'static str, AdapterStatus)> {
        let adapters: Vec<(&'static str, Arc<dyn Adapter>)> = self
            .by_id
            .iter()
            .map(|(id, e)| (*id, e.adapter.clone()))
            .collect();
        let (tx, rx) = std::sync::mpsc::channel();
        // If the OS refuses the thread (extremely rare), `tx` drops, `rx`
        // disconnects, and adapters simply stay "not yet probed" until a
        // manual reprobe — the app still runs on its native paths.
        let _ = std::thread::Builder::new()
            .name("valenx-adapter-probe".into())
            .spawn(move || {
                for (id, adapter) in adapters {
                    // Stop early if the receiver was dropped (app closing).
                    if tx.send((id, classify(&*adapter))).is_err() {
                        break;
                    }
                }
            });
        rx
    }

    /// Apply one streamed probe result from
    /// [`AdapterRegistry::spawn_probe_all`]. Returns `true` if that adapter
    /// is now `Ready`. Unknown ids are ignored.
    pub fn apply_probe_result(&mut self, id: &str, status: AdapterStatus) -> bool {
        let ready = status.is_ready();
        if let Some(e) = self.by_id.get_mut(id) {
            e.status = status;
        }
        ready
    }

    /// Look up an adapter by its [`crate::adapter::AdapterInfo::id`].
    pub fn get(&self, id: &str) -> Option<&AdapterEntry> {
        self.by_id.get(id)
    }

    /// Total number of registered adapters (regardless of probe state).
    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    /// `true` if no adapters have been registered yet.
    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }

    /// Iterate every registered adapter with its status.
    pub fn iter(&self) -> impl Iterator<Item = (&&'static str, &AdapterEntry)> {
        self.by_id.iter()
    }

    /// All `Ready` adapters for a given physics.
    pub fn ready_for_physics(&self, p: Physics) -> Vec<&AdapterEntry> {
        self.by_id
            .values()
            .filter(|e| e.status.is_ready() && e.adapter.info().physics.contains(&p))
            .collect()
    }

    /// All `Ready` adapters advertising a given capability.
    pub fn ready_for_capability(&self, c: Capability) -> Vec<&AdapterEntry> {
        self.by_id
            .values()
            .filter(|e| e.status.is_ready() && e.adapter.capabilities().capabilities.contains(&c))
            .collect()
    }

    /// Counts — useful for the Home-screen Status strip.
    pub fn counts(&self) -> StatusCounts {
        let mut c = StatusCounts::default();
        for e in self.by_id.values() {
            match &e.status {
                AdapterStatus::Ready { .. } => c.ready += 1,
                AdapterStatus::Missing { .. } => c.missing += 1,
                AdapterStatus::Outdated { .. } => c.outdated += 1,
                AdapterStatus::Broken { .. } => c.broken += 1,
                AdapterStatus::Disabled => c.disabled += 1,
            }
        }
        c
    }
}

fn classify(adapter: &dyn Adapter) -> AdapterStatus {
    match adapter.probe() {
        Ok(report) if report.ok => {
            let info = adapter.info();
            if let Some(found) = &report.found_version {
                if !info.version_range.contains(found) {
                    return AdapterStatus::Outdated {
                        expected: format!(
                            "{}..{}",
                            info.version_range.min_inclusive, info.version_range.max_exclusive
                        ),
                        found: found.to_string(),
                    };
                }
            }
            AdapterStatus::Ready { report }
        }
        Ok(_) => AdapterStatus::Missing {
            hint: "probe returned ok = false".into(),
        },
        Err(AdapterError::ToolNotInstalled { hint, .. }) => AdapterStatus::Missing { hint },
        Err(AdapterError::ToolVersionMismatch {
            expected, found, ..
        }) => AdapterStatus::Outdated {
            expected,
            found: found.to_string(),
        },
        Err(e) => AdapterStatus::Broken {
            error: e.to_string(),
        },
    }
}

/// Summary for the Home-screen Status strip.
#[derive(Clone, Copy, Debug, Default)]
pub struct StatusCounts {
    pub ready: usize,
    pub missing: usize,
    pub outdated: usize,
    pub broken: usize,
    pub disabled: usize,
}

impl StatusCounts {
    /// True if nothing is ready — the "nothing works" home state
    /// described in DESIGN.md § 10.
    pub fn nothing_ready(&self) -> bool {
        self.ready == 0 && (self.missing + self.outdated + self.broken) > 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_registry() {
        let r = AdapterRegistry::new();
        assert_eq!(r.len(), 0);
        assert!(r.is_empty());
        assert_eq!(r.counts().ready, 0);
    }

    #[test]
    fn spawn_probe_all_on_empty_registry_completes() {
        // No adapters → the background thread sends nothing and drops the
        // sender, so the receiver disconnects. This exercises the
        // thread+channel plumbing without needing a mock Adapter.
        let r = AdapterRegistry::new();
        let rx = r.spawn_probe_all();
        // recv() blocks only until the (empty) probe thread finishes.
        assert!(
            rx.recv().is_err(),
            "empty probe should disconnect, not yield a result"
        );
    }

    #[test]
    fn apply_probe_result_ignores_unknown_id_and_reports_readiness() {
        let mut r = AdapterRegistry::new();
        // Unknown id is a no-op (no panic); return value reflects the
        // status' readiness.
        assert!(!r.apply_probe_result("nope", AdapterStatus::Missing { hint: "x".into() }));
        assert_eq!(r.counts().ready, 0);
    }

    #[test]
    fn status_counts_nothing_ready() {
        let mut c = StatusCounts {
            missing: 3,
            ..StatusCounts::default()
        };
        assert!(c.nothing_ready());
        c.ready = 1;
        assert!(!c.nothing_ready());
    }
}

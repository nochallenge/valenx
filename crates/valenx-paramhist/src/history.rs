//! History data model — `Vec<HistEntry>` where each entry lists its
//! upstream dependencies by *index*.
//!
//! Indices reference positions in the [`History::entries`] Vec. The
//! invariant is that `dependencies[k] < self_index` for every entry
//! whose dependencies were assigned at insertion time — this keeps
//! the Vec topologically sorted by default. Insert/move ops below
//! preserve that invariant by validating each change.

use std::collections::{BTreeSet, HashSet, VecDeque};

use serde::{Deserialize, Serialize};

use crate::error::ParamHistError;

/// One history entry — name + dependency list + dirty flag.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HistEntry {
    /// Free-form name (e.g. `"sketch1"`, `"extrude1"`).
    pub name: String,
    /// Indices of upstream entries this one depends on.
    pub dependencies: Vec<usize>,
    /// Has been (re)evaluated since the last upstream change.
    pub up_to_date: bool,
}

impl HistEntry {
    /// Convenience constructor.
    pub fn new(name: impl Into<String>, dependencies: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            dependencies,
            up_to_date: false,
        }
    }
}

/// Result of rebuilding a single entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RebuildResult {
    /// Index of the rebuilt entry.
    pub index: usize,
    /// Name copy for readability.
    pub name: String,
}

/// The history list.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct History {
    /// Entries in insertion order.
    pub entries: Vec<HistEntry>,
}

impl History {
    /// New empty history.
    pub fn new() -> Self {
        Self::default()
    }

    /// Push an entry at the tail (the standard FreeCAD/NaroCAD path).
    pub fn push(&mut self, entry: HistEntry) -> Result<usize, ParamHistError> {
        let idx = self.entries.len();
        for d in &entry.dependencies {
            if *d >= idx {
                return Err(ParamHistError::InvalidMove(format!(
                    "push: dependency {d} >= new index {idx}"
                )));
            }
        }
        self.entries.push(entry);
        Ok(idx)
    }

    /// Insert before position `idx`. Dependencies inside `entry`
    /// must reference indices `< idx`. Downstream dependency
    /// references are bumped by 1 automatically.
    pub fn insert_before(&mut self, idx: usize, entry: HistEntry) -> Result<(), ParamHistError> {
        if idx > self.entries.len() {
            return Err(ParamHistError::IndexOutOfRange {
                idx,
                limit: self.entries.len() + 1,
            });
        }
        for d in &entry.dependencies {
            if *d >= idx {
                return Err(ParamHistError::InvalidMove(format!(
                    "insert_before: dependency {d} >= insertion point {idx}"
                )));
            }
        }
        // Bump references in entries with original positions >= idx.
        for e in self.entries.iter_mut().skip(idx) {
            for d in e.dependencies.iter_mut() {
                if *d >= idx {
                    *d += 1;
                }
            }
        }
        self.entries.insert(idx, entry);
        Ok(())
    }

    /// Move entry from `from` to position `to` (in the new
    /// arrangement). Fails if the move would create a back-edge.
    pub fn move_to(&mut self, from: usize, to: usize) -> Result<(), ParamHistError> {
        if from >= self.entries.len() {
            return Err(ParamHistError::IndexOutOfRange {
                idx: from,
                limit: self.entries.len(),
            });
        }
        if to >= self.entries.len() {
            return Err(ParamHistError::IndexOutOfRange {
                idx: to,
                limit: self.entries.len(),
            });
        }
        if from == to {
            return Ok(());
        }

        // Compute the index permutation: pull `from` out, slot at
        // `to`.
        let mut perm: Vec<usize> = (0..self.entries.len()).collect();
        let v = perm.remove(from);
        perm.insert(to, v);

        // Inverse permutation — new_pos[old_idx] = new_idx.
        let mut new_pos = vec![0usize; self.entries.len()];
        for (new_i, old_i) in perm.iter().enumerate() {
            new_pos[*old_i] = new_i;
        }

        // Check: every dependency must still point to a smaller
        // index after the move.
        for (old_i, e) in self.entries.iter().enumerate() {
            let ni = new_pos[old_i];
            for d in &e.dependencies {
                if new_pos[*d] >= ni {
                    return Err(ParamHistError::InvalidMove(format!(
                        "move would put dep {d} after dependent {old_i}"
                    )));
                }
            }
        }

        // Apply permutation: build new entries vector with
        // remapped dependencies.
        let mut new_entries = Vec::with_capacity(self.entries.len());
        for old_i in &perm {
            let e = &self.entries[*old_i];
            let remapped: Vec<usize> = e.dependencies.iter().map(|d| new_pos[*d]).collect();
            new_entries.push(HistEntry {
                name: e.name.clone(),
                dependencies: remapped,
                up_to_date: e.up_to_date,
            });
        }
        self.entries = new_entries;
        Ok(())
    }
}

/// Kahn's algorithm — produce a topological ordering of the index
/// list. Returns the order if acyclic; otherwise [`ParamHistError::Cycle`].
pub fn topological_sort(entries: &[HistEntry]) -> Result<Vec<usize>, ParamHistError> {
    let n = entries.len();
    // In-degree = number of DISTINCT dependencies. Counting duplicates would
    // never let the relax step below — which decrements once per dependency
    // node (via `contains`) — drive the in-degree to zero, so a node with a
    // repeated dependency would stall and the function would falsely report a
    // cycle for an acyclic graph.
    let mut in_deg = vec![0usize; n];
    for (i, e) in entries.iter().enumerate() {
        in_deg[i] = e
            .dependencies
            .iter()
            .copied()
            .collect::<std::collections::HashSet<usize>>()
            .len();
    }
    let mut q: VecDeque<usize> = in_deg
        .iter()
        .enumerate()
        .filter(|(_, k)| **k == 0)
        .map(|(i, _)| i)
        .collect();
    let mut out = Vec::with_capacity(n);
    while let Some(i) = q.pop_front() {
        out.push(i);
        for (j, e) in entries.iter().enumerate() {
            if e.dependencies.contains(&i) {
                in_deg[j] -= 1;
                if in_deg[j] == 0 {
                    q.push_back(j);
                }
            }
        }
    }
    if out.len() != n {
        for (i, k) in in_deg.iter().enumerate() {
            if *k > 0 {
                return Err(ParamHistError::Cycle(i));
            }
        }
    }
    Ok(out)
}

/// Compute the dirty set — all entries transitively *downstream* of
/// `changed`, plus `changed` itself.
pub fn dirty_set(entries: &[HistEntry], changed: usize) -> HashSet<usize> {
    let mut dirty = HashSet::new();
    dirty.insert(changed);
    let mut grew = true;
    while grew {
        grew = false;
        for (i, e) in entries.iter().enumerate() {
            if !dirty.contains(&i) && e.dependencies.iter().any(|d| dirty.contains(d)) {
                dirty.insert(i);
                grew = true;
            }
        }
    }
    dirty
}

/// Partial rebuild — re-evaluate only the dirty entries, in
/// topological order. Returns the rebuilt list.
pub fn partial_rebuild(
    entries: &mut [HistEntry],
    dirty: &HashSet<usize>,
) -> Result<Vec<RebuildResult>, ParamHistError> {
    // Sub-graph in topological order.
    let order = topological_sort(entries)?;
    let mut out = Vec::new();
    let dirty_set: BTreeSet<usize> = dirty.iter().copied().collect();
    for i in order {
        if dirty_set.contains(&i) {
            entries[i].up_to_date = true;
            out.push(RebuildResult {
                index: i,
                name: entries[i].name.clone(),
            });
        }
    }
    Ok(out)
}

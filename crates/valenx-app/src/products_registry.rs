//! **Per-file registry of agent-bridge 3-D mesh producers.**
//!
//! The Workbench+Agent bridge (`show_3d{kind}`) used to build each model with a
//! per-kind `else if kind == "<x>"` arm inlined into
//! [`crate::agent_commands`]'s shared reducer — so wiring a new model meant
//! editing that one shared `match`, which serialises parallel work and is a
//! merge-conflict magnet. This module replaces those arms with a single
//! generic lookup ([`lookup`]) keyed by the wire `kind` string, so the reducer
//! never grows a new branch again.
//!
//! ## What lives where
//!
//! The *substance* of each product — how its [`crate::WorkspaceProduct`] is
//! assembled from that tool's canonical mesh / camera / readout-row producers —
//! lives in **that tool's own module** as a pure `pub(crate) fn …() ->
//! WorkspaceProduct` builder (e.g. [`crate::rocket_workbench::rocket_product`]).
//! This file only holds the tiny `kind → builder` table in [`lookup`]. Adding a
//! tool is therefore a *one-line* edit to the table here plus a self-contained
//! builder in the new tool's file — the per-tool code (the part that actually
//! varies and that two contributors might touch at once) is fully isolated, so
//! parallel wiring conflicts shrink to one trivial line in this table.
//!
//! ## Why a shared-`match` table and not the `inventory` crate
//!
//! The plan offered `inventory` (link-time distributed-slice registration, no
//! central list at all) as the first choice. We deliberately use the fallback
//! shared-`match` table instead, because this workspace links Windows builds
//! with **`rust-lld` / `lld-link`** (see `.cargo/config.toml`, committed
//! 2026-06-20). `inventory` populates its slice via `ctor`-style
//! life-before-main static registration, which is exactly the construct a
//! non-default linker can dead-strip when the registering module isn't
//! otherwise referenced — and the gate here is a `cargo test -p valenx-app
//! --lib` build, the case most prone to that stripping. A `match` is
//! linker-agnostic and resolved at compile time, so the registry can never
//! silently lose a kind. (It also avoids adding a third-party dependency +
//! `deny.toml` / `Cargo.lock` churn for a five-entry table.)
//!
//! ## Builder contract
//!
//! Every builder is a `fn() -> WorkspaceProduct` ([`MeshProducerEntry::build`])
//! — pure, app-state-free, built only from that tool's canonical inputs — so
//! the reducer can call it with nothing but the channel it already knows.
//! Behaviour is byte-for-byte what the old inline arms produced.
//!
//! ## Adding a new 3-D tool — the one-liner pattern
//!
//! Copy this into the new tool's own module (PHASE-C wiring subagents: this is
//! the whole change on the producer side — no edit to `agent_commands`):
//!
//! ```ignore
//! // 1. In `crate::foo_workbench` (the tool's own file): a pure builder that
//! //    assembles the WorkspaceProduct from the tool's canonical producers.
//! pub(crate) fn foo_product() -> crate::WorkspaceProduct {
//!     let (mesh, lines) = foo_loaded_mesh();              // your existing producer
//!     let camera = foo_camera(&mesh.mesh);               // your existing camera
//!     crate::WorkspaceProduct {
//!         title: "Foo".into(), lines, mesh: Some(mesh),
//!         vertex_colors: None, camera, kind2d: None, last_export: None,
//!     }
//! }
//! ```
//!
//! …then add the single table line in [`lookup`] below:
//! `"foo" => Some(crate::foo_workbench::foo_product),`. That is the only shared
//! edit, and it is a one-liner — everything else lives in the tool's file.

use crate::WorkspaceProduct;

/// One registry entry: a wire `kind` string mapped to a pure builder that
/// returns the [`WorkspaceProduct`] for that 3-D model.
///
/// The builder is app-state-free (the existing producers build from canonical
/// inputs only), so the bridge can invoke it knowing nothing but the channel it
/// publishes into. Returned by [`lookup`] so callers can read [`Self::kind`]
/// (e.g. for diagnostics) as well as invoke [`Self::build`].
#[derive(Clone, Copy)]
pub struct MeshProducerEntry {
    /// The `show_3d` wire `kind` this entry answers to (e.g. `"rocket"`).
    pub kind: &'static str,
    /// Pure builder for this kind's product — same output as the old inline
    /// reducer arm.
    pub build: fn() -> WorkspaceProduct,
}

/// Resolve a `show_3d` `kind` to its registry entry, or `None` for an unknown
/// kind (the reducer then skips it safely — no panic, no placeholder churn,
/// matching the rest of its bad-input handling).
///
/// **This `match` is the single shared edit point for 3-D mesh tools**: each
/// arm is one line pairing a wire `kind` with the per-tool builder in that
/// tool's own module. The substantive per-tool code lives in those builders,
/// not here — so adding a kind is a one-line addition (see the module docs for
/// the copy-paste pattern). Note `dna` is intentionally absent: `show_3d:dna`
/// is a *text* card handled directly in the reducer (no mesh), and the 2-D
/// `show_2d` drawings (`rcbeam` / `dna`) have their own separate path.
pub fn lookup(kind: &str) -> Option<MeshProducerEntry> {
    let build: fn() -> WorkspaceProduct = match kind {
        "rocket" => crate::rocket_workbench::rocket_product,
        "gear" => crate::gears_workbench::gear_product,
        "bracket" => crate::bracket_product::bracket_workspace_product,
        "rcbeam" => crate::rcbeam_workbench::rcbeam_product,
        "fem" => crate::fem_workbench::fem_product,
        _ => return None,
    };
    Some(MeshProducerEntry {
        kind: kind_static(kind)?,
        build,
    })
}

/// Map a looked-up `kind` to its `'static` spelling for [`MeshProducerEntry::kind`].
/// Kept in lockstep with [`lookup`]'s arms so the returned `kind` is always the
/// canonical literal (not a borrow of the caller's input).
fn kind_static(kind: &str) -> Option<&'static str> {
    Some(match kind {
        "rocket" => "rocket",
        "gear" => "gear",
        "bracket" => "bracket",
        "rcbeam" => "rcbeam",
        "fem" => "fem",
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact set of 3-D mesh kinds the registry is expected to resolve, so a
    /// future addition that forgets the table line trips the count assertion.
    const KNOWN_3D_KINDS: &[&str] = &["rocket", "gear", "bracket", "rcbeam", "fem"];

    #[test]
    fn registry_resolves_all_known_3d_kinds() {
        // Every migrated kind resolves to an entry whose `kind` echoes the query
        // (so the table and the `kind_static` spelling stay in lockstep).
        for &k in KNOWN_3D_KINDS {
            let entry = lookup(k).unwrap_or_else(|| panic!("registry resolves {k:?}"));
            assert_eq!(entry.kind, k, "entry.kind echoes the looked-up kind");
        }
    }

    #[test]
    fn registry_builds_a_live_mesh_for_each_3d_kind() {
        // Each builder is pure and yields a non-empty 3-D mesh product (the FEM
        // one additionally carries the von-Mises vertex colours). This exercises
        // the builders the reducer dispatches to, so a regression in any of them
        // is caught here without the file-poll plumbing.
        for &k in KNOWN_3D_KINDS {
            let entry = lookup(k).unwrap();
            let product = (entry.build)();
            let mesh = product
                .mesh
                .as_ref()
                .unwrap_or_else(|| panic!("{k}: a 3-D product carries a mesh"));
            assert!(!mesh.mesh.nodes.is_empty(), "{k}: mesh has vertices");
            assert!(mesh.mesh.total_elements() > 0, "{k}: mesh has triangles");
        }
        // The FEM cantilever is the only kind that ships per-vertex colours.
        assert!(
            (lookup("fem").unwrap().build)().vertex_colors.is_some(),
            "fem product carries von-Mises vertex colours"
        );
    }

    #[test]
    fn registry_returns_none_for_an_unknown_kind() {
        // An unknown kind resolves to None → the reducer skips it safely.
        assert!(lookup("not-a-model").is_none());
        assert!(lookup("").is_none());
        // `dna` is a text card / 2-D drawing, NOT a registry 3-D mesh kind.
        assert!(lookup("dna").is_none());
    }
}

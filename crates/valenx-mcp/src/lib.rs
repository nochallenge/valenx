//! # valenx-mcp
//!
//! MCP (Model Context Protocol) server that exposes Valenx
//! capabilities to LLM clients. Speaks JSON-RPC over stdio.
//!
//! Currently advertised tools:
//! - `dock` — run a native docking job, stream events
//! - `dry_run` — pre-flight a docking job
//! - `describe` — return JSON Schema for inputs and outputs
//! - `list_adapters` — return the full adapter catalog
//! - the **generative-design** tool set ([`design`]) — `create_sketch`,
//!   `add_sketch_line`, `add_sketch_circle`, `add_constraint`, `pad`,
//!   `pocket`, `revolve`, `fillet`, `boolean`, `evaluate_design`,
//!   `export_design`, `reset_design` — so an external LLM can compose,
//!   measure, and iterate a parametric CAD part. Valenx ships no ML
//!   model; the LLM driving these tools is the generative part.
//!
//! This is a deliberately small surface — adapters land tool-by-tool
//! as their native crates ship (oxDNA, RNA folding, FoldSeek, etc.).
//!
//! # Path sandboxing
//!
//! Every `receptor_path` / `ligand_path` / `output_path` arriving via
//! the protocol is rejected unless it resolves under
//! `$VALENX_MCP_SANDBOX_DIR` (or `<tempdir>/valenx-mcp` when that env
//! var is unset). The check is implemented in [`sandbox::sandbox_check`];
//! callers should treat the returned absolute path as the canonical
//! form and pass nothing else to file IO. This prevents a misbehaving
//! MCP client from reading or overwriting arbitrary files the server
//! has access to.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod design;
pub mod sandbox;
pub mod server;
pub mod tools;

pub use sandbox::{sandbox_check, sandbox_root, SANDBOX_ENV};
pub use server::serve_stdio;

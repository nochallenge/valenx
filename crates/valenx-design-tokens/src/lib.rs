//! # valenx-design-tokens
//!
//! Single source of truth for Valenx design tokens.
//!
//! The canonical definitions live in **`tokens.json`** at this crate's
//! root. A build script reads that file, validates it, and emits
//! Rust `const`s — the full API below is generated, not authored.
//!
//! Use these consts everywhere in UI code:
//!
//! ```ignore
//! use valenx_design_tokens::color::{surface, text, accent};
//! use valenx_design_tokens::space;
//!
//! let bg        = surface::S1;        // canonical surface-1 colour
//! let fg        = text::T1;           // canonical text-1 colour
//! let highlight = accent::PRIMARY;    // accent by role
//! let padding   = space::S4;          // 16 px gutter
//! ```
//!
//! See [RFC 0006](../rfcs/0006-token-system.md) for the full design
//! and pipeline spec.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

include!(concat!(env!("OUT_DIR"), "/generated.rs"));

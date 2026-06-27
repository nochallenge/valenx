//! # valenx-meshgen
//!
//! Bring-your-own mesh-LLM pipeline. valenx ships this runtime + codec; the
//! user supplies the model weights. See
//! `docs/design/2026-06-22-meshgen-byo-llm-design.md`.
#![forbid(unsafe_code)]

pub mod backend;
pub mod codec;
pub mod error;
pub mod model;

pub use backend::{MeshLlm, MockBackend};
pub use error::MeshGenError;
pub use model::{Device, GenParams, ModelConfig, ModelProfile};

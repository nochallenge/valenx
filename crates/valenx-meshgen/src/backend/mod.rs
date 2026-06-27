//! The MeshLlm runtime trait + a weightless MockBackend for tests.

use std::ops::ControlFlow;

use crate::error::MeshGenError;
use crate::model::{Device, GenParams, ModelConfig};

/// A loadable mesh-generation language model. Backends are swappable; callers
/// never know which concrete type they hold.
pub trait MeshLlm: Send {
    /// Load a model described by `config` onto `device`.
    fn load(config: &ModelConfig, device: Device) -> Result<Self, MeshGenError>
    where
        Self: Sized;

    /// Stream a completion for `prompt`. `on_token` is called once per decoded
    /// token chunk; returning `ControlFlow::Break` cancels promptly.
    fn generate(
        &mut self,
        prompt: &str,
        params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> ControlFlow<()>,
    ) -> Result<(), MeshGenError>;
}

/// A weightless backend that replays a fixed token list. For tests + for the
/// "no model installed" demo path. Does not implement `MeshLlm::load` in a
/// meaningful way (there is nothing to load); construct it with `new`.
pub struct MockBackend {
    tokens: Vec<String>,
}

impl MockBackend {
    /// Build a mock that will replay `tokens` in order.
    pub fn new<I, S>(tokens: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self {
            tokens: tokens.into_iter().map(Into::into).collect(),
        }
    }

    /// Replay the canned tokens through `on_token`, honoring cancellation.
    pub fn generate(
        &mut self,
        _prompt: &str,
        _params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> ControlFlow<()>,
    ) -> Result<(), MeshGenError> {
        for tok in &self.tokens {
            if on_token(tok).is_break() {
                break;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codec::decode;
    use crate::model::ModelProfile;

    #[test]
    fn mock_backend_streams_tokens_that_decode_to_a_mesh() {
        // The "model" emits a tetrahedron, split across arbitrary token chunks.
        let tokens = vec![
            "v 0 0 0\nv 127 0 0\n",
            "v 0 127 0\nv 0 0 127\n",
            "f 1 2 3\nf 1 2 4\n",
            "f 1 3 4\nf 2 3 4\n",
        ];
        let mut backend = MockBackend::new(tokens);

        let mut acc = String::new();
        backend
            .generate("a tetrahedron", &GenParams::default(), &mut |tok| {
                acc.push_str(tok);
                ControlFlow::Continue(())
            })
            .unwrap();

        let mesh = decode(&acc, &ModelProfile::LLAMA_MESH);
        assert_eq!(mesh.nodes.len(), 4);
    }

    #[test]
    fn mock_backend_honors_cancellation() {
        let tokens = vec!["a", "b", "c", "d", "e"];
        let mut backend = MockBackend::new(tokens);
        let mut seen = 0;
        backend
            .generate("x", &GenParams::default(), &mut |_tok| {
                seen += 1;
                if seen == 2 {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            })
            .unwrap();
        assert_eq!(seen, 2, "generation stops at the first Break");
    }
}

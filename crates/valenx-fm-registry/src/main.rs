//! Demo: list the gated foundation-model registry and probe each one.
//!
//! In any environment without staged weights + a GPU runtime (CI, this repo)
//! every model reports `BLOCKED` — the registry documents the integration point
//! honestly and never fabricates a prediction.
//!
//! Run with: `cargo run -p valenx-fm-registry`

use valenx_fm_registry::{probe, registry};

fn main() {
    println!("=== valenx-fm-registry: gated foundation models ===");
    println!("(this tool reports readiness only — it never runs inference)\n");

    let mut blocked = 0usize;
    for m in registry() {
        let status = probe(&m);
        if status.is_blocked() {
            blocked += 1;
        }
        println!("• {} [{}]", m.display_name, m.task.as_str());
        println!("    license : {}", m.license);
        println!("    upstream: {}", m.upstream_url);
        println!("    weights : ${}", m.weights_env);
        println!("    status  : {}", status.message(&m));
        println!();
    }

    let total = registry().len();
    println!("{blocked}/{total} models BLOCKED in this environment (no staged weights / GPU).");
    println!(
        "Stage weights and set the matching $VALENX_*_WEIGHTS path to make a model probe ready; \
         actually running it remains the upstream tool's job."
    );
}

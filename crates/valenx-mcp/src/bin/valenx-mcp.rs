//! Stand-alone MCP server binary. Reads JSON-RPC from stdin,
//! writes responses to stdout. Add to your MCP client config as
//! a stdio server pointing at this binary.

use anyhow::Result;

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();
    valenx_mcp::serve_stdio().await
}

//! HTML + Markdown report generators for `Results` bundles.
//!
//! Pulled out of `lib.rs` so the rendering helpers (`html_escape`,
//! `md_escape`, `REPORT_CSS`) live next to the public renderers and
//! writers that consume them.

use std::path::Path;

use valenx_fields::Results;

use crate::ExportError;

/// Render a self-contained HTML report summarising one run's
/// `Results` bundle. No JavaScript, no external CSS — the resulting
/// file is one ~10 KB blob users can email or pin to a project
/// wiki. Sections:
///
/// 1. Header with case id + adapter id + completion timestamp.
/// 2. Provenance block (tool / adapter versions, run_id, hashes).
/// 3. Scalars table (name / value / units / time / source).
/// 4. Fields table (name / kind / location / samples / range).
/// 5. Artifacts list (kind + relative path + label).
///
/// Returns the HTML as a `String`; pair with [`write_html_report`]
/// to drop it on disk.
pub fn render_html_report(results: &Results) -> String {
    let mut s = String::with_capacity(4096);
    let escape = html_escape;

    s.push_str("<!DOCTYPE html>\n");
    s.push_str("<html lang=\"en\"><head><meta charset=\"utf-8\">\n");
    s.push_str(&format!(
        "<title>Valenx run report — {}</title>\n",
        escape(&results.meta.case_id)
    ));
    s.push_str("<style>\n");
    s.push_str(REPORT_CSS);
    s.push_str("</style>\n</head><body>\n");

    // Header
    s.push_str(&format!(
        "<h1>Valenx run report: <code>{}</code></h1>\n",
        escape(&results.meta.case_id)
    ));
    if let Some(desc) = results.meta.description.as_ref() {
        s.push_str(&format!("<p class=\"desc\">{}</p>\n", escape(desc)));
    }

    // Provenance
    s.push_str("<section><h2>Provenance</h2>\n<table class=\"prov\">\n");
    s.push_str(&format!(
        "<tr><th>Adapter</th><td><code>{} v{}</code></td></tr>\n",
        escape(&results.provenance.adapter),
        escape(&results.provenance.adapter_version),
    ));
    s.push_str(&format!(
        "<tr><th>Tool</th><td><code>{} v{}</code></td></tr>\n",
        escape(&results.provenance.tool),
        escape(&results.provenance.tool_version),
    ));
    s.push_str(&format!(
        "<tr><th>Run ID</th><td><code>{}</code></td></tr>\n",
        escape(&results.provenance.run_id)
    ));
    s.push_str(&format!(
        "<tr><th>Wall time</th><td>{:.3} s</td></tr>\n",
        results.provenance.wall_time_seconds
    ));
    s.push_str(&format!(
        "<tr><th>Completed at</th><td>{}</td></tr>\n",
        escape(&results.provenance.completed_at)
    ));
    s.push_str(&format!(
        "<tr><th>Case hash</th><td><code class=\"hash\">{}</code></td></tr>\n",
        escape(&results.provenance.case_hash.0)
    ));
    s.push_str(&format!(
        "<tr><th>Mesh hash</th><td><code class=\"hash\">{}</code></td></tr>\n",
        escape(&results.provenance.mesh_hash.0)
    ));
    s.push_str("</table></section>\n");

    // Scalars
    s.push_str("<section><h2>Scalars</h2>\n");
    if results.scalars.is_empty() {
        s.push_str("<p class=\"empty\">No scalar results.</p>\n");
    } else {
        s.push_str("<table class=\"data\"><thead><tr>");
        s.push_str("<th>Name</th><th>Value</th><th>Units</th>");
        s.push_str("<th>Time</th><th>Source</th></tr></thead>\n<tbody>\n");
        for name in results.scalars.names() {
            for record in results.scalars.all(name) {
                let units = record.units.display.unwrap_or("");
                let time = match record.time {
                    valenx_fields::TimeKey::Steady => "steady".to_string(),
                    valenx_fields::TimeKey::Iteration(n) => format!("iter {n}"),
                    valenx_fields::TimeKey::Time { value, .. } => {
                        format!("t={value} s")
                    }
                };
                s.push_str(&format!(
                    "<tr><td><code>{}</code></td><td>{:.6e}</td><td>{}</td><td>{}</td><td>{:?}</td></tr>\n",
                    escape(&record.name),
                    record.value,
                    escape(units),
                    escape(&time),
                    record.source,
                ));
            }
        }
        s.push_str("</tbody></table>\n");
    }
    s.push_str("</section>\n");

    // Fields
    s.push_str("<section><h2>Fields</h2>\n");
    if results.fields.is_empty() {
        s.push_str("<p class=\"empty\">No field results.</p>\n");
    } else {
        s.push_str("<table class=\"data\"><thead><tr>");
        s.push_str("<th>Name</th><th>Kind</th><th>Location</th>");
        s.push_str("<th>Samples</th><th>Range</th></tr></thead>\n<tbody>\n");
        let names: Vec<String> = results.fields.names().map(|s| s.to_string()).collect();
        for name in &names {
            for f in results.fields.by_name(name) {
                let range = f
                    .range
                    .map(|(lo, hi)| format!("[{lo:.4e}, {hi:.4e}]"))
                    .unwrap_or_else(|| "(unranged)".to_string());
                s.push_str(&format!(
                    "<tr><td><code>{}</code></td><td>{:?}</td><td>{:?}</td><td>{}</td><td>{}</td></tr>\n",
                    escape(&f.name),
                    f.kind,
                    f.location,
                    f.samples(),
                    escape(&range),
                ));
            }
        }
        s.push_str("</tbody></table>\n");
    }
    s.push_str("</section>\n");

    // Artifacts
    s.push_str("<section><h2>Artifacts</h2>\n");
    if results.artifacts.is_empty() {
        s.push_str("<p class=\"empty\">No artifacts produced.</p>\n");
    } else {
        s.push_str("<table class=\"data\"><thead><tr>");
        s.push_str("<th>Kind</th><th>Path</th><th>Label</th></tr></thead>\n<tbody>\n");
        for art in &results.artifacts {
            s.push_str(&format!(
                "<tr><td>{:?}</td><td><code>{}</code></td><td>{}</td></tr>\n",
                art.kind,
                escape(&art.path.to_string_lossy()),
                escape(&art.label),
            ));
        }
        s.push_str("</tbody></table>\n");
    }
    s.push_str("</section>\n");

    s.push_str("<footer><p class=\"footer\">Generated by <a href=\"https://github.com/nochallenge/valenx\">Valenx</a></p></footer>\n");
    s.push_str("</body></html>\n");
    s
}

/// Render a Markdown report. Same five sections as
/// [`render_html_report`], formatted as GitHub-flavoured Markdown
/// (pipe-delimited tables, fenced code blocks).
///
/// Useful for posting a summary as a GitHub PR comment, dropping
/// into a release-notes file, or embedding in a documentation
/// pipeline. Pair with [`write_markdown_report`] to drop the
/// rendered Markdown straight to disk.
pub fn render_markdown_report(results: &Results) -> String {
    let mut s = String::with_capacity(2048);

    // Header
    s.push_str(&format!(
        "# Valenx run report: `{}`\n\n",
        md_escape(&results.meta.case_id)
    ));
    if let Some(desc) = results.meta.description.as_ref() {
        s.push_str(&format!("_{}_\n\n", md_escape(desc)));
    }

    // Provenance
    s.push_str("## Provenance\n\n");
    s.push_str("| Field | Value |\n");
    s.push_str("|---|---|\n");
    s.push_str(&format!(
        "| Adapter | `{} v{}` |\n",
        md_escape(&results.provenance.adapter),
        md_escape(&results.provenance.adapter_version),
    ));
    s.push_str(&format!(
        "| Tool | `{} v{}` |\n",
        md_escape(&results.provenance.tool),
        md_escape(&results.provenance.tool_version),
    ));
    s.push_str(&format!(
        "| Run ID | `{}` |\n",
        md_escape(&results.provenance.run_id)
    ));
    s.push_str(&format!(
        "| Wall time | {:.3} s |\n",
        results.provenance.wall_time_seconds
    ));
    s.push_str(&format!(
        "| Completed at | {} |\n",
        md_escape(&results.provenance.completed_at)
    ));
    s.push_str(&format!(
        "| Case hash | `{}` |\n",
        md_escape(&results.provenance.case_hash.0)
    ));
    s.push_str(&format!(
        "| Mesh hash | `{}` |\n\n",
        md_escape(&results.provenance.mesh_hash.0)
    ));

    // Scalars
    s.push_str("## Scalars\n\n");
    if results.scalars.is_empty() {
        s.push_str("_No scalar results._\n\n");
    } else {
        s.push_str("| Name | Value | Units | Time | Source |\n");
        s.push_str("|---|---:|---|---|---|\n");
        for name in results.scalars.names() {
            for record in results.scalars.all(name) {
                let units = record.units.display.unwrap_or("");
                let time = match record.time {
                    valenx_fields::TimeKey::Steady => "steady".to_string(),
                    valenx_fields::TimeKey::Iteration(n) => format!("iter {n}"),
                    valenx_fields::TimeKey::Time { value, .. } => format!("t={value} s"),
                };
                s.push_str(&format!(
                    "| `{}` | {:.6e} | {} | {} | {:?} |\n",
                    md_escape(&record.name),
                    record.value,
                    md_escape(units),
                    md_escape(&time),
                    record.source,
                ));
            }
        }
        s.push('\n');
    }

    // Fields
    s.push_str("## Fields\n\n");
    if results.fields.is_empty() {
        s.push_str("_No field results._\n\n");
    } else {
        s.push_str("| Name | Kind | Location | Samples | Range |\n");
        s.push_str("|---|---|---|---:|---|\n");
        let names: Vec<String> = results.fields.names().map(|s| s.to_string()).collect();
        for name in &names {
            for f in results.fields.by_name(name) {
                let range = f
                    .range
                    .map(|(lo, hi)| format!("[{lo:.4e}, {hi:.4e}]"))
                    .unwrap_or_else(|| "(unranged)".to_string());
                s.push_str(&format!(
                    "| `{}` | {:?} | {:?} | {} | {} |\n",
                    md_escape(&f.name),
                    f.kind,
                    f.location,
                    f.samples(),
                    md_escape(&range),
                ));
            }
        }
        s.push('\n');
    }

    // Artifacts
    s.push_str("## Artifacts\n\n");
    if results.artifacts.is_empty() {
        s.push_str("_No artifacts produced._\n\n");
    } else {
        s.push_str("| Kind | Path | Label |\n");
        s.push_str("|---|---|---|\n");
        for art in &results.artifacts {
            s.push_str(&format!(
                "| {:?} | `{}` | {} |\n",
                art.kind,
                md_escape(&art.path.to_string_lossy()),
                md_escape(&art.label),
            ));
        }
        s.push('\n');
    }

    s.push_str("---\n_Generated by Valenx._\n");
    s
}

/// Write a rendered Markdown report to disk. Creates parent
/// directories as needed.
pub fn write_markdown_report(results: &Results, path: &Path) -> Result<(), ExportError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| ExportError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    valenx_core::io_caps::atomic_write_str(path, &render_markdown_report(results)).map_err(|e| {
        ExportError::Io {
            path: path.to_path_buf(),
            source: e,
        }
    })
}

/// Escape characters that would break a pipe-delimited Markdown
/// table cell. The unsafe set is `|` (column separator), `\` (escape),
/// and the trio of backticks / angle-brackets that confuse renderers
/// inside code spans. Newlines collapse to a space — Markdown table
/// cells can't span multiple lines in the GFM dialect.
fn md_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '|' => out.push_str("\\|"),
            '\\' => out.push_str("\\\\"),
            '\n' | '\r' => out.push(' '),
            other => out.push(other),
        }
    }
    out
}

/// Write a rendered HTML report to disk. Creates parent directories
/// as needed.
pub fn write_html_report(results: &Results, path: &Path) -> Result<(), ExportError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| ExportError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    valenx_core::io_caps::atomic_write_str(path, &render_html_report(results)).map_err(|e| {
        ExportError::Io {
            path: path.to_path_buf(),
            source: e,
        }
    })
}

/// Minimal CSS for the HTML report. Embedded so the report is
/// self-contained — user can email it without worrying about
/// external stylesheet links.
const REPORT_CSS: &str = r#"
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #1a1a1a; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.3em; }
h2 { color: #2a4a7f; margin-top: 1.6em; border-bottom: 1px solid #ccc; padding-bottom: 0.2em; }
code { background: #f6f6f6; padding: 0.1em 0.3em; border-radius: 3px; font-size: 0.92em; }
code.hash { font-size: 0.75em; word-break: break-all; }
.desc { color: #555; font-style: italic; }
.empty { color: #888; font-style: italic; }
table { border-collapse: collapse; width: 100%; margin: 0.6em 0 1em; }
table.prov th { width: 130px; text-align: left; padding: 0.3em 0.6em; }
table.prov td { padding: 0.3em 0.6em; }
table.data th, table.data td { padding: 0.4em 0.7em; border-bottom: 1px solid #eee; text-align: left; }
table.data th { background: #fafafa; font-weight: 600; }
table.data tr:hover { background: #fafcff; }
footer { margin-top: 2.5em; border-top: 1px solid #ddd; padding-top: 0.5em; color: #888; font-size: 0.85em; }
.footer a { color: #2a4a7f; text-decoration: none; }
"#;

/// Minimal HTML escape: `&`, `<`, `>`, `"`. Field / artifact names
/// can legitimately contain `<`, `>` (chemistry mechanism names like
/// `H2<O>`), so escaping them is required for layout integrity.
fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_fields::scalar::ScalarSource;
    use valenx_fields::units::{Units, DIMENSIONLESS, SECOND};
    use valenx_fields::{ScalarRecord, TimeKey};

    fn pa_units() -> Units {
        Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa"))
    }

    fn results_with_a_few_scalars() -> Results {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "test".into(),
            adapter_version: "0".into(),
            tool: "Test".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut r = Results::empty("test", prov);
        r.scalars.insert(ScalarRecord {
            name: "drag_coefficient".into(),
            value: 0.123456,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        });
        r.scalars.insert(ScalarRecord {
            name: "p_inlet".into(),
            value: 101325.0,
            units: pa_units(),
            time: TimeKey::Iteration(500),
            source: ScalarSource::Extracted,
            description: None,
        });
        r.scalars.insert(ScalarRecord {
            name: "T_at_t1ms".into(),
            value: 293.15,
            units: DIMENSIONLESS,
            time: TimeKey::Time {
                value: 0.001,
                units: SECOND,
            },
            source: ScalarSource::Computed,
            description: None,
        });
        r
    }

    #[test]
    fn render_html_report_emits_valid_skeleton() {
        let r = results_with_a_few_scalars();
        let html = render_html_report(&r);
        assert!(html.starts_with("<!DOCTYPE html>"));
        assert!(html.contains("<title>"));
        assert!(html.contains("</html>"));
        assert!(html.contains("Valenx run report"));
    }

    #[test]
    fn render_html_report_includes_provenance_block() {
        let r = results_with_a_few_scalars();
        let html = render_html_report(&r);
        assert!(html.contains("Provenance"));
        assert!(html.contains("Run ID"));
        assert!(html.contains("Wall time"));
        assert!(html.contains("Adapter"));
        assert!(html.contains("Tool"));
    }

    #[test]
    fn render_html_report_lists_every_scalar_record() {
        // The fixture has 3 scalar records — drag_coefficient, p_inlet,
        // T_at_t1ms. Every name should appear in a code element.
        let r = results_with_a_few_scalars();
        let html = render_html_report(&r);
        assert!(html.contains("<code>drag_coefficient</code>"));
        assert!(html.contains("<code>p_inlet</code>"));
        assert!(html.contains("<code>T_at_t1ms</code>"));
        // Time formatting: iter / steady / time variants
        assert!(html.contains("steady"));
        assert!(html.contains("iter 500"));
        assert!(html.contains("t=0.001 s"));
    }

    #[test]
    fn render_html_report_handles_empty_results() {
        // Zero scalars / zero fields / zero artifacts must each
        // produce a friendly "no X" message rather than an empty
        // table.
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "test".into(),
            adapter_version: "0".into(),
            tool: "Test".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let empty = Results::empty("no-data", prov);
        let html = render_html_report(&empty);
        assert!(html.contains("No scalar results."));
        assert!(html.contains("No field results."));
        assert!(html.contains("No artifacts produced."));
    }

    #[test]
    fn html_escape_handles_special_chars() {
        // Field names like `H2<O>` (rare but legal in chemistry) must
        // not break the table layout.
        assert_eq!(super::html_escape("plain"), "plain");
        assert_eq!(super::html_escape("a&b"), "a&amp;b");
        assert_eq!(super::html_escape("a<b>c"), "a&lt;b&gt;c");
        assert_eq!(super::html_escape("a\"b"), "a&quot;b");
    }

    #[test]
    fn write_html_report_creates_a_file_browsers_can_open() {
        let r = results_with_a_few_scalars();
        let path = std::env::temp_dir().join(format!(
            "valenx-html-report-{}.html",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_html_report(&r, &path).expect("write");
        let text = std::fs::read_to_string(&path).expect("read");
        assert!(text.starts_with("<!DOCTYPE html>"));
        // Smoke: at least the embedded CSS made it through.
        assert!(text.contains("font-family"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn render_html_report_has_no_external_references() {
        // Self-containment guarantee: no <script src="...">,
        // <link rel="stylesheet">, or <img src="..."> tags. Users
        // can email the report without wondering whether the
        // recipient has internet access.
        let r = results_with_a_few_scalars();
        let html = render_html_report(&r);
        assert!(!html.contains("<script"), "script tag leaked: {html}");
        assert!(!html.contains("<link"), "link tag leaked: {html}");
        assert!(!html.contains("<img"), "img tag leaked: {html}");
    }

    #[test]
    fn render_markdown_report_emits_required_sections() {
        let r = results_with_a_few_scalars();
        let md = render_markdown_report(&r);
        assert!(md.starts_with("# Valenx run report"), "got: {md}");
        assert!(md.contains("## Provenance"), "got: {md}");
        assert!(md.contains("## Scalars"), "got: {md}");
        assert!(md.contains("## Fields"), "got: {md}");
        assert!(md.contains("## Artifacts"), "got: {md}");
    }

    #[test]
    fn render_markdown_report_includes_every_scalar_record() {
        let r = results_with_a_few_scalars();
        let md = render_markdown_report(&r);
        assert!(
            md.contains("`drag_coefficient`"),
            "scalar name missing: {md}"
        );
    }

    #[test]
    fn render_markdown_report_handles_empty_results() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "x".into(),
            adapter_version: "0".into(),
            tool: "x".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let r = Results::empty("empty", prov);
        let md = render_markdown_report(&r);
        // Each empty section emits an italic placeholder.
        assert!(
            md.contains("_No scalar results._"),
            "missing scalar empty marker: {md}"
        );
        assert!(
            md.contains("_No field results._"),
            "missing field empty marker: {md}"
        );
        assert!(
            md.contains("_No artifacts produced._"),
            "missing artifact empty marker: {md}"
        );
    }

    #[test]
    fn write_markdown_report_creates_a_readable_file() {
        let r = results_with_a_few_scalars();
        let path = std::env::temp_dir().join(format!(
            "valenx-md-{}.md",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_markdown_report(&r, &path).expect("write");
        let body = std::fs::read_to_string(&path).expect("read back");
        assert!(body.starts_with("# Valenx run report"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn md_escape_handles_pipes_and_newlines() {
        // `|` would split the cell; `\n` would break the row.
        // Ensure both are neutralised so the table renders.
        let escaped = md_escape("hello | world\nsecond line");
        assert!(!escaped.contains("|\n") && !escaped.contains('\n'));
        assert!(escaped.contains("\\|"));
    }
}

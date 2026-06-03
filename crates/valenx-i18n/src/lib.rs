//! # valenx-i18n
//!
//! Lightweight string catalogue for Valenx's user-visible text.
//!
//! ## What this is
//!
//! A `LocaleCatalogue` wraps a flat `BTreeMap<String, String>` of
//! `key = value` pairs loaded from `.ftl`-shaped files. Lookups
//! return the matched string, or the key itself (with a tracing
//! warning) when the catalogue doesn't have an entry. A simple
//! `{ $name }` placeholder substitution covers the dynamic-value
//! case — see [`LocaleCatalogue::format_with`].
//!
//! ## What this isn't (yet)
//!
//! Not a full fluent-rs implementation. Plurals, gender, complex
//! ICU MessageFormat, and bidi handling are all out of scope for
//! v0.1.0. The catalogue API is intentionally narrow so a v0.2.0
//! swap to fluent-rs is a one-crate change — every call site
//! routes through [`LocaleCatalogue::lookup`] / [`LocaleCatalogue::format_with`], never
//! a macro that bakes in our internal storage.
//!
//! ## File format
//!
//! ```text
//! # comment line — ignored
//! ribbon.run = Run
//! dialog.about.title = About Valenx
//! error.tool-not-installed = Tool `{ $tool }` is not installed.
//! ```
//!
//! Keys: anything before the first `=`, trimmed.
//! Values: anything after, trimmed. Multi-line values are not
//! supported — split into multiple keys.
//!
//! Lines starting with `#` and blank lines are skipped. Lines
//! without an `=` produce a tracing warning + are skipped.

#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

/// Round-18 L2 cap on per-`.ftl` reads. Realistic locale files top
/// out at a few hundred KiB even for fully translated apps; 4 MiB
/// is generous while rejecting a hostile multi-GB file dropped into
/// the locale dir before the line scanner allocates the merged
/// catalogue.
pub const MAX_FTL_FILE_BYTES: u64 = 4 * 1024 * 1024;

/// Round-19 L1 cap on the number of `.ftl` files we'll walk under a
/// single locale dir. A realistic locale dir has 1–20 files (one
/// per top-level UI module). 1000 is well past anything sane; the
/// cap exists so a poisoned locale dir with millions of placeholder
/// `.ftl` files can't allocate an unbounded `Vec<_>` of read-dir
/// entries before any of them get parsed.
pub const MAX_LOCALE_FILES: usize = 1_000;

/// Read `path` to a `String` with [`MAX_FTL_FILE_BYTES`] enforced
/// both by stat AND by a bounded `take()` on the read — the same
/// belt-and-braces pattern the rest of the workspace uses to close
/// the TOCTOU window between metadata and open.
fn read_capped_ftl(path: &Path) -> std::io::Result<String> {
    let meta = std::fs::metadata(path)?;
    if meta.len() > MAX_FTL_FILE_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "ftl file {} exceeds {}-byte cap (actual: {})",
                path.display(),
                MAX_FTL_FILE_BYTES,
                meta.len()
            ),
        ));
    }
    let mut buf = Vec::new();
    std::fs::File::open(path)?
        .take(MAX_FTL_FILE_BYTES)
        .read_to_end(&mut buf)?;
    String::from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// A loaded set of key→string mappings for one locale.
#[derive(Clone, Debug, Default)]
pub struct LocaleCatalogue {
    locale: String,
    entries: BTreeMap<String, String>,
}

impl LocaleCatalogue {
    /// Build a catalogue from in-memory `.ftl` text. Returns an
    /// empty catalogue when the input is empty.
    pub fn from_str(locale: impl Into<String>, text: &str) -> Self {
        let mut entries: BTreeMap<String, String> = BTreeMap::new();
        for (line_no, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let Some((key, value)) = line.split_once('=') else {
                tracing::warn!(
                    target: "valenx-i18n",
                    line = line_no + 1,
                    text = line,
                    "skipping line without `=`"
                );
                continue;
            };
            let key = key.trim().to_string();
            let value = value.trim().to_string();
            if key.is_empty() {
                tracing::warn!(
                    target: "valenx-i18n",
                    line = line_no + 1,
                    "skipping line with empty key"
                );
                continue;
            }
            entries.insert(key, value);
        }
        Self {
            locale: locale.into(),
            entries,
        }
    }

    /// Load every `*.ftl` under `dir` and merge them into one
    /// catalogue. Useful when the locale is split across multiple
    /// files (`app.ftl`, `errors.ftl`, …) for editor ergonomics.
    pub fn from_dir(locale: impl Into<String>, dir: &Path) -> std::io::Result<Self> {
        let locale = locale.into();
        let mut merged: BTreeMap<String, String> = BTreeMap::new();
        if !dir.is_dir() {
            return Ok(Self {
                locale,
                entries: merged,
            });
        }
        // Round-19 L1: cap the read_dir iteration at MAX_LOCALE_FILES.
        // A poisoned locale dir with millions of `.ftl` files would
        // otherwise let the merge loop run unbounded, spending
        // forever in `from_str` before any UI rendering happens.
        for entry in std::fs::read_dir(dir)?.take(MAX_LOCALE_FILES) {
            let entry = entry?;
            let path = entry.path();
            let ext_match = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("ftl"))
                .unwrap_or(false);
            if !ext_match {
                continue;
            }
            let text = read_capped_ftl(&path)?;
            let part = Self::from_str(&locale, &text);
            for (k, v) in part.entries {
                merged.insert(k, v);
            }
        }
        Ok(Self {
            locale,
            entries: merged,
        })
    }

    /// Locale identifier (`"en-US"`, `"de-DE"`, …) — informational,
    /// not a lookup key.
    pub fn locale(&self) -> &str {
        &self.locale
    }

    /// Total entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when the catalogue has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Look up a key. Returns the matched string, or `None` when
    /// the key isn't present. Most call sites should use
    /// [`Self::lookup`] which falls back to the key itself.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(String::as_str)
    }

    /// Look up a key with fallback. When the key is missing the
    /// caller gets the literal key string back + a tracing warning.
    /// Lets the UI render *something* even when a translation is
    /// late — and the warning surfaces the gap during dev.
    pub fn lookup<'a>(&'a self, key: &'a str) -> &'a str {
        match self.entries.get(key) {
            Some(s) => s.as_str(),
            None => {
                tracing::warn!(
                    target: "valenx-i18n",
                    locale = self.locale,
                    key,
                    "missing translation key"
                );
                key
            }
        }
    }

    /// Look up a key + interpolate `{ $name }` placeholders. Each
    /// placeholder must appear in `args`; unfilled placeholders
    /// remain in the output (matching fluent-rs) so a missing arg
    /// is visible at runtime rather than silently swallowed.
    ///
    /// `args` is a slice of `(name, value)` pairs. Order doesn't
    /// matter — placeholders are matched by name.
    pub fn format_with<'a>(&'a self, key: &'a str, args: &[(&str, &str)]) -> String {
        let template = self.lookup(key);
        format_placeholders(template, args)
    }

    /// Insert a key explicitly. Useful for building a pseudo-locale
    /// catalogue programmatically (wrap every key in brackets to
    /// make hard-coded strings visually obvious).
    pub fn insert(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.entries.insert(key.into(), value.into());
    }

    /// Iterate `(key, value)` pairs in stable (BTreeMap) order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &str)> {
        self.entries.iter().map(|(k, v)| (k.as_str(), v.as_str()))
    }

    /// Build a pseudo-locale catalogue by wrapping every value in
    /// `⟦…⟧` brackets. Used in dev builds to make hard-coded
    /// strings visually obvious — anything not wrapped didn't
    /// route through the catalogue.
    pub fn to_pseudo(&self) -> Self {
        let entries = self
            .entries
            .iter()
            .map(|(k, v)| (k.clone(), format!("⟦{v}⟧")))
            .collect();
        Self {
            locale: format!("{}+pseudo", self.locale),
            entries,
        }
    }
}

/// Build a `LocaleCatalogue` from the en-US baseline `.ftl` file
/// shipped at `crates/valenx-i18n/locales/en-US/app.ftl`.
///
/// The string is `include_str!`-baked at compile time so the
/// caller doesn't need to ship the file alongside the binary or
/// resolve it at runtime. Future locales can ride a similar
/// helper (e.g. `embedded_de_de()`).
pub fn embedded_en_us() -> LocaleCatalogue {
    LocaleCatalogue::from_str("en-US", include_str!("../locales/en-US/app.ftl"))
}

/// Pure helper that interpolates `{ $name }` placeholders. Public
/// so tests can hit it without going through a `LocaleCatalogue`.
///
/// The byte cursor walks `template.as_bytes()` because the placeholder
/// markers are all single-byte ASCII (`{`, ` `, `$`, `}`), but when we
/// emit a non-placeholder character we MUST advance by the full UTF-8
/// codepoint or we'd splatter mojibake into the output (every byte
/// above 0x7F would become its own `char` instead of being decoded as
/// part of a multi-byte sequence — round-2 fixed the same bug in the
/// crash-reporter sanitisers; this is the matching i18n fix).
pub fn format_placeholders(template: &str, args: &[(&str, &str)]) -> String {
    let mut out = String::with_capacity(template.len());
    let bytes = template.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        // Look for the literal sequence `{` ` ` `$`.
        if i + 3 <= bytes.len() && bytes[i] == b'{' && bytes[i + 1] == b' ' && bytes[i + 2] == b'$'
        {
            // Find the closing ` }`.
            if let Some(end_rel) = find_close(&bytes[i + 3..]) {
                let name_end = i + 3 + end_rel;
                let name = std::str::from_utf8(&bytes[i + 3..name_end])
                    .unwrap_or("")
                    .trim();
                if let Some((_, value)) = args.iter().find(|(k, _)| *k == name) {
                    out.push_str(value);
                } else {
                    // Unfilled placeholder — preserve it so the gap
                    // is visible at runtime.
                    out.push_str(&template[i..name_end + 2]);
                }
                i = name_end + 2; // skip past " }"
                continue;
            }
        }
        // Not a placeholder start — copy the next UTF-8 codepoint
        // verbatim. `template[i..].chars().next()` is guaranteed to
        // return Some(c) because we've already bounds-checked `i <
        // bytes.len()` and `template` is &str (= valid UTF-8).
        let ch = template[i..]
            .chars()
            .next()
            .expect("i is in-bounds inside a &str");
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

/// Find the byte offset of the closing `" }"` (space + closing
/// brace) sequence. Returns `None` if not found before EOF.
fn find_close(bytes: &[u8]) -> Option<usize> {
    let mut i = 0;
    while i + 1 < bytes.len() {
        if bytes[i] == b' ' && bytes[i + 1] == b'}' {
            return Some(i);
        }
        i += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline_text() -> &'static str {
        r#"
# Header comment.
ribbon.run = Run
ribbon.prepare = Prepare

dialog.about.title = About Valenx
error.tool-not-installed = Tool `{ $tool }` is not installed.
error.run-failed = Run failed with exit code { $code }.

# Trailing comment is fine.
"#
    }

    #[test]
    fn from_str_skips_comments_and_blank_lines() {
        let cat = LocaleCatalogue::from_str("en-US", baseline_text());
        // baseline_text has 5 keys: ribbon.run, ribbon.prepare,
        // dialog.about.title, error.tool-not-installed,
        // error.run-failed. The header / inline / trailing comments
        // and blank lines drop.
        assert_eq!(cat.len(), 5);
        assert!(cat.get("ribbon.run").is_some());
    }

    #[test]
    fn from_str_handles_empty_input() {
        let cat = LocaleCatalogue::from_str("en-US", "");
        assert!(cat.is_empty());
    }

    #[test]
    fn lookup_falls_back_to_key_for_unknown() {
        let cat = LocaleCatalogue::from_str("en-US", baseline_text());
        assert_eq!(cat.lookup("ribbon.run"), "Run");
        // Missing key returns the key itself.
        assert_eq!(cat.lookup("ribbon.does-not-exist"), "ribbon.does-not-exist");
    }

    #[test]
    fn format_with_substitutes_placeholders() {
        let cat = LocaleCatalogue::from_str("en-US", baseline_text());
        let out = cat.format_with("error.tool-not-installed", &[("tool", "openfoam")]);
        assert_eq!(out, "Tool `openfoam` is not installed.");
    }

    #[test]
    fn format_with_handles_multiple_placeholders() {
        let cat = LocaleCatalogue::from_str("en-US", "msg = { $a } and { $b } at { $a } again");
        let out = cat.format_with("msg", &[("a", "first"), ("b", "second")]);
        assert_eq!(out, "first and second at first again");
    }

    #[test]
    fn format_with_preserves_unfilled_placeholders() {
        let cat = LocaleCatalogue::from_str("en-US", "msg = { $name } missing");
        let out = cat.format_with("msg", &[]);
        // Unfilled placeholder visible — easier to debug than a
        // silent empty cell.
        assert_eq!(out, "{ $name } missing");
    }

    #[test]
    fn format_with_treats_int_as_string() {
        // The format helper takes &str — callers convert. This test
        // just confirms the path round-trips.
        let cat = LocaleCatalogue::from_str("en-US", "msg = exit { $code }");
        let code_str = "42";
        let out = cat.format_with("msg", &[("code", code_str)]);
        assert_eq!(out, "exit 42");
    }

    /// Round-3 fix: `format_placeholders` previously walked the
    /// template byte-by-byte and emitted each byte as a `char` via
    /// `out.push(bytes[i] as char)`, which splits multi-byte UTF-8
    /// codepoints into mojibake. This test pins the codepoint-aware
    /// behaviour.
    #[test]
    fn format_placeholders_preserves_multibyte_unicode_in_template() {
        // German "Größe" — 'ö' is U+00F6 (two-byte UTF-8 0xC3 0xB6),
        // 'ß' is U+00DF (two-byte UTF-8 0xC3 0x9F). The byte-wise
        // implementation would produce "Grössse" — wait, worse than
        // that: 0xC3 / 0xB6 each become their own Latin-1 char and
        // the displayed string becomes "GrÃ¶ÃŸe".
        let template = "Hallo Größe { $what }";
        let out = format_placeholders(template, &[("what", "World")]);
        assert_eq!(out, "Hallo Größe World");
        // The literal mojibake substring should NOT appear.
        assert!(!out.contains('\u{00C3}'), "mojibake present: {out}");
    }

    #[test]
    fn format_placeholders_preserves_multibyte_unicode_with_no_placeholders() {
        // Sanity: even with no substitutions the multibyte chars are
        // preserved.
        let template = "日本語 → emoji 🎉";
        let out = format_placeholders(template, &[]);
        assert_eq!(out, "日本語 → emoji 🎉");
    }

    #[test]
    fn pseudo_locale_wraps_every_value() {
        let cat = LocaleCatalogue::from_str("en-US", "ribbon.run = Run\nribbon.prepare = Prepare");
        let pseudo = cat.to_pseudo();
        assert_eq!(pseudo.locale(), "en-US+pseudo");
        assert_eq!(pseudo.lookup("ribbon.run"), "⟦Run⟧");
        assert_eq!(pseudo.lookup("ribbon.prepare"), "⟦Prepare⟧");
    }

    #[test]
    fn from_str_skips_lines_without_equals() {
        let cat = LocaleCatalogue::from_str(
            "en-US",
            "ribbon.run = Run\nthis is not a valid line\nribbon.prepare = Prepare",
        );
        assert_eq!(cat.len(), 2);
    }

    #[test]
    fn from_str_skips_empty_keys() {
        let cat = LocaleCatalogue::from_str("en-US", "= no key here\n  = also no key");
        assert!(cat.is_empty());
    }

    #[test]
    fn iter_returns_pairs_in_stable_order() {
        let cat = LocaleCatalogue::from_str("en-US", "z = z\na = a\nm = m");
        let keys: Vec<&str> = cat.iter().map(|(k, _)| k).collect();
        // BTreeMap order — alphabetical.
        assert_eq!(keys, vec!["a", "m", "z"]);
    }

    #[test]
    fn from_dir_returns_empty_for_missing_dir() {
        let nope = std::env::temp_dir().join("valenx-i18n-does-not-exist-banana");
        let _ = std::fs::remove_dir_all(&nope);
        let cat = LocaleCatalogue::from_dir("en-US", &nope).expect("load");
        assert!(cat.is_empty());
    }

    #[test]
    fn from_dir_merges_multiple_ftl_files() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-i18n-merge-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("app.ftl"), "ribbon.run = Run\n").unwrap();
        std::fs::write(tmp.join("errors.ftl"), "error.io = IO\n").unwrap();
        // Skipped — wrong extension.
        std::fs::write(tmp.join("readme.md"), "ribbon.skip = nope\n").unwrap();

        let cat = LocaleCatalogue::from_dir("en-US", &tmp).expect("load");
        assert_eq!(cat.len(), 2);
        assert_eq!(cat.lookup("ribbon.run"), "Run");
        assert_eq!(cat.lookup("error.io"), "IO");
        assert_eq!(cat.lookup("ribbon.skip"), "ribbon.skip");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-18 L2 RED→GREEN: a `.ftl` larger than
    /// `MAX_FTL_FILE_BYTES` must be rejected by `from_dir` before the
    /// catalogue allocates the parsed entries.
    #[test]
    fn from_dir_rejects_oversize_ftl_file() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx-i18n-oversize-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let big = tmp.join("big.ftl");
        {
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&big)
                .unwrap();
            // Sparse-file trick — metadata.len() reports past-the-cap
            // while consuming ~0 disk blocks.
            f.seek(SeekFrom::Start(MAX_FTL_FILE_BYTES + 1)).unwrap();
            f.write_all(b"x").unwrap();
        }
        let r = LocaleCatalogue::from_dir("en-US", &tmp);
        let _ = std::fs::remove_dir_all(&tmp);
        let err = r.expect_err("must reject oversize .ftl");
        assert_eq!(
            err.kind(),
            std::io::ErrorKind::InvalidData,
            "expected InvalidData (size cap), got: {err}"
        );
    }

    #[test]
    fn embedded_en_us_round_trips_a_known_key() {
        let cat = embedded_en_us();
        // The baseline file ships these keys — surfaces a stale
        // baseline as a precise failure rather than a generic
        // "missing translation" warning.
        assert_eq!(cat.lookup("ribbon.run"), "Run");
        assert_eq!(cat.lookup("dialog.about.title"), "About Valenx");
    }

    #[test]
    fn baseline_locale_loads_without_warnings() {
        // The shipped en-US baseline should parse cleanly. Use a
        // path relative to CARGO_MANIFEST_DIR so the test works
        // regardless of the cwd.
        let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("locales")
            .join("en-US");
        let cat = LocaleCatalogue::from_dir("en-US", &dir).expect("load");
        assert!(!cat.is_empty(), "en-US baseline is empty — file missing?");
        // Spot-check every namespace surfaces at least one key.
        for prefix in [
            "ribbon.", "browser.", "status.", "dialog.", "palette.", "tooltip.", "error.",
        ] {
            let any = cat.iter().any(|(k, _)| k.starts_with(prefix));
            assert!(any, "no `{prefix}*` keys in baseline catalogue");
        }
    }
}

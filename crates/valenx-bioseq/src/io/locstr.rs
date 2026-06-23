//! Parser and writer for GenBank / EMBL feature-table location strings.
//!
//! Location strings follow the INSDC feature-table specification. The
//! subset supported here covers the full set of operators that occur in
//! real flat files:
//!
//! - `467` — a single base
//! - `12..78` — a simple range
//! - `<12..78` / `12..>78` — fuzzy (partial) ends (the `<`/`>` markers
//!   are accepted and dropped; the numeric value is kept)
//! - `complement(12..78)` — the reverse strand
//! - `join(12..78,134..202)` — a multi-segment feature (one contiguous
//!   biological entity)
//! - `order(12..78,134..202)` — a multi-segment listing with no joining
//!   claim
//! - `100^101` — the phosphodiester bond between two adjacent bases
//! - `complement(join(...))` and `join(complement(...),...)` nesting
//!
//! Cross-record references (`J00194.1:1..100`) are detected and
//! reported as a typed [`BioseqError::CrossRecordLocation`] so the
//! caller can fetch the referenced record and re-parse the location,
//! or skip the feature.
//!
//! Coordinates in the flat file are 1-based and inclusive; the parser
//! converts to the crate's 0-based half-open [`Span`] convention.

use crate::error::{BioseqError, Result};
use crate::record::{Location, Span, Strand};

/// Maximum operator-nesting depth the location parser will descend
/// before rejecting the input.
///
/// `parse_inner` recurses through `complement(`/`join(`/`order(`, so an
/// adversarial location string with thousands of nested operators
/// (`complement(complement(…))`, `join(join(…))`) would overflow the
/// call stack and *abort the process* — a stack overflow is not a
/// catchable panic. Real feature locations nest only a handful deep
/// (`complement(join(…))` is two levels; three is already exotic), so a
/// 1 000-level cap rejects only pathological / malicious input.
///
/// The cap fires well before the stack is exhausted: measured in a
/// debug build this parser overflows an 8 MiB stack (the common default
/// main-thread size) between ~5 000 and ~8 000 nested frames, so 1 000
/// keeps a ~5× margin — and release frames are smaller still.
const MAX_LOCATION_DEPTH: usize = 1_000;

/// Parses an INSDC location string into a [`Location`].
///
/// Returns [`BioseqError::CrossRecordLocation`] for a cross-record
/// reference (`accession:1..100`) and [`BioseqError::Parse`] for
/// syntactically broken input (including operator nesting deeper than
/// `MAX_LOCATION_DEPTH`).
pub fn parse_location(s: &str) -> Result<Location> {
    let trimmed = s.trim();
    parse_inner(trimmed, Strand::Forward, OpHint::None, 0)
}

/// Hint about which enclosing operator we are recursing under, so that
/// `order(...)` inside `complement(...)` survives the recursion as
/// `Order`, not `Join`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum OpHint {
    None,
    Join,
    Order,
}

/// Recursive worker. `strand` is the strand inherited from any
/// enclosing `complement(...)`; `hint` is the enclosing list operator;
/// `depth` is the current operator-nesting level, checked against
/// [`MAX_LOCATION_DEPTH`] so a maliciously deep input is rejected with
/// an error rather than overflowing the stack.
fn parse_inner(s: &str, strand: Strand, hint: OpHint, depth: usize) -> Result<Location> {
    if depth > MAX_LOCATION_DEPTH {
        return Err(BioseqError::parse(
            "location",
            format!("operator nesting too deep (max {MAX_LOCATION_DEPTH})"),
        ));
    }
    let s = s.trim();
    if let Some(rest) = strip_call(s, "complement") {
        // complement flips the strand for everything inside.
        return parse_inner(rest, strand.flip(), hint, depth + 1);
    }
    if let Some(rest) = strip_call(s, "join") {
        return parse_list(rest, strand, true, depth + 1);
    }
    if let Some(rest) = strip_call(s, "order") {
        return parse_list(rest, strand, false, depth + 1);
    }
    // A bare span — possibly between-bases. `hint` is not used here
    // because a bare span has no list-vs-`order` distinction to carry.
    let _ = hint;
    parse_span_or_between(s, strand)
}

/// Parses the contents of a `join(...)` (`is_join=true`) or
/// `order(...)` (`is_join=false`) call. `depth` is forwarded to the
/// recursive [`parse_inner`] calls for the nesting-depth guard.
fn parse_list(rest: &str, strand: Strand, is_join: bool, depth: usize) -> Result<Location> {
    let mut spans: Vec<Span> = Vec::new();
    // For order(...) we still allow nested join, complement, between.
    for part in split_top_level(rest) {
        let inner = parse_inner(
            &part,
            strand,
            if is_join { OpHint::Join } else { OpHint::Order },
            depth,
        )?;
        match inner {
            Location::Single(s) => spans.push(s),
            Location::Join(v) | Location::Order(v) => spans.extend(v),
            Location::Between { position, strand } => {
                // A between-bases entry inside a list becomes a
                // zero-length span at the bond position.
                spans.push(Span::with_strand(position, position, strand));
            }
        }
    }
    if spans.is_empty() {
        return Err(BioseqError::parse("location", "empty location list"));
    }
    // For a reverse-strand join, INSDC lists segments in reverse
    // genomic order; reverse so extraction concatenates correctly.
    if strand == Strand::Reverse {
        spans.reverse();
    }
    if spans.len() == 1 && is_join {
        // A one-element join is just the bare span.
        Ok(Location::Single(spans.into_iter().next().unwrap()))
    } else if is_join {
        Ok(Location::Join(spans))
    } else {
        Ok(Location::Order(spans))
    }
}

/// If `s` is `name(...)`, returns the inside of the parentheses.
fn strip_call<'a>(s: &'a str, name: &str) -> Option<&'a str> {
    let s = s.trim();
    let prefix = format!("{name}(");
    if s.starts_with(&prefix) && s.ends_with(')') {
        Some(&s[prefix.len()..s.len() - 1])
    } else {
        None
    }
}

/// Splits a comma-separated list, respecting nested parentheses.
fn split_top_level(s: &str) -> Vec<String> {
    let mut parts = Vec::new();
    let mut depth = 0i32;
    let mut cur = String::new();
    for ch in s.chars() {
        match ch {
            '(' => {
                depth += 1;
                cur.push(ch);
            }
            ')' => {
                depth -= 1;
                cur.push(ch);
            }
            ',' if depth == 0 => {
                parts.push(cur.trim().to_string());
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        parts.push(cur.trim().to_string());
    }
    parts
}

/// Parses a bare span; if the form is `n^m`, returns a
/// [`Location::Between`].
fn parse_span_or_between(s: &str, strand: Strand) -> Result<Location> {
    let s = s.trim();
    if s.contains(':') {
        // Cross-record reference. The accession is the text before
        // the colon; the remainder is the (sub-)location in the
        // referenced record.
        let (acc, _rest) = s.split_once(':').unwrap();
        return Err(BioseqError::cross_record_location(acc.trim(), s));
    }
    if let Some((a, b)) = s.split_once('^') {
        // n^n+1 — bond before base `n+1` (1-based). We accept any
        // adjacent pair n^n+1, including the circular wrap n^1 written
        // by some tools (we map that to position 0, i.e. the bond
        // before the first base / after the last).
        let na = parse_coord(a)?;
        let nb = parse_coord(b)?;
        if na == 0 {
            return Err(BioseqError::parse(
                "location",
                format!("invalid between-bases coordinate `{s}` (1-based)"),
            ));
        }
        // Conventional form: `100^101` -> bond between bases 100 and 101.
        // In 0-based half-open terms that is the boundary at index 100.
        // We accept the wrap form `n^1` and map it to position 0.
        let position = if na.checked_add(1) == Some(nb) {
            na // bond before 1-based base nb -> 0-based index na
        } else if nb == 1 && na > 1 {
            0 // n^1 wrap — bond at the origin
        } else {
            return Err(BioseqError::parse(
                "location",
                format!("`{s}` is not a valid between-bases form (need n^n+1)"),
            ));
        };
        return Ok(Location::Between { position, strand });
    }
    Ok(Location::Single(parse_span(s, strand)?))
}

/// Parses a single span: `467`, `12..78`, `<12..>78`.
fn parse_span(s: &str, strand: Strand) -> Result<Span> {
    let s = s.trim();
    if let Some((a, b)) = s.split_once("..") {
        let start = parse_coord(a)?;
        let end = parse_coord(b)?;
        if start == 0 || end == 0 || end < start {
            return Err(BioseqError::parse(
                "location",
                format!("invalid range `{s}` (1-based, start<=end)"),
            ));
        }
        // 1-based inclusive -> 0-based half-open.
        Ok(Span::with_strand(start - 1, end, strand))
    } else {
        // Single base.
        let p = parse_coord(s)?;
        if p == 0 {
            return Err(BioseqError::parse(
                "location",
                format!("invalid single-base coordinate `{s}`"),
            ));
        }
        Ok(Span::with_strand(p - 1, p, strand))
    }
}

/// Parses a coordinate, stripping any `<` / `>` fuzzy marker.
fn parse_coord(s: &str) -> Result<usize> {
    let cleaned: String = s.trim().chars().filter(|c| c.is_ascii_digit()).collect();
    cleaned
        .parse::<usize>()
        .map_err(|_| BioseqError::parse("location", format!("not a coordinate: `{s}`")))
}

/// Renders a [`Location`] back into an INSDC location string (the
/// inverse of [`parse_location`], 0-based half-open → 1-based
/// inclusive).
pub fn write_location(loc: &Location) -> String {
    // Render each span 1-based inclusive.
    let render_span = |s: &Span| -> String {
        if s.end == s.start + 1 {
            format!("{}", s.start + 1)
        } else {
            format!("{}..{}", s.start + 1, s.end)
        }
    };
    match loc {
        Location::Between { position, strand } => {
            // 0-based position p -> bond between bases p and p+1, both
            // 1-based. (i.e. `p^p+1`)
            let inner = format!("{}^{}", position, position + 1);
            if *strand == Strand::Reverse {
                format!("complement({inner})")
            } else {
                inner
            }
        }
        Location::Single(s) => {
            let reverse = s.strand == Strand::Reverse;
            let inner = render_span(s);
            if reverse {
                format!("complement({inner})")
            } else {
                inner
            }
        }
        Location::Join(spans) => write_list(spans, "join", render_span),
        Location::Order(spans) => write_list(spans, "order", render_span),
    }
}

/// Renders a list-valued location (join or order).
fn write_list<F>(spans: &[Span], op: &str, render: F) -> String
where
    F: Fn(&Span) -> String,
{
    let reverse = spans.iter().any(|s| s.strand == Strand::Reverse);
    // For reverse strand, segments were stored reverse-genomic;
    // emit them in ascending genomic order inside the call.
    let mut ordered: Vec<Span> = spans.to_vec();
    if reverse {
        ordered.sort_by_key(|s| s.start);
    }
    let joined = ordered.iter().map(&render).collect::<Vec<_>>().join(",");
    if reverse {
        format!("complement({op}({joined}))")
    } else {
        format!("{op}({joined})")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    #[test]
    fn single_base() {
        let loc = parse_location("467").unwrap();
        assert_eq!(loc, Location::Single(Span::new(466, 467)));
    }

    #[test]
    fn simple_range() {
        let loc = parse_location("12..78").unwrap();
        // 1-based inclusive 12..78 -> 0-based half-open 11..78.
        assert_eq!(loc, Location::Single(Span::new(11, 78)));
        assert_eq!(loc.total_len(), 67);
    }

    #[test]
    fn fuzzy_ends_are_accepted() {
        let loc = parse_location("<12..>78").unwrap();
        assert_eq!(loc, Location::Single(Span::new(11, 78)));
    }

    #[test]
    fn complement_flips_strand() {
        let loc = parse_location("complement(12..78)").unwrap();
        match loc {
            Location::Single(s) => {
                assert_eq!(s.strand, Strand::Reverse);
                assert_eq!((s.start, s.end), (11, 78));
            }
            _ => panic!("expected single"),
        }
    }

    #[test]
    fn join_multi_segment() {
        let loc = parse_location("join(1..3,10..12)").unwrap();
        let spans = loc.spans();
        assert_eq!(spans.len(), 2);
        assert_eq!((spans[0].start, spans[0].end), (0, 3));
        assert_eq!((spans[1].start, spans[1].end), (9, 12));
        assert!(matches!(loc, Location::Join(_)));
    }

    #[test]
    fn complement_of_join() {
        let loc = parse_location("complement(join(1..3,10..12))").unwrap();
        assert!(loc.is_reverse());
        // Reverse-strand join is stored reverse-genomic.
        let spans = loc.spans();
        assert_eq!(spans.len(), 2);
        assert_eq!(spans[0].start, 9);
        assert_eq!(spans[1].start, 0);
    }

    #[test]
    fn order_operator_parses_as_distinct_variant() {
        let loc = parse_location("order(1..3,10..12)").unwrap();
        assert!(
            matches!(loc, Location::Order(_)),
            "expected Order, got {loc:?}"
        );
        let spans = loc.spans();
        assert_eq!(spans.len(), 2);
        assert_eq!((spans[0].start, spans[0].end), (0, 3));
        assert_eq!((spans[1].start, spans[1].end), (9, 12));
    }

    #[test]
    fn order_under_complement() {
        let loc = parse_location("complement(order(1..3,10..12))").unwrap();
        assert!(loc.is_reverse());
        assert!(matches!(loc, Location::Order(_)));
    }

    #[test]
    fn between_bases_parses() {
        let loc = parse_location("100^101").unwrap();
        match loc {
            Location::Between { position, strand } => {
                assert_eq!(position, 100);
                assert_eq!(strand, Strand::Forward);
            }
            other => panic!("expected Between, got {other:?}"),
        }
        // Zero-length feature.
        assert_eq!(loc.total_len(), 0);
    }

    #[test]
    fn between_bases_complement() {
        let loc = parse_location("complement(100^101)").unwrap();
        match loc {
            Location::Between { position, strand } => {
                assert_eq!(position, 100);
                assert_eq!(strand, Strand::Reverse);
            }
            other => panic!("expected Between, got {other:?}"),
        }
    }

    #[test]
    fn between_bases_wrap_form() {
        // `5^1` is the circular wrap — bond at the origin.
        let loc = parse_location("5^1").unwrap();
        match loc {
            Location::Between { position, .. } => assert_eq!(position, 0),
            other => panic!("expected Between, got {other:?}"),
        }
    }

    #[test]
    fn between_bases_bad_form_errs() {
        // 100^102 — not adjacent.
        assert!(parse_location("100^102").is_err());
    }

    #[test]
    fn deeply_nested_location_is_rejected_not_a_stack_overflow() {
        // A pathological / malicious location string of deeply-nested
        // `complement(...)`. `parse_inner` recurses once per operator,
        // so an unbounded parser overflows the stack and *aborts the
        // process* (a stack overflow is NOT a catchable panic). `cap + 1`
        // nesting levels drive the recursion just past the guard, which
        // must return a clean `Err` reporting the nesting is too deep.
        // The cap is far below the stack-overflow floor, so these frames
        // fit comfortably on the default test stack.
        let depth = MAX_LOCATION_DEPTH + 1;
        let mut s = String::with_capacity(depth * 12);
        for _ in 0..depth {
            s.push_str("complement(");
        }
        s.push_str("1..1");
        for _ in 0..depth {
            s.push(')');
        }
        let err = parse_location(&s).expect_err("deep nesting must be rejected");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);
        assert!(
            err.to_string().contains("nesting too deep"),
            "expected a depth-limit error, got: {err}"
        );

        // The `join(...)` path (which also routes through `parse_list`)
        // is bounded by the same counter.
        let mut j = String::with_capacity(depth * 8);
        for _ in 0..depth {
            j.push_str("join(");
        }
        j.push_str("1..1,2..2");
        for _ in 0..depth {
            j.push(')');
        }
        let err = parse_location(&j).expect_err("deep join nesting must be rejected");
        assert!(
            err.to_string().contains("nesting too deep"),
            "expected a depth-limit error, got: {err}"
        );
    }

    #[test]
    fn realistic_nesting_still_parses() {
        // The deepest nesting that occurs in real flat files is two or
        // three operators — these must be unaffected by the depth guard.
        assert!(parse_location("complement(join(1..3,10..12))").is_ok());
        assert!(parse_location("join(complement(1..3),order(10..12,20..22))").is_ok());
    }

    #[test]
    fn cross_record_reference_is_typed_error() {
        let err = parse_location("J00194.1:1..10").unwrap_err();
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
        match err {
            BioseqError::CrossRecordLocation { accession, raw } => {
                assert_eq!(accession, "J00194.1");
                assert_eq!(raw, "J00194.1:1..10");
            }
            _ => panic!("wrong variant"),
        }
        // Inside a join, the cross-record reference still surfaces.
        let err = parse_location("join(1..10,J00194.1:1..10)").unwrap_err();
        assert!(matches!(err, BioseqError::CrossRecordLocation { .. }));
    }

    #[test]
    fn write_roundtrip_simple() {
        for s in ["467", "12..78"] {
            let loc = parse_location(s).unwrap();
            assert_eq!(write_location(&loc), s);
        }
    }

    #[test]
    fn write_roundtrip_complement_and_join() {
        let loc = parse_location("complement(12..78)").unwrap();
        assert_eq!(write_location(&loc), "complement(12..78)");
        let loc = parse_location("join(1..3,10..12)").unwrap();
        assert_eq!(write_location(&loc), "join(1..3,10..12)");
        let loc = parse_location("complement(join(1..3,10..12))").unwrap();
        assert_eq!(write_location(&loc), "complement(join(1..3,10..12))");
    }

    #[test]
    fn write_roundtrip_order_and_between() {
        let loc = parse_location("order(1..3,10..12)").unwrap();
        assert_eq!(write_location(&loc), "order(1..3,10..12)");
        let loc = parse_location("complement(order(1..3,10..12))").unwrap();
        assert_eq!(write_location(&loc), "complement(order(1..3,10..12))");
        let loc = parse_location("100^101").unwrap();
        assert_eq!(write_location(&loc), "100^101");
        let loc = parse_location("complement(100^101)").unwrap();
        assert_eq!(write_location(&loc), "complement(100^101)");
    }
}

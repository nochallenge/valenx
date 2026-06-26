//! # valenx-waveform
//!
//! In-house **digital waveform capture**: a self-contained
//! [Value Change Dump (VCD)][vcd] parser and signal-inspection API for HDL
//! simulation and logic-analysis traces — Valenx's digital-oscilloscope /
//! logic-analyzer input.
//!
//! ## What
//!
//! A VCD file is the lingua franca dumped by virtually every Verilog /
//! SystemVerilog / VHDL simulator (and by many logic analyzers): a small
//! declaration header that names the signals, followed by a body of
//! time-ordered value changes. [`Waveform::parse`] turns such a source
//! into a [`Waveform`] — a list of [`Signal`]s (each with a `name` and a
//! bit `width`), the overall [`Waveform::time_range`], and, per signal, the
//! full sequence of [`ValueChange`]s (a `time` paired with the new
//! [`SignalValue`]).
//!
//! Once parsed you can [`Waveform::signal`] / [`Waveform::signal_by_name`]
//! to find a wire and [`Signal::transitions`] to walk its edges, or
//! [`Signal::value_at`] to sample the value held at any time.
//!
//! ## Why in-house (not `wellen`)
//!
//! The obvious off-the-shelf option is [`wellen`](https://crates.io/crates/wellen),
//! the fast VCD/FST library behind the *surfer* viewer. Every currently
//! maintained `wellen` release (0.20–0.25) declares `rust-version = 1.90`,
//! which is **above the Valenx workspace MSRV of 1.88** — adopting it would
//! force an MSRV bump or pin us to the stale 0.19 line. VCD itself is a
//! simple, fully specified text format, so this crate parses it directly in
//! a few hundred dependency-light lines that build cleanly on the pinned
//! 1.88 toolchain. FST (the binary, compressed sibling format) is **not yet
//! supported** here and returns [`WaveformError::UnsupportedFormat`].
//!
//! ## Scope
//!
//! Research / educational grade. This is a pragmatic VCD reader covering the
//! constructs real simulators emit — `$var` / `$scope` / `$upscope` /
//! `$timescale` / `$enddefinitions`, scalar and vector value changes, and
//! the `$dumpvars` / `$dumpall` blocks — not an exhaustive IEEE-1364
//! conformance parser. It is *not* a clinical, safety- or production-grade
//! tool.
//!
//! ```
//! use valenx_waveform::{SignalValue, Waveform};
//!
//! // A 1-bit clock toggling and a 2-bit counter, dumped at #0/#5/#10.
//! let src = "\
//! $timescale 1ns $end
//! $scope module top $end
//! $var wire 1 ! clk $end
//! $var wire 2 # cnt [1:0] $end
//! $upscope $end
//! $enddefinitions $end
//! #0
//! 0!
//! b00 #
//! #5
//! 1!
//! b01 #
//! #10
//! 0!
//! b10 #
//! ";
//! let wf = Waveform::parse(src).unwrap();
//! assert_eq!(wf.signals().len(), 2);
//! assert_eq!(wf.time_range(), Some((0, 10)));
//!
//! let clk = wf.signal_by_name("clk").unwrap();
//! assert_eq!(clk.width, 1);
//! // clk: 0 @0, 1 @5, 0 @10
//! assert_eq!(clk.transitions().len(), 3);
//! assert_eq!(clk.value_at(7), Some(&SignalValue::bits("1")));
//! ```
//!
//! [vcd]: https://en.wikipedia.org/wiki/Value_change_dump

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors that can arise while parsing a VCD source.
///
/// Parsing never panics: every malformed or unsupported input is reported as
/// one of these variants instead.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum WaveformError {
    /// The source contained no `$var` declarations, so there are no signals
    /// to capture (an empty or header-less dump).
    #[error("no signals declared in VCD source (missing $var declarations)")]
    NoSignals,

    /// A `$var` declaration was syntactically malformed (wrong number of
    /// fields, or a non-numeric width). The offending text is included.
    #[error("malformed $var declaration: {0:?}")]
    MalformedVar(String),

    /// A `#<time>` timestamp could not be parsed as an unsigned integer.
    #[error("invalid timestamp: {0:?}")]
    InvalidTimestamp(String),

    /// A value-change line referenced an identifier code that was never
    /// declared by a `$var`.
    #[error("value change for unknown signal identifier {0:?}")]
    UnknownIdentifier(String),

    /// A value-change line was malformed (e.g. an empty scalar change, or a
    /// vector change missing its identifier).
    #[error("malformed value change: {0:?}")]
    MalformedValueChange(String),

    /// The input looks like a binary / non-VCD waveform format (for example
    /// FST), which this in-house reader does not handle.
    #[error("unsupported waveform format (only textual VCD is supported, not FST/binary)")]
    UnsupportedFormat,
}

/// Errors from loading a waveform off disk via [`Waveform::load`].
///
/// Separates I/O failures (file missing / unreadable) from the parse-level
/// [`WaveformError`]s, so callers can react to each. The I/O cause is kept as
/// a string so this type stays `Clone` + `PartialEq` like the rest of the
/// crate's public errors.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum LoadError {
    /// The file could not be read (does not exist, permission denied, …). The
    /// underlying [`std::io::Error`] message is captured as text.
    #[error("could not read waveform file: {0}")]
    Io(String),

    /// The file was read but its contents could not be parsed as VCD.
    #[error(transparent)]
    Parse(WaveformError),
}

/// The logic value a signal holds after a value change.
///
/// VCD distinguishes single-bit *scalar* changes from multi-bit *vector*
/// changes; both ultimately resolve to a 4-state bit string over the
/// alphabet `0 1 x z` (the conventional `X`/`Z` are normalised to lower
/// case). Real-valued (`r`/`R`) changes are kept verbatim as a [`Self::Real`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SignalValue {
    /// A 4-state bit string, most-significant bit first (e.g. `"10xz"`). A
    /// 1-bit signal carries a one-character string.
    Bits(String),
    /// A real-valued change, stored as the textual mantissa exactly as it
    /// appeared in the dump (e.g. `"3.14"`).
    Real(String),
}

impl SignalValue {
    /// Construct a [`Self::Bits`] value from a bit string, normalising the
    /// 4-state characters to lower case (`X` → `x`, `Z` → `z`).
    ///
    /// ```
    /// use valenx_waveform::SignalValue;
    /// assert_eq!(SignalValue::bits("1X0Z"), SignalValue::Bits("1x0z".to_string()));
    /// ```
    pub fn bits(s: &str) -> Self {
        let normalised: String = s
            .chars()
            .map(|c| match c {
                'X' => 'x',
                'Z' => 'z',
                other => other,
            })
            .collect();
        SignalValue::Bits(normalised)
    }

    /// The raw textual payload of this value (the bit string, or the real
    /// mantissa).
    pub fn as_str(&self) -> &str {
        match self {
            SignalValue::Bits(s) | SignalValue::Real(s) => s,
        }
    }
}

/// A single timestamped change of one signal's value.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValueChange {
    /// Simulation time of the change, in the dump's timescale units.
    pub time: u64,
    /// The value the signal takes from `time` onward (until the next change).
    pub value: SignalValue,
}

/// One declared waveform signal (a wire / variable) and its value history.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Signal {
    /// The signal's name as declared by `$var`, with its hierarchical scope
    /// prefixed and dotted (e.g. `"top.cnt"`). Any bit-select suffix
    /// (`[1:0]`) declared on the `$var` line is preserved.
    pub name: String,
    /// Bit width: `1` for a scalar wire, `N` for an `N`-bit vector.
    pub width: u32,
    /// The VCD identifier code used in the body to reference this signal
    /// (e.g. `"!"`, `"#"`). Exposed for traceability / debugging.
    pub id_code: String,
    /// The value changes for this signal, in ascending time order.
    changes: Vec<ValueChange>,
}

impl Signal {
    /// All value changes (transitions) for this signal, in ascending time
    /// order.
    pub fn transitions(&self) -> &[ValueChange] {
        &self.changes
    }

    /// The value held by the signal at `time` — that is, the value set by the
    /// most recent change at or before `time`.
    ///
    /// Returns `None` if the signal has no change at or before `time` (its
    /// value is undefined there). The lookup is a binary search over the
    /// time-ordered changes.
    pub fn value_at(&self, time: u64) -> Option<&SignalValue> {
        // Find the last change whose time is <= `time`.
        let idx = self.changes.partition_point(|c| c.time <= time);
        if idx == 0 {
            None
        } else {
            Some(&self.changes[idx - 1].value)
        }
    }

    /// The signal's first (earliest) recorded value, if any.
    pub fn initial_value(&self) -> Option<&SignalValue> {
        self.changes.first().map(|c| &c.value)
    }
}

/// A parsed digital waveform: a set of [`Signal`]s plus the overall time span.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Waveform {
    /// The timescale string as declared by `$timescale` (e.g. `"1ns"`), if
    /// present.
    pub timescale: Option<String>,
    /// The declared signals, in declaration order.
    signals: Vec<Signal>,
}

impl Waveform {
    /// Parse a VCD source held in a string.
    ///
    /// # Errors
    ///
    /// Returns a [`WaveformError`] for an empty / header-less dump
    /// ([`WaveformError::NoSignals`]), a malformed declaration or value
    /// change, an unparsable timestamp, or input that looks like a binary
    /// (non-VCD) format ([`WaveformError::UnsupportedFormat`]). It never
    /// panics.
    pub fn parse(source: &str) -> Result<Self, WaveformError> {
        // Reject obviously-binary input early: FST and other binary wave
        // formats are not textual VCD. A NUL byte never appears in a VCD.
        if source.as_bytes().contains(&0) {
            return Err(WaveformError::UnsupportedFormat);
        }

        let mut timescale: Option<String> = None;
        // id_code -> index into `signals` (multiple names may share a code
        // when the same net is dumped under several aliases; we keep the
        // first, but every alias resolves to the same change stream).
        let mut by_code: HashMap<String, usize> = HashMap::new();
        let mut signals: Vec<Signal> = Vec::new();

        let mut current_time: u64 = 0;
        let mut scope_stack: Vec<String> = Vec::new();
        let mut in_header = true;

        let mut tokens = Tokenizer::new(source);
        while let Some(tok) = tokens.next() {
            if in_header && tok.starts_with('$') {
                match tok {
                    "$timescale" => {
                        // Collect everything up to $end as the timescale text.
                        let body = tokens.take_until_end();
                        if !body.is_empty() {
                            timescale = Some(body.join(" "));
                        }
                    }
                    "$scope" => {
                        // $scope <type> <name> $end — the name is the 2nd word.
                        let body = tokens.take_until_end();
                        if let Some(name) = body.get(1) {
                            scope_stack.push((*name).to_string());
                        }
                    }
                    "$upscope" => {
                        let _ = tokens.take_until_end();
                        scope_stack.pop();
                    }
                    "$var" => {
                        let body = tokens.take_until_end();
                        let signal = parse_var(&body, &scope_stack)?;
                        let code = signal.id_code.clone();
                        // Only register a new entry per distinct id code.
                        by_code.entry(code).or_insert_with(|| {
                            signals.push(signal);
                            signals.len() - 1
                        });
                    }
                    "$enddefinitions" => {
                        let _ = tokens.take_until_end();
                        in_header = false;
                        if signals.is_empty() {
                            return Err(WaveformError::NoSignals);
                        }
                    }
                    // $date, $version, $comment, $dumpvars, $dumpall, $end …
                    // headers we don't need: skip to their $end.
                    _ => {
                        let _ = tokens.take_until_end();
                    }
                }
                continue;
            }

            // Body (and any value changes that appear inside $dumpvars even
            // before $enddefinitions in malformed-but-common dumps).
            if let Some(stripped) = tok.strip_prefix('#') {
                current_time = stripped
                    .parse::<u64>()
                    .map_err(|_| WaveformError::InvalidTimestamp(tok.to_string()))?;
            } else if tok == "$dumpvars"
                || tok == "$dumpall"
                || tok == "$dumpon"
                || tok == "$dumpoff"
                || tok == "$end"
                || tok == "$comment"
            {
                // Section markers in the body: nothing to record. ($comment
                // bodies are skipped token-by-token; they rarely collide with
                // value-change syntax in practice.)
                continue;
            } else if tok.starts_with('$') {
                // An unexpected header keyword after $enddefinitions: skip its
                // body defensively.
                let _ = tokens.take_until_end();
            } else {
                apply_value_change(tok, &mut tokens, current_time, &by_code, &mut signals)?;
            }
        }

        if signals.is_empty() {
            return Err(WaveformError::NoSignals);
        }

        Ok(Waveform { timescale, signals })
    }

    /// Load and parse a VCD waveform from a file on disk.
    ///
    /// A thin convenience over [`Waveform::parse`]: it reads the whole file
    /// into a string and parses it. The file extension is *not* inspected,
    /// but a binary (e.g. FST) payload is detected by content and reported as
    /// [`WaveformError::UnsupportedFormat`].
    ///
    /// # Errors
    ///
    /// Returns [`LoadError::Io`] if the file cannot be read, or
    /// [`LoadError::Parse`] wrapping any [`WaveformError`] from parsing.
    pub fn load<P: AsRef<std::path::Path>>(path: P) -> Result<Self, LoadError> {
        let text = std::fs::read_to_string(path).map_err(|e| LoadError::Io(e.to_string()))?;
        Waveform::parse(&text).map_err(LoadError::Parse)
    }

    /// All declared signals, in declaration order.
    pub fn signals(&self) -> &[Signal] {
        &self.signals
    }

    /// Look up a signal by its index in declaration order.
    pub fn signal(&self, index: usize) -> Option<&Signal> {
        self.signals.get(index)
    }

    /// Look up a signal by name.
    ///
    /// Matches either the fully-qualified dotted name (`"top.cnt"`) or the
    /// bare leaf name (`"cnt"`); the leaf match ignores any `[..]` bit-select
    /// suffix. Returns the first signal that matches.
    pub fn signal_by_name(&self, name: &str) -> Option<&Signal> {
        self.signals.iter().find(|s| {
            if s.name == name {
                return true;
            }
            let leaf = s.name.rsplit('.').next().unwrap_or(&s.name);
            // Strip a trailing bit-select like " [1:0]" or "[1:0]".
            let leaf_base = leaf.split('[').next().unwrap_or(leaf).trim();
            leaf_base == name
        })
    }

    /// The `(min_time, max_time)` spanned by any value change, or `None` if
    /// the waveform recorded no changes at all.
    pub fn time_range(&self) -> Option<(u64, u64)> {
        let mut min = u64::MAX;
        let mut max = u64::MIN;
        let mut seen = false;
        for s in &self.signals {
            if let (Some(first), Some(last)) = (s.changes.first(), s.changes.last()) {
                seen = true;
                min = min.min(first.time);
                max = max.max(last.time);
            }
        }
        if seen {
            Some((min, max))
        } else {
            None
        }
    }
}

/// Apply one value-change token (and, for vector/real changes, the following
/// identifier token) to the right signal.
fn apply_value_change(
    tok: &str,
    tokens: &mut Tokenizer<'_>,
    time: u64,
    by_code: &HashMap<String, usize>,
    signals: &mut [Signal],
) -> Result<(), WaveformError> {
    let first = tok.chars().next().expect("tokenizer never yields empty");
    let (value, id_code): (SignalValue, String) = match first {
        // Scalar change: a single value char immediately followed by the id
        // code, no space — e.g. `0!`, `1#`, `x$`, `z%`.
        '0' | '1' | 'x' | 'X' | 'z' | 'Z' => {
            let mut chars = tok.chars();
            let v = chars.next().unwrap();
            let code: String = chars.as_str().to_string();
            if code.is_empty() {
                return Err(WaveformError::MalformedValueChange(tok.to_string()));
            }
            (SignalValue::bits(&v.to_string()), code)
        }
        // Vector change: `b<bits>` then a space then the id code.
        'b' | 'B' => {
            let bits = &tok[1..];
            let code = tokens
                .next()
                .ok_or_else(|| WaveformError::MalformedValueChange(tok.to_string()))?;
            (SignalValue::bits(bits), code.to_string())
        }
        // Real change: `r<value>` then a space then the id code.
        'r' | 'R' => {
            let val = &tok[1..];
            let code = tokens
                .next()
                .ok_or_else(|| WaveformError::MalformedValueChange(tok.to_string()))?;
            (SignalValue::Real(val.to_string()), code.to_string())
        }
        _ => return Err(WaveformError::MalformedValueChange(tok.to_string())),
    };

    let idx = *by_code
        .get(&id_code)
        .ok_or(WaveformError::UnknownIdentifier(id_code.clone()))?;
    signals[idx].changes.push(ValueChange { time, value });
    Ok(())
}

/// Parse a `$var` declaration body (the words between `$var` and `$end`) into
/// a [`Signal`]. Expected form: `<type> <width> <id> <name> [bit-select]`.
fn parse_var(body: &[&str], scope: &[String]) -> Result<Signal, WaveformError> {
    // Minimum: type, width, id, name.
    if body.len() < 4 {
        return Err(WaveformError::MalformedVar(body.join(" ")));
    }
    let width: u32 = body[1]
        .parse()
        .map_err(|_| WaveformError::MalformedVar(body.join(" ")))?;
    let id_code = body[2].to_string();
    // The name is field 3; any further fields are a bit-select (e.g. `[1:0]`)
    // which we re-attach so the declared name round-trips.
    let leaf = body[3..].join(" ");

    let name = if scope.is_empty() {
        leaf
    } else {
        format!("{}.{}", scope.join("."), leaf)
    };

    Ok(Signal {
        name,
        width,
        id_code,
        changes: Vec::new(),
    })
}

/// A minimal whitespace tokenizer over a VCD source that also understands the
/// `… $end` keyword-terminated regions of the header.
struct Tokenizer<'a> {
    iter: std::str::SplitWhitespace<'a>,
}

impl<'a> Tokenizer<'a> {
    fn new(source: &'a str) -> Self {
        Tokenizer {
            iter: source.split_whitespace(),
        }
    }

    /// The next whitespace-delimited token, or `None` at end of input.
    fn next(&mut self) -> Option<&'a str> {
        self.iter.next()
    }

    /// Consume and collect tokens up to (and consuming) the next `$end`,
    /// returning the tokens in between. Stops at end of input if `$end` is
    /// absent.
    fn take_until_end(&mut self) -> Vec<&'a str> {
        let mut out = Vec::new();
        for t in self.iter.by_ref() {
            if t == "$end" {
                break;
            }
            out.push(t);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A clock + 2-bit counter dumped at #0/#5/#10 — the canonical fixture.
    const CLK_CNT: &str = "\
$date Mon $end
$version valenx $end
$timescale 1ns $end
$scope module top $end
$var wire 1 ! clk $end
$var reg 2 # cnt [1:0] $end
$upscope $end
$enddefinitions $end
$dumpvars
0!
b00 #
$end
#0
0!
b00 #
#5
1!
b01 #
#10
0!
b10 #
";

    #[test]
    fn parses_signal_names_and_widths() {
        let wf = Waveform::parse(CLK_CNT).unwrap();
        assert_eq!(wf.timescale.as_deref(), Some("1ns"));
        assert_eq!(wf.signals().len(), 2);

        let clk = wf.signal_by_name("clk").unwrap();
        assert_eq!(clk.name, "top.clk");
        assert_eq!(clk.width, 1);

        let cnt = wf.signal_by_name("cnt").unwrap();
        assert_eq!(cnt.name, "top.cnt [1:0]");
        assert_eq!(cnt.width, 2);
        // Also reachable by fully-qualified name.
        assert!(wf.signal_by_name("top.cnt [1:0]").is_some());
    }

    #[test]
    fn captures_transitions_at_right_times() {
        let wf = Waveform::parse(CLK_CNT).unwrap();
        let clk = wf.signal_by_name("clk").unwrap();
        // $dumpvars(@0) + #0 + #5 + #10 => four recorded changes for clk.
        let t = clk.transitions();
        assert_eq!(t.len(), 4);
        assert_eq!(
            t[0],
            ValueChange {
                time: 0,
                value: SignalValue::bits("0")
            }
        );
        assert_eq!(
            t[1],
            ValueChange {
                time: 0,
                value: SignalValue::bits("0")
            }
        );
        assert_eq!(
            t[2],
            ValueChange {
                time: 5,
                value: SignalValue::bits("1")
            }
        );
        assert_eq!(
            t[3],
            ValueChange {
                time: 10,
                value: SignalValue::bits("0")
            }
        );

        let cnt = wf.signal_by_name("cnt").unwrap();
        let ct = cnt.transitions();
        assert_eq!(ct.last().unwrap().value, SignalValue::bits("10"));
        assert_eq!(ct.last().unwrap().time, 10);
    }

    #[test]
    fn value_at_samples_held_value() {
        let wf = Waveform::parse(CLK_CNT).unwrap();
        let clk = wf.signal_by_name("clk").unwrap();
        // Held value between edges.
        assert_eq!(clk.value_at(0), Some(&SignalValue::bits("0")));
        assert_eq!(clk.value_at(3), Some(&SignalValue::bits("0")));
        assert_eq!(clk.value_at(5), Some(&SignalValue::bits("1")));
        assert_eq!(clk.value_at(7), Some(&SignalValue::bits("1")));
        assert_eq!(clk.value_at(100), Some(&SignalValue::bits("0")));

        let cnt = wf.signal_by_name("cnt").unwrap();
        assert_eq!(cnt.value_at(5), Some(&SignalValue::bits("01")));
        assert_eq!(cnt.value_at(9), Some(&SignalValue::bits("01")));
        assert_eq!(cnt.value_at(10), Some(&SignalValue::bits("10")));
    }

    #[test]
    fn time_range_spans_all_changes() {
        let wf = Waveform::parse(CLK_CNT).unwrap();
        assert_eq!(wf.time_range(), Some((0, 10)));
    }

    #[test]
    fn handles_x_z_and_real_values() {
        let src = "\
$timescale 1ps $end
$var wire 1 ! a $end
$var real 64 @ v $end
$enddefinitions $end
#0
x!
r3.14 @
#1
z!
r2.5 @
";
        let wf = Waveform::parse(src).unwrap();
        let a = wf.signal_by_name("a").unwrap();
        assert_eq!(a.value_at(0), Some(&SignalValue::bits("x")));
        assert_eq!(a.value_at(1), Some(&SignalValue::bits("z")));
        let v = wf.signal_by_name("v").unwrap();
        assert_eq!(v.value_at(0), Some(&SignalValue::Real("3.14".to_string())));
        assert_eq!(v.value_at(1), Some(&SignalValue::Real("2.5".to_string())));
    }

    #[test]
    fn empty_source_errors_no_signals() {
        assert_eq!(Waveform::parse(""), Err(WaveformError::NoSignals));
        assert_eq!(
            Waveform::parse("   \n\t  \n"),
            Err(WaveformError::NoSignals)
        );
    }

    #[test]
    fn header_without_vars_errors() {
        let src = "$timescale 1ns $end\n$enddefinitions $end\n#0\n";
        assert_eq!(Waveform::parse(src), Err(WaveformError::NoSignals));
    }

    #[test]
    fn malformed_var_errors_cleanly() {
        // Width is not a number.
        let src = "$var wire xx ! clk $end\n$enddefinitions $end\n";
        match Waveform::parse(src) {
            Err(WaveformError::MalformedVar(_)) => {}
            other => panic!("expected MalformedVar, got {other:?}"),
        }
        // Too few fields.
        let src2 = "$var wire 1 ! $end\n$enddefinitions $end\n";
        assert!(matches!(
            Waveform::parse(src2),
            Err(WaveformError::MalformedVar(_))
        ));
    }

    #[test]
    fn bad_timestamp_errors_cleanly() {
        let src = "$var wire 1 ! a $end\n$enddefinitions $end\n#notanumber\n0!\n";
        assert!(matches!(
            Waveform::parse(src),
            Err(WaveformError::InvalidTimestamp(_))
        ));
    }

    #[test]
    fn unknown_identifier_errors_cleanly() {
        let src = "$var wire 1 ! a $end\n$enddefinitions $end\n#0\n0?\n";
        assert_eq!(
            Waveform::parse(src),
            Err(WaveformError::UnknownIdentifier("?".to_string()))
        );
    }

    #[test]
    fn binary_input_rejected_as_unsupported() {
        // A leading NUL (as a real FST/binary blob would contain) is rejected
        // without panicking.
        let src = "\u{0}\u{1}FST\u{0}garbage";
        assert_eq!(Waveform::parse(src), Err(WaveformError::UnsupportedFormat));
    }

    #[test]
    fn load_from_path_round_trips() {
        let mut path = std::env::temp_dir();
        path.push(format!("valenx_waveform_test_{}.vcd", std::process::id()));
        std::fs::write(&path, CLK_CNT).unwrap();

        let wf = Waveform::load(&path).unwrap();
        assert_eq!(wf.signals().len(), 2);
        assert_eq!(wf.time_range(), Some((0, 10)));

        let _ = std::fs::remove_file(&path);

        // Missing file -> Io error, not a panic.
        let missing = std::env::temp_dir().join("valenx_waveform_does_not_exist.vcd");
        assert!(matches!(Waveform::load(missing), Err(LoadError::Io(_))));
    }

    #[test]
    fn signals_are_serializable() {
        let wf = Waveform::parse(CLK_CNT).unwrap();
        let json = serde_json::to_string(&wf).unwrap();
        let round: Waveform = serde_json::from_str(&json).unwrap();
        assert_eq!(wf, round);
    }
}

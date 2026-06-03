//! `valenx-sam-info <path|->` — summary stats on a SAM file.
//!
//! Counts records, mapped vs unmapped, and header lines. `--format
//! json` emits a JSON envelope for CI consumption. BAM (binary) is
//! out of scope — convert with `samtools view -h` first.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: valenx-sam-info [--format text|json] <path|->");
        return ExitCode::from(2);
    }
    let mut format = "text";
    let mut path: Option<String> = None;
    let mut iter = args.iter();
    while let Some(a) = iter.next() {
        match a.as_str() {
            "--format" => match iter.next().map(|s| s.as_str()) {
                Some("text") => format = "text",
                Some("json") => format = "json",
                _ => {
                    eprintln!("--format takes `text` or `json`");
                    return ExitCode::from(2);
                }
            },
            other => path = Some(other.to_string()),
        }
    }
    let path = match path {
        Some(p) => p,
        None => {
            eprintln!("usage: valenx-sam-info [--format text|json] <path|->");
            return ExitCode::from(2);
        }
    };
    let body = match read_input(&path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("read failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let sam = match valenx_bio::format::sam::read_str(&body) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("parse failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mapped = sam.records.iter().filter(|r| r.is_mapped()).count();
    let unmapped = sam.records.len() - mapped;
    if format == "json" {
        let v = serde_json::json!({
            "header_lines": sam.header.len(),
            "records": sam.records.len(),
            "mapped": mapped,
            "unmapped": unmapped,
        });
        println!("{}", serde_json::to_string_pretty(&v).unwrap());
    } else {
        println!("header lines: {}", sam.header.len());
        println!("records: {}", sam.records.len());
        println!("mapped: {mapped}");
        println!("unmapped: {unmapped}");
    }
    ExitCode::SUCCESS
}

fn read_input(path: &str) -> std::io::Result<String> {
    // Round-21 M4: see valenx_fastq for the cap rationale.
    if path == "-" {
        valenx_core::io_caps::read_capped_stdin_to_string(valenx_core::io_caps::MAX_BIO_CLI_BYTES)
    } else {
        valenx_core::io_caps::read_capped_to_string(
            std::path::Path::new(path),
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        )
    }
}

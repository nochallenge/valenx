//! `valenx-fastq inspect <path|->` — quick stats on a FASTQ file.
//!
//! Reads from a path or `-` (stdin). Prints record count, total
//! bases, length range, and Phred quality range. `--format json`
//! emits a flat JSON envelope for CI consumption.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(|s| s.as_str()) != Some("inspect") {
        eprintln!("usage: valenx-fastq inspect [--format text|json] <path|->");
        return ExitCode::from(2);
    }
    let mut format = "text";
    let mut path: Option<String> = None;
    let mut iter = args.iter().skip(1);
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
            eprintln!("usage: valenx-fastq inspect [--format text|json] <path|->");
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
    let recs = match valenx_bio::format::fastq::read_str(&body) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("parse failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let total: usize = recs.iter().map(|r| r.len()).sum();
    let min_len = recs.iter().map(|r| r.len()).min().unwrap_or(0);
    let max_len = recs.iter().map(|r| r.len()).max().unwrap_or(0);
    let min_q = recs.iter().filter_map(|r| r.min_quality()).min();
    let max_q = recs.iter().filter_map(|r| r.max_quality()).max();
    if format == "json" {
        let v = serde_json::json!({
            "records": recs.len(),
            "total_bases": total,
            "min_length": min_len,
            "max_length": max_len,
            "min_quality": min_q,
            "max_quality": max_q,
        });
        println!("{}", serde_json::to_string_pretty(&v).unwrap());
    } else {
        println!("records: {}", recs.len());
        println!("total bases: {total}");
        println!("length range: {min_len}..={max_len}");
        match (min_q, max_q) {
            (Some(a), Some(b)) => println!("quality (Phred): {a}..={b}"),
            _ => println!("quality (Phred): n/a"),
        }
    }
    ExitCode::SUCCESS
}

fn read_input(path: &str) -> std::io::Result<String> {
    // Round-21 M4: bound at MAX_BIO_CLI_BYTES so a piped
    // `cat /dev/zero | valenx-fastq inspect -` produces a clean
    // error rather than OOM.
    if path == "-" {
        valenx_core::io_caps::read_capped_stdin_to_string(valenx_core::io_caps::MAX_BIO_CLI_BYTES)
    } else {
        valenx_core::io_caps::read_capped_to_string(
            std::path::Path::new(path),
            valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
        )
    }
}

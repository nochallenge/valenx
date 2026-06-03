//! `valenx-vcf-info <path|->` — summary stats on a VCF file.
//!
//! Counts header lines, samples, records, and PASS/fail + no-ALT
//! splits. `--format json` emits a JSON envelope for CI consumption.
//! BCF (binary) and BGZF-compressed VCF are out of scope — convert
//! with `bcftools view -O v` first.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: valenx-vcf-info [--format text|json] <path|->");
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
            eprintln!("usage: valenx-vcf-info [--format text|json] <path|->");
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
    let v = match valenx_bio::format::vcf::read_str(&body) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("parse failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let pass = v.records.iter().filter(|r| r.is_pass()).count();
    let fail = v.records.len() - pass;
    let no_alt = v.records.iter().filter(|r| !r.has_alt()).count();
    if format == "json" {
        let env = serde_json::json!({
            "header_lines": v.header.len(),
            "samples": v.samples.len(),
            "records": v.records.len(),
            "pass": pass,
            "fail": fail,
            "no_alt": no_alt,
        });
        println!("{}", serde_json::to_string_pretty(&env).unwrap());
    } else {
        println!("header lines: {}", v.header.len());
        println!("samples: {}", v.samples.len());
        println!("records: {}", v.records.len());
        println!("pass: {pass}");
        println!("fail: {fail}");
        println!("no-alt: {no_alt}");
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

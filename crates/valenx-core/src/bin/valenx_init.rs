//! `valenx-init` — scaffold a fresh `.valenx` project skeleton.
//!
//! Emits the minimal file tree the project loader requires, with
//! a usable `case.toml` chosen from a small template library
//! (cfd / fea / chemistry / empty / one entry per registered
//! adapter). Saves users from copy-pasting the boilerplate out of
//! QUICKSTART.md.
//!
//! Usage:
//!
//! ```text
//! valenx-init <dir> [--template cfd|fea|chemistry|empty|...] [--name <project-name>]
//! valenx-init help
//! ```
//!
//! Exit codes:
//! - 0: project written
//! - 1: target dir already populated / IO failure
//! - 2: invalid CLI usage
//!
//! Templates + rendering live in [`valenx_core::init_templates`] —
//! this binary is a thin CLI shell over that library so the desktop
//! GUI can drive the same scaffolding from its "New case from
//! adapter" flow.
//!
//! No `clap` dep — argument parsing is a single match expression.

use std::path::PathBuf;
use std::process::ExitCode;

use valenx_core::init_templates::{render_template_list, scaffold_project, Template, USAGE};

#[derive(Debug, PartialEq, Eq)]
enum ParsedArgs {
    Help,
    Version,
    ListTemplates,
    Init {
        dir: PathBuf,
        template: Template,
        name: Option<String>,
    },
    Invalid(String),
}

fn parse_args(args: &[String]) -> ParsedArgs {
    if args.is_empty() {
        return ParsedArgs::Invalid("missing target directory".into());
    }
    if matches!(args[0].as_str(), "help" | "-h" | "--help") {
        return ParsedArgs::Help;
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        return ParsedArgs::Version;
    }
    if matches!(
        args[0].as_str(),
        "--list-templates" | "-l" | "list-templates"
    ) {
        return ParsedArgs::ListTemplates;
    }
    let dir = PathBuf::from(&args[0]);
    let mut template = Template::Empty;
    let mut name: Option<String> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--template" | "-t" => {
                i += 1;
                let Some(t) = args.get(i) else {
                    return ParsedArgs::Invalid("--template needs a value".into());
                };
                let Some(parsed) = Template::from_str(t) else {
                    return ParsedArgs::Invalid(format!(
                        "unknown template `{t}` — run `valenx-init --list-templates` \
                         to see what's available, or `valenx-init help` for aliases"
                    ));
                };
                template = parsed;
            }
            "--name" | "-n" => {
                i += 1;
                let Some(v) = args.get(i) else {
                    return ParsedArgs::Invalid("--name needs a value".into());
                };
                name = Some(v.clone());
            }
            other => {
                return ParsedArgs::Invalid(format!("unknown argument `{other}`"));
            }
        }
        i += 1;
    }
    ParsedArgs::Init {
        dir,
        template,
        name,
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        ParsedArgs::Help => {
            print!("{USAGE}");
            ExitCode::from(0)
        }
        ParsedArgs::Version => {
            println!("valenx-init v{}", env!("CARGO_PKG_VERSION"));
            ExitCode::from(0)
        }
        ParsedArgs::ListTemplates => {
            print!("{}", render_template_list());
            ExitCode::from(0)
        }
        ParsedArgs::Init {
            dir,
            template,
            name,
        } => match scaffold_project(&dir, template, name.as_deref()) {
            Ok(()) => {
                println!("Initialised valenx project at {}", dir.display());
                // Next-step hints — chain `valenx-validate` to confirm
                // the project loads cleanly, and point at the case
                // file for hand-editing.
                println!();
                println!("Next steps:");
                println!(
                    "  valenx-validate {}    # confirm the project loads",
                    dir.display()
                );
                println!(
                    "  $EDITOR {}        # tweak the starter case",
                    dir.join("cases").display()
                );
                ExitCode::from(0)
            }
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::from(1)
            }
        },
        ParsedArgs::Invalid(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_args_help_variants() {
        assert_eq!(parse_args(&["help".into()]), ParsedArgs::Help);
        assert_eq!(parse_args(&["-h".into()]), ParsedArgs::Help);
        assert_eq!(parse_args(&["--help".into()]), ParsedArgs::Help);
    }

    #[test]
    fn parse_args_version_variants() {
        assert_eq!(parse_args(&["-V".into()]), ParsedArgs::Version);
        assert_eq!(parse_args(&["--version".into()]), ParsedArgs::Version);
    }

    #[test]
    fn parse_args_list_templates_variants() {
        assert_eq!(
            parse_args(&["--list-templates".into()]),
            ParsedArgs::ListTemplates
        );
        assert_eq!(parse_args(&["-l".into()]), ParsedArgs::ListTemplates);
        assert_eq!(
            parse_args(&["list-templates".into()]),
            ParsedArgs::ListTemplates
        );
    }

    #[test]
    fn parse_args_dir_only_uses_empty_template() {
        match parse_args(&["./out".into()]) {
            ParsedArgs::Init {
                dir,
                template,
                name,
            } => {
                assert_eq!(dir, PathBuf::from("./out"));
                assert_eq!(template, Template::Empty);
                assert!(name.is_none());
            }
            other => panic!("wrong: {other:?}"),
        }
    }

    #[test]
    fn parse_args_template_aliases_normalise() {
        match parse_args(&["out".into(), "--template".into(), "openfoam".into()]) {
            ParsedArgs::Init { template, .. } => assert_eq!(template, Template::Cfd),
            other => panic!("wrong: {other:?}"),
        }
        match parse_args(&["out".into(), "--template".into(), "structural".into()]) {
            ParsedArgs::Init { template, .. } => assert_eq!(template, Template::Fea),
            other => panic!("wrong: {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_template_is_invalid() {
        match parse_args(&["out".into(), "--template".into(), "rocket-science".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown template")),
            other => panic!("wrong: {other:?}"),
        }
    }

    #[test]
    fn parse_args_missing_value_for_template() {
        match parse_args(&["out".into(), "--template".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("--template needs")),
            other => panic!("wrong: {other:?}"),
        }
    }

    #[test]
    fn parse_args_name_overrides_dir_name() {
        match parse_args(&["ignore-me".into(), "--name".into(), "real-name".into()]) {
            ParsedArgs::Init { name, .. } => assert_eq!(name.as_deref(), Some("real-name")),
            other => panic!("wrong: {other:?}"),
        }
    }

    #[test]
    fn parse_args_unknown_arg_is_invalid() {
        match parse_args(&["out".into(), "--bogus".into()]) {
            ParsedArgs::Invalid(msg) => assert!(msg.contains("unknown argument")),
            other => panic!("wrong: {other:?}"),
        }
    }
}

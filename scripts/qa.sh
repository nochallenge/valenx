#!/usr/bin/env bash
# =============================================================================
# Valenx local QA harness — the one-command "is the project healthy?" runner.
#
#   ./scripts/qa.sh            run the whole safe validation suite
#   ./scripts/qa.sh --tests    only the scoped per-crate test runs
#   ./scripts/qa.sh --gates    only the workspace check / clippy / doc / deny gates
#   ./scripts/qa.sh --help     this message
#
# WHY THIS SCRIPT IS SCOPED — READ docs/QA.md.
# A blanket `cargo test --workspace` is FORBIDDEN: the `valenx-app` library's
# UI-coupled unit tests call `rfd::FileDialog`, which pops a native OS file
# dialog that blocks forever in a headless run — and once wedged a machine.
# This harness therefore runs `cargo test` ONLY:
#   * `-p <crate>` for the 20 pure computational crates (none link `rfd`);
#   * `-p valenx-app headless_ui_tests` — a NAME FILTER that selects only the
#     windowless egui-logic tests and excludes every file-dialog test;
#   * `-p valenx-app --test pipeline_e2e` — ONE integration-test file, the
#     cross-crate end-to-end suite, which never touches `rfd`.
# It NEVER runs `cargo test --workspace`, unfiltered `valenx-app` tests,
# `cargo run`, `cargo bench`, or launches the app binary.
# =============================================================================

set -u
cd "$(dirname "$0")/.."   # repo root, regardless of where we were invoked

# --- the 20 pure computational crates: safe to `cargo test -p` ---------------
# These are pure-Rust algorithm libraries — no `rfd`, no GUI, no subprocess.
PURE_CRATES=(
  valenx-bioseq
  valenx-align
  valenx-phylo
  valenx-popgen
  valenx-rnastruct
  valenx-md
  valenx-cheminf
  valenx-biostruct
  valenx-qchem
  valenx-genomics
  valenx-sysbio
  valenx-dock-screen
  valenx-genediting
  valenx-structpredict
  valenx-rnadesign
  valenx-aero
  valenx-cfd-native
  valenx-fem
  valenx-pathtrace
  valenx-render-bridge
)

GREEN=$'\033[32m'; RED=$'\033[31m'; BOLD=$'\033[1m'; DIM=$'\033[2m'; RST=$'\033[0m'
FAILURES=()
START=$(date +%s)

run() {  # run <label> <cmd...>
  local label="$1"; shift
  printf '%s>>> %s%s\n' "$BOLD" "$label" "$RST"
  if "$@"; then
    printf '%s    PASS%s  %s\n' "$GREEN" "$RST" "$label"
  else
    printf '%s    FAIL%s  %s\n' "$RED" "$RST" "$label"
    FAILURES+=("$label")
  fi
  echo
}

do_tests() {
  echo "${BOLD}== Scoped per-crate tests (the 20 pure computational crates) ==${RST}"
  echo "${DIM}   one crate at a time — never 'cargo test --workspace'${RST}"
  echo
  for crate in "${PURE_CRATES[@]}"; do
    run "cargo test -p $crate" cargo test -p "$crate"
  done
  echo "${BOLD}== Workbench UI-logic tests (name-filtered — no file dialogs) ==${RST}"
  echo
  # The `headless_ui_tests` filter selects ONLY the windowless egui-logic
  # tests; valenx-app's lib unit tests (which open rfd dialogs) are excluded.
  run "cargo test -p valenx-app headless_ui_tests" \
      cargo test -p valenx-app headless_ui_tests
  echo "${BOLD}== Cross-crate end-to-end pipeline tests (single test file) ==${RST}"
  echo
  # `--test pipeline_e2e` compiles + runs ONLY that one integration file.
  run "cargo test -p valenx-app --test pipeline_e2e" \
      cargo test -p valenx-app --test pipeline_e2e
}

do_gates() {
  echo "${BOLD}== Workspace gates (check / clippy / doc / deny — build-only, no run) ==${RST}"
  echo
  run "cargo check --workspace"  cargo check --workspace
  run "cargo clippy --workspace --all-targets -- -D warnings" \
      cargo clippy --workspace --all-targets -- -D warnings
  # NOTE: ~5 pre-existing rustdoc warnings in the untouched `valenx-solvespace-3d`
  # crate are a known baseline — `cargo doc` is not failed on them here.
  run "cargo doc --workspace --no-deps" cargo doc --workspace --no-deps
  # Round-8 L18: cargo-deny check delivers the supply-chain audit
  # POLICIES.md / SECURITY.md / CHANGELOG.md claim. Pre-fix the docs
  # advertised `cargo deny check` but `qa.sh` didn't run it; now it
  # does. Skip silently when the tool isn't installed so contributor
  # boxes without `cargo-deny` aren't blocked from running the rest
  # of the suite (CI installs it explicitly).
  if command -v cargo-deny >/dev/null 2>&1; then
    run "cargo deny check"  cargo deny check
  else
    printf '%s>>> cargo deny check%s  %s(skipped — `cargo install cargo-deny` to enable)%s\n\n' \
      "$BOLD" "$RST" "$DIM" "$RST"
  fi
}

case "${1:-}" in
  --help|-h)
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    exit 0 ;;
  --tests)  do_tests ;;
  --gates)  do_gates ;;
  "")       do_tests; do_gates ;;
  *)
    echo "unknown option: $1 (try --help)" >&2
    exit 2 ;;
esac

ELAPSED=$(( $(date +%s) - START ))
echo "============================================================"
if [ ${#FAILURES[@]} -eq 0 ]; then
  printf '%sALL QA STEPS PASSED%s  (%ss)\n' "$GREEN$BOLD" "$RST" "$ELAPSED"
  exit 0
else
  printf '%s%d QA STEP(S) FAILED%s  (%ss)\n' "$RED$BOLD" "${#FAILURES[@]}" "$RST" "$ELAPSED"
  for f in "${FAILURES[@]}"; do echo "  - $f"; done
  exit 1
fi

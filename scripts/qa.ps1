<#
.SYNOPSIS
    Valenx local QA harness — the one-command "is the project healthy?" runner.

.DESCRIPTION
    Runs the entire SAFE validation suite: scoped per-crate tests for the 20
    pure computational crates, the name-filtered workbench UI-logic tests, the
    cross-crate end-to-end pipeline tests, then the workspace check / clippy /
    doc / deny gates.

    WHY THIS SCRIPT IS SCOPED — READ docs/QA.md.
    A blanket `cargo test --workspace` is FORBIDDEN: the `valenx-app` library's
    UI-coupled unit tests call `rfd::FileDialog`, which pops a native OS file
    dialog that blocks forever in a headless run — and once wedged a machine.
    This harness therefore runs `cargo test` ONLY:
      * `-p <crate>` for the 20 pure computational crates (none link `rfd`);
      * `-p valenx-app headless_ui_tests` — a NAME FILTER selecting only the
        windowless egui-logic tests, excluding every file-dialog test;
      * `-p valenx-app --test pipeline_e2e` — ONE integration-test file.
    It NEVER runs `cargo test --workspace`, unfiltered `valenx-app` tests,
    `cargo run`, `cargo bench`, or launches the app binary.

.PARAMETER Tests
    Run only the scoped per-crate test runs.

.PARAMETER Gates
    Run only the workspace check / clippy / doc / deny gates.

.EXAMPLE
    .\scripts\qa.ps1
    Run the whole safe validation suite.

.EXAMPLE
    .\scripts\qa.ps1 -Tests
    Only the scoped per-crate tests.
#>
[CmdletBinding()]
param(
    [switch]$Tests,
    [switch]$Gates
)

$ErrorActionPreference = 'Continue'
# Repo root, regardless of where the script was invoked from.
Set-Location (Join-Path $PSScriptRoot '..')

# --- the 20 pure computational crates: safe to `cargo test -p` ---------------
# These are pure-Rust algorithm libraries — no `rfd`, no GUI, no subprocess.
$PureCrates = @(
    'valenx-bioseq'
    'valenx-align'
    'valenx-phylo'
    'valenx-popgen'
    'valenx-rnastruct'
    'valenx-md'
    'valenx-cheminf'
    'valenx-biostruct'
    'valenx-qchem'
    'valenx-genomics'
    'valenx-sysbio'
    'valenx-dock-screen'
    'valenx-genediting'
    'valenx-structpredict'
    'valenx-rnadesign'
    'valenx-aero'
    'valenx-cfd-native'
    'valenx-fem'
    'valenx-pathtrace'
    'valenx-render-bridge'
)

$script:Failures = @()
$start = Get-Date

function Invoke-Step {
    param([string]$Label, [string[]]$CargoArgs)
    Write-Host ">>> $Label" -ForegroundColor White
    & cargo @CargoArgs
    if ($LASTEXITCODE -eq 0) {
        Write-Host "    PASS  $Label" -ForegroundColor Green
    } else {
        Write-Host "    FAIL  $Label" -ForegroundColor Red
        $script:Failures += $Label
    }
    Write-Host ''
}

function Invoke-Tests {
    Write-Host '== Scoped per-crate tests (the 20 pure computational crates) ==' -ForegroundColor White
    Write-Host '   one crate at a time - never "cargo test --workspace"' -ForegroundColor DarkGray
    Write-Host ''
    foreach ($crate in $PureCrates) {
        Invoke-Step "cargo test -p $crate" @('test', '-p', $crate)
    }
    Write-Host '== Workbench UI-logic tests (name-filtered - no file dialogs) ==' -ForegroundColor White
    Write-Host ''
    # The `headless_ui_tests` filter selects ONLY the windowless egui-logic
    # tests; valenx-app's lib unit tests (which open rfd dialogs) are excluded.
    Invoke-Step 'cargo test -p valenx-app headless_ui_tests' `
        @('test', '-p', 'valenx-app', 'headless_ui_tests')
    Write-Host '== Cross-crate end-to-end pipeline tests (single test file) ==' -ForegroundColor White
    Write-Host ''
    # `--test pipeline_e2e` compiles + runs ONLY that one integration file.
    Invoke-Step 'cargo test -p valenx-app --test pipeline_e2e' `
        @('test', '-p', 'valenx-app', '--test', 'pipeline_e2e')
}

function Invoke-Gates {
    Write-Host '== Workspace gates (check / clippy / doc / deny - build-only, no run) ==' -ForegroundColor White
    Write-Host ''
    Invoke-Step 'cargo check --workspace' @('check', '--workspace')
    Invoke-Step 'cargo clippy --workspace --all-targets -- -D warnings' `
        @('clippy', '--workspace', '--all-targets', '--', '-D', 'warnings')
    # NOTE: ~5 pre-existing rustdoc warnings in the untouched
    # `valenx-solvespace-3d` crate are a known baseline.
    Invoke-Step 'cargo doc --workspace --no-deps' @('doc', '--workspace', '--no-deps')
    # Round-8 L18 (PowerShell mirror): cargo-deny check delivers the
    # supply-chain audit POLICIES.md / SECURITY.md / CHANGELOG.md
    # claim. Skipped silently when cargo-deny isn't installed so
    # contributor boxes without it aren't blocked from running the
    # rest of the suite (CI installs it explicitly).
    if (Get-Command cargo-deny -ErrorAction SilentlyContinue) {
        Invoke-Step 'cargo deny check' @('deny', 'check')
    } else {
        Write-Host '>>> cargo deny check' -ForegroundColor White
        Write-Host '    (skipped - "cargo install cargo-deny" to enable)' -ForegroundColor DarkGray
        Write-Host ''
    }
}

if ($Tests) {
    Invoke-Tests
} elseif ($Gates) {
    Invoke-Gates
} else {
    Invoke-Tests
    Invoke-Gates
}

$elapsed = [int]((Get-Date) - $start).TotalSeconds
Write-Host '============================================================'
if ($script:Failures.Count -eq 0) {
    Write-Host "ALL QA STEPS PASSED  (${elapsed}s)" -ForegroundColor Green
    exit 0
} else {
    Write-Host "$($script:Failures.Count) QA STEP(S) FAILED  (${elapsed}s)" -ForegroundColor Red
    foreach ($f in $script:Failures) { Write-Host "  - $f" }
    exit 1
}

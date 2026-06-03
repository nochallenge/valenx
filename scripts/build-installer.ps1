#!/usr/bin/env pwsh
#
# Build the Valenx Windows .msi installer locally.
#
# This is the one-command equivalent of the `Windows .msi` job in
# `.github/workflows/release.yml`. Prereqs:
#
#   cargo install cargo-wix --locked
#   # WiX Toolset 3.x on PATH (candle.exe + light.exe). cargo-wix
#   # auto-discovers it from the registry; if you don't have it,
#   # cargo-wix prints a download link.
#
# Output: `target/wix/Valenx-<version>-x86_64.msi`.

$ErrorActionPreference = "Stop"

$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot   = Split-Path -Parent $scriptRoot
Set-Location $repoRoot

Write-Host "[1/2] Building valenx-app release binary..."
cargo build --release -p valenx-app

Write-Host "[2/2] Packaging valenx.exe into MSI via cargo-wix..."
cargo wix -p valenx-app --no-build --nocapture

Write-Host ""
Write-Host "Done. Installer(s):"
Get-ChildItem (Join-Path $repoRoot "target/wix/*.msi") | ForEach-Object {
    $sizeKb = [Math]::Round($_.Length / 1024, 1)
    Write-Host "  $($_.FullName)  ($sizeKb KB)"
}

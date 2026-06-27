<#
.SYNOPSIS
  Relaunch valenx and capture its window to a PNG — the one-call version of the
  kill -> build -> launch -> (open product) -> (tile) -> (maximize) -> capture
  scaffolding, so AI/dev demos are a single command instead of inline PowerShell.

.EXAMPLE
  scripts/valenx-shot.ps1 -Out shot.png
  scripts/valenx-shot.ps1 -Product thermo -Maximize
  scripts/valenx-shot.ps1 -Product fem -Tile -Maximize -Out grid.png
  scripts/valenx-shot.ps1 -NoBuild -Out quick.png

.NOTES
  Companion to scripts/valenx-drive.ps1 (bridge commands) and the headless
  `valenx --self-test` verification harness. Bridge env + global cmd/feed files
  match those tools. Window capture is a screen-region copy of the live valenx
  window (GPU content renders fine; PrintWindow does not for wgpu).
#>
param(
  [string]$Product = "",      # optional workbench id to open via the bridge (e.g. thermo, fem, rocket)
  [string]$Out = "$env:TEMP\valenx_shot.png",
  [switch]$Tile,              # send {"cmd":"tile"} to grid the open panels before capture
  [switch]$Maximize,          # maximize the window before capture
  [switch]$Show3d,            # open the 3-D viewport (show the designed model) before capture
  [switch]$NoBuild,           # skip the cargo rebuild (use the current binary)
  [int]$WaitSec = 8           # GPU-init wait after launch
)
$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$base = $env:TEMP
$env:VALENX_ASSISTANT_INBOX = "$base/valenx_chat_inbox.jsonl"
$env:VALENX_ASSISTANT_FEED  = "$base/valenx_chat_feed.jsonl"
$cmdFile = "$base/valenx_chat_cmd.jsonl"
$enc = New-Object System.Text.UTF8Encoding($false)
function Send-Cmd($json) { [System.IO.File]::AppendAllText($cmdFile, $json + "`n", $enc) }

Get-Process valenx -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Milliseconds 1000

if (-not $NoBuild) {
  Push-Location $repo
  cargo build -p valenx-app --bin valenx
  $code = $LASTEXITCODE
  Pop-Location
  if ($code -ne 0) { Write-Output "BUILD FAILED=$code"; exit 1 }
}

Start-Process -FilePath "$repo\target\debug\valenx.exe"
Start-Sleep -Seconds $WaitSec

Add-Type -AssemblyName System.Drawing
Add-Type @"
using System;using System.Runtime.InteropServices;
public class VShot{
[DllImport("user32.dll")]public static extern bool ShowWindow(IntPtr h,int n);
[DllImport("user32.dll")]public static extern bool GetWindowRect(IntPtr h,out RECT r);
[DllImport("user32.dll")]public static extern bool SetForegroundWindow(IntPtr h);
public struct RECT{public int Left,Top,Right,Bottom;}}
"@
$p = Get-Process valenx -ErrorAction SilentlyContinue | Where-Object { $_.MainWindowHandle -ne 0 } | Select-Object -First 1
if ($null -eq $p) { Write-Output "NO WINDOW"; exit 1 }
$h = $p.MainWindowHandle

if ($Product -ne "") {
  Send-Cmd ('{"cmd":"new_tab","name":"' + $Product + '","workbench":"' + $Product + '"}')
  Start-Sleep -Seconds 3
}
if ($Tile)     { Send-Cmd '{"cmd":"tile"}'; Start-Sleep -Seconds 2 }
if ($Show3d)   { Send-Cmd '{"cmd":"show_3d"}'; Start-Sleep -Seconds 2 }   # render the designed model, not just panels
if ($Maximize) { [VShot]::ShowWindow($h, 3) | Out-Null }
[VShot]::SetForegroundWindow($h) | Out-Null
Start-Sleep -Milliseconds 1000

$r = New-Object VShot+RECT
[VShot]::GetWindowRect($h, [ref]$r) | Out-Null
$w = $r.Right - $r.Left; $ht = $r.Bottom - $r.Top
$bmp = New-Object System.Drawing.Bitmap $w, $ht
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.CopyFromScreen($r.Left, $r.Top, 0, 0, $bmp.Size)
$bmp.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose(); $bmp.Dispose()
Write-Output "SHOT $Out  ${w}x${ht}"

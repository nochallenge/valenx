<#
.SYNOPSIS
    Drive a running valenx instance over its AI file-bridge: append ONE command
    to the global command channel, then print the new tail of the global feed.

.DESCRIPTION
    valenx polls a global newline-delimited-JSON command file each frame and
    applies every newly-appended line through the same vetted code path a user
    click would (see crates/valenx-app/src/agent_commands.rs and
    docs/AI-DRIVING.md). This helper is the one-liner front-end to that channel
    so an external agent never has to hand-roll BOM-free UTF-8 appends or guess
    where the files live.

    It:
      (a) resolves the global command + feed files the SAME way valenx does
          (honouring VALENX_ASSISTANT_INBOX / VALENX_ASSISTANT_FEED, else the
          per-OS state dir %APPDATA%\valenx),
      (b) appends your command as a single BOM-FREE UTF-8 line
          ([System.IO.File]::AppendAllText with a no-BOM UTF8Encoding -- egui's
          JSON line parser will choke on a UTF-16/BOM-prefixed line),
      (c) sleeps briefly so the ~1 Hz poll can pick it up, then
      (d) prints the new lines that appeared in the global feed
          (valenx_chat_feed.jsonl) -- the acks / readouts valenx posts back.

    The global channel honours the FULL drive vocabulary (tab ops, set_control,
    run_command, read_readout, list_controls, list_commands, camera, sketch,
    2-D, note, animate, plus the new_unit bootstrap) against the ACTIVE tab --
    no per-unit bootstrap needed.

.PARAMETER Command
    Shorthand command name, e.g. set_control / open_workbench / run_command /
    read_readout / list_controls / list_commands / note / new_tab / focus_tab /
    rename_tab / close_tab / set_view / frame_all / new_unit. The remaining
    positional Rest values are mapped to that command's fields (see EXAMPLES).
    Numbers and booleans in Rest are emitted as JSON numbers/bools; everything
    else as a JSON string.

.PARAMETER Rest
    Positional arguments for the shorthand command (field values, in order).

.PARAMETER Raw
    A complete JSON command object to append verbatim, e.g.
    -Raw '{"cmd":"open_workbench","id":"thermo"}'. Bypasses the shorthand
    mapping -- use this for any command/field the shorthand doesn't cover.

.PARAMETER Wait
    Seconds to sleep before reading the feed tail (default 1.5; the poll runs
    at ~1 Hz, so allow at least ~1 s for the ack to appear).

.PARAMETER Tail
    How many trailing feed lines to print after the command (default 8).

.EXAMPLE
    ./valenx-drive.ps1 set_control "Temperature [K]" 350
    # appends {"cmd":"set_control","name":"Temperature [K]","value":350}

.EXAMPLE
    ./valenx-drive.ps1 open_workbench thermo
    # appends {"cmd":"open_workbench","id":"thermo"}

.EXAMPLE
    ./valenx-drive.ps1 -Raw '{"cmd":"open_workbench","id":"thermo"}'

.EXAMPLE
    ./valenx-drive.ps1 read_readout thermo      # ask a workbench for its result
    ./valenx-drive.ps1 list_controls thermo     # discover settable captions
    ./valenx-drive.ps1 list_commands            # discover runnable command ids
#>
[CmdletBinding(DefaultParameterSetName = 'Shorthand')]
param(
    [Parameter(ParameterSetName = 'Shorthand', Position = 0)]
    [string]$Command,

    [Parameter(ParameterSetName = 'Shorthand', Position = 1, ValueFromRemainingArguments = $true)]
    [string[]]$Rest,

    [Parameter(ParameterSetName = 'Raw', Mandatory = $true)]
    [string]$Raw,

    [double]$Wait = 1.5,
    [int]$Tail = 8
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# ---- Resolve the global command + feed files exactly as valenx does ----------
# Command file: <base-dir>/valenx_chat_cmd.jsonl, where <base-dir> is the parent
# of VALENX_ASSISTANT_INBOX (the chat-channel dir) or the state dir.
function Get-ValenxStateDir {
    if ($env:APPDATA) { return (Join-Path $env:APPDATA 'valenx') }
    return (Join-Path $env:TEMP 'valenx')
}

function Get-CmdBaseDir {
    if ($env:VALENX_ASSISTANT_INBOX) {
        $d = Split-Path -Parent $env:VALENX_ASSISTANT_INBOX
        if ($d) { return $d }
    }
    return (Get-ValenxStateDir)
}

function Get-GlobalFeedPath {
    if ($env:VALENX_ASSISTANT_FEED) { return $env:VALENX_ASSISTANT_FEED }
    return (Join-Path (Get-ValenxStateDir) 'assistant_feed.jsonl')
}

# Quote-and-JSON-escape a string field. Compress so it is one line.
function Quote([string]$s) { return (ConvertTo-Json $s -Compress) }

# A bare token that looks like a number/bool becomes a JSON number/bool; quoted
# or non-numeric tokens stay strings. Mirrors AgentValue's untagged decoding.
function ConvertTo-JsonScalar([string]$s) {
    $n = 0.0
    if ([double]::TryParse($s, [ref]$n)) { return $s }
    if ($s -ieq 'true') { return 'true' }
    if ($s -ieq 'false') { return 'false' }
    return (Quote $s)
}

$baseDir = Get-CmdBaseDir
$cmdPath = Join-Path $baseDir 'valenx_chat_cmd.jsonl'
$feedPath = Get-GlobalFeedPath

# ---- Build the JSON command line ---------------------------------------------
if ($PSCmdlet.ParameterSetName -eq 'Raw') {
    $line = $Raw.Trim()
}
else {
    if (-not $Command) {
        throw "Provide a shorthand command (e.g. set_control ...) or -Raw '<json>'. See -? for help."
    }
    if (-not $Rest) { $Rest = @() }

    # Map the common shorthands to their {cmd,...} object. Anything not listed
    # here should be sent with -Raw.
    switch -Regex ($Command) {
        '^set_control$' {
            if ($Rest.Count -lt 2) { throw 'set_control <caption> <value> [workbench]' }
            $line = '{"cmd":"set_control","name":' + (Quote $Rest[0]) + ',"value":' + (ConvertTo-JsonScalar $Rest[1])
            if ($Rest.Count -ge 3) { $line += ',"workbench":' + (Quote $Rest[2]) }
            $line += '}'
        }
        '^open_workbench$' {
            if ($Rest.Count -lt 1) { throw 'open_workbench <id>' }
            $line = '{"cmd":"open_workbench","id":' + (Quote $Rest[0]) + '}'
        }
        '^run_command$' {
            if ($Rest.Count -lt 1) { throw 'run_command <id>' }
            $line = '{"cmd":"run_command","id":' + (Quote $Rest[0]) + '}'
        }
        '^read_readout$' {
            $line = '{"cmd":"read_readout"'
            if ($Rest.Count -ge 1) { $line += ',"workbench":' + (Quote $Rest[0]) }
            $line += '}'
        }
        '^list_controls$' {
            $line = '{"cmd":"list_controls"'
            if ($Rest.Count -ge 1) { $line += ',"workbench":' + (Quote $Rest[0]) }
            $line += '}'
        }
        '^list_commands$' { $line = '{"cmd":"list_commands"}' }
        '^note$' {
            if ($Rest.Count -lt 1) { throw 'note <text> [kind]' }
            $line = '{"cmd":"note","text":' + (Quote $Rest[0])
            if ($Rest.Count -ge 2) { $line += ',"kind":' + (Quote $Rest[1]) }
            $line += '}'
        }
        '^new_tab$' {
            if ($Rest.Count -lt 1) { throw 'new_tab <name> [workbench]' }
            $line = '{"cmd":"new_tab","name":' + (Quote $Rest[0])
            if ($Rest.Count -ge 2) { $line += ',"workbench":' + (Quote $Rest[1]) }
            $line += '}'
        }
        '^focus_tab$' {
            if ($Rest.Count -lt 1) { throw 'focus_tab <name>' }
            $line = '{"cmd":"focus_tab","name":' + (Quote $Rest[0]) + '}'
        }
        '^rename_tab$' {
            if ($Rest.Count -lt 1) { throw 'rename_tab <name>' }
            $line = '{"cmd":"rename_tab","name":' + (Quote $Rest[0]) + '}'
        }
        '^close_tab$' {
            $line = '{"cmd":"close_tab"'
            if ($Rest.Count -ge 1) { $line += ',"name":' + (Quote $Rest[0]) }
            $line += '}'
        }
        '^set_view$' {
            if ($Rest.Count -lt 1) { throw 'set_view <front|back|left|right|top|bottom|iso>' }
            $line = '{"cmd":"set_view","dir":' + (Quote $Rest[0]) + '}'
        }
        '^frame_all$' { $line = '{"cmd":"frame_all"}' }
        '^new_unit$' {
            $line = '{"cmd":"new_unit"'
            if ($Rest.Count -ge 1) { $line += ',"kind":' + (Quote $Rest[0]) }
            if ($Rest.Count -ge 2) { $line += ',"title":' + (Quote $Rest[1]) }
            $line += '}'
        }
        default {
            throw "Unknown shorthand '$Command'. Send it with -Raw '<json>' (see docs/AI-DRIVING.md for the full vocabulary)."
        }
    }
}

# Validate the JSON before writing so a malformed line never lands on the channel.
try { $null = $line | ConvertFrom-Json } catch { throw "Built an invalid JSON command: $line`n$_" }

# ---- Snapshot the feed length, append BOM-FREE UTF-8, then print the tail -----
$feedBefore = 0
if (Test-Path -LiteralPath $feedPath) {
    $feedBefore = @(Get-Content -LiteralPath $feedPath -ErrorAction SilentlyContinue).Count
}

# Ensure the command dir exists, then append one BOM-free UTF-8 line.
$dir = Split-Path -Parent $cmdPath
if ($dir -and -not (Test-Path -LiteralPath $dir)) {
    New-Item -ItemType Directory -Path $dir -Force | Out-Null
}
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::AppendAllText($cmdPath, $line + "`n", $utf8NoBom)

Write-Host ">> $line" -ForegroundColor Cyan
Write-Host "   cmd : $cmdPath" -ForegroundColor DarkGray
Write-Host "   feed: $feedPath" -ForegroundColor DarkGray

Start-Sleep -Seconds $Wait

if (Test-Path -LiteralPath $feedPath) {
    $all = @(Get-Content -LiteralPath $feedPath -ErrorAction SilentlyContinue)
    if ($all.Count -gt $feedBefore) {
        $new = $all[$feedBefore..($all.Count - 1)]
    }
    else {
        $new = @()
    }
    if ($new.Count -eq 0) {
        # Nothing new yet -- print the existing tail so the caller still sees state.
        $start = [Math]::Max(0, $all.Count - $Tail)
        if ($all.Count -gt 0) { $new = $all[$start..($all.Count - 1)] } else { $new = @() }
        Write-Host "(no new feed lines yet -- showing last $Tail)" -ForegroundColor Yellow
    }
    foreach ($entry in $new) {
        try {
            $o = $entry | ConvertFrom-Json
            $kind = if ($o.PSObject.Properties.Name -contains 'kind') { $o.kind } else { '' }
            $title = if ($o.PSObject.Properties.Name -contains 'title') { $o.title } else { '' }
            $detail = if ($o.PSObject.Properties.Name -contains 'detail') { $o.detail } else { '' }
            $color = switch ($kind) {
                'warn' { 'Red' }
                'result' { 'Green' }
                'ship' { 'Cyan' }
                default { 'Gray' }
            }
            Write-Host ("[{0}] {1}: {2}" -f $kind, $title, $detail) -ForegroundColor $color
        }
        catch {
            Write-Host $entry
        }
    }
}
else {
    Write-Host "(global feed not found at $feedPath -- is valenx running with the bridge active?)" -ForegroundColor Yellow
}

<#
.SYNOPSIS
    Screen EVERY valenx product/workbench over the agent file-bridge and write a
    pass/fail report.

.DESCRIPTION
    Drives a LIVE, already-running valenx instance through the GLOBAL command
    bridge documented in docs/AI-DRIVING.md + crates/valenx-app/src/agent_commands.rs.
    For each of the 56 products in docs/PRODUCTS.md it appends a sequence of
    BOM-free JSON-line commands to `<base>/valenx_chat_cmd.jsonl`. There are TWO
    drive paths:

      MAPPED (the 13 products with a registered bridge run-id):
        new_tab{name,workbench}  ->  open the product on a fresh tab
        list_controls{workbench} ->  enumerate its settable control captions
        run_command{id}          ->  fire the solve via its bridge run-id
        read_readout{workbench}  ->  read the computed result back
        close_tab{name}          ->  close the product tab

      GENERIC (the ~43 products with NO mapped run-id — uses the new generic
      accessibility-name bridge commands so ANY workbench can be run + read with
      no per-workbench wiring; see crates/valenx-app/src/agent_commands.rs
      InvokeNamed / ListButtons / ReadText):
        new_tab{name,workbench}  ->  open the product on a fresh tab
        list_controls{workbench} ->  enumerate its settable control captions
        list_buttons{}           ->  enumerate the active panel's clickable
                                     button captions
        invoke_named{name}       ->  click the PRIMARY ACTION button (the one
                                     whose caption matches a run-verb, in the
                                     priority order documented at $RunVerbs)
        read_text{}              ->  read the active panel's visible text back
        close_tab{name}          ->  close the product tab

    and harvests the NEW lines valenx appends to the GLOBAL feed
    (`<base>/valenx_chat_feed.jsonl`). Every command's ack/warn/result lands in
    the global feed (channel 0) per `apply_global`, so no per-unit bootstrap is
    needed. The harness classifies each product PASS / PARTIAL / FAIL and writes
    a markdown report (noting the action button used for the generic path).

    This script ONLY drives the bridge — it does NOT launch valenx. Start valenx
    first (ideally with $env:VALENX_ASSISTANT_INBOX / $env:VALENX_ASSISTANT_FEED
    pointing at a known temp location) and run this against it.

.PARAMETER Wait
    Milliseconds to wait between commands (default 500). The bridge polls ~once
    per second, so 500ms+ per command gives valenx time to apply and ack each.

.PARAMETER Base
    Bridge base directory (where valenx_chat_cmd.jsonl / valenx_chat_feed.jsonl
    live). Defaults to the directory of $env:VALENX_ASSISTANT_INBOX if set, else
    $env:TEMP — the same resolution valenx uses.

.EXAMPLE
    ./scripts/valenx-screen.ps1
    ./scripts/valenx-screen.ps1 -Wait 800
    ./scripts/valenx-screen.ps1 -Base "C:/Users/me/AppData/Local/Temp"
#>
[CmdletBinding()]
param(
    [int]$Wait = 500,
    [string]$Base
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# ---------------------------------------------------------------------------
# 0. Resolve the bridge base directory the SAME way valenx does
#    (docs/AI-DRIVING.md ss.1): the dir of $VALENX_ASSISTANT_INBOX if set,
#    otherwise $env:TEMP. The global command/feed files live directly in it.
# ---------------------------------------------------------------------------
if ([string]::IsNullOrWhiteSpace($Base)) {
    if (-not [string]::IsNullOrWhiteSpace($env:VALENX_ASSISTANT_INBOX)) {
        $Base = Split-Path -Parent $env:VALENX_ASSISTANT_INBOX
    } else {
        $Base = $env:TEMP
    }
}
if ([string]::IsNullOrWhiteSpace($Base) -or -not (Test-Path -LiteralPath $Base)) {
    Write-Error "Bridge base directory '$Base' does not exist. Pass -Base <dir> pointing at the folder holding valenx_chat_cmd.jsonl."
    exit 1
}

$CmdFile  = Join-Path $Base 'valenx_chat_cmd.jsonl'
# Global feed path: $VALENX_ASSISTANT_FEED if set, else <base>/valenx_chat_feed.jsonl.
if (-not [string]::IsNullOrWhiteSpace($env:VALENX_ASSISTANT_FEED)) {
    $FeedFile = $env:VALENX_ASSISTANT_FEED
} else {
    $FeedFile = Join-Path $Base 'valenx_chat_feed.jsonl'
}
$ReportFile = Join-Path $env:TEMP 'valenx_screen_report.md'

# ---------------------------------------------------------------------------
# 1. Verify valenx is actually running. No process -> nothing to drive.
# ---------------------------------------------------------------------------
$proc = $null
try { $proc = Get-Process -Name 'valenx' -ErrorAction SilentlyContinue } catch { $proc = $null }
if (-not $proc) {
    Write-Error "valenx is not running (Get-Process valenx found nothing). Launch valenx first, then re-run this harness against it."
    exit 1
}

Write-Host "valenx-screen: driving live valenx (pid $($proc[0].Id))"
Write-Host "  bridge cmd : $CmdFile"
Write-Host "  bridge feed: $FeedFile"
Write-Host "  report     : $ReportFile"
Write-Host ""

# ---------------------------------------------------------------------------
# 2. Bridge helpers.
#    - Append a command as ONE BOM-free UTF-8 JSON line (docs ss.1 'Write rules').
#    - Read the feed line count + the feed entries (tolerant of partial lines).
# ---------------------------------------------------------------------------
$Utf8NoBom = New-Object System.Text.UTF8Encoding($false)

function Send-Cmd {
    param([hashtable]$Obj)
    # ConvertTo-Json with -Compress emits a single line; depth covers nested
    # arrays (none here, but safe). Append exactly one '`n'-terminated line.
    $json = ($Obj | ConvertTo-Json -Compress -Depth 6)
    [System.IO.File]::AppendAllText($CmdFile, $json + "`n", $Utf8NoBom)
    Start-Sleep -Milliseconds $Wait
}

function Get-FeedLines {
    # Return the feed as an array of raw, non-empty trimmed lines. A missing or
    # transiently-locked feed yields an empty array (never throws).
    if (-not (Test-Path -LiteralPath $FeedFile)) { return @() }
    try {
        $raw = [System.IO.File]::ReadAllText($FeedFile, $Utf8NoBom)
    } catch {
        try { $raw = Get-Content -LiteralPath $FeedFile -Raw -ErrorAction Stop } catch { return @() }
    }
    if ([string]::IsNullOrEmpty($raw)) { return @() }
    return @($raw -split "`r?`n" | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
}

function Get-FeedCount {
    return @(Get-FeedLines).Count
}

function ConvertFrom-FeedLine {
    # Parse one feed JSON line into {title,detail,kind}; $null on garbage.
    param([string]$Line)
    try {
        $o = $Line | ConvertFrom-Json -ErrorAction Stop
        return [pscustomobject]@{
            title  = ([string]$o.title)
            detail = ([string]$o.detail)
            kind   = ([string]$o.kind)
        }
    } catch {
        return $null
    }
}

# ---------------------------------------------------------------------------
# 2b. Open-settle + ack-polling helpers.
#
#     Heavy workbenches (fem / cfd / genetics / uq / reactdyn ...) take longer
#     than $Wait to finish loading; the FIRST screen FALSE-FAILed several of
#     them purely because `list_controls` fired before the panel existed. So
#     after opening a tab we POLL the feed for the matching command's ack rather
#     than guessing a fixed sleep: append the command, then watch for the new
#     feed line whose detail matches an expected pattern, up to a timeout.
#
#     Wait-ForAck returns $true once a NEW feed line (appended after $FromCount)
#     matches $Pattern, else $false on timeout. It never throws.
# ---------------------------------------------------------------------------
function Wait-ForAck {
    param(
        [int]$FromCount,           # feed line count BEFORE the command was sent
        [string]$Pattern,          # regex the awaited ack's detail must match
        [int]$TimeoutMs = 2500,    # give heavy workbenches time to load
        [int]$PollMs = 150
    )
    $deadline = [DateTime]::UtcNow.AddMilliseconds($TimeoutMs)
    while ([DateTime]::UtcNow -lt $deadline) {
        $lines = @(Get-FeedLines)
        if ($lines.Count -gt $FromCount) {
            for ($i = $FromCount; $i -lt $lines.Count; $i++) {
                $e = ConvertFrom-FeedLine $lines[$i]
                if ($null -ne $e -and $e.detail -match $Pattern) { return $true }
            }
        }
        Start-Sleep -Milliseconds $PollMs
    }
    return $false
}

# Run-verb priority list: when a product has no mapped run-id we pick its
# PRIMARY ACTION button from the list_buttons captions by matching these verbs
# (case-insensitive substring), in THIS priority order. The first verb that
# matches any caption wins; among captions matching the same verb the first as
# returned by list_buttons wins. Pure-view buttons (below) are skipped first so
# e.g. "Show 3-D…" never out-ranks a real action.
$RunVerbs = @(
    'Analyze', 'Compute', 'Run', 'Solve', 'Calculate', 'Simulate',
    'Build', 'Generate', 'Eval', 'Train', 'Route', 'Play', 'Apply'
)
# Pure-view / non-action captions to exclude from action-button selection.
$ViewVerbs = @('Show', 'Preview', 'Reset', 'View', 'Frame', 'Zoom', 'Orbit',
    'Export', 'Save', 'Load', 'Open', 'Close', 'Clear', 'Copy', 'Help', 'About')

function Select-ActionButton {
    # From the list_buttons captions of one product, return the caption of the
    # primary action button (verbatim, INCLUDING any leading "▶ " U+25B6), or
    # $null if none matches a run-verb. View-only captions are excluded first.
    param([string[]]$Captions)
    $caps = @($Captions | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    if ($caps.Count -eq 0) { return $null }
    # Drop pure-view buttons so they never win.
    $actionable = @($caps | Where-Object {
        $c = $_
        $isView = $false
        foreach ($v in $ViewVerbs) {
            if ($c -match [regex]::Escape($v)) { $isView = $true; break }
        }
        -not $isView
    })
    if ($actionable.Count -eq 0) { $actionable = $caps }
    foreach ($verb in $RunVerbs) {
        foreach ($c in $actionable) {
            if ($c -match [regex]::Escape($verb)) { return $c }
        }
    }
    return $null
}

function Get-ButtonCaptions {
    # Pull the captions list out of a "buttons (N): a, b, c" feed result line.
    # Returns @() when none / the "(none …)" form.
    param([object[]]$Entries)
    foreach ($e in $Entries) {
        if ($null -eq $e) { continue }
        if ($e.detail -match '^buttons \(\d+\):\s*(.*)$') {
            $rest = $Matches[1].Trim()
            if ($rest -eq '' -or $rest -like '(none*') { return @() }
            return @($rest -split ',\s*' | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
        }
    }
    return @()
}

# ---------------------------------------------------------------------------
# 3. The product list (docs/PRODUCTS.md 'id' column, all 56) with group, and
#    the best-effort id -> RunCommand-id MAP.
#
#    The run-ids are sourced from the bridge run-dispatch in
#    crates/valenx-app/src/agent_commands.rs (AgentCommand::RunCommand): a
#    special-case block routes a set of BRIDGE-ONLY ids (kept out of the
#    user-facing palette) to each workbench's run path:
#        thermo.compute quantum.run optics.compute acoustics.compute
#        waveform.parse topopt.run nodegraph.eval bondgraph.solve
#        surrogate.train brep.build missionsim.run missionplanner.route
#        morphogenesis.play
#    Products with no registered bridge run-id are recorded "no run-id" (the
#    harness still opens them, lists controls, and tries a readout).
#
#    NOTE on run-id gating: every bridge run helper first checks the run lands
#    on the matching ACTIVE tab; since new_tab makes the product tab active,
#    the mapped run fires against the right workbench.
# ---------------------------------------------------------------------------
$Products = @(
    # id                group              run-id
    @('rocket',         'Aerospace',       $null),
    @('engine',         'Aerospace',       $null),
    @('astro',          'Aerospace',       $null),
    @('aero',           'Aerospace',       $null),
    @('gasdynamics',    'Aerospace',       $null),
    @('rotor',          'Aerospace',       $null),
    @('uas',            'Aerospace',       $null),

    @('blackhole',      'Astrophysics',    $null),

    @('cfd',            'Simulation',      $null),
    @('fem',            'Simulation',      $null),
    @('topopt',         'Simulation',      'topopt.run'),
    @('nodegraph',      'Simulation',      'nodegraph.eval'),
    @('bondgraph',      'Simulation',      'bondgraph.solve'),
    @('surrogate',      'Simulation',      'surrogate.train'),
    @('reactdyn',       'Simulation',      $null),
    @('thermo',         'Simulation',      'thermo.compute'),
    @('quantum',        'Simulation',      'quantum.run'),
    @('fields',         'Simulation',      $null),
    @('fluids',         'Simulation',      $null),
    @('ocean',          'Simulation',      $null),
    @('rom',            'Simulation',      $null),
    @('uq',             'Simulation',      $null),
    @('missionsim',     'Simulation',      'missionsim.run'),
    @('missionplanner', 'Simulation',      'missionplanner.route'),
    @('survivability',  'Simulation',      $null),
    @('cosim',          'Simulation',      $null),
    @('mbd',            'Simulation',      $null),
    @('optics',         'Simulation',      'optics.compute'),
    @('acoustics',      'Simulation',      'acoustics.compute'),
    @('waveform',       'Simulation',      'waveform.parse'),

    @('cad',            'CAD & mesh',      $null),
    @('brep',           'CAD & mesh',      'brep.build'),
    @('mesh',           'CAD & mesh',      $null),
    @('sheetmetal',     'CAD & mesh',      $null),
    @('reverse',        'CAD & mesh',      $null),
    @('photogrammetry', 'CAD & mesh',      $null),
    @('draft2d',        'CAD & mesh',      $null),
    @('render',         'CAD & mesh',      $null),
    @('animate',        'CAD & mesh',      $null),

    @('springs',        'Machine design',  $null),
    @('gears',          'Machine design',  $null),
    @('fasteners',      'Machine design',  $null),
    @('frames',         'Machine design',  $null),
    @('collision',      'Machine design',  $null),

    @('piping',         'Civil & AEC',     $null),
    @('hvac',           'Civil & AEC',     $null),
    @('reinforcement',  'Civil & AEC',     $null),
    @('interior',       'Civil & AEC',     $null),
    @('geomatics',      'Civil & AEC',     $null),

    @('genetics',       'Life sciences',   $null),
    @('neuro',          'Life sciences',   $null),
    @('variant',        'Life sciences',   $null),
    @('ppi',            'Life sciences',   $null),
    @('morphogenesis',  'Life sciences',   'morphogenesis.play'),

    @('sensors',        'Sensors',         $null),
    @('autonomy',       'Sensors',         $null)
)

$runIdCount = @($Products | Where-Object { $null -ne $_[2] }).Count
$genericCount = $Products.Count - $runIdCount
Write-Host ("Screening {0} products; {1} via mapped run-id, {2} via the generic list_buttons/invoke_named/read_text path." -f `
    $Products.Count, $runIdCount, $genericCount)
Write-Host ""

# ---------------------------------------------------------------------------
# 4. Classification of the NEW feed lines harvested for one product.
#
#    Feed lines are {title,detail,kind} (assistant_workbench::append_feed_note).
#    Result lines carry kind="result"; failures carry kind="warn".
#
#    MAPPED path (a run-id fired) — same as before:
#      PASS    : a non-error 'result' readout/run-result line
#                ("<wb> readout:" or "ran <id>") is present.
#      PARTIAL : opened + has controls but no readout/run-result line.
#      FAIL    : no new feed lines, an error ack, or valenx died.
#
#    GENERIC path (no run-id; list_buttons -> invoke_named -> read_text):
#      PASS    : invoke_named queued a click AND read_text returned substantive
#                non-error content (a "text: …" result with real characters).
#      PARTIAL : opened + has controls/buttons but either no action button
#                matched a run-verb, or read_text came back empty / "(no
#                readable text…)".
#      FAIL    : no new feed lines at all, an invoke_named/command error ack, or
#                valenx died.
# ---------------------------------------------------------------------------
function Get-ControlCount {
    param([object[]]$Entries)
    # Find the "controls (N): ..." result line and pull N. -1 if none seen.
    foreach ($e in $Entries) {
        if ($null -eq $e) { continue }
        if ($e.detail -match 'controls \((\d+)\):') { return [int]$Matches[1] }
        if ($e.detail -match 'no settable controls for workbench') { return 0 }
    }
    return -1
}

function Test-IsReadout {
    param([object]$Entry)
    if ($null -eq $Entry) { return $false }
    if ($Entry.kind -ne 'result') { return $false }
    $d = $Entry.detail
    if ([string]::IsNullOrWhiteSpace($d)) { return $false }
    # A genuine readout/run-result line.
    if ($d -match ' readout: ') { return $true }
    if ($d -match '^ran ')       { return $true }
    return $false
}

function Test-IsSubstantiveText {
    # True if a read_text 'result' line returned real, non-empty panel text
    # (the "text: …" form with substantive content), false for the empty
    # "read_text: (no readable text…)" sentinel or a too-short payload.
    param([object]$Entry)
    if ($null -eq $Entry) { return $false }
    if ($Entry.kind -ne 'result') { return $false }
    $d = $Entry.detail
    if ([string]::IsNullOrWhiteSpace($d)) { return $false }
    if ($d -like '*no readable text*') { return $false }
    if ($d -notmatch '^text:\s*(.+)$') { return $false }
    $body = $Matches[1].Trim()
    # Require some real content beyond a stray bullet/space.
    return ($body.Length -ge 3)
}

function Test-IsQueuedClick {
    # True if invoke_named acked that it queued a click (the generic-path action
    # fired successfully).
    param([object]$Entry)
    if ($null -eq $Entry) { return $false }
    if ($Entry.kind -ne 'result') { return $false }
    return ([string]$Entry.detail -match '^invoke_named: queued click')
}

function Test-IsError {
    param([object]$Entry)
    if ($null -eq $Entry) { return $false }
    $d = $Entry.detail
    if ([string]::IsNullOrWhiteSpace($d)) { return $false }
    if ($Entry.kind -eq 'warn') { return $true }
    # Belt-and-suspenders: catch known error phrases even if kind drifted.
    foreach ($p in @(
        'unknown command id',
        'unknown workbench id',
        'not run yet',
        'no readout',
        'has no readout wired',
        'is not the ')) {
        if ($d -like "*$p*") { return $true }
    }
    return $false
}

# Count buttons enumerated by a list_buttons "buttons (N): …" feed line (-1 if
# the line was never seen).
function Get-ButtonCount {
    param([object[]]$Entries)
    foreach ($e in $Entries) {
        if ($null -eq $e) { continue }
        if ($e.detail -match '^buttons \((\d+)\):') { return [int]$Matches[1] }
    }
    return -1
}

# Verdict for the MAPPED path (a bridge run-id fired): PASS on a readout/run
# result, else PARTIAL if it opened with controls, else FAIL.
function Get-VerdictMapped {
    param([object[]]$Entries)
    $entries = @($Entries)
    if ($entries.Count -eq 0) {
        return @{ Verdict = 'FAIL'; Snippet = '(no feed response)' }
    }

    $hasReadout = $false
    foreach ($e in $entries) { if (Test-IsReadout $e) { $hasReadout = $true; break } }

    $ctrl = Get-ControlCount -Entries $entries

    # Prefer a readout/run-result line for the snippet; else the last warn; else
    # the last line.
    $snipEntry = $null
    foreach ($e in $entries) { if (Test-IsReadout $e) { $snipEntry = $e } }
    if ($null -eq $snipEntry) {
        for ($i = $entries.Count - 1; $i -ge 0; $i--) {
            if ((Test-IsError $entries[$i])) { $snipEntry = $entries[$i]; break }
        }
    }
    if ($null -eq $snipEntry) { $snipEntry = $entries[$entries.Count - 1] }
    $snippet = if ($null -ne $snipEntry) { [string]$snipEntry.detail } else { '' }

    if ($hasReadout) {
        return @{ Verdict = 'PASS'; Snippet = $snippet }
    }

    # No readout. If we at least opened and saw controls, it's PARTIAL
    # (viewer / not-run / AI-drive gap). Otherwise FAIL.
    if ($ctrl -ge 0) {
        return @{ Verdict = 'PARTIAL'; Snippet = $snippet }
    }

    return @{ Verdict = 'FAIL'; Snippet = $snippet }
}

# Verdict for the GENERIC path (no run-id): list_buttons -> invoke_named ->
# read_text. PASS when a click was queued AND read_text returned substantive
# content; PARTIAL when it opened (controls/buttons) but no action button
# matched or read_text was empty; FAIL on no response / an error ack.
function Get-VerdictGeneric {
    param(
        [object[]]$Entries,
        [string]$ActionButton   # caption we invoked, or '' if none matched
    )
    $entries = @($Entries)
    if ($entries.Count -eq 0) {
        return @{ Verdict = 'FAIL'; Snippet = '(no feed response)' }
    }

    $clicked  = $false
    $textOk   = $false
    foreach ($e in $entries) {
        if (Test-IsQueuedClick $e)     { $clicked = $true }
        if (Test-IsSubstantiveText $e) { $textOk  = $true }
    }
    $ctrl = Get-ControlCount -Entries $entries
    $btn  = Get-ButtonCount  -Entries $entries

    # Snippet: prefer the read_text content, else the last error, else last line.
    $snipEntry = $null
    foreach ($e in $entries) { if (Test-IsSubstantiveText $e) { $snipEntry = $e } }
    if ($null -eq $snipEntry) {
        for ($i = $entries.Count - 1; $i -ge 0; $i--) {
            if ((Test-IsError $entries[$i])) { $snipEntry = $entries[$i]; break }
        }
    }
    if ($null -eq $snipEntry) { $snipEntry = $entries[$entries.Count - 1] }
    $snippet = if ($null -ne $snipEntry) { [string]$snipEntry.detail } else { '' }

    # PASS: the action button fired and the panel read back real text.
    if ($clicked -and $textOk) {
        return @{ Verdict = 'PASS'; Snippet = $snippet }
    }

    # PARTIAL: opened (saw controls or buttons) but either no action button
    # matched a run-verb, or read_text came back empty.
    if ($ctrl -ge 0 -or $btn -ge 0) {
        return @{ Verdict = 'PARTIAL'; Snippet = $snippet }
    }

    return @{ Verdict = 'FAIL'; Snippet = $snippet }
}

# ---------------------------------------------------------------------------
# 5. Drive every product. Robust: a per-product try/catch keeps one bad
#    product from aborting the whole run; a vanished valenx process is detected
#    and recorded as FAIL for the rest.
# ---------------------------------------------------------------------------
$Results = New-Object System.Collections.Generic.List[object]
$idx = 0
foreach ($p in $Products) {
    $idx++
    $id    = [string]$p[0]
    $group = [string]$p[1]
    $runId = $p[2]   # may be $null
    $runIdLabel = if ($null -ne $runId) { [string]$runId } else { 'no run-id' }

    Write-Host ("[{0,2}/{1}] {2,-16} ({3})" -f $idx, $Products.Count, $id, $group) -NoNewline

    # Detect a dead valenx before driving the next product.
    $alive = $false
    try { $alive = [bool](Get-Process -Name 'valenx' -ErrorAction SilentlyContinue) } catch { $alive = $false }
    if (-not $alive) {
        Write-Host "  -> FAIL (valenx process gone)"
        $Results.Add([pscustomobject]@{
            Id = $id; Group = $group; Controls = '-'; RunId = $runIdLabel;
            Action = '-'; Verdict = 'FAIL'; Snippet = 'valenx process exited mid-run'
        })
        continue
    }

    $before = Get-FeedCount
    $controlCount = '-'
    $verdict = 'FAIL'
    $snippet = ''
    $actionBtn = ''   # caption invoked on the generic path ('-' for mapped)

    try {
        # open the product on a fresh tab named after the id, then WAIT for the
        # tab to settle. Heavy workbenches (fem/cfd/genetics/uq…) take longer
        # than $Wait to load; poll for the list_controls ack rather than guess.
        $preList = Get-FeedCount
        Send-Cmd @{ cmd = 'new_tab'; name = $id; workbench = $id }
        # enumerate its settable control captions
        Send-Cmd @{ cmd = 'list_controls'; workbench = $id }
        # Block (up to ~2.5s) until the controls ack lands — the real settle.
        [void](Wait-ForAck -FromCount $preList `
            -Pattern '(controls \(\d+\):|no settable controls for workbench)' `
            -TimeoutMs ([Math]::Max(2500, 2 * $Wait)))

        if ($null -ne $runId) {
            # ---- MAPPED path: fire the registered bridge run-id. ----
            Send-Cmd @{ cmd = 'run_command'; id = $runId }
            Send-Cmd @{ cmd = 'read_readout'; workbench = $id }
            $actionBtn = '-'   # ran by run-id, not a named button
        } else {
            # ---- GENERIC path: list_buttons -> invoke_named -> read_text. ----
            # 1. Enumerate clickable button captions and wait for that ack so we
            #    can read them back BEFORE choosing the action button.
            $preBtns = Get-FeedCount
            Send-Cmd @{ cmd = 'list_buttons' }
            [void](Wait-ForAck -FromCount $preBtns -Pattern '^buttons \(\d+\):' `
                -TimeoutMs ([Math]::Max(2500, 2 * $Wait)))
            $btnLines = @(Get-FeedLines)
            $btnEntries = @()
            if ($btnLines.Count -gt $preBtns) {
                foreach ($l in $btnLines[$preBtns..($btnLines.Count - 1)]) {
                    $e = ConvertFrom-FeedLine $l
                    if ($null -ne $e) { $btnEntries += $e }
                }
            }
            $captions = Get-ButtonCaptions -Entries $btnEntries
            # 2. Pick the primary action button by run-verb priority and click it
            #    by its VERBATIM caption (keeps any leading "▶ " U+25B6).
            $pick = Select-ActionButton -Captions $captions
            if ($null -ne $pick) {
                $actionBtn = [string]$pick
                Send-Cmd @{ cmd = 'invoke_named'; name = $actionBtn }
            } else {
                $actionBtn = '(none matched)'
            }
            # 3. Read the panel text back to self-verify the result.
            Send-Cmd @{ cmd = 'read_text' }
        }

        # Harvest the NEW feed lines for this product (everything appended since
        # 'before'). One extra settle in case the last ack is still flushing.
        Start-Sleep -Milliseconds $Wait
        $allLines = @(Get-FeedLines)
        $newRaw = @()
        if ($allLines.Count -gt $before) {
            $newRaw = $allLines[$before..($allLines.Count - 1)]
        }
        $entries = @()
        foreach ($l in $newRaw) {
            $e = ConvertFrom-FeedLine $l
            if ($null -ne $e) { $entries += $e }
        }

        $controlCount = Get-ControlCount -Entries $entries
        if ($null -ne $runId) {
            $vc = Get-VerdictMapped -Entries $entries
        } else {
            $vc = Get-VerdictGeneric -Entries $entries -ActionButton $actionBtn
        }
        $verdict = $vc.Verdict
        $snippet = $vc.Snippet

        # close the product tab (routes through the confirm modal; harmless if
        # it leaves a confirm open — we never assert it closed).
        Send-Cmd @{ cmd = 'close_tab'; name = $id }
    } catch {
        $verdict = 'FAIL'
        $snippet = "harness error: $($_.Exception.Message)"
    }

    $ctrlDisplay = if ($controlCount -ge 0) { [string]$controlCount } else { 'n/a' }
    if ([string]::IsNullOrEmpty($actionBtn)) { $actionBtn = '-' }
    Write-Host ("  -> {0} (controls {1}; action {2})" -f $verdict, $ctrlDisplay, $actionBtn)

    $Results.Add([pscustomobject]@{
        Id = $id; Group = $group; Controls = $ctrlDisplay; RunId = $runIdLabel;
        Action = $actionBtn; Verdict = $verdict; Snippet = $snippet
    })
}

# ---------------------------------------------------------------------------
# 6. Summary counts.
# ---------------------------------------------------------------------------
$pass    = @($Results | Where-Object { $_.Verdict -eq 'PASS' })
$partial = @($Results | Where-Object { $_.Verdict -eq 'PARTIAL' })
$fail    = @($Results | Where-Object { $_.Verdict -eq 'FAIL' })

# ---------------------------------------------------------------------------
# 7. Write the markdown report.
# ---------------------------------------------------------------------------
function Format-Cell {
    param([string]$Text)
    if ($null -eq $Text) { return '' }
    # One-line cell: collapse newlines, escape pipes, clamp length.
    $t = ($Text -replace '\r?\n', ' ') -replace '\|', '\|'
    if ($t.Length -gt 140) { $t = $t.Substring(0, 137) + '...' }
    return $t
}

$sb = New-Object System.Text.StringBuilder
[void]$sb.AppendLine('# valenx product/workbench screen')
[void]$sb.AppendLine('')
[void]$sb.AppendLine(('Generated: {0}' -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss')))
[void]$sb.AppendLine(('Bridge base: `{0}`' -f $Base))
[void]$sb.AppendLine(('Feed: `{0}`' -f $FeedFile))
[void]$sb.AppendLine(('Products screened: {0} (run-id mapped: {1}, generic drive: {2})' -f $Results.Count, $runIdCount, $genericCount))
[void]$sb.AppendLine('')
[void]$sb.AppendLine('Action column = the in-panel button invoked via the generic')
[void]$sb.AppendLine('`invoke_named` path (`-` for the run-id-mapped products).')
[void]$sb.AppendLine('')
[void]$sb.AppendLine('| id | group | controls# | run-id | action button | verdict | snippet |')
[void]$sb.AppendLine('|---|---|---|---|---|---|---|')
foreach ($r in $Results) {
    [void]$sb.AppendLine(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} |' -f `
        $r.Id, $r.Group, $r.Controls, $r.RunId, (Format-Cell $r.Action), `
        $r.Verdict, (Format-Cell $r.Snippet)))
}
[void]$sb.AppendLine('')
[void]$sb.AppendLine('## Summary')
[void]$sb.AppendLine('')
[void]$sb.AppendLine(('- PASS: {0}' -f $pass.Count))
[void]$sb.AppendLine(('- PARTIAL: {0}' -f $partial.Count))
[void]$sb.AppendLine(('- FAIL: {0}' -f $fail.Count))
[void]$sb.AppendLine(('- TOTAL: {0}' -f $Results.Count))
[void]$sb.AppendLine('')
[void]$sb.AppendLine('### FAIL ids')
[void]$sb.AppendLine('')
if ($fail.Count -gt 0) {
    [void]$sb.AppendLine(('`' + (($fail | ForEach-Object { $_.Id }) -join '`, `') + '`'))
} else {
    [void]$sb.AppendLine('_none_')
}
[void]$sb.AppendLine('')
[void]$sb.AppendLine('### PARTIAL ids')
[void]$sb.AppendLine('')
if ($partial.Count -gt 0) {
    [void]$sb.AppendLine(('`' + (($partial | ForEach-Object { $_.Id }) -join '`, `') + '`'))
} else {
    [void]$sb.AppendLine('_none_')
}

[System.IO.File]::WriteAllText($ReportFile, $sb.ToString(), $Utf8NoBom)

# ---------------------------------------------------------------------------
# 8. Print summary counts to stdout.
# ---------------------------------------------------------------------------
Write-Host ''
Write-Host '================ valenx-screen summary ================'
Write-Host ("PASS    : {0}" -f $pass.Count)
Write-Host ("PARTIAL : {0}" -f $partial.Count)
Write-Host ("FAIL    : {0}" -f $fail.Count)
Write-Host ("TOTAL   : {0}" -f $Results.Count)
Write-Host '-------------------------------------------------------'
if ($fail.Count -gt 0) {
    Write-Host ("FAIL ids   : {0}" -f (($fail | ForEach-Object { $_.Id }) -join ', '))
}
if ($partial.Count -gt 0) {
    Write-Host ("PARTIAL ids: {0}" -f (($partial | ForEach-Object { $_.Id }) -join ', '))
}
Write-Host '-------------------------------------------------------'
Write-Host ("Report written: {0}" -f $ReportFile)

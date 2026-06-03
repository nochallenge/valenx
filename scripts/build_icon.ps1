#!/usr/bin/env pwsh
# Build the multi-resolution Valenx Windows .ico file.
#
# Produces crates/valenx-app/wix/valenx.ico containing PNG frames at
# 16, 32, 48 and 256 pixels — the four sizes Windows Explorer + Start
# Menu + taskbar pick from. The output is committed to the repo so
# subsequent builds don't have to regenerate it.
#
# The art is a stylized "V" in white on the Valenx accent blue
# (#4B9EFF, lifted from `crates/valenx-design-tokens/tokens.json`).
#
# Run manually whenever the design changes:
#
#   pwsh scripts/build_icon.ps1
#
# Requires .NET System.Drawing (ships with Windows PowerShell on every
# Windows host — no extra installs).

$ErrorActionPreference = "Stop"

Add-Type -AssemblyName System.Drawing

$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot   = Split-Path -Parent $scriptRoot
$outDir     = Join-Path $repoRoot "crates/valenx-app/wix"
$outFile    = Join-Path $outDir "valenx.ico"

if (-not (Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}

# Brand: accent.primary from valenx-design-tokens (#4B9EFF).
$bgColor = [System.Drawing.Color]::FromArgb(255, 0x4B, 0x9E, 0xFF)
$fgColor = [System.Drawing.Color]::White

function New-VFrame([int]$size) {
    $bmp = New-Object System.Drawing.Bitmap($size, $size, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
    $g   = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode     = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
    $g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
    $g.PixelOffsetMode   = [System.Drawing.Drawing2D.PixelOffsetMode]::HighQuality

    # Rounded-square background: classic Windows app-icon look.
    $radius = [Math]::Max(2, [int]($size / 6))
    $path = New-Object System.Drawing.Drawing2D.GraphicsPath
    $path.AddArc(0,                0,                $radius * 2, $radius * 2, 180, 90)
    $path.AddArc($size - $radius * 2, 0,             $radius * 2, $radius * 2, 270, 90)
    $path.AddArc($size - $radius * 2, $size - $radius * 2, $radius * 2, $radius * 2, 0,   90)
    $path.AddArc(0,                $size - $radius * 2, $radius * 2, $radius * 2, 90,  90)
    $path.CloseFigure()
    $bgBrush = New-Object System.Drawing.SolidBrush($bgColor)
    $g.FillPath($bgBrush, $path)
    $bgBrush.Dispose()

    # "V" mark — two thick strokes meeting at the bottom-center.
    $pad        = [Math]::Max(1, [int]($size * 0.18))
    $stroke     = [Math]::Max(2, [int]($size * 0.18))
    $apexY      = $size - $pad
    $leftX      = $pad
    $rightX     = $size - $pad
    $topY       = $pad
    $center     = $size / 2.0

    $pen = New-Object System.Drawing.Pen($fgColor, $stroke)
    $pen.StartCap = [System.Drawing.Drawing2D.LineCap]::Round
    $pen.EndCap   = [System.Drawing.Drawing2D.LineCap]::Round
    $pen.LineJoin = [System.Drawing.Drawing2D.LineJoin]::Round
    $g.DrawLine($pen, $leftX,  $topY, $center, $apexY)
    $g.DrawLine($pen, $rightX, $topY, $center, $apexY)
    $pen.Dispose()

    $g.Dispose()
    $path.Dispose()
    return $bmp
}

function Get-PngBytes($bmp) {
    $ms = New-Object System.IO.MemoryStream
    $bmp.Save($ms, [System.Drawing.Imaging.ImageFormat]::Png)
    $bytes = $ms.ToArray()
    $ms.Dispose()
    return ,$bytes
}

# Build the four PNG frames.
$sizes  = @(16, 32, 48, 256)
$frames = @()
foreach ($s in $sizes) {
    $bmp   = New-VFrame -size $s
    $png   = Get-PngBytes -bmp $bmp
    $bmp.Dispose()
    $frames += [PSCustomObject]@{ Size = $s; Png = $png }
}

# ICO container — see https://en.wikipedia.org/wiki/ICO_(file_format).
#   ICONDIR  (6 bytes):  reserved=0, type=1 (icon), count
#   ICONDIRENTRY (16 bytes per image): width, height, colors=0, reserved=0,
#                                       planes=1, bpp=32, byteCount, offset
#   Image data (PNG-encoded for each frame — supported by Vista+; the
#   four sizes we ship are all routinely PNG-packed in modern .ico files).
$count  = $frames.Count
$header = New-Object byte[] (6 + 16 * $count)
$header[0] = 0x00; $header[1] = 0x00              # reserved
$header[2] = 0x01; $header[3] = 0x00              # type = 1 (icon)
$header[4] = [byte]($count -band 0xFF)
$header[5] = [byte](($count -shr 8) -band 0xFF)

# Image data starts immediately after the directory.
$dataOffset = 6 + 16 * $count
$imageBlob  = New-Object System.IO.MemoryStream

for ($i = 0; $i -lt $count; $i++) {
    $f = $frames[$i]
    # Width / height of 0 means 256 in ICO.
    $w = if ($f.Size -eq 256) { 0 } else { $f.Size }
    $h = $w
    $bytes  = $f.Png.Length
    $entry  = 6 + 16 * $i

    $header[$entry + 0]  = [byte]$w
    $header[$entry + 1]  = [byte]$h
    $header[$entry + 2]  = 0            # color count (0 = >=8bpp)
    $header[$entry + 3]  = 0            # reserved
    $header[$entry + 4]  = 1; $header[$entry + 5] = 0     # planes = 1
    $header[$entry + 6]  = 32; $header[$entry + 7] = 0    # bpp = 32
    $header[$entry + 8]  = [byte](($bytes        ) -band 0xFF)
    $header[$entry + 9]  = [byte](($bytes -shr 8 ) -band 0xFF)
    $header[$entry + 10] = [byte](($bytes -shr 16) -band 0xFF)
    $header[$entry + 11] = [byte](($bytes -shr 24) -band 0xFF)
    $header[$entry + 12] = [byte](($dataOffset        ) -band 0xFF)
    $header[$entry + 13] = [byte](($dataOffset -shr 8 ) -band 0xFF)
    $header[$entry + 14] = [byte](($dataOffset -shr 16) -band 0xFF)
    $header[$entry + 15] = [byte](($dataOffset -shr 24) -band 0xFF)

    $imageBlob.Write($f.Png, 0, $f.Png.Length)
    $dataOffset += $bytes
}

$fs = [System.IO.File]::Create($outFile)
$fs.Write($header, 0, $header.Length)
$ib = $imageBlob.ToArray()
$fs.Write($ib, 0, $ib.Length)
$fs.Dispose()
$imageBlob.Dispose()

$info = Get-Item $outFile
Write-Host "Wrote $($info.FullName) ($($info.Length) bytes, $count frames: $($sizes -join ', '))"

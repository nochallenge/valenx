#!/usr/bin/env python3
"""
Build the multi-resolution macOS .icns + the Linux 256x256 PNG fallback
for the Valenx app icon.

The brand mark mirrors crates/valenx-app/wix/valenx.ico and
packaging/linux/valenx.svg: a white "V" on the Valenx accent blue
(#4B9EFF, lifted from crates/valenx-design-tokens/tokens.json),
rendered onto a rounded square whose corner radius = size/6 and whose
padding + stroke width are 18% of the viewport.

Outputs (always regenerated, idempotent):
  packaging/linux/valenx-256.png   — 256x256 RGBA PNG, fallback for
                                      Linux launchers that don't render
                                      SVG (themed at /usr/share/icons/
                                      hicolor/256x256/apps/valenx.png).
  packaging/macos/valenx.icns      — Apple Icon Image; multi-resolution
                                      pack consumed by Finder, the Dock,
                                      Launchpad, and Cmd-Tab. Picked up
                                      by cargo-bundle via
                                      [package.metadata.bundle].icon.

Usage:
  python scripts/build_icns.py

Requires:
  Python 3.8+
  Pillow (any modern version)

The script is permanent — re-run it whenever the brand colour or the
"V" geometry changes. Commit the resulting valenx.icns + valenx-256.png
into git so CI doesn't have to regenerate (and so signing / notarisation
operate on stable bytes).
"""

from __future__ import annotations

import struct
import sys
from io import BytesIO
from pathlib import Path

from PIL import Image, ImageDraw

# ---------------------------------------------------------------------------
# Geometry — kept identical to scripts/build_icon.ps1 so all three
# platforms ship the same brand mark.
# ---------------------------------------------------------------------------

BRAND_BG = (0x4B, 0x9E, 0xFF, 0xFF)    # accent.primary
BRAND_FG = (0xFF, 0xFF, 0xFF, 0xFF)    # white


def render_v(size: int) -> Image.Image:
    """Render the Valenx "V" mark at <size>x<size>, returning RGBA."""
    # 4x supersampling for crisp antialiasing — render large, downsample
    # to <size> with LANCZOS. The rounded-rect mask path is otherwise
    # noticeably staircased at small sizes.
    ss = 4
    big = size * ss
    img = Image.new("RGBA", (big, big), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    radius = max(2, big // 6)
    # Rounded square background.
    draw.rounded_rectangle((0, 0, big - 1, big - 1), radius=radius, fill=BRAND_BG)

    # "V" — two strokes meeting at the apex (bottom-center).
    pad = max(1, int(big * 0.18))
    stroke = max(2, int(big * 0.18))
    apex_x = big // 2
    apex_y = big - pad
    left_x = pad
    right_x = big - pad
    top_y = pad

    # joint="curve" gives a clean meeting point at the apex
    draw.line(
        [(left_x, top_y), (apex_x, apex_y), (right_x, top_y)],
        fill=BRAND_FG, width=stroke, joint="curve",
    )

    # Cap the two outer ends so they're rounded, matching the .ico's
    # Pen.StartCap/EndCap=Round. PIL doesn't expose per-segment caps,
    # so we draw filled circles at the top of each leg.
    half = stroke // 2
    draw.ellipse((left_x - half, top_y - half, left_x + half, top_y + half), fill=BRAND_FG)
    draw.ellipse((right_x - half, top_y - half, right_x + half, top_y + half), fill=BRAND_FG)
    # And one at the apex so the join can't show seams either.
    draw.ellipse((apex_x - half, apex_y - half, apex_x + half, apex_y + half), fill=BRAND_FG)

    return img.resize((size, size), Image.LANCZOS)


# ---------------------------------------------------------------------------
# .icns writer — see
# https://en.wikipedia.org/wiki/Apple_Icon_Image_format and
# https://github.com/dolanor/icns/blob/master/src/icns.rs for the
# canonical type-code table. We pack PNG payloads (modern .icns
# format; supported by every macOS version cargo-bundle targets, which
# pins osx_minimum_system_version = "11.0").
# ---------------------------------------------------------------------------

# (Apple type code, pixel size). The reader picks the closest match for
# whatever surface (Finder list, Dock at 1x/2x, Launchpad, Cmd-Tab) needs
# the icon. We include both 1x and 2x retina pairs covering 16-1024 so
# any zoom level renders crisply.
# Apple readers accept duplicate sizes across type codes — ic11 is
# technically 16x16@2x and icp5 is 32x32, but PNG bytes are identical so
# both encode from the 32x32 source. Same for ic12 (32@2x) ↔ icp6 (64),
# ic13 (128@2x) ↔ ic08 (256), ic14 (256@2x) ↔ ic09 (512).
ICNS_VARIANTS = [
    (b"icp4",   16),   # 16x16
    (b"icp5",   32),   # 32x32 (also 16@2x)
    (b"icp6",   64),   # 64x64 (also 32@2x)
    (b"ic07",  128),   # 128x128
    (b"ic08",  256),   # 256x256 (also 128@2x)
    (b"ic09",  512),   # 512x512
    (b"ic10", 1024),   # 1024x1024 (also 512@2x — Retina)
    (b"ic11",   32),   # 16x16@2x
    (b"ic12",   64),   # 32x32@2x
    (b"ic13",  256),   # 128x128@2x
    (b"ic14",  512),   # 256x256@2x
]


def png_bytes(img: Image.Image) -> bytes:
    buf = BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def build_icns(out_path: Path) -> int:
    # Render each unique size once, then re-use across the type codes
    # that share a pixel dimension.
    sizes = sorted({s for _, s in ICNS_VARIANTS})
    rendered: dict[int, bytes] = {s: png_bytes(render_v(s)) for s in sizes}

    chunks = []
    for code, size in ICNS_VARIANTS:
        payload = rendered[size]
        # Per-chunk header: 4-byte type + 4-byte BE length (incl. header).
        chunk_header = code + struct.pack(">I", len(payload) + 8)
        chunks.append(chunk_header + payload)

    body = b"".join(chunks)
    # File header: magic "icns" + total file length (BE u32).
    header = b"icns" + struct.pack(">I", len(body) + 8)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(header + body)
    return out_path.stat().st_size


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------

def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    png_path = repo_root / "packaging" / "linux" / "valenx-256.png"
    icns_path = repo_root / "packaging" / "macos" / "valenx.icns"

    png = render_v(256)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    png.save(png_path, format="PNG", optimize=True)
    print(f"Wrote {png_path} ({png_path.stat().st_size} bytes, 256x256 RGBA)")

    icns_size = build_icns(icns_path)
    print(f"Wrote {icns_path} ({icns_size} bytes, {len(ICNS_VARIANTS)} chunks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

# noozra.com download snippet

A self-contained, paste-ready block for **noozra.com** to point visitors at the
latest Valenx release. Two flavours below — pick whichever the site uses.

The links target the **GitHub Releases "latest"** redirect
(`/releases/latest`), so they keep working as new tags ship — no need to edit
the version in the URL each release.

---

## Option A — HTML (drop into a page)

```html
<!-- Valenx download buttons — paste into noozra.com -->
<div class="valenx-download">
  <h2>Download Valenx</h2>
  <p>
    Native desktop simulation suite. <strong>Pre-alpha</strong> (0.1.0-alpha.1) —
    early developer build. Requires working GPU drivers
    (Vulkan/DX12 on Windows, Vulkan/GL on Linux).
  </p>

  <p class="valenx-buttons">
    <a class="valenx-btn"
       href="https://github.com/nochallenge/valenx/releases/latest"
       rel="noopener">
      ⬇ Download for Windows (.zip)
    </a>
    <a class="valenx-btn"
       href="https://github.com/nochallenge/valenx/releases/latest"
       rel="noopener">
      ⬇ Download for Linux (.tar.gz)
    </a>
    <a class="valenx-btn valenx-btn--muted"
       href="https://github.com/nochallenge/valenx/releases/latest"
       rel="noopener">
      Download for macOS (experimental)
    </a>
  </p>

  <p class="valenx-note">
    Windows &amp; Linux are tested-build. <strong>macOS is experimental</strong> —
    it compiles in CI but the GUI is not yet verified on macOS.
    All builds are <a href="https://github.com/nochallenge/valenx/releases">on GitHub Releases</a>.
  </p>
</div>

<style>
  .valenx-download { max-width: 640px; font-family: system-ui, sans-serif; }
  .valenx-buttons { display: flex; flex-wrap: wrap; gap: .6rem; }
  .valenx-btn {
    display: inline-block; padding: .7rem 1.1rem; border-radius: 8px;
    background: #2563eb; color: #fff; text-decoration: none; font-weight: 600;
  }
  .valenx-btn:hover { background: #1d4ed8; }
  .valenx-btn--muted { background: #6b7280; }
  .valenx-btn--muted:hover { background: #4b5563; }
  .valenx-note { font-size: .9rem; color: #555; margin-top: .8rem; }
</style>
```

---

## Option B — Markdown (for a markdown-rendered site)

```markdown
## Download Valenx

Native desktop simulation suite. **Pre-alpha** (0.1.0-alpha.1) — early developer
build. Requires working GPU drivers (Vulkan/DX12 on Windows, Vulkan/GL on Linux).

- **[⬇ Download for Windows (.zip)](https://github.com/nochallenge/valenx/releases/latest)**
- **[⬇ Download for Linux (.tar.gz)](https://github.com/nochallenge/valenx/releases/latest)**
- [Download for macOS (experimental)](https://github.com/nochallenge/valenx/releases/latest)

Windows & Linux are tested-build. **macOS is experimental** — it compiles in CI
but the GUI is not yet verified on macOS. All builds are on
[GitHub Releases](https://github.com/nochallenge/valenx/releases).
```

---

### Per-OS run steps (optional — include if the site has room)

- **Windows:** unzip, run `valenx.exe` (SmartScreen may warn — *More info → Run anyway*; the pre-alpha binary is unsigned).
- **Linux:** `tar -xzf valenx-*.tar.gz`, then `chmod +x valenx && ./valenx`.
- **macOS (experimental):** `tar -xzf valenx-*.tar.gz`, then `xattr -dr com.apple.quarantine valenx && chmod +x valenx && ./valenx`.

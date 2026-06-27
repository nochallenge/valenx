Valenx — AI-controlled design and simulation
=============================================

A native, open-source desktop simulation suite written in Rust. One app
spanning engineering (aerospace, CFD, FEA, EM, thermal, CAD/CAM), chemistry
& materials, computational biology, neuroengineering, gravitational physics,
and more. No browser, no subscription, no cloud — your data stays on your
machine.

Status: 0.1.0-alpha.1 — PRE-ALPHA. Early developer build; expect sharp
edges, please file issues.


How to run
----------

Windows:
  1. Unzip this archive.
  2. Run `valenx.exe`.
  (Windows SmartScreen may warn because the binary is unsigned in this
   pre-alpha build. Choose "More info" -> "Run anyway".)

Linux:
  1. Extract:  tar -xzf valenx-<ver>-x86_64-unknown-linux-gnu.tar.gz
  2. Make it executable and run:
        chmod +x valenx
        ./valenx

macOS (EXPERIMENTAL — GUI not yet verified on macOS):
  1. Extract:  tar -xzf valenx-<ver>-aarch64-apple-darwin.tar.gz
  2. Clear the quarantine flag, then run:
        xattr -dr com.apple.quarantine valenx
        chmod +x valenx
        ./valenx


Requirements
------------

Working GPU drivers are required:
  - Windows: Vulkan or DX12
  - Linux:   Vulkan or OpenGL

This is the DEFAULT build: it is self-contained apart from the OS's
standard GPU drivers. The optional external-tool CAD adapters
(OpenCASCADE / EGL / --all-features) are NOT included.


License
-------

Dual-licensed under MIT OR Apache-2.0. See the bundled LICENSE-MIT and
LICENSE-APACHE files.

Project: https://github.com/nochallenge/valenx

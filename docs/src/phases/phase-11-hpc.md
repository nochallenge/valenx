# Phase 11 — HPC / cluster execution

**Status:** ⚪ Planned.

## Goal

Let the same `.valenx` project execute on a cluster without changes,
via SSH / Slurm / PBS. Desktop remains the authoring surface; HPC
becomes the scale surface.

## Capability inventory

- Remote run-target abstraction in `valenx-core` (SSH, Slurm, PBS,
  LSF, cloud queues).
- Case transfer: diff-based sync to the cluster's project storage.
- Live residual + log streaming over SSH.
- Remote viewport: headless render on the cluster, pixel stream back
  (for ultra-large cases that can't fit in local memory).
- Result cache with pinning rules so interesting runs don't age out.
- Credential management with OS keychain.

## Acceptance checklist

- [ ] Submit a Slurm job from the desktop app without touching a
      terminal.
- [ ] Cancel a running remote job from the app.
- [ ] Receive residuals live while the run is on a remote node.
- [ ] Download only the results you asked for; leave heavy
      intermediates on the cluster.

## Leads into

[Phase 12 — Optimisation + adjoint workflows](./phase-12-optimisation.md).

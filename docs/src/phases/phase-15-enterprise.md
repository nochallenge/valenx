# Phase 15 — Enterprise deployment

**Status:** ⚪ Planned.

## Goal

Drop Valenx into organisations without friction — MSI / PKG
installers, SSO, centralised configuration, audit logs.

## Capability inventory

- MSI installer (Windows) + PKG (macOS) + signed DEB / RPM + Snap
  / Flatpak for Linux.
- MDM-friendly: per-machine preferences via config channel
  (ConfigProfile on macOS, Group Policy on Windows, dconf/NixOS on
  Linux).
- SSO: SAML 2.0 + OIDC for access to shared result archives.
- Audit log of every run, edit, and export, with provenance chain.
- Offline install bundle with every OSS tool pre-packaged.
- Governance: per-user data-retention policy, IP / license posture
  dashboard.

## Acceptance checklist

- [ ] Silent install works on a fleet of 1000 Windows workstations.
- [ ] SSO login round-trips with an org IdP.
- [ ] Audit log survives a client crash.
- [ ] Offline mode: no external network required after install.

## Leads into

[Phase 16 — Stewardship + long-term governance](./phase-16-stewardship.md).

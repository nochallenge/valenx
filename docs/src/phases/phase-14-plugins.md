# Phase 14 — Plugin marketplace

**Status:** ⚪ Planned.

## Goal

Let third parties ship adapters, pattern solvers, and UX extensions
without forking Valenx.

## Capability inventory

- WIT-based plugin API (per [RFC 0003](../../../rfcs/0003-plugin-api.md)).
- Marketplace UX: install, enable, update, sandbox.
- Signature verification: signed manifests, first-run consent prompt.
- Permission model: filesystem scope, subprocess allowlist, network
  policy — enforced by the host, not by the plugin.
- In-app review / rating for plugins.
- Publisher portal: how a plugin author gets onto the marketplace.

## Acceptance checklist

- [ ] A third-party Python kinetics package ships as a Valenx plugin.
- [ ] Plugin cannot read / write outside its declared scope.
- [ ] Revoke a plugin's permissions without restarting the app.
- [ ] Marketplace search returns ≤ 200 ms on broadband.

## Leads into

[Phase 15 — Enterprise deployment](./phase-15-enterprise.md).

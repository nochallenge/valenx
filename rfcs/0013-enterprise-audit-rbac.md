# RFC 0013 — Enterprise audit log + role-based access control

| | |
|---|---|
| **Status** | Draft |
| **Phase target** | 15 (enterprise) |
| **Type** | Architecture |
| **Author** | Maintainers |
| **Discussion** | open in PR — not yet merged |
| **Related** | RFC 0009 (HPC executors), RFC 0014 (auth — not yet drafted) |

## Summary

Define the minimal feature set and integration points for two
enterprise prerequisites:

- **Audit log** — append-only record of every project mutation,
  every run kicked off, every settings change, every adapter
  registration. Enough detail that an auditor can reconstruct what
  happened, when, and who did it.
- **Role-based access control (RBAC)** — three roles (`viewer` /
  `runner` / `admin`) gating which workflow actions a user can
  perform on a given project.

Out of scope for this RFC:

- The auth backend (LDAP / OIDC / SAML / SSO). Covered by a
  follow-up.
- Multi-tenancy proper (separate user accounts on the same Valenx
  install). RBAC pretends-multi-tenancy is good enough for the
  primary use case (a single-machine Valenx install used by a
  small team via an existing auth-aware filesystem).
- Compliance certifications (SOC 2, HIPAA, FedRAMP).

## Motivation

Replacing commercial tools at organisations means meeting their
audit and access requirements. The two most common enterprise
asks:

1. **"Who deleted that case from the project?"** — needs an audit
   log that lists project mutations with user + timestamp.
2. **"This user can view results but can't kick off new runs"** —
   needs role-based gates on the workflow actions.

Without these, Valenx isn't deployable past a single user / single
machine into any organisation with an IT department.

## Design

### Audit log

#### Format

JSONL append-only file at `<state_dir>/audit.log.jsonl`. Each line:

```json
{
  "timestamp": "2026-04-25T14:32:11.123Z",
  "actor": {"id": "alice@example.com", "session_id": "..."},
  "action": "run_case",
  "target": {
    "kind": "case",
    "project": "/projects/airfoil-study",
    "case": "cfd-steady"
  },
  "context": {
    "adapter": "openfoam",
    "workdir": "/tmp/valenx-run-cfd-steady-1234567890"
  }
}
```

Mandatory fields: `timestamp` (ISO 8601 UTC), `actor`, `action`,
`target`. Free-form `context` for action-specific detail.

#### Action vocabulary

A closed enum, extended carefully (each new action is a tiny RFC):

| action | trigger | target |
|---|---|---|
| `project.open` | File → Open project | project path |
| `project.save` | File → Save project | project path |
| `case.create` | new case added | case |
| `case.delete` | case removed | case |
| `case.modify` | case.toml edited via UI | case + before/after diff |
| `run.start` | Run / Run from prepared / Run-selected | case + adapter |
| `run.cancel` | Cancel button | case |
| `run.complete` | run finished (success or failure) | case + exit code |
| `prepare.start` | Prepare-only action | case + adapter |
| `settings.modify` | Settings dialog change | which setting + before/after |
| `adapter.register` | adapter loaded into registry | adapter id + status |
| `plugin.load` | RFC 0003 plugin loaded | plugin id |

#### Append-only guarantees

The audit module opens the log file with `O_APPEND` (Unix) /
`FILE_APPEND_DATA` (Windows) so writes can't tear or overlap with
concurrent processes. The log is never rotated by the app — it's
the operator's responsibility (with logrotate / scheduled task) so
the audit trail's continuity isn't on Valenx's hands.

#### Integrity

Each entry's `prev_hash` field carries the SHA-256 of the previous
entry's serialisation. Tampering is detectable by replaying the
chain. Optional cryptographic signing (per-install key) is a
follow-up.

### RBAC

#### Roles

| role | read project | edit case.toml | run / prepare | manage adapters / settings |
|---|---|---|---|---|
| `viewer` | yes | no | no | no |
| `runner` | yes | yes | yes | no |
| `admin` | yes | yes | yes | yes |

#### Where roles live

- **Per-user assignment** — `<state_dir>/rbac.json` maps user id →
  role. Default role `runner` for unassigned users (sensible single-
  user-machine default; admins can override).
- **Per-project override** — optional `[rbac]` block in the project's
  `project.toml` specifying per-user roles for that project, taking
  precedence over the global file.

#### Enforcement points

The Adapter trait does not gain an auth parameter — that's
infrastructure leaking into physics. Instead, the app's workflow
methods (`run_selected_case`, `prepare_selected_case`, etc.)
gate-check at the entry. A new `AppPermissions` struct on
`ValenxApp` wraps the role lookup; every workflow method calls
`self.permissions.require("run.case")?` before doing anything.

A failed permission check:

- Returns a structured error to `last_error` ("You don't have
  permission to run cases. Ask the admin for the `runner` role.").
- Emits an audit-log entry with the attempted action + denied flag.
- Doesn't perform any side effect (no workdir created, no log
  line written to the run pipe).

### Identity

For this RFC: `actor.id` is just whatever string the host OS
reports as the current user (`whoami` on Unix, `GetUserNameW` on
Windows). No SSO, no token validation, no impersonation prevention.

Real auth (RFC 0014, not yet drafted) replaces this with proper
sessions + OIDC / SAML / SSO + token rotation.

## Drawbacks

- Audit + RBAC add real surface area to every workflow method.
  Even when nobody uses them (single-user installs), the
  permission-check call is in every code path.
- Append-only file integrity is best-effort on consumer hosts.
  Filesystem corruption / partial writes / power loss between
  the write and the fsync can break the chain. Write-ahead-log
  semantics are out of scope.
- Per-project `[rbac]` overrides can conflict with the global
  `rbac.json`. Conflict resolution (project wins) needs clear
  documentation; users will get confused otherwise.

## Alternatives considered

- **Punt audit + RBAC entirely to the OS** — file permissions
  control project access; a wrapper script logs invocations.
  Works for the simplest setups; falls down for "viewer can read
  but not run" because at the OS layer a read can already be
  abused into a run (just `cat case.toml | python ...`).
- **Embed in Adapter trait** — every adapter checks permissions.
  Rejected: 11 adapters × 3 roles × N actions = lots of duplicate
  check logic, all bound to drift.
- **Separate auth daemon** — a sidecar process (think PAM /
  polkit). Too heavy for v0; might revisit for the SSO RFC.

## Migration path

Phase 15.0 (this RFC):
- New crate `valenx-audit` with the JSONL writer + integrity-chain
  helpers.
- New crate `valenx-rbac` with role enum + `AppPermissions` +
  `rbac.json` parser.
- App workflow methods grow `permissions.require(...)?` checks;
  audit emit calls land at the same points.

Phase 15.1: per-project `[rbac]` block in `project.toml`.

Phase 15.2: `valenx audit` CLI for reading + verifying the log
chain.

Phase 15.3: hooks for external SIEM ingestion (syslog / journald
forwarder).

## Open questions

1. **What is "the actor" when Valenx runs a scheduled batch run?**
   Multi-user audit gets weird when there's no human at the keyboard.
   Probably `actor.id = "system:scheduler"` with a session id
   referencing the cron / systemd unit.
2. **How long do we keep the log?** Operator decision, but the app
   should warn when the log file passes some size threshold so it
   doesn't grow unboundedly on a long-running install.
3. **GDPR / right-to-erasure.** EU privacy law might force redaction
   of an actor.id from old entries. Tension with the integrity
   chain (redaction breaks the hash). Solution probably: cryptographic
   redaction (replace name with a fixed marker, leave hash valid).
4. **Plugin sandbox.** RFC 0003 says nothing about whether a plugin
   can write its own audit entries. We need a "plugins emit through
   a host-controlled API only" rule so a malicious plugin can't
   forge an entry attributed to a different user.

## References

- RFC 0003 — plugin API (interacts with audit / RBAC for plugin-
  initiated actions).
- RFC 0009 — HPC executors (remote runs need their actor.id
  forwarded across the wire).
- NIST SP 800-92 — guide to computer security log management.
- OWASP Audit Logging Cheat Sheet.
- AWS CloudTrail event schema (reference for the action vocabulary).

# Releasing

How to cut a tagged release of Valenx. Operationalises lane D of
[NEXT_PHASE.md](./NEXT_PHASE.md).

> **Maintainers only.** This doc assumes you have push access to
> `master` and the maintainer role required to draft GitHub
> Releases. If you're a contributor wanting to ship something,
> open an RFC or coordinate with a maintainer first.

## Overview

A release ships:

| Artefact | Built by | Signed? |
|---|---|---|
| `valenx_<version>_amd64.deb` | `cargo deb` (Linux runner) | No (distros sign repo indices) |
| `valenx-<version>-1.x86_64.rpm` | `cargo generate-rpm` (Linux runner) | No |
| `Valenx-<version>.dmg` (notarised) | `cargo bundle` + `codesign` + `notarytool` + `hdiutil` (macOS runner) | Yes — Apple Developer ID |
| `valenx-<version>-x64.msi` | `cargo wix` + `signtool` (Windows runner) | Yes — Authenticode |

All four are produced by `.github/workflows/release.yml`, which is
currently **`workflow_dispatch`-only** (manually triggered from the
GitHub Actions tab). The pre-alpha workflow was originally tag-push
triggered but was switched to manual-only as part of the round-1/2
CI rework — see `docs/CI.md` for the rationale. Releases today:
maintainer dispatches `release.yml` with the desired tag input, the
workflow builds + uploads the four artefacts to a draft GitHub
Release, the maintainer smoke-tests one artefact per platform, and
publishes.

## Pre-alpha shortcut: skip the certs

For `v0.1.0-alpha.1` the project ships **unsigned binaries** by
default. The release workflow detects when the cert secrets are
unset and falls through to an unsigned build; users see one
Gatekeeper / SmartScreen warning the first time they run the app
and bypass it via right-click → Open (macOS) or "More info → Run
anyway" (Windows).

This is the conventional shape for OSS pre-alphas — rustup, uv,
ripgrep, and many others started this way. Skipping signing
saves ~$99/yr (Apple) + ~$100-300/yr (Authenticode) until the
user base justifies the cost.

To ship unsigned: just **don't add** the secrets in the next
section. Manual `workflow_dispatch` → workflow runs → unsigned
`.dmg` / `.msi` upload to the GitHub Release with a CI warning
noting the unsigned status. No further action needed.

The rest of this document covers the **signed path** that lands
in v0.2.0 (or earlier if a sponsor covers the certs).

## One-time setup (v0.2.0+ — signed path)

### Apple Developer ID (macOS)

1. Enrol an organisation account at
   <https://developer.apple.com/programs/enroll/>. Lead time
   ~1 week (DUNS verification, payment).
2. Generate a Developer ID Application certificate at Xcode →
   Settings → Accounts → Manage Certificates → "+" → "Developer
   ID Application".
3. Export to a `.p12` (right-click in Keychain Access).
4. Base64-encode and add as repo secret:
   ```bash
   base64 -i developer_id.p12 | pbcopy
   ```
   - `APPLE_DEVELOPER_ID_CERT` — the base64 blob
   - `APPLE_DEVELOPER_ID_CERT_PASSWORD` — the .p12 password
   - `APPLE_ID` — the Apple ID email
   - `APPLE_TEAM_ID` — visible in Apple Developer → Membership
   - `APPLE_APP_PASSWORD` — generated at appleid.apple.com →
     Sign-In and Security → App-Specific Passwords. Used by
     `notarytool`.

### Authenticode (Windows)

1. Buy an OV or EV code-signing certificate. OV is sufficient for
   open-source distribution; EV gets you immediate SmartScreen
   reputation. Lead time ~1 week (identity verification).
   Recommended issuers: DigiCert, Sectigo, GlobalSign.
2. Receive the cert as `.pfx` (or convert from `.p7b`).
3. Base64-encode and add as repo secret:
   ```bash
   base64 -i cert.pfx | clip
   ```
   - `WINDOWS_CERT` — base64 blob
   - `WINDOWS_CERT_PASSWORD` — .pfx password

### Linux

No certs needed at the moment. If we later want signed packages
for Debian / Fedora repos, add a GPG key as
`LINUX_GPG_PRIVATE_KEY` + `LINUX_GPG_KEY_PASSWORD` and extend the
workflow to call `dpkg-sig` / `rpm --addsign`.

## Cutting a release

### 1. Update the changelog

Move the items from `## [Unreleased]` to a new
`## [<version>] — <date>` block in
[CHANGELOG.md](./CHANGELOG.md). Keep the structure — Added /
Fixed / Changed / Quality.

Commit on a release-prep branch:

```bash
git checkout -b release/v0.1.0-alpha.1
$EDITOR CHANGELOG.md
git commit -am "chore(release): prep v0.1.0-alpha.1 changelog"
git push -u origin release/v0.1.0-alpha.1
gh pr create --title "chore(release): prep v0.1.0-alpha.1" \
  --body "Move Unreleased → 0.1.0-alpha.1 in CHANGELOG."
```

Merge after one approval.

### 2. Update Cargo.toml workspace version

Bump `version` in the root `Cargo.toml` to match the tag you
plan to dispatch. Commit on master (release-prep PRs include
this).

### 3. Dispatch the release workflow

The release pipeline is a manual `workflow_dispatch` — there's
no auto-tag-then-build trigger. Pass the tag as an input;
unsigned builds run with `draft=true` and land in a draft
GitHub Release for manual smoke-testing.

```bash
gh workflow run release.yml \
  -f tag=v0.1.0-alpha.1 \
  -f draft=true
```

### 4. Watch the release workflow

```bash
# Pick up the run id of the dispatch you just triggered.
RUN_ID=$(gh run list --workflow=release.yml --limit 1 \
  --json databaseId --jq '.[0].databaseId')
gh run watch "$RUN_ID"
```

The four build jobs (`linux-deb` / `linux-rpm` / `macos` /
`windows`) run in parallel and feed `publish` (which assembles
the GitHub Release). Total wall-clock time: ~25–35 minutes
depending on cache hits.

### 5. Smoke test the draft Release

The workflow leaves the Release in **draft** state. Download one
artefact per platform and verify:

| Platform | Smoke test |
|---|---|
| Linux .deb | `sudo dpkg -i valenx_*.deb && valenx --version && valenx-validate --version` |
| Linux .rpm | `sudo rpm -i valenx-*.rpm && valenx --version` |
| macOS .dmg | `hdiutil attach Valenx-*.dmg && open /Volumes/Valenx/Valenx.app` (no Gatekeeper warning if signed) |
| Windows .msi | Double-click → install → Start menu → Valenx (no SmartScreen warning if signed) |

If any smoke test fails, **don't publish** — investigate, fix,
delete the tag, repeat.

### 6. Publish

In the GitHub UI: open the draft Release → review the auto-
generated changelog → "Publish release". This triggers the
"Released" event consumers might be subscribed to (Homebrew tap
auto-update, distro packagers, etc.).

### 7. Post-release

- Open a PR bumping `version` in `Cargo.toml` to the next
  development version (e.g. `0.1.0-alpha.2-dev`).
- Announce in the project's discussion / matrix / wherever the
  community lives.
- File any "found during smoke test" follow-ups as issues against
  the next version.

## Hotfix releases

If a critical bug surfaces after publish:

1. Branch from the tag: `git checkout -b hotfix/v0.1.0-alpha.1.1 v0.1.0-alpha.1`.
2. Cherry-pick the fix.
3. Bump the patch identifier (`0.1.0-alpha.1.1`).
4. Trigger the release workflow with the new tag string —
   `release.yml` is `workflow_dispatch`-only (round-4 / round-8
   change), no tag push fires it. From the CLI:
   ```bash
   gh workflow run release.yml -f tag=v0.1.0-alpha.1.1
   ```
   Or from the GitHub UI: Actions → "Release" → "Run workflow",
   choose the hotfix branch, fill in the tag input.
5. Re-merge into master after the release ships so master doesn't
   regress.

## Failure modes

### Apple notarisation fails

`xcrun notarytool submit` returns a UUID. If status comes back
`Invalid`, fetch the log:

```bash
xcrun notarytool log <UUID> \
  --apple-id "$APPLE_ID" --team-id "$APPLE_TEAM_ID" \
  --password "$APPLE_APP_PASSWORD"
```

The two most common rejections:
- **"signed with a non-Developer ID certificate"** — make sure
  the keychain only has the Developer ID cert imported. Ad-hoc
  signing won't notarise.
- **"executable lacks a hardened runtime"** — `cargo bundle` may
  not opt in by default; the workflow passes `--options runtime`
  to `codesign` to fix this.

### Authenticode timestamp fails

DigiCert's timestamp server is the default but occasionally
returns 503. The workflow uses a single timestamp URL — if it's
flaky, retry the run. For a permanent fix, fall back to
`http://timestamp.sectigo.com` with a second `signtool` attempt.

### Tag was wrong

`release.yml` is `workflow_dispatch`-only (round-4 / round-8
change), so a wrong tag does NOT auto-publish — the workflow only
runs when you trigger it with `gh workflow run release.yml -f
tag=...`. Just re-trigger with the corrected tag string:

```bash
gh workflow run release.yml -f tag=v0.1.0-alpha.1
```

If you DID push a literal git tag that you want gone, delete it
locally and on the remote:

```bash
git tag -d v0.1.0-alpha.1
git push --delete origin v0.1.0-alpha.1
```

The workflow's `concurrency` group handles in-flight cancellation
when you re-trigger.

## What's not in scope yet

- **Nightly channel.** Tags only for now. Nightly builds land
  with v0.2.0+ — they need a separate workflow that uploads to
  GitHub Pages or an S3 bucket rather than a Release.
- **Auto-update.** No in-app updater. Users download fresh
  installers per release.
- **Reproducible builds.** `Cargo.lock` is committed; runner
  binaries differ between platforms (timestamps, build paths).
  Phase 14 / 15 territory.
- **Distro repositories.** No Debian / Fedora / Homebrew /
  winget submission workflow. Once 0.1.0 ships and stabilises,
  the standalone artefacts are enough — repository submissions
  follow.

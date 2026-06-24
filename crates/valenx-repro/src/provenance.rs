//! Content-addressed **provenance / lineage** for a single computation.
//!
//! Where [`crate::ReproBundle`] packages a whole study, a [`Manifest`] is the
//! leaner, self-verifying lineage record of **one** computation: which tool (at
//! which version) consumed which content-addressed inputs, under which parameters,
//! to produce which content-addressed outputs. It carries a single stable
//! [`digest`](Manifest::digest) over all of that, and a [`verify`](Manifest::verify)
//! that **re-derives** the digest from the recorded fields and **fails loud** on
//! any mismatch — so a tampered manifest is caught, not trusted.
//!
//! The digest is canonical (inputs and outputs are sorted by name before
//! hashing), so it is independent of the order entries were added: the same
//! lineage always yields the same digest.
//!
//! ```
//! use valenx_repro::provenance::{HashedItem, Manifest};
//!
//! // Hash two inputs and one output by their bytes, then record the lineage.
//! let manifest = Manifest::builder("translate", "0.1.0")
//!     .param("genetic_code", "1")
//!     .input(HashedItem::from_bytes("mstn.fasta", b"ATGCAA"))
//!     .output(HashedItem::from_bytes("protein.txt", b"MQ"))
//!     .build()
//!     .unwrap();
//!
//! // The digest is a SHA-256 in lower-case hex.
//! assert_eq!(manifest.digest().len(), 64);
//! // A freshly built manifest verifies.
//! manifest.verify().unwrap();
//! ```
//!
//! ## Honest scope
//!
//! Content addressing proves a manifest is **internally consistent and
//! unchanged** — that the recorded outputs were produced from the recorded
//! inputs+params by the recorded tool, and that nobody has edited the record
//! since. It does **not** prove the computation was *correct*, nor does it
//! re-execute anything. Re-running the tool to confirm the outputs is the
//! caller's job; this records and checks the lineage.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::error::ReproError;

/// Lower-case hex SHA-256 of `bytes`.
fn sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    to_hex(&hasher.finalize())
}

/// Render bytes as lower-case hex.
fn to_hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        let _ = write!(s, "{b:02x}");
    }
    s
}

/// A named, content-addressed item: a name paired with the SHA-256 of some
/// bytes. Used for both inputs and outputs of a [`Manifest`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct HashedItem {
    /// File name or logical identifier.
    pub name: String,
    /// Lower-case hex SHA-256 of the item's bytes.
    pub hash: String,
}

impl HashedItem {
    /// Build a hashed item by hashing `bytes` with SHA-256.
    pub fn from_bytes(name: impl Into<String>, bytes: &[u8]) -> Self {
        Self {
            name: name.into(),
            hash: sha256_hex(bytes),
        }
    }

    /// Build a hashed item from a name and an already-computed lower-case hex
    /// SHA-256. The hash is **not** re-validated here; use this when the bytes
    /// are large or already digested elsewhere.
    pub fn new(name: impl Into<String>, hash: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            hash: hash.into(),
        }
    }
}

/// A named parameter `(name, value)` recorded for a computation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParamEntry {
    /// Parameter name.
    pub name: String,
    /// Parameter value, rendered as text.
    pub value: String,
}

/// A self-verifying, content-addressed provenance record for **one** computation.
///
/// Build it with [`Manifest::builder`]; read off [`Manifest::digest`] (a stable
/// SHA-256 over tool + version + params + sorted inputs + sorted outputs) and
/// call [`Manifest::verify`] to re-derive the digest and confirm nothing has been
/// tampered with.
///
/// The stored [`digest`](Self::digest) is computed at build time and serialised
/// with the record. `verify()` recomputes it from the other fields and rejects
/// any mismatch — so an attacker who edits an input hash (or a parameter, or the
/// tool version) without also forging the digest is caught.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Manifest {
    /// The tool / function that ran (e.g. `"translate"`).
    pub tool: String,
    /// The tool's version string (e.g. `"0.1.0"`).
    pub version: String,
    /// Content-addressed inputs.
    pub inputs: Vec<HashedItem>,
    /// Content-addressed outputs.
    pub outputs: Vec<HashedItem>,
    /// The SHA-256 over the canonicalised parameter set (lower-case hex). This
    /// is recorded explicitly so a manifest can prove which parameters were used
    /// even if the full parameter list is summarised or dropped downstream.
    pub params_hash: String,
    /// The recorded parameters (kept so the manifest is self-describing and
    /// `params_hash` can be re-derived during [`Manifest::verify`]).
    pub params: Vec<ParamEntry>,
    /// The stable SHA-256 digest over the whole record (lower-case hex). Set at
    /// build time; recomputed and checked by [`Manifest::verify`].
    pub digest: String,
}

impl Manifest {
    /// Start building a manifest for `tool` at `version`.
    pub fn builder(tool: impl Into<String>, version: impl Into<String>) -> ManifestBuilder {
        ManifestBuilder {
            tool: tool.into(),
            version: version.into(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            params: Vec::new(),
        }
    }

    /// The stable content-addressed digest (lower-case hex SHA-256).
    pub fn digest(&self) -> &str {
        &self.digest
    }

    /// The SHA-256 over the parameter set (lower-case hex).
    pub fn params_hash(&self) -> &str {
        &self.params_hash
    }

    /// Re-derive the digest from the recorded fields and confirm it matches the
    /// stored one — the **fail-loud** integrity check.
    ///
    /// This recomputes both [`params_hash`](Self::params_hash) (from `params`)
    /// and the full [`digest`](Self::digest) (from tool + version + the
    /// recomputed params hash + sorted inputs + sorted outputs), and returns an
    /// error on the *first* mismatch.
    ///
    /// # Errors
    ///
    /// - [`ReproError::ParamsHashMismatch`] if the stored `params_hash` does not
    ///   match the hash recomputed from `params`.
    /// - [`ReproError::DigestMismatch`] if the stored `digest` does not match the
    ///   digest recomputed from the (now-validated) fields.
    pub fn verify(&self) -> Result<(), ReproError> {
        let expected_params = hash_params(&self.params);
        if expected_params != self.params_hash {
            return Err(ReproError::ParamsHashMismatch {
                expected: expected_params,
                found: self.params_hash.clone(),
            });
        }
        let expected_digest = compute_digest(
            &self.tool,
            &self.version,
            &expected_params,
            &self.inputs,
            &self.outputs,
        );
        if expected_digest != self.digest {
            return Err(ReproError::DigestMismatch {
                expected: expected_digest,
                found: self.digest.clone(),
            });
        }
        Ok(())
    }

    /// A short human-readable summary: tool, counts, and the digest.
    pub fn summary(&self) -> String {
        format!(
            "Provenance manifest: {tool} {version}\n  inputs: {ni}  outputs: {no}  params: {np}\n  params_hash (sha256): {ph}\n  digest (sha256): {d}",
            tool = self.tool,
            version = self.version,
            ni = self.inputs.len(),
            no = self.outputs.len(),
            np = self.params.len(),
            ph = self.params_hash,
            d = self.digest,
        )
    }
}

/// Builder for a [`Manifest`]. Add inputs, outputs and params, then
/// [`build`](ManifestBuilder::build) to compute the hashes and seal the record.
#[derive(Debug, Clone)]
pub struct ManifestBuilder {
    tool: String,
    version: String,
    inputs: Vec<HashedItem>,
    outputs: Vec<HashedItem>,
    params: Vec<ParamEntry>,
}

impl ManifestBuilder {
    /// Add a content-addressed input.
    pub fn input(mut self, item: HashedItem) -> Self {
        self.inputs.push(item);
        self
    }

    /// Add a content-addressed output.
    pub fn output(mut self, item: HashedItem) -> Self {
        self.outputs.push(item);
        self
    }

    /// Record a parameter `(name, value)`.
    pub fn param(mut self, name: impl Into<String>, value: impl Into<String>) -> Self {
        self.params.push(ParamEntry {
            name: name.into(),
            value: value.into(),
        });
        self
    }

    /// Seal the manifest: validate, compute `params_hash` and the full `digest`,
    /// and return the immutable [`Manifest`].
    ///
    /// # Errors
    ///
    /// - [`ReproError::Empty`] if `tool` or `version` is empty / whitespace.
    /// - [`ReproError::DuplicateName`] if two inputs, or two outputs, share a
    ///   name (a name must address one item within its role, so the canonical,
    ///   sorted hash is unambiguous).
    pub fn build(self) -> Result<Manifest, ReproError> {
        if self.tool.trim().is_empty() {
            return Err(ReproError::Empty { field: "tool" });
        }
        if self.version.trim().is_empty() {
            return Err(ReproError::Empty { field: "version" });
        }
        check_unique_names(&self.inputs, "input")?;
        check_unique_names(&self.outputs, "output")?;

        let params_hash = hash_params(&self.params);
        let digest = compute_digest(
            &self.tool,
            &self.version,
            &params_hash,
            &self.inputs,
            &self.outputs,
        );
        Ok(Manifest {
            tool: self.tool,
            version: self.version,
            inputs: self.inputs,
            outputs: self.outputs,
            params_hash,
            params: self.params,
            digest,
        })
    }
}

/// Reject duplicate names within a single role (inputs or outputs).
fn check_unique_names(items: &[HashedItem], role: &'static str) -> Result<(), ReproError> {
    for (i, a) in items.iter().enumerate() {
        for b in &items[i + 1..] {
            if a.name == b.name {
                return Err(ReproError::DuplicateName {
                    role,
                    name: a.name.clone(),
                });
            }
        }
    }
    Ok(())
}

/// SHA-256 over the canonicalised parameter set (params sorted by name, then by
/// value to break a duplicate-name tie deterministically). NUL-delimited so the
/// boundaries are unambiguous.
fn hash_params(params: &[ParamEntry]) -> String {
    let mut sorted: Vec<&ParamEntry> = params.iter().collect();
    sorted.sort_by(|a, b| a.name.cmp(&b.name).then(a.value.cmp(&b.value)));
    let mut buf = String::new();
    for p in sorted {
        buf.push_str("param\0");
        buf.push_str(&p.name);
        buf.push('\0');
        buf.push_str(&p.value);
        buf.push('\0');
    }
    sha256_hex(buf.as_bytes())
}

/// The canonical digest over the whole record. Inputs and outputs are sorted by
/// name so insertion order does not matter; every field is NUL-delimited and
/// role-tagged so no field boundary is ambiguous.
fn compute_digest(
    tool: &str,
    version: &str,
    params_hash: &str,
    inputs: &[HashedItem],
    outputs: &[HashedItem],
) -> String {
    let mut buf = String::new();
    buf.push_str("tool\0");
    buf.push_str(tool);
    buf.push_str("\0version\0");
    buf.push_str(version);
    buf.push_str("\0params\0");
    buf.push_str(params_hash);

    let mut ins: Vec<&HashedItem> = inputs.iter().collect();
    ins.sort_by(|a, b| a.name.cmp(&b.name).then(a.hash.cmp(&b.hash)));
    for item in ins {
        buf.push_str("\0in\0");
        buf.push_str(&item.name);
        buf.push('\0');
        buf.push_str(&item.hash);
    }

    let mut outs: Vec<&HashedItem> = outputs.iter().collect();
    outs.sort_by(|a, b| a.name.cmp(&b.name).then(a.hash.cmp(&b.hash)));
    for item in outs {
        buf.push_str("\0out\0");
        buf.push_str(&item.name);
        buf.push('\0');
        buf.push_str(&item.hash);
    }

    sha256_hex(buf.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Pinned hashes for the fixed lineage in `digest_is_pinned_for_a_fixed_input`.
    const PINNED_PARAMS_HASH: &str =
        "02a1f1a4418020d936a81163bdd7c6c933cb978e2f3b0e4356aa53c4be5840d9";
    const PINNED_DIGEST: &str = "b95205b14a1bd1d3f443e8ee171c9b0f035cbce1614cab8989ab12a82edcae06";

    fn sample() -> Manifest {
        Manifest::builder("translate", "0.1.0")
            .param("genetic_code", "1")
            .param("frame", "+1")
            .input(HashedItem::from_bytes("mstn.fasta", b"ATGCAA"))
            .output(HashedItem::from_bytes("protein.txt", b"MQ"))
            .build()
            .unwrap()
    }

    #[test]
    fn item_hashes_known_value() {
        // SHA-256("abc") is a well-known constant.
        let it = HashedItem::from_bytes("x", b"abc");
        assert_eq!(
            it.hash,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    /// Benchmark-pin: a fixed lineage hashes to a known, stable digest. If this
    /// value ever changes, the canonical hashing scheme changed — that is a
    /// breaking change to provenance and must be deliberate.
    #[test]
    fn digest_is_pinned_for_a_fixed_input() {
        let m = Manifest::builder("translate", "0.1.0")
            .param("genetic_code", "1")
            .input(HashedItem::from_bytes("in", b"abc"))
            .output(HashedItem::from_bytes("out", b"MQ"))
            .build()
            .unwrap();
        // Pinned values computed once from the canonical scheme in this module.
        assert_eq!(m.params_hash, PINNED_PARAMS_HASH);
        assert_eq!(m.digest, PINNED_DIGEST);
        // The pinned manifest must of course verify.
        m.verify().unwrap();
    }

    #[test]
    fn digest_is_deterministic_and_64_hex() {
        assert_eq!(sample().digest(), sample().digest());
        assert_eq!(sample().digest().len(), 64);
        assert!(sample().digest().chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn digest_is_insertion_order_independent() {
        let a = Manifest::builder("t", "1")
            .input(HashedItem::from_bytes("a", b"1"))
            .input(HashedItem::from_bytes("b", b"2"))
            .output(HashedItem::from_bytes("y", b"9"))
            .output(HashedItem::from_bytes("x", b"8"))
            .build()
            .unwrap();
        let b = Manifest::builder("t", "1")
            .input(HashedItem::from_bytes("b", b"2"))
            .input(HashedItem::from_bytes("a", b"1"))
            .output(HashedItem::from_bytes("x", b"8"))
            .output(HashedItem::from_bytes("y", b"9"))
            .build()
            .unwrap();
        assert_eq!(a.digest(), b.digest());
        assert_eq!(a.params_hash(), b.params_hash());
    }

    #[test]
    fn fresh_manifest_verifies() {
        sample().verify().unwrap();
    }

    #[test]
    fn tampered_input_hash_fails_verify_loud() {
        let mut m = sample();
        // Flip one input hash without re-sealing the digest.
        m.inputs[0].hash = "deadbeef".to_string();
        let err = m.verify().unwrap_err();
        assert_eq!(err.code(), "digest-mismatch");
    }

    #[test]
    fn tampered_param_fails_verify_loud() {
        let mut m = sample();
        // Edit a parameter value; params_hash no longer matches `params`.
        m.params[0].value = "99".to_string();
        let err = m.verify().unwrap_err();
        assert_eq!(err.code(), "params-hash-mismatch");
    }

    #[test]
    fn tampered_tool_fails_verify_loud() {
        let mut m = sample();
        m.tool = "exfiltrate".to_string();
        let err = m.verify().unwrap_err();
        assert_eq!(err.code(), "digest-mismatch");
    }

    #[test]
    fn forged_digest_is_rejected() {
        let mut m = sample();
        m.digest = "0".repeat(64);
        assert!(m.verify().is_err());
    }

    #[test]
    fn rejects_empty_tool_or_version() {
        assert_eq!(
            Manifest::builder("   ", "1").build().unwrap_err().code(),
            "empty"
        );
        assert_eq!(
            Manifest::builder("t", "").build().unwrap_err().code(),
            "empty"
        );
    }

    #[test]
    fn rejects_duplicate_input_name() {
        let err = Manifest::builder("t", "1")
            .input(HashedItem::from_bytes("dup", b"1"))
            .input(HashedItem::from_bytes("dup", b"2"))
            .build()
            .unwrap_err();
        assert_eq!(err.code(), "duplicate-name");
    }

    #[test]
    fn summary_mentions_counts_and_digest() {
        let s = sample().summary();
        assert!(s.contains(sample().digest()));
        assert!(s.contains("inputs: 1"));
        assert!(s.contains("outputs: 1"));
    }
}

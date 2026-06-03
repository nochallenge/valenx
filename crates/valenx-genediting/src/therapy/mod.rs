//! Gene therapy & safety — Group E (features 24–26).
//!
//! The delivery and safety layer of the crate:
//!
//! - [`cassette`] — feature 24: gene-therapy expression-cassette
//!   design (promoter + transgene + polyA + regulatory elements) with
//!   AAV (~4.7 kb) / lentivirus payload-size checks.
//! - [`delivery`] — feature 25: informational delivery-vector
//!   planning — AAV serotype tropism notes and lipid-nanoparticle
//!   payload limits.
//! - [`safety`] — feature 26: safety-screen aggregation — an
//!   off-target tally from [`valenx_genomics`], rule-based
//!   genotoxicity flags and essential-gene-proximity warnings.
//! - [`safety_db`] — commercial-depth curated reference-gene
//!   catalogues (~110 essential genes, ~110 cancer drivers, the three
//!   canonical safe-harbor loci); [`safety::safety_screen`] consumes
//!   it to produce a per-edit risk verdict.
//!
//! ## v1 scope
//!
//! The cassette designer tracks sized, named elements (it does not
//! assemble the literal vector backbone); the delivery module is an
//! informational reference, not a biodistribution model; the safety
//! aggregator applies transparent rule-based flags, not a trained
//! genotoxicity predictor. See each module's note.

pub mod cassette;
pub mod delivery;
pub mod safety;
pub mod safety_db;

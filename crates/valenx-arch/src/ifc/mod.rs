//! IFC4 export — a minimal-but-real ISO-10303-21 writer.
//!
//! IFC files use the STEP Part 21 syntax (`#1=IFCPROJECT(...)`) with
//! the IFC4 schema. The full IFC4 schema has ~1500 entity types; v1
//! covers the canonical project / site / building / storey hierarchy
//! plus the nine arch entity kinds.
//!
//! ## Output guarantees
//!
//! Every file we emit:
//! - starts with the literal `ISO-10303-21;` line, followed by a
//!   well-formed `HEADER;…ENDSEC;`,
//! - declares the `IFC4` schema in `FILE_SCHEMA(('IFC4'));`,
//! - emits one `DATA;…ENDSEC;` block with `#N=ENTITY(...)` lines,
//! - ends with `END-ISO-10303-21;`.
//!
//! Every IFC entity that requires an `IfcGloballyUniqueId` carries
//! one generated via [`writer::ifc_guid_v4`].

pub mod writer;

pub use writer::{
    emit_pset, emit_rel_space_boundary, emit_rel_voids_element, ifc_guid_v4, write_cable,
    write_chimney, write_conduit, write_covering, write_curtain_wall, write_document, write_duct,
    write_footing, write_furnishing, write_mep_equipment, write_opening_for_door,
    write_opening_for_window, write_pile, write_pipe, write_railing, write_ramp, IfcWriter,
    PropValue,
};

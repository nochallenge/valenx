//! Assembly export hook — bridge a [`KicadBoard`] into the Phase 6
//! assembly machinery.

use std::path::Path;

use valenx_assembly::Assembly;
use valenx_assembly::part::Part;
use valenx_cad::Solid;

use crate::board::KicadBoard;
use crate::error::KicadError;
use crate::tessellate::pcb_to_solid;

/// Build an [`Assembly`] from the board + an externally supplied
/// list of (ref_designator, 3D solid) tuples for the components.
///
/// The board itself becomes the first (fixed) part; each matching
/// component is added at its placed position. Components without
/// a matching solid in `components_3d` are skipped silently.
pub fn build_assembly(
    board: &KicadBoard,
    components_3d: &[(String, Solid)],
) -> Result<Assembly, KicadError> {
    let mut asm = Assembly::new();

    let board_solid = pcb_to_solid(board)?;
    let mut board_part = Part::new(0, "PCB", board_solid);
    board_part.fixed = true;
    asm.add_part(board_part);

    for comp in &board.components {
        if let Some((_, solid)) = components_3d
            .iter()
            .find(|(name, _)| name == &comp.ref_designator)
        {
            let part = Part::new(
                asm.parts.len(),
                comp.ref_designator.clone(),
                solid.clone(),
            );
            asm.add_part(part);
        }
    }

    Ok(asm)
}

/// Stub for STEP-assembly export. v1 builds the assembly + delegates
/// the STEP write to the Phase 6 assembly export machinery (which
/// itself is gated until truck-stepio gains assembly support). For
/// now we surface the deferral via a typed error so the caller can
/// fall back to per-part export.
pub fn export_step_with_components(
    board: &KicadBoard,
    components_3d: &[(String, Solid)],
    _path: impl AsRef<Path>,
) -> Result<(), KicadError> {
    let _asm = build_assembly(board, components_3d)?;
    Err(KicadError::NotImplemented(
        "STEP assembly export delegates to Phase 6 — gated until truck-stepio assembly support lands",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_assembly_includes_board_as_fixed_first_part() {
        let b = KicadBoard::demo_devboard();
        let asm = build_assembly(&b, &[]).unwrap();
        assert!(asm.parts[0].fixed);
        assert_eq!(asm.parts[0].name, "PCB");
    }

    #[test]
    fn build_assembly_skips_components_without_matching_solid() {
        let b = KicadBoard::demo_devboard();
        let asm = build_assembly(&b, &[]).unwrap();
        // Only the board — components have no provided solids.
        assert_eq!(asm.parts.len(), 1);
    }

    #[test]
    fn step_export_surfaces_deferred_error() {
        let b = KicadBoard::demo_devboard();
        let err = export_step_with_components(
            &b,
            &[],
            std::env::temp_dir().join("kicad_test.step"),
        )
        .unwrap_err();
        assert!(matches!(err, KicadError::NotImplemented(_)));
    }
}

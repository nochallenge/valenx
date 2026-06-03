//! Bridge — connect two edge loops with a strip of quads. Both
//! loops must have the same number of vertices.

use crate::error::BlenderOpError;
use crate::mesh::Mesh;

/// Add a quad strip connecting two edge loops.
///
/// `loop1` and `loop2` are vertex-id rings of the same length. The
/// stitched quads use the order of `loop1[i] -> loop1[i+1] ->
/// loop2[i+1] -> loop2[i]`.
pub fn edge_loops(
    mesh: &Mesh,
    loop1: &[usize],
    loop2: &[usize],
) -> Result<Mesh, BlenderOpError> {
    if loop1.len() != loop2.len() {
        return Err(BlenderOpError::Topology(format!(
            "loops differ in length ({} vs {})",
            loop1.len(),
            loop2.len()
        )));
    }
    if loop1.len() < 2 {
        return Err(BlenderOpError::Topology(
            "loops must have >= 2 vertices".into(),
        ));
    }
    let nv = mesh.vertices.len();
    for &v in loop1.iter().chain(loop2.iter()) {
        if v >= nv {
            return Err(BlenderOpError::IndexOutOfRange {
                kind: "vertex",
                idx: v,
                limit: nv,
            });
        }
    }
    let mut out = mesh.clone();
    let k = loop1.len();
    for i in 0..k {
        let a = loop1[i];
        let b = loop1[(i + 1) % k];
        let c = loop2[(i + 1) % k];
        let d = loop2[i];
        out.faces.push(vec![a, b, c, d]);
    }
    Ok(out)
}

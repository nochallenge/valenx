//! Procedural 3-D geometry of a car — a lower body, an inset raked **cabin**
//! (greenhouse) and four **wheels** — for the viewport and renders. The body
//! and cabin are extruded side-profiles; the wheels are cylinders on a lateral
//! axis. Built the same way as the launch-vehicle mesh: simple primitives
//! assembled into a recognizable silhouette.
//!
//! The car lies along **+X** (nose at `x = 0`, tail at `x = length`), is
//! centred on **Y**, and sits on the ground at **z = 0** with the wheels'
//! bottoms on the ground. Triangles are emitted double-sided so the render
//! reads correctly regardless of back-face culling.
//!
//! Honest scope: an *illustrative* procedural body (a clean coupe silhouette),
//! not a surfaced CAD model — exactly what the rocket mesh is for the launch
//! vehicle. A NURBS/surface body is a later step in the suite.

use std::f64::consts::TAU;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Outer dimensions of a car (metres).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarDimensions {
    /// Overall length (m).
    pub length: f64,
    /// Overall width (m).
    pub width: f64,
    /// Overall height to the roof (m).
    pub height: f64,
    /// Wheelbase — front-to-rear axle distance (m).
    pub wheelbase: f64,
    /// Track — left-to-right wheel-centre distance (m).
    pub track: f64,
    /// Wheel (tire) radius (m).
    pub wheel_radius: f64,
    /// Wheel (tire) width (m).
    pub wheel_width: f64,
    /// Ground clearance under the body floor (m).
    pub ground_clearance: f64,
}

impl Default for CarDimensions {
    /// A sporty coupe.
    fn default() -> Self {
        Self {
            length: 4.4,
            width: 1.9,
            height: 1.30,
            wheelbase: 2.6,
            track: 1.6,
            wheel_radius: 0.34,
            wheel_width: 0.26,
            ground_clearance: 0.12,
        }
    }
}

impl CarDimensions {
    /// A large luxury-EV sedan with Tesla Model-S-class proportions — a long
    /// wheelbase, a low wide stance and (paired with [`car_body_mesh`]) the
    /// fastback roofline.
    pub fn tesla() -> Self {
        Self {
            length: 4.97,
            width: 1.96,
            height: 1.45,
            wheelbase: 2.96,
            track: 1.66,
            wheel_radius: 0.355,
            wheel_width: 0.265,
            ground_clearance: 0.12,
        }
    }
}

/// Push a quad `(a, b, c, d)` as two triangles, both windings (double-sided).
fn push_quad(tris: &mut Vec<usize>, a: usize, b: usize, c: usize, d: usize) {
    tris.extend_from_slice(&[a, b, c, a, c, d, a, c, b, a, d, c]);
}

/// Extrude a closed `(x, z)` side profile across the full width `±y_half`,
/// double-sided (the two profile faces are the car's sides; the swept band is
/// its top/bottom/front/rear).
fn extrude(
    profile: &[(f64, f64)],
    y_half: f64,
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
) {
    let n = profile.len();
    let base = nodes.len();
    for &(x, z) in profile {
        nodes.push(Vector3::new(x, -y_half, z));
        nodes.push(Vector3::new(x, y_half, z));
    }
    // Swept lateral band around the closed loop.
    for i in 0..n {
        let j = (i + 1) % n;
        let (l0, r0) = (base + 2 * i, base + 2 * i + 1);
        let (l1, r1) = (base + 2 * j, base + 2 * j + 1);
        push_quad(tris, l0, r0, r1, l1);
    }
    // Two end caps (the flat sides), fan-triangulated from vertex 0,
    // double-sided.
    let (left0, right0) = (base, base + 1);
    for i in 1..n - 1 {
        let (la, lb) = (base + 2 * i, base + 2 * (i + 1));
        tris.extend_from_slice(&[left0, la, lb, left0, lb, la]);
        let (ra, rb) = (base + 2 * i + 1, base + 2 * (i + 1) + 1);
        tris.extend_from_slice(&[right0, ra, rb, right0, rb, ra]);
    }
}

/// Append a wheel — a capped cylinder whose axis is the lateral **Y** — at
/// `(cx, cy, cz)`, double-sided.
#[allow(clippy::too_many_arguments)]
fn add_wheel(
    cx: f64,
    cy: f64,
    cz: f64,
    radius: f64,
    half_width: f64,
    seg: usize,
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
) {
    let base = nodes.len();
    for &y in &[cy - half_width, cy + half_width] {
        for k in 0..seg {
            let a = k as f64 / seg as f64 * TAU;
            nodes.push(Vector3::new(
                cx + radius * a.cos(),
                y,
                cz + radius * a.sin(),
            ));
        }
    }
    let (r0, r1) = (base, base + seg);
    for k in 0..seg {
        let (a0, a1) = (r0 + k, r0 + (k + 1) % seg);
        let (b0, b1) = (r1 + k, r1 + (k + 1) % seg);
        push_quad(tris, a0, a1, b1, b0);
    }
    let c0 = nodes.len();
    nodes.push(Vector3::new(cx, cy - half_width, cz));
    let c1 = nodes.len();
    nodes.push(Vector3::new(cx, cy + half_width, cz));
    for k in 0..seg {
        let (a0, a1) = (r0 + k, r0 + (k + 1) % seg);
        tris.extend_from_slice(&[c0, a0, a1, c0, a1, a0]);
        let (b0, b1) = (r1 + k, r1 + (k + 1) % seg);
        tris.extend_from_slice(&[c1, b0, b1, c1, b1, b0]);
    }
}

/// Build the car **body + cabin** (no wheels) as one mesh.
pub fn car_body_mesh(dims: &CarDimensions) -> Mesh {
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let l = dims.length;
    let w_half = dims.width / 2.0;
    let floor = dims.ground_clearance + dims.wheel_radius * 0.15;
    let belt = floor + dims.height * 0.40;
    let roof = dims.height;
    let hood = floor + dims.height * 0.15;
    let tail = floor + dims.height * 0.34;
    let nose = floor + dims.height * 0.09; // low Tesla-style nose

    // Lower body side profile (closed loop): a low nose, a long flat hood, the
    // beltline, and a short tapered tail (subtle ducktail). Traced top edge
    // (nose → tail) then bottom edge (rear → front) to close the loop.
    let body = [
        (0.015 * l, nose),  // nose tip (low)
        (0.10 * l, hood),   // hood front
        (0.30 * l, belt),   // cowl / windshield base
        (0.82 * l, belt),   // rear beltline
        (l, tail),          // tail top
        (l, floor + 0.02),  // rear-bottom
        (0.015 * l, floor), // front-bottom (long front overhang)
    ];
    extrude(&body, w_half * 0.97, &mut nodes, &mut tris);

    // Inset fastback cabin (greenhouse): a fast windshield, a gently crowned
    // roof, and a long sloping backlight to a short decklid — the Model-S/3
    // silhouette. Extra points keep the roofline reading smooth, not faceted.
    let dz = roof - belt;
    let cabin = [
        (0.30 * l, belt),             // windshield base
        (0.40 * l, belt + 0.66 * dz), // windshield mid
        (0.485 * l, roof),            // roof front
        (0.575 * l, roof),            // roof crown (rear)
        (0.72 * l, belt + 0.50 * dz), // backlight mid (long fastback)
        (0.85 * l, belt + 0.04 * dz), // backlight base / decklid
    ];
    extrude(&cabin, w_half * 0.82, &mut nodes, &mut tris);

    finish(nodes, tris, "valenx-car-body")
}

/// Build the four **wheels** as one mesh.
pub fn car_wheel_mesh(dims: &CarDimensions) -> Mesh {
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let seg = 28usize;
    let x_f = 0.5 * (dims.length - dims.wheelbase);
    let x_r = x_f + dims.wheelbase;
    let ty = dims.track / 2.0;
    let hw = dims.wheel_width / 2.0;
    let r = dims.wheel_radius;
    for &x in &[x_f, x_r] {
        for &sy in &[-1.0_f64, 1.0] {
            add_wheel(x, sy * ty, r, r, hw, seg, &mut nodes, &mut tris);
        }
    }
    finish(nodes, tris, "valenx-car-wheels")
}

/// Build the whole car (body + cabin + wheels) as one mesh — for the viewport.
pub fn car_mesh(dims: &CarDimensions) -> Mesh {
    let body = car_body_mesh(dims);
    let wheels = car_wheel_mesh(dims);
    let mut nodes = body.nodes.clone();
    let mut tris: Vec<usize> = body.element_blocks[0]
        .connectivity
        .iter()
        .map(|&i| i as usize)
        .collect();
    let base = nodes.len();
    nodes.extend(wheels.nodes.iter().copied());
    tris.extend(
        wheels.element_blocks[0]
            .connectivity
            .iter()
            .map(|&i| i as usize + base),
    );
    finish(nodes, tris, "valenx-car")
}

/// Assemble a triangle soup into a single-block triangle mesh.
fn finish(nodes: Vec<Vector3<f64>>, tris: Vec<usize>, id: &str) -> Mesh {
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new(id);
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn car_mesh_is_nonempty_triangulated_and_in_bounds() {
        let d = CarDimensions::default();
        let m = car_mesh(&d);
        assert!(m.nodes.len() > 100, "node count {}", m.nodes.len());
        let conn = &m.element_blocks[0].connectivity;
        assert!(!conn.is_empty() && conn.len() % 3 == 0);
        assert!(conn.iter().all(|&i| (i as usize) < m.nodes.len()));
        // The car spans roughly 0..length in x and sits on the ground (z≈0
        // at the wheel bottoms, up to ~height at the roof).
        let zmin = m.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = m
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmin < 0.02, "wheels reach the ground, zmin={zmin}");
        assert!(
            (d.height - 0.05..d.height + 0.05).contains(&zmax),
            "roof near height, zmax={zmax}"
        );
        let xmax = m
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((d.length - 0.05..d.length + 0.05).contains(&xmax), "length");
    }

    #[test]
    fn tesla_preset_has_sedan_proportions_and_builds() {
        let d = CarDimensions::tesla();
        assert!((4.7..5.1).contains(&d.length), "Tesla length {}", d.length);
        assert!((1.85..2.0).contains(&d.width), "width {}", d.width);
        assert!((1.35..1.55).contains(&d.height), "height {}", d.height);
        assert!(
            d.wheelbase > 2.8 && d.wheelbase < d.length,
            "wheelbase {}",
            d.wheelbase
        );
        let m = car_mesh(&d);
        assert!(m.nodes.len() > 100 && !m.element_blocks[0].connectivity.is_empty());
        // The body spans the full length and the roof reaches the stated height.
        let xmax = m
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((d.length - 0.05..d.length + 0.05).contains(&xmax), "length");
        let zmax = m
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((d.height - 0.05..d.height + 0.05).contains(&zmax), "roof");
    }

    #[test]
    fn body_and_wheels_split_cleanly() {
        let d = CarDimensions::default();
        let body = car_body_mesh(&d);
        let wheels = car_wheel_mesh(&d);
        assert!(!body.nodes.is_empty() && !wheels.nodes.is_empty());
        // The whole car has at least as many nodes as body + wheels.
        let whole = car_mesh(&d);
        assert_eq!(whole.nodes.len(), body.nodes.len() + wheels.nodes.len());
    }

    /// Render a path-traced car PNG to TEMP — run with
    /// `cargo test -p valenx-vehicle dump_car_png -- --ignored --nocapture`.
    #[test]
    #[ignore = "writes a path-traced car PNG to TEMP"]
    fn dump_car_png() {
        use valenx_pathtrace::{render, vec3, PtCamera, PtMaterial, RenderParams, SceneBuilder};
        let d = CarDimensions::tesla();
        let res = 480u32;
        let camera = PtCamera::look_at(
            vec3(
                (d.length * 1.15) as f32,
                (-d.width * 2.0) as f32,
                (d.height * 1.6) as f32,
            ),
            vec3((d.length * 0.5) as f32, 0.0, (d.height * 0.42) as f32),
            vec3(0.0, 0.0, 1.0),
            40f32.to_radians(),
            res,
            res,
        );
        let mut b = SceneBuilder::new(camera);
        let paint = b.add_material(PtMaterial::metal([0.74, 0.07, 0.08], 0.16)); // glossy red
        let tire = b.add_material(PtMaterial::metal([0.05, 0.05, 0.06], 0.55)); // dark rubber
        let ground = b.add_material(PtMaterial::diffuse([0.42, 0.42, 0.45]));
        let key = b.add_material(PtMaterial::emissive([13.0, 12.4, 11.4]));
        let fill = b.add_material(PtMaterial::emissive([2.4, 2.4, 2.7]));
        b.add_mesh(&car_body_mesh(&d), paint);
        b.add_mesh(&car_wheel_mesh(&d), tire);
        let g = 40.0;
        b.add_quad(
            vec3(-g, -g, 0.0),
            vec3(g, -g, 0.0),
            vec3(g, g, 0.0),
            vec3(-g, g, 0.0),
            ground,
        );
        let (lx, ly, lz) = (12.0, 12.0, 16.0);
        b.add_quad(
            vec3(-lx, -ly, lz),
            vec3(lx, -ly, lz),
            vec3(lx, ly, lz),
            vec3(-lx, ly, lz),
            key,
        );
        b.add_quad(
            vec3(-30.0, 18.0, -2.0),
            vec3(30.0, 18.0, -2.0),
            vec3(30.0, 18.0, 22.0),
            vec3(-30.0, 18.0, 22.0),
            fill,
        );
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 200,
            max_depth: 6,
            seed: 0x5eed,
            exposure: 1.1,
        };
        let ldr = render(&scene, &params).expect("render").to_ldr(1.1);
        let path = std::env::temp_dir().join("valenx_car.png");
        let file = std::fs::File::create(&path).expect("create png");
        let mut enc = png::Encoder::new(
            std::io::BufWriter::new(file),
            ldr.width as u32,
            ldr.height as u32,
        );
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut w = enc.write_header().expect("header");
        w.write_image_data(&ldr.pixels).expect("data");
        w.finish().expect("finish");
        println!("WROTE {}", path.display());
    }
}

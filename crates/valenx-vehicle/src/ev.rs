//! Parametric **electric-vehicle assembly** — a Tesla-class "skateboard" EV
//! generated as a real **component tree**: structure, a thousands-cell battery
//! pack, dual drive units, four-corner suspension, braking, wheels, body
//! panels, glass, interior and lighting. Every component is a named, placed
//! part in a category, so the assembly has a genuine **bill-of-materials
//! count** — in the thousands, dominated (as a real EV is) by the battery
//! cells.
//!
//! The car lies along **+X** (nose at `x = 0`), centred on **Y**, on the ground
//! at **z = 0**.
//!
//! Honest scope: the components are **procedural / parametric placeholders**
//! sized realistically (a cell is a small cylinder, a motor a cylinder block,
//! a panel a thin box) — a genuine multi-thousand-component parametric
//! assembly, **not** production surfaced (Class-A) CAD. It captures the
//! *architecture and part count* of a modern EV, not its manufacturing detail.

use std::f64::consts::TAU;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::geometry::CarDimensions;

/// Component category — drives the bill-of-materials breakdown and which parts
/// a render includes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Category {
    /// Frame rails, crossmembers, subframes, floor.
    Structure,
    /// Battery pack: enclosure, modules and cells.
    Battery,
    /// Drive units: motors, inverters, gearboxes, halfshafts.
    Drivetrain,
    /// Suspension arms, springs, dampers, uprights, hubs.
    Suspension,
    /// Brake discs, calipers, pads.
    Braking,
    /// Wheels (rims + tires).
    Wheels,
    /// Exterior body panels.
    Body,
    /// Windows / glazing.
    Glass,
    /// Seats, dash, console, steering.
    Interior,
    /// Head- and tail-lights.
    Lighting,
}

impl Category {
    /// Short human label.
    pub fn label(self) -> &'static str {
        match self {
            Category::Structure => "structure",
            Category::Battery => "battery",
            Category::Drivetrain => "drivetrain",
            Category::Suspension => "suspension",
            Category::Braking => "braking",
            Category::Wheels => "wheels",
            Category::Body => "body",
            Category::Glass => "glass",
            Category::Interior => "interior",
            Category::Lighting => "lighting",
        }
    }
}

/// The primitive a component is built from (axis: 0 = X, 1 = Y, 2 = Z).
#[derive(Clone, Copy)]
enum Shape {
    Box { lx: f64, ly: f64, lz: f64 },
    Cyl { radius: f64, length: f64, axis: u8 },
}

/// One placed component in the assembly.
#[derive(Clone)]
pub struct Component {
    /// Human name (e.g. `battery-cell-1234`).
    pub name: String,
    /// Category for the BOM breakdown.
    pub category: Category,
    shape: Shape,
    center: Vector3<f64>,
}

impl Component {
    fn instantiate(&self, nodes: &mut Vec<Vector3<f64>>, tris: &mut Vec<usize>) {
        match self.shape {
            Shape::Box { lx, ly, lz } => add_box(self.center, lx, ly, lz, nodes, tris),
            Shape::Cyl {
                radius,
                length,
                axis,
            } => add_cyl(self.center, radius, length, axis, 8, nodes, tris),
        }
    }
}

/// Inputs to [`build_ev`].
#[derive(Debug, Clone, Copy)]
pub struct EvSpec {
    /// Overall car dimensions.
    pub dimensions: CarDimensions,
    /// Battery cells along the car's length.
    pub battery_cells_long: usize,
    /// Battery cells across the car's width.
    pub battery_cells_wide: usize,
}

impl Default for EvSpec {
    /// A large premium EV sedan with a ~2000-cell floor pack.
    fn default() -> Self {
        Self {
            dimensions: CarDimensions {
                length: 4.95,
                width: 1.96,
                height: 1.44,
                wheelbase: 2.96,
                track: 1.66,
                wheel_radius: 0.36,
                wheel_width: 0.28,
                ground_clearance: 0.14,
            },
            battery_cells_long: 64,
            battery_cells_wide: 32,
        }
    }
}

/// A generated electric-vehicle component assembly.
pub struct EvAssembly {
    /// Every placed component.
    pub components: Vec<Component>,
    /// The dimensions it was built from.
    pub dimensions: CarDimensions,
}

impl EvAssembly {
    /// Total component (BOM) count.
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Component counts per category, in the [`Category`] declaration order
    /// that appears in the assembly.
    pub fn count_by_category(&self) -> Vec<(Category, usize)> {
        let cats = [
            Category::Structure,
            Category::Battery,
            Category::Drivetrain,
            Category::Suspension,
            Category::Braking,
            Category::Wheels,
            Category::Body,
            Category::Glass,
            Category::Interior,
            Category::Lighting,
        ];
        cats.iter()
            .filter_map(|&c| {
                let n = self.components.iter().filter(|x| x.category == c).count();
                (n > 0).then_some((c, n))
            })
            .collect()
    }

    /// Tessellate the components whose category passes `include` into one mesh.
    pub fn to_mesh(&self, include: impl Fn(Category) -> bool) -> Mesh {
        let mut nodes: Vec<Vector3<f64>> = Vec::new();
        let mut tris: Vec<usize> = Vec::new();
        for c in self.components.iter().filter(|c| include(c.category)) {
            c.instantiate(&mut nodes, &mut tris);
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = tris.iter().map(|&i| i as u32).collect();
        let mut mesh = Mesh::new("valenx-ev");
        mesh.nodes = nodes;
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        mesh
    }

    /// The exterior mesh (body, glass, wheels, lighting).
    pub fn exterior_mesh(&self) -> Mesh {
        self.to_mesh(|c| {
            matches!(
                c,
                Category::Body | Category::Glass | Category::Wheels | Category::Lighting
            )
        })
    }

    /// The cutaway mesh — the internals with the body/glass/interior removed,
    /// so the skateboard chassis, battery pack and drivetrain are visible.
    pub fn cutaway_mesh(&self) -> Mesh {
        self.to_mesh(|c| !matches!(c, Category::Body | Category::Glass | Category::Interior))
    }
}

/// Generate a Tesla-class EV with a sensible default spec.
pub fn tesla_class() -> EvAssembly {
    build_ev(&EvSpec::default())
}

/// Build the parametric EV component assembly.
pub fn build_ev(spec: &EvSpec) -> EvAssembly {
    let d = spec.dimensions;
    let mut cs: Vec<Component> = Vec::new();
    let mut push = |name: String, cat: Category, shape: Shape, center: Vector3<f64>| {
        cs.push(Component {
            name,
            category: cat,
            shape,
            center,
        });
    };
    let l = d.length;
    let w = d.width;
    let x_f = 0.5 * (l - d.wheelbase);
    let x_r = x_f + d.wheelbase;
    let floor = d.ground_clearance + 0.06;

    // ── Structure: skateboard frame. ──────────────────────────────────────
    for (i, sy) in [-1.0, 1.0].into_iter().enumerate() {
        push(
            format!("frame-rail-{i}"),
            Category::Structure,
            Shape::Box {
                lx: l * 0.92,
                ly: 0.08,
                lz: 0.10,
            },
            Vector3::new(l * 0.5, sy * w * 0.42, floor),
        );
    }
    for i in 0..6 {
        let x = l * (0.15 + 0.14 * i as f64);
        push(
            format!("crossmember-{i}"),
            Category::Structure,
            Shape::Box {
                lx: 0.06,
                ly: w * 0.86,
                lz: 0.08,
            },
            Vector3::new(x, 0.0, floor),
        );
    }
    for (name, x) in [("front-subframe", x_f), ("rear-subframe", x_r)] {
        push(
            name.to_string(),
            Category::Structure,
            Shape::Box {
                lx: 0.5,
                ly: w * 0.8,
                lz: 0.12,
            },
            Vector3::new(x, 0.0, floor),
        );
    }
    push(
        "floor-pan".to_string(),
        Category::Structure,
        Shape::Box {
            lx: d.wheelbase,
            ly: w * 0.86,
            lz: 0.02,
        },
        Vector3::new(l * 0.5, 0.0, floor + 0.02),
    );

    // ── Battery: enclosure + the thousands-cell floor pack. ───────────────
    push(
        "battery-enclosure".to_string(),
        Category::Battery,
        Shape::Box {
            lx: d.wheelbase * 0.98,
            ly: w * 0.84,
            lz: 0.14,
        },
        Vector3::new(l * 0.5, 0.0, floor + 0.12),
    );
    let (nx, ny) = (
        spec.battery_cells_long.max(1),
        spec.battery_cells_wide.max(1),
    );
    let x0 = x_f + 0.10;
    let x1 = x_r - 0.10;
    let y_span = w * 0.74;
    let cell_z = floor + 0.12;
    for ix in 0..nx {
        for iy in 0..ny {
            let fx = if nx > 1 {
                ix as f64 / (nx - 1) as f64
            } else {
                0.5
            };
            let fy = if ny > 1 {
                iy as f64 / (ny - 1) as f64
            } else {
                0.5
            };
            let x = x0 + fx * (x1 - x0);
            let y = -0.5 * y_span + fy * y_span;
            push(
                format!("battery-cell-{ix}-{iy}"),
                Category::Battery,
                Shape::Cyl {
                    radius: 0.023,
                    length: 0.10,
                    axis: 2,
                },
                Vector3::new(x, y, cell_z),
            );
        }
    }
    // Module straps every 8 cells along the length.
    for m in 0..(nx / 8).max(1) {
        let x = x0 + (m as f64 + 0.5) * (x1 - x0) / (nx as f64 / 8.0).max(1.0);
        push(
            format!("battery-module-strap-{m}"),
            Category::Battery,
            Shape::Box {
                lx: 0.02,
                ly: y_span,
                lz: 0.02,
            },
            Vector3::new(x, 0.0, cell_z + 0.06),
        );
    }

    // ── Drivetrain: front + rear drive units. ─────────────────────────────
    for (tag, x) in [("front", x_f), ("rear", x_r)] {
        push(
            format!("{tag}-motor"),
            Category::Drivetrain,
            Shape::Cyl {
                radius: 0.13,
                length: 0.34,
                axis: 1,
            },
            Vector3::new(x, 0.0, floor + 0.10),
        );
        push(
            format!("{tag}-inverter"),
            Category::Drivetrain,
            Shape::Box {
                lx: 0.22,
                ly: 0.26,
                lz: 0.14,
            },
            Vector3::new(x, 0.0, floor + 0.30),
        );
        push(
            format!("{tag}-gearbox"),
            Category::Drivetrain,
            Shape::Box {
                lx: 0.18,
                ly: 0.20,
                lz: 0.18,
            },
            Vector3::new(x + 0.18, 0.0, floor + 0.10),
        );
        for sy in [-1.0, 1.0] {
            push(
                format!("{tag}-halfshaft"),
                Category::Drivetrain,
                Shape::Cyl {
                    radius: 0.025,
                    length: d.track * 0.42,
                    axis: 1,
                },
                Vector3::new(x, sy * d.track * 0.25, floor + 0.10),
            );
        }
    }

    // ── Suspension, braking, wheels — four corners. ───────────────────────
    let ty = d.track / 2.0;
    let r = d.wheel_radius;
    for (ax_tag, x) in [("front", x_f), ("rear", x_r)] {
        for sy in [-1.0, 1.0] {
            let side = if sy < 0.0 { "L" } else { "R" };
            let corner = format!("{ax_tag}-{side}");
            let hub = Vector3::new(x, sy * ty, r);
            for (part, dz, lx, ly, lz) in [
                ("upper-arm", 0.18, 0.30, 0.06, 0.04),
                ("lower-arm", 0.02, 0.34, 0.06, 0.04),
                ("upright", 0.10, 0.08, 0.06, 0.22),
                ("tie-rod", 0.06, 0.24, 0.03, 0.03),
            ] {
                push(
                    format!("{corner}-{part}"),
                    Category::Suspension,
                    Shape::Box { lx, ly, lz },
                    Vector3::new(x, sy * (ty - 0.18), r + dz),
                );
            }
            push(
                format!("{corner}-spring"),
                Category::Suspension,
                Shape::Cyl {
                    radius: 0.06,
                    length: 0.30,
                    axis: 2,
                },
                Vector3::new(x, sy * (ty - 0.22), r + 0.18),
            );
            push(
                format!("{corner}-damper"),
                Category::Suspension,
                Shape::Cyl {
                    radius: 0.03,
                    length: 0.34,
                    axis: 2,
                },
                Vector3::new(x, sy * (ty - 0.22), r + 0.18),
            );
            push(
                format!("{corner}-hub"),
                Category::Suspension,
                Shape::Cyl {
                    radius: 0.05,
                    length: 0.08,
                    axis: 1,
                },
                hub,
            );
            // Braking.
            push(
                format!("{corner}-brake-disc"),
                Category::Braking,
                Shape::Cyl {
                    radius: r * 0.62,
                    length: 0.03,
                    axis: 1,
                },
                hub,
            );
            push(
                format!("{corner}-caliper"),
                Category::Braking,
                Shape::Box {
                    lx: 0.10,
                    ly: 0.06,
                    lz: 0.14,
                },
                Vector3::new(x - 0.10, sy * ty, r),
            );
            push(
                format!("{corner}-brake-pads"),
                Category::Braking,
                Shape::Box {
                    lx: 0.06,
                    ly: 0.04,
                    lz: 0.10,
                },
                Vector3::new(x - 0.08, sy * ty, r),
            );
            // Wheel (rim + tire).
            push(
                format!("{corner}-rim"),
                Category::Wheels,
                Shape::Cyl {
                    radius: r * 0.62,
                    length: d.wheel_width,
                    axis: 1,
                },
                hub,
            );
            push(
                format!("{corner}-tire"),
                Category::Wheels,
                Shape::Cyl {
                    radius: r,
                    length: d.wheel_width,
                    axis: 1,
                },
                hub,
            );
        }
    }

    // ── Body panels. ──────────────────────────────────────────────────────
    let belt = floor + d.height * 0.34;
    for (name, x, lx, lz, z) in [
        ("hood", 0.16 * l, 0.7, 0.05, belt),
        ("front-trunk-lid", 0.10 * l, 0.4, 0.05, belt - 0.02),
        ("roof", 0.55 * l, 1.5, 0.05, d.height - 0.02),
        ("trunk-lid", 0.88 * l, 0.7, 0.05, belt + 0.06),
        ("front-bumper", 0.02 * l, 0.18, 0.4, floor + 0.25),
        ("rear-bumper", 0.99 * l, 0.18, 0.4, floor + 0.25),
    ] {
        push(
            name.to_string(),
            Category::Body,
            Shape::Box {
                lx,
                ly: w * 0.94,
                lz,
            },
            Vector3::new(x, 0.0, z),
        );
    }
    for (i, sy) in [-1.0, 1.0].into_iter().enumerate() {
        for (part, x, lx) in [
            ("front-door", 0.40 * l, 0.7),
            ("rear-door", 0.62 * l, 0.7),
            ("front-fender", 0.20 * l, 0.5),
            ("rear-quarter", 0.85 * l, 0.5),
        ] {
            push(
                format!("{part}-{}", if i == 0 { "L" } else { "R" }),
                Category::Body,
                Shape::Box {
                    lx,
                    ly: 0.04,
                    lz: 0.5,
                },
                Vector3::new(x, sy * w * 0.46, belt - 0.05),
            );
        }
    }

    // ── Glass. ────────────────────────────────────────────────────────────
    for (name, x, lx) in [
        ("windshield", 0.42 * l, 0.45),
        ("rear-window", 0.74 * l, 0.40),
    ] {
        push(
            name.to_string(),
            Category::Glass,
            Shape::Box {
                lx,
                ly: w * 0.82,
                lz: 0.02,
            },
            Vector3::new(x, 0.0, d.height - 0.10),
        );
    }
    for (i, sy) in [-1.0, 1.0].into_iter().enumerate() {
        for (part, x) in [
            ("front-side-window", 0.46 * l),
            ("rear-side-window", 0.66 * l),
        ] {
            push(
                format!("{part}-{}", if i == 0 { "L" } else { "R" }),
                Category::Glass,
                Shape::Box {
                    lx: 0.4,
                    ly: 0.02,
                    lz: 0.3,
                },
                Vector3::new(x, sy * w * 0.45, d.height - 0.18),
            );
        }
    }

    // ── Interior. ─────────────────────────────────────────────────────────
    for (i, (x, sy)) in [
        (0.46 * l, -1.0),
        (0.46 * l, 1.0),
        (0.66 * l, -1.0),
        (0.66 * l, 0.0),
        (0.66 * l, 1.0),
    ]
    .into_iter()
    .enumerate()
    {
        push(
            format!("seat-{i}"),
            Category::Interior,
            Shape::Box {
                lx: 0.5,
                ly: 0.45,
                lz: 0.6,
            },
            Vector3::new(x, sy * w * 0.24, belt - 0.05),
        );
    }
    push(
        "dashboard".to_string(),
        Category::Interior,
        Shape::Box {
            lx: 0.18,
            ly: w * 0.8,
            lz: 0.2,
        },
        Vector3::new(0.36 * l, 0.0, belt),
    );
    push(
        "center-console".to_string(),
        Category::Interior,
        Shape::Box {
            lx: 0.9,
            ly: 0.2,
            lz: 0.15,
        },
        Vector3::new(0.52 * l, 0.0, belt - 0.10),
    );
    push(
        "steering-wheel".to_string(),
        Category::Interior,
        Shape::Cyl {
            radius: 0.18,
            length: 0.05,
            axis: 0,
        },
        Vector3::new(0.40 * l, -w * 0.24, belt + 0.05),
    );

    // ── Lighting. ─────────────────────────────────────────────────────────
    for (i, sy) in [-1.0, 1.0].into_iter().enumerate() {
        let side = if i == 0 { "L" } else { "R" };
        push(
            format!("headlight-{side}"),
            Category::Lighting,
            Shape::Box {
                lx: 0.1,
                ly: 0.3,
                lz: 0.12,
            },
            Vector3::new(0.04 * l, sy * w * 0.32, belt - 0.06),
        );
        push(
            format!("taillight-{side}"),
            Category::Lighting,
            Shape::Box {
                lx: 0.08,
                ly: 0.34,
                lz: 0.12,
            },
            Vector3::new(0.98 * l, sy * w * 0.32, belt - 0.02),
        );
    }

    EvAssembly {
        components: cs,
        dimensions: d,
    }
}

/// Append a double-sided axis-aligned box centred at `c`.
fn add_box(
    c: Vector3<f64>,
    lx: f64,
    ly: f64,
    lz: f64,
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
) {
    let (hx, hy, hz) = (lx * 0.5, ly * 0.5, lz * 0.5);
    let base = nodes.len();
    for &sx in &[-1.0, 1.0] {
        for &sy in &[-1.0, 1.0] {
            for &sz in &[-1.0, 1.0] {
                nodes.push(c + Vector3::new(sx * hx, sy * hy, sz * hz));
            }
        }
    }
    // Corner index = (sx,sy,sz) → 4*ix+2*iy+iz with 0=-,1=+.
    let faces = [
        [0, 1, 3, 2], // -x
        [4, 6, 7, 5], // +x
        [0, 4, 5, 1], // -y
        [2, 3, 7, 6], // +y
        [0, 2, 6, 4], // -z
        [1, 5, 7, 3], // +z
    ];
    for f in faces {
        let (a, b, cc, d) = (base + f[0], base + f[1], base + f[2], base + f[3]);
        tris.extend_from_slice(&[a, b, cc, a, cc, d, a, cc, b, a, d, cc]);
    }
}

/// Append a double-sided capped cylinder of `length` along `axis`
/// (0 = X, 1 = Y, 2 = Z), centred at `c`.
fn add_cyl(
    c: Vector3<f64>,
    radius: f64,
    length: f64,
    axis: u8,
    seg: usize,
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
) {
    let h = length * 0.5;
    // Build a ring in the plane perpendicular to `axis`.
    let ring = |t: f64, end: f64| -> Vector3<f64> {
        let (cu, su) = (radius * t.cos(), radius * t.sin());
        match axis {
            0 => c + Vector3::new(end, cu, su),
            1 => c + Vector3::new(cu, end, su),
            _ => c + Vector3::new(cu, su, end),
        }
    };
    let base = nodes.len();
    for &end in &[-h, h] {
        for k in 0..seg {
            nodes.push(ring(k as f64 / seg as f64 * TAU, end));
        }
    }
    let (r0, r1) = (base, base + seg);
    for k in 0..seg {
        let (a0, a1) = (r0 + k, r0 + (k + 1) % seg);
        let (b0, b1) = (r1 + k, r1 + (k + 1) % seg);
        tris.extend_from_slice(&[a0, a1, b1, a0, b1, b0, a0, b1, a1, a0, b0, b1]);
    }
    let c0 = nodes.len();
    nodes.push(match axis {
        0 => c + Vector3::new(-h, 0.0, 0.0),
        1 => c + Vector3::new(0.0, -h, 0.0),
        _ => c + Vector3::new(0.0, 0.0, -h),
    });
    let c1 = nodes.len();
    nodes.push(match axis {
        0 => c + Vector3::new(h, 0.0, 0.0),
        1 => c + Vector3::new(0.0, h, 0.0),
        _ => c + Vector3::new(0.0, 0.0, h),
    });
    for k in 0..seg {
        let (a0, a1) = (r0 + k, r0 + (k + 1) % seg);
        tris.extend_from_slice(&[c0, a0, a1, c0, a1, a0]);
        let (b0, b1) = (r1 + k, r1 + (k + 1) % seg);
        tris.extend_from_slice(&[c1, b0, b1, c1, b1, b0]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tesla_class_has_thousands_of_components() {
        let ev = tesla_class();
        assert!(
            ev.component_count() > 2_000,
            "BOM count {} should be in the thousands",
            ev.component_count()
        );
        // The battery dominates the count, as in a real EV.
        let bat = ev
            .components
            .iter()
            .filter(|c| c.category == Category::Battery)
            .count();
        assert!(bat > ev.component_count() / 2, "battery dominates ({bat})");
        // Every architecture category is represented.
        assert_eq!(ev.count_by_category().len(), 10);
    }

    #[test]
    fn meshes_build_and_cutaway_excludes_the_body() {
        let ev = tesla_class();
        let full = ev.to_mesh(|_| true);
        assert!(!full.nodes.is_empty());
        let cut = ev.cutaway_mesh();
        let ext = ev.exterior_mesh();
        // The cutaway (with the thousands-cell pack) has far more vertices
        // than the exterior shell.
        assert!(cut.nodes.len() > ext.nodes.len());
        assert!(!ext.element_blocks[0].connectivity.is_empty());
    }

    /// Render a path-traced EV cutaway PNG to TEMP — run with
    /// `cargo test -p valenx-vehicle dump_ev_png -- --ignored --nocapture`.
    #[test]
    #[ignore = "writes a path-traced EV cutaway PNG to TEMP"]
    fn dump_ev_png() {
        use valenx_pathtrace::{render, vec3, PtCamera, PtMaterial, RenderParams, SceneBuilder};
        let ev = tesla_class();
        let d = ev.dimensions;
        let res = 460u32;
        let camera = PtCamera::look_at(
            vec3(
                (d.length * 1.05) as f32,
                (-d.width * 1.9) as f32,
                (d.height * 2.0) as f32,
            ),
            vec3((d.length * 0.5) as f32, 0.0, (d.height * 0.3) as f32),
            vec3(0.0, 0.0, 1.0),
            42f32.to_radians(),
            res,
            res,
        );
        let mut b = SceneBuilder::new(camera);
        let steel = b.add_material(PtMaterial::metal([0.55, 0.56, 0.60], 0.30));
        let cell = b.add_material(PtMaterial::metal([0.80, 0.78, 0.40], 0.35)); // brass-y cells
        let ground = b.add_material(PtMaterial::diffuse([0.40, 0.40, 0.43]));
        let key = b.add_material(PtMaterial::emissive([13.0, 12.5, 11.5]));
        let fill = b.add_material(PtMaterial::emissive([2.6, 2.6, 2.9]));
        // Battery cells in one colour, everything else (chassis/drivetrain/
        // suspension/wheels) in steel.
        b.add_mesh(&ev.to_mesh(|c| c == Category::Battery), cell);
        b.add_mesh(
            &ev.to_mesh(|c| {
                matches!(
                    c,
                    Category::Structure
                        | Category::Drivetrain
                        | Category::Suspension
                        | Category::Braking
                        | Category::Wheels
                )
            }),
            steel,
        );
        let g = 40.0;
        b.add_quad(
            vec3(-g, -g, 0.0),
            vec3(g, -g, 0.0),
            vec3(g, g, 0.0),
            vec3(-g, g, 0.0),
            ground,
        );
        let (lx, ly, lz) = (14.0, 14.0, 18.0);
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
            samples_per_pixel: 96,
            max_depth: 5,
            seed: 0x5eed,
            exposure: 1.1,
        };
        let ldr = render(&scene, &params).expect("render").to_ldr(1.1);
        let path = std::env::temp_dir().join("valenx_ev_cutaway.png");
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
        println!(
            "WROTE {} ({} components)",
            path.display(),
            ev.component_count()
        );
    }
}

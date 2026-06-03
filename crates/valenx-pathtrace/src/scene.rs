//! The path-tracer scene — triangles, materials, the camera, and the
//! HDR environment — plus the [`crate::bvh::Bvh`] built over the
//! triangles.
//!
//! A [`Scene`] is the immutable input to [`crate::tracer::render`]: it
//! is assembled once (with [`SceneBuilder`]) and then read-only during
//! the render so the work parallelises trivially.

use valenx_render_bridge::environment::EnvironmentMap;
use valenx_render_bridge::material::Material;

use crate::bvh::Bvh;
use crate::geometry::Triangle;
use crate::light_tree::LightTree;
use crate::math::{vec3, Vec3};

/// A surface material as the path tracer sees it.
///
/// It wraps the shared [`valenx_render_bridge::Material`] (the
/// Cook-Torrance PBR parameters — base colour, metallic, roughness,
/// IOR) and adds an **emission** term so a triangle can also be a light
/// source. A path tracer needs no separate light list: emitters are
/// just triangles whose material glows, and the integrator picks them
/// up both by accident (a path that wanders onto one) and on purpose
/// (next-event estimation, [`crate::tracer`]).
///
/// An optional [`Subsurface`] block makes the material a **subsurface
/// scatterer** (skin, marble, wax). When set the renderer's SSS path
/// ([`crate::sss`]) takes over for that material — the surface still
/// has its Cook-Torrance specular lobe for the Fresnel-reflected part
/// of the incident energy, but the diffuse lobe is replaced by a
/// physically-based subsurface model.
#[derive(Clone, Debug)]
pub struct PtMaterial {
    /// The Cook-Torrance / Lambert reflectance parameters.
    pub pbr: Material,
    /// Linear-RGB radiance this surface emits, in W·m⁻²·sr⁻¹. The zero
    /// vector means a non-emitting surface.
    pub emission: Vec3,
    /// Optional subsurface-scattering parameters. `None` means the
    /// material is a pure surface scatterer (the default).
    pub subsurface: Option<Subsurface>,
}

/// Per-channel subsurface-scattering parameters — the input to the
/// [`crate::sss`] random-walk SSS model.
///
/// Two parameterisations are common: the artist-friendly
/// `(subsurface_color, scale)` (Disney / Burley) and the physics-side
/// `(scattering, absorption)` (PBRT v4). This struct holds the
/// physics-side coefficients directly; the helpful constructor
/// [`Subsurface::from_color_scale`] maps the artist parameters to them
/// via the standard `scattering = 1/(mean-free-path · scale)` /
/// `absorption = scattering · (1 − color)` rule.
#[derive(Clone, Copy, Debug)]
pub struct Subsurface {
    /// Per-channel **scattering** coefficient `σ_s`, in inverse world
    /// units — how often a random-walk step scatters inside the
    /// medium.
    pub scattering: Vec3,
    /// Per-channel **absorption** coefficient `σ_a`, in inverse world
    /// units — how often a random-walk step is absorbed.
    pub absorption: Vec3,
    /// Henyey-Greenstein phase asymmetry `g ∈ (−1, 1)`. `0` is
    /// isotropic; small positive (≈ 0.7 for skin) is forward-peaked.
    pub g: f32,
}

impl Subsurface {
    /// Build the physics parameters from the artist-friendly
    /// `(subsurface_color, scale)` pair.
    ///
    /// `color` is the **single-scattering albedo** (the colour the
    /// material *looks* in a thick blob — skin is pinkish-red), each
    /// channel in `[0, 1]`. `scale` is the inverse mean-free path in
    /// world units (a smaller `scale` makes a more translucent
    /// material; `scale = 100` would be a very thin layer, `1` a
    /// thick one).
    ///
    /// The mapping is the standard PBRT one:
    ///
    /// ```text
    ///   σ_s = scale · color
    ///   σ_a = scale · (1 − color)
    /// ```
    ///
    /// so the extinction `σ_t = σ_a + σ_s = scale` is colour-neutral
    /// and the albedo `σ_s / σ_t = color` is exactly the artist input.
    pub fn from_color_scale(color: [f32; 3], scale: f32) -> Subsurface {
        let c = Vec3::from_array(color);
        let s = scale.max(1e-6);
        Subsurface {
            scattering: c.scale(s),
            absorption: vec3(
                (1.0 - c.x).max(0.0),
                (1.0 - c.y).max(0.0),
                (1.0 - c.z).max(0.0),
            )
            .scale(s),
            g: 0.0,
        }
    }

    /// Per-channel extinction `σ_t = σ_a + σ_s`.
    #[inline]
    pub fn extinction(&self) -> Vec3 {
        self.scattering.add(self.absorption)
    }

    /// Per-channel single-scattering albedo `σ_s / σ_t`, the
    /// per-channel survival probability of a random-walk step.
    #[inline]
    pub fn albedo(&self) -> Vec3 {
        let ext = self.extinction();
        Vec3 {
            x: if ext.x > 0.0 { self.scattering.x / ext.x } else { 0.0 },
            y: if ext.y > 0.0 { self.scattering.y / ext.y } else { 0.0 },
            z: if ext.z > 0.0 { self.scattering.z / ext.z } else { 0.0 },
        }
    }
}

impl PtMaterial {
    /// A non-emitting material wrapping the given PBR parameters.
    pub fn surface(pbr: Material) -> PtMaterial {
        PtMaterial {
            pbr,
            emission: Vec3::ZERO,
            subsurface: None,
        }
    }

    /// A matte (Lambertian) diffuse surface of the given linear-RGB
    /// albedo — roughness 1, non-metallic.
    pub fn diffuse(albedo: [f32; 3]) -> PtMaterial {
        let mut m = Material::matte("pt-diffuse", albedo);
        m.roughness = 1.0;
        m.metallic = 0.0;
        PtMaterial::surface(m)
    }

    /// A glossy metal of the given specular tint and roughness.
    pub fn metal(tint: [f32; 3], roughness: f32) -> PtMaterial {
        let mut m = Material::polished_metal("pt-metal", tint);
        m.roughness = roughness.clamp(0.0, 1.0);
        PtMaterial::surface(m)
    }

    /// An emissive (light-source) material radiating `radiance`.
    ///
    /// The reflectance is set to near-black so an emitter does not also
    /// bounce a bright diffuse lobe — a physical area light is its
    /// emission, nothing more.
    pub fn emissive(radiance: [f32; 3]) -> PtMaterial {
        let mut m = Material::matte("pt-light", [0.0, 0.0, 0.0]);
        m.roughness = 1.0;
        m.metallic = 0.0;
        PtMaterial {
            pbr: m,
            emission: Vec3::from_array(radiance),
            subsurface: None,
        }
    }

    /// A **subsurface scatterer** with the artist-friendly
    /// `(subsurface_color, scale)` parameterisation — see
    /// [`Subsurface::from_color_scale`].
    ///
    /// `color` is what the bulk of the material reads as in a thick
    /// blob (the single-scattering albedo); `scale` is the inverse
    /// mean-free path in world units (small → very translucent, large
    /// → very thin SSS layer). The surface still gets its
    /// Cook-Torrance specular lobe from the same `Material` defaults
    /// (so a marble surface is glossy *and* SSS-translucent, as it
    /// should be).
    pub fn subsurface(color: [f32; 3], scale: f32) -> PtMaterial {
        let mut m = Material::matte("pt-subsurface", color);
        m.roughness = 1.0;
        m.metallic = 0.0;
        PtMaterial {
            pbr: m,
            emission: Vec3::ZERO,
            subsurface: Some(Subsurface::from_color_scale(color, scale)),
        }
    }

    /// True if this material emits any light.
    #[inline]
    pub fn is_emitter(&self) -> bool {
        self.emission.max_component() > 0.0
    }

    /// True if this material uses the subsurface-scattering shading
    /// path.
    #[inline]
    pub fn is_subsurface(&self) -> bool {
        self.subsurface.is_some()
    }
}

/// A pinhole camera reduced to what ray generation needs.
///
/// Built from a `valenx_render_bridge::Camera` (or directly) into a
/// precomputed ray-generation frame: the eye, the bottom-left corner of
/// the image plane, and the per-pixel horizontal / vertical spans. A
/// primary ray for pixel `(px, py)` is then a couple of multiply-adds
/// — see [`crate::tracer`].
#[derive(Clone, Debug)]
pub struct PtCamera {
    /// Eye position.
    pub eye: Vec3,
    /// World-space lower-left corner of the image plane (at unit
    /// distance from the eye).
    pub lower_left: Vec3,
    /// World-space vector spanning the full image width.
    pub horizontal: Vec3,
    /// World-space vector spanning the full image height.
    pub vertical: Vec3,
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
}

impl PtCamera {
    /// Build a camera from an eye, a look-at target, an up hint, the
    /// vertical field of view (radians), and the output resolution.
    ///
    /// The frame is the textbook look-at construction: `w` points from
    /// the target back to the eye, `u` is `up × w`, `v` is `w × u`. The
    /// image plane sits one unit in front of the eye, sized by
    /// `tan(fov/2)` and the aspect ratio.
    pub fn look_at(
        eye: Vec3,
        target: Vec3,
        up: Vec3,
        fov_v_rad: f32,
        width: u32,
        height: u32,
    ) -> PtCamera {
        let aspect = width.max(1) as f32 / height.max(1) as f32;
        let half_h = (fov_v_rad * 0.5).tan();
        let half_w = aspect * half_h;
        // Camera basis.
        let w = eye.sub(target).normalized().unwrap_or(Vec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        });
        let u = up.cross(w).normalized().unwrap_or(Vec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        });
        let v = w.cross(u);
        let horizontal = u.scale(2.0 * half_w);
        let vertical = v.scale(2.0 * half_h);
        // Lower-left = eye − half-width·u − half-height·v − w.
        let lower_left = eye
            .sub(horizontal.scale(0.5))
            .sub(vertical.scale(0.5))
            .sub(w);
        PtCamera {
            eye,
            lower_left,
            horizontal,
            vertical,
            width,
            height,
        }
    }
}

/// The complete, immutable scene a render reads.
///
/// Assemble one with [`SceneBuilder`]; pass it to
/// [`crate::tracer::render`].
#[derive(Clone)]
pub struct Scene {
    /// Every shading triangle in the scene.
    pub triangles: Vec<Triangle>,
    /// The material table; a triangle's `material` field indexes here.
    pub materials: Vec<PtMaterial>,
    /// Indices (into `triangles`) of every emitting triangle —
    /// precomputed so next-event estimation can sample a light without
    /// scanning the whole scene each bounce.
    pub emitters: Vec<u32>,
    /// The bounding-volume hierarchy over `triangles`.
    pub bvh: Bvh,
    /// The **light importance tree** over the emitter triangles. Used
    /// by next-event estimation in place of the uniform emitter pick —
    /// a sampling distribution proportional to each cluster's
    /// power × geometric importance toward the shading point, so a
    /// many-light scene's variance does not climb linearly with the
    /// light count. See [`crate::light_tree`].
    pub light_tree: LightTree,
    /// The camera.
    pub camera: PtCamera,
    /// The HDR environment — sampled as background radiance and as the
    /// IBL light for any ray that escapes the geometry.
    pub environment: EnvironmentMap,
}

impl Scene {
    /// The total emitted area-weighted power is not tracked; this is a
    /// convenience for tests / callers — the number of emitting
    /// triangles.
    pub fn emitter_count(&self) -> usize {
        self.emitters.len()
    }
}

/// Incremental builder for a [`Scene`].
///
/// Materials are registered first (each returns a stable index);
/// triangles reference a material by that index. [`SceneBuilder::build`]
/// finalises the scene — it computes the emitter list and builds the
/// BVH.
pub struct SceneBuilder {
    triangles: Vec<Triangle>,
    materials: Vec<PtMaterial>,
    camera: PtCamera,
    environment: EnvironmentMap,
}

impl SceneBuilder {
    /// Start a new scene with the given camera. The environment
    /// defaults to a dim uniform grey (so a scene with no explicit
    /// environment still renders something rather than pure black);
    /// override it with [`SceneBuilder::environment`].
    pub fn new(camera: PtCamera) -> SceneBuilder {
        SceneBuilder {
            triangles: Vec::new(),
            materials: Vec::new(),
            camera,
            environment: EnvironmentMap::uniform([0.05, 0.05, 0.05]),
        }
    }

    /// Set the HDR environment map.
    pub fn environment(mut self, env: EnvironmentMap) -> SceneBuilder {
        self.environment = env;
        self
    }

    /// Register a material; returns its index for triangles to use.
    pub fn add_material(&mut self, m: PtMaterial) -> usize {
        self.materials.push(m);
        self.materials.len() - 1
    }

    /// Add a triangle with explicit per-vertex normals.
    pub fn add_triangle(&mut self, verts: [Vec3; 3], normals: [Vec3; 3], material: usize) {
        self.triangles
            .push(Triangle::new(verts, normals, material));
    }

    /// Add a flat-shaded triangle (geometric normal on every vertex).
    pub fn add_flat_triangle(&mut self, verts: [Vec3; 3], material: usize) {
        self.triangles.push(Triangle::flat(verts, material));
    }

    /// Add an axis-aligned quad as two flat triangles. Vertices are
    /// given in order around the rectangle (`a → b → c → d`).
    pub fn add_quad(&mut self, a: Vec3, b: Vec3, c: Vec3, d: Vec3, material: usize) {
        self.triangles.push(Triangle::flat([a, b, c], material));
        self.triangles.push(Triangle::flat([a, c, d], material));
    }

    /// Append a [`valenx_mesh::Mesh`]'s `Tri3` blocks as triangles, all
    /// assigned `material`.
    ///
    /// Per-vertex normals are area-weighted from the surrounding faces
    /// so a coarse mesh still shades smoothly. Non-`Tri3` element
    /// blocks are ignored. Returns the number of triangles added.
    pub fn add_mesh(&mut self, mesh: &valenx_mesh::Mesh, material: usize) -> usize {
        use valenx_mesh::element::ElementType;
        // Area-weighted vertex normals.
        let mut vnormals = vec![Vec3::ZERO; mesh.nodes.len()];
        for block in &mesh.element_blocks {
            if block.element_type != ElementType::Tri3 {
                continue;
            }
            for tri in block.connectivity.chunks_exact(3) {
                let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                if a >= mesh.nodes.len() || b >= mesh.nodes.len() || c >= mesh.nodes.len() {
                    continue;
                }
                let va = node_to_vec3(mesh.nodes[a]);
                let vb = node_to_vec3(mesh.nodes[b]);
                let vc = node_to_vec3(mesh.nodes[c]);
                // The un-normalised cross product is already weighted
                // by twice the triangle area.
                let face = vb.sub(va).cross(vc.sub(va));
                vnormals[a] = vnormals[a].add(face);
                vnormals[b] = vnormals[b].add(face);
                vnormals[c] = vnormals[c].add(face);
            }
        }

        let mut added = 0;
        for block in &mesh.element_blocks {
            if block.element_type != ElementType::Tri3 {
                continue;
            }
            for tri in block.connectivity.chunks_exact(3) {
                let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                if a >= mesh.nodes.len() || b >= mesh.nodes.len() || c >= mesh.nodes.len() {
                    continue;
                }
                let va = node_to_vec3(mesh.nodes[a]);
                let vb = node_to_vec3(mesh.nodes[b]);
                let vc = node_to_vec3(mesh.nodes[c]);
                let geo = vb.sub(va).cross(vc.sub(va));
                let na = vnormals[a].normalized().unwrap_or(geo);
                let nb = vnormals[b].normalized().unwrap_or(geo);
                let nc = vnormals[c].normalized().unwrap_or(geo);
                self.triangles
                    .push(Triangle::new([va, vb, vc], [na, nb, nc], material));
                added += 1;
            }
        }
        added
    }

    /// Finalise the scene: collect the emitter list and build the BVH
    /// + light tree.
    pub fn build(self) -> Scene {
        // An emitter is any triangle whose material glows.
        let emitters: Vec<u32> = self
            .triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| {
                self.materials
                    .get(t.material)
                    .map(|m| m.is_emitter())
                    .unwrap_or(false)
            })
            .map(|(i, _)| i as u32)
            .collect();
        let bvh = Bvh::build(&self.triangles);
        let light_tree = LightTree::build(&self.triangles, &self.materials, &emitters);
        Scene {
            triangles: self.triangles,
            materials: self.materials,
            emitters,
            bvh,
            light_tree,
            camera: self.camera,
            environment: self.environment,
        }
    }
}

/// Convert a `valenx-mesh` node (an `f64` `nalgebra` vector) to the
/// renderer's `f32` [`Vec3`].
#[inline]
fn node_to_vec3(n: nalgebra::Vector3<f64>) -> Vec3 {
    Vec3 {
        x: n.x as f32,
        y: n.y as f32,
        z: n.z as f32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    fn test_camera() -> PtCamera {
        PtCamera::look_at(
            vec3(0.0, 0.0, 5.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            60f32.to_radians(),
            64,
            64,
        )
    }

    #[test]
    fn camera_frame_spans_the_image_plane() {
        let cam = test_camera();
        // The image-plane corner plus the full spans is the opposite
        // corner; its centre should sit on the eye→target axis.
        let centre = cam
            .lower_left
            .add(cam.horizontal.scale(0.5))
            .add(cam.vertical.scale(0.5));
        // Centre is one unit in front of the eye toward −Z (target is
        // the origin, eye at +5Z).
        let forward = centre.sub(cam.eye).normalized().unwrap();
        assert!(forward.z < -0.99, "camera should look toward −Z");
    }

    #[test]
    fn material_indices_are_stable() {
        let mut b = SceneBuilder::new(test_camera());
        let red = b.add_material(PtMaterial::diffuse([0.8, 0.1, 0.1]));
        let green = b.add_material(PtMaterial::diffuse([0.1, 0.8, 0.1]));
        assert_eq!(red, 0);
        assert_eq!(green, 1);
    }

    #[test]
    fn builder_collects_emitters() {
        let mut b = SceneBuilder::new(test_camera());
        let diffuse = b.add_material(PtMaterial::diffuse([0.7; 3]));
        let light = b.add_material(PtMaterial::emissive([5.0, 5.0, 5.0]));
        // One emitting triangle, one not.
        b.add_flat_triangle(
            [vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0)],
            diffuse,
        );
        b.add_flat_triangle(
            [vec3(0.0, 0.0, 2.0), vec3(1.0, 0.0, 2.0), vec3(0.0, 1.0, 2.0)],
            light,
        );
        let scene = b.build();
        assert_eq!(scene.emitter_count(), 1, "exactly one emitter triangle");
        // The emitter index must point at the glowing triangle.
        let e = scene.emitters[0] as usize;
        assert!(scene.materials[scene.triangles[e].material].is_emitter());
    }

    #[test]
    fn add_quad_produces_two_triangles() {
        let mut b = SceneBuilder::new(test_camera());
        let m = b.add_material(PtMaterial::diffuse([0.5; 3]));
        b.add_quad(
            vec3(0.0, 0.0, 0.0),
            vec3(1.0, 0.0, 0.0),
            vec3(1.0, 1.0, 0.0),
            vec3(0.0, 1.0, 0.0),
            m,
        );
        let scene = b.build();
        assert_eq!(scene.triangles.len(), 2, "a quad is two triangles");
    }

    #[test]
    fn add_mesh_imports_tri3_blocks_with_smooth_normals() {
        use nalgebra::Vector3 as NV;
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = valenx_mesh::Mesh::new("quad");
        mesh.nodes.push(NV::new(0.0, 0.0, 0.0));
        mesh.nodes.push(NV::new(1.0, 0.0, 0.0));
        mesh.nodes.push(NV::new(1.0, 1.0, 0.0));
        mesh.nodes.push(NV::new(0.0, 1.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 1, 2, 0, 2, 3],
        });
        let mut b = SceneBuilder::new(test_camera());
        let m = b.add_material(PtMaterial::diffuse([0.6; 3]));
        let added = b.add_mesh(&mesh, m);
        assert_eq!(added, 2, "two Tri3 elements imported");
        let scene = b.build();
        // The flat quad's interpolated normals should all point +Z.
        for t in &scene.triangles {
            assert!(t.n0.z > 0.99, "imported normal should face +Z");
        }
    }

    #[test]
    fn emissive_material_carries_no_diffuse_lobe() {
        // An emitter should not also bounce a bright diffuse colour.
        let light = PtMaterial::emissive([10.0, 10.0, 10.0]);
        assert!(light.is_emitter());
        assert_eq!(
            light.pbr.diffuse_color,
            [0.0, 0.0, 0.0],
            "emitter must not double as a diffuse reflector"
        );
    }
}

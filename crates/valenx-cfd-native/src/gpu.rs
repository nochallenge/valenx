//! Optional **GPU compute** acceleration for the dominant per-iteration
//! grid relaxation — the pressure-Poisson weighted-Jacobi sweep.
//!
//! # Why this module exists
//!
//! Profiling a SIMPLE run shows the cost is overwhelmingly the inner
//! pressure-correction Poisson solve ([`crate::linsolve`],
//! [`crate::multigrid`]): every outer iteration relaxes the five-point
//! stencil over the whole grid many times. On a fine grid that inner
//! sweep is the hot loop. It is also *embarrassingly parallel* in its
//! **weighted-Jacobi** form (unlike Gauss-Seidel/SOR, a Jacobi sweep
//! reads only the previous iterate, so every cell updates independently)
//! — exactly the shape a GPU eats for breakfast.
//!
//! This module moves that one sweep onto the GPU via a WGSL compute
//! shader, dispatched through `wgpu`. It is a faithful port of
//! [`crate::multigrid::weighted_jacobi_sweep`]: same five-point stencil,
//! same `xₙ₊₁ = xₙ + ω·(b + Σ aₙᵦ·xₙ,ₙᵦ)/aP − ω·xₙ` update, same
//! homogeneous-Neumann boundary handling (out-of-domain neighbours drop
//! out of the sum), same degenerate-cell guard (`aP ≈ 0` ⇒ leave alone).
//!
//! # This is an *optional, additive* path — not a replacement
//!
//! The CPU solvers in [`crate::linsolve`] / [`crate::multigrid`] remain
//! the default and only solver entry point. This module is gated behind
//! the `gpu` Cargo feature and is consulted only when a caller asks for
//! it *and* a GPU adapter is actually present. [`GpuJacobi::new`]
//! returns `None` when no adapter can be acquired (CI, headless boxes
//! with no compute device, a software-only environment), so the caller
//! transparently falls back to the CPU sweep. Nothing here can make the
//! solver *less* correct or *less* portable — it can only make the inner
//! relaxation faster where the hardware allows.
//!
//! # Precision
//!
//! WGSL's portable storage type is 32-bit `f32` (64-bit `f64` is not a
//! core WGSL type and is unavailable on most adapters), so the GPU sweep
//! runs in single precision. The CPU comparison reference
//! ([`jacobi_sweep_f32_cpu`]) therefore mirrors the sweep in `f32` too,
//! so a GPU-vs-CPU equality test compares like with like to within
//! single-precision round-off. The production `f64` solver is untouched;
//! this single-precision smoother is a perfectly standard choice for the
//! *inner* Poisson relaxation, whose result SIMPLE only ever consumes
//! approximately.

use crate::linsolve::PoissonCoeffs;

/// The WGSL compute shader: one weighted-Jacobi sweep of the five-point
/// pressure-Poisson stencil.
///
/// One invocation per pressure cell. Each cell reads its six per-cell
/// coefficients (`aP, aE, aW, aN, aS, b`, interleaved six-to-a-cell in a
/// single `coeffs` buffer) and the four in-domain neighbour values from
/// the *input* buffer, computes the weighted-Jacobi update, and writes it
/// to the *output* buffer. Because the read source (`x_in`) and the write
/// target (`x_out`) are distinct buffers, every cell is independent — the
/// defining property that makes Jacobi (not Gauss-Seidel) the
/// GPU-friendly smoother.
///
/// The six coefficients share one storage buffer (rather than six
/// separate buffers) so the bind group needs only **three** storage
/// bindings — within the conservative downlevel limit of four
/// (`max_storage_buffers_per_shader_stage`), keeping the kernel portable
/// to low-end and integrated adapters.
///
/// Boundary handling matches the CPU sweep exactly: an out-of-domain
/// neighbour contributes nothing (its coefficient was assembled to zero,
/// and the bounds check skips the read), i.e. a homogeneous-Neumann edge
/// on the pressure correction. A degenerate cell (`|aP| < 1e-30`) is left
/// unchanged.
pub const JACOBI_SWEEP_WGSL: &str = r#"
// Grid + relaxation parameters. std140-friendly layout (vec4 padding).
struct Params {
    nx    : u32,
    ny    : u32,
    omega : f32,
    _pad  : f32,
};

// Six coefficients per cell, interleaved: [aP, aE, aW, aN, aS, b].
const STRIDE : u32 = 6u;

@group(0) @binding(0) var<uniform>             params : Params;
@group(0) @binding(1) var<storage, read>       coeffs : array<f32>;
@group(0) @binding(2) var<storage, read>       x_in   : array<f32>;
@group(0) @binding(3) var<storage, read_write> x_out  : array<f32>;

@compute @workgroup_size(8, 8, 1)
fn jacobi_sweep(@builtin(global_invocation_id) gid : vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    if (i >= params.nx || j >= params.ny) {
        return;
    }
    let idx = i + j * params.nx;
    let c   = idx * STRIDE;

    let ap_c = coeffs[c + 0u];   // aP
    let old  = x_in[idx];
    // Degenerate cell (no stencil): leave it exactly as it was so the
    // output buffer is a complete, consistent next iterate.
    if (abs(ap_c) < 1e-30) {
        x_out[idx] = old;
        return;
    }

    var sum = coeffs[c + 5u];    // b
    if (i + 1u < params.nx) { sum = sum + coeffs[c + 1u] * x_in[idx + 1u]; }        // aE
    if (i >= 1u)            { sum = sum + coeffs[c + 2u] * x_in[idx - 1u]; }        // aW
    if (j + 1u < params.ny) { sum = sum + coeffs[c + 3u] * x_in[idx + params.nx]; } // aN
    if (j >= 1u)            { sum = sum + coeffs[c + 4u] * x_in[idx - params.nx]; } // aS

    let jac = sum / ap_c;
    x_out[idx] = old + params.omega * (jac - old);
}
"#;

/// A single weighted-Jacobi sweep of the five-point Poisson stencil, in
/// `f32`, on the CPU.
///
/// This is the *reference* against which the GPU sweep is validated: it
/// is the exact arithmetic [`JACOBI_SWEEP_WGSL`] performs, run in the
/// same single precision, so the two agree to within `f32` round-off
/// rather than to within the larger `f32`-vs-`f64` gap. It is also a
/// self-contained, always-available CPU smoother that the tests use to
/// confirm a field relaxes toward the analytic Poisson solution.
///
/// The coefficient fields and `x_in`/`x_out` are flat `nx·ny` slices in
/// row-major (x-fastest) order, matching [`crate::grid::Field`]. `x_in`
/// and `x_out` must be distinct slices (Jacobi reads the old iterate
/// while writing the new one).
///
/// # Panics
///
/// Panics if any slice length is not `nx·ny` — a programmer error in
/// buffer setup, not a runtime input.
#[allow(clippy::too_many_arguments)]
pub fn jacobi_sweep_f32_cpu(
    nx: usize,
    ny: usize,
    omega: f32,
    ap: &[f32],
    ae: &[f32],
    aw: &[f32],
    an: &[f32],
    a_s: &[f32],
    b: &[f32],
    x_in: &[f32],
    x_out: &mut [f32],
) {
    let n = nx * ny;
    for (label, len) in [
        ("ap", ap.len()),
        ("ae", ae.len()),
        ("aw", aw.len()),
        ("an", an.len()),
        ("a_s", a_s.len()),
        ("b", b.len()),
        ("x_in", x_in.len()),
        ("x_out", x_out.len()),
    ] {
        assert_eq!(len, n, "buffer `{label}` must have length nx*ny = {n}");
    }

    for j in 0..ny {
        for i in 0..nx {
            let idx = i + j * nx;
            let ap_c = ap[idx];
            let old = x_in[idx];
            if ap_c.abs() < 1e-30 {
                x_out[idx] = old;
                continue;
            }
            let mut sum = b[idx];
            if i + 1 < nx {
                sum += ae[idx] * x_in[idx + 1];
            }
            if i >= 1 {
                sum += aw[idx] * x_in[idx - 1];
            }
            if j + 1 < ny {
                sum += an[idx] * x_in[idx + nx];
            }
            if j >= 1 {
                sum += a_s[idx] * x_in[idx - nx];
            }
            let jac = sum / ap_c;
            x_out[idx] = old + omega * (jac - old);
        }
    }
}

/// Pack a [`PoissonCoeffs`] (whose fields are `f64` [`crate::grid::Field`])
/// into the flat `f32` coefficient slices the sweep consumes.
///
/// Returns `(ap, ae, aw, an, a_s, b)`, each `nx·ny` long, row-major.
fn coeffs_to_f32(coeffs: &PoissonCoeffs) -> CoeffArrays {
    let cast = |f: &crate::grid::Field| f.data.iter().map(|&v| v as f32).collect::<Vec<f32>>();
    CoeffArrays {
        nx: coeffs.nx,
        ny: coeffs.ny,
        ap: cast(&coeffs.ap),
        ae: cast(&coeffs.ae),
        aw: cast(&coeffs.aw),
        an: cast(&coeffs.an),
        a_s: cast(&coeffs.as_),
        b: cast(&coeffs.b),
    }
}

/// Flat `f32` view of a [`PoissonCoeffs`], ready for either the CPU
/// reference sweep or upload into GPU storage buffers.
#[derive(Clone, Debug)]
pub struct CoeffArrays {
    /// Cells along x.
    pub nx: usize,
    /// Cells along y.
    pub ny: usize,
    /// Diagonal `aP`, `nx·ny` row-major.
    pub ap: Vec<f32>,
    /// East coefficient `aE`.
    pub ae: Vec<f32>,
    /// West coefficient `aW`.
    pub aw: Vec<f32>,
    /// North coefficient `aN`.
    pub an: Vec<f32>,
    /// South coefficient `aS`.
    pub a_s: Vec<f32>,
    /// Source `b`.
    pub b: Vec<f32>,
}

impl CoeffArrays {
    /// Build the flat `f32` arrays from a [`PoissonCoeffs`].
    pub fn from_coeffs(coeffs: &PoissonCoeffs) -> CoeffArrays {
        coeffs_to_f32(coeffs)
    }

    /// Run one weighted-Jacobi sweep on the CPU (the GPU reference),
    /// reading `x_in` and writing `x_out`.
    pub fn jacobi_sweep_cpu(&self, omega: f32, x_in: &[f32], x_out: &mut [f32]) {
        jacobi_sweep_f32_cpu(
            self.nx, self.ny, omega, &self.ap, &self.ae, &self.aw, &self.an, &self.a_s, &self.b,
            x_in, x_out,
        );
    }

    /// The six coefficients interleaved six-to-a-cell — `[aP, aE, aW, aN,
    /// aS, b]` per cell, `6·nx·ny` long — the single packed layout the
    /// GPU `coeffs` storage buffer consumes.
    fn interleaved(&self) -> Vec<f32> {
        let n = self.nx * self.ny;
        let mut out = Vec::with_capacity(6 * n);
        for k in 0..n {
            out.push(self.ap[k]);
            out.push(self.ae[k]);
            out.push(self.aw[k]);
            out.push(self.an[k]);
            out.push(self.a_s[k]);
            out.push(self.b[k]);
        }
        out
    }
}

// --------------------------------------------------------------------
// The wgpu compute path.
// --------------------------------------------------------------------

/// Uniform block mirroring the WGSL `Params` struct. 16-byte aligned.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct ParamsUbo {
    nx: u32,
    ny: u32,
    omega: f32,
    _pad: f32,
}

/// A constructed GPU weighted-Jacobi smoother: a `wgpu` device, queue,
/// compute pipeline, and the storage buffers for one grid.
///
/// Acquiring this is the expensive part (adapter request, device
/// creation, shader compilation, buffer allocation); once built, each
/// [`GpuJacobi::sweep`] is a cheap buffer write + dispatch + readback.
/// The smoother is sized to one `nx × ny` grid at construction.
///
/// Construct with [`GpuJacobi::new`], which returns `None` when no GPU
/// adapter is available — the caller then falls back to the CPU sweep.
pub struct GpuJacobi {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    params_buf: wgpu::Buffer,
    /// The six per-cell coefficients interleaved into one storage buffer
    /// (`[aP,aE,aW,aN,aS,b]` per cell) so the bind group stays within the
    /// downlevel 4-storage-buffer limit.
    coeffs_buf: wgpu::Buffer,
    x_in_buf: wgpu::Buffer,
    x_out_buf: wgpu::Buffer,
    readback_buf: wgpu::Buffer,
    nx: usize,
    ny: usize,
}

impl GpuJacobi {
    /// Try to build a GPU smoother for the given coefficients.
    ///
    /// Requests a headless `wgpu` adapter (any backend, low-power
    /// preference — this is a compute job, no surface), creates a device
    /// and queue, compiles [`JACOBI_SWEEP_WGSL`], and uploads the
    /// coefficient buffers. Returns `None` if **no adapter is available**
    /// (the standard headless / CI / software-only case) so the caller
    /// can fall back to the CPU sweep without an error.
    pub fn new(coeffs: &PoissonCoeffs) -> Option<GpuJacobi> {
        let arrays = CoeffArrays::from_coeffs(coeffs);
        pollster::block_on(Self::new_async(&arrays))
    }

    async fn new_async(arrays: &CoeffArrays) -> Option<GpuJacobi> {
        use wgpu::util::DeviceExt;

        let nx = arrays.nx;
        let ny = arrays.ny;
        let n = nx * ny;
        if n == 0 {
            return None;
        }

        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        // No surface: this is a pure compute job.
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::LowPower,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await?;

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("valenx-cfd-jacobi"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::downlevel_defaults(),
                },
                None,
            )
            .await
            .ok()?;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("jacobi-sweep-wgsl"),
            source: wgpu::ShaderSource::Wgsl(JACOBI_SWEEP_WGSL.into()),
        });

        let storage_size = (n * std::mem::size_of::<f32>()) as wgpu::BufferAddress;

        let storage = |label: &str, contents: &[f32]| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: bytemuck::cast_slice(contents),
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
            })
        };

        // The six coefficient fields interleaved into one buffer keeps the
        // storage-binding count at three (coeffs, x_in, x_out) — under the
        // conservative downlevel limit of four.
        let coeffs_buf = storage("coeffs", &arrays.interleaved());
        // x_in / x_out start zeroed; the caller seeds x_in each sweep.
        let zeros = vec![0.0f32; n];
        let x_in_buf = storage("x_in", &zeros);
        let x_out_buf = storage("x_out", &zeros);

        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::bytes_of(&ParamsUbo {
                nx: nx as u32,
                ny: ny as u32,
                omega: 1.0,
                _pad: 0.0,
            }),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: storage_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Bind-group layout: binding 0 a uniform; 1..=7 read-only
        // storage; 8 read-write storage (the output).
        let entry = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("jacobi-bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                entry(1, true),  // coeffs (read-only)
                entry(2, true),  // x_in   (read-only)
                entry(3, false), // x_out  (read-write)
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("jacobi-pl"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("jacobi-pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "jacobi_sweep",
            compilation_options: Default::default(),
        });

        Some(GpuJacobi {
            device,
            queue,
            pipeline,
            bind_group_layout,
            params_buf,
            coeffs_buf,
            x_in_buf,
            x_out_buf,
            readback_buf,
            nx,
            ny,
        })
    }

    /// The grid dimensions this smoother was built for.
    pub fn dims(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Run one weighted-Jacobi sweep on the GPU and return the new
    /// iterate.
    ///
    /// `x_in` (length `nx·ny`) is uploaded as the previous iterate; the
    /// shader writes the updated field, which is read back and returned.
    /// The math is identical to [`CoeffArrays::jacobi_sweep_cpu`] /
    /// [`crate::multigrid::weighted_jacobi_sweep`].
    ///
    /// # Panics
    ///
    /// Panics if `x_in.len() != nx·ny`.
    pub fn sweep(&self, omega: f32, x_in: &[f32]) -> Vec<f32> {
        let n = self.nx * self.ny;
        assert_eq!(x_in.len(), n, "x_in must have length nx*ny = {n}");
        pollster::block_on(self.sweep_async(omega, x_in))
    }

    async fn sweep_async(&self, omega: f32, x_in: &[f32]) -> Vec<f32> {
        let n = self.nx * self.ny;

        // Refresh the params (omega may change between sweeps) and the
        // input iterate.
        self.queue.write_buffer(
            &self.params_buf,
            0,
            bytemuck::bytes_of(&ParamsUbo {
                nx: self.nx as u32,
                ny: self.ny as u32,
                omega,
                _pad: 0.0,
            }),
        );
        self.queue
            .write_buffer(&self.x_in_buf, 0, bytemuck::cast_slice(x_in));

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("jacobi-bg"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.coeffs_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.x_in_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: self.x_out_buf.as_entire_binding(),
                },
            ],
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("jacobi-encoder"),
            });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("jacobi-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            // 8×8 workgroups cover the grid; ceil-divide so the tail
            // cells are dispatched (the shader bounds-checks the excess).
            let gx = (self.nx as u32).div_ceil(8);
            let gy = (self.ny as u32).div_ceil(8);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        // Copy the result into the mappable readback buffer.
        encoder.copy_buffer_to_buffer(
            &self.x_out_buf,
            0,
            &self.readback_buf,
            0,
            (n * std::mem::size_of::<f32>()) as wgpu::BufferAddress,
        );
        self.queue.submit(Some(encoder.finish()));

        // Map, wait, read, unmap.
        let slice = self.readback_buf.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.device.poll(wgpu::Maintain::Wait);
        // The map result is delivered through the poll above.
        let _ = rx.recv();
        let out: Vec<f32> = {
            let data = slice.get_mapped_range();
            bytemuck::cast_slice::<u8, f32>(&data).to_vec()
        };
        self.readback_buf.unmap();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Field;

    /// Build the five-point Laplacian coefficients on an `n × n` grid
    /// with cell size `h`, homogeneous-Neumann boundaries (every interior
    /// neighbour gets `1/h²`, `aP` is their sum). Mirrors the helper in
    /// `linsolve`/`multigrid`.
    fn laplacian(n: usize, h: f64) -> PoissonCoeffs {
        let mut c = PoissonCoeffs::zeros(n, n);
        let w = 1.0 / (h * h);
        for j in 0..n {
            for i in 0..n {
                let mut ap = 0.0;
                if i + 1 < n {
                    c.ae.set(i, j, w);
                    ap += w;
                }
                if i > 0 {
                    c.aw.set(i, j, w);
                    ap += w;
                }
                if j + 1 < n {
                    c.an.set(i, j, w);
                    ap += w;
                }
                if j > 0 {
                    c.as_.set(i, j, w);
                    ap += w;
                }
                c.ap.set(i, j, ap);
            }
        }
        c
    }

    /// A smooth, exactly zero-mean target field φ*(x,y)=cos(πx)cos(πy),
    /// cell-centred, and the source b = A·φ* that the Laplacian maps it
    /// to. Returns `(target, coeffs_with_b)`.
    fn manufactured_problem(n: usize) -> (Vec<f32>, PoissonCoeffs) {
        let h = 1.0 / n as f64;
        let mut c = laplacian(n, h);
        let mut target = Field::zeros(n, n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * h;
                let y = (j as f64 + 0.5) * h;
                target.set(
                    i,
                    j,
                    (std::f64::consts::PI * x).cos() * (std::f64::consts::PI * y).cos(),
                );
            }
        }
        let mean: f64 = target.data.iter().sum::<f64>() / target.data.len() as f64;
        for v in target.data.iter_mut() {
            *v -= mean;
        }
        // b(i,j) = aP·t − Σ a_nb·t_nb.
        for j in 0..n {
            for i in 0..n {
                let mut nb = 0.0;
                if i + 1 < n {
                    nb += c.ae.at(i, j) * target.at(i + 1, j);
                }
                if i > 0 {
                    nb += c.aw.at(i, j) * target.at(i - 1, j);
                }
                if j + 1 < n {
                    nb += c.an.at(i, j) * target.at(i, j + 1);
                }
                if j > 0 {
                    nb += c.as_.at(i, j) * target.at(i, j - 1);
                }
                c.b.set(i, j, c.ap.at(i, j) * target.at(i, j) - nb);
            }
        }
        let target_f32 = target.data.iter().map(|&v| v as f32).collect();
        (target_f32, c)
    }

    /// Subtract the mean (the singular all-Neumann gauge) in place.
    fn pin_mean(x: &mut [f32]) {
        let mean = x.iter().sum::<f32>() / x.len() as f32;
        for v in x.iter_mut() {
            *v -= mean;
        }
    }

    #[test]
    fn cpu_reference_relaxes_toward_the_analytic_solution() {
        // The CPU f32 Jacobi reference is the heart of the correctness
        // gate: iterate it and the field must converge to the known
        // manufactured solution φ* = cos(πx)cos(πy).
        let n = 16;
        let (target, c) = manufactured_problem(n);
        let arrays = CoeffArrays::from_coeffs(&c);

        let mut a = vec![0.0f32; n * n];
        let mut b = vec![0.0f32; n * n];
        // Weighted Jacobi on the 5-point Laplacian needs ω ≤ 1; ω=0.8 is
        // a standard smoother choice. Many sweeps to actually converge.
        let omega = 0.8f32;
        let initial_err = max_err(&a, &target);
        for k in 0..20_000 {
            if k % 2 == 0 {
                arrays.jacobi_sweep_cpu(omega, &a, &mut b);
                pin_mean(&mut b);
            } else {
                arrays.jacobi_sweep_cpu(omega, &b, &mut a);
                pin_mean(&mut a);
            }
        }
        // After an even number of sweeps the latest iterate is in `a`.
        let err = max_err(&a, &target);
        assert!(
            err < 2e-2,
            "Jacobi did not relax to the analytic solution: err {err} (started {initial_err})"
        );
        assert!(
            err < initial_err,
            "relaxation must reduce the error: {initial_err} -> {err}"
        );
    }

    #[test]
    fn cpu_reference_matches_multigrid_jacobi_one_sweep() {
        // One f32 CPU sweep must reproduce the production f64
        // `weighted_jacobi_sweep` to single-precision tolerance — proof
        // the reference really is the same stencil/update as the solver.
        let n = 12;
        let (_t, c) = manufactured_problem(n);
        let arrays = CoeffArrays::from_coeffs(&c);

        // A non-trivial starting iterate.
        let mut x0 = Field::zeros(n, n);
        for j in 0..n {
            for i in 0..n {
                x0.set(i, j, ((i * 7 + j * 13) % 11) as f64 * 0.1 - 0.5);
            }
        }
        let omega = 0.8;

        // Production f64 sweep.
        let mut f64_out = x0.clone();
        crate::multigrid::weighted_jacobi_sweep(&c, &mut f64_out, omega);

        // f32 reference sweep on the same input.
        let x0_f32: Vec<f32> = x0.data.iter().map(|&v| v as f32).collect();
        let mut f32_out = vec![0.0f32; n * n];
        arrays.jacobi_sweep_cpu(omega as f32, &x0_f32, &mut f32_out);

        let mut max_d = 0.0f32;
        for k in 0..n * n {
            max_d = max_d.max((f32_out[k] - f64_out.data[k] as f32).abs());
        }
        assert!(
            max_d < 1e-4,
            "f32 reference sweep diverges from the f64 production sweep: {max_d}"
        );
    }

    #[test]
    fn gpu_pipeline_constructs_or_skips_gracefully() {
        // Constructing the smoother must never panic: it either yields a
        // working pipeline (adapter present) or returns None (headless /
        // CI). Either outcome passes — we are gating that the wgpu
        // plumbing assembles, not that a GPU exists.
        let n = 8;
        let (_t, c) = manufactured_problem(n);
        match GpuJacobi::new(&c) {
            Some(gpu) => {
                assert_eq!(gpu.dims(), (n, n));
            }
            None => {
                eprintln!("no GPU adapter available — pipeline-construction skipped (ok)");
            }
        }
    }

    #[test]
    fn gpu_sweep_matches_cpu_reference_when_adapter_present() {
        // The core GPU validation: if an adapter exists, one GPU sweep
        // must equal one CPU f32 reference sweep to f32 round-off. If no
        // adapter, skip gracefully (do NOT fail).
        let n = 24;
        let (_t, c) = manufactured_problem(n);
        let arrays = CoeffArrays::from_coeffs(&c);

        let Some(gpu) = GpuJacobi::new(&c) else {
            eprintln!("no GPU adapter — GPU-vs-CPU comparison skipped (ok)");
            return;
        };

        // A non-trivial input iterate.
        let mut x_in = vec![0.0f32; n * n];
        for (k, v) in x_in.iter_mut().enumerate() {
            *v = ((k % 13) as f32) * 0.07 - 0.3;
        }
        let omega = 0.85f32;

        let gpu_out = gpu.sweep(omega, &x_in);
        let mut cpu_out = vec![0.0f32; n * n];
        arrays.jacobi_sweep_cpu(omega, &x_in, &mut cpu_out);

        assert_eq!(gpu_out.len(), cpu_out.len());
        let mut max_d = 0.0f32;
        for k in 0..n * n {
            max_d = max_d.max((gpu_out[k] - cpu_out[k]).abs());
        }
        assert!(
            max_d < 1e-3,
            "GPU sweep disagrees with CPU reference by {max_d} (> f32 tolerance)"
        );
    }

    #[test]
    fn gpu_multi_sweep_relaxes_toward_solution_when_adapter_present() {
        // End-to-end on the GPU: iterating the GPU sweep must drive the
        // field toward the manufactured analytic solution, exactly as the
        // CPU does — confirms the GPU path is a real solver, not just a
        // single-sweep match. Adapter-gated.
        let n = 16;
        let (target, c) = manufactured_problem(n);

        let Some(gpu) = GpuJacobi::new(&c) else {
            eprintln!("no GPU adapter — GPU relaxation skipped (ok)");
            return;
        };

        let omega = 0.8f32;
        let mut x = vec![0.0f32; n * n];
        let initial_err = max_err(&x, &target);
        for _ in 0..4000 {
            x = gpu.sweep(omega, &x);
            pin_mean(&mut x);
        }
        let err = max_err(&x, &target);
        assert!(
            err < initial_err && err < 5e-2,
            "GPU relaxation did not converge: {initial_err} -> {err}"
        );
    }

    fn max_err(a: &[f32], b: &[f32]) -> f32 {
        a.iter()
            .zip(b.iter())
            .fold(0.0f32, |m, (&x, &y)| m.max((x - y).abs()))
    }
}

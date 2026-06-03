//! Headless CPU preview of the viewport ground grid.
//!
//! Renders the SAME analytic grid the GPU `grid_fs` shader draws —
//! using valenx-viz's already-tested `ray_from_screen` /
//! `intersect_ground_y0` for per-pixel world XZ, and a CPU port of the
//! shader's coverage / axis-line / distance-fade math — straight to a
//! PNG. Lets the grid's appearance be verified and tuned without a GPU
//! or the GUI (the GPU pipeline itself is validated separately at
//! pipeline-creation time).
//!
//! Run: `cargo run -p valenx-viz --example grid_preview`
//! Output: `target/grid_preview.png`

use valenx_viz::{intersect_ground_y0, nice_grid_spacing, ray_from_screen, OrbitCamera};

const W: u32 = 1200;
const H: u32 = 800;

fn frac(x: f32) -> f32 {
    x - x.floor()
}

/// Mirror of `grid_fs`'s `line_a` for one axis: coverage of a grid line
/// at the given coord, anti-aliased over one screen-pixel (`fw`).
fn line_cov(coord_over_s: f32, fw_over_s: f32) -> f32 {
    let d = (frac(coord_over_s - 0.5) - 0.5).abs() / fw_over_s.max(1e-6);
    1.0 - d.min(1.0)
}

fn world_xz(cam: &OrbitCamera, x: f32, y: f32) -> Option<[f32; 2]> {
    let r = ray_from_screen(cam, W as f32, H as f32, [x, y]);
    intersect_ground_y0(&r).map(|p| [p.x, p.z])
}

fn main() {
    let cam = OrbitCamera::default(); // az 45°, el 30°, dist 10, target origin
    let eye = cam.eye();
    let cam_xz = [eye.x, eye.z];
    let minor = nice_grid_spacing(cam.distance);
    let major = minor * 10.0;
    let fade_dist = cam.distance * 14.0;
    let bg = [0.10f32, 0.11, 0.13];

    let mut px = vec![0u8; (W * H * 4) as usize];
    for y in 0..H {
        for x in 0..W {
            let (xf, yf) = (x as f32 + 0.5, y as f32 + 0.5);
            let mut col = bg;
            if let Some(c) = world_xz(&cam, xf, yf) {
                // Screen-space derivatives of world XZ == WGSL fwidth().
                let cdx = world_xz(&cam, xf + 1.0, yf).unwrap_or(c);
                let cdy = world_xz(&cam, xf, yf + 1.0).unwrap_or(c);
                let fw = [
                    (cdx[0] - c[0]).abs() + (cdy[0] - c[0]).abs(),
                    (cdx[1] - c[1]).abs() + (cdy[1] - c[1]).abs(),
                ];
                let cov = |s: f32| {
                    line_cov(c[0] / s, fw[0] / s).max(line_cov(c[1] / s, fw[1] / s))
                };
                let mut a = (cov(minor) * 0.45).max(cov(major) * 0.85);
                let mut gcol = [0.33f32, 0.35, 0.40];
                let aw = [1.5 * fw[0], 1.5 * fw[1]];
                if c[1].abs() < aw[1] {
                    gcol = [0.85, 0.27, 0.27]; // X axis (z=0) red
                    a = a.max(0.9);
                }
                if c[0].abs() < aw[0] {
                    gcol = [0.32, 0.48, 0.95]; // Z axis (x=0) blue
                    a = a.max(0.9);
                }
                let dist = ((c[0] - cam_xz[0]).powi(2) + (c[1] - cam_xz[1]).powi(2)).sqrt();
                a *= (1.0 - dist / fade_dist).clamp(0.0, 1.0);
                for k in 0..3 {
                    col[k] = bg[k] * (1.0 - a) + gcol[k] * a;
                }
            }
            let i = ((y * W + x) * 4) as usize;
            px[i] = (col[0] * 255.0) as u8;
            px[i + 1] = (col[1] * 255.0) as u8;
            px[i + 2] = (col[2] * 255.0) as u8;
            px[i + 3] = 255;
        }
    }

    let path = "target/grid_preview.png";
    let file = std::fs::File::create(path).expect("create png");
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), W, H);
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header()
        .expect("png header")
        .write_image_data(&px)
        .expect("png data");
    println!("wrote {path} ({W}x{H}); spacing minor={minor} major={major}");
}

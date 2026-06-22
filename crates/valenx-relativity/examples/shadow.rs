//! Render black-hole shadows.
//!
//! Prints an ASCII preview of a Schwarzschild and a near-extremal Kerr shadow
//! (edge-on), and writes full-resolution PPM images next to the working
//! directory. Run with:
//!
//! ```text
//! cargo run -p valenx-relativity --example shadow
//! ```

use std::f64::consts::FRAC_PI_2;
use std::fs::File;
use std::io::{BufWriter, Write};

use valenx_relativity::{kerr, render_shadow, schwarzschild, KerrNewman, ShadowImage};

fn ascii(img: &ShadowImage) {
    for row in 0..img.height {
        let mut line = String::with_capacity(img.width);
        for col in 0..img.width {
            line.push(if img.is_shadow(col, row) { '#' } else { ' ' });
        }
        println!("{line}");
    }
}

fn write_ppm(img: &ShadowImage, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(w, "P6\n{} {}\n255\n", img.width, img.height)?;
    for row in 0..img.height {
        for col in 0..img.width {
            let px: [u8; 3] = if img.is_shadow(col, row) {
                [0, 0, 0] // shadow: black
            } else {
                [30, 30, 60] // sky: dark blue
            };
            w.write_all(&px)?;
        }
    }
    w.flush()
}

fn render(label: &str, bh: &KerrNewman, ppm: &str) {
    println!("\n=== {label} (edge-on) ===");
    // ASCII preview: wide aspect to offset tall character cells.
    let preview = render_shadow(bh, 300.0, FRAC_PI_2, 8.0, 70, 33).expect("ascii render");
    ascii(&preview);
    println!(
        "shadow covers {:.1}% of the frame",
        preview.shadow_fraction() * 100.0
    );
    // Full-resolution image to disk.
    let full = render_shadow(bh, 300.0, FRAC_PI_2, 8.0, 160, 160).expect("ppm render");
    match write_ppm(&full, ppm) {
        Ok(()) => println!("wrote {ppm}"),
        Err(e) => eprintln!("could not write {ppm}: {e}"),
    }
}

fn main() {
    let dir = std::env::temp_dir();
    let p1 = dir.join("schwarzschild_shadow.ppm");
    let p2 = dir.join("kerr_shadow.ppm");
    render(
        "Schwarzschild (a=0)",
        &schwarzschild(1.0),
        &p1.to_string_lossy(),
    );
    render("Kerr (a=0.99 M)", &kerr(1.0, 0.99), &p2.to_string_lossy());
}

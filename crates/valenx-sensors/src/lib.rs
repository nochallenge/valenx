//! # valenx-sensors — simulated sensor models + a small autonomy harness
//!
//! An **in-house, dependency-light** library of *simulated sensors* for
//! autonomy work — perception, state estimation, control, and reinforcement
//! learning — plus a tiny step-based harness that flies a kinematic vehicle past
//! those sensors and emits a synchronised [`SensorFrame`] per tick. Each sensor
//! is a parameterised struct with a sample/measure function, validates its
//! configuration at construction (**fail loud, never `NaN`**), and is pinned in
//! tests against analytic ground truth.
//!
//! ## The sensors
//!
//! * **LiDAR** ([`Lidar`]) — ray-casts a regular angular beam pattern against an
//!   analytic [`Scene`] of planes, spheres, and triangles (the triangle test is
//!   Möller–Trumbore), returning a range per beam over the field of view, or a
//!   **no-return** when a beam misses or the hit lies outside the range window.
//! * **Camera** ([`Camera`]) — a pinhole projection through an intrinsic matrix
//!   `K` with **Brown–Conrady** radial-tangential lens distortion; the inverse
//!   ([`Camera::undistort`]) lifts a pixel back to a normalised ray and
//!   round-trips with [`Camera::project`].
//! * **IMU** ([`Imu`]) — a strap-down accelerometer + gyroscope: the
//!   accelerometer reports **specific force** `Rᵀ·(a − g)` (so it reads `+1 g`
//!   at rest, `0` in free fall), the gyroscope the body angular rate, each with
//!   a constant bias and zero-mean Gaussian noise.
//! * **GPS** ([`Gps`]) — converts a local **ENU** offset about a geodetic datum
//!   to a WGS-84 latitude/longitude/altitude fix, with Gaussian position noise.
//!
//! ## The autonomy harness
//!
//! [`Harness`] carries a minimal kinematic [`VehicleState`] (position, velocity,
//! orientation, body rate in a local ENU world), an analytic [`Scene`], and a
//! set of attached sensors. Each [`Harness::step`] applies a [`Command`]
//! (body-frame acceleration + angular rate), integrates the state, and returns a
//! [`SensorFrame`] with every sensor's reading at the new time — the step API an
//! autonomy stack or RL environment polls. (The kinematic state is intentionally
//! distinct from `valenx-vehicle`'s force-based `Car`, which computes
//! *performance* numbers and has no time-stepped pose to attach sensors to.)
//!
//! ## Determinism — no `rand`, seeded SplitMix64
//!
//! Reproducibility is a first-class requirement: a sensor simulation that gives
//! different readings on every run is impossible to regression-test. This crate
//! therefore takes **no `rand` dependency**. All noise comes from a tiny
//! in-crate [`SplitMix64`] PRNG (the same deterministic, seeded generator used
//! in `valenx-uq`), with standard-normal draws via Box–Muller. Each sensor owns
//! its own seeded generator, so given the same seeds every reading and every
//! frame stream is bit-for-bit identical across runs and machines. The PRNG is
//! **not** used for any security purpose.
//!
//! ## Honesty / scope caveats
//!
//! These are **analytic / graphics-grade models, not hardware-calibrated
//! devices.** They reproduce the *geometry and first-order error structure* of
//! each sensor — exactly what you need to build and V&V an autonomy pipeline in
//! the loop (defense thrusts **M9 autonomy V&V** and **M1 UAS**) — but they are
//! deliberately not high-fidelity device models:
//!
//! * **LiDAR** is pure ranging geometry: no beam divergence, intensity /
//!   reflectivity, multi-echo, motion distortion during a sweep, or atmospheric
//!   attenuation, and the world is analytic surfaces, not a loaded mesh with an
//!   acceleration structure.
//! * **Camera** models intrinsics + Brown–Conrady distortion only — no
//!   photometry, no noise/blur/rolling-shutter, no image rendering; it maps 3-D
//!   points to pixels and back.
//! * **IMU** is the standard *ideal + constant-bias + white-noise* model — no
//!   scale-factor or axis-misalignment error, no random-walk bias drift, no
//!   temperature or `g`-sensitivity, no quantisation.
//! * **GPS** is geometry plus a Gaussian error — no satellite geometry / DOP, no
//!   multipath, no atmospheric delay, no clock or correlated-error model. The
//!   WGS-84 conversion itself is exact to floating point.
//! * **The harness** is *kinematic*, not a dynamics simulator: it integrates
//!   commanded body acceleration/rate directly (no mass, tire, contact, or
//!   aero forces — use `valenx-vehicle` / `valenx-mbd` for those).
//!
//! Calibrating any of these against a specific real device (or validating an
//! autonomy stack for deployment) is out of scope; this crate is the
//! **research/educational-grade, reproducible testbed** that work would build on.
//!
//! ## Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_sensors::{
//!     Harness, VehicleState, Command, Lidar, LidarConfig, Imu, ImuConfig,
//!     Scene, Plane,
//! };
//!
//! // A world with a wall 5 m ahead (its normal faces the sensor).
//! let mut scene = Scene::new();
//! scene.push_plane(Plane::new(Vector3::new(5.0, 0.0, 0.0), Vector3::new(-1.0, 0.0, 0.0)).unwrap());
//!
//! // A single forward LiDAR beam + a default IMU on a vehicle at the origin.
//! let lidar = Lidar::new(
//!     LidarConfig { azimuth_steps: 1, elevation_steps: 1, h_fov: 0.0, v_fov: 0.0,
//!                   min_range: 0.0, max_range: 100.0, range_noise_std: 0.0 },
//!     0,
//! ).unwrap();
//! let imu = Imu::new(ImuConfig::default(), 0).unwrap();
//!
//! let mut harness = Harness::new(VehicleState::default(), scene)
//!     .with_lidar(lidar)
//!     .with_imu(imu);
//!
//! // Step the sim 0.1 s while coasting and read the synchronised frame.
//! let frame = harness.step(&Command::coast(), 0.1).unwrap();
//!
//! // The forward beam ranges the wall at exactly 5 m...
//! let scan = frame.lidar.unwrap();
//! assert!((scan.beams[0].range.unwrap() - 5.0).abs() < 1e-9);
//! // ...and the at-rest IMU reads +1 g on its up (+z) axis.
//! assert!((frame.imu.unwrap().specific_force.z - 9.806_65).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]

pub mod camera;
pub mod gps;
pub mod harness;
pub mod imu;
pub mod lidar;
pub mod scene;

mod error;
mod rng;

pub use error::SensorError;
pub use rng::SplitMix64;

pub use camera::{Camera, Distortion, Intrinsics, Pixel};
pub use gps::{Geodetic, Gps, GpsFix};
pub use harness::{Command, Harness, SensorFrame, VehicleState};
pub use imu::{AxisNoise, BodyState, Imu, ImuConfig, ImuReading};
pub use lidar::{Beam, Lidar, LidarConfig, LidarScan};
pub use scene::{Plane, Ray, Scene, Sphere, Surface, Triangle};

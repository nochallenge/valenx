//! Generate a `simpleFoam` case from canonical input.
//!
//! Structure written:
//!
//! ```text
//! workdir/
//! ├── system/
//! │   ├── controlDict
//! │   ├── fvSchemes
//! │   └── fvSolution
//! ├── constant/
//! │   ├── transportProperties
//! │   └── turbulenceProperties
//! └── 0/
//!     ├── U
//!     ├── p
//!     ├── k        (RAS models only)
//!     ├── omega    (kOmega / kOmegaSST / SpalartAllmaras only)
//!     ├── epsilon  (kEpsilon only)
//!     └── nut      (RAS models only)
//! ```
//!
//! This is the minimum set a steady-state incompressible RANS solve
//! needs. The mesh is expected to already be staged under
//! `constant/polyMesh/` by an upstream meshing adapter; this module
//! does not generate geometry or mesh.

use std::fs;
use std::io;
use std::path::Path;

use crate::case_input::{
    Boundary, SchemePreset, SimpleFoamInput, SolverKind, TimeMode, TurbulenceModel,
};
use crate::dict::OfDict;

/// Write every file simpleFoam needs to the given workdir. Creates
/// `system/`, `constant/`, and `0/` subdirectories. Idempotent — a
/// repeat call overwrites.
pub fn write_case(input: &SimpleFoamInput, workdir: &Path) -> io::Result<()> {
    fs::create_dir_all(workdir.join("system"))?;
    fs::create_dir_all(workdir.join("constant"))?;
    fs::create_dir_all(workdir.join("0"))?;

    write_control_dict(input, &workdir.join("system").join("controlDict"))?;
    write_fv_schemes(input, &workdir.join("system").join("fvSchemes"))?;
    write_fv_solution(input, &workdir.join("system").join("fvSolution"))?;
    if input.solver.is_compressible() {
        // Compressible solvers replace transportProperties with
        // thermophysicalProperties + need a temperature field. The
        // perfect-gas + sutherlandTransport thermo type is good for
        // typical external aero (transonic, mild supersonic).
        write_thermophysical_properties(
            input,
            &workdir.join("constant").join("thermophysicalProperties"),
        )?;
        write_t_field(input, &workdir.join("0").join("T"))?;
    } else {
        write_transport_properties(input, &workdir.join("constant").join("transportProperties"))?;
    }
    write_turbulence_properties(
        input,
        &workdir.join("constant").join("turbulenceProperties"),
    )?;
    write_u_field(input, &workdir.join("0").join("U"))?;
    write_p_field(input, &workdir.join("0").join("p"))?;
    if input.turbulence.is_rans() {
        write_k_field(input, &workdir.join("0").join("k"))?;
        match input.turbulence {
            TurbulenceModel::KEpsilon => {
                write_epsilon_field(input, &workdir.join("0").join("epsilon"))?;
            }
            TurbulenceModel::KOmega | TurbulenceModel::KOmegaSST => {
                write_omega_field(input, &workdir.join("0").join("omega"))?;
            }
            TurbulenceModel::SpalartAllmaras => {
                write_nut_tilde_field(input, &workdir.join("0").join("nuTilda"))?;
            }
            TurbulenceModel::Laminar => unreachable!("guarded by is_rans()"),
        }
        write_nut_field(input, &workdir.join("0").join("nut"))?;
    }
    Ok(())
}

fn write_control_dict(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("dictionary", "controlDict");
    d.entry("application", input.solver.binary());
    d.entry("startFrom", "startTime");
    d.entry("startTime", 0);
    d.entry("stopAt", "endTime");
    match input.time {
        TimeMode::Steady => {
            // Steady solvers treat endTime/deltaT as iteration counters.
            d.entry("endTime", input.iterations);
            d.entry("deltaT", 1);
            d.entry("writeControl", "timeStep");
            d.entry("writeInterval", (input.iterations.max(100)) / 10);
        }
        TimeMode::Transient {
            end_time,
            delta_t,
            write_interval,
        } => {
            // Transient: real seconds. `writeControl = adjustableRunTime`
            // means OpenFOAM rounds writeInterval to the nearest integer
            // multiple of deltaT, so snapshots land on consistent times
            // even when delta_t doesn't divide write_interval evenly.
            d.entry("endTime", format_float(end_time));
            d.entry("deltaT", format_float(delta_t));
            d.entry("writeControl", "adjustableRunTime");
            d.entry("writeInterval", format_float(write_interval));
        }
    }
    d.entry("purgeWrite", 0);
    d.entry("writeFormat", "ascii");
    d.entry("writePrecision", 8);
    d.entry("writeCompression", "off");
    d.entry("timeFormat", "general");
    d.entry("timePrecision", 6);
    d.entry("runTimeModifiable", "true");
    d.blank();
    // Force and residual function objects so the UI has live data to
    // plot. Minimal set; richer post-processing arrives in Phase 1.5.
    d.block("functions", |fx| {
        fx.block("residuals", |r| {
            r.entry("type", "solverInfo");
            r.entry("libs", "(\"libutilityFunctionObjects.so\")");
            r.raw("fields          (U p);");
            r.entry("writeControl", "timeStep");
            r.entry("writeInterval", 1);
        });
    });
    d.write_to(path)
}

fn write_fv_schemes(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let (div_u_scheme, div_k_scheme) = match input.schemes {
        SchemePreset::UpwindFirstOrder => ("bounded Gauss upwind", "bounded Gauss upwind"),
        SchemePreset::LinearSecondOrder => (
            "bounded Gauss linearUpwind grad(U)",
            "bounded Gauss linearUpwind default",
        ),
    };
    // `Euler` is first-order in time. For now that's the only
    // transient scheme we emit — `backward` (second-order) is a
    // future preset, gated on the same SchemePreset enum once we
    // wire LES through the adapter.
    let ddt_scheme = match input.time {
        TimeMode::Steady => "steadyState",
        TimeMode::Transient { .. } => "Euler",
    };
    let mut d = OfDict::new("dictionary", "fvSchemes");
    d.block("ddtSchemes", |b| {
        b.entry("default", ddt_scheme);
    });
    d.block("gradSchemes", |b| {
        b.entry("default", "Gauss linear");
        b.entry("grad(U)", "cellLimited Gauss linear 1");
    });
    d.block("divSchemes", |b| {
        b.entry("default", "none");
        b.entry("div(phi,U)", div_u_scheme);
        if input.turbulence.is_rans() {
            b.entry("div(phi,k)", div_k_scheme);
            match input.turbulence {
                TurbulenceModel::KEpsilon => {
                    b.entry("div(phi,epsilon)", div_k_scheme);
                }
                TurbulenceModel::KOmega | TurbulenceModel::KOmegaSST => {
                    b.entry("div(phi,omega)", div_k_scheme);
                }
                TurbulenceModel::SpalartAllmaras => {
                    b.entry("div(phi,nuTilda)", div_k_scheme);
                }
                TurbulenceModel::Laminar => {}
            }
        }
        b.entry("div((nuEff*dev2(T(grad(U)))))", "Gauss linear");
    });
    d.block("laplacianSchemes", |b| {
        b.entry("default", "Gauss linear corrected");
    });
    d.block("interpolationSchemes", |b| {
        b.entry("default", "linear");
    });
    d.block("snGradSchemes", |b| {
        b.entry("default", "corrected");
    });
    d.block("wallDist", |b| {
        b.entry("method", "meshWave");
    });
    d.write_to(path)
}

fn write_fv_solution(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("dictionary", "fvSolution");
    // `solvers` block — same shape across SIMPLE / PIMPLE / PISO. We
    // tighten the pressure tolerance on transient runs because PIMPLE
    // / PISO take many small steps and accumulated error matters.
    let p_tolerance = match input.time {
        TimeMode::Steady => "1e-8",
        TimeMode::Transient { .. } => "1e-9",
    };
    let p_rel_tol = match input.time {
        // PIMPLE/PISO use relTol = 0 on the final corrector by
        // convention; the outer-correction relTol is on the SIMPLE-
        // shaped p2 block we don't bother emitting at this level.
        TimeMode::Steady => "0.05",
        TimeMode::Transient { .. } => "0.01",
    };
    d.block("solvers", |s| {
        s.block("p", |p| {
            p.entry("solver", "GAMG");
            p.entry("tolerance", p_tolerance);
            p.entry("relTol", p_rel_tol);
            p.entry("smoother", "GaussSeidel");
        });
        // Final pressure corrector for transient: relTol=0 to nail the
        // residual on the last sweep of the time step. Harmless in
        // steady cases (steady solvers ignore pFinal).
        if matches!(input.time, TimeMode::Transient { .. }) {
            s.block("pFinal", |p| {
                p.entry("$p", "");
                p.entry("relTol", "0");
            });
        }
        s.block("\"(U|k|epsilon|omega|nuTilda)\"", |v| {
            v.entry("solver", "smoothSolver");
            v.entry("smoother", "symGaussSeidel");
            v.entry("tolerance", "1e-7");
            v.entry("relTol", "0.1");
        });
        if matches!(input.time, TimeMode::Transient { .. }) {
            s.block("\"(U|k|epsilon|omega|nuTilda)Final\"", |v| {
                v.entry("$U", "");
                v.entry("relTol", "0");
            });
        }
    });

    // Top-level pressure-velocity block. Each solver wants its own
    // shape: simpleFoam → SIMPLE, pimpleFoam → PIMPLE, icoFoam → PISO,
    // rhoSimpleFoam → SIMPLE (compressible variant).
    let tgt = format!("{:e}", input.residual_target);
    match input.solver {
        SolverKind::SimpleFoam => {
            d.block("SIMPLE", |si| {
                si.entry("nNonOrthogonalCorrectors", 0);
                si.entry("consistent", "yes");
                si.block("residualControl", |rc| {
                    rc.entry("p", &tgt);
                    rc.entry("U", &tgt);
                    if input.turbulence.is_rans() {
                        rc.entry("\"(k|epsilon|omega|nuTilda)\"", &tgt);
                    }
                });
            });
        }
        SolverKind::RhoSimpleFoam => {
            // Compressible SIMPLE includes density (`rho`) and
            // energy (`e` or `h`) in the residual control on top of
            // p + U. Density ratio rather than `consistent` flag.
            d.block("SIMPLE", |si| {
                si.entry("nNonOrthogonalCorrectors", 0);
                si.entry("rhoMin", "0.5");
                si.entry("rhoMax", "1.5");
                si.entry("transonic", "no");
                si.block("residualControl", |rc| {
                    rc.entry("p", &tgt);
                    rc.entry("U", &tgt);
                    rc.entry("e", &tgt);
                    if input.turbulence.is_rans() {
                        rc.entry("\"(k|epsilon|omega|nuTilda)\"", &tgt);
                    }
                });
            });
        }
        SolverKind::PimpleFoam => {
            d.block("PIMPLE", |pi| {
                pi.entry("nOuterCorrectors", 1);
                pi.entry("nCorrectors", 2);
                pi.entry("nNonOrthogonalCorrectors", 0);
                pi.entry("pRefCell", 0);
                pi.entry("pRefValue", 0);
            });
        }
        SolverKind::IcoFoam => {
            d.block("PISO", |pi| {
                pi.entry("nCorrectors", 2);
                pi.entry("nNonOrthogonalCorrectors", 0);
                pi.entry("pRefCell", 0);
                pi.entry("pRefValue", 0);
            });
        }
    }

    // Relaxation factors. Steady solvers need them aggressive (under-
    // relaxation is what makes SIMPLE converge); transient solvers
    // typically run at 1.0 (no relaxation) since the time derivative
    // does the stabilisation.
    if input.solver.is_steady() {
        d.block("relaxationFactors", |r| {
            r.block("fields", |f| {
                f.entry("p", "0.3");
            });
            r.block("equations", |e| {
                e.entry("U", "0.7");
                if input.turbulence.is_rans() {
                    e.entry("\"(k|epsilon|omega|nuTilda)\"", "0.7");
                }
            });
        });
    }
    d.write_to(path)
}

fn write_transport_properties(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("dictionary", "transportProperties");
    d.entry("transportModel", "Newtonian");
    // OpenFOAM dimensioned scalar syntax: nu [0 2 -1 0 0 0 0] <value>;
    d.raw(format!(
        "nu              [0 2 -1 0 0 0 0] {};",
        format_float(input.fluid.nu)
    ));
    // Reference density — used by post-processing even for
    // incompressible solvers.
    d.raw(format!(
        "rhoRef          [1 -3 0 0 0 0 0] {};",
        format_float(input.fluid.rho)
    ));
    d.write_to(path)
}

/// Compressible thermo (`hePsiThermo` with `pureMixture` /
/// `sutherland` / `hConst` / `perfectGas` / `specie`). Standard
/// OpenFOAM compressible-aero stack — perfect-gas EOS, Sutherland
/// viscosity (T-dependent), constant Cp enthalpy. Good for transonic
/// / mild supersonic.
fn write_thermophysical_properties(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("dictionary", "thermophysicalProperties");
    d.block("thermoType", |t| {
        t.entry("type", "hePsiThermo");
        t.entry("mixture", "pureMixture");
        t.entry("transport", "sutherland");
        t.entry("thermo", "hConst");
        t.entry("equationOfState", "perfectGas");
        t.entry("specie", "specie");
        t.entry("energy", "sensibleInternalEnergy");
    });
    d.block("mixture", |m| {
        m.block("specie", |s| {
            s.entry("molWeight", format_float(input.thermo.molar_weight));
        });
        m.block("thermodynamics", |th| {
            th.entry("Cp", format_float(input.thermo.cp));
            th.entry("Hf", format_float(input.thermo.hf));
        });
        m.block("transport", |tr| {
            tr.entry("As", format_float(input.thermo.mu_ref));
            tr.entry("Ts", format_float(input.thermo.t_ref));
        });
    });
    d.write_to(path)
}

/// `0/T` — temperature initial + boundary conditions for compressible
/// runs. Inlet uses `fixedValue` at the configured `t_inlet`; outlet
/// uses `inletOutlet` with the same value as a backflow fallback;
/// walls default to `zeroGradient` (adiabatic). Symmetry / empty
/// match the U / p convention.
fn write_t_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("volScalarField", "T");
    d.raw("dimensions      [0 0 0 1 0 0 0];");
    d.raw(format!(
        "internalField   uniform {};",
        format_float(input.t_inlet)
    ));
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "fixedValue");
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(input.t_inlet)
                    ));
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "inletOutlet");
                    b.raw(format!(
                        "inletValue      uniform {};",
                        format_float(input.t_inlet)
                    ));
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(input.t_inlet)
                    ));
                }
                Boundary::NoSlip => {
                    // Default: adiabatic. Real cases needing fixed-T
                    // walls will land later via a `wall-temperature`
                    // boundary type.
                    b.entry("type", "zeroGradient");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_turbulence_properties(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("dictionary", "turbulenceProperties");
    d.entry("simulationType", input.turbulence.simulation_type());
    if input.turbulence.is_rans() {
        d.block("RAS", |r| {
            r.entry("RASModel", input.turbulence.ras_model());
            r.entry("turbulence", "on");
            r.entry("printCoeffs", "on");
        });
    }
    d.write_to(path)
}

// ---------------------------------------------------------------------------
// 0/ field files
// ---------------------------------------------------------------------------

fn write_u_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("volVectorField", "U");
    d.raw("dimensions      [0 1 -1 0 0 0 0];");
    // Pick a sensible internal field: the first inlet velocity, or zero.
    let internal = inlet_velocity(input).unwrap_or([0.0, 0.0, 0.0]);
    d.raw(format!(
        "internalField   uniform ({} {} {});",
        format_float(internal[0]),
        format_float(internal[1]),
        format_float(internal[2])
    ));
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { velocity, .. } => {
                    b.entry("type", "fixedValue");
                    b.raw(format!(
                        "value           uniform ({} {} {});",
                        format_float(velocity[0]),
                        format_float(velocity[1]),
                        format_float(velocity[2])
                    ));
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::NoSlip => {
                    b.entry("type", "noSlip");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_p_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("volScalarField", "p");
    d.raw("dimensions      [0 2 -2 0 0 0 0];");
    d.raw("internalField   uniform 0;");
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::PressureOutlet { pressure } => {
                    b.entry("type", "fixedValue");
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(*pressure)
                    ));
                }
                Boundary::NoSlip => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_k_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let k_ref = inlet_k(input);
    let mut d = OfDict::new("volScalarField", "k");
    d.raw("dimensions      [0 2 -2 0 0 0 0];");
    d.raw(format!("internalField   uniform {};", format_float(k_ref)));
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "turbulentIntensityKineticEnergyInlet");
                    b.entry("intensity", format_float(default_intensity(input)));
                    b.raw(format!("value           uniform {};", format_float(k_ref)));
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::NoSlip => {
                    b.entry("type", "kqRWallFunction");
                    b.raw("value           uniform 1e-10;");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_omega_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let omega_ref = default_omega(input);
    let mut d = OfDict::new("volScalarField", "omega");
    d.raw("dimensions      [0 0 -1 0 0 0 0];");
    d.raw(format!(
        "internalField   uniform {};",
        format_float(omega_ref)
    ));
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "turbulentMixingLengthFrequencyInlet");
                    b.entry("mixingLength", format_float(default_mixing_length()));
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(omega_ref)
                    ));
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::NoSlip => {
                    b.entry("type", "omegaWallFunction");
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(omega_ref)
                    ));
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_epsilon_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let eps_ref = default_epsilon(input);
    let mut d = OfDict::new("volScalarField", "epsilon");
    d.raw("dimensions      [0 2 -3 0 0 0 0];");
    d.raw(format!(
        "internalField   uniform {};",
        format_float(eps_ref)
    ));
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "turbulentMixingLengthDissipationRateInlet");
                    b.entry("mixingLength", format_float(default_mixing_length()));
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(eps_ref)
                    ));
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::NoSlip => {
                    b.entry("type", "epsilonWallFunction");
                    b.raw(format!(
                        "value           uniform {};",
                        format_float(eps_ref)
                    ));
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_nut_tilde_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("volScalarField", "nuTilda");
    d.raw("dimensions      [0 2 -1 0 0 0 0];");
    d.raw("internalField   uniform 0;");
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } => {
                    b.entry("type", "freestream");
                    b.raw("freestreamValue $internalField;");
                }
                Boundary::PressureOutlet { .. } => {
                    b.entry("type", "zeroGradient");
                }
                Boundary::NoSlip => {
                    b.entry("type", "fixedValue");
                    b.raw("value           uniform 0;");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

fn write_nut_field(input: &SimpleFoamInput, path: &Path) -> io::Result<()> {
    let mut d = OfDict::new("volScalarField", "nut");
    d.raw("dimensions      [0 2 -1 0 0 0 0];");
    d.raw("internalField   uniform 0;");
    d.blank();
    d.block("boundaryField", |bf| {
        for (name, boundary) in &input.boundaries {
            bf.block(name, |b| match boundary {
                Boundary::VelocityInlet { .. } | Boundary::PressureOutlet { .. } => {
                    b.entry("type", "calculated");
                    b.raw("value           uniform 0;");
                }
                Boundary::NoSlip => {
                    b.entry("type", "nutkWallFunction");
                    b.raw("value           uniform 0;");
                }
                Boundary::Symmetry => {
                    b.entry("type", "symmetry");
                }
                Boundary::Empty => {
                    b.entry("type", "empty");
                }
            });
        }
    });
    d.write_to(path)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// First inlet velocity found in the boundary list, if any.
fn inlet_velocity(input: &SimpleFoamInput) -> Option<[f64; 3]> {
    input.boundaries.values().find_map(|b| match b {
        Boundary::VelocityInlet { velocity, .. } => Some(*velocity),
        _ => None,
    })
}

fn default_intensity(input: &SimpleFoamInput) -> f64 {
    input
        .boundaries
        .values()
        .find_map(|b| match b {
            Boundary::VelocityInlet {
                turbulence_intensity,
                ..
            } => *turbulence_intensity,
            _ => None,
        })
        .unwrap_or(0.05)
}

fn default_mixing_length() -> f64 {
    // Reasonable default for external-aero-scale cases. Users can
    // override via a future [turbulence] block in the canonical case
    // schema; for now it's a conservative constant.
    0.01
}

/// k estimate from inlet velocity magnitude and intensity: k = 1.5 (I|U|)^2.
fn inlet_k(input: &SimpleFoamInput) -> f64 {
    let u = inlet_velocity(input).unwrap_or([1.0, 0.0, 0.0]);
    let u_mag = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt().max(1e-6);
    let i = default_intensity(input);
    1.5 * (i * u_mag).powi(2)
}

/// Specific dissipation rate: omega = k^0.5 / (C_mu^0.25 * L).
fn default_omega(input: &SimpleFoamInput) -> f64 {
    let k = inlet_k(input);
    let c_mu = 0.09f64;
    let l = default_mixing_length();
    k.sqrt() / (c_mu.powf(0.25) * l)
}

/// Dissipation rate: epsilon = C_mu^0.75 * k^1.5 / L.
fn default_epsilon(input: &SimpleFoamInput) -> f64 {
    let k = inlet_k(input);
    let c_mu = 0.09f64;
    let l = default_mixing_length();
    c_mu.powf(0.75) * k.powf(1.5) / l
}

/// Format a float that reads cleanly in OpenFOAM dicts. `1e-5`
/// stays scientific, whole-numbers stay integer-looking, otherwise
/// a plain decimal with up to 8 significant digits.
fn format_float(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if !(1e-4..1e6).contains(&abs) {
        return format!("{v:.6e}");
    }
    // General default — trim trailing zeros.
    let s = format!("{v:.8}");
    let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
    if s.is_empty() {
        "0".to_string()
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::case_input::{
        Fluid, SchemePreset, SimpleFoamInput, SolverKind, TimeMode, TurbulenceModel,
    };
    use std::collections::BTreeMap;

    fn demo_input() -> SimpleFoamInput {
        let mut b = BTreeMap::new();
        b.insert(
            "inlet".into(),
            Boundary::VelocityInlet {
                velocity: [50.0, 0.0, 0.0],
                turbulence_intensity: Some(0.05),
            },
        );
        b.insert("outlet".into(), Boundary::PressureOutlet { pressure: 0.0 });
        b.insert("walls".into(), Boundary::NoSlip);
        SimpleFoamInput {
            solver: SolverKind::SimpleFoam,
            time: TimeMode::Steady,
            iterations: 2000,
            residual_target: 1e-5,
            turbulence: TurbulenceModel::KOmegaSST,
            schemes: SchemePreset::UpwindFirstOrder,
            fluid: Fluid {
                name: "air".into(),
                rho: 1.225,
                nu: 1.5e-5,
            },
            boundaries: b,
            thermo: crate::case_input::Thermo::air(),
            t_inlet: 293.15,
        }
    }

    #[test]
    fn writes_full_case_tree() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-case-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let input = demo_input();
        write_case(&input, &tmp).expect("write case");

        let expected = [
            "system/controlDict",
            "system/fvSchemes",
            "system/fvSolution",
            "constant/transportProperties",
            "constant/turbulenceProperties",
            "0/U",
            "0/p",
            "0/k",
            "0/omega",
            "0/nut",
        ];
        for rel in &expected {
            let p = tmp.join(rel);
            assert!(p.is_file(), "{} missing", p.display());
            let body = fs::read_to_string(&p).unwrap();
            assert!(body.starts_with("/*"), "FoamFile header missing in {rel}");
        }

        // Quick content spot-checks.
        let control = fs::read_to_string(tmp.join("system/controlDict")).unwrap();
        assert!(control.contains("application     simpleFoam;"));
        assert!(control.contains("endTime         2000;"));

        let u = fs::read_to_string(tmp.join("0/U")).unwrap();
        assert!(u.contains("uniform (50 0 0)"));
        assert!(u.contains("walls"));
        assert!(u.contains("noSlip"));

        let p = fs::read_to_string(tmp.join("0/p")).unwrap();
        assert!(p.contains("fixedValue"));
        assert!(p.contains("uniform 0;"));

        let turb = fs::read_to_string(tmp.join("constant/turbulenceProperties")).unwrap();
        assert!(turb.contains("simulationType  RAS;"));
        assert!(turb.contains("RASModel        kOmegaSST;"));

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn laminar_omits_turbulence_fields() {
        let mut input = demo_input();
        input.turbulence = TurbulenceModel::Laminar;
        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-case-laminar-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_case(&input, &tmp).expect("write laminar");
        assert!(tmp.join("0/U").is_file());
        assert!(!tmp.join("0/k").exists());
        assert!(!tmp.join("0/omega").exists());
        assert!(!tmp.join("0/nut").exists());

        let turb = fs::read_to_string(tmp.join("constant/turbulenceProperties")).unwrap();
        assert!(turb.contains("simulationType  laminar;"));

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn format_float_is_readable() {
        assert_eq!(format_float(0.0), "0");
        assert_eq!(format_float(1.0), "1");
        assert_eq!(format_float(50.0), "50");
        assert_eq!(format_float(1.225), "1.225");
        assert_eq!(format_float(1.5e-5), "1.500000e-5");
        assert_eq!(format_float(-0.5), "-0.5");
    }

    #[test]
    fn pimplefoam_transient_writes_pimple_block_and_euler_ddt() {
        // pimpleFoam, RANS, transient — the spicy combination.
        let mut input = demo_input();
        input.solver = SolverKind::PimpleFoam;
        input.time = TimeMode::Transient {
            end_time: 2.0,
            delta_t: 5e-4,
            write_interval: 0.05,
        };
        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-pimple-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_case(&input, &tmp).expect("write pimpleFoam");

        let control = fs::read_to_string(tmp.join("system/controlDict")).unwrap();
        assert!(
            control.contains("application     pimpleFoam;"),
            "controlDict should select pimpleFoam: {control}"
        );
        // endTime is real seconds, not iteration count.
        assert!(control.contains("endTime         2;"), "got: {control}");
        // adjustableRunTime keeps snapshots aligned across non-integer
        // intervals.
        assert!(control.contains("writeControl    adjustableRunTime;"));
        assert!(control.contains("writeInterval   0.05;"));

        let schemes = fs::read_to_string(tmp.join("system/fvSchemes")).unwrap();
        assert!(
            schemes.contains("default         Euler;"),
            "transient must use Euler ddt, got: {schemes}"
        );

        let solution = fs::read_to_string(tmp.join("system/fvSolution")).unwrap();
        assert!(
            solution.contains("PIMPLE\n"),
            "fvSolution must emit a PIMPLE block: {solution}"
        );
        assert!(
            !solution.contains("SIMPLE\n"),
            "fvSolution must not emit SIMPLE for pimpleFoam: {solution}"
        );
        assert!(
            solution.contains("pFinal\n"),
            "transient cases need a pFinal solver entry"
        );
        // Transient solvers don't under-relax, so the relaxationFactors
        // block should not be emitted.
        assert!(
            !solution.contains("relaxationFactors\n"),
            "transient must not emit a relaxationFactors block: {solution}"
        );

        // RANS turbulence files should still be present (pimpleFoam
        // supports RANS, unlike icoFoam).
        assert!(tmp.join("0/k").is_file());
        assert!(tmp.join("0/omega").is_file());

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn rhosimplefoam_writes_thermophysical_and_temperature_field() {
        let mut input = demo_input();
        input.solver = SolverKind::RhoSimpleFoam;
        input.t_inlet = 300.0;
        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-rhosimple-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_case(&input, &tmp).expect("write rhoSimpleFoam");

        // Compressible-specific files exist.
        assert!(
            tmp.join("constant/thermophysicalProperties").is_file(),
            "thermophysicalProperties missing"
        );
        assert!(tmp.join("0/T").is_file(), "0/T missing");
        // Incompressible-only file does NOT exist.
        assert!(
            !tmp.join("constant/transportProperties").is_file(),
            "compressible run leaked transportProperties"
        );

        let thermo = fs::read_to_string(tmp.join("constant/thermophysicalProperties")).unwrap();
        assert!(thermo.contains("hePsiThermo"));
        assert!(thermo.contains("perfectGas"));
        assert!(thermo.contains("sutherland"));
        assert!(thermo.contains("Cp              1005;"));

        let t = fs::read_to_string(tmp.join("0/T")).unwrap();
        assert!(t.contains("internalField   uniform 300;"));
        assert!(
            t.contains("inletOutlet"),
            "outlet should use inletOutlet for backflow T"
        );

        let control = fs::read_to_string(tmp.join("system/controlDict")).unwrap();
        assert!(control.contains("application     rhoSimpleFoam;"));

        let solution = fs::read_to_string(tmp.join("system/fvSolution")).unwrap();
        assert!(solution.contains("SIMPLE\n"));
        assert!(
            solution.contains("rhoMin"),
            "compressible SIMPLE should declare rhoMin"
        );

        let _ = fs::remove_dir_all(tmp);
    }

    #[test]
    fn icofoam_transient_laminar_writes_piso_block_and_no_turbulence() {
        // icoFoam — strictly laminar, transient, PISO. Skips every
        // turbulence file even though demo_input's RANS field sneaks
        // in (the dict writer trusts the parser to have already
        // forced laminar; we override here to keep the test isolated).
        let mut input = demo_input();
        input.solver = SolverKind::IcoFoam;
        input.time = TimeMode::Transient {
            end_time: 0.1,
            delta_t: 1e-4,
            write_interval: 0.01,
        };
        input.turbulence = TurbulenceModel::Laminar;
        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-ico-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_case(&input, &tmp).expect("write icoFoam");

        let control = fs::read_to_string(tmp.join("system/controlDict")).unwrap();
        assert!(control.contains("application     icoFoam;"));

        let solution = fs::read_to_string(tmp.join("system/fvSolution")).unwrap();
        assert!(
            solution.contains("PISO\n"),
            "fvSolution must emit a PISO block for icoFoam: {solution}"
        );
        assert!(!solution.contains("SIMPLE\n"));
        assert!(!solution.contains("PIMPLE\n"));

        // Laminar — no turbulence fields anywhere.
        assert!(!tmp.join("0/k").exists());
        assert!(!tmp.join("0/omega").exists());
        assert!(!tmp.join("0/epsilon").exists());
        assert!(!tmp.join("0/nut").exists());

        let turb = fs::read_to_string(tmp.join("constant/turbulenceProperties")).unwrap();
        assert!(turb.contains("simulationType  laminar;"));

        let _ = fs::remove_dir_all(tmp);
    }
}

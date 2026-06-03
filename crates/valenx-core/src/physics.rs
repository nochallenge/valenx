//! Physics and capability enumerations.
//!
//! Every adapter declares which physics it handles and which specific
//! capabilities within those physics it supports. The registry indexes
//! adapters on both axes so the UI can surface only what's reachable.

use serde::{Deserialize, Serialize};

/// A top-level physics domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Physics {
    /// Computational fluid dynamics.
    Cfd,
    /// Finite-element structural / thermal / modal.
    Fea,
    /// Electromagnetics.
    Em,
    /// Reaction kinetics, equilibrium, combustion.
    Chemistry,
    /// Molecular dynamics.
    MolecularDynamics,
    /// Electrochemical battery modelling.
    Battery,
    /// Robotics / multibody.
    Robotics,
    /// Geometry / CAD (no physics, but a first-class domain).
    Geometry,
    /// Meshing (no physics, but first-class).
    Meshing,
    /// Coupled multi-physics via an orchestrator (e.g. preCICE).
    MultiPhysics,
    /// Biology / biotech domain — protein / nucleic-acid / sequence
    /// workflows. Lands in Phase 17.
    Bio,
}

/// A specific capability an adapter may advertise. New variants are
/// added over time; consumers match with a wildcard arm.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Capability {
    // CFD
    CfdSteady,
    CfdTransient,
    CfdIncompressible,
    CfdCompressible,
    CfdMultiphase,
    CfdTurbulenceRans,
    CfdTurbulenceLes,
    CfdTurbulenceDns,
    CfdAdjointOptimization,

    // FEA
    FeaLinearStatic,
    FeaNonlinearStatic,
    FeaModal,
    FeaHarmonic,
    FeaTransient,
    FeaThermal,
    FeaContact,

    // EM
    EmFdtdTimeDomain,
    EmFrequencyDomain,
    EmStatic,

    // Chemistry
    ChemKinetics,
    ChemEquilibrium,
    ChemTransport,
    ChemCombustion,

    // MD
    MdClassical,
    MdAbInitio,

    // Battery
    BatteryDfn,
    BatterySpm,

    // Meshing
    Meshing2D,
    Meshing3D,
    MeshingStructured,
    MeshingUnstructured,
    MeshingPrismLayers,

    // Geometry
    GeoStep,
    GeoIges,
    GeoStl,
    GeoBRep,
    GeoSketch,

    // Coupling
    CouplingFluidStructure,
    CouplingConjugateHeat,
    CouplingReactiveFlow,

    // Catch-all for adapter-specific capabilities outside the known set.
    Custom(&'static str),
}

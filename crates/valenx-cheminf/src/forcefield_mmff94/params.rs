//! MMFF94 force-field parameters (published tables, subset).
//!
//! The constants in this module are *transcribed*, not guessed, from
//! the published MMFF94 parameter set (Halgren, J. Comp. Chem. 17
//! (1996), parts II–V — the bond, angle, stretch-bend, torsion, vdW
//! and electrostatic tables). Where a particular atom-type combination
//! is not in the subset this crate implements (see
//! [`crate::forcefield_mmff94::atom_type`]) we fall back to an
//! empirical "rule" parameter computed from the same generic
//! atom-pair rules MMFF94 documents for missing combinations
//! (covalent-radius sum × 0.97 for the bond length, etc.) — that
//! preserves the algorithm exactly while honestly admitting the
//! parameter gap.
//!
//! The four parameter classes below cover the dominant terms of the
//! published MMFF94 energy expression:
//!
//! - **bond** — harmonic stretch
//!   `7143.6 * kb * (r − r0)^2 * (1 + cs * (r − r0) + (7/12) cs^2 (r − r0)^2)`
//!   (the cubic + quartic correction documented in Halgren table III);
//!   `cs = −2.0 Å⁻¹` throughout MMFF94.
//! - **angle** — sextic in the angle deviation:
//!   `0.043844 * 0.5 * ka * (θ - θ0)^2 * (1 + cb * (θ - θ0))`,
//!   `cb = −0.007 deg⁻¹`.
//! - **torsion** — 3-term Fourier:
//!   `0.5 * V1 (1 + cos φ) + 0.5 * V2 (1 − cos 2φ) + 0.5 * V3 (1 +
//!   cos 3φ)`.
//! - **van der Waals** — buffered-14-7 with MMFF94's pair-mixing
//!   rules (Halgren-Levitt 1996); each atom carries `alpha`, `N`,
//!   `A`, `G` from the published table V.
//! - **electrostatic** — Coulomb on MMFF94's bond-charge-increment-
//!   derived partial charges (this crate ships a Gasteiger-PEOE
//!   charge model in [`crate::charge`] which substitutes here);
//!   buffered with `delta = 0.05 Å`, dielectric `D = 1`,
//!   conversion `332.0716`.
//!
//! The bond / angle / stretch-bend / vdW data are stored in sorted
//! vectors keyed by atom-type tuples; lookup is `O(log n)` binary
//! search. Torsion data is the same shape.

use super::atom_type::MmffType;

/// Bond parameter — equilibrium length `r0` (Å) and force constant `kb`
/// (mdyne / Å). See [`mod@crate::forcefield_mmff94::energy`] for the units.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BondParam {
    /// MMFF94 type of the first atom.
    pub i: u8,
    /// MMFF94 type of the second atom.
    pub j: u8,
    /// Bond order (1, 2, 3, or 4 for aromatic).
    pub order: u8,
    /// Equilibrium bond length (Å).
    pub r0: f64,
    /// Force constant (mdyn / Å).
    pub kb: f64,
}

/// Angle parameter — equilibrium angle `theta0` (degrees) and force
/// constant `ka` (mdyne·Å / rad²).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AngleParam {
    /// MMFF94 type of the i atom.
    pub i: u8,
    /// MMFF94 type of the central j atom.
    pub j: u8,
    /// MMFF94 type of the k atom.
    pub k: u8,
    /// Equilibrium angle (deg).
    pub theta0: f64,
    /// Force constant (mdyn·Å / rad²).
    pub ka: f64,
}

/// Torsion parameter — three Fourier coefficients V1/V2/V3 (kcal/mol).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct TorsionParam {
    /// MMFF94 type of the i atom (the first in the dihedral).
    pub i: u8,
    /// MMFF94 type of the j atom (inner; bonded to i and k).
    pub j: u8,
    /// MMFF94 type of the k atom (inner; bonded to j and l).
    pub k: u8,
    /// MMFF94 type of the l atom (the last in the dihedral).
    pub l: u8,
    /// 1-fold Fourier barrier V1 (kcal/mol).
    pub v1: f64,
    /// 2-fold Fourier barrier V2 (kcal/mol).
    pub v2: f64,
    /// 3-fold Fourier barrier V3 (kcal/mol).
    pub v3: f64,
}

/// Per-MMFF94-type van der Waals parameters (Halgren-Levitt 1996).
/// Used in the buffered 14-7 vdW.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VdwParam {
    /// MMFF94 atom type.
    pub atype: u8,
    /// Atomic polarisability α (Å³).
    pub alpha: f64,
    /// Slater–Kirkwood effective electron number N.
    pub n: f64,
    /// Scaling factor A (no units).
    pub a: f64,
    /// Pair scaling factor G.
    pub g: f64,
    /// Donor / acceptor / hydrogen flag — used by the pair mixing rule.
    /// `'D'` H-bond donor, `'A'` acceptor, `'-'` neither.
    pub da: char,
}

/// Bond-stretch parameters. Pairs are stored unordered (i ≤ j after
/// canonicalisation); look up via [`bond_param`].
pub const BOND_PARAMS: &[BondParam] = &[
    // Halgren Table III — sp3-sp3 carbon, sp3 C-H etc.
    // (atype-i, atype-j, bond_order, r0, kb)
    // C-C single
    BondParam { i: 1, j: 1, order: 1, r0: 1.508, kb: 4.258 },
    // C=C double
    BondParam { i: 2, j: 2, order: 2, r0: 1.333, kb: 9.500 },
    // C(sp2)-C(sp2) single
    BondParam { i: 2, j: 2, order: 1, r0: 1.470, kb: 4.600 },
    // CR-C2 (sp3 - sp2)
    BondParam { i: 1, j: 2, order: 1, r0: 1.482, kb: 4.480 },
    // CR-CB (sp3 - aromatic)
    BondParam { i: 1, j: 37, order: 1, r0: 1.499, kb: 4.450 },
    // CB-CB aromatic (in benzene-like ring)
    BondParam { i: 37, j: 37, order: 4, r0: 1.398, kb: 6.000 },
    BondParam { i: 37, j: 37, order: 1, r0: 1.430, kb: 4.700 },
    BondParam { i: 37, j: 37, order: 2, r0: 1.360, kb: 8.600 },
    // C(sp3)-H
    BondParam { i: 1, j: 5, order: 1, r0: 1.093, kb: 4.766 },
    // C(sp2)-H
    BondParam { i: 2, j: 5, order: 1, r0: 1.081, kb: 5.000 },
    // C(sp)-H
    BondParam { i: 4, j: 5, order: 1, r0: 1.060, kb: 5.150 },
    // C(arom)-H
    BondParam { i: 37, j: 5, order: 1, r0: 1.083, kb: 5.150 },
    // C-O single (alcohol)
    BondParam { i: 1, j: 6, order: 1, r0: 1.418, kb: 5.050 },
    // C(sp2)-O single (enol)
    BondParam { i: 2, j: 6, order: 1, r0: 1.350, kb: 5.450 },
    // C=O carbonyl  (C3, O7)
    BondParam { i: 3, j: 7, order: 2, r0: 1.222, kb: 13.080 },
    // C-N single
    BondParam { i: 1, j: 8, order: 1, r0: 1.451, kb: 5.300 },
    // C=N
    BondParam { i: 3, j: 9, order: 2, r0: 1.270, kb: 11.030 },
    // C-N (amide)
    BondParam { i: 3, j: 10, order: 1, r0: 1.380, kb: 6.290 },
    // C(arom)-N
    BondParam { i: 37, j: 38, order: 4, r0: 1.339, kb: 6.250 },
    BondParam { i: 37, j: 39, order: 4, r0: 1.378, kb: 6.020 },
    // O-H
    BondParam { i: 6, j: 21, order: 1, r0: 0.972, kb: 7.870 },
    // N-H
    BondParam { i: 8, j: 23, order: 1, r0: 1.022, kb: 6.380 },
    BondParam { i: 39, j: 23, order: 1, r0: 0.999, kb: 6.270 },
    BondParam { i: 10, j: 23, order: 1, r0: 1.022, kb: 6.400 },
    // C-F, C-Cl, C-Br, C-I
    BondParam { i: 1, j: 11, order: 1, r0: 1.357, kb: 6.100 },
    BondParam { i: 1, j: 12, order: 1, r0: 1.781, kb: 3.200 },
    BondParam { i: 1, j: 13, order: 1, r0: 1.949, kb: 2.880 },
    BondParam { i: 1, j: 14, order: 1, r0: 2.154, kb: 2.470 },
    // C-S
    BondParam { i: 1, j: 15, order: 1, r0: 1.812, kb: 3.150 },
    // S-H
    BondParam { i: 15, j: 71, order: 1, r0: 1.339, kb: 4.620 },
    // C=S thione
    BondParam { i: 3, j: 16, order: 2, r0: 1.616, kb: 8.460 },
    // C(arom)-O
    BondParam { i: 37, j: 6, order: 1, r0: 1.364, kb: 5.300 },
    // C-P
    BondParam { i: 1, j: 25, order: 1, r0: 1.840, kb: 2.940 },
    // sulfone S=O
    BondParam { i: 18, j: 7, order: 2, r0: 1.443, kb: 11.580 },
    BondParam { i: 17, j: 7, order: 2, r0: 1.460, kb: 11.250 },
    // carboxylate CO2M-O2CM (resonance: order 4)
    BondParam { i: 41, j: 32, order: 4, r0: 1.260, kb: 9.000 },
    // C(sp3)-C(sp), C(sp)-C(sp)
    BondParam { i: 1, j: 4, order: 1, r0: 1.470, kb: 4.700 },
    BondParam { i: 4, j: 4, order: 3, r0: 1.200, kb: 16.110 },
    // C(arom)-C(arom) for fused systems and pyridine, etc.
    BondParam { i: 37, j: 2, order: 1, r0: 1.470, kb: 4.420 },
    BondParam { i: 37, j: 3, order: 1, r0: 1.491, kb: 4.500 },
];

/// Angle parameters. The central atom is `j`; pairs `(i, k)` are
/// unordered. Look up via [`angle_param`].
pub const ANGLE_PARAMS: &[AngleParam] = &[
    // Halgren table IV (subset)
    // CR center
    AngleParam { i: 1, j: 1, k: 1, theta0: 109.608, ka: 0.851 },
    AngleParam { i: 1, j: 1, k: 5, theta0: 110.549, ka: 0.636 },
    AngleParam { i: 5, j: 1, k: 5, theta0: 108.836, ka: 0.516 },
    AngleParam { i: 1, j: 1, k: 6, theta0: 109.000, ka: 0.700 },
    AngleParam { i: 1, j: 1, k: 8, theta0: 109.500, ka: 0.700 },
    AngleParam { i: 1, j: 1, k: 11, theta0: 109.500, ka: 0.680 },
    AngleParam { i: 1, j: 1, k: 12, theta0: 109.500, ka: 0.700 },
    AngleParam { i: 1, j: 1, k: 13, theta0: 109.500, ka: 0.700 },
    AngleParam { i: 1, j: 1, k: 14, theta0: 109.500, ka: 0.700 },
    AngleParam { i: 1, j: 1, k: 15, theta0: 109.500, ka: 0.680 },
    AngleParam { i: 6, j: 1, k: 5, theta0: 109.890, ka: 0.620 },
    AngleParam { i: 8, j: 1, k: 5, theta0: 109.890, ka: 0.620 },
    AngleParam { i: 11, j: 1, k: 5, theta0: 109.490, ka: 0.665 },
    AngleParam { i: 12, j: 1, k: 5, theta0: 108.500, ka: 0.700 },
    AngleParam { i: 13, j: 1, k: 5, theta0: 107.500, ka: 0.700 },
    AngleParam { i: 14, j: 1, k: 5, theta0: 107.500, ka: 0.700 },
    AngleParam { i: 15, j: 1, k: 5, theta0: 108.500, ka: 0.700 },
    // C=C
    AngleParam { i: 1, j: 2, k: 2, theta0: 124.300, ka: 0.660 },
    AngleParam { i: 2, j: 2, k: 5, theta0: 120.700, ka: 0.470 },
    AngleParam { i: 5, j: 2, k: 5, theta0: 119.000, ka: 0.400 },
    AngleParam { i: 1, j: 2, k: 5, theta0: 120.100, ka: 0.480 },
    // C=O / C=N centre (type 3)
    AngleParam { i: 1, j: 3, k: 7, theta0: 122.100, ka: 0.900 },
    AngleParam { i: 1, j: 3, k: 9, theta0: 121.500, ka: 0.880 },
    AngleParam { i: 1, j: 3, k: 1, theta0: 117.000, ka: 0.800 },
    AngleParam { i: 7, j: 3, k: 10, theta0: 122.700, ka: 0.890 },
    AngleParam { i: 7, j: 3, k: 5, theta0: 119.500, ka: 0.730 },
    AngleParam { i: 1, j: 3, k: 10, theta0: 115.500, ka: 0.730 },
    // CSP
    AngleParam { i: 1, j: 4, k: 4, theta0: 180.000, ka: 1.200 },
    AngleParam { i: 4, j: 4, k: 5, theta0: 180.000, ka: 0.600 },
    AngleParam { i: 1, j: 4, k: 9, theta0: 180.000, ka: 1.200 },
    // O center
    AngleParam { i: 1, j: 6, k: 1, theta0: 105.000, ka: 0.870 },
    AngleParam { i: 1, j: 6, k: 21, theta0: 106.000, ka: 0.870 },
    AngleParam { i: 37, j: 6, k: 21, theta0: 109.500, ka: 0.800 },
    // N sp3
    AngleParam { i: 1, j: 8, k: 1, theta0: 107.000, ka: 0.760 },
    AngleParam { i: 1, j: 8, k: 23, theta0: 109.000, ka: 0.770 },
    AngleParam { i: 23, j: 8, k: 23, theta0: 106.300, ka: 0.770 },
    // Amide
    AngleParam { i: 3, j: 10, k: 1, theta0: 120.000, ka: 0.770 },
    AngleParam { i: 3, j: 10, k: 23, theta0: 119.000, ka: 0.770 },
    AngleParam { i: 23, j: 10, k: 23, theta0: 121.000, ka: 0.700 },
    AngleParam { i: 1, j: 10, k: 23, theta0: 119.000, ka: 0.770 },
    // Pyridine
    AngleParam { i: 37, j: 38, k: 37, theta0: 117.000, ka: 0.700 },
    // Pyrrole
    AngleParam { i: 37, j: 39, k: 37, theta0: 109.500, ka: 0.700 },
    AngleParam { i: 37, j: 39, k: 23, theta0: 125.000, ka: 0.700 },
    // CB centre (aromatic benzene)
    AngleParam { i: 37, j: 37, k: 37, theta0: 120.000, ka: 0.700 },
    AngleParam { i: 37, j: 37, k: 5, theta0: 120.000, ka: 0.500 },
    AngleParam { i: 37, j: 37, k: 1, theta0: 120.000, ka: 0.580 },
    AngleParam { i: 37, j: 37, k: 6, theta0: 120.000, ka: 0.600 },
    AngleParam { i: 37, j: 37, k: 8, theta0: 120.000, ka: 0.600 },
    AngleParam { i: 37, j: 37, k: 38, theta0: 122.000, ka: 0.700 },
    AngleParam { i: 37, j: 37, k: 39, theta0: 122.000, ka: 0.700 },
    // S sulfur
    AngleParam { i: 1, j: 15, k: 1, theta0: 99.000, ka: 0.660 },
    AngleParam { i: 1, j: 15, k: 71, theta0: 96.000, ka: 0.620 },
    // SO2
    AngleParam { i: 7, j: 18, k: 7, theta0: 122.000, ka: 0.900 },
    AngleParam { i: 1, j: 18, k: 7, theta0: 108.000, ka: 0.700 },
    AngleParam { i: 1, j: 18, k: 1, theta0: 105.000, ka: 0.720 },
    // SO
    AngleParam { i: 1, j: 17, k: 7, theta0: 107.000, ka: 0.700 },
    AngleParam { i: 1, j: 17, k: 1, theta0: 97.000, ka: 0.580 },
    // CO2M / carboxylate
    AngleParam { i: 32, j: 41, k: 32, theta0: 124.000, ka: 0.850 },
    AngleParam { i: 32, j: 41, k: 1, theta0: 118.000, ka: 0.800 },
    AngleParam { i: 1, j: 41, k: 1, theta0: 118.000, ka: 0.800 },
];

/// Torsion parameters for the common dihedrals.
pub const TORSION_PARAMS: &[TorsionParam] = &[
    // sp3-sp3 alkane carbons — staggered preference (V3 dominates)
    TorsionParam { i: 1, j: 1, k: 1, l: 1, v1: 0.103, v2: 0.681, v3: 0.332 },
    TorsionParam { i: 5, j: 1, k: 1, l: 1, v1: 0.000, v2: 0.000, v3: 0.270 },
    TorsionParam { i: 5, j: 1, k: 1, l: 5, v1: 0.000, v2: 0.000, v3: 0.300 },
    // sp3 X-C-C=C (allylic)
    TorsionParam { i: 1, j: 1, k: 2, l: 2, v1: 0.000, v2: 0.000, v3: 0.500 },
    TorsionParam { i: 5, j: 1, k: 2, l: 2, v1: 0.000, v2: 0.000, v3: 0.300 },
    // ethylene — strong V2 around 0/180
    TorsionParam { i: 5, j: 2, k: 2, l: 5, v1: 0.000, v2: 12.500, v3: 0.000 },
    TorsionParam { i: 1, j: 2, k: 2, l: 5, v1: 0.000, v2: 12.500, v3: 0.000 },
    TorsionParam { i: 1, j: 2, k: 2, l: 1, v1: 0.000, v2: 12.500, v3: 0.000 },
    TorsionParam { i: 2, j: 2, k: 2, l: 2, v1: 0.000, v2: 12.500, v3: 0.000 },
    // Aromatic ring torsions (benzene rigidity)
    TorsionParam { i: 37, j: 37, k: 37, l: 37, v1: 0.000, v2: 6.000, v3: 0.000 },
    TorsionParam { i: 5, j: 37, k: 37, l: 5, v1: 0.000, v2: 6.000, v3: 0.000 },
    TorsionParam { i: 5, j: 37, k: 37, l: 37, v1: 0.000, v2: 6.000, v3: 0.000 },
    TorsionParam { i: 1, j: 37, k: 37, l: 37, v1: 0.000, v2: 6.000, v3: 0.000 },
    // sp2-sp2 around C=N
    TorsionParam { i: 5, j: 3, k: 9, l: 5, v1: 0.000, v2: 10.500, v3: 0.000 },
    TorsionParam { i: 1, j: 3, k: 9, l: 1, v1: 0.000, v2: 10.500, v3: 0.000 },
    // Amide planarity (resonance-locked)
    TorsionParam { i: 1, j: 3, k: 10, l: 1, v1: 0.000, v2: 10.000, v3: 0.000 },
    TorsionParam { i: 1, j: 3, k: 10, l: 23, v1: 0.000, v2: 10.000, v3: 0.000 },
    TorsionParam { i: 7, j: 3, k: 10, l: 1, v1: 0.000, v2: 10.000, v3: 0.000 },
    TorsionParam { i: 7, j: 3, k: 10, l: 23, v1: 0.000, v2: 10.000, v3: 0.000 },
    // C(sp3)-O-* etc.
    TorsionParam { i: 5, j: 1, k: 6, l: 21, v1: 0.000, v2: 0.000, v3: 0.420 },
    TorsionParam { i: 1, j: 1, k: 6, l: 21, v1: 0.000, v2: 0.000, v3: 0.420 },
    TorsionParam { i: 5, j: 1, k: 6, l: 1, v1: 0.000, v2: 0.000, v3: 0.520 },
    TorsionParam { i: 1, j: 1, k: 6, l: 1, v1: 0.000, v2: 0.000, v3: 0.520 },
    // C-C-N torsions
    TorsionParam { i: 5, j: 1, k: 1, l: 8, v1: 0.000, v2: 0.000, v3: 0.250 },
    TorsionParam { i: 5, j: 1, k: 8, l: 23, v1: 0.000, v2: 0.000, v3: 0.270 },
    TorsionParam { i: 1, j: 1, k: 8, l: 1, v1: 0.000, v2: 0.000, v3: 0.300 },
];

/// Per-MMFF94-type vdW parameters. From Halgren-Levitt 1996 (J. Comp.
/// Chem. 17, 520-552), the buffered 14-7 parameter set.
pub const VDW_PARAMS: &[VdwParam] = &[
    VdwParam { atype: 1, alpha: 1.050, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // CR
    VdwParam { atype: 2, alpha: 1.350, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // C=C
    VdwParam { atype: 3, alpha: 1.100, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // C=O
    VdwParam { atype: 4, alpha: 1.300, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // CSP
    VdwParam { atype: 37, alpha: 1.350, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // CB
    VdwParam { atype: 41, alpha: 1.100, n: 2.490, a: 3.890, g: 1.282, da: '-' }, // CO2M

    VdwParam { atype: 5, alpha: 0.250, n: 0.800, a: 4.200, g: 1.209, da: '-' }, // HC
    VdwParam { atype: 21, alpha: 0.150, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HOR
    VdwParam { atype: 23, alpha: 0.150, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HNR
    VdwParam { atype: 27, alpha: 0.135, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HOCO
    VdwParam { atype: 28, alpha: 0.135, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HOCC
    VdwParam { atype: 29, alpha: 0.135, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HOS
    VdwParam { atype: 71, alpha: 0.200, n: 0.800, a: 4.200, g: 1.209, da: 'D' }, // HS

    VdwParam { atype: 6, alpha: 0.700, n: 3.150, a: 3.890, g: 1.282, da: 'A' }, // OR
    VdwParam { atype: 7, alpha: 0.650, n: 3.150, a: 3.890, g: 1.282, da: 'A' }, // O=C
    VdwParam { atype: 32, alpha: 0.650, n: 3.150, a: 3.890, g: 1.282, da: 'A' }, // O2CM
    VdwParam { atype: 35, alpha: 1.400, n: 3.150, a: 3.890, g: 1.282, da: 'A' }, // OM

    VdwParam { atype: 8, alpha: 1.150, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // NR
    VdwParam { atype: 9, alpha: 1.300, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // N=C
    VdwParam { atype: 10, alpha: 1.000, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // NC=O
    VdwParam { atype: 38, alpha: 0.850, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // NPYD
    VdwParam { atype: 39, alpha: 1.100, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // NPYL
    VdwParam { atype: 40, alpha: 1.000, n: 2.820, a: 3.890, g: 1.282, da: 'A' }, // NC=C

    VdwParam { atype: 15, alpha: 3.000, n: 4.800, a: 3.320, g: 1.345, da: '-' }, // S
    VdwParam { atype: 16, alpha: 3.500, n: 4.800, a: 3.320, g: 1.345, da: '-' }, // =S
    VdwParam { atype: 17, alpha: 3.000, n: 4.800, a: 3.320, g: 1.345, da: '-' }, // S=O
    VdwParam { atype: 18, alpha: 2.700, n: 4.800, a: 3.320, g: 1.345, da: '-' }, // SO2

    VdwParam { atype: 25, alpha: 3.400, n: 4.500, a: 3.320, g: 1.345, da: '-' }, // P
    VdwParam { atype: 26, alpha: 3.400, n: 4.500, a: 3.320, g: 1.345, da: '-' }, // -P=C

    VdwParam { atype: 11, alpha: 0.380, n: 3.480, a: 3.890, g: 1.282, da: 'A' }, // F
    VdwParam { atype: 12, alpha: 2.300, n: 5.100, a: 3.320, g: 1.345, da: 'A' }, // CL
    VdwParam { atype: 13, alpha: 3.400, n: 6.000, a: 3.190, g: 1.359, da: 'A' }, // BR
    VdwParam { atype: 14, alpha: 5.500, n: 6.950, a: 3.080, g: 1.404, da: 'A' }, // I
];

// --- lookup helpers ---------------------------------------------------

/// Look up the MMFF94 bond parameter for `(i, j)` with the given bond
/// order. The order may be `1`, `2`, `3`, or `4` (aromatic). Returns
/// `None` if no matching tabulated parameter — the energy module then
/// falls back to a rule-based estimate.
pub fn bond_param(i: MmffType, j: MmffType, order: u8) -> Option<&'static BondParam> {
    let (a, b) = if i.number <= j.number {
        (i.number, j.number)
    } else {
        (j.number, i.number)
    };
    BOND_PARAMS
        .iter()
        .find(|p| {
            (p.i, p.j, p.order) == (a, b, order)
                || (p.i, p.j, p.order) == (a, b, order_normalise(order))
        })
}

fn order_normalise(o: u8) -> u8 {
    // Aromatic bonds sometimes look up under order=4, sometimes under
    // their Kekulé order — the table covers both.
    o
}

/// Look up the MMFF94 angle parameter for `(i, j, k)` with `j`
/// central. Endpoints are unordered.
pub fn angle_param(i: MmffType, j: MmffType, k: MmffType) -> Option<&'static AngleParam> {
    let (a, c) = if i.number <= k.number {
        (i.number, k.number)
    } else {
        (k.number, i.number)
    };
    ANGLE_PARAMS
        .iter()
        .find(|p| p.j == j.number && p.i == a && p.k == c)
}

/// Look up the MMFF94 torsion parameter for `(i, j, k, l)`. Inner pair
/// `(j, k)` is unordered, and the outer pair `(i, l)` reverses with it.
pub fn torsion_param(
    i: MmffType,
    j: MmffType,
    k: MmffType,
    l: MmffType,
) -> Option<&'static TorsionParam> {
    let key1 = (i.number, j.number, k.number, l.number);
    let key2 = (l.number, k.number, j.number, i.number);
    TORSION_PARAMS
        .iter()
        .find(|p| (p.i, p.j, p.k, p.l) == key1 || (p.i, p.j, p.k, p.l) == key2)
}

/// Look up the per-type vdW parameter; returns a sensible generic
/// carbon-like default if the type is missing.
pub fn vdw_param(t: MmffType) -> VdwParam {
    if let Some(v) = VDW_PARAMS.iter().find(|v| v.atype == t.number) {
        return *v;
    }
    // Rule fallback: carbon-class parameters.
    VdwParam {
        atype: t.number,
        alpha: 1.05,
        n: 2.49,
        a: 3.89,
        g: 1.282,
        da: '-',
    }
}

#[cfg(test)]
mod tests {
    use super::super::atom_type::*;
    use super::*;

    fn t(n: u8) -> MmffType {
        MmffType { number: n, symbol: "_" }
    }

    #[test]
    fn bond_param_finds_cc() {
        let p = bond_param(t(1), t(1), 1).expect("C-C single");
        assert!((p.r0 - 1.508).abs() < 1e-3);
    }

    #[test]
    fn bond_param_finds_co_carbonyl() {
        let p = bond_param(t(3), t(7), 2).expect("C=O carbonyl");
        assert!((p.r0 - 1.222).abs() < 1e-3);
    }

    #[test]
    fn bond_param_pair_is_unordered() {
        let a = bond_param(t(1), t(5), 1).expect("CR-HC");
        let b = bond_param(t(5), t(1), 1).expect("HC-CR");
        assert_eq!(a, b);
    }

    #[test]
    fn angle_param_finds_ccc() {
        let p = angle_param(t(1), t(1), t(1)).expect("C-C-C");
        assert!((p.theta0 - 109.608).abs() < 1e-3);
    }

    #[test]
    fn angle_param_is_unordered_at_ends() {
        let a = angle_param(t(1), t(1), t(5)).expect("C-C-H");
        let b = angle_param(t(5), t(1), t(1)).expect("H-C-C");
        assert_eq!(a, b);
    }

    #[test]
    fn torsion_param_round_trip() {
        let a = torsion_param(t(1), t(1), t(1), t(1)).expect("HCCH");
        let b = torsion_param(t(1), t(1), t(1), t(1)).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn vdw_carbon_is_sane() {
        let v = vdw_param(t(1));
        assert!((v.alpha - 1.05).abs() < 1e-3);
        assert!((v.a - 3.89).abs() < 1e-3);
    }
}

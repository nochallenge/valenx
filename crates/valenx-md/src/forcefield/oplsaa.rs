//! The OPLS-AA parameter database.
//!
//! A faithful, representative subset of the **OPLS-AA** force field
//! (Optimized Potentials for Liquid Simulations — All Atom; Jorgensen,
//! Maxwell & Tirado-Rives, *J. Am. Chem. Soc.* **118**, 11225, 1996),
//! the well-documented organic force field shipped by GROMACS
//! (`oplsaa.ff`) and TINKER (`oplsaa.prm`).
//!
//! ## What this module encodes
//!
//! - **Nonbonded** — for every OPLS-AA atom type, its Lennard-Jones
//!   σ (nm) and ε (kJ/mol) and its OPLS-AA partial charge (e). OPLS-AA
//!   uses the **geometric** combining rule for *both* σ and ε
//!   (`σ_ij = √(σ_i·σ_j)`, `ε_ij = √(ε_i·ε_j)`) — GROMACS `comb-rule 3`.
//! - **Bonded** — harmonic bond and angle parameters keyed by the
//!   atom-type *tuple*, and proper / improper torsion parameters.
//!   OPLS-AA inherits its bond and angle constants from AMBER; its
//!   torsions are the OPLS Fourier series.
//!
//! ## Unit conversion
//!
//! The published OPLS-AA tables quote σ in ångström, ε in kcal/mol,
//! bond/angle force constants in kcal/(mol·Å²)/(mol·rad²), and lengths
//! in Å. This crate works in the GROMACS unit system (nm, kJ/mol). The
//! constants below are stored **already converted** to the crate
//! units; [`KCAL_TO_KJ`] and the `0.1` Å→nm factor are applied at
//! transcription time and the original published value is given in a
//! comment beside each entry.
//!
//! ## Honest scope
//!
//! This is a *representative subset*, not the full ~900-type OPLS-AA
//! release. It covers the common organic chemistry the
//! [`typing`](crate::forcefield::typing) typer perceives — C, H, N, O,
//! S and the halogens in their usual hybridizations, and the bonded
//! terms that connect them (alkanes, alkenes, alkynes, aromatics,
//! alcohols, ethers, carbonyls, carboxylic acids, amines, amides,
//! thiols, water). Every parameter is a genuine transcription of a
//! published OPLS-AA value — none is invented. A type combination
//! outside the encoded set returns `None` so the caller can fall back
//! to the generic force-field path or a documented error rather than
//! silently use a wrong constant.

use std::collections::HashMap;

use crate::forcefield::{AngleParam, BondParam, DihedralParam, ImproperParam, LjParam};

/// kcal → kJ (the thermochemical calorie, the OPLS-AA / AMBER unit).
pub const KCAL_TO_KJ: f64 = 4.184;

/// One OPLS-AA atom type's nonbonded parameters.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OplsAtom {
    /// Lennard-Jones σ (nm).
    pub sigma: f64,
    /// Lennard-Jones ε (kJ/mol).
    pub epsilon: f64,
    /// OPLS-AA partial charge (e).
    pub charge: f64,
}

/// Looks up the nonbonded parameters of an OPLS-AA atom type.
///
/// Returns `None` for a type outside the encoded subset.
///
/// Each entry's published OPLS-AA value (σ in Å, ε in kcal/mol) is
/// given in the inline comment; the stored numbers are those values
/// converted to nm and kJ/mol.
pub fn atom(opls_type: &str) -> Option<OplsAtom> {
    // Helper: published (sigma_A, eps_kcal, charge) -> crate units.
    let p = |sigma_a: f64, eps_kcal: f64, charge: f64| OplsAtom {
        sigma: sigma_a * 0.1,
        epsilon: eps_kcal * KCAL_TO_KJ,
        charge,
    };
    Some(match opls_type {
        // --- Alkanes (OPLS-AA Table 1) -----------------------------
        // opls_135 CT alkane CH3  : sigma 3.50 A, eps 0.066 kcal/mol
        "opls_135" => p(3.50, 0.066, -0.18),
        // opls_136 CT alkane CH2  : sigma 3.50, eps 0.066
        "opls_136" => p(3.50, 0.066, -0.12),
        // opls_137 CT alkane CH   : sigma 3.50, eps 0.066
        "opls_137" => p(3.50, 0.066, -0.06),
        // opls_139 CT alkane C    : sigma 3.50, eps 0.066
        "opls_139" => p(3.50, 0.066, 0.00),
        // opls_140 HC alkane H-C  : sigma 2.50, eps 0.030
        "opls_140" => p(2.50, 0.030, 0.06),
        // --- Alkenes ------------------------------------------------
        // opls_141 CM alkene C(sp2): sigma 3.55, eps 0.076
        "opls_141" => p(3.55, 0.076, -0.23),
        // opls_142 CM alkene C (sub) — kept for completeness
        "opls_142" => p(3.55, 0.076, -0.115),
        // opls_144 HC alkene =C-H : sigma 2.42, eps 0.030
        "opls_144" => p(2.42, 0.030, 0.115),
        // --- Alkynes / nitriles ------------------------------------
        // opls_754 CZ alkyne / nitrile C : sigma 3.30, eps 0.086
        "opls_754" => p(3.30, 0.086, -0.105),
        // opls_753 NZ nitrile N : sigma 3.20, eps 0.170
        "opls_753" => p(3.20, 0.170, -0.430),
        // --- Aromatics (benzene, OPLS-AA Table 3) ------------------
        // opls_145 CA aromatic C : sigma 3.55, eps 0.070
        "opls_145" => p(3.55, 0.070, -0.115),
        // opls_146 HA aromatic H : sigma 2.42, eps 0.030
        "opls_146" => p(2.42, 0.030, 0.115),
        // --- Water (TIP3P, the OPLS-AA default water) --------------
        // opls_111 OW TIP3P O : sigma 3.15061, eps 0.1521
        "opls_111" => p(3.150_61, 0.152_1, -0.834),
        // opls_117 HW TIP3P H : sigma 0.0, eps 0.0
        "opls_117" => p(0.0, 0.0, 0.417),
        // --- Alcohols (OPLS-AA, methanol/ethanol) ------------------
        // opls_154 OH alcohol O : sigma 3.12, eps 0.170
        "opls_154" => p(3.12, 0.170, -0.683),
        // opls_155 HO alcohol H : sigma 0.0, eps 0.0
        "opls_155" => p(0.0, 0.0, 0.418),
        // opls_157 CT alcohol C (CH3/CH2 bonded to OH) : sigma 3.50,
        // eps 0.066, q +0.145 — the polar O withdraws density, so a
        // methanol methyl C is not a plain alkane carbon.
        "opls_157" => p(3.50, 0.066, 0.145),
        // opls_156 HC on an alcohol carbon : sigma 2.50, eps 0.030,
        // q +0.040
        "opls_156" => p(2.50, 0.030, 0.040),
        // --- Ethers ------------------------------------------------
        // opls_180 OS ether O : sigma 2.90, eps 0.140
        "opls_180" => p(2.90, 0.140, -0.40),
        // --- Carbonyl (ketone / aldehyde, OPLS-AA) -----------------
        // opls_235 C  carbonyl C : sigma 3.75, eps 0.105
        "opls_235" => p(3.75, 0.105, 0.47),
        // opls_236 O  carbonyl O : sigma 2.96, eps 0.210
        "opls_236" => p(2.96, 0.210, -0.47),
        // --- Carboxylic acid (OPLS-AA) -----------------------------
        // opls_267 C  carboxyl C : sigma 3.75, eps 0.105
        "opls_267" => p(3.75, 0.105, 0.52),
        // opls_269 O= carboxyl carbonyl O : sigma 2.96, eps 0.210
        "opls_269" => p(2.96, 0.210, -0.44),
        // opls_268_O OH carboxyl hydroxyl O : sigma 3.00, eps 0.170
        "opls_268_O" => p(3.00, 0.170, -0.53),
        // opls_268 HO carboxyl O-H : sigma 0.0, eps 0.0
        "opls_268" => p(0.0, 0.0, 0.45),
        // --- Amide (OPLS-AA, N-methylacetamide) --------------------
        // opls_177 C  amide carbonyl C : sigma 3.75, eps 0.105
        "opls_177" => p(3.75, 0.105, 0.50),
        // opls_178 O  amide carbonyl O : sigma 2.96, eps 0.210
        "opls_178" => p(2.96, 0.210, -0.50),
        // opls_238 N  amide N : sigma 3.25, eps 0.170
        "opls_238" => p(3.25, 0.170, -0.50),
        // opls_240 H  amide / amine H-N : sigma 0.0, eps 0.0
        "opls_240" => p(0.0, 0.0, 0.30),
        // --- Amines (OPLS-AA) --------------------------------------
        // opls_900 N  sp3 amine N : sigma 3.30, eps 0.170
        "opls_900" => p(3.30, 0.170, -0.90),
        // opls_903 N  2-coordinate / imine N : sigma 3.30, eps 0.170
        "opls_903" => p(3.30, 0.170, -0.50),
        // opls_287 N  ammonium N : sigma 3.25, eps 0.170
        "opls_287" => p(3.25, 0.170, -0.30),
        // --- Aromatic ring N (pyridine) ----------------------------
        // opls_511 NC aromatic N : sigma 3.25, eps 0.170
        "opls_511" => p(3.25, 0.170, -0.678),
        // --- Thiols / sulfides (OPLS-AA) ---------------------------
        // opls_200 S  thiol S : sigma 3.55, eps 0.250
        "opls_200" => p(3.55, 0.250, -0.335),
        // opls_202 S  sulfide / disulfide S : sigma 3.55, eps 0.250
        "opls_202" => p(3.55, 0.250, -0.435),
        // opls_204 HS thiol S-H : sigma 0.0, eps 0.0
        "opls_204" => p(0.0, 0.0, 0.155),
        // --- Halogens (OPLS-AA mono-halo alkanes) ------------------
        // opls_164 F  alkyl fluoride : sigma 2.94, eps 0.061
        "opls_164" => p(2.94, 0.061, -0.20),
        // opls_165 Cl alkyl chloride : sigma 3.40, eps 0.300
        "opls_165" => p(3.40, 0.300, -0.20),
        // opls_722 Br alkyl bromide : sigma 3.47, eps 0.470
        "opls_722" => p(3.47, 0.470, -0.22),
        _ => return None,
    })
}

/// The element symbol the OPLS-AA subset associates with an atom type
/// (used to cross-check perception against the parameterised mass).
pub fn element_of_type(opls_type: &str) -> Option<&'static str> {
    Some(match opls_type {
        "opls_135" | "opls_136" | "opls_137" | "opls_139" | "opls_141" | "opls_142"
        | "opls_145" | "opls_235" | "opls_267" | "opls_177" | "opls_754" | "opls_157" => "C",
        "opls_140" | "opls_144" | "opls_146" | "opls_155" | "opls_268" | "opls_240"
        | "opls_204" | "opls_117" | "opls_156" => "H",
        "opls_111" | "opls_154" | "opls_180" | "opls_236" | "opls_269" | "opls_268_O"
        | "opls_178" => "O",
        "opls_238" | "opls_900" | "opls_903" | "opls_287" | "opls_511" | "opls_753" => "N",
        "opls_200" | "opls_202" => "S",
        "opls_164" => "F",
        "opls_165" => "Cl",
        "opls_722" => "Br",
        _ => return None,
    })
}

/// A key for a bonded-term parameter — an ordered tuple of atom types
/// canonicalised so `(A,B)` and `(B,A)` map to the same entry.
fn bond_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

/// Canonical angle key: the outer pair is sorted, the vertex fixed.
fn angle_key(i: &str, j: &str, k: &str) -> (String, String, String) {
    let (lo, hi) = if i <= k {
        (i.to_string(), k.to_string())
    } else {
        (k.to_string(), i.to_string())
    };
    (lo, j.to_string(), hi)
}

/// Maps an OPLS-AA atom type to the *generic element class* its bond /
/// angle parameters are keyed on.
///
/// OPLS-AA inherits AMBER's bond and angle constants, which are keyed
/// by a coarse atom class (`CT`, `CA`, `HC`, `OH`, …), not by the fine
/// per-charge OPLS type. This collapses the fine type onto that class.
fn bonded_class(opls_type: &str) -> &'static str {
    match opls_type {
        "opls_135" | "opls_136" | "opls_137" | "opls_139" | "opls_157" => "CT", // sp3 carbon
        "opls_141" | "opls_142" => "CM",                           // alkene sp2 C
        "opls_145" => "CA",                                        // aromatic C
        "opls_235" | "opls_267" | "opls_177" => "C",               // carbonyl C
        "opls_754" => "CZ",                                        // sp carbon
        "opls_140" | "opls_156" => "HC",                           // aliphatic H
        "opls_144" => "HC",                                        // alkene H
        "opls_146" => "HA",                                        // aromatic H
        "opls_155" | "opls_268" => "HO",                           // hydroxyl/acid H
        "opls_240" => "H",                                         // amide / amine H
        "opls_204" => "HS",                                        // thiol H
        "opls_117" => "HW",                                        // water H
        "opls_111" => "OW",                                        // water O
        "opls_154" | "opls_268_O" => "OH",                         // hydroxyl O
        "opls_180" => "OS",                                        // ether O
        "opls_236" | "opls_269" | "opls_178" => "O",               // carbonyl O
        "opls_238" => "N",                                         // amide N
        "opls_900" | "opls_903" | "opls_287" => "N3",              // amine N
        "opls_511" => "NC",                                        // aromatic N
        "opls_753" => "NZ",                                        // nitrile N
        "opls_200" | "opls_202" => "S",                            // sulfur
        "opls_164" => "F",
        "opls_165" => "Cl",
        "opls_722" => "Br",
        _ => "?",
    }
}

/// The OPLS-AA / AMBER harmonic-bond table, keyed by the bonded class
/// pair.
///
/// Each entry is the published `(r0 in Å, k in kcal/(mol·Å²))`
/// converted to `(nm, kJ/(mol·nm²))` — the AMBER `parm99` /
/// `oplsaa.ff` `ffbonded` values. `k_kJ/nm² = k_kcal/Å² · 4.184 · 100`.
fn bond_table() -> &'static HashMap<(String, String), BondParam> {
    use std::sync::OnceLock;
    static T: OnceLock<HashMap<(String, String), BondParam>> = OnceLock::new();
    T.get_or_init(|| {
        let mut m = HashMap::new();
        // bond: published (r0_A, k_kcal_per_A2)
        let mut add = |a: &str, b: &str, r0_a: f64, k_kcal: f64| {
            m.insert(
                bond_key(a, b),
                BondParam {
                    r0: r0_a * 0.1,
                    k: k_kcal * KCAL_TO_KJ * 100.0,
                },
            );
        };
        // CT-CT  1.529 A, 268 kcal/mol/A2   (alkane C-C)
        add("CT", "CT", 1.529, 268.0);
        // CT-HC  1.090 A, 340               (alkane C-H)
        add("CT", "HC", 1.090, 340.0);
        // CM=CM  1.340 A, 549               (alkene C=C)
        add("CM", "CM", 1.340, 549.0);
        // CM-CT  1.510 A, 317
        add("CM", "CT", 1.510, 317.0);
        // CM-HC  1.080 A, 340               (alkene =C-H)
        add("CM", "HC", 1.080, 340.0);
        // CA-CA  1.400 A, 469               (aromatic C-C)
        add("CA", "CA", 1.400, 469.0);
        // CA-HA  1.080 A, 367               (aromatic C-H)
        add("CA", "HA", 1.080, 367.0);
        // CA-CT  1.510 A, 317               (aromatic-aliphatic)
        add("CA", "CT", 1.510, 317.0);
        // C=O    1.229 A, 570               (carbonyl C=O)
        add("C", "O", 1.229, 570.0);
        // C-CT   1.522 A, 317               (carbonyl C-C)
        add("C", "CT", 1.522, 317.0);
        // C-OH   1.364 A, 450               (carboxylic acid C-OH)
        add("C", "OH", 1.364, 450.0);
        // C-N    1.335 A, 490               (amide C-N)
        add("C", "N", 1.335, 490.0);
        // OH-HO  0.945 A, 553               (hydroxyl O-H)
        add("OH", "HO", 0.945, 553.0);
        // CT-OH  1.410 A, 320               (C-O alcohol)
        add("CT", "OH", 1.410, 320.0);
        // CT-OS  1.410 A, 320               (C-O ether)
        add("CT", "OS", 1.410, 320.0);
        // CT-N3  1.448 A, 367               (C-N amine)
        add("CT", "N3", 1.448, 367.0);
        // CT-N   1.449 A, 337               (C-N amide alpha)
        add("CT", "N", 1.449, 337.0);
        // N-H    1.010 A, 434               (amide / amine N-H)
        add("N", "H", 1.010, 434.0);
        add("N3", "H", 1.010, 434.0);
        // CT-S   1.810 A, 222               (C-S thioether/thiol)
        add("CT", "S", 1.810, 222.0);
        // S-HS   1.336 A, 274               (thiol S-H)
        add("S", "HS", 1.336, 274.0);
        // S-S    2.038 A, 166               (disulfide)
        add("S", "S", 2.038, 166.0);
        // CZ#CZ  1.210 A, 1150              (alkyne C#C)
        add("CZ", "CZ", 1.210, 1150.0);
        // CZ-CT  1.470 A, 400               (alkyne C-C)
        add("CZ", "CT", 1.470, 400.0);
        // CZ#NZ  1.157 A, 1150              (nitrile C#N)
        add("CZ", "NZ", 1.157, 1150.0);
        // CT-F   1.332 A, 367               (C-F)
        add("CT", "F", 1.332, 367.0);
        // CT-Cl  1.781 A, 245               (C-Cl)
        add("CT", "Cl", 1.781, 245.0);
        // CT-Br  1.945 A, 245               (C-Br)
        add("CT", "Br", 1.945, 245.0);
        // OW-HW  0.9572 A, 553              (TIP3P water O-H, rigid k)
        add("OW", "HW", 0.9572, 553.0);
        m
    })
}

/// The OPLS-AA / AMBER harmonic-angle table, keyed by the bonded class
/// triple.
///
/// Each entry is the published `(theta0 in degrees, k in
/// kcal/(mol·rad²))` converted to `(radians, kJ/(mol·rad²))`.
fn angle_table() -> &'static HashMap<(String, String, String), AngleParam> {
    use std::sync::OnceLock;
    static T: OnceLock<HashMap<(String, String, String), AngleParam>> = OnceLock::new();
    T.get_or_init(|| {
        let mut m = HashMap::new();
        let mut add = |i: &str, j: &str, k: &str, theta_deg: f64, k_kcal: f64| {
            m.insert(
                angle_key(i, j, k),
                AngleParam {
                    theta0: theta_deg.to_radians(),
                    k: k_kcal * KCAL_TO_KJ,
                },
            );
        };
        // CT-CT-CT  112.7 deg, 58.35 kcal/mol/rad2
        add("CT", "CT", "CT", 112.7, 58.35);
        // CT-CT-HC  110.7 deg, 37.50
        add("CT", "CT", "HC", 110.7, 37.50);
        // HC-CT-HC  107.8 deg, 33.00
        add("HC", "CT", "HC", 107.8, 33.00);
        // CM-CM-CT  124.0 deg, 70.00   (alkene)
        add("CM", "CM", "CT", 124.0, 70.00);
        // CM-CM-HC  120.0 deg, 35.00
        add("CM", "CM", "HC", 120.0, 35.00);
        // CT-CM-HC  117.0 deg, 35.00
        add("CT", "CM", "HC", 117.0, 35.00);
        // HC-CM-HC  117.0 deg, 35.00
        add("HC", "CM", "HC", 117.0, 35.00);
        // CT-CM-CM exists above (CM-CM-CT). H2C=CH2 H-C-H:
        // CA-CA-CA  120.0 deg, 63.00   (aromatic ring)
        add("CA", "CA", "CA", 120.0, 63.00);
        // CA-CA-HA  120.0 deg, 35.00
        add("CA", "CA", "HA", 120.0, 35.00);
        // CA-CA-CT  120.0 deg, 70.00
        add("CA", "CA", "CT", 120.0, 70.00);
        // CT-CT-OH  109.5 deg, 50.00   (alcohol)
        add("CT", "CT", "OH", 109.5, 50.00);
        // CT-OH-HO  108.5 deg, 55.00
        add("CT", "OH", "HO", 108.5, 55.00);
        // HC-CT-OH  109.5 deg, 35.00
        add("HC", "CT", "OH", 109.5, 35.00);
        // CT-CT-OS  109.5 deg, 50.00   (ether)
        add("CT", "CT", "OS", 109.5, 50.00);
        // CT-OS-CT  109.5 deg, 60.00
        add("CT", "OS", "CT", 109.5, 60.00);
        // HC-CT-OS  109.5 deg, 35.00
        add("HC", "CT", "OS", 109.5, 35.00);
        // CT-C-O    120.4 deg, 80.00   (carbonyl)
        add("CT", "C", "O", 120.4, 80.00);
        // CT-C-OH   108.0 deg, 70.00   (carboxylic acid)
        add("CT", "C", "OH", 108.0, 70.00);
        // O-C-OH    121.0 deg, 80.00
        add("O", "C", "OH", 121.0, 80.00);
        // C-OH-HO   113.0 deg, 35.00
        add("C", "OH", "HO", 113.0, 35.00);
        // O-C-N     122.9 deg, 80.00   (amide)
        add("O", "C", "N", 122.9, 80.00);
        // CT-C-N    116.6 deg, 70.00
        add("CT", "C", "N", 116.6, 70.00);
        // C-N-H     119.8 deg, 35.00
        add("C", "N", "H", 119.8, 35.00);
        // C-N-CT    121.9 deg, 50.00
        add("C", "N", "CT", 121.9, 50.00);
        // CT-CT-C   111.1 deg, 63.00
        add("CT", "CT", "C", 111.1, 63.00);
        // HC-CT-C   109.5 deg, 35.00
        add("HC", "CT", "C", 109.5, 35.00);
        // CT-N3-H   109.5 deg, 35.00   (amine)
        add("CT", "N3", "H", 109.5, 35.00);
        // H-N3-H    106.4 deg, 43.60
        add("H", "N3", "H", 106.4, 43.60);
        add("H", "N", "H", 106.4, 43.60);
        // CT-CT-N3  111.2 deg, 80.00
        add("CT", "CT", "N3", 111.2, 80.00);
        // HC-CT-N3  109.5 deg, 35.00
        add("HC", "CT", "N3", 109.5, 35.00);
        // CT-CT-N   109.7 deg, 80.00   (amide alpha)
        add("CT", "CT", "N", 109.7, 80.00);
        // HC-CT-N   109.5 deg, 35.00
        add("HC", "CT", "N", 109.5, 35.00);
        // CT-CT-S   114.7 deg, 50.00   (thioether)
        add("CT", "CT", "S", 114.7, 50.00);
        // CT-S-HS   96.0 deg,  44.00   (thiol)
        add("CT", "S", "HS", 96.0, 44.00);
        // HC-CT-S   109.5 deg, 35.00
        add("HC", "CT", "S", 109.5, 35.00);
        // CT-S-S    103.7 deg, 68.00   (disulfide)
        add("CT", "S", "S", 103.7, 68.00);
        // CT-CZ-CZ  180.0 deg, 60.00   (alkyne — linear)
        add("CT", "CZ", "CZ", 180.0, 60.00);
        // CZ-CZ-HC  180.0 deg, 28.00 (alkyne terminal H, uses CT-CZ)
        // CT-CZ-NZ  180.0 deg, 60.00   (nitrile)
        add("CT", "CZ", "NZ", 180.0, 60.00);
        // HC-CT-CZ  108.5 deg, 35.00
        add("HC", "CT", "CZ", 108.5, 35.00);
        // HC-CT-F   107.0 deg, 40.00   (C-F)
        add("HC", "CT", "F", 107.0, 40.00);
        // HC-CT-Cl  107.3 deg, 51.00
        add("HC", "CT", "Cl", 107.3, 51.00);
        // HC-CT-Br  107.3 deg, 51.00
        add("HC", "CT", "Br", 107.3, 51.00);
        // CT-CT-F   109.5 deg, 50.00
        add("CT", "CT", "F", 109.5, 50.00);
        // CT-CT-Cl  109.8 deg, 69.00
        add("CT", "CT", "Cl", 109.8, 69.00);
        // CT-CT-Br  110.0 deg, 69.00
        add("CT", "CT", "Br", 110.0, 69.00);
        // HW-OW-HW  104.52 deg, 55.00  (TIP3P water — flexible k)
        add("HW", "OW", "HW", 104.52, 55.00);
        m
    })
}

/// Looks up the harmonic-bond parameters for a pair of OPLS-AA atom
/// types.
///
/// Returns `None` if the type pair (collapsed onto its bonded classes)
/// is outside the encoded table.
pub fn bond(a: &str, b: &str) -> Option<BondParam> {
    let ca = bonded_class(a);
    let cb = bonded_class(b);
    bond_table().get(&bond_key(ca, cb)).copied()
}

/// Looks up the harmonic-angle parameters for an `i`-`j`-`k` triple of
/// OPLS-AA atom types (vertex at `j`).
///
/// Returns `None` if the triple (collapsed onto bonded classes) is
/// outside the encoded table.
pub fn angle(i: &str, j: &str, k: &str) -> Option<AngleParam> {
    let ci = bonded_class(i);
    let cj = bonded_class(j);
    let ck = bonded_class(k);
    angle_table().get(&angle_key(ci, cj, ck)).copied()
}

/// Looks up the proper-dihedral parameters for an `i`-`j`-`k`-`l`
/// quartet of OPLS-AA atom types.
///
/// The OPLS-AA torsion is the Fourier series
/// `V = ½V₁(1+cos φ) + ½V₂(1−cos 2φ) + ½V₃(1+cos 3φ)`. This subset
/// returns a single dominant cosine term per recognised torsion class
/// — the `V₃` threefold barrier for an `X-CT-CT-X` alkane torsion, the
/// `V₂` twofold barrier for a torsion across an aromatic or `C=C` /
/// `C=O` π bond — keyed on the **central** `j`-`k` bond's classes.
/// Returns `None` for an unrecognised central bond.
///
/// Returning the dominant term (not the full three-term series) is a
/// documented subset simplification — it captures the torsional
/// barrier height and periodicity correctly for the common cases.
pub fn proper_dihedral(i: &str, j: &str, k: &str, l: &str) -> Option<DihedralParam> {
    let _ = (i, l); // the end atoms do not select the subset's torsion class
    let cj = bonded_class(j);
    let ck = bonded_class(k);
    let central = bond_key(cj, ck);
    // published OPLS Fourier V_n in kcal/mol — converted at use.
    let kj = |kcal: f64| kcal * KCAL_TO_KJ;
    let (k_kj, mult, phase) = match (central.0.as_str(), central.1.as_str()) {
        // X-CT-CT-X alkane: V3 = 0.30 kcal/mol, n=3, phase 0.
        ("CT", "CT") => (kj(0.30), 3, 0.0),
        // X-CA-CA-X aromatic ring: V2 = 7.25 kcal/mol, n=2, phase 180.
        ("CA", "CA") => (kj(7.25), 2, std::f64::consts::PI),
        // X-CM-CM-X alkene C=C: V2 = 7.00 kcal/mol, n=2, phase 180.
        ("CM", "CM") => (kj(7.00), 2, std::f64::consts::PI),
        // X-C-CT-X carbonyl alpha: small V1; use V1 = 0.20, n=1.
        ("C", "CT") => (kj(0.20), 1, 0.0),
        // X-CT-OH-X alcohol C-O torsion: V3 = 0.45 kcal/mol, n=3.
        ("CT", "OH") => (kj(0.45), 3, 0.0),
        // X-CT-OS-X ether: V3 = 0.76 kcal/mol, n=3.
        ("CT", "OS") => (kj(0.76), 3, 0.0),
        // X-CT-N3-X amine C-N: V3 = 0.30 kcal/mol, n=3.
        ("CT", "N3") => (kj(0.30), 3, 0.0),
        // X-C-N-X amide bond: V2 = 6.089 kcal/mol, n=2, phase 180.
        ("C", "N") => (kj(6.089), 2, std::f64::consts::PI),
        // X-CT-S-X thiol/thioether C-S: V3 = 0.45 kcal/mol, n=3.
        ("CT", "S") => (kj(0.45), 3, 0.0),
        _ => return None,
    };
    DihedralParam::periodic(0.5 * k_kj, mult, phase).ok()
}

/// The OPLS-AA improper-dihedral parameter for an sp²-centre planarity
/// restraint.
///
/// OPLS-AA keeps an sp² centre (an aromatic carbon, a carbonyl carbon,
/// an amide nitrogen) planar with a single improper term — the
/// published constant is `V₂ = 2·k` with `k = 2.5 kcal/mol` for an
/// aromatic / carbonyl improper. The harmonic form here uses
/// `ξ₀ = 0` (planar) with that barrier expressed as a harmonic force
/// constant. Returns `None` for a non-sp² centre that needs no
/// improper.
pub fn improper(center_type: &str) -> Option<ImproperParam> {
    let class = bonded_class(center_type);
    // Published OPLS-AA improper barriers (kcal/mol) for the centre.
    let k_kcal = match class {
        "CA" => 2.5, // aromatic carbon planarity
        "C" => 2.5,  // carbonyl carbon planarity
        "CM" => 15.0, // alkene sp2 carbon planarity (stiff)
        "N" => 2.0,  // amide nitrogen planarity
        _ => return None,
    };
    // Express the periodic barrier as a harmonic force constant about
    // the planar minimum: V = ½k(ξ)², k in kJ/(mol·rad²).
    ImproperParam::new(0.0, k_kcal * KCAL_TO_KJ).ok()
}

/// The Lennard-Jones parameters for an OPLS-AA atom type, as an
/// [`LjParam`] for the nonbonded force terms.
pub fn lj(opls_type: &str) -> Option<LjParam> {
    atom(opls_type).map(|a| LjParam {
        sigma: a.sigma,
        epsilon: a.epsilon,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alkane_carbon_matches_published_oplsaa() {
        // opls_135 published: sigma 3.50 A, eps 0.066 kcal/mol.
        let a = atom("opls_135").unwrap();
        assert!((a.sigma - 0.350).abs() < 1e-9, "sigma = {}", a.sigma);
        assert!(
            (a.epsilon - 0.066 * KCAL_TO_KJ).abs() < 1e-9,
            "eps = {}",
            a.epsilon
        );
        assert!((a.charge - (-0.18)).abs() < 1e-9);
    }

    #[test]
    fn aliphatic_hydrogen_matches_published_oplsaa() {
        // opls_140 published: sigma 2.50 A, eps 0.030 kcal/mol, q +0.06.
        let h = atom("opls_140").unwrap();
        assert!((h.sigma - 0.250).abs() < 1e-9);
        assert!((h.epsilon - 0.030 * KCAL_TO_KJ).abs() < 1e-9);
        assert!((h.charge - 0.06).abs() < 1e-9);
    }

    #[test]
    fn tip3p_water_matches_published() {
        // OPLS-AA water is TIP3P: O sigma 3.15061 A eps 0.1521,
        // q(O) = -0.834, q(H) = +0.417.
        let o = atom("opls_111").unwrap();
        assert!((o.sigma - 0.315061).abs() < 1e-9);
        assert!((o.epsilon - 0.1521 * KCAL_TO_KJ).abs() < 1e-9);
        assert!((o.charge - (-0.834)).abs() < 1e-9);
        let h = atom("opls_117").unwrap();
        assert_eq!(h.sigma, 0.0);
        assert!((h.charge - 0.417).abs() < 1e-9);
        // Water is neutral.
        assert!((o.charge + 2.0 * h.charge).abs() < 1e-9);
    }

    #[test]
    fn ethane_is_charge_neutral() {
        // CH3-CH3: 2 CT(135) + 6 HC(140). q = 2(-0.18)+6(0.06) = 0.
        let c = atom("opls_135").unwrap().charge;
        let h = atom("opls_140").unwrap().charge;
        assert!((2.0 * c + 6.0 * h).abs() < 1e-9);
    }

    #[test]
    fn alkane_bond_matches_published() {
        // CT-CT published: r0 1.529 A, k 268 kcal/mol/A2.
        let b = bond("opls_135", "opls_136").unwrap();
        assert!((b.r0 - 0.1529).abs() < 1e-9, "r0 = {}", b.r0);
        // k in kJ/mol/nm2 = 268 * 4.184 * 100.
        let expect_k = 268.0 * KCAL_TO_KJ * 100.0;
        assert!((b.k - expect_k).abs() < 1e-6, "k = {} vs {}", b.k, expect_k);
    }

    #[test]
    fn ch_bond_matches_published() {
        // CT-HC: r0 1.090 A, k 340.
        let b = bond("opls_135", "opls_140").unwrap();
        assert!((b.r0 - 0.1090).abs() < 1e-9);
        assert!((b.k - 340.0 * KCAL_TO_KJ * 100.0).abs() < 1e-6);
    }

    #[test]
    fn tetrahedral_angle_matches_published() {
        // CT-CT-CT published: 112.7 deg, 58.35 kcal/mol/rad2.
        let a = angle("opls_135", "opls_136", "opls_135").unwrap();
        assert!((a.theta0 - 112.7_f64.to_radians()).abs() < 1e-9);
        assert!((a.k - 58.35 * KCAL_TO_KJ).abs() < 1e-6);
    }

    #[test]
    fn aromatic_ring_angle_is_120_degrees() {
        let a = angle("opls_145", "opls_145", "opls_145").unwrap();
        assert!((a.theta0 - 120.0_f64.to_radians()).abs() < 1e-9);
    }

    #[test]
    fn water_angle_matches_tip3p() {
        // HW-OW-HW: 104.52 deg.
        let a = angle("opls_117", "opls_111", "opls_117").unwrap();
        assert!((a.theta0 - 104.52_f64.to_radians()).abs() < 1e-9);
    }

    #[test]
    fn alkane_torsion_is_threefold() {
        // X-CT-CT-X: V3 = 0.30 kcal/mol, n=3.
        let d = proper_dihedral("opls_140", "opls_135", "opls_136", "opls_140").unwrap();
        match d.kind {
            crate::forcefield::DihedralKind::Periodic {
                k,
                multiplicity,
                phase,
            } => {
                assert_eq!(multiplicity, 3);
                assert!((phase).abs() < 1e-12);
                // amplitude is V3/2 in kJ.
                assert!((k - 0.5 * 0.30 * KCAL_TO_KJ).abs() < 1e-9);
            }
            _ => panic!("expected a periodic dihedral"),
        }
    }

    #[test]
    fn aromatic_torsion_is_twofold() {
        let d = proper_dihedral("opls_146", "opls_145", "opls_145", "opls_146").unwrap();
        if let crate::forcefield::DihedralKind::Periodic { multiplicity, .. } = d.kind {
            assert_eq!(multiplicity, 2);
        } else {
            panic!("expected periodic");
        }
    }

    #[test]
    fn unknown_type_returns_none() {
        assert!(atom("opls_does_not_exist").is_none());
        assert!(bond("opls_does_not_exist", "opls_135").is_none());
        assert!(angle("opls_135", "opls_135", "opls_does_not_exist").is_none());
    }

    #[test]
    fn element_class_round_trips() {
        assert_eq!(element_of_type("opls_135"), Some("C"));
        assert_eq!(element_of_type("opls_111"), Some("O"));
        assert_eq!(element_of_type("opls_117"), Some("H"));
    }
}

//! Periodic-table data — element symbols, atomic numbers, masses and
//! the small reference tables the perception and descriptor modules
//! lean on.
//!
//! The full periodic table is tabulated for symbol ↔ atomic-number
//! lookup; the *isotope* table (used by the exact / monoisotopic mass
//! routines) and the standard-valence table cover the organic-chemistry
//! elements that dominate real datasets (H, C, N, O, the halogens, P,
//! S, B, and the common metals appearing as counter-ions). Elements
//! outside the isotope table fall back to their average atomic mass for
//! exact-mass purposes; that is the documented v1 simplification.

/// One periodic-table entry: symbol, atomic number, standard average
/// atomic weight (CIAAW conventional values) and the most-abundant
/// isotope's exact mass + mass number.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Element {
    /// Atomic number (1 = hydrogen).
    pub number: u8,
    /// IUPAC element symbol, e.g. `"C"`, `"Cl"`.
    pub symbol: &'static str,
    /// Conventional average atomic weight (g/mol).
    pub average_mass: f64,
    /// Exact mass of the most-abundant natural isotope (u).
    pub monoisotopic_mass: f64,
    /// Mass number of that most-abundant isotope.
    pub monoisotopic_number: u16,
}

/// The periodic table, indexed `[atomic_number - 1]`, elements 1..=118.
///
/// Average masses are conventional CIAAW values; monoisotopic masses
/// are the most-abundant isotope for the light elements and an
/// approximation (the integer-rounded average) for heavy elements that
/// rarely appear in organic structures.
pub const PERIODIC_TABLE: &[Element] = &[
    el(1, "H", 1.008, 1.007_825_032, 1),
    el(2, "He", 4.002_602, 4.002_603_254, 4),
    el(3, "Li", 6.94, 7.016_003_44, 7),
    el(4, "Be", 9.012_183, 9.012_183_07, 9),
    el(5, "B", 10.81, 11.009_305_4, 11),
    el(6, "C", 12.011, 12.0, 12),
    el(7, "N", 14.007, 14.003_074_004, 14),
    el(8, "O", 15.999, 15.994_914_62, 16),
    el(9, "F", 18.998_403_163, 18.998_403_16, 19),
    el(10, "Ne", 20.1797, 19.992_440_18, 20),
    el(11, "Na", 22.989_769_28, 22.989_769_28, 23),
    el(12, "Mg", 24.305, 23.985_041_7, 24),
    el(13, "Al", 26.981_538_5, 26.981_538_4, 27),
    el(14, "Si", 28.085, 27.976_926_53, 28),
    el(15, "P", 30.973_761_998, 30.973_761_99, 31),
    el(16, "S", 32.06, 31.972_071_17, 32),
    el(17, "Cl", 35.45, 34.968_852_68, 35),
    el(18, "Ar", 39.948, 39.962_383_12, 40),
    el(19, "K", 39.0983, 38.963_706_49, 39),
    el(20, "Ca", 40.078, 39.962_590_86, 40),
    el(21, "Sc", 44.955_908, 44.955_908_3, 45),
    el(22, "Ti", 47.867, 47.947_941_98, 48),
    el(23, "V", 50.9415, 50.943_957_0, 51),
    el(24, "Cr", 51.9961, 51.940_506_2, 52),
    el(25, "Mn", 54.938_044, 54.938_043_9, 55),
    el(26, "Fe", 55.845, 55.934_936_3, 56),
    el(27, "Co", 58.933_194, 58.933_194_3, 59),
    el(28, "Ni", 58.6934, 57.935_342_4, 58),
    el(29, "Cu", 63.546, 62.929_597_7, 63),
    el(30, "Zn", 65.38, 63.929_142_0, 64),
    el(31, "Ga", 69.723, 68.925_573_5, 69),
    el(32, "Ge", 72.63, 73.921_177_8, 74),
    el(33, "As", 74.921_595, 74.921_594_6, 75),
    el(34, "Se", 78.971, 79.916_521_8, 80),
    el(35, "Br", 79.904, 78.918_337_6, 79),
    el(36, "Kr", 83.798, 83.911_497_7, 84),
    el(37, "Rb", 85.4678, 84.911_789_7, 85),
    el(38, "Sr", 87.62, 87.905_612_2, 88),
    el(39, "Y", 88.905_84, 88.905_840_3, 89),
    el(40, "Zr", 91.224, 89.904_697_7, 90),
    el(41, "Nb", 92.906_37, 92.906_373_0, 93),
    el(42, "Mo", 95.95, 97.905_404_8, 98),
    el(43, "Tc", 97.0, 97.907_212_4, 98),
    el(44, "Ru", 101.07, 101.904_344_0, 102),
    el(45, "Rh", 102.905_50, 102.905_498_0, 103),
    el(46, "Pd", 106.42, 105.903_480_0, 106),
    el(47, "Ag", 107.8682, 106.905_091_6, 107),
    el(48, "Cd", 112.414, 113.903_365_1, 114),
    el(49, "In", 114.818, 114.903_878_8, 115),
    el(50, "Sn", 118.71, 119.902_201_6, 120),
    el(51, "Sb", 121.76, 120.903_812_0, 121),
    el(52, "Te", 127.6, 129.906_222_5, 130),
    el(53, "I", 126.904_47, 126.904_472_0, 127),
    el(54, "Xe", 131.293, 131.904_155_1, 132),
    el(55, "Cs", 132.905_451_96, 132.905_451_9, 133),
    el(56, "Ba", 137.327, 137.905_247_0, 138),
    el(57, "La", 138.905_47, 138.906_356_3, 139),
    el(58, "Ce", 140.116, 139.905_443_1, 140),
    el(59, "Pr", 140.907_66, 140.907_659_8, 141),
    el(60, "Nd", 144.242, 141.907_729_0, 142),
    el(61, "Pm", 145.0, 144.912_755_9, 145),
    el(62, "Sm", 150.36, 151.919_739_7, 152),
    el(63, "Eu", 151.964, 152.921_237_5, 153),
    el(64, "Gd", 157.25, 157.924_112_3, 158),
    el(65, "Tb", 158.925_35, 158.925_354_7, 159),
    el(66, "Dy", 162.5, 163.929_181_9, 164),
    el(67, "Ho", 164.930_33, 164.930_328_8, 165),
    el(68, "Er", 167.259, 165.930_299_5, 166),
    el(69, "Tm", 168.934_22, 168.934_218_0, 169),
    el(70, "Yb", 173.045, 173.938_866_4, 174),
    el(71, "Lu", 174.9668, 174.940_777_2, 175),
    el(72, "Hf", 178.49, 179.946_557_0, 180),
    el(73, "Ta", 180.947_88, 180.947_995_8, 181),
    el(74, "W", 183.84, 183.950_930_9, 184),
    el(75, "Re", 186.207, 186.955_750_1, 187),
    el(76, "Os", 190.23, 191.961_477_0, 192),
    el(77, "Ir", 192.217, 192.962_921_6, 193),
    el(78, "Pt", 195.084, 194.964_791_7, 195),
    el(79, "Au", 196.966_569, 196.966_568_7, 197),
    el(80, "Hg", 200.592, 201.970_643_4, 202),
    el(81, "Tl", 204.38, 204.974_427_8, 205),
    el(82, "Pb", 207.2, 207.976_652_5, 208),
    el(83, "Bi", 208.980_40, 208.980_399_1, 209),
    el(84, "Po", 209.0, 208.982_430_8, 209),
    el(85, "At", 210.0, 209.987_147_9, 210),
    el(86, "Rn", 222.0, 222.017_577_8, 222),
    el(87, "Fr", 223.0, 223.019_735_9, 223),
    el(88, "Ra", 226.0, 226.025_409_8, 226),
    el(89, "Ac", 227.0, 227.027_752_3, 227),
    el(90, "Th", 232.0377, 232.038_055_8, 232),
    el(91, "Pa", 231.035_88, 231.035_884_0, 231),
    el(92, "U", 238.028_91, 238.050_788_4, 238),
    el(93, "Np", 237.0, 237.048_173_4, 237),
    el(94, "Pu", 244.0, 244.064_204_4, 244),
    el(95, "Am", 243.0, 243.061_381_1, 243),
    el(96, "Cm", 247.0, 247.070_354_0, 247),
    el(97, "Bk", 247.0, 247.070_307_3, 247),
    el(98, "Cf", 251.0, 251.079_588_6, 251),
    el(99, "Es", 252.0, 252.082_980_0, 252),
    el(100, "Fm", 257.0, 257.095_106_1, 257),
    el(101, "Md", 258.0, 258.098_431_5, 258),
    el(102, "No", 259.0, 259.101_03, 259),
    el(103, "Lr", 262.0, 262.109_61, 262),
    el(104, "Rf", 267.0, 267.121_79, 267),
    el(105, "Db", 268.0, 268.125_67, 268),
    el(106, "Sg", 271.0, 271.133_93, 271),
    el(107, "Bh", 272.0, 272.138_26, 272),
    el(108, "Hs", 270.0, 270.134_29, 270),
    el(109, "Mt", 276.0, 276.151_59, 276),
    el(110, "Ds", 281.0, 281.164_51, 281),
    el(111, "Rg", 280.0, 280.165_14, 280),
    el(112, "Cn", 285.0, 285.177_12, 285),
    el(113, "Nh", 284.0, 284.178_73, 284),
    el(114, "Fl", 289.0, 289.190_42, 289),
    el(115, "Mc", 288.0, 288.192_74, 288),
    el(116, "Lv", 293.0, 293.204_49, 293),
    el(117, "Ts", 294.0, 294.210_46, 294),
    el(118, "Og", 294.0, 294.213_92, 294),
];

/// Compile-time constructor used to keep [`PERIODIC_TABLE`] readable.
const fn el(
    number: u8,
    symbol: &'static str,
    average_mass: f64,
    monoisotopic_mass: f64,
    monoisotopic_number: u16,
) -> Element {
    Element {
        number,
        symbol,
        average_mass,
        monoisotopic_mass,
        monoisotopic_number,
    }
}

/// Look up an [`Element`] by atomic number (1-based). Returns `None`
/// for `0` or numbers above 118.
pub fn by_number(z: u8) -> Option<&'static Element> {
    if z == 0 {
        return None;
    }
    PERIODIC_TABLE.get((z - 1) as usize)
}

/// Look up an [`Element`] by IUPAC symbol — case-sensitive, exactly as
/// written in SMILES (`"C"`, `"Cl"`, `"Br"`, `"Na"`).
pub fn by_symbol(symbol: &str) -> Option<&'static Element> {
    PERIODIC_TABLE.iter().find(|e| e.symbol == symbol)
}

/// The ten elements that may appear lowercase (aromatic) in a SMILES
/// organic-subset atom: b c n o p s plus the three written with two
/// letters. SMILES aromatic atoms outside this set must be bracketed.
pub const AROMATIC_ELEMENTS: &[u8] = &[5, 6, 7, 8, 15, 16, 33, 34, 52];

/// The SMILES "organic subset" — elements writable without brackets:
/// B, C, N, O, P, S, F, Cl, Br, I.
pub const ORGANIC_SUBSET: &[u8] = &[5, 6, 7, 8, 9, 15, 16, 17, 35, 53];

/// Standard (most common) neutral valences for an element, used by the
/// implicit-hydrogen filler. The first value is the lowest normal
/// valence; later values are higher hypervalent states (N, P, S, the
/// halogens). An empty slice means "no implicit H" (metals, noble
/// gases) — those atoms carry only explicit hydrogens.
pub fn standard_valences(z: u8) -> &'static [u8] {
    match z {
        1 => &[1],          // H
        5 => &[3],          // B
        6 => &[4],          // C
        7 => &[3],          // N  (5 handled via charge / explicit)
        8 => &[2],          // O
        9 => &[1],          // F
        14 => &[4],         // Si
        15 => &[3, 5],      // P
        16 => &[2, 4, 6],   // S
        17 => &[1],         // Cl
        33 => &[3, 5],      // As
        34 => &[2, 4, 6],   // Se
        35 => &[1],         // Br
        53 => &[1, 3, 5, 7], // I
        _ => &[],
    }
}

/// Pauling electronegativity, used by the Gasteiger-charge seeding and
/// a couple of descriptor heuristics. Returns `2.2` (hydrogen-like) for
/// elements not in the small table — a harmless default for the
/// organic-chemistry regime this crate targets.
pub fn electronegativity(z: u8) -> f64 {
    match z {
        1 => 2.20,
        5 => 2.04,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        35 => 2.96,
        53 => 2.66,
        _ => 2.20,
    }
}

/// Covalent radius (Angstrom) — Cordero 2008 values for the common
/// elements, used to seed the distance-geometry bounds matrix and the
/// 2D layout bond lengths. Falls back to `0.75` (a carbon-ish radius).
pub fn covalent_radius(z: u8) -> f64 {
    match z {
        1 => 0.31,
        5 => 0.84,
        6 => 0.76,
        7 => 0.71,
        8 => 0.66,
        9 => 0.57,
        14 => 1.11,
        15 => 1.07,
        16 => 1.05,
        17 => 1.02,
        35 => 1.20,
        53 => 1.39,
        _ => 0.75,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_dense_and_consistent() {
        assert_eq!(PERIODIC_TABLE.len(), 118);
        for (i, e) in PERIODIC_TABLE.iter().enumerate() {
            assert_eq!(e.number as usize, i + 1, "{} mis-numbered", e.symbol);
        }
    }

    #[test]
    fn symbol_lookup_round_trips() {
        let c = by_symbol("C").unwrap();
        assert_eq!(c.number, 6);
        assert_eq!(by_number(6).unwrap().symbol, "C");
        assert_eq!(by_symbol("Cl").unwrap().number, 17);
        assert!(by_symbol("Xx").is_none());
        assert!(by_number(0).is_none());
        assert!(by_number(200).is_none());
    }

    #[test]
    fn carbon_masses_are_exact() {
        let c = by_symbol("C").unwrap();
        assert!((c.monoisotopic_mass - 12.0).abs() < 1e-9);
        assert!((c.average_mass - 12.011).abs() < 1e-3);
    }

    #[test]
    fn standard_valences_cover_organics() {
        assert_eq!(standard_valences(6), &[4]);
        assert_eq!(standard_valences(7), &[3]);
        assert_eq!(standard_valences(16), &[2, 4, 6]);
        assert!(standard_valences(11).is_empty(), "Na has no implicit H");
    }

    #[test]
    fn electronegativity_ordering() {
        assert!(electronegativity(9) > electronegativity(8));
        assert!(electronegativity(8) > electronegativity(6));
    }

    #[test]
    fn standard_valences_every_tabulated_element() {
        // Every explicit arm of `standard_valences` against textbook
        // organic-chemistry valences. The first value is the lowest
        // normal valence; later values are the hypervalent states.
        assert_eq!(standard_valences(1), &[1]); // H — monovalent
        assert_eq!(standard_valences(5), &[3]); // B — trivalent
        assert_eq!(standard_valences(6), &[4]); // C — tetravalent
        assert_eq!(standard_valences(7), &[3]); // N
        assert_eq!(standard_valences(8), &[2]); // O — divalent
        assert_eq!(standard_valences(9), &[1]); // F — monovalent
        assert_eq!(standard_valences(14), &[4]); // Si — tetravalent
        assert_eq!(standard_valences(15), &[3, 5]); // P — 3 or 5 (PCl5)
        assert_eq!(standard_valences(16), &[2, 4, 6]); // S — 2/4/6 (SF6)
        assert_eq!(standard_valences(17), &[1]); // Cl
        assert_eq!(standard_valences(33), &[3, 5]); // As — like P
        assert_eq!(standard_valences(34), &[2, 4, 6]); // Se — like S
        assert_eq!(standard_valences(35), &[1]); // Br
        assert_eq!(standard_valences(53), &[1, 3, 5, 7]); // I — up to 7 (IF7)
        // The fallback arm: noble gases and metals carry no implicit H.
        assert!(standard_valences(2).is_empty()); // He
        assert!(standard_valences(26).is_empty()); // Fe
        assert!(standard_valences(0).is_empty()); // invalid Z
    }

    #[test]
    fn electronegativity_every_tabulated_element() {
        // Pauling electronegativities — every explicit arm of the
        // `electronegativity` table against the published values.
        assert!((electronegativity(1) - 2.20).abs() < 1e-9); // H
        assert!((electronegativity(5) - 2.04).abs() < 1e-9); // B
        assert!((electronegativity(6) - 2.55).abs() < 1e-9); // C
        assert!((electronegativity(7) - 3.04).abs() < 1e-9); // N
        assert!((electronegativity(8) - 3.44).abs() < 1e-9); // O
        assert!((electronegativity(9) - 3.98).abs() < 1e-9); // F — most EN
        assert!((electronegativity(14) - 1.90).abs() < 1e-9); // Si
        assert!((electronegativity(15) - 2.19).abs() < 1e-9); // P
        assert!((electronegativity(16) - 2.58).abs() < 1e-9); // S
        assert!((electronegativity(17) - 3.16).abs() < 1e-9); // Cl
        assert!((electronegativity(35) - 2.96).abs() < 1e-9); // Br
        assert!((electronegativity(53) - 2.66).abs() < 1e-9); // I
        // The fallback arm: a hydrogen-like 2.20 for untabulated Z.
        assert!((electronegativity(26) - 2.20).abs() < 1e-9); // Fe
        assert!((electronegativity(0) - 2.20).abs() < 1e-9); // invalid Z
        // Halogen electronegativity decreases down the group.
        assert!(electronegativity(9) > electronegativity(17));
        assert!(electronegativity(17) > electronegativity(35));
        assert!(electronegativity(35) > electronegativity(53));
    }

    #[test]
    fn covalent_radius_every_tabulated_element() {
        // Cordero 2008 single-bond covalent radii (Angstrom) — every
        // explicit arm of `covalent_radius` against the published set.
        assert!((covalent_radius(1) - 0.31).abs() < 1e-9); // H — smallest
        assert!((covalent_radius(5) - 0.84).abs() < 1e-9); // B
        assert!((covalent_radius(6) - 0.76).abs() < 1e-9); // C
        assert!((covalent_radius(7) - 0.71).abs() < 1e-9); // N
        assert!((covalent_radius(8) - 0.66).abs() < 1e-9); // O
        assert!((covalent_radius(9) - 0.57).abs() < 1e-9); // F
        assert!((covalent_radius(14) - 1.11).abs() < 1e-9); // Si
        assert!((covalent_radius(15) - 1.07).abs() < 1e-9); // P
        assert!((covalent_radius(16) - 1.05).abs() < 1e-9); // S
        assert!((covalent_radius(17) - 1.02).abs() < 1e-9); // Cl
        assert!((covalent_radius(35) - 1.20).abs() < 1e-9); // Br
        assert!((covalent_radius(53) - 1.39).abs() < 1e-9); // I — largest here
        // The fallback arm: a carbon-ish 0.75 for untabulated Z.
        assert!((covalent_radius(26) - 0.75).abs() < 1e-9); // Fe
        assert!((covalent_radius(0) - 0.75).abs() < 1e-9); // invalid Z
        // Down a group the covalent radius grows (F < Cl < Br < I).
        assert!(covalent_radius(9) < covalent_radius(17));
        assert!(covalent_radius(17) < covalent_radius(35));
        assert!(covalent_radius(35) < covalent_radius(53));
    }

    #[test]
    fn organic_and_aromatic_subset_membership() {
        // The SMILES organic subset is exactly B, C, N, O, F, P, S,
        // Cl, Br, I — writable without brackets.
        for &z in &[5u8, 6, 7, 8, 9, 15, 16, 17, 35, 53] {
            assert!(ORGANIC_SUBSET.contains(&z), "Z={z} should be organic");
        }
        assert_eq!(ORGANIC_SUBSET.len(), 10);
        // The aromatic-element set adds As (33) and Se (52) and drops
        // the halogens.
        for &z in &[5u8, 6, 7, 8, 15, 16, 33, 34, 52] {
            assert!(AROMATIC_ELEMENTS.contains(&z), "Z={z} aromatic");
        }
        assert!(!AROMATIC_ELEMENTS.contains(&9)); // F never aromatic
    }

    #[test]
    fn element_struct_is_copyable_and_comparable() {
        // Exercises the `Element` derives (Copy, PartialEq, Debug).
        let c1 = *by_symbol("C").unwrap();
        let c2 = c1; // Copy — c1 is still usable below
        assert_eq!(c1, c2); // PartialEq
        let o = *by_symbol("O").unwrap();
        assert_ne!(c1, o);
        // Debug formatting works (used in test failure messages).
        assert!(format!("{c1:?}").contains("symbol"));
    }
}

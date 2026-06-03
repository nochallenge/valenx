//! Built-in Gaussian basis-set libraries — STO-3G, 3-21G, 6-31G and
//! 6-31G*.
//!
//! Each library exposes the *raw* (unnormalised-coefficient) shell
//! definitions for hydrogen through neon. [`BasisSet::build`] normalises
//! them on the way in, so the numbers here are the published contraction
//! coefficients exactly as they appear in the EMSL / Basis Set Exchange
//! `Gaussian94` listings.
//!
//! [`BasisSet::build`]: super::BasisSet::build
//!
//! ## What "real v1" means here
//!
//! - **STO-3G** — minimal basis, the full H–Ne table.
//! - **3-21G** — split-valence, the full H–Ne table.
//! - **6-31G** — split-valence with a tighter core, the full H–Ne
//!   table.
//! - **6-31G\*** — 6-31G plus a single *d* polarisation shell on the
//!   first-row atoms (Li–Ne). On hydrogen 6-31G\* equals 6-31G (the
//!   star adds *p* functions only in 6-31G\*\*); this crate follows
//!   that convention.
//!
//! Pople-style *sp* shells are stored split into a separate *s* and
//! *p* [`Shell`] sharing exponents, which is mathematically identical
//! and keeps the integral loop angular-momentum-pure.

use super::{AngularMomentum, Primitive, Shell};
use crate::error::{QchemError, Result};

/// A raw (pre-normalisation) shell definition: an angular momentum and
/// a list of `(exponent, coefficient)` primitive pairs. The `atom_index`
/// and `centre` of a real [`Shell`] are filled in by
/// [`BasisSet::build`](super::BasisSet::build).
#[derive(Clone, Debug)]
pub struct RawShell {
    /// Angular momentum of the shell.
    pub angular: AngularMomentum,
    /// `(exponent, contraction-coefficient)` primitive pairs.
    pub primitives: Vec<Primitive>,
}

impl RawShell {
    fn new(angular: AngularMomentum, prims: &[(f64, f64)]) -> RawShell {
        RawShell {
            angular,
            primitives: prims
                .iter()
                .map(|&(e, c)| Primitive {
                    exponent: e,
                    coefficient: c,
                })
                .collect(),
        }
    }
}

/// Convert a [`RawShell`] into a placed [`Shell`] — used by tests and by
/// callers building custom basis sets atom-by-atom.
pub fn place_shell(raw: &RawShell, atom_index: usize, centre: [f64; 3]) -> Shell {
    Shell {
        atom_index,
        centre,
        angular: raw.angular,
        primitives: raw.primitives.clone(),
    }
}

/// A built-in basis-set library.
pub trait BasisLibrary {
    /// The library's canonical name.
    fn name(&self) -> &'static str;
    /// The raw shells for element `z`, or `None` when the library has
    /// no definition for it.
    fn shells_for(&self, z: u8) -> Option<Vec<RawShell>>;
}

/// Resolve a basis-set name (case-insensitive, common aliases) to its
/// library.
///
/// # Errors
///
/// Returns [`QchemError::Parse`] for an unknown name.
pub fn resolve(name: &str) -> Result<Box<dyn BasisLibrary>> {
    let key: String = name
        .chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| c.to_ascii_lowercase())
        .collect();
    match key.as_str() {
        "sto-3g" | "sto3g" => Ok(Box::new(sto3g::StoNg)),
        "3-21g" | "321g" => Ok(Box::new(pople_321g::Pople321)),
        "6-31g" | "631g" => Ok(Box::new(pople_631g::Pople631)),
        "6-31g*" | "631g*" | "6-31gd" | "631gd" | "6-31g(d)" => {
            Ok(Box::new(pople_631gs::Pople631s))
        }
        _ => Err(QchemError::parse(
            "basis",
            format!("unknown basis set `{name}` (have sto-3g, 3-21g, 6-31g, 6-31g*)"),
        )),
    }
}

/// All basis-set names this crate ships, for UI population.
pub fn available_names() -> &'static [&'static str] {
    &["sto-3g", "3-21g", "6-31g", "6-31g*"]
}

// =====================================================================
// STO-3G
// =====================================================================

/// The STO-3G minimal basis set (Hehre, Stewart, Pople 1969).
pub mod sto3g {
    use super::*;

    /// STO-3G library marker type.
    pub struct StoNg;

    impl BasisLibrary for StoNg {
        fn name(&self) -> &'static str {
            "sto-3g"
        }
        fn shells_for(&self, z: u8) -> Option<Vec<RawShell>> {
            Some(match z {
                1 => vec![RawShell::new(
                    AngularMomentum::S,
                    &[
                        (3.42525091, 0.15432897),
                        (0.62391373, 0.53532814),
                        (0.16885540, 0.44463454),
                    ],
                )],
                2 => vec![RawShell::new(
                    AngularMomentum::S,
                    &[
                        (6.36242139, 0.15432897),
                        (1.15892294, 0.53532814),
                        (0.31364979, 0.44463454),
                    ],
                )],
                3 => sto3g_first_row(
                    &[
                        (16.119575, 0.15432897),
                        (2.9362007, 0.53532814),
                        (0.79465050, 0.44463454),
                    ],
                    &[
                        (0.63628970, -0.09996723),
                        (0.14786033, 0.39951283),
                        (0.04808870, 0.70011547),
                    ],
                    &[
                        (0.63628970, 0.15591627),
                        (0.14786033, 0.60768372),
                        (0.04808870, 0.39195739),
                    ],
                ),
                4 => sto3g_first_row(
                    &[
                        (30.167871, 0.15432897),
                        (5.4951153, 0.53532814),
                        (1.4871927, 0.44463454),
                    ],
                    &[
                        (1.3148331, -0.09996723),
                        (0.3055389, 0.39951283),
                        (0.09937070, 0.70011547),
                    ],
                    &[
                        (1.3148331, 0.15591627),
                        (0.3055389, 0.60768372),
                        (0.09937070, 0.39195739),
                    ],
                ),
                5 => sto3g_first_row(
                    &[
                        (48.791113, 0.15432897),
                        (8.8873622, 0.53532814),
                        (2.4052670, 0.44463454),
                    ],
                    &[
                        (2.2369561, -0.09996723),
                        (0.5198205, 0.39951283),
                        (0.1690618, 0.70011547),
                    ],
                    &[
                        (2.2369561, 0.15591627),
                        (0.5198205, 0.60768372),
                        (0.1690618, 0.39195739),
                    ],
                ),
                6 => sto3g_first_row(
                    &[
                        (71.616837, 0.15432897),
                        (13.045096, 0.53532814),
                        (3.5305122, 0.44463454),
                    ],
                    &[
                        (2.9412494, -0.09996723),
                        (0.6834831, 0.39951283),
                        (0.2222899, 0.70011547),
                    ],
                    &[
                        (2.9412494, 0.15591627),
                        (0.6834831, 0.60768372),
                        (0.2222899, 0.39195739),
                    ],
                ),
                7 => sto3g_first_row(
                    &[
                        (99.106169, 0.15432897),
                        (18.052312, 0.53532814),
                        (4.8856602, 0.44463454),
                    ],
                    &[
                        (3.7804559, -0.09996723),
                        (0.8784966, 0.39951283),
                        (0.2857144, 0.70011547),
                    ],
                    &[
                        (3.7804559, 0.15591627),
                        (0.8784966, 0.60768372),
                        (0.2857144, 0.39195739),
                    ],
                ),
                8 => sto3g_first_row(
                    &[
                        (130.70932, 0.15432897),
                        (23.808861, 0.53532814),
                        (6.4436083, 0.44463454),
                    ],
                    &[
                        (5.0331513, -0.09996723),
                        (1.1695961, 0.39951283),
                        (0.3803890, 0.70011547),
                    ],
                    &[
                        (5.0331513, 0.15591627),
                        (1.1695961, 0.60768372),
                        (0.3803890, 0.39195739),
                    ],
                ),
                9 => sto3g_first_row(
                    &[
                        (166.67913, 0.15432897),
                        (30.360812, 0.53532814),
                        (8.2168207, 0.44463454),
                    ],
                    &[
                        (6.4648032, -0.09996723),
                        (1.5022812, 0.39951283),
                        (0.4885885, 0.70011547),
                    ],
                    &[
                        (6.4648032, 0.15591627),
                        (1.5022812, 0.60768372),
                        (0.4885885, 0.39195739),
                    ],
                ),
                10 => sto3g_first_row(
                    &[
                        (207.01561, 0.15432897),
                        (37.708151, 0.53532814),
                        (10.205297, 0.44463454),
                    ],
                    &[
                        (8.2463151, -0.09996723),
                        (1.9162662, 0.39951283),
                        (0.6232293, 0.70011547),
                    ],
                    &[
                        (8.2463151, 0.15591627),
                        (1.9162662, 0.60768372),
                        (0.6232293, 0.39195739),
                    ],
                ),
                _ => return None,
            })
        }
    }

    /// Assemble the three STO-3G shells of a first-row atom: a `1s`
    /// core, a `2s` valence shell and a `2p` valence shell. The `2s`
    /// and `2p` share their exponents (the Pople `sp` convention).
    fn sto3g_first_row(
        s1: &[(f64, f64)],
        s2: &[(f64, f64)],
        p2: &[(f64, f64)],
    ) -> Vec<RawShell> {
        vec![
            RawShell::new(AngularMomentum::S, s1),
            RawShell::new(AngularMomentum::S, s2),
            RawShell::new(AngularMomentum::P, p2),
        ]
    }
}

// =====================================================================
// 3-21G
// =====================================================================

/// The 3-21G split-valence basis set (Binkley, Pople, Hehre 1980).
pub mod pople_321g {
    use super::*;

    /// 3-21G library marker type.
    pub struct Pople321;

    impl BasisLibrary for Pople321 {
        fn name(&self) -> &'static str {
            "3-21g"
        }
        fn shells_for(&self, z: u8) -> Option<Vec<RawShell>> {
            Some(match z {
                1 => vec![
                    RawShell::new(
                        AngularMomentum::S,
                        &[
                            (5.44717800, 0.15628498),
                            (0.82454724, 0.90469091),
                        ],
                    ),
                    RawShell::new(AngularMomentum::S, &[(0.18319158, 1.0)]),
                ],
                2 => vec![
                    RawShell::new(
                        AngularMomentum::S,
                        &[
                            (13.6267000, 0.17523000),
                            (1.99935000, 0.89348300),
                        ],
                    ),
                    RawShell::new(AngularMomentum::S, &[(0.38299300, 1.0)]),
                ],
                3 => split_valence_first_row(
                    &[
                        (36.8382000, 0.0696868),
                        (5.4817200, 0.3813460),
                        (1.1117100, 0.6817020),
                    ],
                    &[(0.54020500, -0.2631270), (0.10225500, 1.1433900)],
                    &[(0.54020500, 0.1615460), (0.10225500, 0.9156630)],
                    &[(0.02856500, 1.0)],
                    &[(0.02856500, 1.0)],
                ),
                4 => split_valence_first_row(
                    &[
                        (71.8876000, 0.0644263),
                        (10.7289000, 0.3660960),
                        (2.2220500, 0.6959340),
                    ],
                    &[(1.29548000, -0.4210640), (0.26881000, 1.2240700)],
                    &[(1.29548000, 0.2051320), (0.26881000, 0.8825280)],
                    &[(0.07735000, 1.0)],
                    &[(0.07735000, 1.0)],
                ),
                5 => split_valence_first_row(
                    &[
                        (116.4340000, 0.0629605),
                        (17.4315000, 0.3633040),
                        (3.6801600, 0.6972550),
                    ],
                    &[(2.28187000, -0.3686620), (0.46524800, 1.1994400)],
                    &[(2.28187000, 0.2311520), (0.46524800, 0.8667640)],
                    &[(0.12432800, 1.0)],
                    &[(0.12432800, 1.0)],
                ),
                6 => split_valence_first_row(
                    &[
                        (172.2560000, 0.0617667),
                        (25.9109000, 0.3587940),
                        (5.5335000, 0.7007130),
                    ],
                    &[(3.66498000, -0.3958970), (0.77054500, 1.2158400)],
                    &[(3.66498000, 0.2364600), (0.77054500, 0.8606190)],
                    &[(0.19585700, 1.0)],
                    &[(0.19585700, 1.0)],
                ),
                7 => split_valence_first_row(
                    &[
                        (242.7660000, 0.0598657),
                        (36.4851000, 0.3529550),
                        (7.8144900, 0.7065130),
                    ],
                    &[(5.42522000, -0.4133010), (1.14915000, 1.2244200)],
                    &[(5.42522000, 0.2379720), (1.14915000, 0.8589530)],
                    &[(0.28320500, 1.0)],
                    &[(0.28320500, 1.0)],
                ),
                8 => split_valence_first_row(
                    &[
                        (322.0370000, 0.0592394),
                        (48.4308000, 0.3515000),
                        (10.4206000, 0.7076580),
                    ],
                    &[(7.40294000, -0.4044530), (1.57620000, 1.2215600)],
                    &[(7.40294000, 0.2445860), (1.57620000, 0.8539550)],
                    &[(0.37368400, 1.0)],
                    &[(0.37368400, 1.0)],
                ),
                9 => split_valence_first_row(
                    &[
                        (413.8010000, 0.0585483),
                        (62.2246000, 0.3493080),
                        (13.4340000, 0.7096320),
                    ],
                    &[(9.77759000, -0.4073270), (2.08351000, 1.2231400)],
                    &[(9.77759000, 0.2466800), (2.08351000, 0.8517430)],
                    &[(0.48238300, 1.0)],
                    &[(0.48238300, 1.0)],
                ),
                10 => split_valence_first_row(
                    &[
                        (515.7240000, 0.0581430),
                        (77.6178000, 0.3479510),
                        (16.8121000, 0.7107140),
                    ],
                    &[(12.4830000, -0.4099220), (2.66451000, 1.2243100)],
                    &[(12.4830000, 0.2474600), (2.66451000, 0.8517430)],
                    &[(0.60625000, 1.0)],
                    &[(0.60625000, 1.0)],
                ),
                _ => return None,
            })
        }
    }
}

// =====================================================================
// 6-31G
// =====================================================================

/// The 6-31G split-valence basis set (Hehre, Ditchfield, Pople 1972;
/// Hariharan & Pople 1973).
pub mod pople_631g {
    use super::*;

    /// 6-31G library marker type.
    pub struct Pople631;

    impl BasisLibrary for Pople631 {
        fn name(&self) -> &'static str {
            "6-31g"
        }
        fn shells_for(&self, z: u8) -> Option<Vec<RawShell>> {
            raw_631g(z)
        }
    }
}

// =====================================================================
// 6-31G*
// =====================================================================

/// The 6-31G\* basis set — 6-31G plus a *d* polarisation shell on the
/// first-row atoms (Hariharan & Pople 1973).
pub mod pople_631gs {
    use super::*;

    /// 6-31G\* library marker type.
    pub struct Pople631s;

    /// Polarisation `d`-shell exponents for Li–Ne in 6-31G\*. A single
    /// uncontracted `d` primitive per atom.
    const D_POLARISATION: [(u8, f64); 8] = [
        (3, 0.200),
        (4, 0.400),
        (5, 0.600),
        (6, 0.800),
        (7, 0.800),
        (8, 0.800),
        (9, 0.800),
        (10, 0.800),
    ];

    impl BasisLibrary for Pople631s {
        fn name(&self) -> &'static str {
            "6-31g*"
        }
        fn shells_for(&self, z: u8) -> Option<Vec<RawShell>> {
            let mut shells = raw_631g(z)?;
            // 6-31G* adds d functions only on first-row atoms; on H/He
            // it is identical to 6-31G (the d→p H polarisation is
            // 6-31G**).
            if let Some(&(_, exp)) = D_POLARISATION.iter().find(|&&(zz, _)| zz == z) {
                shells.push(RawShell::new(AngularMomentum::D, &[(exp, 1.0)]));
            }
            Some(shells)
        }
    }
}

/// The raw 6-31G shells for element `z` (`1..=10`). Shared by the 6-31G
/// and 6-31G\* libraries.
fn raw_631g(z: u8) -> Option<Vec<RawShell>> {
    Some(match z {
        1 => vec![
            RawShell::new(
                AngularMomentum::S,
                &[
                    (18.73113700, 0.03349460),
                    (2.82539437, 0.23472695),
                    (0.64012170, 0.81375733),
                ],
            ),
            RawShell::new(AngularMomentum::S, &[(0.16127778, 1.0)]),
        ],
        2 => vec![
            RawShell::new(
                AngularMomentum::S,
                &[
                    (38.42163400, 0.02376600),
                    (5.77803000, 0.15467900),
                    (1.24177400, 0.46963000),
                ],
            ),
            RawShell::new(AngularMomentum::S, &[(0.29796400, 1.0)]),
        ],
        3 => split_valence_first_row(
            &[
                (642.4189200, 0.002142607),
                (96.7985150, 0.016208872),
                (22.0911210, 0.077315575),
                (6.2010703, 0.24578605),
                (1.9351177, 0.47018900),
                (0.63673578, 0.34547043),
            ],
            &[
                (2.324918408, -0.03509174),
                (0.6324372080, -0.1912328),
                (0.07905344372, 1.0839875),
            ],
            &[
                (2.324918408, 0.008941508),
                (0.6324372080, 0.1410094),
                (0.07905344372, 0.9453638),
            ],
            &[(0.03596202837, 1.0)],
            &[(0.03596202837, 1.0)],
        ),
        4 => split_valence_first_row(
            &[
                (1264.5856900, 0.00194475),
                (189.9363680, 0.01483505),
                (43.1590890, 0.07209058),
                (12.0986627, 0.2377154),
                (3.8063232, 0.4698746),
                (1.2728903, 0.3567452),
            ],
            &[
                (3.196446980, -0.1126470),
                (0.747813319, -0.2297614),
                (0.219575629, 1.1868760),
            ],
            &[
                (3.196446980, 0.0559802),
                (0.747813319, 0.2615511),
                (0.219575629, 0.7939842),
            ],
            &[(0.08233573900, 1.0)],
            &[(0.08233573900, 1.0)],
        ),
        5 => split_valence_first_row(
            &[
                (2068.8822000, 0.001866274),
                (310.6495700, 0.014251487),
                (70.6830030, 0.069564797),
                (19.8610803, 0.23292699),
                (6.2993048, 0.46703612),
                (2.1270270, 0.36342178),
            ],
            &[
                (4.727971071, -0.1303938),
                (1.190337771, -0.1307823),
                (0.359292671, 1.1309444),
            ],
            &[
                (4.727971071, 0.07459760),
                (1.190337771, 0.30798472),
                (0.359292671, 0.74358494),
            ],
            &[(0.1280885454, 1.0)],
            &[(0.1280885454, 1.0)],
        ),
        6 => split_valence_first_row(
            &[
                (3047.5248800, 0.00183473713),
                (457.3695180, 0.0140373228),
                (103.9486850, 0.0688426222),
                (29.2101553, 0.232184443),
                (9.2866630, 0.467941348),
                (3.1639270, 0.362311985),
            ],
            &[
                (7.868272350, -0.1193324198),
                (1.881288540, -0.1608541519),
                (0.5442492580, 1.143456438),
            ],
            &[
                (7.868272350, 0.06899906659),
                (1.881288540, 0.3164239609),
                (0.5442492580, 0.7443083963),
            ],
            &[(0.1687144785, 1.0)],
            &[(0.1687144785, 1.0)],
        ),
        7 => split_valence_first_row(
            &[
                (4173.5114600, 0.00183477216),
                (627.4579110, 0.0139946270),
                (142.9020930, 0.0685865513),
                (40.2343293, 0.232271435),
                (12.8202129, 0.469059752),
                (4.3904370, 0.360455199),
            ],
            &[
                (11.626358140, -0.1149611198),
                (2.716280007, -0.1691176588),
                (0.7722183700, 1.145851947),
            ],
            &[
                (11.626358140, 0.06758000564),
                (2.716280007, 0.3239072703),
                (0.7722183700, 0.7408953393),
            ],
            &[(0.2120315449, 1.0)],
            &[(0.2120315449, 1.0)],
        ),
        8 => split_valence_first_row(
            &[
                (5484.6716600, 0.00183107443),
                (825.2349460, 0.0139501724),
                (188.0469580, 0.0684450785),
                (52.9645000, 0.232714336),
                (16.8975704, 0.470192898),
                (5.7996353, 0.358520853),
            ],
            &[
                (15.539616250, -0.1107775495),
                (3.599933586, -0.1480262633),
                (1.013761750, 1.130767015),
            ],
            &[
                (15.539616250, 0.07087426882),
                (3.599933586, 0.3397528394),
                (1.013761750, 0.7271585773),
            ],
            &[(0.2700058230, 1.0)],
            &[(0.2700058230, 1.0)],
        ),
        9 => split_valence_first_row(
            &[
                (7001.7130900, 0.00181961690),
                (1051.3660900, 0.0139160796),
                (239.4285690, 0.0684053246),
                (67.3974453, 0.233185760),
                (21.5195739, 0.471872724),
                (8.2163984, 0.356866300),
            ],
            &[
                (20.847953280, -0.1085069753),
                (4.808308390, -0.1464516584),
                (1.344063900, 1.128688581),
            ],
            &[
                (20.847953280, 0.07162872405),
                (4.808308390, 0.3459126772),
                (1.344063900, 0.7224693995),
            ],
            &[(0.3588151270, 1.0)],
            &[(0.3588151270, 1.0)],
        ),
        10 => split_valence_first_row(
            &[
                (8425.8515300, 0.00188434850),
                (1268.5194000, 0.0143368296),
                (289.6214140, 0.0701096147),
                (81.8590040, 0.237373266),
                (26.2515079, 0.473007126),
                (9.0947204, 0.348401241),
            ],
            &[
                (26.532131000, -0.1071182872),
                (6.101755010, -0.1461638213),
                (1.696271530, 1.127773504),
            ],
            &[
                (26.532131000, 0.07190958856),
                (6.101755010, 0.3495136043),
                (1.696271530, 0.7199405129),
            ],
            &[(0.4485647750, 1.0)],
            &[(0.4485647750, 1.0)],
        ),
        _ => return None,
    })
}

/// Assemble the five split-valence shells of a first-row atom: a `1s`
/// core, an inner valence `s`, an inner valence `p`, an outer valence
/// `s` and an outer valence `p`. The inner `s`/`p` share their
/// exponents, as do the outer pair (the Pople `sp` convention).
fn split_valence_first_row(
    core_s: &[(f64, f64)],
    inner_s: &[(f64, f64)],
    inner_p: &[(f64, f64)],
    outer_s: &[(f64, f64)],
    outer_p: &[(f64, f64)],
) -> Vec<RawShell> {
    vec![
        RawShell::new(AngularMomentum::S, core_s),
        RawShell::new(AngularMomentum::S, inner_s),
        RawShell::new(AngularMomentum::P, inner_p),
        RawShell::new(AngularMomentum::S, outer_s),
        RawShell::new(AngularMomentum::P, outer_p),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_known_names() {
        assert_eq!(resolve("STO-3G").unwrap().name(), "sto-3g");
        assert_eq!(resolve("3-21g").unwrap().name(), "3-21g");
        assert_eq!(resolve("6-31G").unwrap().name(), "6-31g");
        assert_eq!(resolve("6-31g*").unwrap().name(), "6-31g*");
        assert_eq!(resolve("6-31G(d)").unwrap().name(), "6-31g*");
    }

    #[test]
    fn resolve_unknown_errors() {
        assert!(resolve("def2-tzvp").is_err());
    }

    #[test]
    fn sto3g_hydrogen_has_one_s_shell() {
        let shells = sto3g::StoNg.shells_for(1).unwrap();
        assert_eq!(shells.len(), 1);
        assert_eq!(shells[0].angular, AngularMomentum::S);
        assert_eq!(shells[0].primitives.len(), 3);
    }

    #[test]
    fn sto3g_carbon_has_three_shells() {
        let shells = sto3g::StoNg.shells_for(6).unwrap();
        // 1s, 2s, 2p.
        assert_eq!(shells.len(), 3);
        assert_eq!(shells[2].angular, AngularMomentum::P);
    }

    #[test]
    fn pople_321g_hydrogen_split_into_two_shells() {
        let shells = pople_321g::Pople321.shells_for(1).unwrap();
        assert_eq!(shells.len(), 2);
    }

    #[test]
    fn pople_631g_carbon_has_five_shells() {
        let shells = pople_631g::Pople631.shells_for(6).unwrap();
        // core 1s, inner 2s, inner 2p, outer 2s, outer 2p.
        assert_eq!(shells.len(), 5);
        // core s is the 6-primitive contraction.
        assert_eq!(shells[0].primitives.len(), 6);
    }

    #[test]
    fn pople_631gs_carbon_adds_a_d_shell() {
        let plain = pople_631g::Pople631.shells_for(6).unwrap();
        let star = pople_631gs::Pople631s.shells_for(6).unwrap();
        assert_eq!(star.len(), plain.len() + 1);
        assert_eq!(star.last().unwrap().angular, AngularMomentum::D);
    }

    #[test]
    fn pople_631gs_hydrogen_equals_631g() {
        let plain = pople_631g::Pople631.shells_for(1).unwrap();
        let star = pople_631gs::Pople631s.shells_for(1).unwrap();
        assert_eq!(star.len(), plain.len());
    }

    #[test]
    fn all_libraries_cover_h_through_ne() {
        for name in available_names() {
            let lib = resolve(name).unwrap();
            for z in 1..=10u8 {
                assert!(
                    lib.shells_for(z).is_some(),
                    "{name} missing element Z={z}"
                );
            }
            // Beyond neon is out of v1 scope.
            assert!(lib.shells_for(11).is_none(), "{name} should stop at Ne");
        }
    }
}

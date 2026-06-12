//! Panel 14 — **RNA Designer** — an end-to-end RNA / mRNA design
//! workbench over [`valenx_rnastruct`] + [`valenx_rnadesign`] +
//! [`valenx_genediting`].
//!
//! Where the [`super::rnastruct`] panel is a single-purpose folder, this
//! panel is the design *workbench*: five linked sections that take you
//! from a sequence (or a target structure, or a protein) all the way to
//! a synthesis-ready mRNA construct. Every computation is a native crate
//! call — the panel never reimplements folding or design.
//!
//! ## The five sections
//!
//! 1. **Fold a structure** — fold an RNA sequence (exact Zuker for short
//!    RNA, LinearFold for long, or auto by length); the MFE structure,
//!    MFE energy, the LinearPartition ensemble free energy and the
//!    Boltzmann frequency of the MFE structure.
//! 2. **Structure visualization** — the predicted secondary structure
//!    drawn as a real **2-D diagram** with egui's painter (backbone +
//!    bases + base-pair bonds, from the `valenx-rnastruct` naview-class
//!    2-D layout), plus a mountain plot and a base-pair-probability
//!    dot-plot heatmap.
//! 3. **Inverse design** — ensemble-defect inverse folding from a target
//!    dot-bracket, with the constrained-design controls (locked
//!    positions, GC band, forbidden motifs); the designed sequence, its
//!    ensemble defect, and a fold-back confirming it matches the target.
//! 4. **mRNA design (LinearDesign)** — the LinearDesign joint mRNA
//!    optimiser: a protein → an optimised CDS over `MFE + λ·codon-penalty`
//!    with a λ slider, plus a λ-sweep that plots the stability/efficiency
//!    Pareto front.
//! 5. **mRNA construct** — wrap the optimised CDS into a full five-part
//!    construct (5′UTR + Kozak, CDS, 3′UTR, poly-A, cap) via
//!    `valenx-genediting`.
//!
//! ## Threading
//!
//! Inverse folding, LinearDesign on a long protein and the λ-sweep can
//! each take several seconds. They run on a background `std::thread`
//! (the pattern the Aerodynamics workbench uses); the panel polls a
//! handle once per frame and shows a busy spinner so the egui thread
//! never blocks.

use std::thread::JoinHandle;

use eframe::egui;
use egui_plot::{Line, MarkerShape, Plot, PlotPoints, Points};

use valenx_genediting::mrna::codon::ExpressionHost;
use valenx_genediting::mrna::construct::{MrnaConstruct, MrnaConstructBuilder};
use valenx_genediting::mrna::tailcap::{recommend_cap, recommend_poly_a, MrnaUseCase};
use valenx_genediting::mrna::utr::{reference_utr3, reference_utr5};
use valenx_rnadesign::constraints::{lock_entry, DesignConstraintSet};
use valenx_rnadesign::goal::DesignConstraints;
use valenx_rnadesign::inverse::{inverse_fold_constrained, EnsembleDefectParams};
use valenx_rnadesign::lineardesign::{
    linear_design, pareto_sweep, LinearDesignRequest, LinearDesignResult, ParetoPoint,
    DEFAULT_BEAM_SIZE,
};
use valenx_rnastruct::fold::energy::{GAS_CONSTANT, T37_KELVIN};
use valenx_rnastruct::fold::linear::fold_linear;
use valenx_rnastruct::fold::zuker::mfe;
use valenx_rnastruct::layout::{layout, Layout};
use valenx_rnastruct::{linear_partition, RnaSeq, Structure};

use super::common;
use crate::ValenxApp;

// ---------------------------------------------------------------------
// Workbench sections
// ---------------------------------------------------------------------

/// One of the five workbench sections — the top selector.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Section {
    /// Fold an RNA sequence.
    #[default]
    Fold,
    /// Visualise the folded structure (2-D diagram, mountain, dot-plot).
    Visualize,
    /// Ensemble-defect inverse design from a target structure.
    Inverse,
    /// LinearDesign joint mRNA optimisation.
    Mrna,
    /// Wrap the optimised CDS into a full mRNA construct.
    Construct,
}

impl Section {
    /// Every section, in workbench order.
    pub const ALL: [Section; 5] = [
        Section::Fold,
        Section::Visualize,
        Section::Inverse,
        Section::Mrna,
        Section::Construct,
    ];

    /// A short tab label.
    pub fn label(self) -> &'static str {
        match self {
            Section::Fold => "1 · Fold",
            Section::Visualize => "2 · Visualize",
            Section::Inverse => "3 · Inverse design",
            Section::Mrna => "4 · mRNA (LinearDesign)",
            Section::Construct => "5 · mRNA construct",
        }
    }

    /// The section heading.
    pub fn title(self) -> &'static str {
        match self {
            Section::Fold => "Fold a structure",
            Section::Visualize => "Structure visualization",
            Section::Inverse => "Inverse design — ensemble-defect folding",
            Section::Mrna => "mRNA design — the LinearDesign joint optimiser",
            Section::Construct => "mRNA construct assembly",
        }
    }

    /// A one-line subtitle.
    pub fn subtitle(self) -> &'static str {
        match self {
            Section::Fold => "Predict the secondary structure of an RNA sequence.",
            Section::Visualize => {
                "The folded structure as a 2-D diagram, a mountain plot and a base-pair dot-plot."
            }
            Section::Inverse => "Design a sequence that folds to a target structure.",
            Section::Mrna => "Optimise a coding sequence jointly for stability and codon usage.",
            Section::Construct => "Wrap the optimised CDS into a synthesis-ready five-part mRNA.",
        }
    }
}

// ---------------------------------------------------------------------
// Fold engine choice
// ---------------------------------------------------------------------

/// Which folding engine the Fold section uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum FoldEngine {
    /// Pick the engine by sequence length (Zuker ≤ 200 nt, else
    /// LinearFold).
    #[default]
    Auto,
    /// Exact Zuker minimum-free-energy folding (`O(n³)`).
    Zuker,
    /// LinearFold linear-time beam-search folding.
    LinearFold,
}

impl FoldEngine {
    /// A short label.
    pub fn label(self) -> &'static str {
        match self {
            FoldEngine::Auto => "Auto (by length)",
            FoldEngine::Zuker => "Exact Zuker MFE",
            FoldEngine::LinearFold => "LinearFold (linear-time)",
        }
    }

    /// The length at / above which `Auto` switches to LinearFold.
    pub const AUTO_LINEAR_THRESHOLD: usize = 200;

    /// Resolves `Auto` to a concrete engine for a sequence of `len` nt.
    pub fn resolve(self, len: usize) -> FoldEngine {
        match self {
            FoldEngine::Auto => {
                if len >= Self::AUTO_LINEAR_THRESHOLD {
                    FoldEngine::LinearFold
                } else {
                    FoldEngine::Zuker
                }
            }
            other => other,
        }
    }
}

// ---------------------------------------------------------------------
// Fold result
// ---------------------------------------------------------------------

/// The outcome of a Fold-section run — every figure the panel shows,
/// plus the `Structure` the Visualize section consumes.
#[derive(Clone, Debug)]
pub struct FoldOutcome {
    /// The folded RNA sequence (`A C G U`).
    pub sequence: String,
    /// The concrete engine that produced this fold.
    pub engine: FoldEngine,
    /// The MFE secondary structure.
    pub structure: Structure,
    /// The dot-bracket of [`structure`](Self::structure).
    pub dot_bracket: String,
    /// The minimum free energy, kcal/mol.
    pub mfe: f64,
    /// The LinearPartition ensemble free energy `G`, kcal/mol.
    pub ensemble_g: f64,
    /// The Boltzmann frequency of the MFE structure in `[0, 1]`.
    pub mfe_frequency: f64,
    /// Base pairs with probability ≥ the dot-plot threshold, as
    /// `(i, j, p)` — for the Visualize section's dot-plot / heatmap.
    pub bpp: Vec<(usize, usize, f64)>,
    /// The peak base-pair probability across the ensemble.
    pub max_bpp: f64,
}

/// Folds `seq_text` with `engine`, returning a [`FoldOutcome`] or an
/// error message. Extracted from the Run-button closure so the headless
/// UI tests can call it directly.
///
/// The ensemble free energy and the base-pair probabilities always come
/// from **LinearPartition** (it scales to long mRNA); the MFE structure
/// comes from the chosen MFE engine.
pub fn run_fold(seq_text: &str, engine: FoldEngine) -> Result<FoldOutcome, String> {
    let cleaned = common::clean_sequence(seq_text);
    if cleaned.is_empty() {
        return Err("enter an RNA sequence".to_string());
    }
    let rna = RnaSeq::parse(cleaned.as_bytes()).map_err(|e| e.to_string())?;
    let resolved = engine.resolve(rna.len());

    // MFE structure + energy from the chosen engine.
    let (structure, mfe_energy) = match resolved {
        FoldEngine::Zuker => {
            let r = mfe(&rna).map_err(|e| e.to_string())?;
            (r.structure, r.energy)
        }
        FoldEngine::LinearFold => {
            let r = fold_linear(&rna).map_err(|e| e.to_string())?;
            (r.structure, r.energy)
        }
        FoldEngine::Auto => unreachable!("resolve() never returns Auto"),
    };

    // Ensemble characterisation from LinearPartition.
    let lp = linear_partition(&rna).map_err(|e| e.to_string())?;
    let ensemble_g = lp.ensemble_free_energy();

    // MFE Boltzmann frequency: exp(-(E_mfe - G)/RT).
    let rt = GAS_CONSTANT * T37_KELVIN;
    let mfe_frequency = (-(mfe_energy - ensemble_g) / rt).exp().clamp(0.0, 1.0);

    // Significant base pairs for the dot-plot.
    let mut bpp = lp.significant_pairs(0.0);
    bpp.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    let max_bpp = bpp.iter().map(|&(_, _, p)| p).fold(0.0_f64, f64::max);

    let dot_bracket = structure.to_dot_bracket();
    Ok(FoldOutcome {
        sequence: rna.as_str().to_string(),
        engine: resolved,
        structure,
        dot_bracket,
        mfe: mfe_energy,
        ensemble_g,
        mfe_frequency,
        bpp,
        max_bpp,
    })
}

// ---------------------------------------------------------------------
// Inverse-design run
// ---------------------------------------------------------------------

/// The outcome of an Inverse-design run.
#[derive(Clone, Debug)]
pub struct InverseOutcome {
    /// The target dot-bracket the design was steered toward.
    pub target_db: String,
    /// The designed RNA sequence (`A C G U`).
    pub sequence: String,
    /// The achieved ensemble defect (expected mis-(un)paired bases).
    pub ensemble_defect: f64,
    /// The normalised ensemble defect in `[0, 1]`.
    pub normalized_defect: f64,
    /// `true` if the normalised-defect target was reached.
    pub solved: bool,
    /// The base-pair distance between the design's MFE fold and the
    /// target — `0` means the single MFE structure *is* the target.
    pub mfe_distance: usize,
    /// The dot-bracket of the design's own MFE fold (the fold-back).
    pub fold_back_db: String,
    /// Mutation steps accepted / attempted.
    pub accepted: usize,
    /// Mutation steps attempted.
    pub total: usize,
    /// Whether a constraint set was applied.
    pub constrained: bool,
    /// The number of locked positions honoured.
    pub locked_count: usize,
}

/// Runs an ensemble-defect inverse-folding design for `target_db`,
/// honouring the constraint form. Extracted so the headless UI tests can
/// call it directly. Always uses [`inverse_fold_constrained`] — an empty
/// constraint form is just an all-free, full-GC-band set, so the one
/// code path covers both.
pub fn run_inverse(
    target_db: &str,
    gc_min: f64,
    gc_max: f64,
    max_homopolymer: usize,
    forbidden_motifs: &str,
    locked: &[(usize, char)],
) -> Result<InverseOutcome, String> {
    let db = target_db.trim();
    if db.is_empty() {
        return Err("enter a target dot-bracket structure".to_string());
    }
    let target = Structure::from_dot_bracket(db).map_err(|e| e.to_string())?;
    let n = target.len();
    if n == 0 {
        return Err("the target structure is empty".to_string());
    }

    // Build the declarative constraints, then the active constraint set.
    let mut constraints = DesignConstraints::default().with_gc_range(gc_min, gc_max);
    constraints.max_homopolymer = max_homopolymer.max(1);
    for m in split_list(forbidden_motifs) {
        constraints = constraints.forbid_motif(&m);
    }
    // Locked positions are stored as required_subsequences "@idx=base".
    let mut required: Vec<String> = Vec::new();
    for &(idx, base) in locked {
        if idx >= n {
            return Err(format!("locked position {idx} is past the {n}-nt target"));
        }
        required.push(lock_entry(idx, base as u8));
    }
    constraints.required_subsequences = required;
    constraints.validate().map_err(|e| e.to_string())?;

    let set = DesignConstraintSet::new(&constraints, n);
    let locked_count = set.locked_count();
    let design = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default())
        .map_err(|e| format!("[{}] {e}", e.code()))?;

    // Fold the design back to confirm it matches the target.
    let fold_back = mfe(&RnaSeq::parse(&design.sequence).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;

    Ok(InverseOutcome {
        target_db: db.to_string(),
        sequence: design.sequence_str().to_string(),
        ensemble_defect: design.ensemble_defect,
        normalized_defect: design.normalized_defect,
        solved: design.solved,
        mfe_distance: design.mfe_distance,
        fold_back_db: fold_back.structure.to_dot_bracket(),
        accepted: design.accepted_steps,
        total: design.total_steps,
        constrained: locked_count > 0
            || !constraints.forbidden_motifs.is_empty()
            || gc_min > 0.0
            || gc_max < 1.0,
        locked_count,
    })
}

// ---------------------------------------------------------------------
// LinearDesign run
// ---------------------------------------------------------------------

/// The outcome of an mRNA-design (LinearDesign) run.
#[derive(Clone, Debug)]
pub struct MrnaOutcome {
    /// The optimised coding sequence (RNA, `AUG`…stop).
    pub cds: String,
    /// The predicted minimum free energy of the CDS, kcal/mol.
    pub mfe: f64,
    /// The codon adaptation index of the CDS, `(0, 1]`.
    pub cai: f64,
    /// The dot-bracket of the CDS's predicted structure.
    pub dot_bracket: String,
    /// The trade-off weight `λ` used.
    pub lambda: f64,
    /// `true` if the lattice DP was exact (no beam pruning).
    pub exact: bool,
    /// The combined objective `MFE + λ·codon-penalty`.
    pub objective: f64,
    /// The expression host.
    pub host: ExpressionHost,
}

/// Runs the LinearDesign joint optimiser for `protein` at trade-off
/// `lambda`. Extracted so the headless UI tests can call it directly.
pub fn run_linear_design(
    protein: &str,
    host: ExpressionHost,
    lambda: f64,
) -> Result<MrnaOutcome, String> {
    let cleaned = common::clean_sequence(protein);
    if cleaned.is_empty() {
        return Err("paste a protein sequence (one-letter codes)".to_string());
    }
    let req = LinearDesignRequest {
        protein: cleaned.into_bytes(),
        host,
        lambda,
        beam_size: DEFAULT_BEAM_SIZE,
    };
    let r: LinearDesignResult = linear_design(&req).map_err(|e| format!("[{}] {e}", e.code()))?;
    Ok(MrnaOutcome {
        cds: r.cds_str().to_string(),
        mfe: r.mfe,
        cai: r.cai,
        dot_bracket: r.structure.to_dot_bracket(),
        lambda: r.lambda,
        exact: r.exact,
        objective: r.objective,
        host,
    })
}

/// The λ values the Pareto sweep evaluates — from a pure-MFE design at
/// `λ = 0` through a strongly CAI-weighted design.
const SWEEP_LAMBDAS: [f64; 7] = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0];

/// Runs a λ-sweep for `protein`, returning the stability/efficiency
/// Pareto front. Extracted so the headless UI tests can call it directly.
pub fn run_pareto_sweep(protein: &str, host: ExpressionHost) -> Result<Vec<ParetoPoint>, String> {
    let cleaned = common::clean_sequence(protein);
    if cleaned.is_empty() {
        return Err("paste a protein sequence (one-letter codes)".to_string());
    }
    pareto_sweep(cleaned.as_bytes(), host, &SWEEP_LAMBDAS, DEFAULT_BEAM_SIZE)
        .map_err(|e| format!("[{}] {e}", e.code()))
}

// ---------------------------------------------------------------------
// Construct run
// ---------------------------------------------------------------------

/// The outcome of an mRNA-construct assembly.
#[derive(Clone, Debug)]
pub struct ConstructOutcome {
    /// The assembled, validated construct.
    pub construct: MrnaConstruct,
    /// The full transcript body (5′UTR + CDS + 3′UTR + poly-A), RNA.
    pub transcript: String,
    /// GC content of the whole transcript in `[0, 1]`.
    pub gc: f64,
    /// The use case the cap / poly-A were chosen for.
    pub use_case: MrnaUseCase,
}

/// Wraps `cds` into a full five-part mRNA construct via
/// `valenx-genediting` — reference UTRs, a use-case-recommended poly-A
/// length and cap. Extracted so the headless UI tests can call it
/// directly. The CDS is used **verbatim** (no re-codon-optimisation), so
/// a LinearDesign-optimised CDS is preserved exactly.
pub fn run_construct(cds: &str, use_case: MrnaUseCase) -> Result<ConstructOutcome, String> {
    let cleaned = common::clean_sequence(cds);
    if cleaned.is_empty() {
        return Err("no optimised CDS — run the mRNA design step first".to_string());
    }
    let poly_a = recommend_poly_a(use_case);
    let cap = recommend_cap(use_case);
    let construct = MrnaConstructBuilder::new()
        .cap(cap.cap)
        .utr5(reference_utr5())
        .cds(cleaned.as_bytes())
        .utr3(reference_utr3())
        .poly_a(poly_a.length)
        .build()
        .map_err(|e| format!("[{}] {e}", e.code()))?;
    let transcript_bytes = construct.transcript();
    let gc = gc_fraction(&transcript_bytes);
    let transcript = String::from_utf8_lossy(&transcript_bytes).into_owned();
    Ok(ConstructOutcome {
        construct,
        transcript,
        gc,
        use_case,
    })
}

// ---------------------------------------------------------------------
// Background run plumbing
// ---------------------------------------------------------------------

/// Which long-running action a background thread is computing.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Job {
    /// An ensemble-defect inverse-design run.
    Inverse,
    /// A LinearDesign joint-optimisation run.
    LinearDesign,
    /// A λ-sweep Pareto-front run.
    Sweep,
}

/// The outcome a finished background thread carries back.
enum JobResult {
    /// An inverse-design run finished.
    Inverse(Result<InverseOutcome, String>),
    /// A LinearDesign run finished.
    LinearDesign(Result<MrnaOutcome, String>),
    /// A λ-sweep finished.
    Sweep(Result<Vec<ParetoPoint>, String>),
}

/// A live handle to a background design run.
///
/// Round-6 fix: carries an `Arc<AtomicBool> cancelled` flag that the
/// App's `on_exit` flips when the window closes. The current design
/// workers (`run_inverse`, `run_linear_design`, `run_pareto_sweep`)
/// don't yet poll the flag at every inner iteration — that requires
/// plumbing through `valenx-rnadesign` / `valenx-rnastruct` which is
/// out of scope for this round — but the handle is in place so any
/// future cancellable worker can honour it, and so the `spawn_with_cancel`
/// path exposed below lets tests + future codepaths cancel work
/// they HAVE plumbed.
struct BackgroundRun {
    /// Which action is in flight (for the busy label).
    job: Job,
    /// The worker thread; `None` once joined.
    thread: Option<JoinHandle<JobResult>>,
    /// Cooperative cancel flag — flipped by [`BackgroundRun::cancel`]
    /// and visible to the worker via the `Arc` clone returned from
    /// `spawn_with_cancel`.
    cancelled: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl BackgroundRun {
    /// Spawns `f` on a named background thread.
    fn spawn(job: Job, f: impl FnOnce() -> JobResult + Send + 'static) -> BackgroundRun {
        let thread = std::thread::Builder::new()
            .name("valenx-rnadesign-workbench".to_string())
            .spawn(f)
            .expect("spawn rna-designer worker thread");
        BackgroundRun {
            job,
            thread: Some(thread),
            cancelled: std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)),
        }
    }

    /// Spawns `f` with access to a shared cancellation flag. The
    /// closure receives an `Arc<AtomicBool>` it can poll at iteration
    /// boundaries; when the flag is `true` the worker should exit
    /// early with whatever result it can.
    #[cfg(test)]
    fn spawn_with_cancel(
        job: Job,
        f: impl FnOnce(std::sync::Arc<std::sync::atomic::AtomicBool>) -> JobResult + Send + 'static,
    ) -> BackgroundRun {
        let cancelled = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let inner = cancelled.clone();
        let thread = std::thread::Builder::new()
            .name("valenx-rnadesign-workbench".to_string())
            .spawn(move || f(inner))
            .expect("spawn rna-designer worker thread");
        BackgroundRun {
            job,
            thread: Some(thread),
            cancelled,
        }
    }

    /// `true` once the worker thread has finished.
    fn is_finished(&self) -> bool {
        self.thread
            .as_ref()
            .map(|t| t.is_finished())
            .unwrap_or(true)
    }

    /// Flip the cooperative cancel flag. Workers that observe the
    /// flag will exit on their next polling checkpoint. The handle
    /// itself stays in place until the next `poll_run` collects the
    /// finished thread.
    fn cancel(&self) {
        self.cancelled
            .store(true, std::sync::atomic::Ordering::SeqCst);
    }

    /// Takes the result once finished; `None` while in flight.
    fn take(&mut self) -> Option<JobResult> {
        if !self.is_finished() {
            return None;
        }
        let thread = self.thread.take()?;
        Some(thread.join().unwrap_or_else(|_| match self.job {
            Job::Inverse => JobResult::Inverse(Err("the design thread panicked".into())),
            Job::LinearDesign => JobResult::LinearDesign(Err("the design thread panicked".into())),
            Job::Sweep => JobResult::Sweep(Err("the sweep thread panicked".into())),
        }))
    }
}

// ---------------------------------------------------------------------
// Panel state
// ---------------------------------------------------------------------

/// Snapshot of every editable input across the RNA Designer's five
/// sections. Stored per Run so `Ctrl+Z` can rewind the full design
/// state.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RnaDesignerSnapshot {
    pub(crate) fold_seq: String,
    pub(crate) fold_engine: FoldEngine,
    pub(crate) inverse_target: String,
    pub(crate) inverse_gc_min: f64,
    pub(crate) inverse_gc_max: f64,
    pub(crate) inverse_max_homopolymer: usize,
    pub(crate) inverse_forbidden: String,
    pub(crate) inverse_locked: String,
    pub(crate) protein: String,
    pub(crate) host: ExpressionHost,
    pub(crate) lambda: f64,
    pub(crate) use_case: MrnaUseCase,
}

/// Form + result state for the RNA Designer workbench.
pub struct RnaDesignerPanel {
    /// The section currently shown.
    section: Section,

    // --- section 1: fold ---------------------------------------------
    /// The RNA / DNA sequence to fold.
    fold_seq: String,
    /// The folding engine.
    fold_engine: FoldEngine,
    /// The last fold outcome — drives sections 1 and 2.
    fold: Option<FoldOutcome>,
    /// A fold-section error.
    fold_error: Option<String>,

    // --- section 2: visualization ------------------------------------
    /// Which visualization is showing.
    viz: VizKind,
    /// Base-pair-probability significance threshold for the dot-plot.
    bpp_threshold: f64,

    // --- section 3: inverse design -----------------------------------
    /// The target dot-bracket to design toward.
    inverse_target: String,
    /// GC-band lower bound.
    inverse_gc_min: f64,
    /// GC-band upper bound.
    inverse_gc_max: f64,
    /// Homopolymer-run cap.
    inverse_max_homopolymer: usize,
    /// Comma / newline separated forbidden motifs.
    inverse_forbidden: String,
    /// Comma / newline separated locked positions, `idx=base` form.
    inverse_locked: String,
    /// The last inverse-design outcome.
    inverse: Option<InverseOutcome>,
    /// An inverse-section error.
    inverse_error: Option<String>,

    // --- section 4: mRNA design (LinearDesign) -----------------------
    /// The protein to encode.
    protein: String,
    /// The expression host.
    host: ExpressionHost,
    /// The structure ↔ codon trade-off weight λ.
    lambda: f64,
    /// The last single-λ LinearDesign outcome.
    mrna: Option<MrnaOutcome>,
    /// The last λ-sweep Pareto front.
    pareto: Option<Vec<ParetoPoint>>,
    /// An mRNA-section error.
    mrna_error: Option<String>,

    // --- section 5: mRNA construct -----------------------------------
    /// The therapeutic use case (drives cap / poly-A).
    use_case: MrnaUseCase,
    /// The last construct outcome.
    construct: Option<ConstructOutcome>,
    /// A construct-section error.
    construct_error: Option<String>,

    // --- background run ----------------------------------------------
    /// The in-flight background run, if any.
    run: Option<BackgroundRun>,

    /// Undo / redo over every editable input across the five sections.
    history: crate::undo::History<RnaDesignerSnapshot>,
}

/// Which structure visualization the Visualize section shows.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum VizKind {
    /// The 2-D secondary-structure diagram.
    #[default]
    Diagram,
    /// The nesting-depth mountain plot.
    Mountain,
    /// The base-pair-probability dot-plot heatmap.
    DotPlot,
}

impl VizKind {
    /// A short label.
    pub fn label(self) -> &'static str {
        match self {
            VizKind::Diagram => "2-D diagram",
            VizKind::Mountain => "Mountain plot",
            VizKind::DotPlot => "Base-pair dot-plot",
        }
    }
}

impl Default for RnaDesignerPanel {
    fn default() -> Self {
        RnaDesignerPanel {
            section: Section::Fold,
            fold_seq: "GGGAAAUCCUCUUUACCCGGAAGAGGGAAACCC".to_string(),
            fold_engine: FoldEngine::Auto,
            fold: None,
            fold_error: None,
            viz: VizKind::Diagram,
            bpp_threshold: 0.10,
            inverse_target: "(((((((....)))))))".to_string(),
            inverse_gc_min: 0.30,
            inverse_gc_max: 0.70,
            inverse_max_homopolymer: 6,
            inverse_forbidden: String::new(),
            inverse_locked: String::new(),
            inverse: None,
            inverse_error: None,
            protein: "MKVLAGDRENT".to_string(),
            host: ExpressionHost::Human,
            lambda: 1.0,
            mrna: None,
            pareto: None,
            mrna_error: None,
            use_case: MrnaUseCase::Vaccine,
            construct: None,
            construct_error: None,
            run: None,
            history: crate::undo::History::new(),
        }
    }
}

impl RnaDesignerPanel {
    /// Snapshot every editable input across the five sections.
    pub(crate) fn snapshot(&self) -> RnaDesignerSnapshot {
        RnaDesignerSnapshot {
            fold_seq: self.fold_seq.clone(),
            fold_engine: self.fold_engine,
            inverse_target: self.inverse_target.clone(),
            inverse_gc_min: self.inverse_gc_min,
            inverse_gc_max: self.inverse_gc_max,
            inverse_max_homopolymer: self.inverse_max_homopolymer,
            inverse_forbidden: self.inverse_forbidden.clone(),
            inverse_locked: self.inverse_locked.clone(),
            protein: self.protein.clone(),
            host: self.host,
            lambda: self.lambda,
            use_case: self.use_case,
        }
    }
    fn restore(&mut self, s: RnaDesignerSnapshot) {
        self.fold_seq = s.fold_seq;
        self.fold_engine = s.fold_engine;
        self.inverse_target = s.inverse_target;
        self.inverse_gc_min = s.inverse_gc_min;
        self.inverse_gc_max = s.inverse_gc_max;
        self.inverse_max_homopolymer = s.inverse_max_homopolymer;
        self.inverse_forbidden = s.inverse_forbidden;
        self.inverse_locked = s.inverse_locked;
        self.protein = s.protein;
        self.host = s.host;
        self.lambda = s.lambda;
        self.use_case = s.use_case;
    }
    /// Undo the most recent input edit.
    pub fn undo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(prev) = self.history.undo(current) {
            self.restore(prev);
            self.fold_error = None;
            self.inverse_error = None;
            self.mrna_error = None;
            self.construct_error = None;
            true
        } else {
            false
        }
    }
    /// Redo the most recently undone input edit.
    pub fn redo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(next) = self.history.redo(current) {
            self.restore(next);
            self.fold_error = None;
            self.inverse_error = None;
            self.mrna_error = None;
            self.construct_error = None;
            true
        } else {
            false
        }
    }
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
    /// Record the current snapshot before a Run. Callable from the
    /// inverse-design / mRNA spawn-thread paths so a long-running
    /// design is undo-able too.
    pub(crate) fn record_snapshot(&mut self) {
        let snap = self.snapshot();
        self.history.record(snap);
    }

    /// The section currently shown (for tests).
    #[cfg(test)]
    fn current_section(&self) -> Section {
        self.section
    }

    /// `true` while a background run is in flight.
    fn is_running(&self) -> bool {
        self.run.is_some()
    }

    /// Cooperatively cancel any in-flight background run. Used by
    /// the App's `on_exit` so closing the window flips the cancel
    /// flag on the RNA-designer worker; the worker exits on its
    /// next polling checkpoint (workers that don't yet poll the
    /// flag will still finish on their own — the flag is a no-op
    /// guarantee for them, not a hard kill).
    pub fn cancel_run(&self) {
        if let Some(run) = &self.run {
            run.cancel();
        }
    }

    /// Folds the Fold-section sequence, storing the outcome or the error.
    fn do_fold(&mut self) {
        self.fold_error = None;
        match run_fold(&self.fold_seq, self.fold_engine) {
            Ok(outcome) => self.fold = Some(outcome),
            Err(e) => {
                self.fold = None;
                self.fold_error = Some(e);
            }
        }
    }

    /// Parses the locked-position form into `(index, base)` pairs.
    fn parse_locked(&self) -> Result<Vec<(usize, char)>, String> {
        parse_locked_form(&self.inverse_locked)
    }

    /// Kicks off a background inverse-design run for the current form.
    ///
    /// The form is validated **up front** — a malformed locked-position
    /// list, an empty / unbalanced target dot-bracket or an
    /// out-of-range locked index surfaces an error here rather than
    /// spawning a doomed background thread.
    fn start_inverse(&mut self) {
        self.inverse_error = None;
        let locked = match self.parse_locked() {
            Ok(l) => l,
            Err(e) => {
                self.inverse_error = Some(e);
                return;
            }
        };
        // Validate the target dot-bracket synchronously.
        let target = self.inverse_target.clone();
        let n = match validate_target(&target) {
            Ok(n) => n,
            Err(e) => {
                self.inverse_error = Some(e);
                return;
            }
        };
        for &(idx, _) in &locked {
            if idx >= n {
                self.inverse_error =
                    Some(format!("locked position {idx} is past the {n}-nt target"));
                return;
            }
        }
        let (gc_min, gc_max) = (self.inverse_gc_min, self.inverse_gc_max);
        let max_homo = self.inverse_max_homopolymer;
        let forbidden = self.inverse_forbidden.clone();
        self.inverse = None;
        self.run = Some(BackgroundRun::spawn(Job::Inverse, move || {
            JobResult::Inverse(run_inverse(
                &target, gc_min, gc_max, max_homo, &forbidden, &locked,
            ))
        }));
    }

    /// Kicks off a background LinearDesign run for the current form.
    fn start_linear_design(&mut self) {
        self.mrna_error = None;
        let protein = self.protein.clone();
        let host = self.host;
        let lambda = self.lambda;
        self.mrna = None;
        self.run = Some(BackgroundRun::spawn(Job::LinearDesign, move || {
            JobResult::LinearDesign(run_linear_design(&protein, host, lambda))
        }));
    }

    /// Kicks off a background λ-sweep run for the current form.
    fn start_sweep(&mut self) {
        self.mrna_error = None;
        let protein = self.protein.clone();
        let host = self.host;
        self.pareto = None;
        self.run = Some(BackgroundRun::spawn(Job::Sweep, move || {
            JobResult::Sweep(run_pareto_sweep(&protein, host))
        }));
    }

    /// Assembles the mRNA construct from the last LinearDesign CDS.
    fn do_construct(&mut self) {
        self.construct_error = None;
        let Some(mrna) = &self.mrna else {
            self.construct_error =
                Some("run the mRNA design step first to get an optimised CDS".to_string());
            return;
        };
        match run_construct(&mrna.cds, self.use_case) {
            Ok(outcome) => self.construct = Some(outcome),
            Err(e) => {
                self.construct = None;
                self.construct_error = Some(e);
            }
        }
    }

    /// Polls the background run; on completion stores the result on the
    /// matching section.
    fn poll_run(&mut self) {
        let outcome = self.run.as_mut().and_then(BackgroundRun::take);
        if let Some(result) = outcome {
            self.run = None;
            match result {
                JobResult::Inverse(Ok(o)) => self.inverse = Some(o),
                JobResult::Inverse(Err(e)) => self.inverse_error = Some(e),
                JobResult::LinearDesign(Ok(o)) => self.mrna = Some(o),
                JobResult::LinearDesign(Err(e)) => self.mrna_error = Some(e),
                JobResult::Sweep(Ok(o)) => self.pareto = Some(o),
                JobResult::Sweep(Err(e)) => self.mrna_error = Some(e),
            }
        }
    }
}

// ---------------------------------------------------------------------
// Small parsing / formatting helpers
// ---------------------------------------------------------------------

/// Splits a comma / newline / whitespace separated list into uppercased
/// trimmed non-empty tokens.
fn split_list(raw: &str) -> Vec<String> {
    raw.split([',', '\n', '\r', ' ', '\t', ';'])
        .map(|t| t.trim().to_ascii_uppercase())
        .filter(|t| !t.is_empty())
        .collect()
}

/// Parses the locked-position form — `idx=base` tokens, comma / newline
/// separated, e.g. `"0=G, 13=C"` — into `(index, base)` pairs.
fn parse_locked_form(raw: &str) -> Result<Vec<(usize, char)>, String> {
    let mut out = Vec::new();
    for tok in raw
        .split([',', '\n', '\r', ';'])
        .map(str::trim)
        .filter(|t| !t.is_empty())
    {
        let Some(eq) = tok.find('=') else {
            return Err(format!("locked position `{tok}` must be `index=base`"));
        };
        let (idx_str, base_str) = tok.split_at(eq);
        let idx: usize = idx_str
            .trim()
            .parse()
            .map_err(|_| format!("`{}` is not a valid position index", idx_str.trim()))?;
        let base = match base_str[1..].trim().to_ascii_uppercase().chars().next() {
            Some(b @ ('A' | 'C' | 'G' | 'U')) => b,
            Some('T') => 'U',
            _ => return Err(format!("locked base in `{tok}` must be A / C / G / U")),
        };
        out.push((idx, base));
    }
    Ok(out)
}

/// Validates a target dot-bracket string, returning its length on
/// success — used to reject a malformed target before a background
/// inverse-design run is spawned.
fn validate_target(target_db: &str) -> Result<usize, String> {
    let db = target_db.trim();
    if db.is_empty() {
        return Err("enter a target dot-bracket structure".to_string());
    }
    let s = Structure::from_dot_bracket(db).map_err(|e| e.to_string())?;
    if s.is_empty() {
        return Err("the target structure is empty".to_string());
    }
    Ok(s.len())
}

/// GC fraction of an RNA / DNA sequence in `[0, 1]`.
fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Wraps a sequence into fixed-width lines for a mono output box.
fn wrap_seq(seq: &str, width: usize) -> String {
    if width == 0 || seq.is_empty() {
        return seq.to_string();
    }
    // Wrap every `width` *characters*. Engine output is pure ASCII
    // (A C G U / A C G T / dot-bracket), so char count == byte count in
    // practice, but stepping byte offsets by `width` (the old approach)
    // would slice mid-char and panic on any non-ASCII input. Chunking
    // by `char_indices` keeps every slice on a char boundary and is
    // identical to the old behaviour for ASCII.
    let mut out = String::with_capacity(seq.len() + seq.len() / width + 1);
    let mut chars_in_line = 0;
    for ch in seq.chars() {
        if chars_in_line == width {
            out.push('\n');
            chars_in_line = 0;
        }
        out.push(ch);
        chars_in_line += 1;
    }
    out
}

/// `true` if the dot-bracket and sequence are equal (a fold-back match).
fn structures_match(a: &str, b: &str) -> bool {
    a == b
}

// ---------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------

/// Render the RNA Designer workbench panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    // Poll the background run first so a just-finished job shows this
    // frame, and keep repainting while one is in flight.
    app.genetics.rna_designer.poll_run();
    if app.genetics.rna_designer.is_running() {
        ui.ctx().request_repaint();
    }

    draw_section_tabs(app, ui);
    ui.separator();

    let section = app.genetics.rna_designer.section;
    ui.label(egui::RichText::new(section.title()).heading());
    ui.label(egui::RichText::new(section.subtitle()).weak().small());
    ui.add_space(4.0);

    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| match section {
            Section::Fold => draw_fold_section(app, ui),
            Section::Visualize => draw_visualize_section(app, ui),
            Section::Inverse => draw_inverse_section(app, ui),
            Section::Mrna => draw_mrna_section(app, ui),
            Section::Construct => draw_construct_section(app, ui),
        });
}

/// The horizontal section selector — five tabs.
fn draw_section_tabs(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let current = app.genetics.rna_designer.section;
    let (can_undo, can_redo) = (
        app.genetics.rna_designer.can_undo(),
        app.genetics.rna_designer.can_redo(),
    );
    let mut undo_clicked = false;
    let mut redo_clicked = false;
    ui.horizontal_wrapped(|ui| {
        for s in Section::ALL {
            let resp = ui.selectable_label(s == current, s.label());
            if resp.clicked() {
                app.genetics.rna_designer.section = s;
            }
        }
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, can_undo, can_redo);
        undo_clicked = u;
        redo_clicked = r;
    });
    if undo_clicked {
        app.genetics.rna_designer.undo_edit();
    }
    if redo_clicked {
        app.genetics.rna_designer.redo_edit();
    }
}

// ---------------------------------------------------------------------
// Section 1 — Fold
// ---------------------------------------------------------------------

fn draw_fold_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let running = app.genetics.rna_designer.is_running();
    let p = &mut app.genetics.rna_designer;

    common::section(ui, "Input sequence");
    common::seq_input(
        ui,
        "rnad_fold_seq",
        "RNA (or DNA — auto-transcribed):",
        &mut p.fold_seq,
        4,
    );
    let len = common::clean_sequence(&p.fold_seq).len();
    ui.label(
        egui::RichText::new(format!("length: {len} nt"))
            .weak()
            .small(),
    );

    ui.add_space(4.0);
    common::section(ui, "Folding engine");
    ui.horizontal_wrapped(|ui| {
        for e in [FoldEngine::Auto, FoldEngine::Zuker, FoldEngine::LinearFold] {
            ui.selectable_value(&mut p.fold_engine, e, e.label());
        }
    });
    if p.fold_engine == FoldEngine::Auto {
        let resolved = p.fold_engine.resolve(len);
        ui.label(
            egui::RichText::new(format!(
                "Auto → {} (Zuker below {} nt, LinearFold at or above).",
                resolved.label(),
                FoldEngine::AUTO_LINEAR_THRESHOLD,
            ))
            .weak()
            .small(),
        );
    }

    ui.add_space(6.0);
    if common::run_button(ui, "Fold structure") {
        p.record_snapshot();
        p.do_fold();
    }

    common::error_line(ui, &p.fold_error);

    if let Some(fold) = &p.fold {
        ui.separator();
        common::section(ui, "Predicted secondary structure");
        egui::Grid::new("rnad_fold_grid")
            .num_columns(2)
            .spacing([10.0, 3.0])
            .show(ui, |ui| {
                common::kv(ui, "engine used", fold.engine.label());
                common::kv(ui, "length", format!("{} nt", fold.sequence.len()));
                common::kv(ui, "MFE", format!("{:.2} kcal/mol", fold.mfe));
                common::kv(
                    ui,
                    "ensemble ΔG",
                    format!("{:.2} kcal/mol (LinearPartition)", fold.ensemble_g),
                );
                common::kv(
                    ui,
                    "MFE-structure frequency",
                    format!(
                        "{:.1}% of the Boltzmann ensemble",
                        fold.mfe_frequency * 100.0
                    ),
                );
                common::kv(
                    ui,
                    "base pairs",
                    format!("{}", fold.structure.pairs().len()),
                );
            });
        ui.add_space(2.0);
        common::mono_output(
            ui,
            "rnad_fold_db",
            &format!(
                "sequence  {}\nstructure {}",
                fold.sequence, fold.dot_bracket,
            ),
            4,
        );
        ui.add_space(2.0);
        common::ok_line(
            ui,
            "Folded — open the Visualize section for the 2-D diagram.",
        );
        if running {
            ui.label(
                egui::RichText::new("(a background job is running)")
                    .weak()
                    .small(),
            );
        }
    }
}

// ---------------------------------------------------------------------
// Section 2 — Structure visualization
// ---------------------------------------------------------------------

fn draw_visualize_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.rna_designer;

    let Some(fold) = p.fold.clone() else {
        ui.label(
            egui::RichText::new(
                "No folded structure yet — fold a sequence in the Fold section first.",
            )
            .weak(),
        );
        return;
    };

    common::section(ui, "View");
    ui.horizontal_wrapped(|ui| {
        for v in [VizKind::Diagram, VizKind::Mountain, VizKind::DotPlot] {
            ui.selectable_value(&mut p.viz, v, v.label());
        }
    });
    ui.add_space(4.0);

    ui.label(
        egui::RichText::new(format!(
            "{} nt · {} base pair(s) · MFE {:.2} kcal/mol",
            fold.sequence.len(),
            fold.structure.pairs().len(),
            fold.mfe,
        ))
        .weak()
        .small(),
    );
    ui.add_space(4.0);

    match p.viz {
        VizKind::Diagram => draw_structure_diagram(ui, &fold),
        VizKind::Mountain => draw_mountain_plot(ui, &fold),
        VizKind::DotPlot => draw_dot_plot(p, ui, &fold),
    }
}

/// Draws the 2-D secondary-structure diagram with egui's painter, using
/// the `valenx-rnastruct` naview-class layout coordinates: the backbone
/// as a polyline, each base as a small disc, and each base pair as a
/// bond line.
fn draw_structure_diagram(ui: &mut egui::Ui, fold: &FoldOutcome) {
    common::section(ui, "2-D secondary-structure diagram");
    ui.label(
        egui::RichText::new(
            "naview-class layout — helices as ladders, loops as circles. \
             Backbone in grey, base-pair bonds in blue.",
        )
        .weak()
        .small(),
    );
    ui.add_space(4.0);

    let lay: Layout = layout(&fold.structure);
    if lay.is_empty() {
        ui.label(egui::RichText::new("(empty structure)").weak());
        return;
    }

    // Allocate a square-ish canvas.
    let avail = ui.available_width().clamp(180.0, 560.0);
    let canvas = egui::vec2(avail, (avail * 0.78).clamp(160.0, 460.0));
    let (rect, _resp) = ui.allocate_exact_size(canvas, egui::Sense::hover());
    let painter = ui.painter_at(rect);

    // Backdrop.
    painter.rect_filled(rect, 4.0, egui::Color32::from_rgb(22, 24, 30));

    // Map layout-space (the structure's bounding box) into the canvas
    // rect, preserving aspect ratio, with a margin.
    let (min_x, min_y, max_x, max_y) = lay.bounding_box();
    let span_x = (max_x - min_x).max(1.0);
    let span_y = (max_y - min_y).max(1.0);
    let margin = 18.0_f32;
    let inner_w = (rect.width() - 2.0 * margin).max(1.0);
    let inner_h = (rect.height() - 2.0 * margin).max(1.0);
    let scale = (inner_w / span_x as f32).min(inner_h / span_y as f32);
    // Centre the drawing in the rect.
    let draw_w = span_x as f32 * scale;
    let draw_h = span_y as f32 * scale;
    let off_x = rect.left() + (rect.width() - draw_w) * 0.5;
    let off_y = rect.top() + (rect.height() - draw_h) * 0.5;
    // Layout y grows upward; egui y grows downward — flip y.
    let to_screen = |x: f64, y: f64| -> egui::Pos2 {
        egui::pos2(
            off_x + ((x - min_x) as f32) * scale,
            off_y + draw_h - ((y - min_y) as f32) * scale,
        )
    };

    // 1) the backbone — a polyline through every base in 5'→3' order.
    let backbone_stroke = egui::Stroke::new(1.6, egui::Color32::from_rgb(120, 126, 140));
    for w in lay.points.windows(2) {
        painter.line_segment(
            [to_screen(w[0].x, w[0].y), to_screen(w[1].x, w[1].y)],
            backbone_stroke,
        );
    }

    // 2) the base-pair bonds.
    let bond_stroke = egui::Stroke::new(1.8, egui::Color32::from_rgb(95, 165, 240));
    for &(i, j) in &lay.pairs {
        if let (Some(pi), Some(pj)) = (lay.points.get(i), lay.points.get(j)) {
            painter.line_segment([to_screen(pi.x, pi.y), to_screen(pj.x, pj.y)], bond_stroke);
        }
    }

    // 3) the bases — a small disc each, paired bases highlighted. Label
    // the bases when the discs are large enough to be legible.
    let base_radius = (scale * 0.16).clamp(2.0, 7.0);
    let label_bases = base_radius >= 5.0 && fold.sequence.len() <= 90;
    let seq_bytes = fold.sequence.as_bytes();
    for (idx, pt) in lay.points.iter().enumerate() {
        let centre = to_screen(pt.x, pt.y);
        let paired = fold.structure.partner(idx).is_some();
        let fill = if paired {
            egui::Color32::from_rgb(95, 165, 240)
        } else {
            egui::Color32::from_rgb(150, 156, 168)
        };
        painter.circle_filled(centre, base_radius, fill);
        if label_bases {
            if let Some(&b) = seq_bytes.get(idx) {
                painter.text(
                    centre,
                    egui::Align2::CENTER_CENTER,
                    (b as char).to_string(),
                    egui::FontId::monospace(base_radius * 1.5),
                    egui::Color32::from_rgb(18, 20, 26),
                );
            }
        }
    }

    // Mark the 5' end.
    let start = to_screen(lay.points[0].x, lay.points[0].y);
    painter.text(
        start + egui::vec2(0.0, -base_radius - 7.0),
        egui::Align2::CENTER_CENTER,
        "5'",
        egui::FontId::monospace(11.0),
        egui::Color32::from_rgb(200, 205, 215),
    );

    ui.add_space(2.0);
    ui.label(
        egui::RichText::new(format!(
            "{} bases drawn · {} base-pair bond(s){}",
            lay.len(),
            lay.pairs.len(),
            if label_bases {
                ""
            } else {
                " · base letters hidden (zoom in by folding a shorter RNA)"
            },
        ))
        .weak()
        .small(),
    );
}

/// Draws the nesting-depth mountain plot — height at each position is the
/// number of base pairs enclosing it. Rendered with `egui_plot`.
fn draw_mountain_plot(ui: &mut egui::Ui, fold: &FoldOutcome) {
    common::section(ui, "Mountain plot");
    ui.label(
        egui::RichText::new(
            "Height = base pairs enclosing each position. Stems are peaks, \
             loops are plateaus.",
        )
        .weak()
        .small(),
    );
    ui.add_space(4.0);

    let heights = mountain_heights(&fold.structure);
    let points: PlotPoints = heights
        .iter()
        .enumerate()
        .map(|(i, &h)| [i as f64, h])
        .collect();
    let peak = heights.iter().cloned().fold(0.0_f64, f64::max);

    Plot::new("rnad_mountain_plot")
        .height(180.0)
        .show_axes([true, true])
        .allow_zoom(false)
        .allow_drag(false)
        .allow_scroll(false)
        .show(ui, |plot_ui| {
            plot_ui.line(
                Line::new(points)
                    .name("nesting depth")
                    .color(egui::Color32::from_rgb(95, 165, 240)),
            );
        });
    ui.label(
        egui::RichText::new(format!("peak nesting depth: {peak:.0}"))
            .weak()
            .small(),
    );
}

/// Draws the base-pair-probability dot-plot — a heatmap of the
/// LinearPartition base-pair probabilities, with the painter. Each
/// significant pair `(i, j)` is a square whose shade tracks `p(i,j)`.
fn draw_dot_plot(p: &mut RnaDesignerPanel, ui: &mut egui::Ui, fold: &FoldOutcome) {
    common::section(ui, "Base-pair-probability dot-plot");
    ui.label(
        egui::RichText::new(
            "Each cell (i, j) is shaded by the LinearPartition probability \
             that base i pairs base j across the whole Boltzmann ensemble.",
        )
        .weak()
        .small(),
    );
    ui.horizontal(|ui| {
        ui.label("show pairs with p ≥");
        ui.add(
            egui::DragValue::new(&mut p.bpp_threshold)
                .speed(0.01)
                .range(0.0..=1.0),
        );
    });
    ui.add_space(4.0);

    let n = fold.sequence.len();
    if n == 0 {
        ui.label(egui::RichText::new("(empty sequence)").weak());
        return;
    }
    let shown: Vec<&(usize, usize, f64)> = fold
        .bpp
        .iter()
        .filter(|&&(_, _, prob)| prob >= p.bpp_threshold)
        .collect();

    // A square canvas: n × n cells.
    let avail = ui.available_width().clamp(160.0, 460.0);
    let side = avail;
    let (rect, _resp) = ui.allocate_exact_size(egui::vec2(side, side), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 4.0, egui::Color32::from_rgb(22, 24, 30));

    let cell = (side / n as f32).max(0.5);
    // Cell (i, j) — i across, j down.
    let cell_rect = |i: usize, j: usize| -> egui::Rect {
        egui::Rect::from_min_size(
            egui::pos2(rect.left() + i as f32 * cell, rect.top() + j as f32 * cell),
            egui::vec2(cell, cell),
        )
    };
    // Faint diagonal so the matrix reads as a matrix.
    for k in 0..n {
        painter.rect_filled(cell_rect(k, k), 0.0, egui::Color32::from_rgb(40, 43, 52));
    }
    // The significant pairs — shade by probability. Plot both (i, j) and
    // its mirror (j, i) so the symmetric dot-plot is filled.
    for &&(i, j, prob) in &shown {
        let t = prob.clamp(0.0, 1.0) as f32;
        let color = egui::Color32::from_rgb(
            (60.0 + 30.0 * t) as u8,
            (90.0 + 100.0 * t) as u8,
            (130.0 + 110.0 * t) as u8,
        );
        if i < n && j < n {
            painter.rect_filled(cell_rect(i, j), 0.0, color);
            painter.rect_filled(cell_rect(j, i), 0.0, color);
        }
    }

    ui.add_space(2.0);
    ui.label(
        egui::RichText::new(format!(
            "{} of {} ensemble base pair(s) shown at p ≥ {:.2} · peak p = {:.3}",
            shown.len(),
            fold.bpp.len(),
            p.bpp_threshold,
            fold.max_bpp,
        ))
        .weak()
        .small(),
    );

    // The top pairs, as a compact table under the heatmap.
    if !shown.is_empty() {
        common::section(ui, "Strongest base pairs");
        let mut text = String::new();
        for &&(i, j, prob) in shown.iter().take(16) {
            let bar = "#".repeat((prob * 20.0).round() as usize);
            text.push_str(&format!(
                "{:>4}-{:<4}  p={:.3}  {}\n",
                i + 1,
                j + 1,
                prob,
                bar
            ));
        }
        common::mono_output(ui, "rnad_dotplot_table", text.trim_end(), 8);
    }
}

/// The mountain-plot heights of a structure — `heights[k]` is the number
/// of base pairs `(i, j)` with `i < k < j`. An `O(n + pairs)`
/// difference-array accumulation, mirroring `valenx-rnastruct`'s
/// `mountain_plot`.
fn mountain_heights(s: &Structure) -> Vec<f64> {
    let n = s.len();
    let mut delta = vec![0.0_f64; n + 1];
    for bp in s.pairs() {
        delta[bp.i + 1] += 1.0;
        delta[bp.j] -= 1.0;
    }
    let mut heights = vec![0.0; n];
    let mut acc = 0.0;
    for (k, h) in heights.iter_mut().enumerate() {
        acc += delta[k];
        *h = acc;
    }
    heights
}

// ---------------------------------------------------------------------
// Section 3 — Inverse design
// ---------------------------------------------------------------------

fn draw_inverse_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let running = app.genetics.rna_designer.is_running();
    let p = &mut app.genetics.rna_designer;

    common::section(ui, "Target structure");
    ui.label("Target secondary structure (dot-bracket):");
    ui.add(
        egui::TextEdit::multiline(&mut p.inverse_target)
            .id_source("rnad_inv_target")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(2)
            .hint_text("e.g. (((((((....)))))))"),
    );
    let target_len = p
        .inverse_target
        .trim()
        .chars()
        .filter(|c| !c.is_whitespace())
        .count();
    ui.label(
        egui::RichText::new(format!(
            "length: {target_len}  ·  the designed sequence will be {target_len} nt",
        ))
        .weak()
        .small(),
    );

    ui.add_space(4.0);
    common::section(ui, "Design constraints");
    ui.horizontal(|ui| {
        ui.label("GC band");
        ui.add(
            egui::DragValue::new(&mut p.inverse_gc_min)
                .speed(0.01)
                .range(0.0..=1.0)
                .custom_formatter(|v, _| format!("{:.0}%", v * 100.0)),
        );
        ui.label("to");
        ui.add(
            egui::DragValue::new(&mut p.inverse_gc_max)
                .speed(0.01)
                .range(0.0..=1.0)
                .custom_formatter(|v, _| format!("{:.0}%", v * 100.0)),
        );
    });
    ui.horizontal(|ui| {
        ui.label("max homopolymer run");
        ui.add(
            egui::DragValue::new(&mut p.inverse_max_homopolymer)
                .speed(1.0)
                .range(1..=30),
        );
        ui.label("nt");
    });
    ui.label("Forbidden motifs (comma / newline separated):");
    ui.add(
        egui::TextEdit::singleline(&mut p.inverse_forbidden)
            .id_source("rnad_inv_forbidden")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .hint_text("e.g. GGGG, AAUAAA"),
    );
    ui.label("Locked positions (index=base, comma separated):");
    ui.add(
        egui::TextEdit::singleline(&mut p.inverse_locked)
            .id_source("rnad_inv_locked")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .hint_text("e.g. 0=G, 17=C"),
    );

    ui.add_space(6.0);
    if running {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label("Designing — ensemble-defect inverse folding…");
        });
    } else if common::run_button(ui, "Run inverse design") {
        p.record_snapshot();
        p.start_inverse();
    }

    common::error_line(ui, &p.inverse_error);

    if let Some(inv) = &p.inverse {
        ui.separator();
        common::section(ui, "Designed sequence");
        let matches = structures_match(&inv.fold_back_db, &inv.target_db);
        egui::Grid::new("rnad_inv_grid")
            .num_columns(2)
            .spacing([10.0, 3.0])
            .show(ui, |ui| {
                common::kv(ui, "length", format!("{} nt", inv.sequence.len()));
                common::kv(
                    ui,
                    "ensemble defect",
                    format!(
                        "{:.3} ({:.4} normalised)",
                        inv.ensemble_defect, inv.normalized_defect
                    ),
                );
                common::kv(
                    ui,
                    "defect target",
                    if inv.solved { "reached" } else { "not reached" },
                );
                common::kv(
                    ui,
                    "MFE-fold distance",
                    format!("{} base pair(s) from the target", inv.mfe_distance),
                );
                common::kv(
                    ui,
                    "mutations",
                    format!("{} accepted / {} attempted", inv.accepted, inv.total),
                );
                if inv.constrained {
                    common::kv(
                        ui,
                        "constraints",
                        format!(
                            "{} locked position(s), GC / motif band honoured",
                            inv.locked_count
                        ),
                    );
                }
            });
        ui.add_space(2.0);
        common::mono_output(
            ui,
            "rnad_inv_seq",
            &format!(
                "designed   {}\nfold-back  {}\ntarget     {}",
                inv.sequence, inv.fold_back_db, inv.target_db,
            ),
            5,
        );
        ui.add_space(2.0);
        if matches {
            common::ok_line(ui, "✔ The design's MFE structure is exactly the target.");
        } else {
            ui.colored_label(
                egui::Color32::from_rgb(235, 195, 110),
                format!(
                    "▲ The MFE fold differs from the target by {} pair(s) — the \
                     ensemble defect is the principled objective and stays low.",
                    inv.mfe_distance,
                ),
            );
        }
    }
}

// ---------------------------------------------------------------------
// Section 4 — mRNA design (LinearDesign)
// ---------------------------------------------------------------------

fn draw_mrna_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let running = app.genetics.rna_designer.is_running();
    let p = &mut app.genetics.rna_designer;

    common::section(ui, "Protein to encode");
    common::seq_input(
        ui,
        "rnad_mrna_protein",
        "Protein sequence (one-letter amino-acid codes):",
        &mut p.protein,
        3,
    );
    ui.label(
        egui::RichText::new(format!(
            "{} residue(s) — a leading M and a trailing stop are added if absent",
            common::clean_sequence(&p.protein).len(),
        ))
        .weak()
        .small(),
    );

    ui.add_space(4.0);
    common::section(ui, "Host & trade-off");
    ui.horizontal(|ui| {
        ui.label("Host organism:");
        egui::ComboBox::from_id_source("rnad_mrna_host")
            .selected_text(host_label(p.host))
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut p.host,
                    ExpressionHost::Human,
                    host_label(ExpressionHost::Human),
                );
                ui.selectable_value(
                    &mut p.host,
                    ExpressionHost::EColi,
                    host_label(ExpressionHost::EColi),
                );
            });
    });
    ui.horizontal(|ui| {
        ui.label("λ (structure ↔ codon trade-off):");
        ui.add(
            egui::Slider::new(&mut p.lambda, 0.0..=16.0)
                .step_by(0.1)
                .fixed_decimals(1),
        );
    });
    ui.label(
        egui::RichText::new(
            "λ = 0 → the most stable mRNA (pure minimum free energy); \
             large λ → the most codon-optimal mRNA (pure CAI).",
        )
        .weak()
        .small(),
    );

    ui.add_space(6.0);
    if running {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label("Optimising — the LinearDesign lattice DP…");
        });
    } else {
        ui.horizontal(|ui| {
            if ui
                .add_sized([0.0, 24.0], egui::Button::new("Run LinearDesign"))
                .clicked()
            {
                p.start_linear_design();
            }
            if ui
                .add_sized([0.0, 24.0], egui::Button::new("Sweep λ → Pareto front"))
                .clicked()
            {
                p.start_sweep();
            }
        });
    }

    common::error_line(ui, &p.mrna_error);

    // Single-λ result.
    if let Some(m) = &p.mrna {
        ui.separator();
        common::section(ui, "Optimised coding sequence");
        egui::Grid::new("rnad_mrna_grid")
            .num_columns(2)
            .spacing([10.0, 3.0])
            .show(ui, |ui| {
                common::kv(ui, "host", host_label(m.host));
                common::kv(ui, "λ used", format!("{:.2}", m.lambda));
                common::kv(ui, "CDS length", format!("{} nt", m.cds.len()));
                common::kv(ui, "MFE", format!("{:.2} kcal/mol", m.mfe));
                common::kv(ui, "CAI", format!("{:.3}", m.cai));
                common::kv(
                    ui,
                    "combined objective",
                    format!("{:.3}  (MFE + λ·codon-penalty)", m.objective),
                );
                common::kv(
                    ui,
                    "lattice DP",
                    if m.exact {
                        "exact (no beam pruning)"
                    } else {
                        "beam-pruned (near-optimal)"
                    },
                );
            });
        ui.add_space(2.0);
        common::mono_output(
            ui,
            "rnad_mrna_cds",
            &format!(
                "{}\n\nstructure:\n{}",
                wrap_seq(&m.cds, 60),
                wrap_seq(&m.dot_bracket, 60),
            ),
            8,
        );
        ui.add_space(2.0);
        common::ok_line(
            ui,
            "CDS ready — the mRNA construct section wraps it into a full transcript.",
        );
    }

    // λ-sweep Pareto front.
    if let Some(front) = &p.pareto {
        ui.separator();
        common::section(ui, "Stability / efficiency Pareto front");
        ui.label(
            egui::RichText::new(
                "Each point is a CDS designed at one λ. As λ rises CAI climbs \
                 (more codon-optimal) and MFE rises (less stable) — the classic \
                 trade-off.",
            )
            .weak()
            .small(),
        );
        ui.add_space(4.0);
        draw_pareto_plot(ui, front);
        common::section(ui, "Sweep points");
        let mut text = String::from("    λ        MFE (kcal/mol)    CAI\n");
        for pt in front {
            text.push_str(&format!(
                "{:>6.2}   {:>14.2}   {:>8.3}\n",
                pt.lambda, pt.mfe, pt.cai,
            ));
        }
        common::mono_output(ui, "rnad_pareto_table", text.trim_end(), 9);
    }
}

/// Draws the λ-sweep Pareto front — CAI (x) against MFE (y) — with
/// `egui_plot`: a connected line plus a marker per λ.
fn draw_pareto_plot(ui: &mut egui::Ui, front: &[ParetoPoint]) {
    if front.is_empty() {
        ui.label(egui::RichText::new("(no sweep points)").weak());
        return;
    }
    let line_pts: PlotPoints = front.iter().map(|p| [p.cai, p.mfe]).collect();
    let marker_pts: PlotPoints = front.iter().map(|p| [p.cai, p.mfe]).collect();

    Plot::new("rnad_pareto_plot")
        .height(200.0)
        .show_axes([true, true])
        .x_axis_label("CAI (codon adaptation index) →")
        .y_axis_label("MFE (kcal/mol)")
        .allow_zoom(false)
        .allow_drag(false)
        .allow_scroll(false)
        .show(ui, |plot_ui| {
            plot_ui.line(
                Line::new(line_pts)
                    .name("Pareto front")
                    .color(egui::Color32::from_rgb(120, 170, 240)),
            );
            plot_ui.points(
                Points::new(marker_pts)
                    .name("λ point")
                    .radius(4.0)
                    .shape(MarkerShape::Circle)
                    .color(egui::Color32::from_rgb(240, 190, 110)),
            );
        });
}

// ---------------------------------------------------------------------
// Section 5 — mRNA construct
// ---------------------------------------------------------------------

fn draw_construct_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.rna_designer;

    let has_cds = p.mrna.is_some();
    common::section(ui, "Optimised CDS");
    match &p.mrna {
        Some(m) => {
            ui.label(
                egui::RichText::new(format!(
                    "Using the LinearDesign CDS — {} nt, MFE {:.2} kcal/mol, CAI {:.3}.",
                    m.cds.len(),
                    m.mfe,
                    m.cai,
                ))
                .weak()
                .small(),
            );
        }
        None => {
            ui.label(
                egui::RichText::new(
                    "No optimised CDS yet — run the mRNA design (LinearDesign) section first.",
                )
                .weak(),
            );
        }
    }

    ui.add_space(4.0);
    common::section(ui, "Construct settings");
    ui.horizontal(|ui| {
        ui.label("Use case:");
        egui::ComboBox::from_id_source("rnad_construct_usecase")
            .selected_text(p.use_case.name())
            .show_ui(ui, |ui| {
                for uc in [
                    MrnaUseCase::Vaccine,
                    MrnaUseCase::ProteinReplacement,
                    MrnaUseCase::TransientEditing,
                ] {
                    ui.selectable_value(&mut p.use_case, uc, uc.name());
                }
            });
    });
    ui.label(
        egui::RichText::new(
            "The use case picks the poly-A length and cap chemistry; reference \
             5′/3′ UTRs (Kozak context, stability elements) wrap the CDS.",
        )
        .weak()
        .small(),
    );

    ui.add_space(6.0);
    let assemble = ui
        .add_enabled(
            has_cds,
            egui::Button::new("Assemble mRNA construct")
                .min_size(egui::vec2(ui.available_width(), 24.0)),
        )
        .clicked();
    if assemble {
        p.do_construct();
    }

    common::error_line(ui, &p.construct_error);

    if let Some(c) = &p.construct {
        let con = &c.construct;
        ui.separator();
        common::section(ui, "Assembled construct");
        egui::Grid::new("rnad_construct_grid")
            .num_columns(2)
            .spacing([10.0, 3.0])
            .show(ui, |ui| {
                common::kv(ui, "cap", con.cap.name());
                common::kv(ui, "5′UTR", format!("{} nt", con.utr5.len()));
                common::kv(
                    ui,
                    "CDS",
                    format!("{} nt ({} codons)", con.cds.len(), con.codon_count()),
                );
                common::kv(ui, "3′UTR", format!("{} nt", con.utr3.len()));
                common::kv(ui, "poly-A tail", format!("{} nt", con.poly_a_len));
                common::kv(
                    ui,
                    "total transcript",
                    format!("{} nt (excl. cap)", con.len()),
                );
                common::kv(ui, "GC content", format!("{:.1}%", c.gc * 100.0));
                common::kv(ui, "CDS start", format!("position {}", con.cds_start()));
            });
        ui.add_space(2.0);
        common::section(ui, "Construct map (5′ → 3′)");
        let map = format!(
            "[cap: {}]\n5'UTR   1..{}\nCDS     {}..{}\n3'UTR   {}..{}\npoly-A  {}..{}",
            con.cap.name(),
            con.utr5.len(),
            con.utr5.len() + 1,
            con.utr5.len() + con.cds.len(),
            con.utr5.len() + con.cds.len() + 1,
            con.utr5.len() + con.cds.len() + con.utr3.len(),
            con.utr5.len() + con.cds.len() + con.utr3.len() + 1,
            con.len(),
        );
        common::mono_output(ui, "rnad_construct_map", &map, 6);
        ui.add_space(2.0);
        common::section(ui, "Full transcript (5′UTR + CDS + 3′UTR + poly-A)");
        common::mono_output(
            ui,
            "rnad_construct_transcript",
            &wrap_seq(&c.transcript, 60),
            8,
        );
        ui.add_space(2.0);
        common::ok_line(
            ui,
            "✔ A validated five-part mRNA construct — the synthesis-ready transcript.",
        );
    }
}

/// A human-readable expression-host label.
fn host_label(h: ExpressionHost) -> &'static str {
    match h {
        ExpressionHost::Human => "Homo sapiens (human)",
        ExpressionHost::EColi => "Escherichia coli",
    }
}

// ---------------------------------------------------------------------
// Tests — pure (non-UI) logic
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sections_have_text() {
        assert_eq!(Section::ALL.len(), 5);
        for s in Section::ALL {
            assert!(!s.label().is_empty());
            assert!(!s.title().is_empty());
            assert!(!s.subtitle().is_empty());
        }
    }

    #[test]
    fn fold_engine_auto_resolves_by_length() {
        // Short → Zuker, long → LinearFold.
        assert_eq!(FoldEngine::Auto.resolve(50), FoldEngine::Zuker);
        assert_eq!(FoldEngine::Auto.resolve(5000), FoldEngine::LinearFold);
        assert_eq!(
            FoldEngine::Auto.resolve(FoldEngine::AUTO_LINEAR_THRESHOLD),
            FoldEngine::LinearFold,
        );
        // An explicit engine is never overridden.
        assert_eq!(FoldEngine::Zuker.resolve(99_999), FoldEngine::Zuker);
        assert_eq!(FoldEngine::LinearFold.resolve(1), FoldEngine::LinearFold);
    }

    #[test]
    fn default_panel_starts_at_fold() {
        let p = RnaDesignerPanel::default();
        assert_eq!(p.current_section(), Section::Fold);
        assert!(p.fold.is_none());
        assert!(!p.is_running());
    }

    #[test]
    fn background_run_cancel_signals_a_polling_worker() {
        // Round-6 RED→GREEN: BackgroundRun now carries an
        // Arc<AtomicBool> cancel flag. A worker that polls the flag
        // exits early when `cancel()` is called. We exercise this
        // via the cfg(test)-only `spawn_with_cancel` helper.
        use std::sync::atomic::Ordering;
        let run = BackgroundRun::spawn_with_cancel(Job::Inverse, |cancelled| {
            // Polled-busy loop: stop the moment cancel is observed.
            // Bounded so the test doesn't hang if cancel never fires.
            for _ in 0..10_000 {
                if cancelled.load(Ordering::SeqCst) {
                    return JobResult::Inverse(Err("cancelled by test".into()));
                }
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
            JobResult::Inverse(Err("test polled-loop expired".into()))
        });
        // Cancel a few ms in to confirm the flag actually flips
        // mid-flight.
        std::thread::sleep(std::time::Duration::from_millis(20));
        run.cancel();
        // Give the worker a bounded window to exit.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
        while !run.is_finished() && std::time::Instant::now() < deadline {
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(run.is_finished(), "worker did not honour cancel within 2s");
    }

    #[test]
    fn panel_cancel_run_is_safe_with_no_active_run() {
        // Defensive: cancelling when no run is in flight must be a
        // no-op so on_exit can call cancel_run unconditionally.
        let p = RnaDesignerPanel::default();
        assert!(!p.is_running());
        p.cancel_run(); // must not panic
    }

    #[test]
    fn split_list_parses_mixed_separators() {
        assert_eq!(
            split_list("gggg, aauaaa\nccc"),
            vec!["GGGG", "AAUAAA", "CCC"],
        );
        assert!(split_list("   ").is_empty());
    }

    #[test]
    fn parse_locked_form_reads_index_base_pairs() {
        let got = parse_locked_form("0=G, 13=c, 5=T").unwrap();
        // T folds to U; case is normalised.
        assert_eq!(got, vec![(0, 'G'), (13, 'C'), (5, 'U')]);
        assert!(parse_locked_form("").unwrap().is_empty());
    }

    #[test]
    fn parse_locked_form_rejects_bad_input() {
        // Missing '='.
        assert!(parse_locked_form("0G").is_err());
        // Bad base.
        assert!(parse_locked_form("0=Z").is_err());
        // Bad index.
        assert!(parse_locked_form("x=G").is_err());
    }

    #[test]
    fn gc_fraction_is_correct() {
        assert!((gc_fraction(b"GGCC") - 1.0).abs() < 1e-9);
        assert!((gc_fraction(b"AAUU") - 0.0).abs() < 1e-9);
        assert!((gc_fraction(b"ACGU") - 0.5).abs() < 1e-9);
        assert_eq!(gc_fraction(b""), 0.0);
    }

    #[test]
    fn wrap_seq_breaks_into_lines() {
        assert_eq!(wrap_seq("ACGUACGUAC", 4), "ACGU\nACGU\nAC");
        assert_eq!(wrap_seq("ACGU", 80), "ACGU");
        assert_eq!(wrap_seq("", 10), "");
    }

    #[test]
    fn wrap_seq_multibyte_does_not_panic() {
        // R32 L2: wrap_seq sliced `&seq[i..end]` on BYTE offsets stepped
        // by `width`, relying on an unguarded "pure ASCII" invariant.
        // Engine output is ACGU/dot-bracket (ASCII) today, but a
        // multibyte char would land `end` mid-char → panic. Now
        // char-aware: wraps every `width` *characters*, never panics.
        assert_eq!(
            wrap_seq("\u{20AC}\u{20AC}\u{20AC}", 1),
            "\u{20AC}\n\u{20AC}\n\u{20AC}"
        );
        // ASCII behaviour is unchanged (chars == bytes).
        assert_eq!(wrap_seq("ACGUAC", 2), "AC\nGU\nAC");
    }

    #[test]
    fn mountain_heights_match_nesting() {
        // A hairpin's loop sits at the peak depth.
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let h = mountain_heights(&s);
        assert_eq!(h.len(), 10);
        assert_eq!(h[5], 3.0, "the loop should be at nesting depth 3");
        assert_eq!(h[0], 0.0, "the 5' end is at ground level");
        assert_eq!(*h.last().unwrap(), 0.0, "the 3' end is at ground level");
        // An unpaired structure is flat.
        let flat = mountain_heights(&Structure::empty(6));
        assert!(flat.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn structures_match_compares_dot_brackets() {
        assert!(structures_match("((..))", "((..))"));
        assert!(!structures_match("((..))", "(....)"));
    }

    #[test]
    fn validate_target_accepts_and_rejects() {
        // A balanced dot-bracket validates and reports its length.
        assert_eq!(validate_target("(((....)))").unwrap(), 10);
        // Whitespace is trimmed.
        assert_eq!(validate_target("  ((..))  ").unwrap(), 6);
        // An empty target is rejected.
        assert!(validate_target("   ").is_err());
        // An unbalanced dot-bracket is rejected.
        assert!(validate_target("(((....").is_err());
    }
}

/// Headless egui UI-logic tests for the RNA Designer workbench panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    /// Draws the panel once into a headless egui context.
    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    /// A fresh app with the RNA Designer panel active.
    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::RnaDesigner;
        app
    }

    /// Drives whatever background run is in flight to completion — kicks
    /// the poll loop until the panel is idle again.
    fn drain_background(app: &mut ValenxApp) {
        for _ in 0..12_000 {
            app.genetics.rna_designer.poll_run();
            if !app.genetics.rna_designer.is_running() {
                return;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        panic!("a background run did not finish in time");
    }

    // --- draw-every-section, every state -----------------------------

    #[test]
    fn draws_every_section_fresh_without_panic() {
        // Every section renders from a fresh panel (no results yet).
        for section in Section::ALL {
            let mut app = app_with_panel();
            app.genetics.rna_designer.section = section;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_every_visualization_without_panic() {
        // The Visualize section's three views all render once a fold
        // exists.
        let mut app = app_with_panel();
        app.genetics.rna_designer.do_fold();
        assert!(app.genetics.rna_designer.fold.is_some());
        for viz in [VizKind::Diagram, VizKind::Mountain, VizKind::DotPlot] {
            app.genetics.rna_designer.section = Section::Visualize;
            app.genetics.rna_designer.viz = viz;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_error_states_without_panic() {
        // A per-section error line on every section must still render.
        let mut app = app_with_panel();
        app.genetics.rna_designer.fold_error = Some("enter an RNA sequence".into());
        app.genetics.rna_designer.inverse_error = Some("bad target".into());
        app.genetics.rna_designer.mrna_error = Some("empty protein".into());
        app.genetics.rna_designer.construct_error = Some("no CDS".into());
        for section in Section::ALL {
            app.genetics.rna_designer.section = section;
            draw_headless(&mut app);
        }
    }

    // --- section 1: Fold ---------------------------------------------

    #[test]
    fn fold_run_produces_a_structure() {
        // The default sequence folds via the real valenx-rnastruct API.
        let mut app = app_with_panel();
        app.genetics.rna_designer.do_fold();
        let p = &app.genetics.rna_designer;
        assert!(p.fold_error.is_none(), "fold errored: {:?}", p.fold_error);
        let fold = p.fold.as_ref().expect("fold produced no outcome");
        // The dot-bracket is the sequence's length.
        assert_eq!(fold.dot_bracket.len(), fold.sequence.len());
        // The MFE structure carries the bulk of a real ensemble.
        assert!((0.0..=1.0).contains(&fold.mfe_frequency));
        // The structure is pseudoknot-free (MFE folders are).
        assert!(fold.structure.is_nested());
    }

    #[test]
    fn fold_each_engine_runs() {
        // Every explicit engine folds the default sequence.
        for engine in [FoldEngine::Zuker, FoldEngine::LinearFold] {
            let outcome =
                run_fold("GGGAAAUCCUCUUUACCCGGAAGAGGGAAACCC", engine).expect("fold should succeed");
            assert_eq!(outcome.engine, engine);
            assert_eq!(outcome.dot_bracket.len(), outcome.sequence.len());
        }
        // Auto resolves and folds.
        let auto = run_fold("GGGGAAAACCCC", FoldEngine::Auto).unwrap();
        assert_eq!(auto.engine, FoldEngine::Zuker);
    }

    #[test]
    fn fold_surfaces_error_on_empty_input() {
        let mut app = app_with_panel();
        app.genetics.rna_designer.fold_seq = "   ".to_string();
        app.genetics.rna_designer.do_fold();
        assert!(app.genetics.rna_designer.fold_error.is_some());
        assert!(app.genetics.rna_designer.fold.is_none());
    }

    #[test]
    fn fold_surfaces_error_on_malformed_input() {
        // Non-nucleotide letters are malformed RNA.
        let err = run_fold("ZZZXXXQQQ", FoldEngine::Zuker);
        assert!(err.is_err(), "malformed RNA should error");
    }

    // --- section 3: Inverse design -----------------------------------

    #[test]
    fn inverse_design_runs_and_folds_to_target() {
        // The real ensemble-defect designer on a background thread.
        let mut app = app_with_panel();
        app.genetics.rna_designer.section = Section::Inverse;
        app.genetics.rna_designer.start_inverse();
        drain_background(&mut app);
        let p = &app.genetics.rna_designer;
        assert!(
            p.inverse_error.is_none(),
            "inverse errored: {:?}",
            p.inverse_error
        );
        let inv = p.inverse.as_ref().expect("inverse produced no outcome");
        // The designed sequence is the target's length.
        assert_eq!(inv.sequence.len(), inv.target_db.len());
        // A clean hairpin design should reach a low normalised defect.
        assert!(
            inv.normalized_defect < 0.3,
            "normalised defect too high: {}",
            inv.normalized_defect,
        );
        // Post-run the Inverse section renders its result.
        draw_headless(&mut app);
    }

    #[test]
    fn inverse_design_honours_locked_positions() {
        // Lock the outer pair of a 14-nt hairpin and confirm it holds.
        let outcome = run_inverse("(((((....)))))", 0.0, 1.0, 6, "", &[(0, 'G'), (13, 'C')])
            .expect("constrained inverse design should succeed");
        assert!(outcome.sequence.starts_with('G'), "locked pos 0 not held");
        assert!(outcome.sequence.ends_with('C'), "locked pos 13 not held");
        assert_eq!(outcome.locked_count, 2);
        assert!(outcome.constrained);
    }

    #[test]
    fn inverse_design_surfaces_error_on_bad_target() {
        // An empty target — start_inverse must error, not spawn a run.
        let mut app = app_with_panel();
        app.genetics.rna_designer.inverse_target = "   ".to_string();
        app.genetics.rna_designer.start_inverse();
        assert!(app.genetics.rna_designer.inverse_error.is_some());
        assert!(!app.genetics.rna_designer.is_running());
        // A malformed locked-position form is also caught up front.
        let mut app = app_with_panel();
        app.genetics.rna_designer.inverse_locked = "not-a-pair".to_string();
        app.genetics.rna_designer.start_inverse();
        assert!(app.genetics.rna_designer.inverse_error.is_some());
        assert!(!app.genetics.rna_designer.is_running());
    }

    #[test]
    fn inverse_design_rejects_unbalanced_dot_bracket() {
        // An unbalanced dot-bracket is not a valid structure.
        let err = run_inverse("(((....", 0.0, 1.0, 6, "", &[]);
        assert!(err.is_err(), "unbalanced dot-bracket should error");
    }

    // --- section 4: mRNA design (LinearDesign) -----------------------

    #[test]
    fn linear_design_runs_and_produces_a_cds() {
        // The real LinearDesign optimiser on a background thread.
        let mut app = app_with_panel();
        app.genetics.rna_designer.section = Section::Mrna;
        app.genetics.rna_designer.start_linear_design();
        drain_background(&mut app);
        let p = &app.genetics.rna_designer;
        assert!(
            p.mrna_error.is_none(),
            "linear-design errored: {:?}",
            p.mrna_error
        );
        let m = p.mrna.as_ref().expect("linear-design produced no CDS");
        // The CDS is a valid AUG…stop coding sequence.
        assert!(m.cds.len() % 3 == 0, "CDS length not a multiple of 3");
        assert!(m.cds.starts_with("AUG"), "CDS does not start with AUG");
        assert!((0.0..=1.0).contains(&m.cai), "CAI out of range: {}", m.cai);
        assert_eq!(m.dot_bracket.len(), m.cds.len());
        draw_headless(&mut app);
    }

    #[test]
    fn linear_design_lambda_extremes_differ_in_cai() {
        // λ = 0 maximises stability; a large λ maximises CAI. The
        // CAI-weighted design must have a CAI at least as high.
        let stable = run_linear_design("MKVLAGD", ExpressionHost::Human, 0.0)
            .expect("λ=0 design should succeed");
        let optimal = run_linear_design("MKVLAGD", ExpressionHost::Human, 16.0)
            .expect("λ=16 design should succeed");
        assert!(
            optimal.cai >= stable.cai - 1e-6,
            "a higher λ should not lower the CAI ({} vs {})",
            optimal.cai,
            stable.cai,
        );
    }

    #[test]
    fn pareto_sweep_runs_and_is_monotone() {
        // The λ-sweep on a background thread; the front should be
        // CAI-monotone (parametric-optimisation theory).
        let mut app = app_with_panel();
        app.genetics.rna_designer.section = Section::Mrna;
        app.genetics.rna_designer.start_sweep();
        drain_background(&mut app);
        let p = &app.genetics.rna_designer;
        assert!(p.mrna_error.is_none(), "sweep errored: {:?}", p.mrna_error);
        let front = p.pareto.as_ref().expect("sweep produced no front");
        assert!(front.len() >= 2, "the sweep should have multiple points");
        // CAI is non-decreasing in λ.
        for w in front.windows(2) {
            assert!(
                w[1].cai >= w[0].cai - 1e-6,
                "CAI not monotone across the λ sweep: {} then {}",
                w[0].cai,
                w[1].cai,
            );
        }
        draw_headless(&mut app);
    }

    #[test]
    fn linear_design_surfaces_error_on_empty_protein() {
        let mut app = app_with_panel();
        app.genetics.rna_designer.protein = "   ".to_string();
        app.genetics.rna_designer.start_linear_design();
        drain_background(&mut app);
        assert!(app.genetics.rna_designer.mrna_error.is_some());
        assert!(app.genetics.rna_designer.mrna.is_none());
    }

    #[test]
    fn linear_design_rejects_bad_residue() {
        // 'Z' is not a standard amino-acid code.
        let err = run_linear_design("MKZLA", ExpressionHost::Human, 1.0);
        assert!(err.is_err(), "a bad residue should error");
    }

    // --- section 5: mRNA construct -----------------------------------

    #[test]
    fn construct_assembly_wraps_a_linear_design_cds() {
        // Run LinearDesign, then assemble the construct from its CDS.
        let mut app = app_with_panel();
        app.genetics.rna_designer.start_linear_design();
        drain_background(&mut app);
        assert!(app.genetics.rna_designer.mrna.is_some());
        app.genetics.rna_designer.section = Section::Construct;
        app.genetics.rna_designer.do_construct();
        let p = &app.genetics.rna_designer;
        assert!(
            p.construct_error.is_none(),
            "construct errored: {:?}",
            p.construct_error,
        );
        let c = p.construct.as_ref().expect("construct produced no outcome");
        // The construct is a real five-part transcript.
        assert!(c.construct.codon_count() >= 1);
        assert!(!c.construct.utr5.is_empty(), "no 5'UTR");
        assert!(!c.construct.utr3.is_empty(), "no 3'UTR");
        assert!(c.construct.poly_a_len > 0, "no poly-A tail");
        // The transcript concatenates the parts.
        assert_eq!(c.transcript.len(), c.construct.len());
        draw_headless(&mut app);
    }

    #[test]
    fn construct_run_helper_builds_a_valid_construct() {
        // run_construct on a hand-built valid CDS.
        let outcome = run_construct("AUGGCCCUGUAA", MrnaUseCase::Vaccine)
            .expect("a valid CDS should assemble");
        assert_eq!(outcome.construct.codon_count(), 4);
        // The CDS is preserved verbatim inside the transcript.
        assert!(String::from_utf8_lossy(&outcome.construct.cds).contains("AUGGCCCUGUAA"),);
    }

    #[test]
    fn construct_surfaces_error_without_a_cds() {
        // No LinearDesign result → do_construct must error gracefully.
        let mut app = app_with_panel();
        app.genetics.rna_designer.section = Section::Construct;
        app.genetics.rna_designer.do_construct();
        assert!(app.genetics.rna_designer.construct_error.is_some());
        assert!(app.genetics.rna_designer.construct.is_none());
    }

    #[test]
    fn construct_run_helper_rejects_an_invalid_cds() {
        // A CDS with no start codon is rejected by the builder.
        let err = run_construct("GCCGCCUAA", MrnaUseCase::Vaccine);
        assert!(err.is_err(), "a CDS without AUG should error");
    }

    // --- post-run / populated states ---------------------------------

    #[test]
    fn draws_every_section_populated_without_panic() {
        // Run all four actions, then draw every section with its result
        // populated — the post-run render path.
        let mut app = app_with_panel();
        app.genetics.rna_designer.do_fold();
        app.genetics.rna_designer.start_inverse();
        drain_background(&mut app);
        app.genetics.rna_designer.start_linear_design();
        drain_background(&mut app);
        app.genetics.rna_designer.start_sweep();
        drain_background(&mut app);
        app.genetics.rna_designer.do_construct();
        assert!(app.genetics.rna_designer.fold.is_some());
        assert!(app.genetics.rna_designer.inverse.is_some());
        assert!(app.genetics.rna_designer.mrna.is_some());
        assert!(app.genetics.rna_designer.pareto.is_some());
        assert!(app.genetics.rna_designer.construct.is_some());
        for section in Section::ALL {
            app.genetics.rna_designer.section = section;
            draw_headless(&mut app);
        }
    }
}

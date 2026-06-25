//! The right-side **Protein-Interaction (PPI / interactome) Workbench** panel
//! — a native front-end over the in-house [`valenx_ppi`] crate (Valenx's
//! sequence-first coevolution PPI / interactome engine).
//!
//! A protein-protein interaction network is a graph whose **nodes are
//! proteins** and whose **edges are (scored) interactions**. valenx-ppi infers
//! those edges *from sequence*: given a **paired** multiple-sequence alignment
//! of two chains' orthologues, it scores how strongly the two chains coevolve
//! (APC-corrected mutual information over inter-chain alignment columns), folds
//! that into one comparable `[0, 1]` [`valenx_ppi::PpiScore`], and an
//! all-vs-all [`valenx_ppi::interactome_screen`] ranks every host × pathogen
//! pair into a [`valenx_ppi::RankedInteractions`] table. A real screen pairs
//! deep, well-curated orthologue alignments; building those needs external
//! sequence databases that cannot exist in the headless CI environment, so this
//! workbench drives the **real** [`valenx_ppi::interactome_screen`] /
//! [`valenx_ppi::predict_contacts`] over a fully-native, fully-transparent
//! **demo interactome**: a small named host × pathogen panel whose paired
//! orthologue rows are generated *deterministically* (no RNG) so a known set of
//! protein pairs genuinely coevolve (the planted "true" interactions) while the
//! rest are conserved noise.
//!
//! ```text
//!   hosts:     GUARD  KINASE  RECEPTOR        (the central GUARD is the hub)
//!   pathogens: EFF-A  EFF-B   EFF-C   EFF-D
//!
//!   GUARD coevolves with every pathogen effector  -> GUARD is the hub
//!   KINASE / RECEPTOR coevolve with one effector each (or none)
//! ```
//!
//! From valenx-ppi's scored pairs the workbench assembles the **interaction
//! network**: a node per protein, and an edge for every host × pathogen pair
//! whose real [`PpiScore::value`](valenx_ppi::PpiScore) clears the user's
//! highlight threshold. Over that real network it then computes two standard,
//! analytically-checkable graph metrics **in the workbench** (the underlying
//! numbers — every edge weight — are valenx-ppi's, never invented):
//!
//! * **Degree centrality** — the fraction of other nodes each protein is
//!   linked to (`deg / (N - 1)`). On a hub-and-spoke interactome the hub has
//!   the maximum degree centrality; the PIN asserts exactly that.
//! * **Shortest path** — the BFS hop distance between two selected proteins
//!   over the thresholded edge graph (∞ / "unreachable" when disconnected).
//!   The PIN asserts the hop count on a known graph and that a disconnected
//!   pair yields no path (an in-panel notice, not a panic).
//!
//! The user picks the analysis (degree / betweenness / shortest path), the
//! demo network size, and the edge-highlight threshold, clicks **Run**, and the
//! panel renders the network with a deterministic circular layout (positions
//! derived from node index — no `Math.random`), node size + colour keyed to
//! centrality, the chosen shortest path (or above-threshold edges) highlighted,
//! and a readout: the top-centrality proteins, the edge count, the path length,
//! and valenx-ppi's standing **review-required** banner.
//!
//! Mirrors the other workbenches (`cosim_workbench`, `uq_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_ppi_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"ppi"` (aliases
//! `"interactome"` / `"network"`; see [`crate::project_tabs::TabKind`]). Every
//! numeric control is `.labelled_by` an accessible caption so the panel is
//! AI-drivable by name.
//!
//! Honesty: valenx-ppi is a **research / educational coevolution heuristic**
//! that *ranks candidate interactions for a human to triage — it never emits a
//! validated "interacts" verdict* ([`PpiScore::requires_review`] is always
//! `true`, and so is [`RankedInteractions::requires_review`]). Plain
//! APC-corrected MI does not separate direct from transitive couplings the way
//! a full DCA model does, and the planted toy interactome here proves the
//! ranking + the graph metrics are wired up and behave monotonically — **not**
//! that the method achieves any accuracy on real proteomes. Every predicted
//! edge needs orthogonal evidence and wet-lab validation. A degenerate input
//! (zero nodes, or a shortest path between two disconnected proteins) surfaces
//! an in-panel notice / `Err` — **not** a panic. The tests pin the
//! hub-has-max-degree-centrality and the shortest-path hop-count PINs and the
//! degenerate-handling.

use eframe::egui;
use valenx_align::msa::Msa;
use valenx_ppi::{interactome_screen, predict_contacts, PairedMsa, ScreenEntry};

use crate::agent_commands::AgentValue;
use crate::ValenxApp;

// ---------------------------------------------------------------------------
// The named demo interactome (host × pathogen orthologue panel)
//
// We generate a paired-orthologue MSA per protein DETERMINISTICALLY: a protein
// is described by a per-column "channel" pattern over `DEPTH` organisms. A
// host and a pathogen *coevolve* (a real coevolution signal valenx-ppi will
// pick up) when they share a varying channel on a column; otherwise their
// columns are conserved (zero mutual information). No RNG — the whole panel is
// a fixed function of the protein names below.
// ---------------------------------------------------------------------------

/// Organisms (rows) in every orthologue alignment. Deep enough that mutual
/// information is meaningful (valenx-ppi's floor is
/// [`valenx_ppi::MIN_PAIRED_DEPTH`] = 3; we use far more) yet tiny + fixed.
const DEPTH: usize = 12;

/// Alignment width (columns) of every protein's orthologue MSA.
const WIDTH: usize = 6;

/// Host protein names, in node order. The first — `GUARD` — is the planted
/// **hub**: it coevolves with *every* pathogen effector, so on the interaction
/// network it has the maximum degree.
const HOST_NAMES: [&str; 3] = ["GUARD", "KINASE", "RECEPTOR"];

/// Pathogen ("effector") protein names, in node order.
const PATHOGEN_NAMES: [&str; 4] = ["EFF-A", "EFF-B", "EFF-C", "EFF-D"];

/// The two varying two-state channels a coevolving pair can share. Channel 0
/// flips on `k % 2`, channel 1 on `(k / 2) % 2`, so the two are independent
/// patterns down the `DEPTH` organisms.
const HOST_STATES: [[u8; 2]; 2] = [[b'A', b'T'], [b'D', b'E']];
/// The pathogen residues coupled to [`HOST_STATES`] (a fixed 1:1 substitution,
/// so a shared channel is a *perfect* coevolution signal valenx-ppi recovers).
const PATH_STATES: [[u8; 2]; 2] = [[b'C', b'G'], [b'K', b'R']];

/// Value of the varying channel `chan` for organism `k` (host side).
fn host_channel(chan: usize, k: usize) -> u8 {
    let s = if chan == 0 { k % 2 } else { (k / 2) % 2 };
    HOST_STATES[chan][s]
}

/// Value of the varying channel `chan` for organism `k` (pathogen side) —
/// the residue coupled 1:1 to [`host_channel`].
fn path_channel(chan: usize, k: usize) -> u8 {
    let s = if chan == 0 { k % 2 } else { (k / 2) % 2 };
    PATH_STATES[chan][s]
}

/// The planted coevolution map: which `(host i, pathogen j)` pairs genuinely
/// coevolve, and on which shared channel. `GUARD` (host 0) couples to every
/// effector (so it is the hub); `KINASE` (host 1) couples to `EFF-A` only;
/// `RECEPTOR` (host 2) couples to `EFF-C` only. A pair not listed here shares
/// no varying channel -> conserved columns -> ~zero coevolution.
fn coupled_channel(host: usize, pathogen: usize) -> Option<usize> {
    match (host, pathogen) {
        // GUARD <-> every effector (the hub). Alternate the two channels so
        // the columns carrying the signal differ across partners.
        (0, j) => Some(j % 2),
        // KINASE <-> EFF-A.
        (1, 0) => Some(0),
        // RECEPTOR <-> EFF-C.
        (2, 2) => Some(1),
        _ => None,
    }
}

/// Build host protein `i`'s orthologue [`Msa`]. Column `c` carries host
/// channel `c` (varying) iff *some* pathogen couples to this host on channel
/// `c`; otherwise it is conserved (`'M'`). Deterministic.
fn host_msa(i: usize, n_pathogen: usize) -> Result<Msa, String> {
    // Which channels does this host use? (the set of channels it couples on)
    let uses_channel = |chan: usize| (0..n_pathogen).any(|j| coupled_channel(i, j) == Some(chan));
    let mut rows: Vec<Vec<u8>> = Vec::with_capacity(DEPTH);
    for k in 0..DEPTH {
        let mut row = Vec::with_capacity(WIDTH);
        for c in 0..WIDTH {
            if c < 2 && uses_channel(c) {
                row.push(host_channel(c, k));
            } else {
                row.push(b'M'); // conserved
            }
        }
        rows.push(row);
    }
    Msa::new(rows).map_err(|e| format!("host MSA {} ({}) invalid: {e}", i, HOST_NAMES[i]))
}

/// Build pathogen protein `j`'s orthologue [`Msa`]. Column `c` carries the
/// pathogen residue coupled to host channel `c` iff *some* host couples to
/// this pathogen on channel `c`; otherwise conserved (`'K'`). Deterministic.
fn pathogen_msa(j: usize, n_host: usize) -> Result<Msa, String> {
    let uses_channel = |chan: usize| (0..n_host).any(|i| coupled_channel(i, j) == Some(chan));
    let mut rows: Vec<Vec<u8>> = Vec::with_capacity(DEPTH);
    for k in 0..DEPTH {
        let mut row = Vec::with_capacity(WIDTH);
        for c in 0..WIDTH {
            if c < 2 && uses_channel(c) {
                row.push(path_channel(c, k));
            } else {
                row.push(b'K'); // conserved
            }
        }
        rows.push(row);
    }
    Msa::new(rows).map_err(|e| format!("pathogen MSA {} ({}) invalid: {e}", j, PATHOGEN_NAMES[j]))
}

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Which network analysis the workbench computes over the PPI interaction
/// graph. A thin UI enum (so it derives the traits the egui `selectable_value`
/// widget needs and the label text lives here).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum PpiAnalysis {
    /// Per-node **degree centrality** (`deg / (N - 1)`) over the thresholded
    /// edge graph. The hub of a hub-and-spoke interactome scores highest.
    #[default]
    DegreeCentrality,
    /// Per-node **betweenness centrality** — the share of all-pairs shortest
    /// paths that pass through each node (Brandes' algorithm on the
    /// unweighted thresholded graph). A bridging hub scores highest.
    Betweenness,
    /// **Shortest path** (BFS hop distance) between the two selected node
    /// indices over the thresholded edge graph.
    ShortestPath,
}

impl PpiAnalysis {
    /// Human-readable label for the combo box / status line.
    fn label(self) -> &'static str {
        match self {
            PpiAnalysis::DegreeCentrality => "Degree centrality",
            PpiAnalysis::Betweenness => "Betweenness centrality",
            PpiAnalysis::ShortestPath => "Shortest path (BFS hops)",
        }
    }

    /// Parse an analysis name (for the agent `SetControl` bridge) into a
    /// [`PpiAnalysis`]. Case-insensitive; accepts the short menu words. Fail-loud
    /// on an unrecognised name so a typo is a `warn` note, not a silent no-op.
    fn from_name(s: &str) -> Result<PpiAnalysis, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "degree" | "degree centrality" | "degreecentrality" => {
                Ok(PpiAnalysis::DegreeCentrality)
            }
            "betweenness" | "betweenness centrality" => Ok(PpiAnalysis::Betweenness),
            "shortest" | "shortest path" | "shortestpath" | "path" => Ok(PpiAnalysis::ShortestPath),
            other => Err(format!(
                "unknown analysis '{other}' (expected 'degree', 'betweenness', or \
                 'shortest path')"
            )),
        }
    }
}

/// Editable interactome inputs shown in the workbench.
#[derive(Clone, Copy, Debug)]
pub struct PpiParams {
    /// Number of host proteins to include from the demo panel (`1..=`
    /// [`HOST_NAMES`]`.len()`). The pathogen count is fixed to the full
    /// effector panel so the hub keeps its spokes.
    pub n_hosts: usize,
    /// Number of pathogen effectors to include (`1..=`[`PATHOGEN_NAMES`]`.len()`).
    pub n_pathogens: usize,
    /// Which analysis to compute over the resulting interaction network.
    pub analysis: PpiAnalysis,
    /// Edge-highlight / graph threshold on [`PpiScore::value`](valenx_ppi::PpiScore):
    /// a host × pathogen pair becomes a network **edge** (and is drawn
    /// highlighted) iff its real fused score is `>= threshold`. In `[0, 1]`.
    pub threshold: f64,
    /// Source node index (into the combined node list `hosts ++ pathogens`)
    /// for the [`PpiAnalysis::ShortestPath`] query.
    pub path_from: usize,
    /// Target node index for the shortest-path query.
    pub path_to: usize,
}

impl Default for PpiParams {
    fn default() -> Self {
        Self {
            n_hosts: HOST_NAMES.len(),
            n_pathogens: PATHOGEN_NAMES.len(),
            analysis: PpiAnalysis::DegreeCentrality,
            // A mid threshold: planted coevolving pairs clear it, conserved
            // noise pairs do not, so the network is the planted one.
            threshold: 0.3,
            // GUARD (node 0) -> EFF-A (the first pathogen node) — a direct edge.
            path_from: 0,
            path_to: HOST_NAMES.len(), // first pathogen index in the combined list
        }
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// One node of the computed interaction network: a protein, its layout
/// position, and its centrality score for the selected analysis.
#[derive(Clone, Debug)]
pub struct PpiNode {
    /// Protein label (e.g. `"GUARD"`, `"EFF-A"`).
    pub name: String,
    /// `true` for a host protein, `false` for a pathogen effector (drives the
    /// two-colour split in the layout).
    pub is_host: bool,
    /// Number of incident above-threshold edges.
    pub degree: usize,
    /// The selected centrality measure for this node, in `[0, 1]`
    /// (degree centrality `deg/(N-1)`, or normalised betweenness). `0` for the
    /// shortest-path analysis (which has no per-node centrality).
    pub centrality: f64,
}

/// One above-threshold network edge: the two node indices it joins and the
/// real valenx-ppi fused score that produced it.
#[derive(Clone, Copy, Debug)]
pub struct PpiEdge {
    /// Index of the host endpoint in the combined node list.
    pub a: usize,
    /// Index of the pathogen endpoint in the combined node list.
    pub b: usize,
    /// The real [`PpiScore::value`](valenx_ppi::PpiScore) for this pair.
    pub score: f64,
}

/// The computed interactome + selected analysis. Edge scores come straight
/// from valenx-ppi; the graph metrics are standard derivations over them.
#[derive(Default, Clone)]
pub struct PpiResult {
    /// Network nodes (`hosts` first, then `pathogens`), in combined-index
    /// order — the index space the shortest-path endpoints use.
    pub nodes: Vec<PpiNode>,
    /// Above-threshold edges.
    pub edges: Vec<PpiEdge>,
    /// Total host × pathogen pairs scored by the screen (incl. below-threshold).
    pub pairs_scored: usize,
    /// The single highest fused score over all pairs (the strongest candidate
    /// interaction).
    pub top_score: f64,
    /// `(host, pathogen)` node labels of that top pair.
    pub top_pair: (String, String),
    /// Number of predicted interface contacts ([`predict_contacts`]) for the
    /// top pair — a glimpse of *where* the interface is.
    pub top_contacts: usize,
    /// For [`PpiAnalysis::ShortestPath`]: the node-index path from
    /// `path_from` to `path_to` (inclusive endpoints), or `None` when the two
    /// are disconnected over the thresholded graph.
    pub path: Option<Vec<usize>>,
    /// Which analysis produced this result.
    pub analysis: PpiAnalysis,
}

impl PpiResult {
    /// The number of network edges.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// The top-`k` nodes by centrality (descending; ties broken by node
    /// index for determinism).
    pub fn top_centrality(&self, k: usize) -> Vec<&PpiNode> {
        let mut idx: Vec<usize> = (0..self.nodes.len()).collect();
        idx.sort_by(|&i, &j| {
            self.nodes[j]
                .centrality
                .partial_cmp(&self.nodes[i].centrality)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(i.cmp(&j))
        });
        idx.into_iter().take(k).map(|i| &self.nodes[i]).collect()
    }

    /// The shortest-path hop length (`path.len() - 1`), or `None` if no path /
    /// not a shortest-path analysis.
    pub fn path_hops(&self) -> Option<usize> {
        self.path.as_ref().map(|p| p.len().saturating_sub(1))
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the PPI / interactome workbench.
#[derive(Default)]
pub struct PpiWorkbenchState {
    /// User-editable parameters.
    pub params: PpiParams,
    /// Last successful result (populated after a successful Run).
    pub result: Option<PpiResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

// ---------------------------------------------------------------------------
// Run the interactome through the REAL valenx-ppi engine
// ---------------------------------------------------------------------------

impl PpiWorkbenchState {
    /// Build the demo interactome, screen it through the **real**
    /// [`interactome_screen`], assemble the interaction network from the real
    /// scored pairs, and compute the selected graph analysis — fail-loud.
    ///
    /// Every failure path returns an `Err(String)` — never a panic, never an
    /// invented number. Degenerate inputs (zero hosts/pathogens, a
    /// shortest-path endpoint out of range) are rejected up front; a
    /// disconnected shortest-path pair returns a result whose
    /// [`PpiResult::path`] is `None` (rendered as an in-panel notice, not an
    /// error).
    pub fn run(&self) -> Result<PpiResult, String> {
        let p = &self.params;

        // --- Degenerate-input guards (fail loud, no panic) ------------------
        let max_hosts = HOST_NAMES.len();
        let max_path = PATHOGEN_NAMES.len();
        if p.n_hosts == 0 || p.n_hosts > max_hosts {
            return Err(format!(
                "number of hosts must be 1..={max_hosts} (got {})",
                p.n_hosts
            ));
        }
        if p.n_pathogens == 0 || p.n_pathogens > max_path {
            return Err(format!(
                "number of pathogens must be 1..={max_path} (got {})",
                p.n_pathogens
            ));
        }
        if !p.threshold.is_finite() || !(0.0..=1.0).contains(&p.threshold) {
            return Err(format!(
                "edge threshold must be finite and in [0, 1] (got {})",
                p.threshold
            ));
        }
        let n_nodes = p.n_hosts + p.n_pathogens;
        if p.analysis == PpiAnalysis::ShortestPath
            && (p.path_from >= n_nodes || p.path_to >= n_nodes)
        {
            return Err(format!(
                "shortest-path endpoints must be node indices in 0..{n_nodes} \
                 (got from={}, to={})",
                p.path_from, p.path_to
            ));
        }

        // --- Build the demo orthologue panel (host + pathogen MSAs) ---------
        let mut host: Vec<ScreenEntry> = Vec::with_capacity(p.n_hosts);
        for (i, &name) in HOST_NAMES.iter().enumerate().take(p.n_hosts) {
            host.push(ScreenEntry::new(name, host_msa(i, p.n_pathogens)?));
        }
        let mut pathogen: Vec<ScreenEntry> = Vec::with_capacity(p.n_pathogens);
        for (j, &name) in PATHOGEN_NAMES.iter().enumerate().take(p.n_pathogens) {
            pathogen.push(ScreenEntry::new(name, pathogen_msa(j, p.n_hosts)?));
        }

        // --- Screen the interactome through the REAL valenx-ppi engine ------
        let screen = interactome_screen(&host, &pathogen)
            .map_err(|e| format!("interactome screen failed: {e}"))?;

        // --- Nodes: hosts first, then pathogens (combined index space) ------
        let mut nodes: Vec<PpiNode> = Vec::with_capacity(n_nodes);
        for &name in HOST_NAMES.iter().take(p.n_hosts) {
            nodes.push(PpiNode {
                name: name.to_string(),
                is_host: true,
                degree: 0,
                centrality: 0.0,
            });
        }
        for &name in PATHOGEN_NAMES.iter().take(p.n_pathogens) {
            nodes.push(PpiNode {
                name: name.to_string(),
                is_host: false,
                degree: 0,
                centrality: 0.0,
            });
        }

        // --- Edges: every above-threshold real PPI score becomes an edge ----
        let mut edges: Vec<PpiEdge> = Vec::new();
        let mut top_score = f64::NEG_INFINITY;
        let mut top_pair = (String::new(), String::new());
        for inter in &screen.ranked {
            let val = inter.score.value;
            if val > top_score {
                top_score = val;
                top_pair = (
                    HOST_NAMES[inter.host].to_string(),
                    PATHOGEN_NAMES[inter.pathogen].to_string(),
                );
            }
            if val >= p.threshold {
                let a = inter.host; // host node index
                let b = p.n_hosts + inter.pathogen; // pathogen node index
                edges.push(PpiEdge { a, b, score: val });
            }
        }
        if top_score == f64::NEG_INFINITY {
            top_score = 0.0;
        }

        // Adjacency (undirected) over the thresholded edges for the metrics.
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
        for e in &edges {
            adj[e.a].push(e.b);
            adj[e.b].push(e.a);
            nodes[e.a].degree += 1;
            nodes[e.b].degree += 1;
        }

        // --- The selected analysis over the REAL interaction network --------
        let mut path: Option<Vec<usize>> = None;
        match p.analysis {
            PpiAnalysis::DegreeCentrality => {
                let denom = (n_nodes.saturating_sub(1)).max(1) as f64;
                for node in &mut nodes {
                    node.centrality = node.degree as f64 / denom;
                }
            }
            PpiAnalysis::Betweenness => {
                let bc = betweenness_centrality(&adj);
                for (node, &b) in nodes.iter_mut().zip(&bc) {
                    node.centrality = b;
                }
            }
            PpiAnalysis::ShortestPath => {
                path = bfs_path(&adj, p.path_from, p.path_to);
            }
        }

        // --- Where is the interface? contacts for the top pair --------------
        // (predict_contacts on the strongest pair — a glimpse of the interface;
        // count the predicted contacts, fail-loud on a bad pairing.)
        let top_contacts = if !top_pair.0.is_empty() {
            let hi = HOST_NAMES
                .iter()
                .position(|&n| n == top_pair.0)
                .unwrap_or(0);
            let pj = PATHOGEN_NAMES
                .iter()
                .position(|&n| n == top_pair.1)
                .unwrap_or(0);
            let paired = PairedMsa::new(host_msa(hi, p.n_pathogens)?, pathogen_msa(pj, p.n_hosts)?)
                .map_err(|e| format!("re-pairing the top pair failed: {e}"))?;
            let coev = predict_contacts(&paired)
                .map_err(|e| format!("contact prediction for the top pair failed: {e}"))?;
            // Count contacts whose APC-corrected score is positive (a coupling
            // the model actually flags) — a glimpse of the interface size.
            coev.ranked.iter().filter(|c| c.score > 0.0).count()
        } else {
            0
        };

        Ok(PpiResult {
            nodes,
            edges,
            pairs_scored: screen.ranked.len(),
            top_score,
            top_pair,
            top_contacts,
            path,
            analysis: p.analysis,
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. The captions match exactly what
    /// the workbench form draws (and what each control is `labelled_by`).
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "# host proteins",
            "# pathogen proteins",
            "analysis",
            "edge threshold",
            "path from (node #)",
            "path to (node #)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the wrong
    /// type returns `Err(String)` (the bridge turns it into a `warn` feed note) —
    /// never a panic, and no field is written on error. The numeric captions read
    /// [`AgentValue::as_i64`] / [`AgentValue::as_f64`]; the `analysis` enum caption
    /// reads [`AgentValue::as_str`].
    pub fn agent_set(&mut self, name: &str, value: &AgentValue) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            "# host proteins" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("# host proteins must be >= 0, got {n}"));
                }
                p.n_hosts = n as usize;
            }
            "# pathogen proteins" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("# pathogen proteins must be >= 0, got {n}"));
                }
                p.n_pathogens = n as usize;
            }
            "analysis" => p.analysis = PpiAnalysis::from_name(value.as_str()?)?,
            "edge threshold" => p.threshold = value.as_f64()?,
            "path from (node #)" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("path from (node #) must be >= 0, got {n}"));
                }
                p.path_from = n as usize;
            }
            "path to (node #)" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("path to (node #) must be >= 0, got {n}"));
                }
                p.path_to = n as usize;
            }
            other => return Err(format!("unknown PPI control: {other:?}")),
        }
        Ok(())
    }
}

/// BFS shortest path (fewest hops) between `from` and `to` over an undirected
/// adjacency list. Returns the node-index path inclusive of both endpoints, or
/// `None` if `to` is unreachable from `from`. `from == to` is a zero-hop path
/// of one node.
fn bfs_path(adj: &[Vec<usize>], from: usize, to: usize) -> Option<Vec<usize>> {
    if from >= adj.len() || to >= adj.len() {
        return None;
    }
    if from == to {
        return Some(vec![from]);
    }
    let mut prev: Vec<Option<usize>> = vec![None; adj.len()];
    let mut seen = vec![false; adj.len()];
    let mut queue = std::collections::VecDeque::new();
    seen[from] = true;
    queue.push_back(from);
    while let Some(u) = queue.pop_front() {
        for &v in &adj[u] {
            if !seen[v] {
                seen[v] = true;
                prev[v] = Some(u);
                if v == to {
                    // Reconstruct the path back to `from`.
                    let mut path = vec![to];
                    let mut cur = to;
                    while let Some(p) = prev[cur] {
                        path.push(p);
                        cur = p;
                    }
                    path.reverse();
                    return Some(path);
                }
                queue.push_back(v);
            }
        }
    }
    None
}

/// Brandes' betweenness centrality on an unweighted, undirected graph,
/// normalised to `[0, 1]` by the `(N-1)(N-2)` undirected pair count. `0` for
/// every node when the graph is too small to have a through-path.
fn betweenness_centrality(adj: &[Vec<usize>]) -> Vec<f64> {
    let n = adj.len();
    let mut bc = vec![0.0f64; n];
    if n < 3 {
        return bc;
    }
    for s in 0..n {
        // Single-source shortest-path counts (Brandes).
        let mut stack: Vec<usize> = Vec::new();
        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut sigma = vec![0.0f64; n];
        let mut dist = vec![-1i64; n];
        sigma[s] = 1.0;
        dist[s] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(s);
        while let Some(v) = queue.pop_front() {
            stack.push(v);
            for &w in &adj[v] {
                if dist[w] < 0 {
                    dist[w] = dist[v] + 1;
                    queue.push_back(w);
                }
                if dist[w] == dist[v] + 1 {
                    sigma[w] += sigma[v];
                    preds[w].push(v);
                }
            }
        }
        // Accumulation.
        let mut delta = vec![0.0f64; n];
        while let Some(w) = stack.pop() {
            for &v in &preds[w] {
                if sigma[w] != 0.0 {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if w != s {
                bc[w] += delta[w];
            }
        }
    }
    // Undirected: each pair counted twice; normalise by (n-1)(n-2).
    let norm = ((n - 1) * (n - 2)) as f64;
    if norm > 0.0 {
        for b in &mut bc {
            *b /= norm;
        }
    }
    bc
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the PPI / interactome workbench. A no-op unless toggled on via
/// View → Protein interaction (PPI).
///
/// Mirrors [`crate::cosim_workbench::draw_cosim_workbench`].
pub fn draw_ppi_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_ppi_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_ppi_workbench",
        "Protein interaction (PPI / interactome)",
        ppi_workbench_body,
    );
    if close {
        app.show_ppi_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn ppi_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Protein-protein interaction network \u{2014} a small named host \u{00D7} pathogen \
             demo interactome (deterministic coevolving orthologue alignments) screened through \
             the REAL in-house valenx-ppi engine (APC-corrected mutual-information coevolution \
             -> a fused [0,1] PPI score per pair -> a ranked all-vs-all screen). Nodes = \
             proteins, edges = above-threshold scored interactions; degree / betweenness \
             centrality + BFS shortest path are computed over valenx-ppi's real scored network. \
             [research / educational coevolution HEURISTIC \u{2014} ranks candidates, NEVER a \
             verdict; every result requires human + wet-lab review]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.ppi;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Demo interactome").strong());
        egui::Grid::new("ppi_network_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("# host proteins");
                ui.add(
                    egui::DragValue::new(&mut p.n_hosts)
                        .speed(1)
                        .range(1..=HOST_NAMES.len()),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "How many host proteins to include from the demo panel (GUARD, KINASE, \
                     RECEPTOR). GUARD is the planted hub (coevolves with every effector). \
                     Must be 1..=3.",
                );
                ui.end_row();

                let lbl = ui.label("# pathogen proteins");
                ui.add(
                    egui::DragValue::new(&mut p.n_pathogens)
                        .speed(1)
                        .range(1..=PATHOGEN_NAMES.len()),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "How many pathogen effectors to include (EFF-A..EFF-D). The total node \
                     count is #hosts + #pathogens. Must be 1..=4.",
                );
                ui.end_row();

                let lbl = ui.label("analysis");
                egui::ComboBox::from_id_source("ppi_analysis_combo")
                    .selected_text(p.analysis.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut p.analysis,
                            PpiAnalysis::DegreeCentrality,
                            PpiAnalysis::DegreeCentrality.label(),
                        );
                        ui.selectable_value(
                            &mut p.analysis,
                            PpiAnalysis::Betweenness,
                            PpiAnalysis::Betweenness.label(),
                        );
                        ui.selectable_value(
                            &mut p.analysis,
                            PpiAnalysis::ShortestPath,
                            PpiAnalysis::ShortestPath.label(),
                        );
                    })
                    .response
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Which graph analysis to compute over the real PPI network: per-node \
                         degree or betweenness centrality, or the BFS shortest path between two \
                         selected proteins.",
                    );
                ui.end_row();

                let lbl = ui.label("edge threshold");
                ui.add(
                    egui::DragValue::new(&mut p.threshold)
                        .speed(0.01)
                        .range(0.0..=1.0)
                        .max_decimals(3),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "A host \u{00D7} pathogen pair becomes a network edge (and is highlighted) \
                     when its real valenx-ppi fused score is >= this value. In [0, 1]. \
                     Nodes' degree / centrality are over the thresholded graph.",
                );
                ui.end_row();
            });

        // Shortest-path endpoints. Shown always (greyed unless the analysis is
        // ShortestPath) so the form layout + accessible names stay stable
        // (mirrors cosim_workbench's enabled-ui pattern).
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Shortest-path endpoints").strong());
        egui::Grid::new("ppi_path_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let is_path = p.analysis == PpiAnalysis::ShortestPath;
                let max_node = (p.n_hosts + p.n_pathogens).saturating_sub(1);
                ui.add_enabled_ui(is_path, |ui| {
                    let lbl = ui.label("path from (node #)");
                    ui.add(
                        egui::DragValue::new(&mut p.path_from)
                            .speed(1)
                            .range(0..=max_node),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Source node index for the shortest-path query. Node order is hosts \
                         first (0..#hosts), then pathogens. Must be a valid node index.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_path, |ui| {
                    let lbl = ui.label("path to (node #)");
                    ui.add(
                        egui::DragValue::new(&mut p.path_to)
                            .speed(1)
                            .range(0..=max_node),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Target node index for the shortest-path query. A pair with no path \
                         over the thresholded edges is reported as unreachable (an in-panel \
                         notice, not an error).",
                    );
                });
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Screen the demo interactome through valenx-ppi, assemble the interaction \
                     network from the real scored pairs, and compute the selected centrality / \
                     shortest-path analysis.",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside the params borrow) --------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.ppi;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_ppi_viz(s, ui);
}

/// Run the screen and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.ppi;
    match s.run() {
        Ok(res) => {
            let analysis = res.analysis.label();
            if res.analysis == PpiAnalysis::ShortestPath {
                let path_txt = match res.path_hops() {
                    Some(h) => format!("{h} hop(s)"),
                    None => "unreachable".to_string(),
                };
                s.status = format!(
                    "\u{2714} {} nodes \u{00B7} {} edges \u{00B7} {} \u{00B7} path {} \u{00B7} \
                     REVIEW REQUIRED",
                    res.nodes.len(),
                    res.n_edges(),
                    analysis,
                    path_txt,
                );
            } else {
                let top = res
                    .top_centrality(1)
                    .first()
                    .map(|n| n.name.clone())
                    .unwrap_or_default();
                s.status = format!(
                    "\u{2714} {} nodes \u{00B7} {} edges \u{00B7} {} \u{00B7} top: {} \u{00B7} \
                     REVIEW REQUIRED",
                    res.nodes.len(),
                    res.n_edges(),
                    analysis,
                    top,
                );
            }
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (painter network graph + readout)
// ---------------------------------------------------------------------------

fn draw_ppi_viz(s: &PpiWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to screen the demo interactome through valenx-ppi and draw the \
                 protein-interaction network",
            )
            .weak(),
        );
        return;
    };

    ui.label(egui::RichText::new("Interaction network").strong());
    ui.label(
        egui::RichText::new(
            "cyan nodes = host proteins \u{00B7} amber nodes = pathogen effectors \u{00B7} node \
             size + brightness = centrality \u{00B7} green edges = shortest path; grey = other \
             above-threshold edges",
        )
        .weak()
        .small(),
    );

    draw_network(res, ui);

    // Readouts grid below the graph.
    ui.add_space(6.0);
    egui::Grid::new("ppi_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "proteins (nodes)", format!("{}", res.nodes.len()));
            row(
                ui,
                "interactions (edges \u{2265} threshold)",
                format!("{} of {} pairs", res.n_edges(), res.pairs_scored),
            );
            row(ui, "analysis", res.analysis.label().to_string());
            row(
                ui,
                "strongest pair",
                format!(
                    "{} \u{2013} {} (score {:.3}, {} contacts)",
                    res.top_pair.0, res.top_pair.1, res.top_score, res.top_contacts
                ),
            );

            if res.analysis == PpiAnalysis::ShortestPath {
                let names: Vec<String> = match &res.path {
                    Some(p) => p
                        .iter()
                        .map(|&i| res.nodes.get(i).map(|n| n.name.clone()).unwrap_or_default())
                        .collect(),
                    None => Vec::new(),
                };
                match res.path_hops() {
                    Some(h) => {
                        row(ui, "shortest path", format!("{h} hop(s)"));
                        row(ui, "route", names.join(" \u{2192} "));
                    }
                    None => {
                        row(
                            ui,
                            "shortest path",
                            "unreachable (disconnected over the thresholded graph)".to_string(),
                        );
                    }
                }
            } else {
                let top: Vec<String> = res
                    .top_centrality(3)
                    .iter()
                    .map(|n| format!("{} ({:.2})", n.name, n.centrality))
                    .collect();
                row(ui, "top-centrality proteins", top.join(", "));
            }
            row(
                ui,
                "status",
                "review required \u{2014} candidates, not a verdict".to_string(),
            );
        });
}

/// Draw the interaction network with the egui painter: a deterministic
/// circular layout (positions derived from node index — no RNG), node size +
/// colour keyed to centrality, the shortest path (or above-threshold edges)
/// highlighted. Self-contained — no graph-drawing dependency.
fn draw_network(res: &PpiResult, ui: &mut egui::Ui) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(460.0, 300.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    let n = res.nodes.len();
    if n == 0 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "empty network (no nodes)",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Deterministic circular layout: node i at angle 2*pi*i/n on a centred
    // circle. Purely a function of the index — reproducible, no Math.random.
    let center = rect.center();
    let radius = (rect.width().min(rect.height()) * 0.5 - 36.0).max(10.0);
    let pos = |i: usize| -> egui::Pos2 {
        if n == 1 {
            return center;
        }
        let theta = std::f32::consts::TAU * (i as f32) / (n as f32) - std::f32::consts::FRAC_PI_2;
        egui::pos2(
            center.x + radius * theta.cos(),
            center.y + radius * theta.sin(),
        )
    };

    // Edges that lie on the shortest path (drawn highlighted). Build a set of
    // unordered index pairs from the path node sequence.
    let path_edges: std::collections::HashSet<(usize, usize)> = match &res.path {
        Some(p) => p
            .windows(2)
            .map(|w| (w[0].min(w[1]), w[0].max(w[1])))
            .collect(),
        None => std::collections::HashSet::new(),
    };

    // --- Edges first (under the nodes) --------------------------------------
    for e in &res.edges {
        let key = (e.a.min(e.b), e.a.max(e.b));
        let on_path = path_edges.contains(&key);
        let color = if on_path {
            egui::Color32::from_rgb(90, 220, 120) // green — the shortest path
        } else {
            // Grey, with opacity scaled by the edge's PPI score.
            let g = (90.0 + 120.0 * e.score.clamp(0.0, 1.0)) as u8;
            egui::Color32::from_rgb(g / 2, g, g)
        };
        let width = if on_path { 2.6 } else { 1.4 };
        painter.line_segment([pos(e.a), pos(e.b)], egui::Stroke::new(width, color));
    }

    // --- Nodes on top -------------------------------------------------------
    for (i, node) in res.nodes.iter().enumerate() {
        let pcenter = pos(i);
        // Radius scales with centrality (degree/betweenness); shortest-path
        // analysis has zero centrality, so use the raw degree for sizing then.
        let mag = if res.analysis == PpiAnalysis::ShortestPath {
            // Normalise degree to [0,1] for sizing in path mode.
            let denom = (n.saturating_sub(1)).max(1) as f64;
            (node.degree as f64 / denom).clamp(0.0, 1.0)
        } else {
            node.centrality.clamp(0.0, 1.0)
        };
        let r = 6.0 + 10.0 * mag as f32;

        // Host = cyan, pathogen = amber; brightness rises with centrality.
        let base = if node.is_host {
            (70u8, 200u8, 210u8)
        } else {
            (230u8, 180u8, 70u8)
        };
        let lift = (0.55 + 0.45 * mag) as f32;
        let color = egui::Color32::from_rgb(
            (base.0 as f32 * lift) as u8,
            (base.1 as f32 * lift) as u8,
            (base.2 as f32 * lift) as u8,
        );

        // Highlight nodes on the shortest path with a green ring.
        let on_path = res.path.as_ref().is_some_and(|p| p.contains(&i));
        if on_path {
            painter.circle_stroke(
                pcenter,
                r + 2.5,
                egui::Stroke::new(2.0, egui::Color32::from_rgb(90, 220, 120)),
            );
        }
        painter.circle_filled(pcenter, r, color);

        // Label just outside the node.
        let label_pos = egui::pos2(pcenter.x, pcenter.y - r - 7.0);
        painter.text(
            label_pos,
            egui::Align2::CENTER_CENTER,
            &node.name,
            egui::FontId::monospace(11.0),
            egui::Color32::from_gray(210),
        );
    }
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring cosim_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = PpiWorkbenchState::default();
        let res = s.run().expect("default PPI run should succeed");
        assert_eq!(
            res.nodes.len(),
            HOST_NAMES.len() + PATHOGEN_NAMES.len(),
            "all hosts + pathogens become nodes"
        );
        assert_eq!(
            res.pairs_scored,
            HOST_NAMES.len() * PATHOGEN_NAMES.len(),
            "every host x pathogen pair is scored by the real screen"
        );
        // The planted coevolving pairs clear the default threshold -> edges.
        assert!(res.n_edges() >= 1, "planted coevolving pairs become edges");
        // valenx-ppi's review flag is reflected: top score is a real [0,1] value.
        assert!((0.0..=1.0).contains(&res.top_score));
        assert_eq!(res.analysis, PpiAnalysis::DegreeCentrality);
    }

    #[test]
    fn hub_has_max_degree_centrality_pin() {
        // PIN (analytic): on the planted hub-and-spoke interactome the hub
        // (GUARD, node 0) coevolves with EVERY pathogen effector, so it has
        // the maximum degree and thus the maximum degree centrality. No other
        // node may strictly exceed it.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::DegreeCentrality;
        let res = s.run().expect("degree-centrality run should succeed");

        // GUARD is node 0.
        assert_eq!(res.nodes[0].name, "GUARD");
        let hub_c = res.nodes[0].centrality;
        let hub_deg = res.nodes[0].degree;

        // The hub links to all pathogens.
        assert_eq!(
            hub_deg,
            PpiParams::default().n_pathogens,
            "the hub GUARD should connect to every pathogen effector"
        );
        // Max centrality overall.
        let max_c = res
            .nodes
            .iter()
            .map(|n| n.centrality)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (hub_c - max_c).abs() < 1e-12,
            "hub centrality {hub_c} must be the maximum {max_c}"
        );
        // And strictly greater than at least one spoke (the network is not a
        // trivial clique).
        assert!(
            res.nodes.iter().skip(1).any(|n| n.centrality < hub_c),
            "the hub must out-rank some spoke"
        );
    }

    #[test]
    fn shortest_path_hop_count_pin() {
        // PIN (analytic): GUARD (node 0) <-> EFF-A (node = n_hosts, the first
        // pathogen) is a DIRECT planted edge, so the shortest path is exactly
        // ONE hop. And GUARD -> KINASE (node 1), two hosts that share no direct
        // edge but both touch EFF-A, is exactly TWO hops (GUARD -> EFF-A ->
        // KINASE).
        let n_hosts = HOST_NAMES.len();

        // (a) one-hop direct edge.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::ShortestPath;
        s.params.path_from = 0; // GUARD
        s.params.path_to = n_hosts; // EFF-A
        let res = s.run().expect("shortest-path run should succeed");
        assert_eq!(
            res.path_hops(),
            Some(1),
            "GUARD -> EFF-A is a direct (1-hop) interaction; got {:?}",
            res.path
        );

        // (b) two-hop path GUARD -> EFF-A -> KINASE.
        let mut s2 = PpiWorkbenchState::default();
        s2.params.analysis = PpiAnalysis::ShortestPath;
        s2.params.path_from = 0; // GUARD
        s2.params.path_to = 1; // KINASE
        let res2 = s2.run().expect("shortest-path run should succeed");
        assert_eq!(
            res2.path_hops(),
            Some(2),
            "GUARD -> KINASE should be 2 hops via EFF-A; got {:?}",
            res2.path
        );
    }

    #[test]
    fn disconnected_pair_has_no_path_no_panic() {
        // A degenerate query: raise the threshold so high that NO pair clears
        // it (no edges at all). Then any two distinct nodes are disconnected,
        // so the shortest path is `None` — surfaced, never a panic.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::ShortestPath;
        s.params.threshold = 1.0; // nothing scores a perfect 1.0 here
        s.params.path_from = 0;
        s.params.path_to = 1;
        let res = s
            .run()
            .expect("run should still succeed (empty graph is valid)");
        assert_eq!(res.n_edges(), 0, "threshold 1.0 admits no edges");
        assert!(
            res.path.is_none(),
            "disconnected nodes must yield no path, got {:?}",
            res.path
        );
        assert_eq!(res.path_hops(), None);
    }

    #[test]
    fn self_path_is_zero_hops() {
        // from == to is a valid zero-hop path of one node (not a panic).
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::ShortestPath;
        s.params.path_from = 0;
        s.params.path_to = 0;
        let res = s.run().expect("self-path run should succeed");
        assert_eq!(res.path_hops(), Some(0));
        assert_eq!(res.path.as_deref(), Some(&[0usize][..]));
    }

    #[test]
    fn higher_threshold_never_adds_edges() {
        // Monotone: raising the edge threshold can only remove edges.
        let edges_at = |t: f64| -> usize {
            let mut s = PpiWorkbenchState::default();
            s.params.threshold = t;
            s.run().expect("run should succeed").n_edges()
        };
        let lo = edges_at(0.1);
        let hi = edges_at(0.6);
        assert!(
            hi <= lo,
            "raising the threshold ({lo} -> {hi}) must not add edges"
        );
    }

    #[test]
    fn betweenness_runs_and_hub_scores_high() {
        // Betweenness on the planted graph: the hub GUARD bridges otherwise
        // disconnected effectors and hosts, so it has the max betweenness too.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::Betweenness;
        let res = s.run().expect("betweenness run should succeed");
        let hub = res.nodes[0].centrality; // GUARD
        let max_c = res
            .nodes
            .iter()
            .map(|n| n.centrality)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (hub - max_c).abs() < 1e-12,
            "hub GUARD betweenness {hub} must be the maximum {max_c}"
        );
        assert!(hub > 0.0, "the bridging hub must have positive betweenness");
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_hosts_returns_err() {
        let mut s = PpiWorkbenchState::default();
        s.params.n_hosts = 0;
        assert!(s.run().is_err(), "0 hosts must return Err, not panic");
    }

    #[test]
    fn zero_pathogens_returns_err() {
        let mut s = PpiWorkbenchState::default();
        s.params.n_pathogens = 0;
        assert!(s.run().is_err(), "0 pathogens must return Err, not panic");
    }

    #[test]
    fn too_many_hosts_returns_err() {
        let mut s = PpiWorkbenchState::default();
        s.params.n_hosts = HOST_NAMES.len() + 5;
        assert!(s.run().is_err(), "out-of-range host count must return Err");
    }

    #[test]
    fn out_of_range_threshold_returns_err() {
        let mut s = PpiWorkbenchState::default();
        s.params.threshold = 1.5;
        assert!(s.run().is_err(), "threshold > 1 must return Err");
        s.params.threshold = f64::NAN;
        assert!(s.run().is_err(), "NaN threshold must return Err");
    }

    #[test]
    fn out_of_range_path_endpoint_returns_err() {
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::ShortestPath;
        s.params.path_to = 999;
        assert!(
            s.run().is_err(),
            "out-of-range path endpoint must return Err"
        );
    }

    #[test]
    fn coevolving_pairs_outrank_noise_edges() {
        // The real valenx-ppi screen must score the planted coevolving pairs
        // above the conserved-noise pairs — so at a mid threshold the edges
        // are (a superset of) the planted interactions and exclude pure noise.
        let mut s = PpiWorkbenchState::default();
        s.params.threshold = 0.3;
        let res = s.run().expect("run should succeed");
        // GUARD (host 0) couples to all 4 effectors -> GUARD has degree 4.
        assert_eq!(res.nodes[0].name, "GUARD");
        assert_eq!(
            res.nodes[0].degree,
            PATHOGEN_NAMES.len(),
            "GUARD should edge to every effector at a mid threshold"
        );
    }

    // ---- agent_set / agent_control_names (the SetControl bridge) ----

    #[test]
    fn agent_set_sets_params_and_rejects_unknown_and_typemismatch() {
        let mut s = PpiWorkbenchState::default();

        // A representative numeric param, verified via state.
        s.agent_set("# host proteins", &AgentValue::Int(2))
            .expect("set # host proteins");
        assert_eq!(s.params.n_hosts, 2);
        // A float param.
        s.agent_set("edge threshold", &AgentValue::Float(0.42))
            .expect("set edge threshold");
        assert!((s.params.threshold - 0.42).abs() < 1e-12);
        // The enum-by-name combo.
        s.agent_set("analysis", &AgentValue::Str("shortest path".into()))
            .expect("set analysis");
        assert_eq!(s.params.analysis, PpiAnalysis::ShortestPath);

        // Unknown caption -> Err (not a panic).
        assert!(s.agent_set("nope", &AgentValue::Int(1)).is_err());
        // Type mismatch: a numeric caption fed a string -> Err.
        assert!(s
            .agent_set("# host proteins", &AgentValue::Str("two".into()))
            .is_err());
        // Type mismatch: the enum caption fed a number -> Err.
        assert!(s.agent_set("analysis", &AgentValue::Int(1)).is_err());
        // An unknown enum name -> Err.
        assert!(s
            .agent_set("analysis", &AgentValue::Str("bogus".into()))
            .is_err());

        // Every advertised control name is settable with a value of its type.
        for name in PpiWorkbenchState::agent_control_names() {
            let v = match *name {
                "analysis" => AgentValue::Str("degree".into()),
                "edge threshold" => AgentValue::Float(0.5),
                _ => AgentValue::Int(1),
            };
            assert!(
                s.agent_set(name, &v).is_ok(),
                "advertised control '{name}' must be settable"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_ppi_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_ppi_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_ppi_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        let res = app.ppi.run().expect("run should succeed");
        app.ppi.result = Some(res);
        app.ppi.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_shortest_path_result_without_panic() {
        // Exercise the shortest-path readout rows + green path highlight.
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        app.ppi.params.analysis = PpiAnalysis::ShortestPath;
        app.ppi.params.path_from = 0;
        app.ppi.params.path_to = 1;
        let res = app.ppi.run().expect("shortest-path run should succeed");
        app.ppi.result = Some(res);
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_unreachable_path_without_panic() {
        // Disconnected (threshold 1.0) shortest-path -> "unreachable" readout.
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        app.ppi.params.analysis = PpiAnalysis::ShortestPath;
        app.ppi.params.threshold = 1.0;
        let res = app.ppi.run().expect("run should succeed");
        assert!(res.path.is_none());
        app.ppi.result = Some(res);
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        // Trigger an error state (0 hosts is fail-loud in run()).
        app.ppi.params.n_hosts = 0;
        let result = app.ppi.run();
        app.ppi.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.ppi.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // The numeric DragValues (# hosts, # pathogens, edge threshold, path
        // from, path to) MUST each carry an accessible name (be labelled_by a
        // caption) so the panel is AI-drivable.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected at least 4 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check the specific captions are present as named accessibility nodes.
        for caption in [
            "# host proteins",
            "# pathogen proteins",
            "analysis",
            "edge threshold",
            "path from (node #)",
            "path to (node #)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption; each `labelled_by` target must RESOLVE to a real named
        // caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_ppi_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        for caption in ["# host proteins", "edge threshold"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn hub_centrality_pin_from_ui_state() {
        // Mirror of the unit pin, exercised from the UI-state struct: the hub
        // GUARD has the maximum degree centrality on the planted interactome.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::DegreeCentrality;
        let res = s.run().expect("degree run");
        let hub = res.nodes[0].centrality;
        let max_c = res
            .nodes
            .iter()
            .map(|n| n.centrality)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((hub - max_c).abs() < 1e-12, "hub is the centrality max");
    }

    #[test]
    fn shortest_path_pin_from_ui_state() {
        // Mirror of the shortest-path pin: GUARD -> EFF-A is one hop.
        let mut s = PpiWorkbenchState::default();
        s.params.analysis = PpiAnalysis::ShortestPath;
        s.params.path_from = 0;
        s.params.path_to = HOST_NAMES.len();
        let res = s.run().expect("path run");
        assert_eq!(res.path_hops(), Some(1));
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // Zero hosts (or an out-of-range threshold) must surface the error
        // in-panel, not panic.
        let mut state = PpiWorkbenchState::default();
        state.params.n_hosts = 0;
        assert!(state.run().is_err(), "0 hosts must produce Err, not panic");
        state.params.n_hosts = HOST_NAMES.len();
        state.params.threshold = -0.5;
        assert!(
            state.run().is_err(),
            "bad threshold must produce Err, not panic"
        );
    }

    #[test]
    fn agent_bridge_ppi_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "ppi" }`:
        //   1. TabKind::from_id("ppi") -> Some(TabKind::Ppi)
        //      (plus the aliases "interactome" / "network")
        //   2. set_workbench_flag(app, "ppi", true) -> show_ppi_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup (canonical + aliases).
        assert_eq!(
            TabKind::from_id("ppi"),
            Some(TabKind::Ppi),
            "\"ppi\" must resolve to TabKind::Ppi"
        );
        assert_eq!(TabKind::from_id("interactome"), Some(TabKind::Ppi));
        assert_eq!(TabKind::from_id("network"), Some(TabKind::Ppi));
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("  Ppi  "), Some(TabKind::Ppi));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_ppi_workbench);
        set_workbench_flag(&mut app, "ppi", true);
        assert!(
            app.show_ppi_workbench,
            "set_workbench_flag(\"ppi\", true) must set the flag"
        );
        set_workbench_flag(&mut app, "ppi", false);
        assert!(!app.show_ppi_workbench);
    }
}

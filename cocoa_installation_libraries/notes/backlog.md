# Cocoa execution backlog

Tracked in Git so unfinished work survives a new clone and a new
maintenance session. The commit that closes a ticket updates its entry
in the same change.

## Contents

- [Open tickets](#open-tickets)
- [Closed tickets](#closed-tickets)

## How to read this backlog

One line beginning exactly `- OPEN` per unfinished ticket in the index
below, each carrying `BUG FIX` or `NEW FUNCTIONALITY`. Never a second
`- OPEN` line inside a ticket. Every ticket has one anchored section
with its summary, status, and what is missing; pointers to files and
commits go in the collapsed technical record at the end of the section.

# Open tickets

Severity sorts first, type second: High bugs, Medium bugs, Medium
features, Low bugs, Low features.

- **High** — the science can be wrong, data can be lost, or a core
  operation halts.
- **Medium** — a concrete problem reasonably likely during normal
  work, below the High boundary.
- **Low** — concrete but improbable edge cases.

## Open ticket index

### Medium

- OPEN **MEDIUM** **BUG** — [Cache-key hardening for in-process reconfiguration](#open-cache-key-hardening)
- OPEN **MEDIUM** **NEW FUNCTIONALITY** — [Cache the Fourier non-Limber band-center corrections](#open-fourier-nonlimber-cache)

### Low

- OPEN **LOW** **NEW FUNCTIONALITY** — [IA x higher-order-bias (gb2) cross terms in cfastpt](#open-cfastpt-gb2-ia-bias)
- OPEN **LOW** **NEW FUNCTIONALITY** — [Contract the dlnxi/dlnw scale-cut cache builds](#open-scuts-dlnxi-contraction)
- OPEN **LOW** **NEW FUNCTIONALITY** — [Finish the Compton-y port (C_gy, C_ys, C_ky, C_yy)](#open-compton-y-port)
- OPEN **LOW** **NEW FUNCTIONALITY** — [Web-halo model (WHM, arXiv:2508.10902) for nonlinear P(k)](#open-web-halo-model)

<a id="open-cache-key-hardening"></a>
## Cache-key hardening for in-process reconfiguration

### High-level summary

The 2026-09-28 crash fixes keyed the ggl pair maps on
tomo.random_ggl; the same review found the remaining static state
that survives an in-process reconfiguration it should not survive.
None of it is reachable from the shipped pytest suites (every model
rebuild redraws Ntable.random), but notebooks driving the C inits
directly can hit each one.

### Current status

**Ticket type: BUG.**

**OPEN.** The list:

- FAST-PT table ownership (pt_cfastpt.c / generic_interface.cpp):
  `owned_tab` detects a python-installed table by pointer identity;
  a free-then-malloc of the same base address (common for equal-size
  blocks) defeats it. Fix: a uint64 generation stamp in the FPTIA
  and FPTbias structs, bumped by set_IA_PS/set_bias_PS, keyed in
  get_FPT_IA/get_FPT_bias. Do NOT key on nuisance.random_ia: that
  redraws per sampled point and would rebuild FAST-PT (~72 ms) every
  evaluation.
- Z1/Z2/N_shear/ZCL1/ZCL2/N_CL (redshift_spline.c): build-once maps
  with no invalidation key; stale when bin COUNTS change in-process.
  Natural keys exist (redshift.random_shear, random_clustering);
  warm in init_ntomo_powerspectra beside the ggl maps.
- C_cl_tomo statics (cosmo2D.c): CLnl/CLlin/lx and the FFT block are
  sized with the current clustering_nbin but rekeyed only on
  Ntable.random; add redshift.random_clustering (the class of bug
  e7a51af fixed elsewhere).
- `test_kmax`'s chiref (redshift_spline.c): built once, cosmology
  dependent, never invalidated; currently dormant (no callers).
- init_binning_fourier/init_binning_real_space/init_probes/
  init_survey/init_bias draw no cache-busting random: a mid-process
  change does not invalidate Ntable.random-keyed tables.
- CONFIRMED AND MEASURED (2026-09-28, found while validating the
  internal ell grid): a bare post-initialize Ntable.random bump -
  init_ntable_ell_internal called with the DEFAULT value, so no
  setting changes - deterministically shifts the lsst_y1 frozen
  fiducial chi2 by 2.5e-4 (3.4317e-2 -> 3.4065e-2). Some state built
  during likelihood initialization is Ntable-coupled but rebuilds to
  different values on a bump; init_ntable_lmax and every accuracy
  init that bumps Ntable.random share the trigger. Reproduction:
  build the frozen lsst_y1 example2 model, call
  ci.init_ntable_ell_internal(nell_internal=192), evaluate. The
  planned sector-ladder cache test would have caught this class;
  finding the stale consumer is the first task of this ticket.

**Severity: MEDIUM.** Latent, notebook-reachable; each fix is one
key plus one warm-up in the established idiom.

### What is already in place

tomo.random_ggl and the init_ntomo_powerspectra warm-up show the
pattern end to end.

### What is missing

The stamps, the keys, the warm-ups, and the sector-ladder cache test
(planned per project) that exercises exactly this class.

<a id="open-fourier-nonlimber-cache"></a>
## Cache the Fourier non-Limber band-center corrections

### High-level summary

adopt_limber_gs = 0 is now the roman_fourier and roman_kl default,
and C_gs_tomo_ells has no caching: every evaluation redoes C_gs_tomo
(the FFTLog split), two 150-multipole Limber batches, and the
band-center interpolation, about +70 ms per evaluation. The
real-space path caches everything behind w_gammat_tomo's key set;
the Fourier path deserves the same.

### Current status

**Ticket type: NEW FUNCTIONALITY.**

**OPEN.** Measured 2026-09-28: roman_fourier full evaluation +26%
with the flag on (22 -> 92 ms cosmolike time).

**Severity: MEDIUM.** A per-evaluation cost on two default
configurations; the fix is the standard static-table-plus-randoms
pattern around C_gs_tomo_ells (and C_gg_tomo_ells for any Fourier
project that switches adopt_limber_gg to 0).

### What is already in place

w_gammat_tomo's cache condition lists the exact key set; the ells
functions are pure of static state, so the cache wraps cleanly.

### What is missing

The cached wrapper and a timing note in the two project READMEs.


<a id="open-cfastpt-gb2-ia-bias"></a>
## IA x higher-order-bias (gb2) cross terms in cfastpt

### High-level summary

The galaxy-intrinsic (gI) part of gamma_t is computed with linear
galaxy bias times the full NLA/TATT alignment spectrum, while the
density part of the same probe keeps the quadratic-bias terms (b2,
bs). The missing sector is the cross between the two expansions: the
spectra <delta^2|E> and <s^2|E> (times b2/2 and bs/2), which upstream
FAST-PT ships as the IA_gb2 module (tables fe, he, F2, G2, S2F2,
S2G2, S2fe, S2he). Add these tables to cfastpt and wire them into the
gI integrand, closing the one-loop consistency gap.

### Current status

**Ticket type: NEW FUNCTIONALITY.**

**OPEN.** Idea stage; gated on upstream FAST-PT resolving issue #29
(see below).

**Severity: LOW.** The truncation matches the community baseline
(DES-Y3 TATT, cosmosis tatt_interface, CCL): no current result is
wrong. The terms are a loop correction inside a correction (of order
b2*A1 relative to the linear-bias gI term); whether LSST-Y1/Roman
gamma_t precision cares is a quantifiable question to answer as part
of this ticket.

### What is already in place

`get_FPT_bias` (pt_cfastpt.c) provides the exact machinery: hardcoded
(alpha, beta, ell) J-tables accumulated through one `J_abl` call,
with exchange pairs collapsed into single doubled rows - an idiom
that is structurally immune to the transcription slip found upstream
(see the technical record). Its 13 rows were verified analytically on
2026-09-25, and the Pd1s2 rows are literally the S2F2 kernel of the
gb2 family, so part of the derivation work already exists.
`get_FPT_IA` provides the TATT (ta/tt/mix) side the new terms couple
to. Upstream FAST-PT 4.0.0 is installed in `.local` and carries the
IA_gb2 module as a cross-check reference (with caveats below).

### What is missing

- Derive the eight gb2 tables independently and write them in the
  `get_FPT_bias` collapsed-row style. Do NOT transcribe upstream:
  FAST-PT issue #28 (confirmed 2026-09-25 by direct Legendre
  derivation) has a typo in IA_gb2_S2G2 - the last row must be
  (-1,1,l=3,1/5), not l=1 - and issue #29 (a duplicated row in
  IA_gb2_he) is unresolved; the row placement is numerically inert,
  but whether the total is -1/3 or -2/3 needs the source derivation.
- Wire the new spectra into the gI integrand of gamma_t with the
  b2/2 and bs/2 coefficients paired to C1/C1delta/C2 per the gb2
  module's source paper.
- Quantify the effect on gamma_t at LSST-Y1 and Roman precision
  (delta^T C^-1 delta against the truncated model) before deciding
  whether any default changes.
- Validate against python FAST-PT once upstream has fixed #28 and
  resolved #29 (a comparison against the current upstream would
  inherit its typo).

<details><summary>Technical record</summary>

- Owners: `external_modules/code/cosmolike_core/cosmolike/
  pt_cfastpt.c` (new table block beside `get_FPT_bias`), `IA.c` /
  `cosmo2D.c` (gI integrand), plus the source-paper coefficient map.
- Upstream references: FAST-PT issues jablazek/FAST-PT#28 and #29
  (both filed 2026-09-25), fastpt/IA/IA_gb2.py at commit b91f6b7.
- Analytic facts established 2026-09-25: F2 and G2 share the
  identical (q1/q2 + q2/q1)(mu/2) term, so the (+-1,-+1) rows of the
  S2F2 and S2G2 tables must coincide; (mu/2)(mu^2 - 1/3) =
  (2/15) P1 + (1/5) P3 fixes those rows, proving #28. With
  alpha = beta = 0 and the same P(k) on both legs, J(l1,l2) equals
  J(l2,l1), which is why #29's duplicated row is numerically
  equivalent to the symmetric pair and only the TOTAL coefficient is
  in question.
- Exposure audit (2026-09-25): no gb2-family table or caller exists
  anywhere in cocoa (cfastpt, PyFAST-PT, likelihoods, notebooks);
  the cosmolike function `gb2(z, ni)` is the b2(z) galaxy bias, a
  name collision only.

</details>

<a id="open-compton-y-port"></a>
## Finish the Compton-y port (C_gy, C_ys, C_ky, C_yy)

### High-level summary

The Compton-y cross spectra were never finished when cosmolike was
ported into cocoa: nothing enables the gy/sy/ky/yy probes, and every
python binding for them was commented out. The pieces now live outside
the build in
`cosmolike_core/future_port_unfinished/cosmo2d_tmp.c` (declarations,
C functions, and the commented wrapper bindings), moved there on
2026-09-26 to clean cosmo2D.c.

### Current status

**Ticket type: NEW FUNCTIONALITY.**

**OPEN.** Idea stage; parked deliberately.

**Severity: LOW.** No analysis in cocoa uses a y probe today.

### What is already in place

The parked code compiles against the pre-batch API it was written
for, and the radial-weight infrastructure the port needs
(`radial_weights.c`) already carries the lensing kernels.

### What is missing

Port the y radial weight, revive the parked functions on the _work
batch design (they follow the retired single-grid per-point pattern),
wire probes and bindings, and validate against the original cosmolike.

<details><summary>Technical record</summary>

- Parking lot: `cosmolike_core/future_port_unfinished/cosmo2d_tmp.c`
  (NOT compiled; header comment lists the contents).
- The like struct already carries the gy/sy/ky/yy probe flags.

</details>

<a id="open-scuts-dlnxi-contraction"></a>

## Contract the dlnxi/dlnw scale-cut cache builds

### High-level summary

The per-cosmology scale-cut refill on roman_real costs ~5.7 s, and
measurement localized it to the dlnxi (and dlnw_ks) cache builds:
Ntable.dCX_dlnk_nlnk = 256 nointerp calls, each a Glpm Legendre
sweep over every integer multipole for all (theta, pair) rows. The
RF layer above is closed-form (6e160ac) and the dC tables below
cost ~44 ms, so this build is the only expensive stage left. Not
chain-critical (scale cuts run offline); worth doing when a
scale-cut study makes 5.7 s per cosmology annoying.

### Current status

Design verified against the code (2026-09-28). The nointerp per-k
cost splits three ways: (1) the dCgrid node fill - 2 x pairs x
N_ell scalar dC_ss_dlnk_tomo_limber calls, each paying cache-key
checks and a bilinear read to evaluate the table at its own ell
nodes at one fixed k; (2) limber_fill_interp expansion to every
integer ell (already AVX2, optimal); (3) the Glpm x cx simd
reduction over LMAX ~ 1e5 (already threaded and vectorized).

### What is missing

1. A perf profile of one nointerp call confirming the split (the
   RF work taught the attribution lesson twice: measure, never
   infer).
2. Step-(1) fix, standalone and small: at fixed k the node fill is
   one fixed-weight blend of two dC table k-rows (the same t for
   every entry, contiguous, SIMD-trivial); needs a _fill-style
   sharing struct for the dC statics - the C_ss/ss_ precedent.
   Expected ~1.5-2x on the build by itself.
3. The contraction, subsuming (2)+(3): both apply the same fixed
   linear map every call - limber_fill_interp's node ->
   integer-ell weights followed by the Ntable-cached Glpm dot.
   Precompute W(theta, node) = sum_ell Glpm(theta, ell)
   w_node(ell) once per Ntable change (one of today's 1e5 sweeps,
   ~200 KB); each k-node then costs 2 x pairs x Ntheta dot
   products of length N_ell. Expected order 50-200x on the build;
   validate to float reassociation (~1e-14 relative) against
   nointerp on a probe grid, nointerp kept as the guarded
   reference. The weights to contract against are exactly
   limber_fill_interp's map (verified in code, no hidden scheme).

<details><summary>Technical record</summary>

- Owner: `cosmolike_core/cosmolike/cosmo2D_scuts.c`
  (`dlnxi_dlnk_pm_tomo_nointerp`, `dlnw_ks_dlnk_tomo_nointerp`,
  their cache builders), the Glpm kernels, the la/ldx coupling.
- Timing evidence (roman_real, 2026-09-28): rf-path per-cosmology
  refill 5.676 s after 6e160ac; dC+dlnxi table stage in isolation
  42.8 ms; the difference is this build.

</details>

<a id="open-web-halo-model"></a>
## Web-halo model (WHM, arXiv:2508.10902) for nonlinear P(k)

### High-level summary

The web-halo model (Stuecker et al., arXiv:2508.10902) is a
parameter-free variant of the halo model: alongside halos it adds
structures collapsed along two axes (filaments) and one axis
(sheets), combined with 1-loop Lagrangian Perturbation Theory in
one consistent framework. Against N-body power spectra it reaches
better than 2% up to k = 0.4, 0.7 and 1.3 h/Mpc at z = 0, 0.8 and
1.5, matching baccoemu and EuclidEmulator2 across their full
w0waCDM + sum m_nu space - where HMcode-2020 needs 12 tuned
parameters and degrades at high z. The authors released WHMcode
integrated into the HMcode implementations for CAMB and CLASS.

### Current status

**Ticket type: NEW FUNCTIONALITY.**

**OPEN.** Longer-term interest declared by the maintainer
(2026-09-28); intake ticket, no implementation started. Sequencing:
after (or alongside the later phases of) the halo.c modernization
campaign, which supplies the unit-test harness and the HMcode-family
porting precedent any WHM work would build on.

**Severity: LOW.** Pure feature; no current pipeline depends on it.

### What is already in place

- The halo.c campaign staging (halo_wrapper.cpp bindings plus
  test_halo.py with frozen references) gives new halo-model code a
  validation harness.
- Original-cosmolike `theory/halo_hmcode.c` (HMCode-2020, Ludlow16
  concentrations, multi-field I0j_y) is available as a reference and
  its port is a staged campaign phase - the natural precedent for a
  WHM sector.
- EuclidEmulator2 already runs inside Cocoa, so one of the paper's
  two accuracy baselines is available in-tree for validation.
- Cocoa ships a pinned custom CAMB fork; WHMcode being released
  inside HMcode for CAMB means an upstream-adoption route exists in
  addition to a native cosmolike-side port.

### What is missing

- Route decision: adopt WHMcode through the custom CAMB fork's
  HMcode (pin bump, Fortran fences per house rules) versus a native
  C port next to the halo_hmcode phase; criteria include whether the
  halo-model sector must stay differentiable/emulable inside
  cosmolike and how the y/gas extensions interact. Maintainer lean
  (2026-09-28): a native cosmolike implementation, including the
  1-loop LPT ingredient itself, "may be advantageous for us ...
  when the time comes" - so scope the native route first and treat
  the CAMB adoption as the comparison baseline, not the target.
- Physics review of the WHM construction (sheet/filament mass
  fractions, 1l-LPT matching) before any port.
- Validation plan: reproduce the paper's accuracy claims against
  EuclidEmulator2 in-tree over the shared parameter space.

<details><summary>Technical record</summary>

- Paper: https://arxiv.org/abs/2508.10902 (WHM; public code WHMcode
  in the CAMB and CLASS HMcode implementations).
- Maintainer intake: 2026-09-28, "longer term - I am also very
  interested in implement this".
- Related staging: halo campaign documents (session scratchpad
  halo/STAGE.md, 2026-09-28) list the halo_hmcode.c inclusion phase
  this ticket would follow.

</details>

# Closed tickets

Grouped by subject and compressed. Nothing here is open work; dated
measurements are kept, since this section is a decision record, not a
README. To reopen a ticket, move its content back under
[Open tickets](#open-tickets) as a full ticket section and add its
`- OPEN` index line.

## Scale-cut refill: coarse dln caches and the 128 default (2026-09-28)

- **Implemented** (cosmolike_core 1f67c0b + one binding commit per
  project): the dlnxi/dlnw_ks caches - the measured owner of the
  scale-cut refill, one full nointerp Legendre build per ln k node -
  run their exact builds on the coarse nodes of the scale-cut k
  knob, and the house cubic spline upsamples every row onto the
  unchanged 256-node grid. The nointerp node fills read the dC
  tables' two bracketing k-rows through the new dCss_/dCks_ sharing
  structs via the file-local limber_krow_blend (SIMDe: contiguous
  loads and one fused multiply-add per lane - no gathers, since the
  weight and row offset are shared by every entry, unlike
  limber_fill_interp's per-multipole indices).
- The knob default moved 0 -> 128 from measurement (roman_real):
  refill 5.628 -> 2.788 s (2.0x; 96 gives 2.107 s = 2.7x; 192 gives
  4.317 s), max |dRF| 2.8e-3 ss / 5.9e-3 ks with ~1e-6 medians - at
  or below the response error the retired quadrature imposed; chi2
  is knob-independent, and the shipped default was verified end to
  end (2.795 s with no knob calls). The RF workers keep their
  scalar cache sampling: measured at 1.1 ms warm, nothing to win -
  the third attribution lesson of the day, this one caught before
  implementing. Suites lsst_y1 57 + roman_real 50 on the shipped
  default.

## RF workers: closed-form cumulative, quadrature retired (2026-09-28)

- **Implemented** (cosmolike_core 6e160ac): RF(kmax) integrates
  |dlnX/dlnk|, and the tabulated response the workers consume is
  piecewise linear in ln k, so both RF integrals are closed form
  (trapezoids; two triangles where the response crosses zero; a
  partial = prefix + the analytically cut last interval, direct
  index, no search). All four workers (RF_xi, RF_w_ks, RF_C_ss,
  RF_C_ks) dropped their per-kmax Gauss-Legendre sweeps for one
  cumulative per row: a per-kmax sweep was an expensive
  approximation of a function whose integral has a closed form.
- Measured: removed Gauss-Legendre error max |dRF| 3.1e-3 (ss) /
  5.4e-3 (ks), medians 8e-17 / 5.5e-5 - the old rule failed at the
  |.| kinks, where the closed form is exact, and the error was
  comparable to the coarse-table dRF, so RF accuracy was never
  better than few x 1e-3 before this change. Honest timing:
  roman_real per-cosmology rf refill 5.754 -> 5.676 s; the sweeps
  were ~80 ms of it (the open dlnxi ticket holds the true owner).
  Suites lsst_y1 57 + roman_real 50.

## Threaded scuts plane upsamples (2026-09-28)

- **Implemented** (cosmolike_core caaba62): the four coarse-grid
  refills upsample their ~120 independent (ln k, ln l) planes under
  omp parallel for (collapse(2) on the ss pair). Bitwise against
  the serial build; interleaved A/B timing (one process, 12
  alternating samples per setting): fully coarse derivative refill
  45.0 ms vs 42.8 ms exact - the upsampling overhead is gone and
  the coarse path is break-even on the already-cheap derivative
  tables.

## sigma2(M): lobe-summed quadrature, the 14.1 cutoff retired (2026-09-28)

- **Implemented** (cosmolike_core 080347a + one binding commit per
  project): the mass-variance quadrature integrates lobe by lobe
  between the zeros of the top-hat window (roots of tan x = x,
  McMahon start + Newton), one low-order Gauss-Legendre rule per
  single-bump segment, stopping when the near-geometric tail bound
  s_j r/(1 - r) drops below 1e-7 of the running total. The node
  cache (nodes, weights folded with 9 j1^2, segment offsets) is
  mass-independent and Ntable-keyed; sigma2_work fills the whole
  table threaded over masses with one p_lin read + fma per node.
  int_for_sigma2 and the fixed x < 14.1 cutoff are DELETED;
  sigma2_nointerp is the one-mass point diagnostic. The cached
  table upsamples from Ntable.N_M_internal = 192 coarse ln M nodes
  (house cubic spline in ln sigma^2); the knob joined
  init_accuracy_boost's ladder.
- Measured: the old cutoff cost up to 4.577e-3 relative in sigma^2
  (worst at M = 2e12 Msun/h; median 5.7e-4 over 1e10..1e16) - the
  cutoff sat just past the window's 4th zero (14.0662) with no
  stated error budget. Table upsampling + interpolation: max
  2.7e-5. Stage-2 point diagnostic bitwise-reproduces the stage-1
  lobe values; OMP 1 vs 4 bitwise identical; suites lsst_y1 57 +
  roman_real 50.

- **Implemented** (cosmolike_core 0172c0e + one binding commit per
  project): C_gk_tomo_limber and C_ks_tomo_limber adopt the internal
  coarse ell grid (exact quadrature on Ntable.N_ell_internal nodes,
  the house cubic spline upsampled onto the unchanged N_ell table,
  workspace static in the Ntable rebuild block, threaded
  evaluation); C_gg alone keeps the exact grid. The four dC_X/dlnk
  scale-cut tables coarsen in 2D through the new
  spline2d_upsample_uniform (basics.c: tensor-product natural
  bicubic, two passes of the house 1D spline, precomputed direct
  indices): the ell axis follows N_ell_internal, the ln k axis gets
  Ntable.dCX_dlnk_nlnk_internal (default 0 = exact; ln k carries
  the BAO wiggles). init_accuracy_boost became the global super
  function - one call scales every sampling knob, internal grids
  and FPT_internal_accuracy_boost included, from first-call
  baselines, so the coarse/dense ratios are boost-invariant -
  recorded in the cosmolike-dev skill with the rule that every new
  sampling knob joins its ladder.
- Validated: lsst_y1 frozen chi2 3.4316970976e-02 bitwise (3x2pt
  untouched by gk/ks); desy1xplanck 6x2pt exact-vs-192 dv max rel
  1.159e-05 (median 2.2e-07), chi2 delta 4.3e-08; suites lsst_y1
  57, roman_real 50, desy1xplanck 45, all green. Scale-cut impact
  measured through RF, the real consumer (the dlnC batch reads
  compute on caller grids and bypass the cached tables): ell at
  192, max |dRF| 3.2e-3 ss / 4.3e-3 ks (medians 2e-16 / 8e-5);
  k at 192-of-256 adds a similar 3.2e-3 / 3.5e-3 - RF integrates
  |dlnC/dlnk| over ln k, so wiggle interpolation errors partly
  cancel. Also fixed in passing: a factor-2 error in the
  spline_coeffs_uniform header derivation (the code was always
  right).
- Ticket-B evidence from the A/B: on desy1xplanck the fresh chi2 at
  internal=192 depends on how many pre-first-evaluation
  Ntable.random bumps ran (8.6002538251e-06 after one init call,
  8.5994603071e-06 after two; the exact path is bump-invariant in
  the same comparison) - the invalidation-sequence dependence the
  cache-hardening ticket hunts.
- Timing post-mortem (three protocols, roman_real): on
  Ntable.random bumps the scuts rebuild is ~5.9 s of Ntable-keyed
  machinery regardless of the knobs; on cosmology-keyed refills it
  is ~5.7 s dominated by the dlnxi cache build (its open ticket
  holds the contraction design); the dC+dlnxi table stage in
  isolation is 42.8 ms exact vs 154 ms coarse-serial vs 45.0 ms
  coarse-threaded (caaba62). Lesson recorded: time the
  invalidation class the workflow actually triggers, and profile
  before attributing.

## Sector-ladder cache-consistency test (2026-09-28)

- **Implemented** (one commit per project, all six):
  `tests/test_cache_consistency.py` walks a deterministic parameter
  ladder in one process - three cosmology-only steps, then IA-only,
  source-photo-z, lens-photo-z and shear-calibration steps, a sector
  the configuration does not sample dropping out (the Roman lenses
  are the source sample; roman_kl also fixes every IA amplitude) -
  evaluating after every step with cobaya's cache bypassed, then
  scrambles every sector at once and returns to the ladder's final
  point: the data vector and chi2 must reproduce bit for bit, and a
  second instance walking the mirrored sector order must land on the
  same vector. Every step must change the vector, each
  shear-calibration step must equal the analytic (1+m_i)(1+m_j)
  block rescale to 1e-12 (the ks cross scales by its source factor;
  the point-mass term rides inside gamma_t), and a no-op update must
  change nothing. Both IA models run: the TATT ladder exercises the
  FAST-PT rebuild machinery NLA never touches.
- All six projects pass (NLA + TATT, forward + mirrored). The target
  is the partial-invalidation class the per-point suites never
  exercise - the class of the cache-hardening ticket's measured
  post-initialize Ntable.random bump defect.

## Internal coarse ell grid for the C_ss/C_gs tables (2026-09-28)

- **Implemented** (cosmolike_core a3d19e2 + one binding commit per
  project): the exact quadrature runs on Ntable.N_ell_internal
  log-spaced nodes (default 192) and a cubic spline in ln l upsamples
  each (component, pair) row onto the unchanged 512-node table at
  cache-build time; init_ntable_ell_internal (bound in all six
  projects, 0 = exact) is the A/B switch. C_gg keeps the exact grid
  at every node, with the BAO-wiggle warning in its header.
- Validated (lsst_y1 frozen fiducial, equal footing): data-vector
  upsampling error max 2.1e-4 relative (median 5.6e-7), chi2 delta
  1.9e-6. Timing (median of 8, 4 threads; the deltas are CAMB-free -
  the theory cost is identical in both settings and cancels):
  lsst_y1 saves ~3 ms per evaluation, roman_real ~19.5 ms
  (1.5037 -> 1.4842 s; ~13% of its cosmolike share), roman_kl
  ~68 ms (2.2294 -> 2.1615 s; its 10 source bins carry the largest
  tables) - the pair-count scaling the ticket predicted. Suites
  lsst_y1 55, roman_real 48 after the six projects
  re-pinned their baryon drift vectors on this build (the
  test_accuracy_baryons in-place rewrite, now documented in the
  maintenance skill, surfaced with the first shear-touching change).
- The A/B also exposed and recorded a pre-existing defect in the
  cache-hardening ticket: a bare post-initialize Ntable.random bump
  shifts the frozen fiducial chi2 by 2.5e-4.
- Follow-up after review (cosmolike_core e1e33f8): the coarse-grid
  workspace - nodes, prefactors, result and coefficient tables, and
  the precomputed fine-node interval indices and offsets - is
  allocated only in the Ntable.random-keyed rebuild block, and the
  GSL spline objects were replaced by the house natural cubic spline
  (spline_coeffs_uniform + direct-index Horner evaluation; both
  grids are uniform in ln l, so there is no search and no
  accelerator state), with the evaluation loop threaded. The lsst_y1
  frozen fiducial chi2 3.4316970976e-02 reproduces exactly.

## cfftlog activity mask: empty radial slots skipped (2026-09-28)

- **Implemented** (cosmolike_core da14c54): cfftlog_ells_p1/_p2 accept a
  per-(row, slot) activity mask (NULL = all active); inactive slots
  skip their forward and inverse FFTs and write exact zeros.
  C_gs_tomo declares the mask its fx assembly implies (lens rows:
  density always, RSD under include_RSD_GS, magnification when some
  gbmag != 0; source rows: slot 2 alone) - 10 of 30 row-slots active
  at the lsst_y1 defaults. C_cl_tomo keeps NULL: its conditional
  SIZE2 already drops the magnification slot when bmag = 0. The
  relevance grew with the 2026-09-28 default switch: non-Limber ggl
  now runs by default in four projects.
- Validated: bitwise (the skipped slots held identically zero
  kernels; lsst_y1 frozen references reproduce to all ten printed
  digits); full-evaluation median 1.677 -> 1.673 s (macOS arm64,
  4 threads, median of 8). Suites: lsst_y1 55, roman_real 48.

## Correctness sweep of the 2026-09-28 review findings (2026-09-28)

- **All confirmed defects fixed in one pass** (cosmolike_core
  c4aa392): the sigma2_nointerp out-of-bounds ar[1] read (the scale
  factor now reaches p_lin), the chi_all fallback j+2 clamp, finite
  f_growth/norm_growfac_all at exactly z = 0 (all four variants,
  z > 0 branch bit-identical), the swapped f_K open/closed labels,
  the swapped zmean_source cache stamps (table rebuilt every call),
  the N_ggl < -1 rebuild guard (an excluded pair (0, 0) rebuilt the
  map per call), the W_RSD bin guard, the stale bs2 on b1-only
  updates, the set_source_sample bump-before-warm-up order, the
  init_probes mixed-case out_of_range, the init_ggl_exclude
  odd-length refusal, the dlnxi normalization swap (1 - p; reported
  2026-09-27), and the RF_C_ss/RF_C_ks vanishing-denominator guard
  (the BB plane under NLA was 0/0 = NaN). Cleanups: misattributed
  fname strings, wrong wrapper error names, dead has_b2 copy and the
  disabled flat-vector block, double semicolons, the gg scalar
  wrapper overload declared.
- Validated: frozen lsst_y1 references reproduce the pivot-build
  values at all compared digits; test_scale_cut_diagnostics 4/4
  (the dlnxi and RF changes live under its assertions); suites
  lsst_y1 55, roman_real 48, green.

## Per-bin pivot for the FKEM separable spectrum (2026-09-28)

- **Implemented the same day it was proposed** (cosmolike_core
  3c9ed69): the separable linear spectrum of the non-Limber split
  is anchored per lens bin at a_piv = 1/(1 + zmean(bin)),
  P_sep = (D(z)/D(a_piv))^2 P_lin(k, a_piv), in the FFTLog terms of
  C_cl_tomo and C_gs_tomo and both use_linear_ps = 1 Limber legs
  together (the cancellation invariant); C_gs_tomo's shared P_lin
  table became per lens bin. COSMO2D_FKEM_PIVOT_Z0 restores the
  z = 0 anchor exactly.
- Validated: the l = 149 cancellation floor unchanged (des_y3 BMAG
  fiducials, 1.126e-3 vs 1.125e-3); lsst_y1 references move by
  delta chi2 = 0.034 (3x2pt) / 0.031 (2x2pt), inside the 0.2 band,
  deliberately not refrozen; at the fiducial Sigma m_nu the
  correction content shifts by 0.05-0.15% of itself at l <= 50 (the
  removed error grows with Sigma m_nu). Suites on the pivot build:
  lsst_y1 55, roman_real 48, all green.

## Non-Limber defaults, regenerated data vectors, and the refreeze (2026-09-28)

- **ggl defaults switched to non-Limber** (`adopt_limber_gs: 0` in
  the yamls, interface bindings, prototype fallbacks, tests and
  READMEs) for lsst_y1, roman_real, roman_fourier and roman_kl: the
  measured Limber cost, delta^T C^-1 delta = 1.86 / 0.49 / 1.30 /
  1.63 at the frozen fiducials, was judged too large to absorb.
  des_y3 and desy1xplanck keep Limber ggl (0.011 / 0.0037).
- **Shipped simulated data vectors regenerated** at the frozen
  fiducial points with the gg growth fix and the new defaults:
  lsst_y1_theory.modelvector (NLA) and lsst_y1_theory_TATT
  .modelvector, roman_real example1.modelvector, roman_fourier
  roman_example.datavector, roman_kl roman_kl_3x2.modelvector and
  roman_kl_2.modelvector (KL-mode vectors: the ggl switch moves
  every 3x2 mode and the probe-mixed third of the mcmc vector, by
  construction of the compression).
- **All six snapshots refrozen** (`--overwrite`, `--baryons`): every
  exact-configuration reference sits at its minimum again
  (0.000000-level; the des_y3 0.168 and desy1xplanck 0.0092 gg-fix
  drifts absorbed into the regenerated frozen fiducial vectors; the
  measured drifts matched the predictions to three digits). lsst_y1
  emul2 references now read the pure emulator-vs-exact difference
  (0.289 / 0.772, was 0.554 / 1.383 against the stale vector).
- **Freeze-generator gap closed**: the per-mask TATT dataset
  descriptors of the `--mask` comparison sweeps (hand-added in
  lsst_y1 712739e and siblings) were not reproduced by
  `generate_frozen_reference.py`, so any refreeze silently dropped
  them. All six generators now carry `TATT_MASK_VARIANTS` +
  `generate_tatt_mask_datasets()` (run inside `--overwrite`, plus
  the incremental `--tatt-masks` mode).

## int_for retirement and the C_gk _work port (2026-09-28)

- **int_for_C_gs/gg/gk_tomo_limber retired** with their GSL
  per-(node, ell) quadratures; every probe keeps one single-ell
  point diagnostic in the C_ss_tomo_limber_nointerp pattern (one
  batch call, one entry read), and C_ks gained the missing one.
  int_for_C_kk_limber stays: kk has no tomography. Validated by a
  rebuilt lsst_y1 reproducing the frozen reference chi2 to every
  printed digit.
- **C_gk moved onto the _work batch API** (C_gk_tomo_limber_work +
  _nointerp_ells/_batch in the gg design with one leg W_k, the
  64-node ladder of the retired scalar kept, RSD gate and one-loop
  FPTbias tables as in gg): table builder, w_gk low-ell loop and
  both wrapper overloads migrated. Validated on the desy1xplanck
  6x2pt fiducial: ss/ggl/gg/ks/kk bitwise unchanged, the 37 changed
  w_gk entries within 4.0e-15 relative, chi2 unchanged at ten
  digits. Closes the gk/ks _work ticket.

## Non-Limber galaxy-galaxy lensing and the gg batch API (2026-09-27)

- **Non-Limber C_gs shipped** (cosmolike_core): `C_gs_tomo`, the
  FKEM split (arXiv:1911.11947) in the house FFTLog design (hoisted
  forward FFT, per-pair early exit, Limber continuation), with lens
  rows on the Limber lens range (magnification foreground kept),
  one combined lensing + NLA source kernel per source bin, and the
  subtracted Limber term in separable growth. Both Limber terms come
  from the batched `C_gs_tomo_limber_linpsopt_nointerp_ells`
  (`use_linear_ps` switch in `C_gs_tomo_limber_work`, bitwise neutral
  at 0). `C_gs_tomo_ells` adds the correction to Fourier band centers
  (interpolated between integer l). `like.adopt_limber_gs` (default
  1) with `init_adopt_limber_gs`, the yaml key in every ggl
  likelihood of the six projects, `w_gammat_tomo` keyed on the flag,
  and one `tests/test_nonlimber_ggl.py` per project. The reference
  code (cluster_chto `cosmo2D_exact.c`) was not ported line by line:
  its gamma_t Legendre kernel has a typo (Pmin[l+1]-Pmax[l] for
  Pmax[l+1]; wrong by orders of magnitude at l = 2), its l = 1 terms
  are uninitialized, its early exit is dead code, and its RSD uses the
  gamma approximation f = Omega_m(z)^0.55.
- **Validation (lsst_y1 3x2pt fiducial).** FFTLog term vs a
  brute-force double integral (scipy spherical_jn, 20000-point chi
  grid, same kernels): 1e-6 to 1e-4 relative for l >= 3 in all 25
  pairs (~1% at l = 2 for the highest lens bin). FFTLog -> Limber at
  l = 149: ~1e-4 for lens-in-front pairs. Early exit vs none within
  the 1% tolerance. Default data vectors bitwise equal to HEAD (full
  precision). Debug build (UBSan, scalar fallbacks) and aggressive
  build agree with the default to 8e-13 and 6e-12 in lens-bin units;
  OMP 1 vs 8 bitwise; repeated calls bitwise.
- **gg batch API.** `C_gg_tomo_limber_work` + `_linpsopt_nointerp_ells`
  / `_nointerp_ells` / `_nointerp_batch` in the ss/gs design, same
  128-point rule as the scalar path; the table, the w_gg low-l loop,
  the C_cl_tomo Limber pair, the Fourier gg data vector and the
  notebook wrapper moved onto it. Batch vs scalar: 2.2e-15 relative;
  w(theta) moves by 3.3e-14. Cost: 33 ms vs 45 ms (1750 l x 5 bins,
  one thread; W_RSD and P_delta dominate). The measurement of the gs
  batch Limber pair: 3.3 ms vs 51 ms scalar (25 pairs x 150 l).
- **Aggressive-mode n(z) fix.** Under COSMOLIKE_AGGRESSIVE_MODE the
  `arma::find(nofz > c*nofz.max())` of set_lens/source_sample
  returned an empty list for a valid column once this work changed the
  code generation of generic_interface.cpp (HEAD happened to escape).
  The bin-support search is two plain loops, with identical results
  in the default build.
- **Non-Limber C_gg growth fix (kept 2026-09-28, references not
  refrozen by the maintainer's decision).** `C_cl_tomo` subtracted
  C_limber(P_lin) with CAMB's p_lin(k, a) while its FFTLog term uses
  the separable D(z1) D(z2) P_lin(k, 0); D(a)^2 differs from
  P_lin(k,a)/P_lin(k,0) by 0.7% (z = 0.3) to 1.6% (z = 1), so the pair
  never cancelled at high l (-0.85% to -1.6% at l = 149; the 1% early
  exit fired near l = 50 on the crossing of the decaying non-Limber
  excess with the offset). The subtracted term uses D(a)^2
  P_lin(k, 0) (batch and scalar integrands); mismatch at l = 149:
  0.02% to 0.17%. w(theta) delta^T C^-1 delta: lsst_y1 1.52, des_y3
  0.168 (Y3) / 0.030 (Y1), desy1xplanck 0.0092, roman_real 0.0083,
  Fourier projects 0. The lsst_y1 3x2pt/2x2pt reference checks fail
  by 1.43 to 1.52 until a refreeze; the other projects stay within
  the 0.2 limit. Each real-space project's tests/README.md explains
  it in an appendix FAQ.
- **gg Limber switch as a yaml key (2026-09-28).** `adopt_limber_gg`
  (`init_adopt_limber_gg`; 0 = non-Limber in the real-space projects,
  1 = Limber in roman_fourier and roman_kl, so no default changed),
  w_gg_tomo keyed on the flag, and `C_gg_tomo_ells` for Fourier band
  centers. `tests/test_nonlimber_gg.py` in every project measures
  non-Limber minus Limber clustering, delta^T C^-1 delta: lsst_y1 148,
  roman_kl 57.4, des_y3 9.26, roman_real 8.58, desy1xplanck 6.52,
  roman_fourier 3.34 (against 0.0037-1.86 for ggl).
- **Two pre-existing bugs found by the full test suites.** (1) The
  python FAST-PT setters `set_IA_PS`/`set_bias_PS` replaced
  `FPTIA/FPTbias.tab` but left `tab_int` aliasing the freed table;
  the next cfastpt rebuild freed it again (malloc abort in
  `get_FPT_IA`, lsst_y1 suite, also in crash reports of 2026-09-26).
  The setters re-alias `tab_int`, and `get_FPT_IA`/`get_FPT_bias`
  rebuild whenever `tab` is not the table they allocated
  (cfastpt -> FAST-PT -> cfastpt reproduces the first chi2 exactly).
  (2) The static ggl pair maps (`test_zoverlap`, `ZL`, `ZS`, `N_ggl`)
  were filled once per process, so a cosmic-shear model built after a
  3x2pt model with another `ggl_exclude` list miscounted the pairs
  (roman_kl: "IP::set_mask: inconsistent mask"). They rebuild on
  `tomo.random_ggl`, refreshed by `init_ntomo_powerspectra` and
  `init_ggl_exclude`, which also warms them single-threaded.
- **Suites (2026-09-28, one project at a time):** roman_real 47
  passed, roman_fourier 42, roman_kl 46, lsst_y1 48 passed and 6
  failed (the frozen 3x2pt/2x2pt references, off by the gg growth
  fix's delta chi2 of 1.43-1.52); des_y3 60 and desy1xplanck 42 on
  the build before the two bug fixes.

## Scale-cut diagnostics for C_ks and w_ks (2026-09-27)

- **The ks family shipped** (cosmolike_core 7ded4a7; desy1xplanck
  694bd8d): `dC_ks_dlnk_tomo_limber_work` (single-node grid
  derivative, per-bin source-support gating, fused normalize mode),
  cached dC/dlnC tables, `RF_C_ks`, the real-space `dlnw_ks_dlnk`
  pipeline (limber_fill_interp gather, CMB beam/pixel filter,
  gamma_t-type Legendre kernel, normalized by w_ks) and `RF_w_ks`,
  wrapper overloads at the doc standard, desy1xplanck bindings (plus
  plain C_ks_tomo_limber and a fixed live w_ks_tomo binding), six
  notebook wrappers, and a 15-cell EXAMPLE_EVALUATE1 section.
- **Validation.** Truth check: the dC_ks/dlnk integral reconstructs
  C_ks to 4.7e-5 relative (integration_accuracy 3; the 1.9e-3 at the
  default is the 64-point C_ks rule itself). Scalar-vs-array
  overloads and repeated calls bit-identical; RF bounded and
  saturating at 1; 6x2pt chi2 unchanged to the last digit; full
  41-test suite green; the notebook executes end to end (98 cells,
  33 figures). The truth check also exposed the inherited ss
  amplitude bug (see the 2026-09-27 bullet under the scuts
  refactor).

## Runtime n(z) photo-z conventions (2026-09)

- **Two runtime flags shipped (2026-09-24/25).** `Ntable.photoz_interpolation_type`
  (0 = cspline default, 1 = linear, 2+ = Steffen monotone; the
  documented-but-unimplemented switch in basics.c became real, inside
  `malloc_gsl_interp`/`malloc_gsl_spline`, whose only callers are the
  photo-z readers) and `Ntable.photoz_zmid_convention` (0 = z column
  read as Z_LOW left bin edges, values at centers z + dz/2, the
  historical behavior; 1 = Z_MID sample points). One shared
  `init_photoz_conventions` in generic_interface.cpp; the n(z) caches
  watch both values through a packed slot, so a runtime flip rebuilds
  the tables (both directions covered by a bit-identical round-trip
  assertion). No preprocessor flag: a compile-time `#ifdef` was
  rejected because the unit test would need two builds. Defaults
  unchanged everywhere; every default chi2 reproduces its frozen
  reference exactly.
- **Wired end to end in all six projects.** pybind binding, the
  likelihood call, yaml keys in every likelihood variant (31 yamls),
  the notebook wrappers of lsst_y1/roman_real (_CONFIG +
  configure()), and the inline des_y3 notebooks.
- **Unit test + figures + README section in all six projects.**
  `test_photoz_conventions.py` evaluates five settings in one process
  and measures each alternative as delta^T C^-1 delta against the
  default (masked inverse covariance; dead-flag floors);
  `generate_photoz_convention_figure.py` makes the per-pair
  delta-xi/delta-Cl figures the tests README records.
- **Measured (2026-09-24/25, frozen cosmic-shear fiducials).**
  Steffen: 4.1e-5 (lsst_y1), 1.3e-6 (roman_real), 2.9e-6
  (roman_fourier), 1.6e-5 (des_y3), 5.7e-5 (desy1xplanck), 1.5e-6
  (roman_kl); linear: 6.3e-4, 1.6e-4, 9.4e-4, 3.1e-4, 1.3e-3,
  2.0e-5; Z_MID: 1.63, 0.94, 3.60, 0.30, 0.42, 2.35. The interpolant
  is far below statistical precision; the half-bin z-column reading
  is the one photo-z convention that matters, as the 2026-09 DES-Y6
  three-code comparison predicted. Full suites green in all six
  projects (2026-09-25): lsst_y1 (49), roman_real (42), roman_fourier
  (41), des_y3 (59), desy1xplanck (41), roman_kl (45).
- **Deferred decisions, deliberately.** Flipping the default to
  Steffen (the measurements support it whenever desired; a deliberate
  documented accuracy change under the full validation protocol when
  taken). The correct Z_LOW-vs-Z_MID setting per survey is a data-
  product fact (which column the tables were exported from), declared
  per analysis in its likelihood yamls. Motivating measurements: the
  DESY6 source tables carry an overflow-like last row against which
  cspline rings at -0.3% to -0.6% of peak; a single-cell spike makes
  cspline undershoot -13.7% of peak while Steffen stays non-negative;
  detector-plus-rebuild designs were rejected because sampled-n(z)
  analyses rebuild the tables per likelihood call.

## FAST-PT two-grid upgrade (2026-09)

- **Two-grid theory block.** The wrapper now computes its FFTLog
  convolutions on a coarse internal grid and upsamples with a cubic
  spline in log k onto the dense output table the likelihood reads
  linearly. `accuracyboost` scales the output table,
  `internal_accuracyboost` the internal grid, and both are rebased so
  1.0 is the converged default; values above 8 are refused with the
  reason. Shipped as PyFAST-PT v1.0; README added in v1.01.
- **Divergence located.** The historical CFASTPT vs FASTPT
  disagreement (up to Δχ² = 21 at LSST-Y1 precision and Δχ² = 87184
  at Roman precision on the shared single 1,100-point grid) lived
  entirely in the density of the interpolated output table and none
  of it in the convolutions.
- **Six projects converged.** The comparison test passes at the
  defaults on lsst_y1, roman_real, desy1xplanck, des_y3, and roman_kl
  (which keeps a density-independent floor ≈ 0.175); roman_fourier
  needs `accuracyboost: 2.0`. Each `tests/README.md` carries the
  contract, the measurement table, and the corner plot; the example
  yamls carry the rebased defaults.

## Batched shear wrappers and the cosmo2D_scuts refactor (2026-09-26)

- **cosmo2D_wrapper on the batch APIs (2026-09-25).** The C_ss and
  C_gs array overloads make one `_nointerp_ells` call and scatter;
  the scalar overloads are one-element calls into the same engines,
  documented as point diagnostics. Measured on lsst_y1 (100
  multipoles, 4 threads): C_ss 13.2 -> 0.5 ms, C_gs 10.8 -> 0.9 ms
  per call. C_ss matches the retired loop at reassociation level
  (<= 4e-15); C_gs moved onto the likelihood's own 64-point rule
  (the legacy scalar used 96 points), each verified at its rule's
  error against a 2000-point reference.
- **cosmo2D_scuts on the _work design (2026-09-26).**
  `dC_ss_dlnk_tomo_limber_work` computes the scale-cut derivative on
  a (ln k, ell) grid (each output one core evaluation at the Limber
  node chi = (l+1/2)/k, amplitude dchida/fK); its normalize flag
  computes C_ss on the same thread team and divides cache-hot rows to
  dlnC. The hidden deriv argument of int_for_C_ss_tomo_limber is
  gone (the dC point value is the C_ell integrand times chi). The
  cached tables start at l = 1, killing the exact-scalar low-l
  branches; `RF_C_ss/RF_xi_tomo_limber_work` run both RF integrals on
  precomputed Gauss-Legendre node arrays with the kmax-independent
  denominators computed once. RF array call (26x26): 108 -> 25 ms.
- **dlnxi redesigned (2026-09-26).** dlnxi_dlnk_pm_tomo_nointerp now
  reads the dC table on its multipole log-grid, gathers to every
  integer multipole with the vectorized limber_fill_interp (exported
  through cosmo2D.h), and Legendre-sums with restrict + SIMD: first
  call 4.04 s -> 0.04 s, steady 15 -> 8 ms, and the 6.1 GB
  per-integer-multipole cache is gone (work arrays are
  Ntable/tomography-keyed statics).
- **Wrapper cleanup (2026-09-26).** Every int_for_* python wrapper
  and binding deleted from cosmo2D_wrapper and all six projects'
  interface.cpp (never used from notebooks, incompatible with the
  _work design); the unfinished Compton-y family parked in
  future_port_unfinished/cosmo2d_tmp.c; per-argument documentation
  across cosmo2D_scuts, cosmo2D_wrapper and cosmo2D_scuts_wrapper.
- **Validation.** Every step compared old-vs-new dumps on the lsst_y1
  frozen TATT fiducial: bit-identical wherever the arithmetic was
  unchanged (the single-work fusion, the statics, y-move and int_for
  removal were all bitwise), reassociation-level (<= 1e-13) for the
  reorganized sums, and lnl-respacing-level shifts only where the
  table grids deliberately changed. Full suites green through the
  arc; `test_scale_cut_diagnostics.py` added to lsst_y1 and
  roman_real. A cosmetic leftover: the notebook wrappers call the
  dlnC function `dlnC_dlss_tomo_limber` (scrambled name), a
  rename-commit candidate.
- **End-to-end rf_xi measured (2026-09-27).** roman_real, 76-value
  kmax grid, idle machine: the retired path (GSL RF integrals over
  the per-integer-multipole dlnxi table) took 22.4 s on the first
  call and 19.7 s on every repeat at the same cosmology — the GSL
  quadrature re-drove the scalar integrand on each call, cached
  tables or not, while the one-time table builds cost only the
  remaining 2.7 s (and 10.3 GB at roman lmax). The _work path takes
  4.2 s, first call and repeat alike (5.3x / 4.7x), on ~40 MB of
  statics.
- **dC/dlnk node amplitude fixed (2026-09-27).** The per-node
  amplitude multiplied the per-a integrand by |dchi/dlnk| = chi
  without converting the measure (|da/dlnk| = fK/dchida), leaving a
  spurious dchida — inherited verbatim from the retired GSL
  implementation and preserved by the old-vs-new migration, which
  never ran a truth-integral check on ss. Caught by the ks family's
  reconstruction test and fixed (cosmolike_core 81319fc): the
  integral of dlnC_ss/dlnk over ln k is now 0.9996-1.0001 (it was
  1.25-1.57, the kernel-weighted mean of dchida). Every ss
  diagnostic (dC, dlnC, RF_C_ss, dlnxi, RF_xi) carried the smooth
  k-dependent tilt; the likelihood never read this path, and chi2
  is bit-identical.

## Notebook derivative crashes (2026-09-26)

- **Root cause found: a deterministic exit.** The retired low-l RF
  branch underflowed k = exp(-(1-t)/t) to 0 at small quadrature t,
  and the dC point diagnostic's k > 0 guard was log_fatal -> exit(1):
  any rf_C_ss call at l <= 20 killed the process — and a jupyter
  kernel dies from exit(1) with no traceback, matching the reported
  symptom. Two aggravators: every log_fatal/critical in the
  diagnostics exits the kernel the same way (for example running
  cells out of order), and the old dlnxi build allocated 6.1 GB and
  took 4 s per cosmology.
- **Fixed by the scuts refactor** (the branch and the point function
  are gone; dlnxi's cache is gone), and **verified with the
  notebooks' own calls**: drivers running the exact derivative cells
  of EXAMPLE_EVALUATE1 (roman_real, and the section newly ported to
  lsst_y1) complete with finite outputs; the regression test
  `test_scale_cut_diagnostics.py` pins the old fatal multipoles
  (rf_C_ss at l = 3) in both suites.

## Two-grid C-FAST-PT internal grid (2026-09-25)

- **Runtime knob shipped.** `Ntable.FPT_internal_accuracy_boost`
  (yaml key `internal_accuracyboost`; setter `init_fpt_internal_boost`,
  which refuses values <= 0) runs the C-FAST-PT FFTLog convolutions
  of `get_FPT_IA` and `get_FPT_bias` on an internal grid of
  ceil(N * boost) points rounded up to even (the FFTLog engine
  requires an even count — discovered by the convergence scan) and
  moves the spectra onto the unchanged output table with
  `spline_coeffs_uniform` plus Horner evaluation on the uniform ln k
  grid (`fpt_regrid`, no binary search); 1.0 recovers the single-grid
  path exactly and is the reference arm of the accuracy tests. Work
  tables (`FPT.tab_int`), regrid scratch and the FFTLog
  configurations persist across cosmologies and rebuild only when
  Ntable changes; the FFTLog padding/extrapolation spans scale with
  the grid ratio (they are ln-k spans, not point counts).
- **Default 0.5 (Nint = 550), chosen with margin.** Measured
  2026-09-25 on lsst_y1 frozen TATT fiducials: delta^T C^-1 delta vs
  the single-grid path <= 1e-10 (shear) and <= 1e-9 (3x2pt) down to
  298 internal points, converged from above at 2200; the default
  keeps a ~100x chi2 margin over the 298-point arm. FPT table cost
  fits 1.55 us * N ln N on the M-series Mac: 11.9 -> 5.4 ms per
  evaluation (6.5 ms saved, CAMB excluded); the TATT premium over
  NLA drops from ~25 ms to ~7 ms.
- **Validated in all six projects (2026-09-25).**
  `internal_accuracyboost` is in every project's high-accuracy
  settings and one-at-a-time accuracy knob scan; full suites green:
  lsst_y1 (49), roman_real (42), roman_fourier (41), des_y3 (59),
  desy1xplanck (41), roman_kl (45). Parity at the default point is
  exact (lsst_y1 chi2 0.267263), and the convergence scan stayed
  bit-identical through every review refactor.
- **Deferred decisions, deliberately.** `perf stat -r 3` on the
  Linux roman benchmark (the citable numbers) and the
  DEBUG/AGGRESSIVE build-mode sweep are maintainer steps.

## Shared test harness (2026-09)

- **cocoa_testing.py.** About ten thousand duplicated lines of unit
  test python across five projects consolidated into
  `external_modules/code/cosmolike_core/cocoa_testing.py`
  (`CocoaTestHarness`: frozen-state integrity, chi2 pipeline, baryon
  checks, worker isolation, race checks, CFASTPT vs FASTPT
  comparison). The projects became data-only shims binding one
  harness instance; the full suites reproduced their digits exactly
  before and after. Shipped in cosmolike_core v4.11.6.

## CFASTPT vs FASTPT: 3x2pt, 2x2pt, and scale-cut masks (2026-09-23)

- **Comparison extended to every probe.** All six projects gained
  the 3x2pt (desy1xplanck: 6x2pt) and 2x2pt sweeps next to the
  cosmic-shear one, passing the 0.2 rule at their contracts; each
  `tests/README.md` carries the dated measurements. The frozen
  configurations fix the one-loop bias amplitudes at zero, so the
  sweeps score the IA tables; lsst_y1's exploratory b2-activated
  variant measured the bias tables separately (max 0.044,
  density-independent floor 0.013).
- **--mask option.** The sweeps rerun under a chosen scale-cut mask
  (frozen contract, lsst_y1's M2-M6, or the all-ones no-cuts mask)
  via frozen TATT dataset variants; under a non-frozen mask the
  chi2 baseline is regenerated from the cfastpt fiducial (zero by
  construction). conftest.py content was deduplicated into
  cocoa_testing.py (the file stays per project for pytest
  discovery; a shim binds the shared hooks).
- **Findings.** The ones-mask growth at default camb/cosmolike
  settings collapses under the pushed settings (cosmolike
  integration accuracy, not the FAST-PT grids); three shipped
  covariances (desy1xplanck 6x2pt, roman_real 3x2pt, roman_fourier
  3x2pt) are not positive definite when fully unmasked; roman_kl's
  frozen 3x2pt mask is already the all-ones mask.
- **Stray-file repairs.** roman_fourier's ones.mask was a
  wrong-length lsst_y1 byte-copy, replaced by the correct
  1,485-row mask (one guarded frozen replacement); roman_fourier's
  and roman_kl's real-space calculate_mask.py strays were replaced
  by Fourier generators that reproduce every shipped mask
  byte-identically.

## Nonlinear P(k): Halofit vs EE2, and the EE2 modifications (2026-09-23)

- **Halofit vs EE2 (advisory checks NL1/NL2, all six projects).**
  Ten shared seeded cosmologies in omegam/ns/As, the EE2 vector as
  each cosmology's fiducial, the difference weighted under the
  --mask scale cuts, reported through the per-project corner
  figures. The two sources are not interchangeable at survey
  precision under the frozen cuts anywhere (medians from 1.4 at
  DES-Y3 shear to 912 at roman_kl 3x2pt, always worst at high
  omegam).
- **Cocoa's EE2 vs the original EE2.** The original (commit
  ff59f66) cannot run inside Cocoa as-is: no get_boost2, and a
  silent overflow beyond 101 redshifts (the likelihoods send ~110)
  - two of the defects the Cocoa modifications fixed. Bridged with
  a compatibility patch (a get_boost2 adapter plus 100-redshift
  chunking), the side-by-side ten-cosmology comparison passes on
  all six projects (max delta chi2 4.3e-4, at roman_fourier).
  Shipped as lsst_y1's test 18, which compiles the original at
  test time (offline, --ignore-installed) and gates at 0.2; the
  modifications, their measured 14x speed-up, and the validation
  table are documented in the euclidemu2 repository README.
- **EE2 OpenMP race tests.** Every project ships a race check with
  the nonlinear P(k) from EE2 (ten_in_a_row_chi2's ee2 flag), all
  passing with exact eight-decimal agreement.

## Installation scripts (2026-09)

- **SWITCH_TO_DEV_MODE for setup_pyfastpt.sh.** The wrapper clone now
  routes through the per-script `devurl()` helper like every other
  CosmoLike-org repository; `set_installation_options.sh` keeps the
  https address.
- **Tag checkout refname fix.** The wrapper's tag checkout uses the
  `-b "${TAG}TMP"` form the other setup scripts use, removing the
  ambiguous-refname failure when a branch and a tag share a name.

## Release train (2026-09-23)

- **Nine repositories released** and pinned in
  `set_installation_options.sh`: cocoa v4.11.7, cosmolike_core
  v4.11.6, lsst_y1 v4.11.2, roman_real v4.11.2, desy1xplanck v4.11.1,
  des_y3 v4.11.1, roman_fourier v4.11.1, roman_kl v4.11.4, PyFAST-PT
  v1.01, BFMT v1.01. Every pin verified against an existing remote
  tag; every repository cycled onto a fresh `bugfix` from its updated
  `main`.

## New projects and the tests/ folder (2026-09-24)

- **Rename script teaches tests/.**
  `projects/copy_and_rename_project.sh` now renames the project
  inside `tests/*.py` and `tests/README.md` with the same sed
  families it applies to the other folders, and deletes the
  untransferable donor pins — `tests/frozen`,
  `tests/manifest_sha256.json`, the measured figures, and
  `tests/__pycache__` — with a comment pointing at
  `tests/generate_frozen_reference.py --overwrite` for regeneration.
  Validated on a copy of `lsst_y1/tests`: zero old-name references
  remain and every renamed test file compiles.
- **FAQ documents tests/.** `projects/README.md`: the folder-structure
  tree and a NOTE now describe the suite and its frozen pins; the
  easy way carries a NOTE on what the script does to `tests/` plus
  the regeneration command; the hard way gained a "Changes in the
  `tests` folder" section (rename seds, pin deletions, regeneration,
  and a tests-scoped leftover-name grep). Pruning survey-specific
  entries (for example the `M2`-`M6` scale-cut datasets in
  `tests/cocoa_test_utils.py`) and re-measuring the README's tables
  and figures stay documented manual steps.
- **The one script runs on Linux and macOS.** The bash-4
  `${var,,}`/`${var^^}` expansions became tr-precomputed case
  variants (bash 3.2, the macOS /bin/bash, suffices), the linux-only
  `rename` tool became `find -depth -print0` + `mv` with bash pattern
  substitution, and a guard aborts with a message when GNU sed (from
  the cocoa environment) is missing, before anything is copied or
  deleted. Validated by a full end-to-end run on macOS 13 /
  bash 3.2.57: the FAQ final-check greps return nothing, every
  renamed python file compiles, and the renamed `tests/` is
  byte-identical to the precomputed-expansion reference; the FAQ's
  hard-way data-file rename command was made portable the same way.
- **Exhaustive leftover scan (2026-09-24).** Scanning EVERY file of a
  freshly created project (all extensions, contents and filenames,
  case-insensitive) caught three survivors the code-extension greps
  missed: the `EXAMPLE_MCMC*.covmat` header line names the sampled
  parameters with the old survey prefix (a silent proposal-matrix
  mismatch), `.gitattributes` kept the LFS pattern of the old
  covariance name (the renamed large file would escape LFS), and the
  top-level `README.md` is donor documentation that renaming would
  only falsify. The script now seds the covmat headers (bodies
  verified byte-identical) and `.gitattributes`, and replaces the
  README with a stub; the FAQ hard way gained the matching step and
  its final check became the exhaustive pair `grep -rli lsst` +
  `find -iname "*lsst*"`, both returning nothing on the validated
  macOS run.
- **Compile-and-run validation (2026-09-24).** A created project was
  taken through the full easy-way cycle on macOS: re-sourcing
  `start_cocoa.sh` auto-created its data and cobaya-likelihood
  symlinks, `compile_xxx.sh` built `cosmolike_xxx_interface.so`, and
  `EXAMPLE_EVALUATE1.yaml` (cosmic shear) and `EXAMPLE_EVALUATE2.yaml`
  (3x2pt) both ran, with the full evaluate output identical to the
  lsst_y1 donor's after name normalization (chi2 0.267263 and
  0.443061, log-posteriors -1076.03 and -1064.55). A failed `cd` in
  the rename loops now aborts instead of letting `find` rename the
  wrong folder.
- **Easy and hard way re-synced (2026-09-24).** Auditing the FAQ hard
  way against the script found drift in both directions: the script
  kept `scripts/EXAMPLE_PLOT_*.py` and `scripts/*.sbatch` (they plot
  emulator chains; now deleted), and the hard-way Final cleanup
  lacked the chains junk and `interface/*.{so,o}` removals (now
  listed). After the sync a script-created project's `scripts/`
  carries exactly the three lifecycle scripts, and both exhaustive
  final checks stay empty.
- **Project symlinks ignored everywhere.** `start_all_projects.sh`
  deleted `cobaya/cobaya/likelihoods/.gitignore` but never recreated
  it (the generated project list was copied only to
  `external_modules/{data,code}`), so every project's likelihood
  symlink sat untracked in the cobaya repository. The copy now
  reaches `cobaya/cobaya/likelihoods/` too, the generated file
  ignores itself (the cobaya repository has no other rule for it),
  and `stop_all_projects.sh` removes it symmetrically. Verified: no
  project link and no generated `.gitignore` shows in any repo's
  `git status`, for existing projects and a new one alike.

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

### High

No open HIGH tickets.

### Medium

- OPEN **MEDIUM** **BUG FIX** — [Find why the dlnxi/Css/dlnk notebook functions crash](#open-notebook-derivative-crash)
- OPEN **MEDIUM** **NEW FUNCTIONALITY** — [Move the cosmo2D_wrapper array overloads to the batched _ells API](#open-cosmo2d-wrapper-ells-api)

### Low

- OPEN **LOW** **NEW FUNCTIONALITY** — [IA x higher-order-bias (gb2) cross terms in cfastpt](#open-cfastpt-gb2-ia-bias)

<a id="open-notebook-derivative-crash"></a>
## Find why the dlnxi/Css/dlnk notebook functions crash

### High-level summary

The example notebooks call the diagnostic functions that compute the
Limber $C_{ss}(\ell)$ spectra and the $d\ln\xi_\pm/d\ln k$ derivatives
of the data vector, and those calls crash on jupyter. Until the crash
is reproduced and located, the notebook diagnostics built on them are
unusable.

### Current status

**Ticket type: BUG FIX.**

**OPEN.** Not yet reproduced under a debugger; only the notebook-level
crash is reported.

**Severity: MEDIUM.** The crash blocks notebook diagnostics, not a
likelihood evaluation or a chain; there is no evidence yet that a
primary result is wrong.

### What is already in place

The functions exist along the whole chain: pybind11 bindings
(`dlnxi_dlnk_pm_tomo_limber` and the `Css` spectrum functions) in each
project's compiled interface, python wrappers in the projects'
notebook-wrapper modules, and calls in the EXAMPLE_EVALUATE notebooks.

### What is missing

Reproduce the crash from one project's example notebook, capture the
traceback (or the signal, if the compiled core aborts), and bisect
whether it enters in the notebook wrapper, the pybind11 layer, or the
C core. Then fix where it enters and rerun the notebook end to end.

<details><summary>Technical record</summary>

- Owners: `projects/<name>/interface/interface.cpp` (the
  `dlnxi_dlnk_pm_tomo_limber` and `Css` bindings),
  `projects/<name>/interface/cosmolike_<name>_notebook_wrappers.py`,
  and the `EXAMPLE_EVALUATE*.ipynb` notebooks that call them.
- A crash inside the compiled core will not show a python traceback;
  run the reproduction under `python -X faulthandler` or gdb.

</details>

<a id="open-cosmo2d-wrapper-ells-api"></a>
## Move the cosmo2D_wrapper array overloads to the batched _ells API

### High-level summary

The python-facing array overloads `C_ss_tomo_limber_cpp(arma::Col l)`
and `C_gs_tomo_limber_cpp(arma::Col l)` in `cosmo2D_wrapper.cpp` still
compute their spectra by calling the scalar
`C_ss_tomo_limber_nointerp` / `C_gs_tomo_limber_nointerp` once per
(multipole, bin-pair), after a serial `init=1` warm-up pass. The
optimized batch entry points `C_ss_tomo_limber_nointerp_ells` and
`C_gs_tomo_limber_nointerp_ells` already exist in `cosmo2D.h` and are
what the likelihood path uses. The ss and gs array overloads should
call the batch API and drop the scalar loop, and the ss and gs scalar
overloads should be rewritten as one-element calls into the array
overloads, so the wrapper layer has no dependency left on the old
non-batched shear cosmo2D API.

### Current status

**Ticket type: NEW FUNCTIONALITY** (performance and single-code-path
consolidation; no wrong numbers are known).

**OPEN.**

**Severity: MEDIUM.** Two parallel implementations of the same spectra
(the wrapper's scalar loop and the batch `_ells` path) can drift apart
and make notebook diagnostics disagree with the likelihood path; the
scalar loop is also the slow way to fill the tables the notebooks ask
for.

### What is already in place

The batch C API with per-argument documentation
(`C_ss_tomo_limber_nointerp_ells` returning EE and BB as
`[NSIZE][nell]`, `C_gs_tomo_limber_nointerp_ells` returning
`[NSIZE][nell]`), the `_batch` thin wrappers over integer multipoles,
and the bin-pair maps `Z1`/`Z2` and `ZL`/`ZS` that translate the
power-spectrum index `nz` into the cube coordinates the python side
expects.

### What is missing

Rewrite the ss and gs array overloads to allocate the `[NSIZE][nell]`
work arrays with `malloc2d`, call the `_ells` batch functions once,
and scatter the rows into the returned `arma::Cube` via `Z1`/`Z2`
(EE and BB) and `ZL`/`ZS`. Rewrite the scalar overloads
`C_ss_tomo_limber_cpp(l, ni, nj)` and
`C_gs_tomo_limber_cpp(l, ni, nj)` as one-element calls into their
array overloads, with a doc warning that the scalar API is a
point diagnostic (it pays the full batch cost per call) and the array
API is the one to use. Validate that the new wrappers return the same
values as the current code at notebook multipole arrays before
deleting the loops, and rerun one EXAMPLE_EVALUATE notebook per
affected project end to end.

<details><summary>Technical record</summary>

- Owner: `external_modules/code/cosmolike_core/cosmolike/`
  (`cosmo2D_wrapper.cpp`; batch API declared in `cosmo2D.h`).
- Current scalar-loop overloads: `C_ss_tomo_limber_cpp(arma::Col l)`
  and `C_gs_tomo_limber_cpp(arma::Col l)` in `cosmo2D_wrapper.cpp`
  (each does a serial `init=1` pass over bin pairs, then an
  `omp collapse(2)` loop of scalar calls).
- Only ss and gs have batch entry points; gg has no `_ells` function
  (the cosmo_nodes batch refactor was never applied to gg), so the gg
  overloads are out of scope until gg gets a batch API.
- Related: the notebook-derivative crash ticket
  (#open-notebook-derivative-crash) lives in the same wrapper layer;
  if the migration touches shared state, retest that reproduction.

</details>

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

# Closed tickets

Grouped by subject and compressed. Nothing here is open work; dated
measurements are kept, since this section is a decision record, not a
README. To reopen a ticket, move its content back under
[Open tickets](#open-tickets) as a full ticket section and add its
`- OPEN` index line.

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

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
- OPEN **MEDIUM** **BUG FIX** — [Teach copy_and_rename_project.sh and the projects FAQ about tests/](#open-new-project-tests)
- OPEN **MEDIUM** **NEW FUNCTIONALITY** — [Two-grid sampling inside cfastpt for faster TATT](#open-cfastpt-two-grid)

### Low

No open LOW tickets.

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

<a id="open-new-project-tests"></a>
## Teach copy_and_rename_project.sh and the projects FAQ about tests/

### High-level summary

Projects now ship `tests/`: a pytest suite bound to the shared
`CocoaTestHarness`, manifest-pinned `frozen/` state, and a
`tests/README.md`. `projects/copy_and_rename_project.sh` never touches
`tests/`, and the FAQ "How do we create a new Cosmolike project?"
(the easy way and the hard way) predates it, so a copied project
carries the donor's tests verbatim: name references, frozen manifests,
and README all point at the old project.

### Current status

**Ticket type: BUG FIX.**

**OPEN.** Verified: the script contains no reference to `tests/`.

**Severity: MEDIUM.** Every new-project creation inherits misleading
tests; existing projects are untouched.

### What is already in place

The script handles the rename of the likelihood, interface, and data
layers; the FAQ documents both creation paths.

### What is missing

Make the rename script rewrite the project name inside `tests/` and
either regenerate or instruct the user to regenerate `frozen/` (its
manifests pin the donor's state and are never hand-edited). Update
both FAQ paths in `projects/README.md` to say what `tests/` is and
what a new project must do with it.

<details><summary>Technical record</summary>

- Owners: `projects/copy_and_rename_project.sh` and
  `projects/README.md` anchors `appendix_projects_new`,
  `appendix_projects_new_easy`, `appendix_projects_new_hard`.
- The harness binding lives in each project's
  `tests/cocoa_test_utils.py`; the shared machinery is
  `external_modules/code/cosmolike_core/cocoa_testing.py`, whose
  module docstring carries a worked shim example.

</details>

<a id="open-cfastpt-two-grid"></a>
## Two-grid sampling inside cfastpt for faster TATT

### High-level summary

The FAST-PT theory block computes its convolutions on a coarse
internal grid and upsamples with a cubic spline in $\log k$ onto the
dense output table the likelihood reads linearly. cfastpt
(`IA_code: 0`) still runs everything on one grid; the same two-tier
strategy could make the C path faster for TATT.

### Current status

**Ticket type: NEW FUNCTIONALITY.**

**OPEN.** Idea stage; no cfastpt change exists.

**Severity: MEDIUM.** Performance work; current cfastpt is correct
and remains the reference implementation.

### What is already in place

The strategy is proven twice in python: the FAST-PT theory block and
the `bfmt` baryon block both separate the computation grid from the
output grid. The upgrade located the accuracy in the output-table
density, and none in the convolutions, so a coarse convolution grid
is safe.

### What is missing

Implement coarse-convolution plus spline-upsampling inside cfastpt,
verify the result against current cfastpt with the existing
comparison machinery below the test tolerance, and measure the
speed-up before adopting. Keep the single-grid path until the
two-grid version matches it.

<details><summary>Technical record</summary>

- Owner: `external_modules/code/cosmolike_core/cfastpt/`.
- Precedents: `external_modules/code/PyFAST-PT/fastpt.py` (the
  two-grid mechanism and its README's Design section) and
  `code/baryon_suppression/bfmt.py` (the nz/nk internal grid).
- The grid-density analysis is the decision record in
  `projects/lsst_y1/tests/README.md`.

</details>

# Closed tickets

Grouped by subject and compressed. Nothing here is open work; dated
measurements are kept, since this section is a decision record, not a
README. To reopen a ticket, move its content back under
[Open tickets](#open-tickets) as a full ticket section and add its
`- OPEN` index line.

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

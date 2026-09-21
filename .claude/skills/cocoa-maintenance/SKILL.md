---
name: cocoa-maintenance
description: >
  Rules for editing Cocoa's README files and bash installation scripts
  (set_installation_options.sh, installation_scripts/*.sh, setup_cocoa.sh,
  compile_cocoa.sh, start/stop scripts). Use this skill whenever the task
  touches Cocoa's README.md, a project README, any *.sh file, the conda yml
  files, release tags, or environmental keys. This repository does NOT contain
  the Cosmolike C code; this skill covers only Cocoa's shell layer and
  documentation.
---

# Cocoa maintenance: README and bash scripts

Follow these rules exactly. When a rule below conflicts with your own idea of
"better code", the rule wins. Do not improvise improvements.

## 0. Absolute rules (read first)

1. **NEVER merge, deduplicate, or refactor the bash scripts.** Every
   `setup_XXX.sh` / `compile_XXX.sh` / `unxv_XXX.sh` intentionally repeats the
   same helper functions and structure. This repetition is a design decision:
   each script must be readable top-to-bottom on its own and must work mostly
   standalone. Do NOT create shared helper libraries, do NOT extract common
   functions into a sourced file, do NOT "simplify" `set_installation_options.sh`
   by looping over keys. If you think a merge would help: it does not. Skip it.
2. **Make minimal diffs.** Change only the lines the task requires. Never
   reformat, re-indent, or re-wrap lines you are not changing.
3. **Never `git push`.** Work on a branch (usually `bugfix`), commit when
   asked, and let a human push and squash-merge.
4. **Never commit in the wrong repository.** `Cocoa/projects/<name>/`,
   `Cocoa/external_modules/code/<name>/`, and `Cocoa/cobaya/` may be separate
   git repositories. Before committing, run `git rev-parse --show-toplevel`
   and confirm you are in the repository you intend.
5. **Verify before claiming.** After editing a script run `bash -n <file>`.
   After editing a README run the render check (Section 3.7). Report the
   actual output; if a check fails, say so.

## 1. Repository layout facts

- `set_installation_options.sh` — every user-facing environmental key:
  `IGNORE_*` (skip a package), `OVERWRITE_*`, URLs, and version pins.
- `installation_scripts/` — one `setup_*.sh` (download), one `compile_*.sh`
  (build/install), and for large datasets one `unxv_*.sh` per package.
- `installation_scripts/flags_impl_unset_keys.sh` — unsets EVERY key and
  function Cocoa defines. Kept in alphabetical order (variables, then
  functions).
- `installation_scripts/flags_derived.sh` — derived keys (for example the
  `OVERWRITE_EXISTING_ALL_PACKAGES` cascade) and `PYTHON3`/`PIP3`.
- `cocoapy311*.yml` (repo root) — conda environment files. `-linux` and
  `-osxarm` are conda-lock files; do not hand-edit locked package lists.
- `cocoa_installation_libraries/docker/` — the whovian-cocoa image recipe.
- `README.md` (repo root) — the main documentation. Project READMEs live in
  each project repository.

## 2. Bash script rules

### 2.1 Scripts are sourced, never executed

Every installation script runs via `source script.sh`. Consequences:

- Use `return`, never `exit` (`exit` kills the user's shell).
- `return 99` means "skipped because the IGNORE key is set".
- `return 55` at the end means "ran completely" — `setup_cocoa.sh` and
  `compile_cocoa.sh` cache success only on 55. Keep this line and its comment.
- Every variable and every function the script defines MUST be unset before
  returning, or it leaks into the user's shell. That is what
  `unset_env_vars` / `unset_env_funcs` / `unset_all` are for.

### 2.2 The standard script skeleton

When creating a new `setup_XXX.sh` or `compile_XXX.sh`, copy an existing one
(`setup_pyspk.sh` and `compile_pyspk.sh` are clean models) and rename. The
skeleton is:

```bash
#!/bin/bash
if [ -n "${IGNORE_XXX_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL ECODEF FOLDER PACKDIR    # every variable the body defines
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error      # every function defined below
  unset -f unset_env_funcs
  cdroot || return 1;
}

unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

# ... local helpers (cdfolder, cpfolder, ...), then the body ...

unset_all || return 1

return 55; # setup/compile_cocoa caches this script only on full completion
```

If you add a function, add it to `unset_env_funcs`. If you add a variable,
add it to `unset_env_vars`. No exceptions.

### 2.3 Error-handling contract (do not get this wrong)

`error()` prints the failure message AND runs the full cleanup (`unset_all`).
Therefore:

- A call to a helper whose failure path goes through `error()` — that is
  `cdfolder`, `cpfolder`, `cpfile`, `gitact*`, `wgetact*`, and similar — uses:

      cdfolder "${PACKDIR:?}" || return 1;

  NOT `|| { unset_all; return 1; }`. The cleanup already ran inside `error`;
  calling `unset_all` again prints `unset_all: command not found`.

- A direct command (git, tar, pip, cp without a helper) that does not go
  through `error()` uses either:

      some_command ... || { error "MESSAGE"; return 1; }

  or, when no message is wanted:

      some_command ... || { unset_all; return 1; }

- `ptop`/`pbottom` calls keep `|| { unset_all; return 1; }` (they do not
  clean up on failure).

Never `raise`/propagate by letting a command fail silently: every failure
path must either call `error` or perform cleanup explicitly.

### 2.4 rm safety

Every `rm` uses `:?`-guarded absolute variables so an unset variable can
never expand to `/`:

```bash
PRJ="${ROOTDIR:?}/projects/${NEW_PROJECT:?}"
rm -rf "${PRJ:?}"
rm -f  "${PRJ:?}"/chains/*.txt 2>/dev/null
```

Never write `rm -rf $VAR` or `rm -rf ${VAR}/` without the `:?` guard.

### 2.5 Environmental keys

Adding a key requires touching, in this order:

1. `set_installation_options.sh` — define it in the matching section, with
   the same comment style as its neighbors. Keys that default OFF are written
   commented: `#export KEY=1`. Keys the user may uncomment to skip a package
   are named `IGNORE_XXX_CODE` / `IGNORE_XXX_DATA`.
2. `installation_scripts/flags_impl_unset_keys.sh` — add `unset -v KEY` in
   case-insensitive alphabetical position, never at the end (functions go
   under `# Functions` as `unset -f`). Verify with Section 2.14's checks.
3. `installation_scripts/flags_derived.sh` — only if the key participates in
   a cascade (for example each `OVERWRITE_EXISTING_XXX_CODE` is set to 1 when
   `OVERWRITE_EXISTING_ALL_PACKAGES` is set).
4. If a new setup/compile script consumes it: add the script name to the
   matching array (`CORE`/`THEORY`/`ML`/`LIKELIHOOD`) in `setup_cocoa.sh`
   and/or `compile_cocoa.sh`. Adding or removing a list entry is safe: the
   cache file rebuilds automatically when the list length changes.
5. If the key's script clones or downloads into `external_modules/code` or
   `external_modules/data`: add the destination folder to `Cocoa/.gitignore`
   (Section 2.12). This step is as mandatory as the other four.

### 2.6 Version pins

- Every external repository is pinned with a trio:
  `XXX_URL`, `XXX_NAME`, and exactly one of
  `XXX_GIT_COMMIT` / `XXX_GIT_BRANCH` / `XXX_GIT_TAG`.
- Precedence when more than one is set: COMMIT beats BRANCH beats TAG.
  An accidentally active `XXX_GIT_BRANCH="main"` silently unpins the package.
- Release rule: every dependency must be pinned; debug and developer keys
  (`COCOA_OUTPUT_VERBOSE`, `COCOA_OUTPUT_DEBUG`, `COSMOLIKE_DEBUG_MODE`,
  `SWITCH_TO_DEV_MODE`) must ship commented out.
- Before writing a tag into this file or the README, confirm it exists:

      git ls-remote --tags <URL> | grep <tag>

### 2.7 The internet invariant (strict)

- `setup_*.sh` and `unxv_*.sh` MAY use the internet (git clone, wget,
  pip index access, emulator model-file downloads).
- `compile_*.sh` MUST NEVER use the internet. They may run on compute nodes
  without network access. Compile-time pip installs use exactly:

      ${PIP3:?} install "${PACKDIR:?}" \
        --no-dependencies \
        --prefix="${ROOTDIR:?}/.local" \
        --no-index \
        --no-build-isolation

- If a package downloads data files on first use (some emulators do),
  trigger that download in its `setup_*.sh`, never at compile or run time.
- If a step needs both installing and downloading, put the whole step in
  setup — do not split it into compile.

### 2.8 pip dependency rules

- Runtime Python dependencies of packages installed with `--no-dependencies`
  are seeded in the `PIPCP` array of
  `installation_scripts/setup_pip_core_packages.sh`, pinned by version.
- `.local` site-packages SHADOWS the conda environment. Before adding any
  package to `PIPCP`, run in the cocoa environment:

      pip install --dry-run --prefer-binary <the full PIPCP list>

  If any already-pinned package (numpy, scipy, scikit-learn, matplotlib,
  pandas, typing-extensions, ...) appears in the "Would install" line, add an
  equal-version guard pin to `PIPCP` instead of letting pip upgrade it.
- Determine a package's real runtime dependencies from its `import`
  statements, not from `requirements.txt` (which mixes in test/docs tools).
- The `PIPCP_HASH` sentinel hashes the array: editing the array automatically
  triggers reinstall on the next setup run. Do not bypass it.

### 2.9 Developer SSH mode

`SWITCH_TO_DEV_MODE=1` switches clones of repositories the team owns from
https to SSH. The mechanism is a small local `devurl()` function copied into
each consuming setup script (yes, copied — see rule 0.1):

```bash
devurl() {
  # SWITCH_TO_DEV_MODE=1: rewrite GitHub https URLs to their ssh form
  local U="${1:?}"
  if [ -n "${SWITCH_TO_DEV_MODE:-}" ]; then
    case "${U}" in
      https://*github.com/*) U="git@github.com:${U#*github.com/}" ;;
    esac
  fi
  echo "${U}"
}
...
URL="${XXX_URL:-"https://github.com/CosmoLike/xxx.git"}"
URL=$(devurl "${URL:?}")
```

Apply it ONLY to repositories the team has write access to (CosmoLike/SBU
organizations, emulator-training repos). Third-party upstream codes are
cloned via plain https with no `devurl`.

### 2.10 Misc script facts

- `sed` is GNU sed (from the conda environment): the house style is
  `sed --in-place --regexp-extended`.
- Lists of items in scripts are written as bash arrays with loops, not as
  long repeated sed lines — but never loop across DIFFERENT scripts (rule 0.1).
- Some project repositories use Git LFS; git operations there need `git-lfs`
  on PATH (the cocoa conda environment provides it).
- After a squash merge, delete the local topic branch with
  `git branch -D <branch>` (`-d` refuses: squash-merged branches never look
  merged to git) and start the next branch fresh from updated `main`.
- **NEVER `git push --force`** (including `--force-with-lease`). When
  a remote topic branch is stale (for example the pre-squash leftover
  of a merged branch), the maintainer deletes it on GitHub and the
  local branch is pushed plain.
- **NEVER `git commit --amend`.** Follow-up fixes go in a new commit.
  The maintainer may push the branch at any moment; amending rewrites
  the local history and the next `git pull` then ends in a merge
  conflict between the amended commit and its pushed twin. The branch
  is squash-merged in the end, so extra small commits cost nothing.

### 2.11 Conservative version policy (Python, NumPy, everything)

Cocoa deliberately runs OLD, known-good versions. Never upgrade a version as
a side effect of another task, and never pick "latest" when adding a pin.

- Python stays at the pinned version (`export PYTHON_VERSION=3.11` in
  `set_installation_options.sh`) until that version is within about a month
  of its upstream end-of-life; only then does Cocoa move to the next minor
  (3.11 → 3.12, and so on). Until that point, do not propose newer Pythons,
  and the migration itself is a maintainer-led task, never a side effect.
- NumPy is capped BELOW 2.0: `COCOA_NUMPY_VERSION` in
  `setup_pip_core_packages.sh` is `1.26.3` (`1.23.5` under
  `COCOA_FORCE_NUMPY_1_23`), and 1.26.x is the permanent maximum. NumPy 2.x
  breaks the compiled stack (CARMA/cosmolike bindings, older scipy/numba
  wheels). That is why `numpy==${COCOA_NUMPY_VERSION}` is repeated inside
  the pip commands in that script: pip must never get a chance to resolve
  `numpy>=2` on its own.
- The same idea applies to every package: when a new dependency is needed,
  pin the oldest version that works with the already-pinned stack, not the
  newest release. Emulator model files are often serialized against one
  exact library version (BCemu loads its trained files only with
  `smt==1.0.0`; other versions are incompatible).
- Version bumps are a deliberate maintainer task of their own (tested on
  Linux and macOS, conda-lock files regenerated) — never bundled into a
  feature or bugfix change.
- ENFORCEMENT — be obsessive about this. The policy only holds if versions
  are locked tightly everywhere a resolver runs, above all in
  `setup_pip_core_packages.sh` and the conda ymls:
  - Every entry you add to `PIPCP`/`PIPCPFR` carries an exact `==` pin.
    Never an unpinned name, never `>=`, never a range: any of those is an
    open door for pip to upgrade the stack later.
  - A pin that exists in both a conda yml and a pip array must be the SAME
    version in both places; they protect each other.
  - Before AND after touching these files, run the Section 2.8 dry-run.
    If any already-pinned package (numpy, scipy, matplotlib,
    typing-extensions, ...) shows up in pip's "Would install" line, that is
    a bug in your change — stop and add the missing guard pin. Never accept
    the upgrade, never "fix" it by bumping the existing pin.
  - The same discipline applies to conda recipe edits: adding a package to
    a yml without an exact version lets the solver move other packages.
    Pin the addition and confirm the solver keeps everything else fixed.
  - If a task appears to REQUIRE upgrading a pinned package, do not do it:
    report the conflict to the maintainer and stop. An upgrade that
    violates the intent of `setup_pip_core_packages.sh` or the conda
    recipes is never an acceptable side effect, even if it makes the task
    "work".

### 2.12 .gitignore review for installed packages (quite important)

Every `setup_*.sh` clone destination and every `unxv_*.sh` data folder
under `external_modules/code/` and `external_modules/data/` MUST have a
matching entry in `Cocoa/.gitignore`.

Why this matters: cloned packages are full git repositories. If one is not
ignored, a user running `git add --all` in Cocoa stages that directory as a
bare gitlink — an accidental submodule with no `.gitmodules` entry — which
commits an unusable pointer and breaks everyone else's clone. Users must
always be able to run `git add --all` safely; the `.gitignore` is what
guarantees that.

Rules:

- When adding a package, add its destination folder to `Cocoa/.gitignore`
  next to its neighbors, using the DEFAULT `XXX_NAME` value (the entry
  cannot follow a user's rename). Example entries from the baryon work:

      external_modules/code/baryon_suppression
      external_modules/code/pyspk
      external_modules/code/bcemu

- Data folders get `external_modules/data/<name>` entries; families may use
  a wildcard (`external_modules/data/roman*`), matching the existing style.
- Verify after running the new setup script, from the repository root:

      git check-ignore -v "external_modules/code/<folder>"
      git status --porcelain | grep external_modules

  The first must print the matching rule; the second must print nothing.
- Review this whenever a script's destination folder is added or renamed —
  a renamed `XXX_NAME` default silently orphans the old `.gitignore` entry.

### 2.13 Science packages: always install the Cocoa way (the getdist pattern)

When a task or a user needs a new scientific package — a cosmology code, an
emulator, a sampler, an analysis tool — NEVER suggest `pip install <pkg>`
into the environment. The answer is always the key mechanism in
`set_installation_options.sh`, with `getdist` as the reference example:

- Keys (set_installation_options.sh):

      #export IGNORE_GETDIST_CODE=1 #dev getdist with code tweaks
      export GETDIST_URL="https://github.com/cmbant/getdist.git"
      export GETDIST_GIT_COMMIT="ff477beea2e7e2231a3de4941bdc3d64bd1f0bb4"
      export GETDIST_NAME="getdist"

- `setup_getdist.sh` clones the pinned commit into
  `external_modules/code/getdist`; `compile_getdist.sh` pip-installs that
  clone offline into `.local` (Section 2.7 flags). Both are wired into the
  runner arrays, and the destination is gitignored (Sections 2.5, 2.12).

Why this is non-negotiable: scientists-as-developers break backward
compatibility ALL the time. A directly pip-installed package floats with
upstream and silently diverges between users and machines; the key
mechanism pins an exact commit, so every Cocoa installation is identical
and a breaking upstream change can never reach users until the maintainer
deliberately moves the pin (Section 2.11).

Consequences:

- Adding a science package = the full recipe: pin-trio keys + IGNORE key,
  setup (and, if offline-installable, compile) scripts, flag bookkeeping,
  runner lists, `.gitignore`, and PIPCP pins for its runtime dependencies.
  There is no shortcut version of this.
- The rule extends to documentation and support answers: when a README or
  a reply explains how to get a package, it shows the key block to
  enable/add — it never tells users to `pip install` the package directly
  (the theory-block READMEs carry an explicit warning about this).

### 2.14 Environment hygiene: the shell must be clean after stop_cocoa.sh

Cocoa is obsessive about leaving the user's shell exactly as it found it:
after `source stop_cocoa.sh`, the environment must be indistinguishable
from before `source start_cocoa.sh`. Hygiene has three mechanisms, and
every variable belongs to exactly one of them:

- per-script `unset_env_vars`/`unset_env_funcs`/`unset_all`
  (Sections 2.1–2.2) clean the names a script defines for itself;
- `flags_impl_unset_keys.sh` (sourced by `stop_cocoa.sh`) UNSETS every
  Cocoa-invented exported key (from `set_installation_options.sh`,
  `flags_derived.sh`, `start_cocoa.sh`) that deliberately survives across
  scripts until stop;
- `flags_save_old.sh` / `flags_recover_old.sh` handle variables that
  ALREADY EXIST in the system and that Cocoa merely modifies (`PATH`,
  `LD_LIBRARY_PATH`, `PYTHONPATH`, `C_INCLUDE_PATH`, `LDFLAGS`,
  `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `CUDA_VISIBLE_DEVICES`, ... — 21 in
  total, mirrored between the two files). `start_cocoa.sh` saves each into
  `OLD_<NAME>` (with the sentinel value `"x"` meaning "was unset");
  `stop_cocoa.sh` restores the saved value (or unsets, on the sentinel)
  and drops the `OLD_` copy. These variables must NEVER go into
  `flags_impl_unset_keys.sh`: a blind unset would destroy the user's own
  pre-Cocoa value instead of recovering it.

Rules:

- When adding a variable, first decide its class. A name Cocoa invents
  (keys, pins, `IGNORE_*`, ...) goes in `flags_impl_unset_keys.sh`. A
  standard system/toolchain name Cocoa modifies goes as a matching pair of
  blocks in `flags_save_old.sh` AND `flags_recover_old.sh`, copying the
  existing `"x"`-sentinel pattern.

- Every new key gets an `unset -v KEY` line in the `# Variables` block;
  every new function gets `unset -f name` in the `# Functions` block. Add
  ALL spellings of a pin family (`_URL`, `_NAME`, `_GIT_COMMIT`,
  `_GIT_BRANCH`, `_GIT_TAG`, `IGNORE_`, `OVERWRITE_EXISTING_`) as
  insurance, even the ones not currently exported.
- Both blocks are kept in CASE-INSENSITIVE alphabetical order, and a new
  key is inserted at its position — NEVER appended at the end. The
  ordering exists so a developer can spot-check for a missing key at a
  glance; an appended block defeats that. (This happened: the
  baryon-emulator key families shipped appended after `WGET_VERSION` and
  were only caught later by the sort check below.)
- Verification, after any key change:

      grep '^unset -v' installation_scripts/flags_impl_unset_keys.sh | sort -cf
      grep '^unset -f' installation_scripts/flags_impl_unset_keys.sh | sort -cf

  Both must print nothing. Then the coverage check:

      comm -23 \
        <(grep -hoE '^(export )?[A-Za-z][A-Za-z0-9_]+=' \
            set_installation_options.sh \
            installation_scripts/flags_derived.sh \
          | sed -E 's/^export //; s/=$//' | sort -u) \
        <(grep -oE '^unset -v [A-Za-z0-9_]+' \
            installation_scripts/flags_impl_unset_keys.sh \
          | awk '{print $3}' | sort -u)

  Every name this prints must belong to the save/recover class — confirm
  it appears in BOTH `flags_save_old.sh` and `flags_recover_old.sh` (as of
  v4.11.5 the output is exactly `OMP_NUM_THREADS` and `OMP_PROC_BIND`,
  the two save/recover variables that `set_installation_options.sh`
  exports directly). A name in neither the unset file nor the
  save/recover pair is a hygiene bug: fix it in the class it belongs to.
- End-to-end smoke test: `source start_cocoa.sh` then
  `source stop_cocoa.sh`, then

      env | grep -E '_GIT_(COMMIT|TAG|BRANCH)=|^IGNORE_|^OVERWRITE_'

  must print nothing.

## 3. README rules

### 3.0 Imitate before you write (the governing rule)

Every correction the maintainer made to freshly written README pages
was already visible in the main README: the numbered Contents list
with `<a name>` anchors, the assumption paragraph followed by
`**Step :one:**:` blocks, NOTE/TIP callouts for asides, appendix
FAQs for detail beyond the main flow, tables and indented
`key: value` blocks for anything enumerable, LaTeX for math. So:
before writing or editing ANY README in this tree, open the main
README and the closest existing page of the same kind, and build the
new page by imitating their shapes section by section. Never compose
a page fresh from prose and let it get patched back into the house
shape rule by rule. When a page and the main README disagree in
shape, the main README wins.

### 3.1 Writing style (mandatory)

- Golden rule: **"less is more, and when you need more — put it into an
  appendix."** The main flow is commands; explanations go to appendix FAQs.
- Concrete actions only. No adjectives, no marketing ("faster", "robust",
  "powerful" — delete on sight). Describe what a command does, literally.
- TIPs in the main flow are pure signposts: "If X happens, see the appendix
  [FAQ: ...]" — zero inline explanation.
- Notes explaining command flags are flag-first, one item per flag:
  `` `--flag value` ``: imperative action. Do not re-enumerate in prose what
  the command block already shows.
- Appendix FAQ titles are humble: "How can users deal with X (Platform): a
  possible cause" — never claim the only cause.
- FAQ bodies: one or two intro sentences ending "If that is the case, follow
  the steps below", then `**Step :one:**:` blocks — one short imperative
  line plus one command each.
- Every command flow uses `**Step :one:**:` blocks — ALWAYS, even
  when there is a single step. A bare command block with a lead-in
  sentence is not a flow.
- Every flow is self-contained: repeat the assumption paragraph
  ("We assume users are in the Conda cocoa environment ...") and the
  `source start_cocoa.sh` step in EVERY flow. Never write "with the
  environment of [some other section]" — the reader landing from
  the table of contents must not chase references.
- Table rows are atomic: one item per row (one test, one knob, one
  file), so an item number maps to exactly one row. Do not condense
  ("1-4", "NLA and TATT" in one cell); expand instead. Checks and
  configuration are different columns: NLA/TATT is configuration
  (write "IA modeling: NLA"), a race test is a check (write "race
  condition (OpenMP threading)", never bare "race").
- Name things by physics, never by internal bookkeeping: "cosmic
  shear", "3x2pt", "6x2pt" — not "example1"/"example2". Example
  numbers may appear only where the reader touches the actual file
  (a `test_example1.py` filename, an `example1.dataset` path).
- A README with more than three sections opens with a numbered
  Contents list linking `<a name>` anchors placed on each heading
  (the main README pattern).
- One topic = one TIP (a tip may contain two numbered solutions; do not
  split into two tips).
- A paragraph must not carry what a display shows better. Settings
  with values, parameter lists, and commands go into indented blocks
  (the main README's mpirun blocks are the model), enumerable facts
  and measured numbers into tables, sequences of actions into
  numbered lists. NO-GO: "(cosmolike accuracyboost 2,
  integration_accuracy 10, lmax 200000, kmax_boltzmann 40; CAMB
  AccuracyBoost 2, k_per_logint 50)" buried inside a sentence. GO: a
  short sentence, then an indented block with one `key: value` per
  line.
- Paragraph size: about six rendered lines is the ceiling. A
  ten-line paragraph that chains a mechanism, a report format, a
  cost estimate, and an activation switch is four paragraphs (or a
  numbered list plus two short paragraphs) pretending to be one.
- Math is LaTeX, settings are code spans. Write $\chi^2 < 0.2$ and
  $\Delta\chi^2$, never the bare words "chi2 < 0.2"; write
  `accuracyboost: 2`, never accuracyboost 2. Inside command blocks
  and quoted terminal output keep everything verbatim. In table
  cells never put a raw | inside math: use \lvert and \rvert for
  absolute values, or the pipe breaks the table.
- A list item is one short clause: the action, nothing else. The
  reasoning ("so the chi2 against it is zero by construction") moves
  to a sentence after the list that refers back to the steps.
- Numbered lists count 1, 2, 3 with no gaps and no grouped items:
  the renderer recounts an ordered list, so "5. -8." followed by
  "9. -14." renders as "6.", "7." and the printed numbers lie.
  Grouped ranges belong in a table.
- READMEs name only what users type or open: commands, files,
  folders, yaml keys. Never internal python function names (naming
  the function that implements a check helps nobody); a pytest
  selector inside a command (`-k nmodels`) is part of the command
  and stays.
- Operational caveats (a pager swallowing output, a platform quirk,
  anything of the form "if the terminal does X, do Y") are
  `> [!NOTE]` or `> [!TIP]` callouts, never body paragraphs.
- Banned words in READMEs (user-vetoed): "suite" (say "the tests")
  and "keypress" (say what is pressed, e.g. "until space or enter is
  pressed"). Identifiers in code stay as the code spells them.
- Titles never repeat the repository name: inside the lsst_y1
  project the tests page is "Unit tests for the likelihoods", never
  "Unit tests for the lsst_y1 likelihoods". The repository already
  scopes its pages.
- On a project tests page the appendix heading is exactly
  `# Appendix <a name="appendix"></a>`; never a descriptive tail
  ("Appendices about the frozen state").
- FAQ titles are short, simple questions ("FAQ: Do the tests keep
  their own data?"); a compound question ("How do the tests keep
  their own copy of configurations and data?") is cut down and the
  detail moves to the body.
- Each check family gets a concept section ("Accuracy checks") that
  says what is measured, and a "Running ..." subsection under it
  that opens directly on the `**Step :one:**:` blocks.
- Write $\Delta\chi^2$, never bare $\chi^2$, for every drift or
  accuracy quantity: the $\Delta$ tells the reader the value should
  be zero. Bare $\chi^2$ stays only where the raw statistic itself
  is meant.
- Coined labels are defined at first use or replaced with everyday
  words: "stored"/"reference" instead of decorative "frozen" (the
  `frozen/` path stays literal, with one sentence at first mention
  saying what the snapshot is); "the fiducial evaluated on its own
  vs after nine other cosmologies" instead of "fresh vs 10th-of-10";
  "accuracy parameter" or "setting" instead of "knob".
- Counts name their members: "cosmic shear, 3x2pt, and 2x2pt", never
  "the three probes"; "NLA and TATT", never "both IA models".
- Sentences lead with the subject and action; fronted contrast
  openers ("Instead of the one frozen fiducial point, this
  check...") are rewritten as direct statements with the contrast at
  the end or dropped.
- Explanations restate their setup in place ("Again, we set the data
  vector at the point itself, so the $\Delta\chi^2$ is zero by
  construction"), never point at another step or section ("the
  vector of step 2").
- No ALL-CAPS emphasis in body text ("GENERATED WITH TATT", "the
  DEFAULT settings"); the sentence's contrast carries the stress.
- No commentary about other projects or past incidents ("which in
  other projects exposed interface breakdowns"): state what the
  check does here.
- Never quote measured values in a README: they go stale with the
  next commit and force a page redo. State the conclusion ("the
  measured values sit far below the 0.2 band") and point at the
  checks that print the numbers; the terminal output on the current
  code is the source of truth. Contract numbers (pass limits such as
  0.2, tolerances such as 1e-4) are not measurements and stay.

### 3.2 Quote blocks must match the source files

Blocks introduced by `[Adapted from <file>]` quote real files. Rules:

- The quoted lines must match the current file content byte-for-byte, except
  where `(...)` marks elided regions.
- When you edit `set_installation_options.sh` (keys, tags, comments), grep
  the README(s) for quotes of the changed lines and update them in the same
  change:

      grep -n "KEY_NAME" README.md

- When quoting keys, show them in the SHIPPED default state (commented if
  the default is off), not in your local working state.

### 3.3 Never present text as images

Code, terminal output, key listings, and data tables are markdown code
blocks or markdown tables — never screenshots. Plots and diagrams may be
images.

### 3.4 Version references: pointers vs boundaries

Two kinds of version mentions exist. Before a release, bump the pointers and
leave the boundaries alone:

- POINTERS (bump every release): `refs/tags/<tag>/` wget URLs,
  `git clone ... --branch <tag>`, the Docker `COCOA_TAG` default and build
  command, sentences of the form "replace the tag `<tag>` in the URL".
- BOUNDARIES (never bump — they are historical facts): sentences like
  "From `v4.11.4` on, the yml files are named `cocoapy311-*`",
  "tags older than `v4.11.4` carry the `cocoapy310` prefix",
  "installed before `v4.11.4`".

After bumping, run `grep -n "v4\." README.md` and justify every remaining
old-version mention as a boundary.

### 3.5 mpirun command conventions

- All Linux MPI commands: `--mca pml ob1 --mca btl vader,tcp,self`
  (never `^ucx`).
- Commands that can span multiple nodes additionally carry
  `--mca btl_tcp_if_exclude lo,docker0,virbr0,ib0`, the `-x` environment
  forwarding list, and `--mca mpi_yield_when_idle 1`.
- Mapping: OpenMP-heavy commands (plain MCMC, hybrid EMUL2) use
  `--map-by numa:pe=${OMP_NUM_THREADS}`; single-thread emulator (EMUL)
  commands use `--map-by slot`.
- Long commands break after `-n X --oversubscribe \` with
  `--mca pml ob1 --mca btl vader,tcp,self \` on its own continuation line.
  The same line-break style applies to sbatch scripts.
- macOS commands stay bare `mpirun` with no mca flags. Never add the Linux
  flags to macOS commands "for consistency" — the omission is deliberate.
- Project READMEs carry ONE standard note ("Running on more than one node",
  four flag-first items). It must stay byte-identical across the main README
  and all project READMEs; verify with `md5` on the extracted note.

### 3.6 Anti-AI prose rules (README, comments, docstrings, commit messages)

These rules apply to every piece of explanatory text: README prose,
Python comments and docstrings, error messages, and commit messages.
The target reader is a physics student who knows no AI-agent language
and may know little Git.

- **Describe the current state, never the editing history.** When a rule
  changes, rewrite the explanation in place. Never write "now", "new
  rule", "previously", "as of", a date attached to a rule, or a
  reference to the request/review that caused the change. Git holds the
  history.
- **Abstractions need real examples.** A broad term (test, key, block,
  freeze) gets one or two concrete repository examples nearby: a real
  file or command, the action, and the visible result.
- **Define terms where first used, keep one name per object.** Never
  rotate synonyms (script/tool/utility for the same file); never write
  "X (called Y in the code)" before the difference is explained.

#### Replace generic praise with checkable facts

Use concrete actors and verbs. "The watcher moves the file" is easier
to check than "a filesystem transition occurs." "`--cycle 2` admits no
more than two tickets" carries more information than "this provides
useful runtime control." The following claims need evidence or
removal:

- `easy`, `simple`, `intuitive`, `automatic`, `safe`, `fast`,
  `lightweight`, `flexible`, `powerful`, `production-ready`,
  `seamless`, `complete`, and `supported`;
- "handles errors", "works out of the box", "uses best practices",
  "improves performance", and "provides a better experience";
- "standard approach", "conventional method", or "robust workflow"
  without the exact algorithm, version, condition, test, or failure
  behavior.

For example, replace "safe automatic shutdown" with the actual
condition: "The watcher stops starting jobs, waits for jobs already
running, and then prints the Ctrl-C countdown." Replace "fast" with a
measured time and machine, or remove the claim.

#### Vocabulary bans

- do not use `thereby`;
- do not use `commendable`, `innovative`, `meticulous`, `intricate`,
  `notable`, or `versatile` as adjectives;
- replace decorative `delve`, `crucial`, `comprehensive`, `notably`,
  `underscores`, `highlights`, `showcases`, `sheds light on`,
  `leveraging`, and `utilize` with the fact they were trying to
  decorate;
- keep `robust`, `robustness`, and other domain terms when they carry
  a precise technical meaning;

#### Hard-zero words

Prohibited in all prose. An exact command, code identifier,
external title, or quotation may contain one only when changing it
would make the reference false; name every such exception.

    commendable    comprehensive   crucial        crucially
    delve          delves          delving        dwelve
    innovative     intricate       meticulous     multifaceted
    notable        notably         pivotal        thereby
    transformative versatile

These forms are prohibited for the same reason:

    showcase   showcases   showcasing
    underscore underscores underscoring
    unveil     unveils     unveiling

#### Hard-zero phrases

    as an AI                     it is worth noting
    certainly, here is           plays a crucial role
    complex interplay            provides valuable insights
    deeper understanding         sheds light on
    evolving landscape           stands as a testament
    future work should explore   taken together, these results
    important implications       the realm of
    it is important to note      valuable insights

Also remove bot closers and chat residue such as "I hope this helps",
"Let me know if you would like", "Here is the revised version",
"In conclusion", and "Overall" when they merely end a generated
answer.

#### Decorative vocabulary to replace

These words and phrases usually hide a simpler verb or a missing
fact. Replace them unless a literal technical meaning is identified:

    actionable insights          landscape
    advance our understanding    leverage
    at the forefront             meaningful
    broader implications         nuanced
    cornerstone                  offer insights
    cutting-edge                 pave the way
    dynamic                      poised to
    ecosystem                    powerful
    elevate                      realm
    empower                      roadmap
    enhance                      seamless
    facilitate                   state-of-the-art
    foster                       streamline
    harness                      synergy
    holistic                     tapestry
    impactful                    unlock
    in order to                  utilize
    in the context of            valuable
    insightful                   vibrant
    journey                      within the context of

Prefer `use` to `utilize` or `leverage`. Prefer `is` to `serves as`.
Prefer the name of the operation to `streamline`, `enhance`, or
`facilitate`. These generated-sounding phrases are also banned unless
they are literal names or technically necessary:

- "It is important to note", "It is worth noting", "This is important
  because", and "A few points should be noted";
- "plays a crucial role", "serves as", "stands as", "acts as", and
  "represents a key step" when the plain verb `is`, `uses`, `checks`,
  or `runs` says the fact;
- `rich`, `vibrant`, `pivotal`, `transformative`, `cutting-edge`,
  `state-of-the-art`, `multifaceted`, `nuanced`, `evolving landscape`,
  `ecosystem`, `realm`, `tapestry`, and `cornerstone` as praise;
- "valuable insight", "meaningful contribution", "important
  implication", "deeper understanding", "comprehensive perspective",
  and "future work should explore";
- `unlock`, `foster`, `enhance`, `streamline`, `empower`,
  `facilitate`, `harness`, and `unveil` when a concrete verb names
  the operation;
- `Moreover`, `Furthermore`, `Additionally`, `Notably`, `Importantly`,
  `Consequently`, `Taken together`, and `Overall` when they merely
  make a paragraph sound connected.

Do not replace one banned phrase with a synonym that performs the same
empty job. Delete the decoration or write the missing fact.

#### Sentence-level patterns

Look for repeated patterns across a changed section: many
medium-length sentences with the same shape; paragraphs sharing an
announce-explain-caveat-summary arc; repeated "This suggests",
"While X, Y", or "By doing X, we Y" openings; repeated endings such
as ", highlighting" or ", underscoring"; abstract nouns where a
person or program could act; roadmap sentences that announce, perform,
and recap a simple step.

Check the em dash. Inside a sentence it is the punctuation a reader
notices first, and a passage that reaches for it repeatedly reads as
machine-written rhythm. Use the mark that carries the join: a colon
when what follows restates the clause before it, a comma for an
appositive, parentheses for a true aside, a full stop when two
statements were welded together. The em dash keeps one job:
separating a label from its body in a list entry, table cell, or
heading. More than one sentence-interrupting dash in a changed
section requires review, of the unspaced `word—word` form as well as
the spaced one.

Check these sentence shapes explicitly (three or more rhetorical uses
in one section require a rewrite as direct statements; `thereby` is
banned outright):

1. `This suggests/indicates/demonstrates/highlights that ...`
2. `While X, Y ...`
3. `By doing X, we Y ...`
4. `X, thereby Y ...`
5. `X, highlighting/indicating/reinforcing Y ...`
6. `Although X, it is important to note Y ...`
7. `Not only X, but also Y ...`
8. `X does not merely do A; it also does B ...`

Check abstract-noun stacking: "The validation of the configuration
enables the identification of the availability of the route" hides
every actor; write "The watcher checks the configuration and reports
whether the route is available." Four or more nearby nouns ending in
`-tion`, `-ment`, `-ity`, `-ance`, `-ence`, `-ness`, or `-ization`
require review. Check hedge stacking: two or more of `may`, `might`,
`could`, `appears`, `potentially`, `somewhat` in one sentence usually
mean the actual condition is unstated; prefer "This happens only when
X" or "This case is not tested." Check causal words: `therefore`,
`thus`, `hence`, `consequently` must connect facts that actually
prove the conclusion.

#### Paragraph-level and formatting patterns

Generated paragraphs repeat a five-part arc: announce, explain
broadly, add one detail, add a balanced caveat, summarize. Short,
uneven paragraphs are normal in a README; repeating the polished arc
across a section buries the action. Do not add a roadmap sentence
before every table or example ("The next section explores...", "We
now examine..."). Do not announce, show, and then recap a command
whose result is already clear; use the space to say where to run it,
what it changes, and what the output means.

Keep three items when there are exactly three real items; do not
reshape two or four facts into three for rhythm. In FAQ appendices,
answer the heading directly; no "What does this mean? The answer lies
in...". Bold lead-ins, colons, semicolons, parentheses, passive
voice, and complete sentences are not AI evidence by themselves;
review repetition and reader cost. Do not rotate through synonyms to
avoid repeating the correct term: if the code calls it a watcher,
keep calling it a watcher.

### 3.7 Verification after every README edit

Run all that apply and report the results:

1. Render check (markdown-it, commonmark preset with tables enabled):
   the file must render; tables produce `<table>`, notes produce
   `<blockquote>`; count them if the edit touched one.
2. Quote sync: every `[Adapted from ...]` block touched by the change
   matches its source file (Section 3.2).
3. Link/anchor check: every `#anchor` used in the TOC or cross-references
   has a matching `<a name="anchor">`.
4. URL check (after tags exist on GitHub):
   `curl -s -o /dev/null -w "%{http_code}" <raw URL>` must print 200.
5. GitHub alert blocks use `> [!NOTE]`, `> [!TIP]`, `> [!Warning]` syntax.

## 4. Release checklist

Before a Cocoa tag is created, verify in `set_installation_options.sh`:

1. Debug/dev keys commented: `COCOA_OUTPUT_VERBOSE`, `COCOA_OUTPUT_DEBUG`,
   `COSMOLIKE_DEBUG_MODE`, `SWITCH_TO_DEV_MODE`.
2. Every dependency pinned (no repository cloning a floating branch); each
   pinned tag verified to exist with `git ls-remote --tags`.
3. Project tags point at the versions containing this release's features.
4. README pointer versions bumped (Section 3.4); quote blocks synced
   (Section 3.2); Docker `COCOA_TAG` default matches.
5. After the tag is pushed: the `refs/tags/<tag>/cocoapy311-linux.yml` raw
   URL returns 200.

## 5. Known traps (each of these happened; do not repeat them)

- Calling `unset_all` after a helper that already called `error` prints
  `unset_all: command not found` (Section 2.3).
- `raise`-style aborts hidden in helpers: in this repo's Python-adjacent
  tooling (cobaya theory blocks), rejecting a sample must not abort the run —
  but that code lives in other repositories; here, the analogue is: never
  let a script `exit`.
- An active `XXX_GIT_BRANCH` silently overrides `XXX_GIT_TAG` (Section 2.6).
- pip "upgrading" numpy/scipy into `.local` shadows the conda pins and
  breaks numba/GPy later (Section 2.8).
- Compiler symlink names: the conda-forge triplet is
  `x86_64-conda-linux-gnu-*`. The old `x86_64-conda_cos6-linux-gnu-*` names
  are retired; any instruction using them creates dangling symlinks.
- Filenames referencing another project (for example a `filename_baryon_pca`
  pointing at a different project's chains folder) survive copy-based project
  creation; grep for the source project's name after any copy.
- Tags must share the `v` prefix style; a tag named `4.X.Y` (no `v`) will be
  missed by `v*` globs in scripts and searches.
- NumPy 2.0: any pip command that can resolve numpy without an explicit
  `numpy==${COCOA_NUMPY_VERSION}` in the same command may pull numpy 2.x
  into `.local` and break the compiled stack. NumPy stays at 1.26.x maximum,
  forever (Section 2.11).
- macOS arm compile failure — mpicxx dies with
  `arm64-apple-darwin20.0.0-clang++: No such file or directory`: conda-forge
  `cxx-compiler` 2.0.0 dropped the `clangxx_osx-arm64` triplet shim that
  the MPI compiler wrappers expect. Fix: pin `c-compiler`, `cxx-compiler`,
  and `fortran-compiler` to 1.11.0 in the `cocoapy311-osxarm*.yml` files.
  Version fixes go into the base/loose ymls AND the lock file is regenerated
  with conda-lock (`-p osx-arm64`); never hand-edit a lock file.
- Key-name mismatch between options file and script: if
  `set_installation_options.sh` exports `XXX_THEORY_URL` but the setup
  script reads `${XXX_URL:-<default>}`, the default URL silently wins, and
  a pinned commit that exists only in the intended fork fails later with
  `fatal: unable to read tree <commit>`. When touching a pin trio, grep the
  consuming script for the exact key names it actually reads.
- Emulators fail silently outside their training box: out-of-range inputs
  (for example baccoemu with `omega_baryon` below its box) return -inf/NaN
  with no exception. The theory block must validate the box and reject the
  sample (Section 7, step 8).
- A guard copied from the neighboring block: the bfmt symlink blocks in
  `start_cocoa.sh`/`stop_cocoa.sh` shipped guarded by `IGNORE_FASTPT_CODE`,
  copy-pasted from the FAST-PT block above them. After copying any guarded
  block, grep the new block for the donor's key.
- Git LFS repos: branch operations fail with "This repository is configured
  for Git LFS but 'git-lfs' was not found on your path" when the cocoa
  conda environment is not active. Activate it (it provides `git-lfs`)
  before any git work in project repositories.
- `~/.condarc` edits: steps that temporarily change conda channel settings
  (the channel allowlist) must save a copy first and restore the file
  byte-for-byte afterwards (verify with `diff`) — never reconstruct it from
  memory.
- Validating cobaya yamls: they may contain the `!defaults` tag, which
  plain `yaml.safe_load` rejects. Use cobaya's own loader (or register a
  tolerant constructor) when checking that a yaml parses.
- Emulator first-use downloads: some packages fetch their model files at
  the first evaluation; on an offline compute node the first MCMC step then
  dies mid-run. Trigger the download in the package's `setup_*.sh`
  (Section 2.7), never leave it to run time.
- Notebooks lag interface upgrades: the EXAMPLE_EVALUATE notebooks call the
  compiled cosmolike bindings directly (their own CAMB run, `set_cosmology`,
  `compute_data_vector_masked`), a path the cobaya pipeline never exercises.
  After any binding change (`init_IA` growing `ia_code`) or C-side
  validation change, sweep every project's `*.ipynb` for the old call and
  run the project's `tests/test_notebook_interface.py`, which rebuilds the
  notebook call sequence against the frozen data.
- Cosmolike grid validation aborts the process, not the call: a
  non-monotone grid handed to `set_cosmology` (for example a chi(z) grid
  concatenated from linspace segments that each include both endpoints,
  duplicating the junction values) hits `log_fatal` + `exit` in basics.c.
  Inside Jupyter that is a dead kernel with no traceback; the message
  ("chi z: not strictly increasing at i=N") only appears on a terminal.
  Reproduce notebook crashes headlessly (extract the cells into a script
  and run it in a subprocess) before guessing. Interior linspace segments
  take `endpoint=False`; the junction value belongs to one segment only.

## 6. Bash style guide (observed across all installation_scripts)

When writing or editing a script, imitate these conventions exactly. They
hold across the ~85 scripts in `installation_scripts/`.

### 6.1 Layout

- `#!/bin/bash` first line (even though scripts are sourced).
- Sections separated by full-width divider comments:

      # ------------------------------------------------------------------------------
      # SECTION TITLE IN WORDS ------------------------------------------------------
      # ------------------------------------------------------------------------------

- Order inside a script: IGNORE-key guard, ROOTDIR guard, flags_check in a
  subshell, unset functions, error/helper functions, body, `unset_all`,
  `return 55`.

### 6.2 Naming

- Environmental keys and script-level variables: UPPER_SNAKE
  (`PACKDIR`, `ECODEF`, `CCIL`, `PLIB`, `URL`, `FOLDER`, `TMP`, `TMP2`).
- `PRINTNAME` holds the banner text passed to `ptop`/`pbottom`.
- Functions: lowercase (`cdfolder`, `cpfolder`, `cpfile`, `gitact0`,
  `wgetact`, `devurl`, `error`, `unset_all`).
- Abbreviation comments are welcome where a name is dense:
  `# E = EXTERNAL, CODE, F=FODLER`.

### 6.3 Quoting and parameter expansion

- Must-exist expansion everywhere a wrong value would be destructive:
  `"${VAR:?}"`. Every `rm` path uses it (Section 2.4).
- Defaults with fallback: `"${XXX_URL:-"https://github.com/..."}"`.
- Optional test: `[ -n "${KEY:-}" ]` / `[ -z "${KEY:-}" ]` — always with
  the `:-` so `set -u` (debug mode) does not break.
- Case conversion via expansion, not `tr`: `${VAR,,}` (lower), `${VAR^^}`.

### 6.4 Output and error text

- Commands are silenced by appending `>>${OUT1:?} 2>>${OUT2:?}`; OUT1/OUT2
  are chosen by the verbosity keys, so never hardcode `/dev/null`.
- User-visible progress uses `ptop "DOING X"` / `pbottom "DOING X"` pairs
  with identical text.
- Error messages come from the `EC<N>` catalog exported in
  `flags_derived.sh` (for example `EC15="GIT CLONE"`,
  `EC34="SYMLINK CREATION FAILED"`). To add a new error string, add the next
  `EC<N>` there and reference it as `"${EC<N>:?}"`; do not inline new
  free-text messages in scripts when a code fits.
- No `set -e` in scripts: every command's failure is handled explicitly with
  `|| { ...; return 1; }`. (`COCOA_OUTPUT_DEBUG=1` turns on strict mode
  externally; scripts must still work without it.)

### 6.5 Idempotence patterns (scripts run twice safely)

- Clone guard: `if [ ! -d "${PACKDIR:?}" ]; then git clone ...; fi`,
  preceded by `if [ -n "${OVERWRITE_EXISTING_XXX_CODE:-}" ]; then rm -rf ...`.
- Symlink guard: `if [[ ! -L "${LINK}" ]]; then ln -s ...; fi` on create,
  `if [[ -L "${LINK}" ]]; then rm -f ...; fi` on remove.
- Expensive pip stages use a sentinel file whose name contains a hash of the
  inputs (see `PIPCP_HASH` in `setup_pip_core_packages.sh`); changing the
  inputs invalidates the sentinel automatically.
- Cleanup `rm` lines append `2>/dev/null` so missing files are not errors.

### 6.6 Loops, arrays, platform

- Arrays: `declare -a NAME=("item1"` ... one item per line ... `)`.
  Script lists in `setup_cocoa.sh`/`compile_cocoa.sh` follow this form.
- Index loops are C-style: `for (( i=0; i<${#TMP[@]}; i++ ))`; use
  glob loops (`for f in ...`) when no index is needed.
- Platform switches: `case "$(uname -s)" in Linux) ... ;; Darwin) ... ;; esac`.
  Linux and macOS variants stay as separate explicit blocks — do not merge
  them with clever conditionals.
- `sed` is GNU sed from the conda environment:
  `sed --in-place --regexp-extended 's@old@new@g'` (note `@` delimiters when
  paths contain `/`).

## 7. Integrating an external code as a Cobaya theory block (the bfmt case)

This is the complete checklist for adding a new theory block, in the order
the pieces were built for `bfmt` (the baryonic feedback block). A theory
block has two kinds of ingredients: the theory-block repository (a Python
`Theory` class for Cobaya) and zero or more external emulator/model codes it
imports. For each numbered item, the bfmt example is named so you can open
the real files and copy their shape.

**Step 1 — keys in `set_installation_options.sh`:**
- One `IGNORE_XXX_CODE` per repository (commented = installed by default).
  bfmt: `IGNORE_BFMT_CODE`, plus `IGNORE_PYSPK_CODE`, `IGNORE_BCEMU_CODE`,
  `IGNORE_FBRE_CODE`, `IGNORE_BACCOEMU_CODE` for its emulators.
- One URL/NAME/pin trio per repository (Section 2.6). Pin third-party codes
  by COMMIT; pin the team-owned theory repo by TAG.
  bfmt: `BFMT_THEORY_URL`, `BFMT_NAME="baryon_suppression"`, `BFMT_GIT_TAG`.

**Step 2 — flag bookkeeping:**
- `flags_impl_unset_keys.sh`: `unset -v` for every new key
  (URL, NAME, GIT_COMMIT, GIT_TAG, GIT_BRANCH, IGNORE, OVERWRITE).
- `flags_derived.sh`: add `OVERWRITE_EXISTING_XXX_CODE=1` inside the
  `OVERWRITE_EXISTING_ALL_PACKAGES` cascade.

**Step 3 — setup scripts (one per repository, Section 2.2 skeleton):**
- Theory repo: clone at the pin into `external_modules/code/${XXX_NAME}`.
  bfmt: `setup_bfmt.sh`.
- Each emulator code: clone at the pin. If the code needs a source patch
  (for example removing an unwanted import that hijacks another module),
  apply it in setup with an idempotent `sed`, placed AFTER the clone block
  so re-runs also patch existing clones.
  bfmt: `setup_bcemu.sh` comments out BCemu's `from .spectra import ...`
  (it imported camb at package load, shadowing Cobaya's path-checked CAMB).
- The internet invariant decides where pip install goes (Section 2.7):
  if the package downloads model files, BOTH the pip install and a download
  trigger go in setup (bfmt: `setup_bcemu.sh`, `setup_baccoemu.sh`);
  if it is fully offline-installable, pip goes in a `compile_XXX.sh`
  (bfmt: `compile_pyspk.sh`, `compile_fbre.sh`).

**Step 4 — wire the runner lists:**
- Add each `setup_XXX.sh` to the matching array in `setup_cocoa.sh` and each
  `compile_XXX.sh` to `compile_cocoa.sh`. A repository whose setup does
  everything has NO compile script and no compile entry.

**Step 5 — Python runtime dependencies:**
- Because compile-time pip uses `--no-dependencies`, every runtime import of
  the new codes must be seeded in `PIPCP`
  (`setup_pip_core_packages.sh`) with a pinned version, after the dry-run
  check of Section 2.8. Find the imports with:

      grep -rhE "^(import|from) [a-zA-Z0-9_]+" <pkg>/*.py | sort -u

  bfmt additions: `pydantic` (pyspk), `smt==1.0.0` + `wget` + `msgpack`
  (BCemu), `swiftemulator` (FBRE), `progressbar2` (baccoemu).

**Step 6 — the symlink into Cobaya (start/stop pair):**
Cobaya finds a theory class at `cobaya/cobaya/theories/<blockname>/`.
Cocoa does not copy files there; `start_cocoa.sh` creates a symlink and
`stop_cocoa.sh` removes it. Add one block to EACH file, guarded by the
block's OWN IGNORE key (a copy-pasted guard from the neighboring block is a
real bug that happened):

```bash
# in start_cocoa.sh
if [[ -z "${IGNORE_BFMT_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
  TMP="${BFMT_NAME:-"baryon_suppression"}"
  TMP2="bfmt"   # the name Cobaya sees: theory block `bfmt`
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}" "${COBTH:?}/${TMP2}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi
  unset -v ECODEF COBTH TMP TMP2
fi

# in stop_cocoa.sh
if [[ -z "${IGNORE_BFMT_CODE}" ]]; then
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
  TMP="bfmt"
  if [[ -L "${COBTH:?}/${TMP}" ]]; then
    rm -f "${COBTH:?}/${TMP:?}"
  fi
  unset -v COBTH TMP
fi
```

The repository must contain `<blockname>.py` defining
`class <blockname>(Theory)` at its top level for the symlinked folder to
work as a Cobaya theory package.

**Step 7 — housekeeping:**
- Add every clone destination (`external_modules/code/<name>`) to
  `Cocoa/.gitignore` and verify with `git check-ignore` (Section 2.12 —
  otherwise `git add --all` stages the clone as an accidental submodule).
- Document the block: keys quote + usage in the theory repository's README,
  and a short section in the consuming project READMEs pointing to it.

**Step 8 — theory-block code expectations (lives in the theory repo, not
here, but reviewers should check):**
- Sample rejection returns `False` from `calculate()` after logging a
  warning. NEVER `raise LoggedError` for a per-sample rejection: LoggedError
  is in Cobaya's `always_stop_exceptions` and kills the whole run.
- Validate sampled parameters and emulator training boxes; out-of-range
  redshifts degrade to S=1 or a clamped boundary value, never silent
  extrapolation.
- Load emulator files ONCE in `initialize()`; per-sample updates go through
  the emulator's public per-call API.
- Compute on a small internal grid and 2D-spline onto the grid the
  likelihood requests (the Cobaya matter-power-interpolator strategy), so
  cost does not scale with the likelihood's interpolation grid.

**Verification for the whole integration:** `bash -n` every touched script;
`source setup_XXX.sh` end-to-end in the cocoa environment; run a cobaya
evaluate with the block enabled and once with it disabled, and compare.

## 8. Python scripts: style and documentation contract

This contract applies to every Python file committed to a Cocoa
repository: tests, maintenance
scripts, generators, drivers. The intended reader is a physics student
who understands C-like control flow but may not know advanced Python
idioms.

### 8.1 Code shape

- Explicit, C-like control flow. One consequential operation per line.
  An `if` block over a clever expression; an explicit loop over a
  comprehension whenever the loop validates, mutates, logs, or does
  more than one transformation.
- Stage operations into named intermediate variables. NO-GO:

      training = torch.from_numpy(rows[idx].astype("float32")).to(dev)

  GO (each step can fail on its own line and be named in an error):

      training_rows = rows[idx].astype("float32")
      training = torch.from_numpy(training_rows)
      training = training.to(device=dev)

- Use named (keyword) arguments whenever the callee accepts names.
- A dictionary or list with three or more entries puts one item per
  line.
- Forbidden when a plainer form exists: the walrus operator `:=`,
  nested comprehensions, chained ternaries, a `lambda` where a named
  function is clearer, starred unpacking that hides value order, and
  monkey patching (replacing a function/method/module attribute at
  run time) anywhere, including tests.
- Never de-vectorize a numerical hot path for readability: vectorized
  numpy stays vectorized, with a comment giving the mathematical
  reason or shape invariant instead.
- Validate inputs BEFORE any mutation, file write, or expensive setup.
  A failure message states: what failed, the observed value, the
  required condition, and the corrective action. "invalid input" is
  never enough. No silent fallback, no silent coercion.
- Keep lines within 90 columns; prefer parentheses over backslash
  continuation; match the indentation style of the file being edited.

### 8.2 Documentation density (the activations.py standard)

Roughly half of every file is explanation, written for someone who
knows the physics but not this code. This sample shows the expected
density; imitate its shape (a constant with the meaning of its value,
a docstring teaching the mechanism and the reason, an Arguments block
in `name = description` form, and inline comments that state
invariants rather than narrate lines):

```python
# The race tests must run multi-threaded: with one thread there is no
# thread scheduling, so an OpenMP race could never show up.
REQUIRED_OMP_THREADS = "4"


def verify_frozen():
    """Fail every test up front when the frozen state was edited.

    Compares the stored manifest with a fresh hash of tests/frozen/ in
    both directions, so an edited file (CHANGED), a deleted file
    (MISSING), and a new file (EXTRA) are all reported. This runs
    before any model is built: a tampered frozen state must not
    produce a plausible-looking chi2.

    Returns:
      nothing when every frozen file matches the manifest.

    Raises:
      AssertionError listing every mismatched path and pointing to
      generate_frozen_reference.py --overwrite for a deliberate
      refresh.
    """
    # expected = the {relative path: sha256 digest} table written at
    # freeze time; it is the definition of "untouched"
    with open(MANIFEST_FILE) as f:
        expected = json.load(f)["files"]
    # actual = the same table computed from the files on disk right
    # now (compute_manifest walks frozen/ and fingerprints each file
    # with SHA-256, skipping only python bytecode caches)
    actual = compute_manifest()
    # collect every discrepancy before raising: a report naming all
    # problem files at once beats failing on the first one
    problems = []
    for rel, digest in expected.items():
        if rel not in actual:
            # the manifest lists it but the file is gone from disk
            problems.append(f"MISSING  {rel}")
        elif actual[rel] != digest:
            # the file exists but at least one byte differs
            problems.append(f"CHANGED  {rel}")
    # both directions matter: a file ADDED to frozen/ is as suspicious
    # as an edited one, so the reverse scan runs too
    for rel in actual:
        if rel not in expected:
            problems.append(f"EXTRA    {rel}")
    if problems:
        raise AssertionError(
            "Frozen test data does not match the manifest:\n  "
            + "\n  ".join(problems))
```

At a call site, say what the called helper produces, not just its
name: the reader should not need to open another file to follow the
flow.

Concretely:

- **Module docstring** teaches the domain first: what the file
  computes, the definition of every non-obvious term it relies on
  (assume the reader has NOT seen the framework before), how the
  pieces relate, and how to run it when it is runnable.
- **Every function, method, and class gets a docstring** with:
  - a first sentence containing a subject and a verb;
  - a paragraph explaining the mechanism and WHY it is built this way
    (the non-obvious decision, the failure it prevents);
  - an `Arguments:` block naming every parameter as
    `name = what it is, units/shape/valid range when relevant`;
  - a `Returns:` block (type, shape, units);
  - a `Raises:` block when the function refuses inputs;
  - side effects (files written, directories changed, state mutated).
- **Comments state reasons, invariants, units, and failure
  boundaries — never the next line.** NO-GO: `# add one to counter`.
  GO: `# count accepted rows only; rejected rows must not shift
  checkpoint indices`.
- **Comments teach; jargon without its lesson is banned.** A comment
  must not lean on project-internal vocabulary — drift, descriptor,
  frozen, fiducial, coverage — unless that comment (or the docstring
  it sits under) spells the concept out in plain words at the point
  of use. Write for a physics student opening the file for the first
  time. Calibrate "knows python" carefully: assume BASIC working
  knowledge with gaps, the way a C programmer reads python — they
  will not remember what an idiom does on the fly, and the comments
  must save them the trip to the python documentation. Whenever a
  line uses a non-obvious python construct (splitlines with
  keepends, a comprehension with a condition, a generator expression
  inside sum, functools.partial, a with-statement's cleanup
  guarantee, star-unpacking, str.format field syntax), the comment
  at that line says in plain words what the construct produces.
  Prefer the plainer construct when it reads better in C terms: an
  explicit loop that a comment can narrate beats a nested one-liner
  that needs a paragraph.
  Corrected failures, kept as calibration: "the same drift condition
  build_point checks" named neither the condition nor what drifts
  (say instead: the drawn point and today's model must name the same
  sampled parameters; a mismatch means the code gained or lost a
  parameter since the freeze). "the new descriptor is the frozen
  descriptor text with only its data_file line replaced" assumed the
  reader knows a descriptor (say instead: the ".dataset" file, the
  small text file listing which data files the likelihood reads, one
  `key = filename` line each, rewritten with one line changed).
  Repeating a lesson where the concept reappears beats defining it
  once far away.
- **Every constant** carries a comment with the meaning of the chosen
  value (why 0.2, why 4 threads, why this list of nine cosmologies).
- Python prose follows the anti-AI rules of Section 3.6 in full.

### 8.3 Scope discipline

- A narrow bug gets a narrow fix: no new registry, framework, or
  validation subsystem where a short direct check repairs the named
  problem.
- Do not turn one task into a repository-wide cleanup; note other
  problem sites and report them instead.
- Tests may be longer than the code they test (they show valid and
  invalid cases), but they follow every rule above, including full
  docstrings.

## 9. Interpreting the project accuracy tests

Every Cosmolike project carries `tests/test_accuracy.py`. Its output
is advisory: no assertion fails on a large delta chi2, because how
much numerical error an analysis tolerates is a judgment call. This
section says how to read the numbers and what to change when they
move.

### 9.1 What the output is

The file runs two kinds of checks on frozen configurations:

- `KNOB` lines — a one-knob-at-a-time scan on the 3x2pt NLA
  configuration: each accuracy knob is raised alone and the resulting
  chi2 is compared with the frozen default-settings reference.
- A1-A6 — all knobs raised at once, one check per probe (cosmic
  shear, 2x2pt, 3x2pt) and IA model (NLA, TATT). The all-knobs set
  compares the default cosmolike accuracyboost 1 against 3, the
  highest value that stays healthy in every project scanned; the
  KNOB scan keeps 5 as a deliberate stress entry.

Every reference is frozen and sits at the chi2 minimum, because the
data vectors are synthetic: the model generated them at the frozen
point. At a minimum the chi2 responds quadratically to a settings
change, so the reported deltas are stable and small. The comfort
target is |delta chi2| below 0.2. A larger delta is not a failure; it
is a number to investigate.

### 9.2 Investigation order when knobs move the chi2

Settle the cheap knobs before blaming the expensive one:

1. cosmolike `accuracyboost` first (cheap at run time);
2. camb `k_per_logint` second;
3. camb `AccuracyBoost` last. It is expensive, and an apparent CAMB
   sensitivity can masquerade as unresolved cheap-knob resolution: in
   roman_kl an apparent +0.80 from camb `AccuracyBoost` collapsed to
   +0.002 once `k_per_logint` reached 50.

`kmax_boltzmann` (cosmolike) and camb `kmax` are one physical cutoff
seen from two sides; move them together, never one alone.

### 9.3 Reading a knob's convergence scan

To decide what a delta means, scan the knob: evaluate the chi2 at
increasing knob values, always against the vector generated at the
default settings. Three shapes occur:

- Monotone rise to a plateau — the plateau is the real numerical
  error of the default. The roman_kl camb `k_per_logint` scan
  plateaus at 0.255 by 25; raising the default removes exactly that
  error.
- Explosion far beyond every other knob — an interface breakdown,
  not a refinement. The desy1xplanck `accuracyboost` scan is fine
  through 3.25, then +0.26 at 4 and +27 at 5. Cap the knob and file
  an interface bug.
- Oscillation with no plateau — a re-phasing interpolation grid, not
  a convergence property. The roman_kl `accuracyboost` scan jumped by
  0.3 to 3 between boosts 1.25 and 6 because the z-node count of the
  power-spectrum tables grew additively, re-phasing the linear-
  interpolation sawtooth at every boost. Do not chase such a shape by
  raising defaults; fix the grid. The cure, now in every project's
  likelihood prototype: the boost refines the z grid dyadically
  (nested nodes, boost-1 grid unchanged, CAMB's 256-redshift cap
  respected by a fixed request grid), after which the same scan is
  monotone below 0.007.

### 9.4 What a measurement changes

- Fix cheap knobs in the example yamls, with a comment stating the
  measured numbers (the scanned values, the plateau, the resulting
  delta).
- Never bake a camb `AccuracyBoost` increase into a yaml without
  re-measuring at the raised cheap knobs; the roman_kl case in
  Section 9.2 shows the expensive knob absorbing blame that belonged
  to a cheap one.
- Document caps in the likelihood yamls, next to the knob they limit,
  so a user raising the knob past the cap finds the warning where
  they type.

### 9.5 Fixing the cosmolike C core

A fix inside external_modules/code/cosmolike_core follows two rules
at once:

- **Never propose cubic interpolation for the z axis of the
  cosmolike hot-path tables** (p_lin, p_nonlin, growth): it is too
  expensive at run time. Interpolation-grid problems are fixed
  through the node count and node placement handed in by the
  likelihood prototypes, not by raising the interpolation order.
- **Smallest possible diff.** Change only the lines the demonstrated
  failure requires; never restructure around the fix. The FFTLog
  padding repair is the model: three constants became three scaled
  expressions, nothing else moved.
- **Didactic to a physics student.** The comment above the changed
  lines must let a student who has never seen the file understand the
  numerical mechanism: what the quantity is, why the old form failed
  (with the measured number), and why the new form is right. A fix a
  student cannot learn from is not finished.

Work on a `bugfix` branch of the core repository (check
`git rev-parse --show-toplevel` and the remote first: the core is a
separate repository pinned by the COSMOLIKE keys), commit, never
push, and verify that results at the unboosted defaults are unchanged
(the frozen project references must not move) before claiming the
fix.

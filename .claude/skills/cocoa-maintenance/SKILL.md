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
   After editing a README run the render check (Section 3.6). Report the
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
   alphabetical position (functions go under `# Functions` as `unset -f`).
3. `installation_scripts/flags_derived.sh` — only if the key participates in
   a cascade (for example each `OVERWRITE_EXISTING_XXX_CODE` is set to 1 when
   `OVERWRITE_EXISTING_ALL_PACKAGES` is set).
4. If a new setup/compile script consumes it: add the script name to the
   matching array (`CORE`/`THEORY`/`ML`/`LIKELIHOOD`) in `setup_cocoa.sh`
   and/or `compile_cocoa.sh`. Adding or removing a list entry is safe: the
   cache file rebuilds automatically when the list length changes.

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

## 3. README rules

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
- One topic = one TIP (a tip may contain two numbered solutions; do not
  split into two tips).

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

### 3.6 Verification after every README edit

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
  `.gitignore`.
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

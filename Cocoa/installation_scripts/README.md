# The installation scripts

This folder holds the scripts that download, build, and register every
external package Cocoa uses. Two runners one folder up, `setup_cocoa.sh` and
`compile_cocoa.sh`, source them in a fixed order during a full installation;
each script can also be sourced on its own to redo one package.

The file names follow five families:

- `setup_*.sh` — download one package (git clone, wget). Setup scripts MAY
  use the internet.
- `compile_*.sh` — build and install what the matching setup script
  downloaded. Compile scripts NEVER use the internet, because they may run
  on compute nodes without network access.
- `unxv_*.sh` — unpack one dataset archive into `external_modules/data`
  (the body runs `tar xf`).
- `flags_*.sh` — bookkeeping: `flags_check.sh` confirms the required
  variables exist, `flags_derived.sh` computes derived keys and the error
  catalog, `flags_save_old.sh` / `flags_recover_old.sh` save and restore
  the system variables Cocoa modifies, and `flags_impl_unset_keys.sh`
  unsets every key Cocoa invents.
- project helpers — `create_link_all_projects.sh` /
  `delete_link_all_projects.sh` create and remove the symlinks (filesystem
  entries that point at another path) exposing each project inside cobaya;
  `start_all_projects.sh` / `stop_all_projects.sh` run every project's own
  start and stop script.

Every script here is SOURCED, never executed. `source script.sh` runs the
file's commands inside the current shell, so the script can create, change,
and delete variables in that shell. `bash script.sh` runs the commands in a
child process instead; every variable it sets vanishes when the child ends.
Cocoa scripts must export paths and keys that outlive them, and only
`source` can do that.

A "key" is an exported shell variable used as a setting. For example,
`export IGNORE_PYSPK_CODE=1` tells every pyspk script to skip itself. The
user-facing keys live in `set_installation_options.sh` (one folder up).

> [!NOTE]
> Because the scripts run inside the user's shell, they end with `return`,
> never `exit`. See [FAQ: Why must these scripts never call exit?](#faq-no-exit).

## Anatomy of one script: setup_pyspk.sh

`setup_pyspk.sh` downloads pyspk, a baryonic feedback model. It is one of
the shortest scripts in the folder, and every other `setup_*.sh` repeats its
structure. Reading it top to bottom teaches the whole family.

### The IGNORE guard

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
#!/bin/bash
# ------------------------------------------------------------------------------
# setup_pyspk.sh: download pyspk (the SP(k) baryonic feedback model).
#
# Sourced, never executed (only `source` keeps its environment changes).
(...)
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_PYSPK_CODE:-}" ]; then
  return 99
fi
```

`[ -n "..." ]` tests "is this string non-empty". The suffix `:-` inside
`${IGNORE_PYSPK_CODE:-}` means "expand to empty when the variable is unset",
so the test also works under bash's strict mode, where reading an unset
variable is an error. When the user has set the IGNORE key, the script
returns 99, the agreed code for "skipped on purpose". The runner treats 99
as neither success nor failure.

### The ROOTDIR guard and flags_check

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
```

`ROOTDIR` is the absolute path of the `Cocoa/` folder; `start_cocoa.sh`
exports it along with every other key. `[ -z ... ]` tests "is this string
empty": when `ROOTDIR` is missing, the user sourced this script directly in
a fresh shell, so the script starts a Cocoa session itself. This guard is
what makes single-package reinstalls work, and it is why a direct source
must happen from the `Cocoa/` folder (the path `start_cocoa.sh` is relative).

The parentheses on the next line create a subshell: a child copy of the
shell. `flags_check.sh` runs there, so the helper function it defines never
reaches the caller; only its return status comes back. The check confirms
that `PYTHON3`, `PIP3`, the compiler variables, `CMAKE`, and `GIT` are all
defined before any real work starts.

### The cleanup trio

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR  
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error cpfile
  unset -f unset_env_funcs
  cdroot || return 1;
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}
```

This is the price of sourcing: every variable (`unset -v`) and every
function (`unset -f`) the script defines stays in the user's shell unless
the script removes it before returning. Each script therefore carries its
own list of names to forget. If a new variable or function is added to the
body, its name must be added to the matching list. `unset_all` runs both
lists and then deletes itself, so after cleanup no trace of the script
remains.

### The error contract

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}
```

`error()` does two jobs at once: it prints which step failed, and it runs
the full cleanup. That single fact fixes the shape of every failure path in
the body:

- a helper whose own failure path calls `error()` (such as `cdfolder`
  below) is called with plain `|| return 1`, because by the time control
  returns, the cleanup already ran;
- a direct command (git, pip, tar) is called with
  `|| { error "MESSAGE"; return 1; }`, because nothing else will report or
  clean up for it;
- `ptop`/`pbottom` (the banner printers, described below) keep
  `|| { unset_all; return 1; }`, because they do not route through
  `error()`.

Mixing these up has a visible symptom: calling `unset_all` after a helper
that already called `error()` prints `unset_all: command not found`, since
the first cleanup deleted the function.

> [!TIP]
> If that message appears while editing a script, see
> [FAQ: What does `unset_all: command not found` mean?](#faq-unset-all-twice).

### The body: names first

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
unset_env_vars || return 1

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

URL="${PYSPK_URL:-"https://github.com/jemme07/pyspk.git"}"

FOLDER="${PYSPK_NAME:-"pyspk"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"
ptop "INSTALLING BARYONIC PYSPK FEEDBACK MODELS" || { unset_all; return 1; }
```

The body starts by clearing its own names (a leftover from an earlier run
must not survive) and then defines every path it will use. Two parameter
expansions carry the safety rules of the whole folder:

- `${VAR:?}` expands to the value of `VAR`, but if `VAR` is unset or empty
  the command fails with an error message instead of expanding. Every
  destructive command uses it: `rm -rf "${PACKDIR:?}"` can never become
  `rm -rf /` through an unset variable, because the expansion aborts first.
- `${PYSPK_URL:-"https://..."}` expands to `PYSPK_URL` when the user set
  it, and to the written default otherwise. The user can point a package at
  a fork by exporting one key, without touching the script.

`ptop` prints a colored banner announcing the step; `pbottom` later prints
the matching "DONE" banner with identical text. Both come from
`flags_derived.sh`.

### Downloading only when needed

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
if [ -n "${OVERWRITE_EXISTING_PYSPK_CODE:-}" ]; then
  rm -rf "${PACKDIR:?}"
fi

if [ ! -d "${PACKDIR:?}" ]; then
  echo "${PACKDIR:?}"

  cdfolder "${ECODEF:?}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
```

The two guards make a second run harmless. Without the OVERWRITE key, an
existing clone is kept and the script falls through to the end; with it,
the old clone is deleted first (`rm -rf` on a `:?`-guarded path) and
downloaded again. `OVERWRITE_EXISTING_PYSPK_CODE` is set for the user by
the `OVERWRITE_EXISTING_ALL_PACKAGES` cascade in `flags_derived.sh`.

The clone line shows two more folder-wide conventions:

- `>>${OUT1:?} 2>>${OUT2:?}` appends the command's normal output to `OUT1`
  and its error output to `OUT2`. Both are chosen by `flags_derived.sh`:

  [Adapted from installation_scripts/flags_derived.sh]

  ```bash
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  ```

  With `COCOA_OUTPUT_VERBOSE` set, they point at the terminal instead.
  Scripts never hardcode `/dev/null`; they route through these variables so
  one key controls the verbosity of every step.

- `"${EC15:?}"` is an entry in the error catalog, a numbered list of
  standard failure messages exported by `flags_derived.sh`:

  [Adapted from installation_scripts/flags_derived.sh]

  ```bash
  export EC15="GIT CLONE"

  export EC16="GIT CHECKOUT"
  ```

  A script that needs a new message adds the next `EC<N>` to the catalog
  and references it as `"${EC<N>:?}"`, instead of inventing free text at
  the call site. The catalog keeps failure output uniform and searchable.

### Pinning: the URL / NAME / version trio

Every external repository is "pinned": fixed to one exact version, so every
Cocoa installation is identical. A pin is three keys in
`set_installation_options.sh`: the URL, the folder NAME, and exactly one of
`GIT_COMMIT` / `GIT_BRANCH` / `GIT_TAG`. Here is the pyspk block in its
shipped state (the IGNORE key ships commented, meaning the package installs
by default):

[Adapted from set_installation_options.sh]

```bash
#export IGNORE_PYSPK_CODE=1
(...)
export PYSPK_URL="https://github.com/jemme07/pyspk.git"
export PYSPK_GIT_COMMIT="50737f9295fee75ef2fe97e4bdd284134ed0d474"
export PYSPK_NAME="pyspk"
```

The setup script applies the pin with an `if`/`elif` chain, and the order
of the chain is the precedence rule: COMMIT beats BRANCH beats TAG. An
accidentally exported `PYSPK_GIT_BRANCH="main"` therefore overrides a TAG
pin and silently unpins the package, which is why exactly one of the three
should ever be active.

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
  if [ -n "${PYSPK_GIT_COMMIT:-}" ]; then

    "${GIT:?}" checkout "${PYSPK_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
(...)
  elif [ -n "${PYSPK_GIT_BRANCH:-}" ]; then
(...)
  elif [ -n "${PYSPK_GIT_TAG:-}" ]; then
(...)
  fi
```

### The ending: return 55

[Adapted from installation_scripts/setup_pyspk.sh]

```bash
cdfolder "${ROOTDIR:?}" || return 1;

pbottom "INSTALLING BARYONIC PYSPK FEEDBACK MODELS" || { unset_all; return 1; }
(...)
unset_all || return 1
(...)
return 55; # why this odd number? Setup_cocoa will cache this installation only
```

The last line of every setup and compile script is `return 55`, the agreed
code for "ran to the very end". The runner records success only on 55, so a
script that stops early (a closed terminal in the middle of a git clone
leaves `PACKDIR` present but corrupted) is run again next time. Any return
value reaching the runner other than 55 or 99 counts as failure.

## The compile side: compile_pyspk.sh

`compile_pyspk.sh` repeats the same skeleton (IGNORE guard, ROOTDIR guard,
flags_check, cleanup trio, `error()`, `return 55`) and differs only in its
body: instead of cloning, it installs the already-downloaded clone into
Cocoa's private Python environment.

[Adapted from installation_scripts/compile_pyspk.sh]

```bash
PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
(...)
rm -rf "${PACKDIR:?}/build/"
rm -rf "${PLIB:?}"/pyspk/
rm -rf "${PLIB:?}"/pyspk-*
(...)
#prevent all compile_XXX.sh from using the internet (run @compute nodes)
#FROM: https://github.com/pypa/pip/issues/12050
#That is why we use --no-dependencies --no-index --no-build-isolation
(
  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install "${PACKDIR:?}" \
      --no-dependencies \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-build-isolation
) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC3:?}"; return 1; }
```

`.local` is a Python virtual environment (a private folder of Python
packages and executables) that `setup_cocoa.sh` creates inside `Cocoa/`.
Installing there, instead of into the conda environment, keeps every Cocoa
package in one deletable folder and lets Cocoa pin versions independently
of conda. The `rm` lines above delete any previous build first, so a
recompile starts clean.

Each pip flag closes one door to the network or to version drift:

- `--no-dependencies`: install only this package; do not resolve or fetch
  its dependencies (they are pre-seeded, pinned, by
  `setup_pip_core_packages.sh`).
- `--prefix="${ROOTDIR:?}/.local"`: install into Cocoa's private
  environment.
- `--no-index`: never contact PyPI.
- `--no-build-isolation`: build with the packages already installed instead
  of downloading a fresh build environment.

> [!TIP]
> For the reasoning behind the offline rule, see
> [FAQ: Why do compile scripts never use the internet?](#faq-compile-offline).

## How a full installation runs

`setup_cocoa.sh` and `compile_cocoa.sh` each hold bash arrays listing the
scripts to run, grouped as `CORE`, `THEORY`, `ML`, `LIKELIHOOD`, and (for
setup) `DATA` and `COSMOLIKE`:

[Adapted from setup_cocoa.sh]

```bash
declare -a THEORY=("setup_hyrec2.sh"
                   "setup_cosmorec.sh"
                   "setup_camb.sh"
(...)
                   "setup_pyspk.sh"
                   "setup_bcemu.sh"
                   "setup_fbre.sh"
                   "setup_baccoemu.sh"
                  )
```

The groups are concatenated into one `SCRIPTS` array, and a loop sources
each entry in a subshell, recording the outcome in a cache: a text file
with one line per script, `1` for finished and `0` for pending
(`.setup_local_packages.txt` in `Cocoa/` for setup,
`.local/compile_local_packages.txt` for compile). The runner rebuilds the
cache from zeros whenever its length stops matching the script list, so
adding or removing a list entry is safe.

[Adapted from setup_cocoa.sh]

```bash
for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
(...)
  [[ "${CACHE[i]}" -eq 1 ]] && continue
(...)
    ( source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" )
(...)
  rc=$?

  if [[ ${rc:?} -eq 55 ]]; then
    CACHE[i]=1
    save_cache
  elif [[ ${rc:?} -eq 99 ]]; then
    continue
```

The three return codes of the anatomy section meet here: `55` marks the
script finished in the cache, `99` (skipped by an IGNORE key) leaves the
cache at `0` without complaint, and anything else reports a failure. The
runners also accept one mode option (`--soft`, `--hard`, `--aggressive`,
and for setup `--extreme` and `--purge`) that resets the cache entries of
whole groups to force a redo.

To redo a single package, skip the runner and source its scripts directly
from the `Cocoa/` folder:

```bash
cd Cocoa
source installation_scripts/setup_pyspk.sh
source installation_scripts/compile_pyspk.sh
```

> [!TIP]
> If a rerun of the runner skips a package you expected it to redo, see
> [FAQ: Why did my second run skip a package?](#faq-cache-skip).

## Environment hygiene

Cocoa's goal: after `source stop_cocoa.sh`, the shell is indistinguishable
from the shell before `source start_cocoa.sh`. Three mechanisms divide the
work, and every variable belongs to exactly one of them.

**1. Per-script unsets.** The cleanup trio of the anatomy section: names a
script defines for its own use live only until its `unset_all`.

**2. The unset catalog.** Keys that deliberately survive across scripts
until stop (everything from `set_installation_options.sh`,
`flags_derived.sh`, and `start_cocoa.sh`) are unset in one place,
`flags_impl_unset_keys.sh`, which `stop_cocoa.sh` sources. The file is kept
in case-insensitive alphabetical order, variables then functions, so a
missing or duplicated key stands out on inspection; a new key is inserted
at its position, never appended.

[Adapted from installation_scripts/flags_impl_unset_keys.sh]

```bash
# Variables
unset -v ACTDR4_GIT_BRANCH
unset -v ACTDR4_GIT_COMMIT
unset -v ACTDR4_GIT_TAG
unset -v ACTDR4_NAME
unset -v ACTDR4_URL
```

**3. Save and recover.** Variables that already exist in the system and
that Cocoa merely modifies (`PATH`, `PYTHONPATH`, `LD_LIBRARY_PATH`,
`OMP_NUM_THREADS`, compiler flags, and their kin) must not be unset: that
would destroy the user's own values. Instead, `flags_save_old.sh` (sourced
by `start_cocoa.sh`) copies each into an `OLD_` twin before Cocoa touches
it, writing the sentinel value `"x"` (a marker chosen to mean "this
variable did not exist") when there is nothing to save:

[Adapted from installation_scripts/flags_save_old.sh]

```bash
if [ -n "${PYTHONPATH}" ]; then
  export OLD_PYTHONPATH=$PYTHONPATH
else
  export OLD_PYTHONPATH="x"
fi
```

`flags_recover_old.sh` (sourced by `stop_cocoa.sh`) reverses it: restore
the saved value, or unset the variable when the twin holds the sentinel,
and drop the `OLD_` copy either way:

[Adapted from installation_scripts/flags_recover_old.sh]

```bash
if [ -n "${OLD_PYTHONPATH}" ]; then
  if [ "${OLD_PYTHONPATH}" != "x" ]; then
    export PYTHONPATH=$OLD_PYTHONPATH 
  else
    unset PYTHONPATH
  fi
  unset OLD_PYTHONPATH
fi
```

The two files mirror each other block for block: a variable added to one
must be added to the other.

## Adding a new package

New scripts are created by copying an existing one (`setup_pyspk.sh` and
`compile_pyspk.sh` are clean models) and renaming the package-specific
parts. Never invent new structure: the repetition of the skeleton across
all the scripts in this folder is a design decision, made so each script
reads top to bottom on its own, and it must not be refactored into shared
helper files. The checklist, with pointers:

1. Keys in `set_installation_options.sh`: an `IGNORE_XXX_CODE` key plus the
   URL / NAME / version-pin trio, in the matching section.
2. `flags_impl_unset_keys.sh`: one `unset -v` line per new key, each
   inserted at its alphabetical position.
3. `flags_derived.sh`: an `OVERWRITE_EXISTING_XXX_CODE=1` line inside the
   `OVERWRITE_EXISTING_ALL_PACKAGES` cascade, when the script honors an
   OVERWRITE key.
4. Runner lists: the new `setup_XXX.sh` in the matching `setup_cocoa.sh`
   array, and the new `compile_XXX.sh` (if one exists) in
   `compile_cocoa.sh`.
5. `Cocoa/.gitignore`: an entry for the clone destination
   (`external_modules/code/<name>`), so `git add --all` never stages the
   cloned repository.

## Common questions

<a name="faq-compile-offline"></a>
### Why do compile scripts never use the internet?

Compilation often happens on the compute nodes of a cluster, and compute
nodes frequently have no outbound network access. A compile script that
downloaded anything would work on a laptop and die on the cluster, in the
middle of a long installation. The split is therefore strict: every
download (git clone, wget, pip index access, and emulator model files that
a package would otherwise fetch on first use) happens in `setup_*.sh` or
`unxv_*.sh`, and `compile_*.sh` only builds what is already on disk. The
pip flags `--no-dependencies --no-index --no-build-isolation` enforce the
rule mechanically: with them, pip cannot contact an index even by accident.

<a name="faq-cache-skip"></a>
### Why did my second run skip a package?

The first run finished that package's script, the script returned 55, and
the runner wrote a `1` on its line in the cache file
(`.setup_local_packages.txt` for setup,
`.local/compile_local_packages.txt` for compile). On the next run the loop
sees the `1` and moves on. That is the intended behavior: it lets an
interrupted installation resume where it stopped instead of repeating hours
of finished work. A script that was interrupted keeps its `0`, because only
the `return 55` on the script's last line reports completion; a half-made
clone or build is redone on the rerun. To force a redo of a finished
package, source its scripts directly (see the single-package commands
above), or pass a mode option such as `--soft` to the runner to reset a
whole group.

<a name="faq-no-exit"></a>
### Why must these scripts never call exit?

A sourced script runs inside the user's shell. `exit` terminates the
current shell process, so an `exit` in a sourced script closes the user's
terminal session (or kills the batch job) instead of stopping the script.
`return` stops only the sourced file and hands a status code to the
caller. The one `exit` in this codebase proves the rule: the runners and
`start_cocoa.sh` begin with a guard that detects being executed rather than
sourced, and only in that case, in a child process where `exit` is
harmless, do they call it:

[Adapted from start_cocoa.sh]

```bash
if [[ ! "${BASH_SOURCE[0]}" != "$0" ]]; then
  FILE="$(basename ${BASH_SOURCE[0]})"
  MSG="\033[0;31m ${FILE} must be sourced (not executed as program)"
  MSG2=", e.g.: \n source ${FILE}\033[0m"
  echo -e "${MSG}${MSG2}"
  unset FILE MSG MSG2
  exit 1
fi
```

`BASH_SOURCE[0]` is the file being read; `$0` is the name of the running
program. They differ when the file is sourced and match when it is
executed, so this branch runs only in the executed (child-process) case.
The guard exists because an executed installation script would run all its
commands and change nothing in the user's shell, a silent no-op that is
harder to debug than this refusal.

<a name="faq-unset-all-twice"></a>
### What does `unset_all: command not found` mean?

The cleanup ran twice. `error()` already calls `unset_all`, and `unset_all`
ends by deleting itself (`unset -f unset_all`). A call site written as
`cdfolder "${PACKDIR:?}" || { unset_all; return 1; }` therefore fails a
second time inside the failure path: `cdfolder` routed through `error()`,
the cleanup ran and erased the function, and the explicit `unset_all` finds
nothing to call. The fix is the contract from the anatomy section: helpers
that route through `error()` get plain `|| return 1`; only direct commands
and `ptop`/`pbottom` name the cleanup themselves.

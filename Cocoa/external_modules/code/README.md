# Table of contents
1. [Integrating an external code as a Cobaya theory block](#appendix_new_theory_block)
    1. [What the theory-block Python class must do](#appendix_theory_block_class)
    2. [Running the verification](#appendix_theory_block_verification)
2. [FAQ: How to switch Cocoa's adopted CAMB/CLASS? (the easy way)](#appendix_new_camb_class)
3. [FAQ: What about installing Polychord and Velocileptors?](#appendix_new_polychord)
4. [FAQ: How to download and compile likelihoods for modern CMB data?](#new_planck_likelihoods)
5. [FAQ: How to switch the default CAMB/CLASS? (the not-so-easy way)](#appendix_new_camb_class_medium)
6. [FAQ: How to add additional patches to the default CAMB/CLASS?](#appendix_new_camb_class_patches)
7. [Understanding CAMB's patches](#appendix_patch_camb)
8. [Understanding CLASS's patches](#appendix_patch_class)

The folder `Cocoa/external_modules/code` is where external code lives. Every code Cocoa downloads and compiles - Boltzmann codes (`CAMB`, `class_public`), samplers (`PolyChordLite`), CMB likelihood codes, machine-learning emulators, and the Cosmolike ecosystem (`cosmolike_core` plus the project interfaces) - is cloned into its own subfolder here by a `setup_XXX.sh` script and, when needed, compiled by a `compile_XXX.sh` script, both located at `Cocoa/installation_scripts`. New codes belong here as well: never hand-installed into the Conda environment, and never copied into the `cobaya` tree.

Installing a new CAMB code in Cocoa requires a few changes to the existing CAMB/CLASS code. Fortunately, Cocoa provides a set of scripts, located at `Cocoa/installation_scripts`, and patches, located at `Cocoa/../cocoa_installation_libraries/XXX_changes` where XXX is the specific code to be patched, that automatically handle the necessary code adjustments.

The shell scripts are split into two categories. The first category of shell scripts downloads the code and applies necessary patches; they are named `setup_XXX.sh`, where XXX is the specific code to be downloaded. For instance, the script `setup_camb.sh` downloads CAMB Boltzmann code from its original GitHub repository. The second category of shell scripts compiles the code; they are named `compile_XXX.sh`. For instance, the script `setup_camb.sh` compiles CAMB.

The patches are always located at `Cocoa/../cocoa_installation_libraries`. These patches enforce that these codes are compiled with the Cocoa-prescribed compilers and linked against the Cocoa-prescribed version of any necessary numerical library. Consistency when compiling and linking code is one of the main advantages of working within the Cocoa framework.

On an advanced note, three scripts are worth mentioning separately: `setup_core_packages.sh`, `unxv_core_packages.sh`, and `compile_core_packages.sh`. The first two scripts manage the installation, while the third manages the compilation of stable numerical libraries that other codes may need. *The vast majority of these codes are provided by the Cocoa Conda environment*, with a few exceptions, including the CUBA integration library. Nevertheless, these scripts provide a unified and simple interface for users to add new required numerical libraries that may not be available on Conda (or that need some unique way to be compiled).
      
## Integrating an external code as a Cobaya theory block <a name="appendix_new_theory_block"></a>

In Cobaya, a *theory block* is a Python class that computes an intermediate physical quantity at every point the sampler visits. It sits between the Boltzmann code and the likelihood: Cobaya runs the Boltzmann code, hands its products to the theory block, and hands the combined products to the likelihood. Integrating a block into Cocoa means teaching Cocoa's installation scripts to download the code, install its dependencies, and place it where Cobaya looks for theories.

> [!NOTE]
> The worked example in this section is `bfmt`, Cocoa's baryonic feedback block. It multiplies the matter power spectrum $P(k, z)$ by a suppression factor $S(k, z)$ predicted by baryonic feedback emulators, so every likelihood downstream sees a feedback-corrected spectrum. Each step below names the real bfmt files, so users can open them and copy their shape.

A theory block brings two kinds of ingredients, and each ingredient is a separate git repository with its own installation scripts.

| Ingredient | What it is | The bfmt case |
|------------|------------|---------------|
| Theory-block repository | the Python package with the `Theory` class that Cobaya imports | `baryon_suppression` |
| Emulator codes (zero or more) | third-party packages the theory class imports at run time to predict the physical quantity | `pyspk`, `bcemu`, `fbre`, and `baccoemu` |

The integration touches a fixed list of Cocoa scripts, summarized in the table below; each numbered step then walks through one row and explains why that row exists.

| File | What users add there |
|------|----------------------|
| `set_installation_options.sh` | one `IGNORE` switch and one URL/NAME/pin trio per repository |
| `installation_scripts/flags_impl_unset_keys.sh` | one `unset -v` line per new environmental key |
| `installation_scripts/flags_derived.sh` | the new `OVERWRITE` key inside the overwrite cascade |
| `installation_scripts/setup_XXX.sh` (new files) | one download script per repository |
| `setup_cocoa.sh` and `compile_cocoa.sh` | the new scripts appended to the runner lists |
| `installation_scripts/setup_pip_core_packages.sh` | a pinned version of every Python package the new codes import |
| `start_cocoa.sh` and `stop_cocoa.sh` | the symlink that lets Cobaya find the block, created and removed |
| `.gitignore` | one ignore rule per clone destination |

We assume users are in the Conda cocoa environment from a previous `conda activate cocoa` command, that the shell is bash, and that the current folder is the cocoa main folder `cocoa/Cocoa`. Steps :one: to :seven: edit installation scripts; the section [Running the verification](#appendix_theory_block_verification) then exercises the result end to end.

**Step :one::** Declare the new environmental keys in `set_installation_options.sh`.

This shell script is Cocoa's single control panel: every package is described by environmental keys, and no installation script hardcodes a URL or a version. A new repository gets, first, an `IGNORE_XXX_CODE` switch. While the switch stays commented out the package is installed; exporting it skips the package everywhere at once (download, compilation, and the symlink of Step :six:).

    [Adapted from Cocoa/set_installation_options.sh shell script]

    #export IGNORE_PYSPK_CODE=1
    #export IGNORE_BCEMU_CODE=1
    #export IGNORE_FBRE_CODE=1
    #export IGNORE_BACCOEMU_CODE=1

    (...)

    #export IGNORE_BFMT_CODE=1 # Baryon Feedback Theory Block

The repository also gets a URL/NAME/pin trio. The *pin* freezes the exact version that the download script checks out, so an upstream commit can never silently change Cocoa's results. Cocoa pins third-party codes by commit hash, because their owners can move or delete tags, and pins the team-owned theory repository by tag, because the team controls when that tag advances.

    [Adapted from Cocoa/set_installation_options.sh shell script]

    export BFMT_THEORY_URL="https://github.com/CosmoLike/cocoa_baryonic_feedback_models_theory.git"
    export BFMT_NAME="baryon_suppression"    # the folder name inside external_modules/code
    export BFMT_GIT_TAG="v1.0"               # team-owned repository: pinned by TAG

    (...)

    export BCEMU_URL="https://github.com/sambit-giri/BCemu.git"
    export BCEMU_GIT_COMMIT="c32577654e1b48b4bdf079c01a26b16f1473598f"  # third party: pinned by COMMIT
    export BCEMU_NAME="bcemu"

**Step :two::** Register every new key in the flag-bookkeeping scripts.

Cocoa's scripts are sourced, meaning they run inside the user's own shell session rather than in a child process. Every key they define would therefore survive in the session after the script ends, unless it is removed explicitly. The script `flags_impl_unset_keys.sh` is that cleanup list; add one `unset -v` line per new key.

    [Adapted from Cocoa/installation_scripts/flags_impl_unset_keys.sh shell script]

    unset -v BFMT_GIT_BRANCH
    unset -v BFMT_GIT_COMMIT
    unset -v BFMT_GIT_TAG
    unset -v BFMT_NAME
    unset -v BFMT_THEORY_URL

The script `flags_derived.sh` computes derived keys. One of them matters here: when users export `OVERWRITE_EXISTING_ALL_PACKAGES`, setup scripts delete and re-clone every package folder instead of keeping what is on disk, and the cascade must learn that the new package exists.

    [Adapted from Cocoa/installation_scripts/flags_derived.sh shell script]

    if [ -n "${OVERWRITE_EXISTING_ALL_PACKAGES}" ]; then
      (...)
      export OVERWRITE_EXISTING_BCEMU_CODE=1
      (...)
      export OVERWRITE_EXISTING_BFMT_CODE=1
    fi

**Step :three::** Write one `setup_XXX.sh` download script per repository.

Do not write these scripts from scratch: copy the neighboring script whose situation matches and rename its keys (for bfmt, the five scripts `setup_bfmt.sh`, `setup_pyspk.sh`, `setup_bcemu.sh`, `setup_fbre.sh`, and `setup_baccoemu.sh` all follow the standard skeleton). A setup script clones the URL into `external_modules/code/${XXX_NAME}` and checks out the pin.

If the external code needs a source correction, the setup script applies it with a `sed` command placed after the clone block, written so that running it twice is harmless (the pattern matches only the unpatched line). The bfmt integration has a real example: BCemu's package `__init__` imported its `spectra` module, which imports `camb` at load time — shadowing the path-checked CAMB that Cobaya manages — so `setup_bcemu.sh` comments that import out.

    [Adapted from Cocoa/installation_scripts/setup_bcemu.sh shell script]

    # Cocoa patch: BCemu's __init__ imports its spectra module, which imports camb

    (...)

    sed --in-place --regexp-extended \
      's@^from \.spectra import@#from .spectra import@' \

    (...)

> [!NOTE]
> Cocoa enforces an *internet invariant*: setup scripts may use the internet (git clone, wget, the pip index), and every other script must run to completion on a cluster node with no outside connection. This invariant decides where a Python package gets pip-installed. A package that downloads model files on first use is installed in its setup script, which also triggers the download right there (`setup_bcemu.sh`, `setup_baccoemu.sh`). A package that installs fully offline is installed in a compile script instead (`compile_pyspk.sh`, `compile_fbre.sh`).

**Step :four::** Append the new scripts to the runner lists.

The scripts `setup_cocoa.sh` and `compile_cocoa.sh` execute a declared list of scripts in order. Add each new `setup_XXX.sh` to the first and each new `compile_XXX.sh` to the second. A repository whose setup script does everything has no compile script and no compile entry — of the five bfmt repositories, only `pyspk` and `fbre` appear in `compile_cocoa.sh`.

    [Adapted from Cocoa/setup_cocoa.sh shell script]

    declare -a TSCRIPTS=(
                         (...)

                         "setup_bfmt.sh"
                         "setup_pyspk.sh"
                         "setup_bcemu.sh"
                         "setup_fbre.sh"
                         "setup_baccoemu.sh"

                         (...)
                        )

**Step :five::** Seed the runtime Python dependencies.

Every pip install in Cocoa runs with `--no-dependencies`. This is deliberate: pip is never allowed to resolve a dependency tree on its own, because a resolved tree could silently upgrade a core package (NumPy, SciPy) underneath every other code in the environment. The price of this safety is manual work: every package the new codes import at run time must already exist in the environment, which means it must be listed, with a pinned version, in the `PIPCP` list inside `setup_pip_core_packages.sh`.

To find what a package imports, list its import statements:

    grep -rhE "^(import|from) [a-zA-Z0-9_]+" external_modules/code/bcemu/BCemu/*.py | sort -u

Before pinning a version, dry-run it (`pip install --dry-run <package>==<version>`) and confirm the plan does not upgrade anything already pinned. The bfmt integration seeded the packages below.

    [Adapted from Cocoa/installation_scripts/setup_pip_core_packages.sh shell script]

    # pydantic (pyspk); smt + wget + msgpack (BCemu; smt==1.0.0 - other
    # versions are incompatible); swiftemulator (FBRE; pulls george, attrs,
    # unyt, SALib, velociraptor); progressbar2 (baccoemu)
    'pydantic==2.11.4'
    'smt==1.0.0'
    'wget==3.2'
    'msgpack==1.1.0'
    'swiftemulator==1.3.1'
    'progressbar2==4.5.0'

**Step :six::** Create the symlink pair that makes Cobaya see the block.

Cobaya discovers a theory block by folder name: a block named `bfmt` must appear at `cobaya/cobaya/theories/bfmt/`, and that folder must contain a file `bfmt.py` defining `class bfmt(Theory)` at its top level. Cocoa never copies files into the `cobaya` tree. Instead, `start_cocoa.sh` creates a symlink from the clone in `external_modules/code` to that path, and `stop_cocoa.sh` removes it, so a stopped shell leaves `cobaya` exactly as it was cloned.

Add one guarded block to each of the two scripts; the guard must be the block's own `IGNORE` key.

    [Adapted from Cocoa/start_cocoa.sh shell script]

    if [[ -z "${IGNORE_BFMT_CODE}" ]]; then
      ECODEF="${ROOTDIR:?}/external_modules/code"
      COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
      TMP="${BFMT_NAME:-"baryon_suppression"}"
      TMP2="bfmt"   # the folder name Cobaya sees: theory block `bfmt`
      if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
        ln -s "${ECODEF:?}/${TMP}" "${COBTH:?}/${TMP2}" \
          >>${OUT1:?} 2>>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
      fi
      unset -v ECODEF COBTH TMP TMP2
    fi

and

    [Adapted from Cocoa/stop_cocoa.sh shell script]

    if [[ -z "${IGNORE_BFMT_CODE}" ]]; then
      COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
      TMP="bfmt"
      if [[ -L "${COBTH:?}/${TMP}" ]]; then
        rm -f "${COBTH:?}/${TMP:?}"
      fi
      unset -v COBTH TMP
    fi

> [!Warning]
> When copying the guard from a neighboring block, change the `IGNORE` key. A guard left pointing at the neighbor's key makes the new block appear and disappear with the wrong switch — this exact bug has happened.

**Step :seven::** Housekeeping.

Add every clone destination to `.gitignore` and verify each rule with `git check-ignore <folder>`. This step protects the repository: a nested git clone that is not ignored would be staged by `git add --all` as a *gitlink* (an accidental submodule reference), which breaks fresh clones of Cocoa.

    [Adapted from Cocoa/.gitignore]

    external_modules/code/baryon_suppression
    external_modules/code/pyspk
    external_modules/code/bcemu
    external_modules/code/fbre
    external_modules/code/baccoemu

Finally, document the block: quote its keys and usage in the theory repository's own README, and add a short pointer section to the README of every project that consumes it.

### What the theory-block Python class must do <a name="appendix_theory_block_class"></a>

The class itself lives in the theory repository, not in Cocoa, but reviewers of a new block check the four rules below. Each rule prevents a specific failure inside a long sampling run.

| Rule | Why |
|------|-----|
| To reject a sample, log a warning and return `False` from `calculate()`; never raise `LoggedError` | `LoggedError` is on Cobaya's list of always-stop exceptions, so raising it kills the entire run instead of rejecting one point |
| Validate sampled parameters against the emulator's training box; outside it, return a defined fallback (for example $S = 1$, or the clamped boundary value) | an emulator is only trustworthy inside its training box; silent extrapolation poisons chains with plausible-looking numbers |
| Load emulator files once, in `initialize()` | loading per sample adds file input on every one of the many likelihood evaluations of an MCMC run |
| Evaluate on a small internal $(k, z)$ grid and interpolate (2D spline) onto the grid the likelihood requests | the block's cost stays fixed no matter how fine the likelihood's grid is — the same strategy Cobaya uses for its matter-power interpolator |

A block that follows these rules can be listed in any Cocoa project's yaml without changes to the project.

### Running the verification <a name="appendix_theory_block_verification"></a>

We assume users are in the Conda cocoa environment from a previous `conda activate cocoa` command, that the shell is bash, and that the current folder is the cocoa main folder `cocoa/Cocoa`.

**Step :one::** Activate the private Python environment by sourcing `start_cocoa.sh`.

    source start_cocoa.sh

**Step :two::** Check the syntax of every touched shell script before running any of them; `bash -n` parses a script without executing it.

    bash -n installation_scripts/setup_bfmt.sh

**Step :three::** Run each new setup script by itself, so a failure is attributable to one script.

    source installation_scripts/setup_bfmt.sh

**Step :four::** Run each new compile script the same way (skip this step for repositories that have no compile script).

    source installation_scripts/compile_pyspk.sh

**Step :five::** Evaluate a likelihood twice, once with the theory block listed in the yaml and once without it, for example via `cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f` on a project that consumes the block. The run without the block must reproduce the result the project had before the integration, and both runs must complete without errors.

## :interrobang: FAQ: How to switch the Cocoa's adopted CAMB/CLASS? (the easy way) <a name="appendix_new_camb_class"></a> 

Swapping the default CAMB/CLASS is simple. Go to Cocoa's main folder and open the file `set_installation_options.sh`. Then, adjust the following environmental keys. 
     
    [Adapted from Cocoa/set_installation_options.sh shell script]
    
    export CAMB_URL="https://github.com/cmbant/CAMB"
    export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
    export CAMB_NAME='CAMB'
     
    export CLASS_URL="https://github.com/lesgourg/class_public.git"
    export CLASS_GIT_COMMIT="8df566c1ff2d0b3e40e106567c435575aea337be"
    export CLASS_NAME="class_public"

As long as your new CAMB/CLASS makefiles are not altered to the extent that Cocoa patch files designed to adjust [CAMB installation](../../../cocoa_installation_libraries/camb_changes) and [CLASS installation](../../../cocoa_installation_libraries/class_changes)  fail, that is all that is needed to change the adopted Boltzmann codes.

What happens if the Cocoa patch files fail? We have three sections dedicated to this problem. Check [FAQ: How to add additional patches to the default CAMB/CLASS?](#appendix_new_camb_class_patches), and [Understanding CAMB's patches](#appendix_patch_camb), and [Understanding CLASS's patches](#appendix_patch_class).

Finally, what if the user wants to skip CAMB or CLASS compilation? In this case, the shell script `Cocoa/set_installation_options.sh` provides the following environmental keys.
    
    [Adapted from Cocoa/set_installation_options.sh shell script]
    #export IGNORE_CAMB_CODE=1
    #export IGNORE_CLASS_CODE=1

## :interrobang: FAQ: What about installing Polychord and Velocileptors? <a name="appendix_new_polychord"></a> 

The shell script `set_installation_options.sh` provides the following keys that manage the download and compilation of Polychord and Velocileptor.

    [Adapted from Cocoa/set_installation_options.sh shell script]
    
    #export IGNORE_POLYCHORD_SAMPLER_CODE=1
    #export IGNORE_VELOCILEPTORS_CODE=1

    (...)

    export POLY_URL="https://github.com/PolyChord/PolyChordLite.git"
    export POLYCHORD_GIT_COMMIT="daba49d1385d065122db76a2b384050f9e95d278"
    export POLY_NAME="PolyChordLite"

    (...)

    export VELOCILEPTORS_URL="https://github.com/sfschen/velocileptors.git"
    export VELOCILEPTORS_GIT_COMMIT="889a0c98895831eb23b250a26162cfb8a93237bd"
    export VELOCILEPTORS_NAME="velocileptors"
    
## :interrobang: FAQ: How to download and compile likelihoods for modern CMB data? <a name="new_planck_likelihoods"></a>

The CMB data sets require specialized likelihoods. Cocoa will download, patch, and compile them as long as these following keys are not set on `Cocoa/set_installation_options.sh`

    [Adapted from Cocoa/set_installation_options.sh shell script]

    # ------------------------------------------------------------------------------
    # The keys below control which packages will be installed and compiled when ----  
    # running setup_cocoa.sh and compile_cocoa.sh. They are mostly helpful when ----
    # debugging cocoa --------------------------------------------------------------
    # ------------------------------------------------------------------------------
    
    (...)
    
    #export IGNORE_PLANCK_LIKELIHOOD_CODE=1
    #export IGNORE_ACTDR4_CODE=1
    #export IGNORE_ACTDR6_CODE=1
    #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1
    #export IGNORE_CAMSPEC_LIKELIHOOD_CODE=1
    #export IGNORE_LIPOP_LIKELIHOOD_CODE=1

Cocoa selects the URL to download the likelihoods (and the version of the likelihood) using the following keys.

    [Adapted from Cocoa/set_installation_options.sh shell script]
    
    # ------------------------------------------------------------------------------
    # PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
    # ------------------------------------------------------------------------------

    (...)
    
    export HILLIPOP_URL="https://github.com/planck-npipe/hillipop.git"
    export HILLIPOP_GIT_COMMIT="cc9cbe31991d4662522241543a46d44d2cdec251"

    export LOLLIPOP_URL="https://github.com/planck-npipe/lollipop.git"
    export LOLLIPOP_GIT_COMMIT="280a9c93d33bc6a058d6bf769ec82d9f7fdbd2b6"

## :interrobang: FAQ: How to switch the default CAMB/CLASS? (the not-so-easy way) <a name="appendix_new_camb_class_medium"></a> 

If users want to create their setup and compile shell scripts to work seamlessly with multiple modified Boltzmann codes, they must follow the steps below. Here, we assume the users want to create a new modified CAMB named CAMBQ (the modified class case is similar).

**Step :one::** Copy the setup and compile shell scripts. 

    cp "${ROOTDIR:?}"/installation_scripts/setup_camb "${ROOTDIR:?}"/installation_scripts/setup_cambq

and

    cp "${ROOTDIR:?}"/installation_scripts/compile_camb "${ROOTDIR:?}"/installation_scripts/compile_cambq

**Step :two::** Modify the name of the environmental variables `CAMB_URL`, `CAMB_NAME`, and `CAMB_GIT_COMMIT` on `setup_cambq.sh` shell script.

    [Adapted from Cocoa/installation_scripts/setup_cambq.sh shell script]
  
    CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

    (...)
    
    #URL="${CAMB_URL:-"https://github.com/cmbant/CAMB"}"   # Original line - commented
    URL="${CAMBQ_URL:-"https://github.com/CAMBQ"}"   
    
    (...)

    # If the default patches work on the modified CAMBQ, there is no need to modify the line below
    #CHANGES="${CCIL:?}/camb_changes"                     # Original line - commented
    CHANGES="${CCIL:?}/cambq_changes"
    
    #FOLDER="${CAMB_NAME:-"CAMB"}"                        # Original line - commented
    FOLDER="${CAMBQ_NAME:-"CAMBQ"}"
    
    (...)
    
    #if [ -n "${CAMB_GIT_COMMIT}" ]; then                 # Original line - commented
    #  "${GIT:?}" checkout "${CAMB_GIT_COMMIT:?}" \       # Original line - commented
    if [ -n "${CAMBQ_GIT_COMMIT}" ]; then                
      "${GIT:?}" checkout "${CAMBQ_GIT_COMMIT:?}" \       
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

**Step :three::** Modify the name of the environmental variable `CAMB_NAME` on `compile_cambq.sh` shell script.

    [Adapted from Cocoa/installation_scripts/setup_cambq.sh shell script]
    
    #FOLDER="${CAMB_NAME:-"CAMB"}"   # Original line - commented
    FOLDER="${CAMBQ_NAME:-"CAMB"}"  

**Step :four::** Add the environmental variables `CAMB_URL`, `CAMB_NAME`, and `CAMB_GIT_COMMIT` to `set_installation_options.sh`. This is optional if `setup_cambq.sh` and `compile_cambq.sh` provide default reasonable values for these variables. 

    [Adapted from Cocoa/set_installation_options.sh shell script]

    export CAMB_URL="https://github.com/cmbant/CAMB"
    export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
    export CAMB_NAME='CAMB'

    # Add and adapt the lines below
    export CAMBQ_URL="https://github.com/CAMBQ"
    export CAMBQ_GIT_COMMIT="XXX"
    export CAMBQ_NAME='CAMBQ'

**Step :five::** Add `setup_cambq.sh` to the list of files run by `setup_cocoa.sh` as shown below.

    [Adapted from Cocoa/setup_cocoa.sh shell script]

    declare -a TSCRIPTS=("setup_core_packages.sh" 
                         
                         (...)
                         
                         "setup_velocileptors.sh"
                         "setup_cambq.sh) 

**Step :six::** Add `compile_cambq.sh` to the list of files run by `compile_cocoa.sh` as shown below.

    [Adapted from Cocoa/compile_cocoa.sh shell script]

    declare -a TSCRIPTS=("compile_camb.sh"
                         
                         (...)
                         
                         "compile_velocileptors.sh"
                         "compile_cambq.sh")                    
   
**Step :seven::** Add the following line to the script `Cocoa/installation_scripts/flags_impl_unset_keys.sh` so these new environmental variables don't pollute the shell session

    [Adapted from Cocoa/installation_scripts/flags_impl_unset_keys.sh shell script]
    
    #add the line below
    unset -v CAMBQ_URL CAMBQ_GIT_COMMIT CAMBQ_NAME

## :interrobang: FAQ: How to add additional patches to the default CAMB/CLASS? <a name="appendix_new_camb_class_patches"></a> 

Adding additional patches to the default CAMB/CLASS is pretty straightforward. Here, we assume the users want to add a patch to modify CAMB (the modified class case is similar). 

**Step :one::** Copy and save the new patch files to `cocoa_installation_libraries/camb_changes`. 

**Step :two::** Modify the `setup_cambq.sh` shell script as shown below

    [Adapted from Cocoa/installation_scripts/setup_camb.sh shell script]

    # Patch CAMB to be compatible w/ COCOA environment --------------------------

    (...)
    
    # T = TMP
    declare -a TFOLDER=("camb/" 
                        "fortran/" 
                        "forutils/"
                        # Add here the subfolder where the file to be patched is located
                        ) # If nonblank, path must include /
  
    # T = TMP
    declare -a TFILE=("_compilers.py" 
                      "Makefile" 
                      "Makefile_compiler"
                      # Add here the file that needs to be patched
                      )

    #T = TMP, P = PATCH
    declare -a TFILEP=("_compilers.patch" 
                       "Makefile.patch" 
                       "Makefile_compiler.patch"
                       # Add here the file that the patch file
                       )
                     
## Understanding CAMB's patches (developers only :bangbang: :scream: ☠️ :bangbang: ) <a name="appendix_patch_camb"></a> 

To start, we show the current list of CAMB patches below.
    
    cocoa_installation_libraries/camb_changes/camb/_compilers.patch 
    cocoa_installation_libraries/camb_changes/fortran/Makefile.patch
    cocoa_installation_libraries/camb_changes/forutils/Makefile_compiler

Below, we explain what these patches do.

**:one: Patch [camb/_compilers.patch](../../../cocoa_installation_libraries/camb_changes/camb/_compilers.patch)**: This patch modifies the Python function`get_gfortran_version` located in the file `camb/_compilers.py`. 
    
    [Adapted from Cocoa/external_modules/code/CAMB/camb/_compilers.py]
    
    def get_gfortran_version(command='gfortran'):
        #ver = call_command(command + " -dumpversion")           # Original line - commented
        
        # add the line below
        ver = call_command("$FORTRAN_COMPILER -dumpversion")
        
        if ver and '.' not in ver:
            #ver = call_command(command + " -dumpfullversion")   # Original line - commented

            # add the line below
            ver = call_command("$FORTRAN_COMPILER -dumpfullversion")
        return ver
    
    (...)
    
**:two: Patch [fortran/Makefile.patch](../../../cocoa_installation_libraries/camb_changes/fortran/Makefile.patch)**: This patch modifies the file `fortran/Makefile`.

    [Adapted from Cocoa/external_modules/code/CAMB/fortran/Makefile]
    
    #Will detect ifort/gfortran or edit for your compiler            # Original line - commented
    #ifneq ($(COMPILER),gfortran)                                    # Original line - commented
    #  ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)     # Original line - commented
    #else                                                            # Original line - commented
    #  ifortErr = 1                                                  # Original line - commented
    #endif                                                           # Original line - commented
    
    # add the line below
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
   
        (...)
    
    else
        #gfortErr = $(shell which gfortran >/dev/null; echo $$?)    # Original line - commented
        
        # add the line below
        gfortErr = 0
        
        ifeq "$(gfortErr)" "0"
          #Gfortran compiler (version 6+):                          # Original line - commented
          #compiler_ver = $(shell gfortran -dumpversion 2>&1)       # Original line - commented
          #COMPILER = gfortran                                      # Original line - commented
          #F90C     = gfortran                                      # Original line - commented
          
          # add the line below
          COMPILER ?= $(FORTRAN_COMPILER)  
          # add the line below
          F90C     ?= $(FORTRAN_COMPILER)                           
          
          (...)
          #FFLAGS+=-march=native                                    # Original line - commented
          (...)
        endif
     endif

**:three: Patch [forutils/Makefile_compiler.patch](../../../cocoa_installation_libraries/camb_changes/forutils/Makefile_compiler.patch)**: This patch modifies the file `forutils/Makefile_compiler`

    [Adapted from Cocoa/external_modules/code/CAMB/forutils/Makefile_compiler]
    
    #ifneq ($(COMPILER),gfortran)                                    # Original line - commented
    #   ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)    # Original line - commented
    #else                                                            # Original line - commented
    #   ifortErr = 1                                                 # Original line - commented
    #endif                                                           # Original line - commented
    
    # add the line below
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
        
      (...)
    
    else
    
    #  major_version = $(shell gfortran -dumpversion 2>&1 | cut -d " " -f 3 | cut -d. -f 1)  # Original line - commented
    #  ifneq ($(shell test $(major_version) -gt 5; echo $$?),0)                              # Original line - commented
    #    $(error gfortran version 6.3 or higher (or ifort 14+) is required)                  # Original line - commented
    #  endif                                                                                 # Original line - commented
    
    #  compiler_ver = $(shell gfortran -dumpversion 2>&1)                                    # Original line - commented
        
      # add the line below
      F90C ?= $(FORTRAN_COMPILER)
        
      (...)
        
    endif
        
## Understanding CLASS's patches (developers only :bangbang: :scream: ☠️ :bangbang:) <a name="appendix_patch_class"></a> 

To start, we show the current list of CLASS patches below.
    
    cocoa_installation_libraries/class_changes/Makefile.patch
    cocoa_installation_libraries/class_changes/python/setup.patch

Below, we explain what these patches do.

**:one: Patch [Makefile.patch](../../../cocoa_installation_libraries/class_changes/Makefile.patch)**: This patch modifies the file `Makefile` 
    
    [Adapted from Cocoa/external_modules/code/class_public/Makefile]
     
    # your C compiler:
    #CC       = gcc                           # Original line - commented
    #CC       = icc                           # Original line - commented
    #CC       = pgcc                          # Original line - commented   
    
    # add the line below
    CC       ?= $(C_COMPILER) 
   
**:two: Patch [python/setup.patch](../../../cocoa_installation_libraries/class_changes/python/setup.patch)**: This patch modifies the file `python/setup.py` 
    
    [Adapted from Cocoa/external_modules/code/class_public/python/setup.py]
    
    #GCCPATH_STRING = sbp.Popen(                           # Original line - commented
    #    ['gcc -print-libgcc-file-name'],                  # Original line - commented
    #    stdout=sbp.PIPE, shell=True).communicate()[0]     # Original line - commented
    
    # add the line below
    GCCPATH_STRING = sbp.check_output(["$C_COMPILER -print-libgcc-file-name"], shell=True)
    
    (...)
    
    #MVEC_STRING = sbp.Popen(                          # Original line - commented
    #    ['gcc', '-lmvec'],                            # Original line - commented
    #    stderr=sbp.PIPE).communicate()[1]             # Original line - commented
    #if b"mvec" not in MVEC_STRING:                    # Original line - commented
    #    liblist += ["mvec","m"]                       # Original line - commented




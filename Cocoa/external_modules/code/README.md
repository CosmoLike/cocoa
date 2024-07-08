# Table of contents
1. [FAQ: How to switch Cocoa's adopted CAMB/CLASS? (the easy way)](#appendix_new_camb_class)
2. [FAQ: What about Polychord and Velocileptors?](#appendix_new_polychord)
3. [FAQ: How to download and compile likelihoods for modern CMB data?](#new_planck_likelihoods)
4. [FAQ: How to switch the default CAMB/CLASS? (the not-so-easy way)](#appendix_new_camb_class_medium)
5. [FAQ: How to add additional patches to the default CAMB/CLASS?](#appendix_new_camb_class_patches)
6. [Understanding CAMB's patches](#appendix_patch_camb)
7. [Understanding CLASS's patches](#appendix_patch_class)

Installing a new CAMB code in Cocoa requires a few changes to the existing CAMB/CLASS code. Fortunately, Cocoa provides a set of scripts, located at `Cocoa/installation_scripts`, and patches, located at `Cocoa/../cocoa_installation_libraries/XXX_changes` where XXX is the specific code to be patched, that automatically handle the necessary code adjustments.

The shell scripts are split into two categories. The first category of shell scripts downloads the code and applies necessary patches; they are named `setup_XXX.sh`, where XXX is the specific code to be downloaded. For instance, the script `setup_camb.sh` downloads CAMB Boltzmann code from its original GitHub repository. The second category of shell scripts compiles the code; they are named `compile_XXX.sh`. For instance, the script `setup_camb.sh` compiles CAMB.

The patches are always located at `Cocoa/../cocoa_installation_libraries`. These patches enforce that these codes are compiled with the Cocoa-prescribed compilers and linked against the Cocoa-prescribed version of any necessary numerical library. Consistency when compiling and linking code is one of the main advantages of working within the Cocoa framework.

On an advanced note, three scripts are worth mentioning separately: `setup_core_packages.sh`, `unxv_core_packages.sh`, and `compile_core_packages.sh`. The first two scripts manage the installation, while the third manages the compilation of stable numerical libraries that other codes may need. *The vast majority of these codes are provided by the Cocoa Conda environment*, with a few exceptions, including the CUBA integration library. Nevertheless, these scripts provide a unified and simple interface for users to add new required numerical libraries that may not be available on Conda (or that need some unique way to be compiled).
      
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
    #export IGNORE_CAMB_COMPILATION=1
    #export IGNORE_CLASS_COMPILATION=1

## :interrobang: FAQ: What about the installation of Polychord and Velocileptors? <a name="appendix_new_polychord"></a> 

The shell script `set_installation_options.sh` provides the following keys that manage the download and compilation of Polychord and Velocileptor.

    [Adapted from Cocoa/set_installation_options.sh shell script]
    
    #export IGNORE_POLYCHORD_COMPILATION=1
    #export IGNORE_VELOCILEPTORS_COMPILATION=1

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
    
    #export IGNORE_PLANCK_COMPILATION=1
    #export IGNORE_ACTDR4_COMPILATION=1
    #export IGNORE_ACTDR6_COMPILATION=1
    #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_COMPILATION=1
    #export IGNORE_CAMSPEC_LIKELIHOOD_COMPILATION=1
    #export IGNORE_LIPOP_LIKELIHOOD_COMPILATION=1

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




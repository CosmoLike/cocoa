# Table of contents
1. [FAQ: How to switch the default CAMB/CLASS? (the easy way)](#appendix_new_camb_class)
2. [FAQ: How to download and compile likelihoods for modern CMB data?](#new_planck_likelihoods)
3. [FAQ: How to switch the default CAMB/CLASS? (the not-so-easy way)](#appendix_new_camb_class_medium)
4. [FAQ: How to add additional patches to the default CAMB/CLASS?](#appendix_new_camb_class_patches)
5. [Understanding CAMB's patches (developers only)](#appendix_patch_camb)
6. [Understanding CLASS's patches (developers only)](#appendix_patch_class)

Installing a new CAMB code in Cocoa requires a few changes to the existing CAMB/CLASS code. Fortunately, Cocoa provides a set of scripts and patches that automatically handle the necessary adjustments.
      
## :interrobang: FAQ: How to switch the default CAMB/CLASS? (the easy way) <a name="appendix_new_camb_class"></a> 

Swapping the default CAMB/CLASS is simple. Go to Cocoa's main folder and open the file `set_installation_options.sh`. Then, adjust the following environmental keys. 
     
    [Adapted from Cocoa/set_installation_options.sh]
    
    export CAMB_URL="https://github.com/cmbant/CAMB"
    export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
    export CAMB_NAME='CAMB'
     
    export CLASS_URL="https://github.com/lesgourg/class_public.git"
    export CLASS_GIT_COMMIT="8df566c1ff2d0b3e40e106567c435575aea337be"
    export CLASS_NAME="class_public"

As long as your new CAMB/CLASS makefiles are not altered to the extent that [Cocoa patch files designed to adjust CAMB installation](https://github.com/CosmoLike/cocoa/tree/dev/cocoa_installation_libraries/camb_changes) fail, that is all what is needed.

What happens if the Cocoa patch files fail? We have three sections dedicated to this problem. Check [FAQ: How to add additional patches to the default CAMB/CLASS?](#appendix_new_camb_class_patches), and [Understanding CAMB's patches (developers only)](#appendix_patch_camb), and [Understanding CLASS's patches (developers only)](#appendix_patch_class).

### :interrobang: FAQ: How to download and compile likelihoods for modern CMB data? <a name="new_planck_likelihoods"></a>

The CMB data sets require specialized likelihoods. Cocoa will download, patch and compile them as long as these following keys are not set on `Cocoa/set_installation_options.sh`

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

Cocoa selects the URL to download the likelihoods (and the version of the likelohood) via the following keys also shown on `Cocoa/set_installation_options.sh`

    [Adapted from Cocoa/set_installation_options.sh shell script]
    
    # ------------------------------------------------------------------------------
    # PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
    # ------------------------------------------------------------------------------

    (...)
    
    export HILLIPOP_URL="https://github.com/planck-npipe/hillipop.git"
    export HILLIPOP_GIT_COMMIT="cc9cbe31991d4662522241543a46d44d2cdec251"

    export LOLLIPOP_URL="https://github.com/planck-npipe/lollipop.git"
    export LOLLIPOP_GIT_COMMIT="280a9c93d33bc6a058d6bf769ec82d9f7fdbd2b6"

## :interrobang: FAQ: How to switch the default CAMB/CLASS? (the not-so-easy way :heavy_exclamation_mark: :scream: :heavy_exclamation_mark:) <a name="appendix_new_camb_class_medium"></a> 

If users want to create their own `setup_camb.sh` and `compile_camb.sh` scripts so they can work seamlessly with multiple modified Boltznman codes, they need to follow the steps below. Here, we assume the users want to create a new modified CAMB named CAMBQ (the modified class case is similar).

**Step :one::** Copy the setup and compile scripts. 

    cp "${ROOTDIR:?}"/installation_scripts/setup_camb "${ROOTDIR:?}"/installation_scripts/setup_cambq

and

    cp "${ROOTDIR:?}"/installation_scripts/compile_camb "${ROOTDIR:?}"/installation_scripts/compile_cambq

**Step :two::** Modify the name of the environmental variables `CAMB_URL`, `CAMB_NAME`, and `CAMB_GIT_COMMIT` on `setup_cambq.sh` shell script.

    [Adapted from Cocoa/installation_scripts/setup_cambq.sh]
  
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

    [Adapted from Cocoa/installation_scripts/setup_cambq.sh]
    
    #FOLDER="${CAMB_NAME:-"CAMB"}"   # Original line - commented
    FOLDER="${CAMBQ_NAME:-"CAMB"}"  

**Step :four::** Add the environmental variables `CAMB_URL`, `CAMB_NAME`, and `CAMB_GIT_COMMIT` to `set_installation_options.sh`. This is optional if `setup_cambq.sh` and `compile_cambq.sh` provide default reasonable values for these variables. 

    [Adapted from Cocoa/set_installation_options.sh]

    export CAMB_URL="https://github.com/cmbant/CAMB"
    export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
    export CAMB_NAME='CAMB'

    # add and adapt the lines below
    export CAMBQ_URL="https://github.com/CAMBQ"
    export CAMBQ_GIT_COMMIT="XXX"
    export CAMBQ_NAME='CAMBQ'

**Step :five::** Add `setup_cambq.sh` to the list of files run by `setup_cocoa.sh` as shown below.

    [Adapted from Cocoa/setup_cocoa.sh]

    declare -a TSCRIPTS=("setup_core_packages.sh" 
                         (...)
                         "setup_velocileptors.sh"
                         "setup_cambq.sh) 

**Step :six::** Add `compile_cambq.sh` to the list of files run by `compile_cocoa.sh` as shown below.

    [Adapted from Cocoa/compile_cocoa.sh]

    declare -a TSCRIPTS=("compile_camb.sh"
                         (...)
                         "compile_velocileptors.sh"
                         "compile_cambq.sh")                    
   
**Step :seven::** Add the following line to the script `Cocoa/installation_scripts/flags_impl_unset_keys.sh` so these new environmental variables don't pollute the shell session

    [Extracted and adapted from Cocoa/installation_scripts/flags_impl_unset_keys.sh]

    (...)
    
    #add the line below
    unset -v CAMBQ_URL CAMBQ_GIT_COMMIT CAMBQ_NAME

## :interrobang: FAQ: How to add additional patches to the default CAMB/CLASS? <a name="appendix_new_camb_class_patches"></a> 

Adding additional patches to the default CAMB/CLASS is pretty straightforward. Here, we assume the users want to add a patch to modify CAMB (the modified class case is similar). 

**Step :one::** Copy and save the new patch files to `cocoa_installation_libraries/camb_changes`. 

**Step :two::** Modify the `setup_cambq.sh` shell script as shown below

    [Extracted and adapted from Cocoa/installation_scripts/setup_camb.sh]

    # Patch CAMB to be compatible w/ COCOA environment --------------------------

    (...)
    
    # T = TMP
    declare -a TFOLDER=("camb/" 
                        "fortran/" 
                        "forutils/"
                        # add the subfolder where the file to be patched is located
                        ) # If nonblank, path must include /
  
    # T = TMP
    declare -a TFILE=("_compilers.py" 
                      "Makefile" 
                      "Makefile_compiler"
                      # add here the file that needs to be patched
                      )

    #T = TMP, P = PATCH
    declare -a TFILEP=("_compilers.patch" 
                       "Makefile.patch" 
                       "Makefile_compiler.patch"
                       # add here the file that the patch file
                       )
                     
## Understanding CAMB's patches (developers only :bangbang: ☠️ :bangbang: ) <a name="appendix_patch_camb"></a> 

To start, we show below the current list of CAMB patches
    
    cocoa_installation_libraries/camb_changes/camb/_compilers.patch 
    cocoa_installation_libraries/camb_changes/fortran/Makefile.patch
    cocoa_installation_libraries/camb_changes/forutils/Makefile_compiler
    
**:one: Patch `camb/_compilers.patch`**: This patch modifies the Python function`get_gfortran_version` located in the file `camb/_compilers.py`. 
    
    [Extracted and adapted from Cocoa/external_modules/code/CAMB/camb/_compilers.py]
    
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
    
**:two: Patch `fortran/Makefile.patch`**: This patch modifies the file `fortran/Makefile`.

    [Extracted and adapted from Cocoa/external_modules/code/CAMB/fortran/Makefile]
    
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

**:three: Patch `forutils/Makefile.patch`**: This patch modifies the file `forutils/Makefile_compiler`

    [Extracted and adapted from Cocoa/external_modules/code/CAMB/forutils/Makefile_compiler]
    
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
        
## Understanding CLASS's patches (developers only :bangbang: ☠️ :bangbang:) <a name="appendix_patch_class"></a> 

To start, we show below the current list of CLASS patches
    
    cocoa_installation_libraries/class_changes/Makefile.patch
    cocoa_installation_libraries/class_changes/python/setup.patch
    
**:one: Patch `Makefile.patch`**: This patch modifies the file `Makefile` 
    
    [Extracted and adapted from Cocoa/external_modules/code/class_public/Makefile]
     
    # your C compiler:
    #CC       = gcc                           # Original line - commented
    #CC       = icc                           # Original line - commented
    #CC       = pgcc                          # Original line - commented   
    
    # add the line below
    CC       ?= $(C_COMPILER) 
   
**:two: Patch python/setup.patch**: This patch modifies the file `python/setup.py` 
    
    [Extracted and adapted from Cocoa/external_modules/code/class_public/python/setup.py]
    
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




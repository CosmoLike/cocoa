# Table of contents
1. [Switching the default CAMB/CLASS (the easy way)](#appendix_new_camb_class)
2. [Adding a new modified CAMB/CLASS (the not-so-easy way)](#appendix_new_camb_class_hard)
3. [Understanding CAMB's patches (developers only)](#appendix_patch_camb)
4. [Understanding CLASS's patches (developers only)](#appendix_patch_class)

Installing a new CAMB code in Cocoa requires a few changes to the existing CAMB/CLASS code. Fortunately, CoCoA provides the following set of scripts and patches that automatically handle all the necessary changes.

    installation_scripts/setup_camb.sh
    installation_scripts/compile_camb.sh
    cocoa_installation_libraries/camb_changes/camb/_compilers.patch 
    cocoa_installation_libraries/camb_changes/fortran/Makefile.patch
    cocoa_installation_libraries/camb_changes/forutils/Makefile_compiler
      
    installation_scripts/setup_class.sh
    installation_scripts/compile_class.sh
    cocoa_installation_libraries/class_changes/Makefile.patch
    cocoa_installation_libraries/class_changes/python/setup.patch
      
## Switching the default CAMB/CLASS (the easy way) <a name="appendix_new_camb_class"></a> 

Swapping the default CAMB/CLASS is simple. Go to Cocoa's main folder and open the file `set_installation_options.sh`. Then, adjust the following environmental keys. 
     
    [Extracted and adapted from set_installation_options.sh]
    
    export CAMB_URL="https://github.com/cmbant/CAMB"
    export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
    export CAMB_NAME='CAMB'
     
    export CLASS_URL="https://github.com/lesgourg/class_public.git"
    export CLASS_GIT_COMMIT="8df566c1ff2d0b3e40e106567c435575aea337be"
    export CLASS_NAME="class_public"

As long as your CAMB/CLASS makefiles were not altered to the extend that the `.patch` files shown above would fail, that is all what is needed to switch CAMB/CLASS. In this case, the scripts below will handle the download and compilation of CAMB/CLASS Boltznman codes.
     
    # setup scripts: they download CAMB/CLASS and apply the appropriate patches.
    installation_scripts/setup_class.sh
    installation_scripts/setup_camb.sh
     
    installation_scripts/compile_class.sh
    installation_scripts/compile_camb.sh
       
## Understanding CAMB's patches (developers only) <a name="appendix_patch_camb"></a> 
    
**Patch `camb/_compilers.patch`**: This patch modifies the Python function`get_gfortran_version` located in the file `camb/_compilers.py`. 
    
    [Extracted and adapted from ${ROOTDIR}/external_modules/code/CAMB/camb/_compilers.py]
    
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
    
**Patch `fortran/Makefile.patch`**: This patch modifies the file `fortran/Makefile`.

    (...)
    [Extracted and adapted from ${ROOTDIR}/external_modules/code/CAMB/fortran/Makefile]
    
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

**Patch `forutils/Makefile.patch`**: This patch modifies the file `forutils/Makefile_compiler`

    [Extracted and adapted from ${ROOTDIR}/external_modules/code/CAMB/forutils/Makefile_compiler]
    
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
        
## Understanding CLASS's patches (developers only) <a name="appendix_patch_class"></a> 

**Patch `Makefile.patch`**: This patch modifies the file `Makefile` 
    
    (...)
    [Extracted and adapted from ${ROOTDIR}/external_modules/code/class_public/Makefile]
     
    # your C compiler:
    #CC       = gcc                           # Original line - commented
    #CC       = icc                           # Original line - commented
    #CC       = pgcc                          # Original line - commented   
    
    # add the line below
    CC       ?= $(C_COMPILER) 
   
**Patch python/setup.patch**: This patch modifies the file `python/setup.py` 
    
    [Extracted and adapted from ${ROOTDIR}/external_modules/code/class_public/python/setup.py]
    
    #GCCPATH_STRING = sbp.Popen(                                     # Original line - commented
    #    ['gcc -print-libgcc-file-name'],                            # Original line - commented
    #    stdout=sbp.PIPE, shell=True).communicate()[0]               # Original line - commented
    
    # add the line below
    GCCPATH_STRING = sbp.check_output(["$C_COMPILER -print-libgcc-file-name"], shell=True)
    
    (...)
    
    #MVEC_STRING = sbp.Popen(                          # Original line - commented
    #    ['gcc', '-lmvec'],                            # Original line - commented
    #    stderr=sbp.PIPE).communicate()[1]             # Original line - commented
    #if b"mvec" not in MVEC_STRING:                    # Original line - commented
    #    liblist += ["mvec","m"]                       # Original line - commented
    

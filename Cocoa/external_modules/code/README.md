# Table of contents
1. [Adding a new modified CAMB](#appendix_new_camb)
2. [Adding a new modified CLASS](#appendix_new_class)
 
## Adding a new modified CAMB <a name="appendix_new_camb"></a> 

Installing a new CAMB code in Cocoa requires a few changes to the existing code that ensures the CAMB scripts use the correct compiler and that Cocoa's shell scripts can compile and link CAMB. In this section, we summarize these steps, which can be helpful for users who maintain modified versions of the Boltzmann code. We assume the user opted for the easier Cocoa Conda installation and located the terminal in the folder where Cocoa was cloned.

**Step 1 of 13**: Activate the Cocoa conda environment, then go to the Cocoa main folder, and start the `(.local)` Python private environment.

    $ conda activate cocoa
    $(cocoa) cd ./cocoa/Cocoa
    $(cocoa) source start_cocoa
    
**Step 2 of 13**: Move the modified CAMB code to the folder `$ROOTDIR/external_modules/code/XXX`. The string `XXX` should be replaced by the new code's adopted name (e.g., CAMBQ). 
    
**Step 3 of 13**: Modify the `get_gfortran_version` Python function located in the file `$ROOTDIR/external_modules/code/XXX/camb/_compilers.py`. `$ROOTDIR` is the absolute path of the main `Cocoa/` folder.
    
    (...)
    
    def get_gfortran_version(command='gfortran'):
        #ver = call_command(command + " -dumpversion")
        ver = call_command("$FORTRAN_COMPILER -dumpversion")
        if ver and '.' not in ver:
            #ver = call_command(command + " -dumpfullversion")
            ver = call_command("$FORTRAN_COMPILER -dumpfullversion")
        return ver
    
    (...)
    
**Step 4 of 13**: Modify the following lines located in the file `$ROOTDIR/external_modules/code/XXX/fortran/Makefile`

    (...)
    
    #Will detect ifort/gfortran or edit for your compiler
    #ifneq ($(COMPILER),gfortran)
    #  ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
    #else
    #  ifortErr = 1
    #endif
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
        (...)
    else
        #gfortErr = $(shell which gfortran >/dev/null; echo $$?)
        gfortErr = 0
        ifeq "$(gfortErr)" "0"
          #Gfortran compiler (version 6+):
          #compiler_ver = $(shell gfortran -dumpversion 2>&1)
          #COMPILER = gfortran
          #F90C     = gfortran
          COMPILER ?= $(FORTRAN_COMPILER)
          F90C     ?= $(FORTRAN_COMPILER)
          (...)
          #FFLAGS+=-march=native
          (...)
        endif
     endif

**Step 5 of 13**: Modify the following lines located in the file `$ROOTDIR/external_modules/code/XXX/forutils/Makefile_compiler`

    (...)
    #ifneq ($(COMPILER),gfortran)
    #  ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
    #else
    #  ifortErr = 1
    #endif
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
        (...)
    else
        #major_version = $(shell gfortran -dumpversion 2>&1 | cut -d " " -f 3 | cut -d. -f 1)
        #ifneq ($(shell test $(major_version) -gt 5; echo $$?),0)
        #  $(error gfortran version 6.3 or higher (or ifort 14+) is required)
        #endif
        #compiler_ver = $(shell gfortran -dumpversion 2>&1)
        F90C ?= $(FORTRAN_COMPILER)
        (...)
     endif

**Step 6 of 13**: Copy the file `$ROOTDIR/installation_scripts/compile_camb` to `$ROOTDIR/installation_scripts/` and rename it to compile_XXX via the command.

    $(cocoa)(.local) cp $ROOTDIR/installation_scripts/compile_camb $ROOTDIR/installation_scripts/compile_XXX

**Step 7 of 13**: Modify the following lines located in the file `$ROOTDIR/installation_scripts/compile_XXX`

    (...)
    
    if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
        #echo 'COMPILING CAMB'
        #cd $ROOTDIR/external_modules/code/CAMB/
        echo 'COMPILING XXX'
        cd $ROOTDIR/external_modules/code/XXX/
        (...)
    fi

**Step 8 of 13**: Modify the following lines located in the file `compile_external_modules` located at $ROOTDIR

    (...)
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # ------------------------ COMPILE EXTERNAL MODULES --------------------------
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    source $ROOTDIR/installation_scripts/compile_camb
    # add the line below
    source $ROOTDIR/installation_scripts/compile_XXX
    (...)

This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CAMB

**Step 9 of 13**: Go to Cocoa main directory and compile the new modified CAMB

    $(cocoa)(.local) cd $ROOTDIR
    $(cocoa)(.local) source ./installation_scripts/compile_XXX

PS: you can also use `source compile_external_modules` instead of `source ./installation_scripts/compile_XXX`

**Step 10 of 13**: Copy the file `$ROOTDIR/installation_scripts/clean_camb` to `$ROOTDIR/installation_scripts/` and rename it to clean_XXX via the command.

    $(cocoa)(.local) cp $ROOTDIR/installation_scripts/clean_camb $ROOTDIR/installation_scripts/clean_XXX

**Step 11 of 13**: Modify he following lines located in the file `$ROOTDIR/installation_scripts/clean_XXX`

    (...) 
    #cd $ROOTDIR/external_modules/code/camb/
    cd $ROOTDIR/external_modules/code/XXX/
    (...)

**Step 12 of 13**: Modify the script `clean_all` located at $ROOTDIR (Cocoa main folder)

    (...)
    source $ROOTDIR/installation_scripts/clean_camb
    # add the line below
    source $ROOTDIR/installation_scripts/clean_XXX
    
    (...)

**Step 13 of 13**: Modify any YAML file that should load the new CAMB, adding the option `path` to the CAMB section

    (...)
    theory:
        camb:
            path: ./external_modules/code/XXX   
            (...)
    (...)
    
## Adding a new modified CLASS <a name="appendix_new_class"></a> 

**Step 1 of 12**: Activate Conda environment, go to the Cocoa main folder and start the private python environment

    $ conda activate cocoa
    $(cocoa) cd ./cocoa/Cocoa
    $(cocoa) source start_cocoa
    
**Step 2 of 12**: Move the Boltzmann code to `./external_modules/code/XXX`

`XXX` should be replaced by whatever name the user adopts to their modified CLASS (e.g., CLASSQ). 

**Step 3 of 12**: Modify the file `./external_modules/code/XXX/Makefile` 
    
    (...)
    
    # your C compiler:
    CC       ?= $(C_COMPILER) 
    #CC       = icc
    #CC       = pgcc
    
    (...)
    
    classy: libclass.a python/classy.pyx python/cclassy.pxd
    #ifdef OMPFLAG
    #cp python/setup.py python/autosetup.py
    #else
    #grep -v "lgomp" python/setup.py > python/autosetup.py
    #endif
    #cd python; export CC=$(CC); $(PYTHON) autosetup.py install || $(PYTHON) autosetup.py install --user
    #rm python/autosetup.py
    
    (...)

**Step 4 of 12**: Modify the file `./external_modules/code/XXX/python/setup.py` 
    
    (...)
    
    #GCCPATH_STRING = sbp.Popen(
    #    ['gcc -print-libgcc-file-name'],
    #    stdout=sbp.PIPE, shell=True).communicate()[0]
    GCCPATH_STRING = sbp.check_output(["$C_COMPILER -print-libgcc-file-name"], shell=True)
    GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()
    
    (...)
    
    #MVEC_STRING = sbp.Popen(
    #    ['gcc', '-lmvec'],
    #    stderr=sbp.PIPE).communicate()[1]
    #if b"mvec" not in MVEC_STRING:
    #    liblist += ["mvec","m"]
    
    (...)
    
**Step 5 of 12**: Copy `./installation_scripts/compile_class` to `./installation_scripts/compile_XXX` via the command.

    $(cocoa)(.local) cp ./installation_scripts/compile_class ./installation_scripts/compile_XXX

**Step 6 of 12**: Modify `./installation_scripts/compile_XXX`

    (...)
    
    if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
        echo 'COMPILING XXX'

        cd $ROOTDIR/external_modules/code/XXX/
        
        (...)
        
**Step 7 of 12**: Modify [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/compile_external_modules)

    (...)
    
    source $ROOTDIR/installation_scripts/compile_class
    
    source $ROOTDIR/installation_scripts/compile_XXX
    
    (...)

(**expert**) This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CLASS

**Step 8 of 12**: Compile CLASS

    $(cocoa)(.local) source ./installation_scripts/compile_XXX

**Step 9 of 12**: Copy `./installation_scripts/clean_class` to `./installation_scripts/clean_XXX` via the command.

    $(cocoa)(.local) cp ./installation_scripts/clean_class ./installation_scripts/clean_XXX

**Step 10 of 12**: Modify `./installation_scripts/clean_XXX`

    (...) 
    
    cd $ROOTDIR/external_modules/code/XXX/
    
    (...)

**Step 11 of 12**: Modify [clean_all](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/clean_all)

    (...)

    source $ROOTDIR/installation_scripts/clean_camb
    
    source $ROOTDIR/installation_scripts/clean_XXX
    
    (...)
    
**Step 12 of 12**: Modify any YAML file that loads the new CLASS, adding the option `path` to the CLASS section

    (...)
    
    theory:
        classy:
            path: ./external_modules/code/XXX   
            (...)
    
    (...)


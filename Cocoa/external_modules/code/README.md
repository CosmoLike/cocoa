# Table of contents
1. [Adding a new modified CAMB](#appendix_new_camb)
2. [Adding a new modified CLASS](#appendix_new_class)
 
## Adding a new modified CAMB <a name="appendix_new_camb"></a> 

Installing a new CAMB code requires a few additional steps to ensure that CAMB scripts use the correct compiler, and Cocoa's shell scripts can compile and link CAMB. This section is helpful for users who possess a modified version of the Boltzmann code. Again, we assume the user opted for the easier *Conda installation* and located the terminal at the folder *where Cocoa was cloned*

**Step 1 of 10**: Activate Conda environment, go to the Cocoa main folder and start the private python environment

    $ conda activate cocoa
    $(cocoa) cd ./cocoa/Cocoa
    $(cocoa) source start_cocoa
    
**Step 2 of 10**: Move the Boltzmann code to `./external_modules/code/XXX`

`XXX` should be replaced by whatever name the user adopts to their modified CAMB (e.g., CAMBQ). 
    
**Step 3 of 10**: Modify the file `./external_modules/code/XXX/camb/_compilers.py` 
    
    (...)
    
    def get_gfortran_version():
        ver = call_command("$FORTRAN_COMPILER -dumpversion")
        if ver and '.' not in ver:
            ver = call_command("$FORTRAN_COMPILER -dumpfullversion")
        return ver
    
    (...)
    
**Step 4 of 10**: Modify the file `./external_modules/code/XXX/fortran/Makefile`

    (...)
    
    #Will detect ifort/gfortran or edit for your compiler
    #ifneq ($(COMPILER),gfortran)
    #ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
    #else
    #ifortErr = 1
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
        COMPILER ?= $(FORTRAN_COMPILER)
        F90C     ?= $(FORTRAN_COMPILER)
        #COMPILER = gfortran
        #F90C     = gfortran
        
        (...)
        
        #FFLAGS+=-march=native
        
        (...)
        
        endif
     endif

**Step 5 of 10**: Modify the file `./external_modules/code/XXX/forutils/Makefile_compiler`

    (...)
    
    #For standalone compiling set the compiler
    #ifneq ($(COMPILER),gfortran)
    #ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
    #else
    #ifortErr = 1
    #endif
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
        (...)
    else
        #major_version = $(shell gfortran -dumpversion 2>&1 | cut -d " " -f 3 | cut -d. -f 1)
        #ifneq ($(shell test $(major_version) -gt 5; echo $$?),0)
        #$(error gfortran version 6.3 or higher (or ifort 14+) is required)
        #endif
        #compiler_ver = $(shell gfortran -dumpversion 2>&1)
        F90C ?= $(FORTRAN_COMPILER)
        
        (...)
     endif

**Step 6 of 10**: Copy `./installation_scripts/compile_camb` to `./installation_scripts/compile_XXX` via the command.

    $(cocoa)(.local) cp ./installation_scripts/compile_camb ./installation_scripts/compile_XXX

**Step 7 of 10**: Modify `./installation_scripts/compile_XXX`

    (...)
    
    if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
        echo 'COMPILING XXX'

        cd $ROOTDIR/external_modules/code/XXX/
        
        (...)

**Step 8 of 10**: Add `source $ROOTDIR/installation_scripts/compile_XXX` command to the shell script [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/compile_external_modules)

    (...)
    
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # ------------------------ COMPILE EXTERNAL MODULES --------------------------
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    source $ROOTDIR/installation_scripts/compile_camb
    
    source $ROOTDIR/installation_scripts/compile_XXX
    
    (...)

(**expert**) This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CAMB

**Step 9 of 10**: Compile CAMB

    $(cocoa)(.local) source ./installation_scripts/compile_XXX
 
**Step 10 of 10**: Modify any YAML file that loads the new CAMB, adding the option `path` to the CAMB section

    (...)
    
    theory:
        camb:
            path: ./external_modules/code/XXX   
            (...)
    
    (...)

## Adding a new modified CLASS <a name="appendix_new_class"></a> 

**Step 1 of 9**: Activate Conda environment, go to the Cocoa main folder and start the private python environment

    $ conda activate cocoa
    $(cocoa) cd ./cocoa/Cocoa
    $(cocoa) source start_cocoa
    
**Step 2 of 9**: Move the Boltzmann code to `./external_modules/code/XXX`

`XXX` should be replaced by whatever name the user adopts to their modified CLASS (e.g., CLASSQ). 

**Step 3 of 9**: Modify the file `./external_modules/code/XXX/Makefile` 
    
    (...)
    
    # your C compiler:
    CC       ?= $(C_COMPILER) 
    #CC       = icc
    #CC       = pgcc
    
    (...)

**Step 4 of 9**: Modify the file `./external_modules/code/XXX/python/setup.py` 
    
    (...)
    
    GCCPATH_STRING = sbp.Popen(
        ['$C_COMPILER -print-libgcc-file-name'],
        stdout=sbp.PIPE, shell=True).communicate()[0]
    GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()
    
    (...)
    
**Step 5 of 9**: Copy `./installation_scripts/compile_class` to `./installation_scripts/compile_XXX` via the command.

    $(cocoa)(.local) cp ./installation_scripts/compile_class ./installation_scripts/compile_XXX

**Step 6 of 9**: Modify `./installation_scripts/compile_XXX`

    (...)
    
    if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
        echo 'COMPILING XXX'

        cd $ROOTDIR/external_modules/code/XXX/
        
        (...)
        
**Step 7 of 9**: Add `source $ROOTDIR/installation_scripts/compile_XXX` command to the shell script [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/compile_external_modules)

    (...)
    
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # ------------------------ COMPILE EXTERNAL MODULES --------------------------
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    source $ROOTDIR/installation_scripts/compile_class
    
    source $ROOTDIR/installation_scripts/compile_XXX
    
    (...)

(**expert**) This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CLASS

**Step 8 of 9**: Compile CLASS

    $(cocoa)(.local) source ./installation_scripts/compile_XXX
    
**Step 9 of 9**: Modify any YAML file that loads the new CLASS, adding the option `path` to the CLASS section

    (...)
    
    theory:
        classy:
            path: ./external_modules/code/XXX   
            (...)
    
    (...)


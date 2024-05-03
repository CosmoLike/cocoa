# Table of contents
1. [Adding a new modified CAMB](#appendix_new_camb)
2. [Adding a new modified CLASS](#appendix_new_class)
 
## Adding a new modified CAMB <a name="appendix_new_camb"></a> 

Installing a new CAMB code in Cocoa requires a few changes to the existing code that ensure the CAMB scripts use the correct compiler and that Cocoa's shell scripts can compile and link CAMB. In this section, we summarize these steps, which can be helpful for users who maintain modified versions of the Boltzmann code. We assume the user performed the Cocoa Conda installation and located the terminal in the folder where Cocoa was cloned.

**Step :one:**: Activate the Cocoa conda environment, go to the Cocoa main folder, and start the `(.local)` Python private environment.

    conda activate cocoa
    cd ./cocoa/Cocoa
    source start_cocoa
    
**Step :two:**: Move (or clone) the modified CAMB code to the folder `$ROOTDIR/external_modules/code/XXX`. The string `XXX` should be replaced by the new code's adopted name (e.g., CAMBQ). 

If the user does not want this modified CAMB to be added to the Cocoa repo, then add its path to the `.gitignore` file located in the `/cocoa/Cocoa` main folder.
    
**Step :three:**: Modify the `get_gfortran_version` Python function located in the file `$ROOTDIR/external_modules/code/XXX/camb/_compilers.py`. Here, `$ROOTDIR` is the absolute path of the main `cocoa/Cocoa/` folder.
    
    (...)
    [Extracted and adapted from $ROOTDIR/external_modules/code/XXX/camb/_compilers.py]
    
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
    
**Step 4️⃣**: Modify the following lines located in the file `$ROOTDIR/external_modules/code/XXX/fortran/Makefile`.

    (...)
    [Extracted and adapted from $ROOTDIR/external_modules/code/CAMB/fortran/Makefile]
    
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

These modifications will force CAMB to adopt the conda gfortran compiler. 

**Step :five:**: Modify the following lines located in the file `$ROOTDIR/external_modules/code/XXX/forutils/Makefile_compiler`

    (...)
    [Extracted and adapted from $ROOTDIR/external_modules/code/CAMB/forutils/Makefile_compiler]
    
    #ifneq ($(COMPILER),gfortran)                                     # Original line - commented
    #  ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)      # Original line - commented
    #else                                                             # Original line - commented
    #  ifortErr = 1                                                   # Original line - commented
    #endif                                                            # Original line - commented
    
    # add the line below
    ifortErr = 1
    
    (...)
    
    ifeq "$(ifortErr)" "0"
        (...)
    else
        #major_version = $(shell gfortran -dumpversion 2>&1 | cut -d " " -f 3 | cut -d. -f 1)  # Original line - commented
        #ifneq ($(shell test $(major_version) -gt 5; echo $$?),0)                              # Original line - commented
        #  $(error gfortran version 6.3 or higher (or ifort 14+) is required)                  # Original line - commented
        #endif                                                                                 # Original line - commented
        #compiler_ver = $(shell gfortran -dumpversion 2>&1)                                    # Original line - commented
        
        # add the line below
        F90C ?= $(FORTRAN_COMPILER)
        
        (...)
     endif

These modifications will force the `forutils` package inside CAMB to adopt the conda gfortran compiler. 

**Step :six:**: Copy the file `$ROOTDIR/installation_scripts/compile_camb` to `$ROOTDIR/installation_scripts/` and rename it to compile_XXX via the command.

    cp $ROOTDIR/installation_scripts/compile_camb $ROOTDIR/installation_scripts/compile_XXX

This script will manage CAMB compilation

**Step :seven:**: Modify the following lines located in the file `$ROOTDIR/installation_scripts/compile_XXX`

    (...)
    [Extracted and adapted from $ROOTDIR/installation_scripts/compile_camb]
    
    if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
        #cd $ROOTDIR/external_modules/code/CAMB/                      # Original line - commented
        
        # add the line below
        cd $ROOTDIR/external_modules/code/XXX/
        
        (...)
    fi

    (...)

**Step :eight:**: Modify the following lines located in the file `compile_external_modules` located at $ROOTDIR

    (...)
    [Extracted and adapted from $ROOTDIR/compile_external_modules]
    
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # ------------------------ COMPILE EXTERNAL MODULES --------------------------
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    
    (...)
    
    source $ROOTDIR/installation_scripts/compile_camb
    
    # add the line below
    source $ROOTDIR/installation_scripts/compile_XXX
    
    (...)

This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CAMB

**Step :nine:**: Go to Cocoa main directory and compile the new modified CAMB

    cd $ROOTDIR
    source ./installation_scripts/compile_XXX

PS: you can also use `source compile_external_modules` instead of `source ./installation_scripts/compile_XXX`

**Step :ten:**: Copy the file `$ROOTDIR/installation_scripts/clean_camb` to `$ROOTDIR/installation_scripts/` and rename it to clean_XXX via the command.

    cp $ROOTDIR/installation_scripts/clean_camb $ROOTDIR/installation_scripts/clean_XXX

**Step :eleven:**: Modify the following lines located in the file `$ROOTDIR/installation_scripts/clean_XXX`

    (...) 
    [Extracted and adapted from $ROOTDIR/installation_scripts/clean_camb]
    
    #cd $ROOTDIR/external_modules/code/camb/                            # Original line - commented
    
    # add the line below
    cd $ROOTDIR/external_modules/code/XXX/
    
    (...)

**Step :twelve:**: Modify the script `clean_all` located at $ROOTDIR (Cocoa main folder)

    (...)
    
    source $ROOTDIR/installation_scripts/clean_camb
    
    # add the line below
    source $ROOTDIR/installation_scripts/clean_XXX
    
    (...)

**Step :thirteen:**: Modify any YAML file that should load the new CAMB, adding the option `path` to the CAMB section

    (...)
    
    theory:
        camb:
            path: ./external_modules/code/XXX   
            (...)
    (...)
    
## Adding a new modified CLASS <a name="appendix_new_class"></a> 

**Step :one:**: Activate Conda environment, go to the Cocoa main folder and start the private python environment

    conda activate cocoa
    cd ./cocoa/Cocoa
    source start_cocoa
    
**Step :two:**: Move the Boltzmann code to `$ROOTDIR/external_modules/code/XXX`

`XXX` should be replaced by whatever name the user adopts to their modified CLASS (e.g., CLASSQ). 

**Step :three:**: Modify the file `$ROOTDIR/external_modules/code/XXX/Makefile` 
    
    (...)
    [Extracted and adapted from $ROOTDIR/external_modules/code/class_public/]
     
    # your C compiler:
    #CC       = gcc                           # Original line - commented
    #CC       = icc                           # Original line - commented
    #CC       = pgcc                          # Original line - commented   
    
    # add the line below
    CC       ?= $(C_COMPILER) 
   
**Step :four:**: Modify the file `$ROOTDIR/external_modules/code/XXX/python/setup.py` 
    
    (...)
    [Extracted and adapted from $ROOTDIR/external_modules/code/class_public//python/setup.py]
    
    #GCCPATH_STRING = sbp.Popen(                                     # Original line - commented
    #    ['gcc -print-libgcc-file-name'],                            # Original line - commented
    #    stdout=sbp.PIPE, shell=True).communicate()[0]               # Original line - commented
    
    # add the line below
    GCCPATH_STRING = sbp.check_output(["$C_COMPILER -print-libgcc-file-name"], shell=True)
    
    (...)
    
    #MVEC_STRING = sbp.Popen(                                        # Original line - commented
    #    ['gcc', '-lmvec'],                                          # Original line - commented
    #    stderr=sbp.PIPE).communicate()[1]                           # Original line - commented
    #if b"mvec" not in MVEC_STRING:                                  # Original line - commented
    #    liblist += ["mvec","m"]                                     # Original line - commented
    
    (...)

**Step :five:**: Move the folder `$ROOTDIR/external_modules/code/XXX/include` to `$ROOTDIR/external_modules/code/XXX/include2`

Cocoa git repository has a restriction on `.gitignore` against adding `include/` folders. So move `/include/` to `/include2`. 

     mv `$ROOTDIR/external_modules/code/XXX/include` `$ROOTDIR/external_modules/code/XXX/include2`

This will not break the code. The script `$ROOTDIR/installation_scripts/compile_class` that compiles CLASS knows that the headers should be located on `/include2` as shown below

    (...)
    [Extracted and adapted from $ROOTDIR/installation_scripts/compile_class]
    #Workaround around Cocoa .gitignore entry on /include
    rm -rf ./include
    cp -r ./include2 ./include
    (...)
     
**Step :six:**: Copy `$ROOTDIR/installation_scripts/compile_class` to `$ROOTDIR/installation_scripts/compile_XXX` via the command.

    cp $ROOTDIR/installation_scripts/compile_class $ROOTDIR/installation_scripts/compile_XXX

**Step :seven:**: Modify `$ROOTDIR/installation_scripts/compile_XXX`

    (...)
    [Extracted from $ROOTDIR/installation_scripts/compile_class]
    
    if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
        echo 'COMPILING XXX'

        cd $ROOTDIR/external_modules/code/XXX/
        
        (...)
    
    (...)
    
**Step :eight:**: Modify [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/compile_external_modules)

    (...)
    [Extracted and adapted from from $ROOTDIR/compile_external_modules]
     
    source $ROOTDIR/installation_scripts/compile_class
     
    # add the line below
    source $ROOTDIR/installation_scripts/compile_XXX
    
    (...)

This step ensures that the command `source compile_external_modules` compiles all modules including the new modified CLASS

**Step :nine:**: Compile CLASS

    source ./installation_scripts/compile_XXX

**Step :ten:**: Copy `$ROOTDIR/installation_scripts/clean_class` to `$ROOTDIR/installation_scripts/clean_XXX` via the command.

    cp $ROOTDIR/installation_scripts/clean_class $ROOTDIR/installation_scripts/clean_XXX

**Step :eleven:**: Modify `$ROOTDIR/installation_scripts/clean_XXX`

    (...) 
    [Extracted and adapted from from $ROOTDIR/installation_scripts/clean_class]

    # ---------------------------------------------------------------------------
    # ---------------------------------------------------------------------------
    #cd $ROOTDIR/external_modules/code/class_public/                # Original line - commented
    
    # add the line below
    cd $ROOTDIR/external_modules/code/XXX/
    # ---------------------------------------------------------------------------
    # ---------------------------------------------------------------------------
   
    (...)

**Step :twelve:**: Modify [clean_all](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/clean_all)

    (...)
    [Extracted and adapted from from $ROOTDIR/clean_all]
    
    source $ROOTDIR/installation_scripts/clean_camb
    
    # add the line below
    source $ROOTDIR/installation_scripts/clean_XXX
    
    (...)
    
**Step :thirteen:**: Modify any YAML file that loads the new CLASS, adding the option `path` to the CLASS section

    (...)
    
    theory:
        classy:
            path: ./external_modules/code/XXX   
            (...)
    
    (...)


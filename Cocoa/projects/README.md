# Table of contents
1. [The Projects Folder](#appendix_projects_folder)
2. [FAQ: How do we download and run Cosmolike projects?](#running_cosmolike_projects)
3. [FAQ: How do we set Weak Lensing YAML files in Cobaya?](#appendix_example_runs)
4. [FAQ: How do we set Slow/Fast decomposition with Cosmolike?](#manual_blocking_cosmolike)
5. [FAQ: How do we create a new Cosmolike project?](#appendix_lsst_y1_new)
   1. [The easy way](#appendix_lsst_y1_new_small)
   2. [The hard way](#appendix_lsst_y1_new_small2)
 
## The Projects Folder <a name="appendix_projects_folder"></a> 

The `projects` folder includes all the projects linked to Cosmolike; they can also help organize general investigations, even if they don't use Cosmolike directly. Projects that utilize Cosmolike need to have more or less the following structure, taken from the [LSST_Y1 project](https://github.com/CosmoLike/cocoa_lsst_y1)

    +-- cocoa_lsst_y1
    |    +-- likelihood
    |    |   +-- _cosmolike_prototype_base.py
    |    |   +-- lsst_3x2pt.py
    |    |   +-- lsst_3x2pt.yaml
    |    |   +-- lsst_2x2pt.py
    |    |   +-- lsst_2x2pt.yaml
    |    |   +-- lsst_clustering.py
    |    |   +-- lsst_clustering.yaml
    |    |   +-- lsst_cosmic_shear.py
    |    |   +-- lsst_cosmic_shear.yaml
    |    +-- scripts
    |    |   +-- compile_lsst_y1
    |    |   +-- start_lsst_y1
    |    |   +-- stop_lsst_y1
    |    +-- data
    |    |   +-- LSST_Y1.dataset
    |    |   +-- datavector.txt
    |    |   +-- covariance.txt
    |    |   +-- nzlens.txt
    |    |   +-- nzsource.txt
    |    |   +-- mask.mask
    |    +-- interface
    |    |   +-- MakefileCosmolike
    |    |   +-- cosmolike_lsst_y1_interface.py
    |    |   +-- interface.cpp
    |    |   +-- interface.hpp
    |    +-- chains
    |    |   +-- README
    |    +-- EXAMPLE_EVALUATE_1.YAML
    |    +-- EXAMPLE_MCMC_1.YAML

> [!Note]
> Projects should be hosted on independent GitHub repositories. By convention, the Cosmolike Organization adds the prefix `cocoa_` to all Cobaya-Cosmolike projects. For instance, the repository `cocoa_XXX` targets project `XXX`. 

## :interrobang: FAQ: How do we download and run Cosmolike projects? <a name="running_cosmolike_projects"></a> 

### Part I: The semi-automatic way 

Cocoa's `set_installation_options.sh` shell script includes instructions to install several Cosmolike projects. To activate them, manipulate the following lines on `set_installation_options.sh` 

     [Adapted from Cocoa/set_installation_options.sh shell script]

     # ------------------------------------------------------------------------------
     # The keys below control which cosmolike projects will be installed and compiled
     # ------------------------------------------------------------------------------
     export IGNORE_COSMOLIKE_LSSTY1_CODE=1
     export IGNORE_COSMOLIKE_DES_Y3_CODE=1
     #export IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE=1
     #export IGNORE_COSMOLIKE_ROMAN_REAL_CODE=1
     (...)
     # ------------------------------------------------------------------------------
     # OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
     # overwrite: delete the existing PACKAGE folder and install it again -----------
     # redownload: delete the compressed file and download data again ---------------
     # These keys are only relevant if you run setup_cocoa multiple times -----------
     # ------------------------------------------------------------------------------
     (...)
     export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible loss of uncommitted work)
                                                # If unset, users must manually delete cosmolike projects
     (...)
     # ------------------------------------------------------------------------------
     # Cosmolike projects below -----------------------------------------------------
     # ------------------------------------------------------------------------------
     (...)
     export ROMAN_REAL_URL="https://github.com/CosmoLike/cocoa_roman_real.git"
     export ROMAN_REAL_NAME="roman_real"
     #BRANCH: if unset, load the latest commit on the specified branch
     #export ROMAN_REAL_BRANCH="dev"
     #COMMIT: if unset, load the specified commit
     export ROMAN_REAL_COMMIT="a5cf62ffcec7b862dda5bf343bf6bb19124bb5d0"
     #TAG: if unset, load the specified TAG
     #export ROMAN_REAL_TAG="v4.0-beta17"
 
Every time `set_installation_options.sh` is edited, users need to reload `(.local)` by rerunning `start_cocoa.sh`. Therefore, to download and compile cosmolike projects users must first run the following commands:

      cd ./cocoa/Cocoa
      source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
      
 and then
 
      source ./installation_scripts/setup_cosmolike_projects.sh   # download all cosmolike projects  
      source ./installation_scripts/compile_all_projects.sh       # compile  all cosmolike project

In case users only want to compile a single cosmolike project (let's say the `roman_real` project)

      source ./projects/roman_real/scripts/compile_roman_real.sh

### Part II: New Projects 

Below, we provide instructions on how to download and install cosmolike projects that have not been previously configured on `set_installation_options.sh` shell script.

**Step :one:**: Activate conda environment and go to the projects folder (`Cocoa/projects`)
    
    conda activate cocoa
    cd ./projects  # here we assume users are located in the Cocoa's main folder that contains the `start_cocoa.sh` script

**Step :two:**  Clone the desired Cosmolike project (assumed here to have with the fictitious name `XXX`):

    git clone https://github.com/CosmoLike/cocoa_lsst_XXX.git XXX # users can clone their fork instead

> [!Note]
> By convention, the *Cosmolike Organization* adds the prefix `cocoa_` to all Cobaya-Cosmolike projects, which must be removed when cloning the repository. Why? The cocoa script `compile_all_projects.sh`, which finds and compiles all Cosmolike projects, requires the prefix to be removed.
 
**Step :three:**: Go back to the Cocoa main folder and activate the private Python environment
    
    cd ../
    source start_cocoa.sh # even if (.local) is loaded, users need to rerun start_cocoa after cloning the project repository
 
> [!Warning]
> Users must run the script `start_cocoa.sh` after cloning the project repository so that Cocoa can reload `(.local)` environment and create appropriate soft links to cobaya.

**Step :four:**: Compile the project (two possibilities)
 
    source ./projects/XXX/scripts/compile_XXX

or

    source ./installation_scripts/compile_all_projects.sh # compile all cosmolike projects
 
**Step :five:**: Select the number of OpenMP cores and run a template YAML file
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

and 

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f

If users want to make a particular Cosmolike project widely available in Cocoa, implement the following changes to Cocoa's configuration scripts:

**Step :one:**: Add the following env keys on `set_installation_options.sh`

    [adapted from Cocoa/set_installation_options.sh shell script]
    
    # ------------------------------------------------------------------------------
    # The keys below control which cosmolike projects will be installed and compiled
    # ------------------------------------------------------------------------------
    (...)
    #export IGNORE_COSMOLIKE_XXX_CODE=1
    (...)
    # ------------------------------------------------------------------------------
    # Cosmolike projects below -------------------------------------------
    # ------------------------------------------------------------------------------
    (...)
    export XXX_URL="https://github.com/.../cocoa_lsst_XXX.git"
    export XXX_NAME="XXX"
    #BRANCH: if unset, load the latest commit on the specified branch
    #export XXX_BRANCH="dev"
    #COMMIT: if unset, load the specified commit
    export XXX_COMMIT="abc123"
    #TAG: if unset, load the specified TAG
    #export XXX_TAG="v4.0-beta17"

**Step :two:**: Add all new defined environment keys to `flags_impl_unset_keys.sh` 
    
    [adapted from Cocoa/installation_scripts/flags_impl_unset_keys.sh]
    unset -v XXX_URL XXX_NAME XXX_COMMIT IGNORE_COSMOLIKE_XXX_CODE XXX_TAG XXX_BRANCH

This will ensure the bash script `stop_cocoa.sh` unsets these keys before unloading Cocoa's `(.local)` environment. Failing to do so will pollute your bash environment.

**Step :three:** Add and adapt the following block of code to the shell script `Cocoa/installation_scripts/setup_cosmolike_projects.sh`.

     [adapted from Cocoa/installation_scripts/setup_cosmolike_projects.sh shell script]
     if [ -z "${IGNORE_COSMOLIKE_XXX_CODE}" ]; then 
         PRINTNAME="XXX"
       
         ptop "GETTING ${PRINTNAME:?}" || return 1

         FOLDER="${XXX_NAME}"
         URL="${XXX_URL}"

         if [ -n "${XXX_COMMIT}" ]; then
             gitact2 "${FOLDER:?}" "${URL:?}" "${XXX_COMMIT:?}"  || return 1
         elif [ -n "${XXX_BRANCH}" ]; then 
             gitact1 "${FOLDER:?}" "${URL:?}" "${XXX_BRANCH:?}" || return 1
        elif [ -n "${XXX_TAG}" ]; then 
           gitact3 "${FOLDER:?}" "${URL:?}" "${XXX_TAG:?}" || return 1
        fi
      
        pbottom "GETTING ${PRINTNAME:?}" || return 1
    fi

> [!NOTE]
> Cocoa contains the scripts
>
>      Cocoa/installation_scripts/setup_cosmolike_projects.sh.  # download cosmolike projects set on set_installation_options.sh
>      Cocoa/installation_scripts/compile_all_projects.sh            # finds and compiles all cosmolike projects located on Cocoa/projects
> 
> designed to download and compile all Cosmolike projects defined in the `Cocoa/projects` folder.

> [!Warning]
> Never delete a folder from `projects` without first running `stop_cocoa.sh`; otherwise, Cocoa will have ill-defined links to these projects.

### :interrobang: FAQ: How do we set Weak Lensing YAML files in Cobaya? <a name="appendix_example_runs"></a>

The CosmoLike pipeline requires $\Omega_m$ and $\Omega_b$ to be provided, but the CAMB Boltzmann code only accepts $\Omega_c h^2$ and $\Omega_b h^2$ in Cobaya. Given that, there are two ways of creating YAML compatible with CAMB and Cosmolike: 

1. CMB parameterization: $\big(\Omega_c h^2,\Omega_b h^2\big)$ as primary MCMC parameters and $\big(\Omega_m,\Omega_b\big)$ as derived quantities.

        omegabh2:
            prior:
                min: 0.01
                max: 0.04
            ref:
                dist: norm
                loc: 0.022383
                scale: 0.005
            proposal: 0.005
            latex: \Omega_\mathrm{b} h^2
        omegach2:
            prior:
                min: 0.06
                max: 0.2
            ref:
                dist: norm
                loc: 0.12011
                scale: 0.03
            proposal: 0.03
            latex: \Omega_\mathrm{c} h^2
        mnu:
            value: 0.06
            latex: m_\\nu
        omegam:
            latex: \Omega_\mathrm{m}
        omegamh2:
            derived: 'lambda omegam, H0: omegam*(H0/100)**2'
            latex: \Omega_\mathrm{m} h^2
        omegab:
            derived: 'lambda omegabh2, H0: omegabh2/((H0/100)**2)'
            latex: \Omega_\mathrm{b}
        omegac:
            derived: 'lambda omegach2, H0: omegach2/((H0/100)**2)'
            latex: \Omega_\mathrm{c}

2. Weak Lensing parameterization: $\big(\Omega_m,\Omega_b\big)$ as primary MCMC parameters and $\big(\Omega_c h^2, \Omega_b h^2\big)$ as derived quantities.

        omegab:
            prior:
                min: 0.03
                max: 0.07
            ref:
                dist: norm
                loc: 0.0495
                scale: 0.004
            proposal: 0.004
            latex: \Omega_\mathrm{b}
            drop: true
        omegam:
            prior:
                min: 0.1
                max: 0.9
            ref:
                dist: norm
                loc: 0.316
                scale: 0.02
            proposal: 0.02
            latex: \Omega_\mathrm{m}
            drop: true
        mnu:
            value: 0.06
            latex: m_\\nu
        omegabh2:
            value: 'lambda omegab, H0: omegab*(H0/100)**2'
            latex: \Omega_\mathrm{b} h^2
        omegach2:
            value: 'lambda omegam, omegab, mnu, H0: (omegam-omegab)*(H0/100)**2-(mnu*(3.046/3)**0.75)/94.0708'
            latex: \Omega_\mathrm{c} h^2
        omegamh2:
            derived: 'lambda omegam, H0: omegam*(H0/100)**2'
            latex: \Omega_\mathrm{m} h^2

> [!Warning]
>Adopting $\big(\Omega_m,\Omega_b\big)$ as main MCMC parameters can create a silent bug in Cobaya. We are unsure if this problem persists in newer Cobaya versions; therefore, users are advised to follow our instructions. The problem occurs when the option `drop: true` is absent in $\big(\Omega_m,\Omega_b\big)$ parameters, and there are no expressions that define the derived $\big(\Omega_c h^2, \Omega_b h^2\big)$ quantities. The bug is silent because the MCMC runs without any warnings, but the CAMB Boltzmann code does not update the cosmological parameters at every MCMC iteration. As a result, the posteriors are flawed, but they may seem reasonable to those unfamiliar with the issue. 

### :interrobang: FAQ: How do we set Slow/Fast decomposition with Cosmolike?  <a name="manual_blocking_cosmolike"></a>

Cosmolike can't cache the intermediate products associated with the previous two evaluations, which are necessary to exploit Cobaya optimizations associated with dragging (`drag: True`). Still, it can cache the last evaluation, allowing the user to take advantage of a slow/fast decomposition of parameters in Cobaya's main Markov Chain Monte Carlo (MCMC) sampler. 

Cobaya cannot automatically handle parameters with different speed hierarchies associated with the same likelihood. Luckily, we can manually impose a speed hierarchy in Cobaya using the `blocking:` option in the YAML file. The main limitation of this method is that parameters of all adopted likelihoods, not only the ones required by Cosmolike, must be manually specified.

Below, we provide an example YAML configuration for an MCMC chain with DES 3x2pt likelihood.

        likelihood: 
            des_y3.des_3x2pt:
            path: ./external_modules/data/des_y3
            data_file: DES_Y1.dataset
         (...)
        sampler:
            mcmc:
                covmat: "./projects/des_y3/EXAMPLE_MCMC22.covmat"
                covmat_params:
                # ---------------------------------------------------------------------
                # Proposal covariance matrix learning
                # ---------------------------------------------------------------------
                learn_proposal: True
                learn_proposal_Rminus1_min: 0.03
                learn_proposal_Rminus1_max: 100.
                learn_proposal_Rminus1_max_early: 200.
                # ---------------------------------------------------------------------
                # Convergence criteria
                # ---------------------------------------------------------------------
                # Maximum number of posterior evaluations
                max_samples: .inf
                # Gelman-Rubin R-1 on means
                Rminus1_stop: 0.015
                # Gelman-Rubin R-1 on std deviations
                Rminus1_cl_stop: 0.17
                Rminus1_cl_level: 0.95
                # ---------------------------------------------------------------------
                # Exploiting Cosmolike speed hierarchy
                # ---------------------------------------------------------------------
                measure_speeds: False 
                drag: False
                oversample_power: 0.2
                oversample_thin: True
                blocking:
                - [1,
                    [
                        As_1e9, H0, omegab, omegam, ns
                    ]
                  ]
                - [4,
                    [
                        DES_DZ_S1, DES_DZ_S2, DES_DZ_S3, DES_DZ_S4, DES_A1_1, DES_A1_2,
                        DES_B1_1, DES_B1_2, DES_B1_3, DES_B1_4, DES_B1_5,
                        DES_DZ_L1, DES_DZ_L2, DES_DZ_L3, DES_DZ_L4, DES_DZ_L5
                    ]
                  ]
                - [25,
                    [
                        DES_M1, DES_M2, DES_M3, DES_M4, DES_PM1, DES_PM2, DES_PM3, DES_PM4, DES_PM5
                    ]
                  ]
                # ---------------------------------------------------------------------
                max_tries: 100000
                burn_in: 0
                Rminus1_single_split: 4


## :interrobang: FAQ: How do we create a new Cosmolike project <a name="appendix_lsst_y1_new"></a> 

Adapting the LSST_Y1 folder to construct a new project involves many small core changes and a few major ones. They are tedious but straightforward. The easier way to apply the minor core changes to the code is via the bash script *transfer_project.sh*.

## The easy way <a name="appendix_lsst_y1_new_small"></a> 

 To properly use the bash script *transfer_project.sh*., users must set the following variables at the beginning of the file:

     OLD_PROJECT="lsst_y1"
     OLD_SURVEY="LSST"

     NEW_PROJECT="des_y3"
     NEW_SURVEY="DES"

After that, type

    conda activate cocoa
    source start_cocoa.sh
    cd ./projects

and

    bash transfer_project.sh

## The hard way <a name="appendix_lsst_y1_new_small2"></a> 

### Create the new project

**Step 1:** Initialize cocoa environments.

    cd ./cocoa/Cocoa

and
	
    conda activate cocoa

and

    source start_cocoa.sh

**Step 2:** Choose a project name (e.g., project `xxx`), and copy the `LSST_Y1` project using the command below
    
    cp -r "${ROOTDIR:?}"/projects/lsst_y1/ "${ROOTDIR:?}"/projects/xxx

**Step 3:** Remove the git repository associated with LSST_Y1 project

    rm -rf "${ROOTDIR:?}"/projects/xxx/.git/

### Changes in the `Cocoa/projects/xxx/interface` folder 

**Step 1:** Change the file `Cocoa/projects/xxx/interface/MakefileCosmolike` following the instructions below (replace all instances of `lsst_y1` with `xxx`

    [adapted from Cocoa/projects/lsst_y1/interface/MakefileCosmolike line ~118]
    shared: cosmolike_lsst_y1_interface.so # delete this line
    shared: cosmolike_xxx_interface.so     # add this line

and

    [adapted from Cocoa/projects/lsst_y1/interface/MakefileCosmolike line ~172]
    cosmolike_lsst_y1_interface.so: $(OBJECTC) $(CSOURCES) interface.cpp # delete this line
    cosmolike_xxx_interface.so: $(OBJECTC) $(CSOURCES) interface.cpp     # add this line
        $(CXX) $(CXXFLAGS) -DCOBAYA_SAMPLER -shared -fPIC -o $@ $(OBJECTC) interface.cpp $(LDFLAGS)
        @rm *.o
    
and
    
     [adapted from Cocoa/projects/lsst_y1/interface/MakefileCosmolike line ~176]
     clean:
        @rm -rf cosmolike_lsst_y1_interface.so cosmolike_lsst_y1_interface.so.dSYM  *.o # delete this line
        @rm -rf cosmolike_xxx_interface.so cosmolike_xxx_interface.so.dSYM  *.o         # add this line

**Step 2:** Change the name of the file `cosmolike_lsst_y1_interface.py` using the command below
       
    mv "${ROOTDIR:?}"/projects/xxx/interface/cosmolike_lsst_y1_interface.py "${ROOTDIR:?}"/projects/xxx/interface/cosmolike_xxx_interface.py

and remove any previously compiled dynamic library 

    rm -f  "${ROOTDIR:?}"/projects/xxx/interface/cosmolike_lsst_y1_interface.so
    
**Step 3** Change the newly created file `cosmolike_xxx_interface.py` following the instructions below

	[adapted from Cocoa/projects/lsst_y1/interface/cosmolike_xxx_interface.py]
    def __bootstrap__():
        (...)
        __file__ = pkg_resources.resource_filename(__name__,'cosmolike_lsst_y1_interface.so') # delete this line
        __file__ = pkg_resources.resource_filename(__name__,'cosmolike_xxx_interface.so')     # add this line
        
**Step 4** Change the file `Cocoa/projects/XXX/interface/interface.cpp` following the instructions below
    
    [adapted from Cocoa/projects/lsst_y1/interface/interface.cpp line ~43]
    PYBIND11_MODULE(cosmolike_lsst_y1_interface, m) # delete this line
    PYBIND11_MODULE(cosmolike_xxx_interface, m)     # add this line
    {
        m.doc() = "CosmoLike Interface for LSST_Y1 3x2pt Module"; # delete this line
        m.doc() = "CosmoLike Interface for XXX 3x2pt Module";     # add this line   
       (...)
    
### Changes in the `Cocoa/projects/xxx/scripts` folder

**Step 1:** Change the name of the file `compile_lsst_y1.sh` using the command below 
    
    mv "${ROOTDIR:?}"/projects/xxx/scripts/compile_lsst_y1.sh "${ROOTDIR:?}"/projects/xxx/scripts/compile_xxx.sh
    
**Step 2:** Change the name of the file `start_lsst_y1.sh` using the command below 
    
    mv "${ROOTDIR:?}"/projects/xxx/scripts/start_lsst_y1.sh "${ROOTDIR:?}"/projects/xxx/scripts/start_xxx.sh
    
**Step 3:** Change the name of the file `stop_lsst_y1.sh` using the command below 
    
    mv "${ROOTDIR:?}"/projects/xxx/scripts/stop_lsst_y1.sh "${ROOTDIR:?}"/projects/xxx/scripts/stop_xxx.sh

**Step 4:** Change the file `compile_xxx.sh` following the instructions below

    [adapted from Cocoa/projects/lsst_y1/scripts/compile_lsst_y1.sh line ~4;43;65]
    if [ -z "${IGNORE_COSMOLIKE_LSSTY1_CODE}" ]; then # delete this line
    if [ -z "${IGNORE_COSMOLIKE_XXX_CODE}" ]; then    # add this line
        (...)
	FOLDER="${LSST_Y1_NAME:-"lsst_y1"}" # delete this line
        FOLDER="${XXX_NAME:-"xxx"}"         # add this line
        (...)
        PRINTNAME="LSST_Y1"  # delete this line
        PRINTNAME="XXX"      # add this line
     
**Step 5:** Change the file `start_xxx.sh` following the instructions below

    [adapted from Cocoa/projects/lsst_y1/scripts/start_lsst_y1.sh line ~4;13]
    if [ -z "${IGNORE_COSMOLIKE_LSSTY1_CODE}" ]; then # delete this line
    if [ -z "${IGNORE_COSMOLIKE_XXX_CODE}" ]; then    # add this line
        (...)
	FOLDER="${LSST_Y1_NAME:-"lsst_y1"}" # delete this line
        FOLDER="${XXX_NAME:-"xxx"}"         # add this line

### Changes in the `Cocoa/projects/xxx/likelihood` folder

**Step 1:** Change the file `_cosmolike_prototype_base.py` following the instructions below

    [adapted from Cocoa/projects/lsst_y1/likelihood/_cosmolike_prototype_base.sh line ~19;21]
    import cosmolike_lsst_y1_interface as ci #delete this line
    import cosmolike_xxx_interface as ci     #add this line

    survey = "LSST"  #delete this line
    survey = "xxx"   # add this line

> [!Tip]
> If the project name `xxx` contains more than the experiment name (e.g., the release year), we suggest assigned `survey` to just the experiment name. For example, if `XXX = DES_Y3`, then assigned `survey = "DES"`.

**Step 2:** Change the file `combo_3x2pt.py` following the instructions below

    [adapted from Cocoa/projects/lsst_y1/likelihood/combo_3x2pt.py line ~1-10;]
    from cobaya.likelihoods.lsst_y1._cosmolike_prototype_base import _cosmolike_prototype_base, survey # delete this line
    from cobaya.likelihoods.xxx._cosmolike_prototype_base import _cosmolike_prototype_base, survey     # add this line
    
    import cosmolike_lsst_y1_interface as ci # delete this line
    import cosmolike_xxx_interface as ci     # add this line

Users should perform similar changes to `combo_2x2pt.py`, `combo_xi_gg.py`, `combo_xi_ggl.py`, and `cosmic_shear.py`.
    
**Step 3:** Change the file `combo_3x2pt.yaml` following the instructions below
   
    [adapted from Cocoa/projects/lsst_y1/likelihood/combo_3x2pt.py lines ~2;30;31]
    data_file: lsst_y1_M1_GGL0.05.dataset #delete this line
    data_file: xxx_yyy.dataset            #add and adapt this line; yyy = the adopted scale cuts
    (...)
    filename_baryon_pca: "./projects/lsst_y1/data/pca.txt" #delete this line
    filename_baryon_pca: "./projects/xxx/data/pca.txt"     #add this line; users will need to recompute pca.txt if they want to use PCA for baryons
    (...)
    print_datavector_file: "./projects/lsst_y1/chains/lsst_y1_theory.modelvector" #delete this line
    print_datavector_file: "./projects/xxx/chains/xxx_theory.modelvector"         #add this line

Users should perform similar changes to `combo_2x2pt.yaml`, `combo_xi_gg.yaml`, `combo_xi_ggl.yaml`, and `cosmic_shear.yaml`.

**Step 4:** Rename the prefix `LSST_` of all lens and source-related parameters, located on `params_lens.yaml` and `params_source.yaml`, as shown below. 

    [adapted from Cocoa/projects/lsst_y1/likelihood/params_lens.yaml lines ~2-12]
    LSST_DZ_L1:  # delete this line
    XX_DZ_L1:    # add this line
        prior:
            dist: norm
            loc: 0.0
            scale: 0.005
        ref:
            dist: norm
            loc: 0.0
            scale: 0.005
        proposal: 0.005
        latex: \Delta z_\mathrm{l,LSST}^1 # delete this line
	latex: \Delta z_\mathrm{l,XXX}^1 # delete this line
     
     (...) # Don't forget to rename ALL parameters
         
### Changes in the `Cocoa/projects/xxx/data` folder

**Step 1:** Rename the `.dataset` file by adapting the command below 

    # yyy = the adopted scale cuts on xxx_yyy.dataset
    mv "${ROOTDIR:?}"/projects/xxx/data/lsst_y1_M1_GGL0.05.dataset $ROOTDIR/projects/xxx/data/xxx_yyy.dataset

> [!Tip]
> There are many datasets and masks associated with different scale cuts in the `lsst_y1/data` folder. Some of them are listed below.
>
>    lsst_y1_M[2-6]_GGL0.05.dataset
>    lsst_y1_M[2-6]_GGLOLAP0.05.mask
>
> We recommend that users delete these files, as they will not be applicable in the new projects. The same goes for data vectors, covariances, and n(z) files.

**Step 2:** Update `xxx_yyy.dataset` file with the names of the new data vector, covariance, n(z), binning, mask...

    [adapted from Cocoa/projects/lsst_y1/likelihood/params_lens.yaml lines ~2-12]
    data_file = XXX_nonlim
    cov_file = cov_XXX
    mask_file = 3x2pt_baseline.mask
    nz_lens_file = lens_XXX.nz
    nz_source_file = source_XXX.nz
    lens_ntomo = 5
    source_ntomo = 5
    n_theta = 26
    theta_min_arcmin = 2.5
    theta_max_arcmin = 900.
    baryon_pca_file = pca.txt

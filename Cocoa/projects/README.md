# Table of contents
1. [The Projects Folder](#appendix_projects_folder)
2. [FAQ: How do we download and run Cosmolike projects?](#running_cosmolike_projects)
3. [FAQ: How do we set Weak Lensing YAML files in Cobaya?](#appendix_example_runs)
4. [FAQ: How do we set Slow/Fast decomposition with Cosmolike?](#manual_blocking_cosmolike)
5. [FAQ: How do we create a new Cosmolike project?](#appendix_lsst_y1_new)
   1. [Minor changes: the easy way](#appendix_lsst_y1_new_small)
   2. [Minor changes: the hard way](#appendix_lsst_y1_new_small2)
   3. [Major changes](#appendix_lsst_y1_new_major)
 
### The Projects Folder <a name="appendix_projects_folder"></a> 

The `projects` folder includes all the projects linked to Cosmolike; they can also help organize general investigations, even if they don't use Cosmolike directly. 

Projects should be hosted on independent GitHub repositories; our convention is to name the repository cocoa_XXX, where XXX is the intended project name. Projects that utilize Cosmolike need to have more or less the following structure, taken from the [LSST_Y1 project](https://github.com/CosmoLike/cocoa_lsst_y1)

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

### :interrobang: FAQ: How do we download and run Cosmolike projects? <a name="running_cosmolike_projects"></a> 

**Step :one:**: Activate conda environment and go to the projects folder (`./projects`)
    
    conda activate cocoa
    cd ./cocoa/Cocoa/projects

**Step :two:**  Clone a Cosmolike project (assumed here to have with the fictitious name `XXX`):

    git clone https://github.com/CosmoLike/cocoa_lsst_XXX.git XXX

By convention, the *Cosmolike Organization* adds the prefix `cocoa_` to all Cobaya-Cosmolike projects, which must be removed when cloning the repository.
 
**Step :three:**: Go back to the Cocoa main folder and activate the private Python environment
    
    cd ../

and

    source start_cocoa.sh
 
> [!Warning]
> Users must run the script `start_cocoa.sh` after cloning the project repository so that Cocoa can reload `(.local)` environment and create appropriate soft links.

**Step :four:**: Compile the project, as shown below (two possibilities)
 
    source ./projects/XXX/scripts/compile_XXX

or

    source ./installation_scripts/compile_all_projects.sh # compile all cosmolike projects
 
**Step :five:**: Select the number of OpenMP cores and run a template YAML file
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

and 

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f

If users want to make a particular Cosmolike project widely available in Cocoa, implement the following steps:

**Step :one:**: Add the following env keys on `set_installation_options.sh`

    [adapted from Cocoa/set_installation_options.sh script]
    
    #export IGNORE_COSMOLIKE_XXX_CODE=1

    (...)
   
    export XXX_URL="https://github.com/.../cocoa_lsst_XXX.git"
    export XXX_NAME="XXX"
    #BRANCH: if unset, load the latest commit on the specified branch
    #export XXX_BRANCH="dev"
    #COMMIT: if unset, load the specified commit
    export XXX_COMMIT="abc123"
    #TAG: if unset, load the specified TAG
    #export XXX_TAG="v4.0-beta17"

**Step :two:**: Adapt and add the new defined keys to `flags_impl_unset_keys.sh` 

    [adapted from Cocoa/installation_scripts/flags_impl_unset_keys.sh]

    unset -v XXX_URL XXX_NAME XXX_COMMIT IGNORE_COSMOLIKE_XXX_CODE XXX_TAG XXX_BRANCH

This will ensure that `stop_cocoa.sh` unsets them before exiting Cocoa.

**Step :three:** Add and adapt the following block to `setup_cosmolike_projects.sh`.

     [adapted from Cocoa/installation_scripts/setup_cosmolike_projects.sh script]
     
     if [ -z "${IGNORE_COSMOLIKE_XXX_CODE}" ]; then 
       ptop "GETTING XXX" || return 1

       if [ -n "${XXX_COMMIT}" ]; then
         gitact2 "${XXX_NAME:?}" "${XXX_URL:?}" "${XXX_COMMIT:?}"  || return 1
       fi

       pbottom "GETTING XXX" || return 1
     fi

> [!NOTE]
> Cocoa contains the scripts
>
>      Cocoa/installation_scripts/setup_cosmolike_projects.sh
>      Cocoa/installation_scripts/compile_all_projects.sh
> 
> designed to download and compile all Cosmolike projects defined in the `Cocoa/projects` folder.

> [!Warning]
> Never delete a folder from `projects` without first running `stop_cocoa.sh`; otherwise, Cocoa will have ill-defined links.

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
>Adopting $\big(\Omega_m,\Omega_b\big)$ as main MCMC parameters can create a silent bug in Cobaya. We are unsure if this problem persists in newer Cobaya versions; therefore, users must follow our instructions. The problem occurs when the option `drop: true` is absent in $\big(\Omega_m,\Omega_b\big)$ parameters, and there are no expressions that define the derived $\big(\Omega_c h^2, \Omega_b h^2\big)$ quantities. The bug is silent because the MCMC runs without any warnings, but the CAMB Boltzmann code does not update the cosmological parameters at every MCMC iteration. As a result, the posteriors are flawed, but they may seem reasonable to those unfamiliar with the issue. 

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

## Minor changes: the easy way <a name="appendix_lsst_y1_new_small"></a> 

 To properly use the bash script *transfer_project.sh*., users must set the following variables at the beginning of the file:

     OLD_PROJECT="lsst_y1"
     OLD_SURVEY="LSST"

     NEW_PROJECT="des_y3"
     NEW_SURVEY="DES"

After that, just type

     $(cocoa)(.local) bash transfer_project.sh

## Minor changes: the hard way <a name="appendix_lsst_y1_new_small2"></a> 

### Create the new project

**Step 1:** Choose a project name (e.g., project XXX), and copy the `LSST_Y1` project using the command below
    
    $(cocoa)(.local) cp $ROOTDIR/projects/lsst_y1/ $ROOTDIR/projects/xxx

**Step 2:** Remove the git repository associated with LSST_Y1 project

    $(cocoa)(.local) rm -rf $ROOTDIR/projects/$NEW_PROJECT/.git/

### Changes in the interface folder 

**Step 1:** Change the file `$ROOTDIR/projects/XXX/interface/MakefileCosmolike` following the instructions below

    (...)
     
    CSOURCES += \
        (...)
        ${ROOTDIR}/external_modules/code/cosmolike/pt_cfastpt.c \
        // add additional files from /external_modules/code/cosmolike that is needed
    
    (...)
    
    OBJECTC += \
        (...)
        ./pt_cfastpt.o \
        // add additional files from /external_modules/code/cosmolike that is needed
    
    (...)
    
    all:  shared
    // change cosmolike_lsst_y1_interface.so to cosmolike_XXX_interface.so in the line below
    shared: cosmolike_lsst_y1_interface.so
    
    (...)
    
    // change cosmolike_lsst_y1_interface.so to cosmolike_XXX_interface.so in the line below
    cosmolike_lsst_y1_interface.so: $(OBJECTC) $(CSOURCES) interface.cpp
        $(CXX) $(CXXFLAGS) -DCOBAYA_SAMPLER -shared -fPIC -o $@ $(OBJECTC) interface.cpp $(LDFLAGS)
        @rm *.o
    
    (...)
    
    // change cosmolike_lsst_y1_interface.so to cosmolike_XXX_interface.so in the line below 
    clean:
        @rm -rf cosmolike_lsst_y1_interface.so cosmolike_lsst_y1_interface.so.dSYM  *.o

**Step 2:** Change the name of the File `$ROOTDIR/projects/XXX/interface/cosmolike_lsst_y1_interface.py` using the command below
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/interface/cosmolike_lsst_y1_interface.py $ROOTDIR/projects/XXX/interface/cosmolike_XXX_interface.py

**Step 3** Changes in the newly created file `$ROOTDIR/projects/XXX/interface/cosmolike_XXX_interface.py` 

    def __bootstrap__():
        (...)
        // change cosmolike_lsst_y1_interface.so to cosmolike_XXX_interface.so in the line below 
        __file__ = pkg_resources.resource_filename(__name__,'cosmolike_lsst_y1_interface.so')
        
**Step 4** Change the file `$ROOTDIR/projects/XXX/interface/interface.cpp` following the instructions below
    
    (...)
    
    // change cosmolike_lsst_y1_interface to cosmolike_XXX_interface in the line below
    PYBIND11_MODULE(cosmolike_lsst_y1_interface, m)
    {
        // change the description below
        m.doc() = "CosmoLike Interface for LSST_Y1 3x2pt Module";
        
       (...)
    }
    
### Changes in the script folder

**Step 1:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/compile_lsst_y1` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/compile_lsst_y1 $ROOTDIR/projects/XXX/scripts/compile_XXX
    
**Step 2:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/start_lsst_y1` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/start_lsst_y1 $ROOTDIR/projects/XXX/scripts/start_XXX
    
**Step 3:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/stop_lsst_y1` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/stop_lsst_y1 $ROOTDIR/projects/XXX/scripts/stop_XXX

**Step 4:** Change the file `$ROOTDIR/projects/XXX/scripts/compile_lsst_y1` following the instructions below

    (...)

    // change $ROOTDIR/projects/lsst_y1/interface to $ROOTDIR/projects/XXX/interface in the line below 
    cd $ROOTDIR/projects/lsst_y1/interface
    
**Step 5:** Change the file `$ROOTDIR/projects/XXX/scripts/start_lsst_y1` following the instructions below

    (...)

    // change $ROOTDIR/projects/lsst_y1/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar LD_LIBRARY_PATH $ROOTDIR/projects/lsst_y1/interface
    
    // change $ROOTDIR/projects/lsst_y1/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar PYTHONPATH $ROOTDIR/projects/lsst_y1/interface

### Changes in the likelihood folder

**Step 1:** Change the file `$ROOTDIR/projects/XXX/likelihood/_cosmolike_prototype_base.py` following the instructions below

    (...) 
    
    // change cosmolike_lsst_y1_interface to cosmolike_XXX_interface in the line below 
    import cosmolike_lsst_y1_interface as ci
    
    (...)
     
    def set_source_related(self, **params_values):
        ci.set_nuisance_shear_calib(
          M = [
            params_values.get(p, None) for p in [
              // change LSST_ to the name of the survey associated w/ XXX)
              "LSST_M"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_shear_photoz(
          bias = [
            params_values.get(p, None) for p in [
              // change LSST_ to the name of the survey associated w/ XXX)
              "LSST_DZ_S"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_ia(
          A1 = [
            params_values.get(p, None) for p in [
              // change LSST_ to the name of the survey associated w/ XXX)
              "LSST_A1_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          A2 = [
            params_values.get(p, None) for p in [
              // change LSST_ to the name of the survey associated w/ XXX)
              "LSST_A2_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          B_TA = [
            params_values.get(p, None) for p in [
              // change LSST_ to the name of the survey associated w/ XXX)
              "LSST_BTA_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
        )
     
    (...)
     
    def set_lens_related(self, **params_values):
        ci.set_nuisance_bias(
            B1 = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "LSST_B1_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ], 
            B2 = [
                  params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "LSST_B2_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ],
            B_MAG = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "LSST_BMAG_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_nuisance_clustering_photoz(
            bias = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "LSST_DZ_L"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_point_mass(
            PMV = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "LSST_PM"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
     
    (...)
     
    def set_baryon_related(self, **params_values):
        // change LSST_ to the name of the survey associated w/ XXX)
        self.baryon_pcs_qs[0] = params_values.get("LSST_BARYON_Q1", None)
        self.baryon_pcs_qs[1] = params_values.get("LSST_BARYON_Q2", None)
        self.baryon_pcs_qs[2] = params_values.get("LSST_BARYON_Q3", None)
        self.baryon_pcs_qs[3] = params_values.get("LSST_BARYON_Q4", None)

If the project name `XXX` contains more than the experiment name, we suggest replacing `LSST_` with just the experiment name. For example, if `XXX = DES_Y3`, then adopt `DES_DZ_L1` for the name of the redshift shift on lens bin 1. The convention adopted must be followed when changing the files `params_des_cosmic_shear.yaml` and `params_des_3x2pt.yaml`. 

**Step 2:** Change the file `$ROOTDIR/projects/XXX/likelihood/lsst_3x2pt.py` following the instructions below
    
    // change lsst_y1 to XXX in the line below
    from cobaya.likelihoods.lsst_y1._cosmolike_prototype_base import _cosmolike_prototype_base
    // change cosmolike_lsst_y1_interface to cosmolike_XXX_interface in the line below
    import cosmolike_lsst_y1_interface as ci
    
**Step 3:** Change the file `$ROOTDIR/projects/XXX/likelihood/lsst_3x2pt.yaml` following the instructions below
   
    (...)
    // change LSST_Y1.dataset to XXX.dataset in the line below (adopted convention: .dataset file name = project name all in CAPS)
    data_file: LSST_Y1.dataset
    (...)
    // change params_lsst_3x2pt to params_XXX_3x2pt in the line below
    params: !defaults [params_lsst_3x2pt]

**Step 4:** Rename the file `params_lsst_3x2pt.yaml` to `params_XXX_3x2pt.yaml`. Also, rename the associated parameter names, 
replacing the `LSST_` prefix as shown below. 

      XXX_DZ_S1:
         prior:
            dist: norm
            loc: 0.0
            scale: 0.005
         ref:
            dist: norm
            loc: 0.0
            scale: 0.005
            proposal: 0.005
         latex: \Delta z_\mathrm{s, XXX}^1
         
Similar changes must be made in `params_XXX_cosmic_shear.yaml`. Note that changes either in the number of lenses or source bins will demand the introduction of new parameters in 
`params_XXX_cosmic_shear.yaml` and `params_XXX_3x2pt.yaml`

### Changes in the data folder

**Step 1:** Rename the `.dataset` file. Our adopted convention is: `.dataset` file name = project name capitalized

     $(cocoa)(.local) cp $ROOTDIR/projects/LSST_Y1/data/LSST_Y1.dataset $ROOTDIR/projects/XXX/data/XXX.dataset
     
**Step 2:** Update `XXX.dataset` file with the names of the new data vector, covariance, n(z), binning, mask...

      data_file = XXX_nonlim
      cov_file = cov_XXX
      mask_file = 3x2pt_baseline.mask
      nz_lens_file = lens_XXX.nz
      nz_source_file = source_XXX.nz
      lensing_overlap_cut = 0.0
      lens_ntomo = 5
      source_ntomo = 5
      n_theta = 26
      IA_model = 4
      theta_min_arcmin = 2.5
      theta_max_arcmin = 900.
      #baryon_pca_file = pca.txt

## Major changes: <a name="appendix_lsst_y1_new_major"></a>  

* Computation of a new covariance matrix using either [CosmoCov](https://github.com/CosmoLike/CosmoCov) or [CosmoCovFourier](https://github.com/CosmoLike/CosmoCov_Fourier)
* Simulation of new `n(z)` for lenses and sources
* Updates to the Cosmolike C++ interface so the appropriate routines can be called from the Python likelihood
* Updates to the Cosmolike Python likelihoods and their associated Yaml files. These include, for example, `/likelihood/lsst_3x2pt.py` and `/likelihood/lsst_3x2pt.yaml`

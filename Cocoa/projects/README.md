# Table of contents
1. [The Projects Folder](#appendix_projects_folder)
2. [Adapting the COCOA_LSST_Y1 repository to a new project](#appendix_lsst_y1_new)
   1. [Minor changes: the easy way](#appendix_lsst_y1_new_small)
   2. [Minor changes: the hard way](#appendix_lsst_y1_new_small2)
   3. [Major changes](#appendix_lsst_y1_new_major)
 
# The Projects Folder <a name="appendix_projects_folder"></a> 

The `projects` folder includes all the projects linked to Cosmolike; they can also help organize general investigations even if they don't use Cosmolike directly. 

Projects should be hosted on independent GitHub repositories; our convention is to name the repository cocoa_XXX, where XXX is the intended project name. Projects that utilize Cosmolike need to have more or less the following structure, taken from a [LSST_Y1 project](https://github.com/CosmoLike/cocoa_lsst_y1)

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

# Adapting the COCOA_LSST_Y1 repository to a new project <a name="appendix_lsst_y1_new"></a> 

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

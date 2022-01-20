# Table of contents
1. [The Projects Folder](#appendix_projects_folder)
2. [Adapting DES-Y3 to a new project](#appendix_des_y3_new)
3. [Minor core changes when adapting DES-Y3 to a new project - the easy way](#appendix_des_y3_new_small)
3. [Minor core changes when adapting DES-Y3 to a new project - the hard way](#appendix_des_y3_new_small2)
 
## The Projects Folder <a name="appendix_projects_folder"></a> 

The `projects` folder includes all the projects that our group is developing. Individual projects must be hosted on independent folders named `cocoa_XXX` where XXX is the project name. The majority of projects we are working on are not public (yet). However, the backbone Cosmolike software is publicly available at `external_modules/code`.

The `cocoa_XXX` folder that host the `XXX` project needs to have more or less the following structure (taken from our private DES-Y3 project)

    +-- cocoa_des_y3
    |    +-- likelihood
    |    |   +-- _cosmolike_prototype_base.py
    |    |   +-- des_3x2pt.py
    |    |   +-- des_3x2pt.yaml
    |    |   +-- des_2x2pt.py
    |    |   +-- des_2x2pt.yaml
    |    |   +-- des_clustering.py
    |    |   +-- des_clustering.yaml
    |    |   +-- des_cosmic_shear.py
    |    |   +-- des_cosmic_shear.yaml
    |    |   +-- des_ggl.py
    |    |   +-- des_ggl.yaml
    |    |   +-- des_xi_ggl.py
    |    |   +-- des_xi_ggl.yaml
    |    +-- scripts
    |    |   +-- compile_des_y3
    |    |   +-- start_des_y3
    |    |   +-- stop_des_y3
    |    +-- data
    |    |   +-- DES.paramnames
    |    |   +-- DES_Y3.dataset
    |    |   +-- datavector.txt
    |    |   +-- covariance.txt
    |    |   +-- nzlens.txt
    |    |   +-- nzsource.txt
    |    |   +-- mask.mask
    |    +-- interface
    |    |   +-- MakefileCosmolike
    |    |   +-- cosmolike_des_y3_interface.py
    |    |   +-- interface.cpp
    |    |   +-- interface.hpp
    |    +-- chains
    |    |   +-- README

## Adapting DES-Y3 to a new project <a name="appendix_des_y3_new"></a> 

(**warning**) The DES-Y3 project is not public yet, but our group will be release the code soon. 

Adapting the DES-Y3 folder to construct a new project involves many small core changes and a few major ones. **All minor core changes have been automatized in the script [transfer_project.sh](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/projects/transfer_project.sh)**.

**The Major changes are:**

* Computation of the covariance matrix using either [CosmoCov](https://github.com/CosmoLike/CosmoCov) or [CosmoCovFourier](https://github.com/CosmoLike/CosmoCov_Fourier) (replacing `/projects/des_y3/data/cov_Y3.txt`)
* Simulation of new `n(z)` for lenses and sources (replacing `/projects/des_y3/data/nz_lens_Y3.txt` and `/projects/des_y3/data/nz_source_Y3.txt`)
* Changes to the Cosmolike C++ interface so the appropriate routines can be called from the Python likelihood (will your project compute `3x2pt`, `6x2pt`, `4x2pt+N` or what?)
* Changes to the Cosmolike Python likelihood so `Cobaya` can call the appropriate routines
* Additional changes in the files located at `/data`, including the `DES_Y3.dataset`
* Changes to the number of lens and source bins in `.dataset` and in `params_XXX_3x2pt` files (and any additional files where the nuisance parameters are listed)

We list below the long list of small core changes so the C - C++ - Python interface can work flawlessly. They are tedious but straightforward. **All core changes have been automatized in the script [transfer_project.sh](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/projects/transfer_project.sh)**

## Minor core changes when adapting DES-Y3 to a new project - the easy way <a name="appendix_des_y3_new_small"></a> 

The easier way to create a new project and apply the many minor core changes to the code is via the bash script [transfer_project.sh](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/projects/transfer_project.sh). To proper use the bash script, users must set the following variables (set at the beginning of the [transfer_project.sh](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/projects/transfer_project.sh) file):

     OLD_PROJECT="des_y3"
     OLD_SURVEY="DES"

     NEW_PROJECT="lsst_y1"
     NEW_SURVEY="LSST"

After that, just type

     $(cocoa)(.local) bash transfer_project.sh

## Minor core changes when adapting DES-Y3 to a new project - the hard way <a name="appendix_des_y3_new_small2"></a> 

### Create the new project

**Step 1:** Choose a project name (e.g., project XXX), and copy the DES-Y3 project using the command below
    
    $(cocoa)(.local) cp $ROOTDIR/projects/des_y3/ $ROOTDIR/projects/xxx

**Step 2:** Remove the git repository associated with DES-Y3 project

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
    // change cosmolike_des_y3_interface.so to cosmolike_XXX_interface.so in the line below
    shared: cosmolike_des_y3_interface.so
    
    (...)
    
    // change cosmolike_des_y3_interface.so to cosmolike_XXX_interface.so in the line below
    cosmolike_des_y3_interface.so: $(OBJECTC) $(CSOURCES) interface.cpp
        $(CXX) $(CXXFLAGS) -DCOBAYA_SAMPLER -shared -fPIC -o $@ $(OBJECTC) interface.cpp $(LDFLAGS)
        @rm *.o
    
    (...)
    
    // change cosmolike_des_y3_interface.so to cosmolike_XXX_interface.so in the line below 
    clean:
        @rm -rf cosmolike_des_y3_interface.so cosmolike_des_y3_interface.so.dSYM  *.o

**Step 2:** Change the name of the File `$ROOTDIR/projects/XXX/interface/cosmolike_des_y3_interface.py` using the command below
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/interface/cosmolike_des_y3_interface.py $ROOTDIR/projects/XXX/interface/cosmolike_XXX_interface.py

**Step 3** Changes in the newly created file `$ROOTDIR/projects/XXX/interface/cosmolike_XXX_interface.py` 

    def __bootstrap__():
        (...)
        // change cosmolike_des_y3_interface.so to cosmolike_XXX_interface.so in the line below 
        __file__ = pkg_resources.resource_filename(__name__,'cosmolike_des_y3_interface.so')
        
**Step 4** Change the file `$ROOTDIR/projects/XXX/interface/interface.cpp` following the instructions below
    
    (...)
    
    // change cosmolike_des_y3_interface to cosmolike_XXX_interface in the line below
    PYBIND11_MODULE(cosmolike_des_y3_interface, m)
    {
        // change the description below
        m.doc() = "CosmoLike Interface for DES-Y3 3x2 Module";
        
       (...)
    }
    
### Changes in the script folder

**Step 1:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/compile_des_y3` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/compile_des_y3 $ROOTDIR/projects/XXX/scripts/compile_XXX
    
**Step 2:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/start_des_y3` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/start_des_y3 $ROOTDIR/projects/XXX/scripts/start_XXX
    
**Step 3:** Change the name of the file `$ROOTDIR/projects/XXX/scripts/stop_des_y3` using the command below 
    
    $(cocoa)(.local) mv $ROOTDIR/projects/XXX/scripts/stop_des_y3 $ROOTDIR/projects/XXX/scripts/stop_XXX

**Step 4:** Change the file `$ROOTDIR/projects/XXX/scripts/compile_des_y3` following the instructions below

    (...)

    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    cd $ROOTDIR/projects/des_y3/interface
    
**Step 5:** Change the file `$ROOTDIR/projects/XXX/scripts/start_des_y3` following the instructions below

    (...)

    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar LD_LIBRARY_PATH $ROOTDIR/projects/des_y3/interface
    
    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar PYTHONPATH $ROOTDIR/projects/des_y3/interface

### Changes in the likelihood folder

**Step 1:** Change the file `$ROOTDIR/projects/XXX/likelihood/_cosmolike_prototype_base.py` following the instructions below

    (...) 
    
    // change cosmolike_des_y3_interface to cosmolike_XXX_interface in the line below 
    import cosmolike_des_y3_interface as ci
    
    (...)
     
    def set_source_related(self, **params_values):
        ci.set_nuisance_shear_calib(
          M = [
            params_values.get(p, None) for p in [
              // change DES_ to the name of the survey associated w/ XXX)
              "DES_M"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_shear_photoz(
          bias = [
            params_values.get(p, None) for p in [
              // change DES_ to the name of the survey associated w/ XXX)
              "DES_DZ_S"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_ia(
          A1 = [
            params_values.get(p, None) for p in [
              // change DES_ to the name of the survey associated w/ XXX)
              "DES_A1_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          A2 = [
            params_values.get(p, None) for p in [
              // change DES_ to the name of the survey associated w/ XXX)
              "DES_A2_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          B_TA = [
            params_values.get(p, None) for p in [
              // change DES_ to the name of the survey associated w/ XXX)
              "DES_BTA_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
        )
     
    (...)
     
    def set_lens_related(self, **params_values):
        ci.set_nuisance_bias(
            B1 = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "DES_B1_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ], 
            B2 = [
                  params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "DES_B2_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ],
            B_MAG = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "DES_BMAG_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_nuisance_clustering_photoz(
            bias = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "DES_DZ_L"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_point_mass(
            PMV = [
                params_values.get(p, None) for p in [
                  // change DES_ to the name of the survey associated w/ XXX)
                  "DES_PM"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
     
    (...)
     
    def set_baryon_related(self, **params_values):
        // change DES_ to the name of the survey associated w/ XXX)
        self.baryon_pcs_qs[0] = params_values.get("DES_BARYON_Q1", None)
        // change DES_ to the name of the survey associated w/ XXX)
        self.baryon_pcs_qs[1] = params_values.get("DES_BARYON_Q2", None)
        // change DES_ to the name of the survey associated w/ XXX)
        self.baryon_pcs_qs[2] = params_values.get("DES_BARYON_Q3", None)
        // change DES_ to the name of the survey associated w/ XXX)
        self.baryon_pcs_qs[3] = params_values.get("DES_BARYON_Q4", None)

(**Warning**) If the project name `XXX` contains more than the experiment name (e.g., `XXX = LSST_Y1`), we suggest to replacing `DES_` with just the experiment name (e.g., `LSST_BARYON_Q1`, `LSST_PM` and `LSST_DZ_L`). The convention adopted must be followed when changing the files `params_des_cosmic_shear.yaml` and `params_des_3x2pt.yaml`. 

**Step 2:** Change the file `$ROOTDIR/projects/XXX/likelihood/des_3x2pt.py` following the instructions below
    
    // change des_y3 to XXX in the line below
    from cobaya.likelihoods.des_y3._cosmolike_prototype_base import _cosmolike_prototype_base
    // change cosmolike_des_y3_interface to cosmolike_XXX_interface in the line below
    import cosmolike_des_y3_interface as ci

(**Warning**) Similar changes must be made in the following files
    
    +-- cocoa_des_y3
    |    +-- likelihood
    |    |   +-- des_2x2pt.py
    |    |   +-- des_clustering.py
    |    |   +-- des_cosmic_shear.py
    |    |   +-- des_ggl.py
    |    |   +-- des_xi_ggl.py
    
**Step 3:** Change the file `$ROOTDIR/projects/XXX/likelihood/des_3x2pt.yaml` following the instructions below
   
    (...)
    // change DES_Y3.dataset to XXX.dataset in the line below (adopted convention: .dataset file name = project name all in CAPS)
    data_file: DES_Y3.dataset
   
    (...)
    // change params_des_3x2pt to params_XXX_3x2pt in the line below
    params: !defaults [params_des_3x2pt]

(**Warning**) Similar changes must be made in the following files
    
    +-- cocoa_des_y3
    |    +-- likelihood
    |    |   +-- des_2x2pt.yaml
    |    |   +-- des_clustering.yaml
    |    |   +-- des_cosmic_shear.yaml
    |    |   +-- des_ggl.yaml
    |    |   +-- des_xi_ggl.yaml

**Step 4:** Change the file `params_des_3x2pt.yaml` following the instructions below

Replace the `DES_` prefix to the name of the survey associated w/ XXX.
    
    DES_DZ_S1:
    	prior:
      	    dist: norm
      	    loc: 0.0
      	    scale: 0.005
    ref:
        dist: norm
        loc: 0.0
        scale: 0.005
    proposal: 0.005
    latex: \Delta z_\mathrm{s, DES}^1

(**Warning**) If the project name `XXX` contains more than the experiment name  (e.g., `XXX = LSST_Y1`), we suggest to replacing `DES_` with just the experiment name (e.g., `LSST_DZ_S1` and `\Delta z_\mathrm{s, LSST}^1`)

(**Warning**) Similar changes must be made in `params_des_cosmic_shear.yaml`

(**Warning**) Changes in either the number of lenses or source bins will require the introduction of new paremeters in 
`params_des_cosmic_shear.yaml` and `params_des_3x2pt.yaml`

### Changes in the data folder

**Step 1** Rename the `.dataset` dile (adopted convention: `.dataset` file name = project name capitalized)

     $(cocoa)(.local) mv $ROOTDIR/projects/XXX/data/DES_Y3.dataset $ROOTDIR/projects/XXX/data/XXX.dataset
     
**Step 2** Update `XXX.dataset` file with the names of the new data vector, covariance, n(z), binning, mask...

     data_file = LSST_Y1_nonlim
     cov_file = cov_lsst_y1
     mask_file = 3x2pt_baseline.mask
     nz_lens_file = lens_LSSTY1.nz
     nz_source_file = source_LSSTY1.nz
     lensing_overlap_cut = 0.0015
     lens_ntomo = 5
     source_ntomo = 4
     n_theta = 20
     IA_model = 6
     theta_min_arcmin = 2.5
     theta_max_arcmin = 250.
     #baryon_pca_file = pca.txt

### Changes in the chain folder

**Step 1:** Remove possible old chains associated with the old `des_y3` project 

     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.1.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.2.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.3.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.4.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.5.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.6.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.7.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.8.txt
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.py
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.yaml
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.input.yaml
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.updated.yaml 
     $(cocoa)(.local) rm $ROOTDIR/projects/$NEW_PROJECT/chains/*.pyc

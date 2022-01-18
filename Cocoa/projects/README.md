# Table of contents
1. [The Projects Folder](#appendix_projects_folder)
2. [Adapting DES-Y3 to a new project](#appendix_des_y3_new)

## The Projects Folder <a name="appendix_projects_folder"></a> 

The `projects` folder includes all the projects that are being developed by our group. Individual projects must be hosted on independent folders named `cocoa_XXX` where XXX is the project name. The majority of projects we are working on are not public (yet), and they are safeguarded on the private repositories listed on `project/clone_all.sh` (the backbone Cosmolike software, however, is publicly available at `external_modules/code`!). You can add your projects there, and the script `setup_cocoa_installation_packages` will try to clone all listed projects. Having inaccessible repositories listed at `project/clone_all.sh` will not cause any errors. 

The `cocoa_XXX` folder that host the `XXX` project needs to have the more or less the following structure (taken from our private DES-Y3 project)

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

(**warning**) The DES-Y3 project is not public yet, but the code will be release in the near future. 

Adapting the DES-Y3 folder to construct a new project involves many small changes and a few big, significant ones. 

The big changes are:

* Computation of the covariance matrix using either [CosmoCov](https://github.com/CosmoLike/CosmoCov) or [CosmoCovFourier](https://github.com/CosmoLike/CosmoCov_Fourier) (replacing `./projects/des_y3/data/cov_Y3.txt`)
* Simulation of new `n(z)` for lenses and sources (replacing `./projects/des_y3/data/nz_lens_Y3.txt` and `./projects/des_y3/data/nz_source_Y3.txt`)
* Changes to the Cosmolike C++ interface so the appropriate routines can be called from the Python likelihood
* Changes to the Cosmolike Python likelihood so the appropriate routines can be called from Cobaya
* Additional changes in the `data` file, including the `DES_Y3.dataset` file

Now we list the long list of small changes so the C - C++ - Python interface can work flawlessly. They are tedious, but straightforward

### Changes in the interface folder 

**Step 1** Choose a project name (e.g., project XXX), and copy the DES-Y3 project using the command below
    
    $ cp ./projects/des_y3/ ./projects/xxx

**Step 2** Change the file `./projects/XXX/interface/MakefileCosmolike` following the instructions below

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

**Step 3** Change the name of the File `./projects/XXX/interface/cosmolike_des_y3_interface.py` using the command below
    
    $ mv ./projects/XXX/interface/cosmolike_des_y3_interface.py ./projects/XXX/interface/cosmolike_XXX_interface.py

**Step 4** Changes in the newly created file `./projects/XXX/interface/cosmolike_XXX_interface.py` 

    def __bootstrap__():
        (...)
        // change cosmolike_des_y3_interface.so to cosmolike_XXX_interface.so in the line below 
        __file__ = pkg_resources.resource_filename(__name__,'cosmolike_des_y3_interface.so')
        
**Step 5** Change the file `./projects/XXX/interface/interface.cpp` following the instructions below
    
    (...)
    
    // change cosmolike_des_y3_interface to cosmolike_XXX_interface in the line below
    PYBIND11_MODULE(cosmolike_des_y3_interface, m)
    {
        // change the description below
        m.doc() = "CosmoLike Interface for DES-Y3 3x2 Module";
        
       (...)
    }
    
### Changes in the script folder

**Step 1** Change the name of the file `./projects/XXX/scripts/compile_des_y3` using the command below 
    
    $ mv ./projects/XXX/scripts/compile_des_y3 ./projects/XXX/scripts/compile_XXX
    
**Step 2** Change the name of the file `./projects/XXX/scripts/start_des_y3` using the command below 
    
    $ mv ./projects/XXX/scripts/start_des_y3 ./projects/XXX/scripts/start_XXX
    
**Step 3** Change the name of the file `./projects/XXX/scripts/stop_des_y3` using the command below 
    
    $ mv ./projects/XXX/scripts/stop_des_y3 ./projects/XXX/scripts/stop_des_y3

**Step 4** Change the file `./projects/XXX/scripts/compile_des_y3` following the instructions below

    (...)

    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    cd $ROOTDIR/projects/des_y3/interface
    
**Step 5** Change the file `./projects/XXX/scripts/start_des_y3` following the instructions below

    (...)

    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar LD_LIBRARY_PATH $ROOTDIR/projects/des_y3/interface
    
    // change $ROOTDIR/projects/des_y3/interface to $ROOTDIR/projects/XXX/interface in the line below 
    addvar PYTHONPATH $ROOTDIR/projects/des_y3/interface

### Changes in the likelihood folder

**Step 1** Change the file `./projects/XXX/likelihood/_cosmolike_prototype_base.py` following the instructions below

    (...) 
    
    // change cosmolike_des_y3_interface to cosmolike_XXX_interface in the line below 
    import cosmolike_des_y3_interface as ci
    
    (...)
     
    def set_source_related(self, **params_values):
        ci.set_nuisance_shear_calib(
          M = [
            params_values.get(p, None) for p in [
              // change DES_ to XXX_ in the line below
              "DES_M"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_shear_photoz(
          bias = [
            params_values.get(p, None) for p in [
              // change DES_ to XXX_ in the line below
              "DES_DZ_S"+str(i+1) for i in range(self.source_ntomo)
            ]
          ]
        )

        ci.set_nuisance_ia(
          A1 = [
            params_values.get(p, None) for p in [
              // change DES_ to XXX_ in the line below
              "DES_A1_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          A2 = [
            params_values.get(p, None) for p in [
              // change DES_ to XXX_ in the line below
              "DES_A2_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
          B_TA = [
            params_values.get(p, None) for p in [
              // change DES_ to XXX_ in the line below
              "DES_BTA_"+str(i+1) for i in range(self.source_ntomo)
            ]
          ],
        )
     
    (...)
     
    def set_lens_related(self, **params_values):
        ci.set_nuisance_bias(
            B1 = [
                params_values.get(p, None) for p in [
                  // change DES_ to XXX_ in the line below
                  "DES_B1_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ], 
	    B2 = [
                params_values.get(p, None) for p in [
                  // change DES_ to XXX_ in the line below
                  "DES_B2_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ],
	    B_MAG = [
                params_values.get(p, None) for p in [
                  // change DES_ to XXX_ in the line below
                  "DES_BMAG_"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_nuisance_clustering_photoz(
            bias = [
                params_values.get(p, None) for p in [
                  // change DES_ to XXX_ in the line below
                  "DES_DZ_L"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
        ci.set_point_mass(
            PMV = [
                params_values.get(p, None) for p in [
                  // change DES_ to XXX_ in the line below
                  "DES_PM"+str(i+1) for i in range(self.lens_ntomo)
                ]
            ]
        )
     
    (...)
     
    def set_baryon_related(self, **params_values):
        // change DES_ to XXX_ in the line below
        self.baryon_pcs_qs[0] = params_values.get("DES_BARYON_Q1", None)
        // change DES_ to XXX_ in the line below
        self.baryon_pcs_qs[1] = params_values.get("DES_BARYON_Q2", None)
        // change DES_ to XXX_ in the line below
        self.baryon_pcs_qs[2] = params_values.get("DES_BARYON_Q3", None)
        // change DES_ to XXX_ in the line below
        self.baryon_pcs_qs[3] = params_values.get("DES_BARYON_Q4", None)

**Step 2** Change the file `./projects/XXX/likelihood/des_3x2pt.py` following the instructions below
    
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
    
**Step 3** Change the file `./projects/XXX/likelihood/des_3x2pt.yaml` following the instructions below
   
    (...)
    // change DES_Y3.dataset to XXX.dataset in the line below
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

# Table of contents
1. [FAQ: How to download modern CMB data?](#new_planck_data)
2. [FAQ: FAQ: How to download new data from current experiments (using git)?](#new_likelihood_and_data)
3. [FAQ: FAQ: How to download new data from current experiments (using wget)?](#new_likelihood_and_data2)
 
Cocoa provides a list of shell scripts, located at `Cocoa/installation_scripts`, that manages the download and installation of Lipop (CMB), Camspec (CMB), SPT (CMB), Simons Observatory (CMB), H0licow (Strong Lensing), and other datasets. They all start with the prefix `unxv_`. 

## :interrobang: FAQ: How to download modern CMB data? <a name="new_planck_data"></a>

Cocoa can download [CamSpec](https://people.ast.cam.ac.uk/~stg20/camspec/index.html), [sroll2](https://web.fe.infn.it/~pagano/low_ell_datasets/sroll2/), [Hillipop](https://github.com/planck-npipe/hillipop.git), and [Lollipop](https://github.com/planck-npipe/lollipop.git) Planck-CMB likelihoods. Cocoa can also download data from some CMB high-resolution ground observatories. These datasets can be quite large, containing files that are several Gigabytes in size, so the shell script  `Cocoa/set_installation_options.sh` contains keys that allow users to skip their download, as shown below.

    [Adapted from Cocoa/set_installation_options.sh shell script] 

    # ------------------------------------------------------------------------------
    # The flags below allow users to skip downloading specific datasets ------------
    # ------------------------------------------------------------------------------
    
    (...)

    # export IGNORE_SETUP_SPT_CMB_DATA=1
    export IGNORE_SETUP_SIMONS_OBSERVATORY_CMB_DATA=1
    # export IGNORE_SETUP_PLANCK_CMB_DATA=1
    export IGNORE_SETUP_CAMSPEC_CMB_DATA=1
    export IGNORE_SETUP_LIPOP_CMB_DATA=1

Cocoa selects the URL to download the data (and its version) using the following keys.

    [Adapted from Cocoa/set_installation_options.sh shell script] 

    # ------------------------------------------------------------------------------
    # PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
    # ------------------------------------------------------------------------------

    (...)
    
    export LIPOP_DATA_URL="https://portal.nersc.gov/cfs/cmb/planck2020/likelihoods"
    export LIPOP_DATA_VERSION=4.2

    export SPT3G_DATA_URL='https://github.com/SouthPoleTelescope/spt3g_y1_dist.git'
    export SPT3G_DATA_GIT_COMMIT="66da8e9e2f325024566fe13245788bf8ede897bc"

    export ACT_DR6_DATA_URL="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"
    export ACT_DR6_DATA_FILE="ACT_dr6_likelihood_v1.2.tgz"

    export SO_DATA_URL="https://portal.nersc.gov/cfs/sobs/users/MFLike_data"
    # Cocoa can download multiple versions of the data (to reproduce existing work)
    # This is only possible because each version is saved in a separate folder
    export SO_DATA_VERSION="v0.7.1 v0.8"

## :interrobang: FAQ: How to download new data from experiments (using git)? <a name="new_likelihood_and_data"></a>

 The script `unxv_github_template.sh` provides a basic template for users to add a new dataset by cloning a git repository. In this script, the main lines that need to be modified are shown below.

    [Adapted from Cocoa/installation_scripts/unxv_git_template.sh shell script] 

    URL="${XXX_DATA_URL:-"https://github.com/XXX"}"
    
    # FOLDER = the directory name of the dataset - git clone $URL $FOLDER
    FOLDER="XXX"

    # Name to be printed on this shell script messages
    PRINTNAME=XXX
  
    (...) 

    # SET XXX_DATA_GIT_COMMIT IN CASE THE USER WANTS TO CHECKOUT A SPECIFIC COMMIT
    if [ -n "${XXX_DATA_GIT_COMMIT}" ]; then

      (...)
      
      ${GIT:?} checkout "${XXX_DATA_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
    fi
  
## :interrobang: FAQ: How to download new data from experiments (using wget)? <a name="new_likelihood_and_data2"></a>

# Table of contents
1. [FAQ: How to download modern CMB data?](#new_planck_data)
2. [FAQ: FAQ: How to download new data using git?](#new_likelihood_and_data)
3. [FAQ: FAQ: How to download new data using wget?](#new_likelihood_and_data2)
 
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

## :interrobang: FAQ: How to download new data using git? <a name="new_likelihood_and_data"></a>

 Suppose the user wants to download a dataset for which Cocoa does not have an already developed shell script at `Cocoa/installation_scripts`. In that case, the script `unxv_github_template.sh` provides a basic template for adding a new dataset by cloning a git repository. The main lines that need to be modified in this script are shown below.

    [Adapted from Cocoa/installation_scripts/unxv_git_template.sh shell script] 
    
    if [ -z "${IGNORE_SETUP_XXX_DATA}" ]; then              # Change the IGNORE_SETUP_XXX__DATA key name

      (...)
                                                            # URL = GitHub URL where data is located
      URL="${XXX_DATA_URL:-"https://github.com/XXX"}"       # Change the string associated with the URL key
    
                                                            # FOLDER = the dataset directory name
      FOLDER="XXX"                                          # Change the string associated with the FOLDER key

                                                            # PRINTNAME = Name to be printed on messages
      PRINTNAME=XXX                                         # Change the string associated with the PRINTNAME key
  
      (...) 
                                                            # XXX_DATA_GIT_COMMIT = commit hash
      if [ -n "${XXX_DATA_GIT_COMMIT}" ]; then              # Change the XXX_DATA_GIT_COMMIT key name

        (...)
      
        ${GIT:?} checkout "${XXX_DATA_GIT_COMMIT:?}" \
          >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
      fi
      
      (...)
    
    fi
  
## :interrobang: FAQ: How to download new data using wget? <a name="new_likelihood_and_data2"></a>

 Suppose the user wants to download a dataset for which Cocoa does not have an already developed shell script at `Cocoa/installation_scripts`. In that case, the script `unxv_wget_template.sh` provides a basic template for adding a new dataset by downloading files from an FTP server. The main lines that need to be modified in this script are shown below.

    [Adapted from Cocoa/installation_scripts/unxv_git_template.sh shell script] 
    
    if [ -z "${IGNORE_SETUP_XXX_DATA}" ]; then              # Change the IGNORE_SETUP_XXX__DATA key name

      (...)
                                                            # URL = FTP URL where data is located
      URL="${XXX_DATA_URL:-"https://website/XXX"}"          # Change the string associated with the URL key

                                                            # FOLDER = the directory name of the dataset
      FOLDER="XXX"                                          # Change the string associated with the FOLDER key

      declare -a FILE=( "filename1"                         # FILE = list the names of the files to be downloaded
                        "filename2"                         # Change the strings associated with the FILE list 
                        "filename3"
                      )

      declare -a EXT=( "tar.gz"                             # EXT = list the extension of the files to be downloaded
                       "tar.gz"                             # Change the strings associated with the EXT list  
                       "tar.gz"
                      )
                                                            # PRINTNAME = Name to be printed on messages
      PRINTNAME=XXX                                         # Change the string associated with the PRINTNAME key

      (...)

      # ------------------------------------------------------------------------------------------
      # Here, you may need to add an option below related to the file extension of your dataset
      # if file extension != "tar.gz" or "tar.xz"
      
      for (( i=0; i<${#FILE[@]}; i++ ));
      do
    
        (...)
    
        if [ "${EXT:?}" == "tar.gz" ]; then
          
          (...)
        
        elif [ "${EXT:?}" == "DATASET FILE EXTENSION" ]; then

          # Add code on how to decompress the file extension of your dataset
        
        else
        
          (...)
 
        fi
    
      done
    
    fi

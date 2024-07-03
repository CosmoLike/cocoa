# Table of contents
1. [FAQ: How to download modern Planck data](#new_planck_data)


### :interrobang: FAQ: How to download modern CMB data? <a name="new_planck_data"></a>

Cocoa can download [CamSpec](https://people.ast.cam.ac.uk/~stg20/camspec/index.html), [sroll2](https://web.fe.infn.it/~pagano/low_ell_datasets/sroll2/), [Hillipop](https://github.com/planck-npipe/hillipop.git), and [Lollipop](https://github.com/planck-npipe/lollipop.git) Planck-CMB likelihoods. CoCoa can also download data from some CMB high-resolution ground observatories. These datasets can be quite large, containing files that are several Gigabytes in size, so the shell script  `Cocoa/set_installation_options.sh` contains keys that allow users to skip their download, as shown below.

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

    

.. _WMAP:

Using WMAP9 likelihood
======================

This utility is only available when the optional python tools are installed using  waf.

If the sources for the ``wmap9`` likelihood are available, one can also use this likelihood code along with the WMAP9 likelihood data file whithin clik. To do so, one must also create likelihood files that refers to the WMAP9 dataset. The tool :program:`prepare_wmap` allows to prepare this files. It needs a parameter file describing the location of the WMAP data the description of the range of ells to use and a few flags. Example parameter files are present in the ``examples`` directory. See for ``wmap_full.par`` that defines likelihood using all of the WMAP9 data for all ells.

.. code-block:: python

    # prepare a file for calling the wmap full likelihood through clik

    res_object = wmap_9_full.clik

    ttmin = 2
    ttmax = 1200
    temin = 2
    temax = 1200
    use_gibbs = 0
    use_lowl_pol = 1

    wmap_data = /THE/PATH/TO/wmap_likelihood_v5
    cl_save = wmap_9_full.cltest

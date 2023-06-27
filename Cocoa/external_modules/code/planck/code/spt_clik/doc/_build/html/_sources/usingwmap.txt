.. _WMAP:

Using WMAP7 likelihood
======================

If the sources for the ``wmap7`` likelihood are available, one can also use this likelihood code along with the WMAP7 likelihood data file whithin clik. To do so, one must also create likelihood files that refers to the WMAP7 dataset. The tool :program:`prepare_wmap` allows to prepare this files. It needs a parameter file describing the location of the WMAP data the description of the range of ells to use and a few flags. Example parameter files are present in the ``examples`` directory. See for ``wmap_full.par`` that defines likelihood using all of the WMAP7 data for all ells.

.. code-block:: python

    # prepare a file for calling the wmap full likelihood through clik

    res_object = wmap_7_full.h5

    ttmin = 2
    ttmax = 1200
    temin = 2
    temax = 1200
    use_gibbs = 0
    use_lowl_pol = 1

    wmap_data = /THE/PATH/TO/likelihood_v4p1
    cl_save = wmap_7_full.cltest

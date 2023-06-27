Creating synthetic likelihoods
==============================

The tool :program:`synthetic_smica` allows to create synthetic likelihood files based on input power spectra and a description of the instrument. The likelihood approximation used is the so-called SMICA likelihood (an offset inverse wishart shape). The program expects a single command line argument, the path to a parameter file. The directory ``examples`` in the clik package contains a few example parameter files.

Here is an example with some explanation of the meaning of the different options.

.. code-block:: python

    # this is a likelihood for the HFI channels 143 and 217
    # TT TE EE
    
    #cl file
    # it can be either the output from CAMB (tot_cls) or the cls corresponding to the has_cl from l=0 
    cl = bestfit_lensedCls.dat
    
    #lmax for all the cls. The file can contain more mode, they will be discarded
    lmax = 1500

    # optional lmin. Do not use mode below lmin
    # lmin = 32

    # list of flags for the presence of each power spectra
    # order is TT EE BB TE TB EB
    has_cl =    1  1  0  1  0  0


    # optional list of mix values this is the gain of each detector. Better leave it to 1
    mixcol = 1 1 1 1
    # or file containing the same info
    # Acmb = 

    # optional file containing the binning matrix
    # bins = 
    # or size for each bin
    binsz = 10


    # number of Temperature channels
    nT = 2
    # number of Polar channels
    nP = 2

    # optional offset matrix file
    # Rq_0 = 

    # optional noise matrix file
    # nQ = 
    # or list of noise level for each channel (T then P)
    # noise is in microK^2
    noise = 0.0005 0.0001 0.001 0.0002 

    # list of full width half max for each channel (T then P)
    fwhm = 9.6 7 9.6 7

    # optional fsky
    fsky = .8
    # or weight for each bin (in a file)
    # wq = 

    # name of the resulting lkl file
    res_object = fake_smica_TE_32_1500_b10_100x143.h5


    # if meanfield is set to 1, no synthetic data is produced
    # meanfield = 1

    # optional seed for fake data
    # seed = 123456 
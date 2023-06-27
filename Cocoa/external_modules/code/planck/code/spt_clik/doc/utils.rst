Playing around with likelihood files
====================================

Computing a log likelihood from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The example codes, :program:`clik_example_C`, :program:`clik_example_f90` and :program:`clik_example_py` allow to compute a 
the log likelihoods for any numbers of files containing Cls andforeground parameters. 

:program:`clik_example_C` *usage:*

.. code-block:: none

    clik_example_C lkl_file.clik [clfile1 ...]

``lkl_file.clik`` is the likelihood file. The ``clfile1 ...`` files must be ascii and contains 
Cls from 0 to the lmax (included) of the likelihood file, followed by the nuisance parameter values in the order shown when 
using :program:`clik_print` or using of the the query function (for example, in c :cfunction:`clik_get_extra_parameter_names`). 

The program :program:`clik_example_py` is only available when the optional python tools are installed either by make or waf.


Printing info about a file
^^^^^^^^^^^^^^^^^^^^^^^^^^

This utility is only available when the optional python tools are installed either by make or waf.

.. program:: clik_print

The tool  :program:`clik_print` displays some information on the content of a likelihood files. The range of modes for each power spectrum, the list of extra parameters, and for each component of the full likelihood, the same info.

*Usage:*

.. code-block:: none

    clik_print somelikelihoodfile.clik

``somelikelihoodfile.clik`` is a likelihood file.

*Example output:*

.. code-block:: none

    $> clik_print ../release/clik_7.4/CAMspec_v6.2TN_2013_02_26.clik/
    ----
    clik version 5869
      CAMspec e61cec87-3a37-43ca-8ed1-edcfcaf5c00a
    Checking likelihood '../release/clik_7.4/CAMspec_v6.2TN_2013_02_26.clik/' on test data. got -3910.03 expected -3910.03 (diff -2.09184e-10)
    ----
    clik lkl file =  ../release/clik_7.4/CAMspec_v6.2TN_2013_02_26.clik/
      number of likelihoods = 1
      lmax ( TT = 2500 )
      number of varying extra parameters 15
        A_ps_100
        A_ps_143
        A_ps_217
        A_cib_143
        A_cib_217
        A_sz
        r_ps
        r_cib
        n_Dl_cib
        cal_100
        cal_143
        cal_217
        xi_sz_cib
        A_ksz
        Bm_1_1

      lkl_0
        lkl_type = CAMspec
        unit = 1
        TT = [50 , 2500]
        number of extra parameters = 15 ('A_ps_100', 'A_ps_143', 'A_ps_217', 'A_cib_143', 'A_cib_217', 'A_sz', 'r_ps', 'r_cib', 'n_Dl_cib', 'cal_100', 'cal_143', 'cal_217', 'xi_sz_cib', 'A_ksz', 'Bm_1_1')



Modifying the content of a likelihood file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This utility is only available when the optional python tools are installed either by make or waf.

The tools :program:`clik_join` and  :program:`clik_disjoin` allow to either join toghether one or more likelihood files in a single one, or cut a likelihood files into as many files as it has components.

.. program:: clik_join

:program:`clik_join` *usage:*

.. code-block:: none

    clik_join lkl_file_1.clik lkl_file_2.clik [lkl_file_3.clik ...] result_lkl_file.clik

``lkl_file_1.clik``, ``lkl_file_2.clik``... are input likelihood files. The resulting file ``result_lkl_file.clik`` defines a likelihood file so that the log likelihood a Cl (+extra parameters) is the sum of the log likelihood of each input files.

.. program:: clik_disjoin

:program:`clik_disjoin` *usage:*

.. code-block:: none

    clik_disjoin lkl_file.clik

The input file is ``lkl_file.clik`` is split in as many likelihood as it has component. Each likelihood is saved in its own file, named ``lkl_file.lkl_X.clik`` where ``X`` is a number between 0 and the number of components. 


Dealing with likelihood files with external data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This utility is only available when the optional python tools are installed either by make or waf.


This is only valid for likelihood files containing only one component and when this component is either a BOPIX or WMAP likelihood. In both cases, the likelihood relies on external data. This data is either included in the file (as a big tarfile) or install somewhere in the file system. the tools :program:`clik_extract_external` and :program:`clik_include_external` allows to go from one choice to the other. It is either, when distribution, to include the external data whithin the file, and more efficient to run with the external data installed somewhere in the file system.

.. program:: clik_extract_external

:program:`clik_extract_external` *usage:*

.. code-block:: none

    clik_extract_external parameterfile

*Example parameter file*

.. code-block:: none

    input_object = wmap_7_full.clik              # input likelihood file. Data is included
    install_path = /data/wmap_likelihood_data    # where to install the data
    res_object = wmap_7_full.external.clik       # output likelihood file. Data is no more included
    
.. program:: clik_include_external

:program:`clik_include_external` *usage:*

.. code-block:: none

    clik_include_external parameterfile

*Example parameter file*

.. code-block:: none

    input_object = wmap_7_full.external.clik   # input likelihood file. Data is installed somewhere
    res_object = wmap_7_full.clik              # output likelihood file. Data is included


Extracting the test Cl from a likelihood file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This utility is only available when the optional python tools are installed either by make or waf.


:program:`clik_get_selfcheck` *usage:*

.. code-block:: none

    clik_get_selfcheck lkl_file.clik clfile

``lkl_file.clik`` is the likelihood file. ``clfile`` is the cl+nuisance parameter array used to compute the selfchek displayed at each initialization of the likelihood. Same format as the one needed for :program:`clik_example_C`


Computing a slice through a log likelihood from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This utility is only available when the optional python tools are installed either by make or waf.

One can quickly compute conditionals through a likelihood along the direction of one of the nuisance parameter using :program:`clik_explore_1d`.

:program:`clik_explore_1d` *usage:*

.. code-block:: none

    clik_explore_1d parfile

``parfile`` is a parameter file similar to:

.. code-block:: python

    # slice 

    #lkl
    input_object = CAMspec_v6.2TN_2013_02_26.clik
    
    #data for the other dimensions. Same format as for clik_example_C. 
    initdata = bestfilcl.camspec
    
    #name of the varying parameter
    parameter = r_cib

    #begin and end values
    beg = -1
    end = 1.5

    #number of computations
    step = 300

    #ascii file that will hold the result as a 2d array, parameter value, lkl value
    res = myresult.txt

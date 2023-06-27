Playing around with likelihood files
====================================

Here is a list of utilities to manipulate the likelihood files.

Printing info about a file
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. program:: clik_print

The tool  :program:`clik_print` displays some information on the content of a likelihood files. The range of modes for each power spectrum, the list of extra parameters, and for each component of the full likelihood, the same info.

*Usage:*

.. code-block:: none

    clik_print somelikelihoodfile.clik

``somelikelihoodfile.clik`` is a likelihood file.

*Example output:*

.. code-block:: none

    Checking likelihood 'fake_smica_TE_2_1700_b10_143.clik' on test data. got -244.097 expected -244.097 (diff 0)
    clik lkl file =  fake_smica_TE_2_1700_b10_143.clik
      number of likelihoods = 1
      lmax ( TT = 1700 EE = 1700 TE = 1700 )
      number of extra parameters = 0 ()
      lkl_0
        lkl_type = smica
        unit = 1
        TT = [2 , 1700] EE = [2 , 1700] TE = [2 , 1700]
        nbins = 507
        number of extra parameters = 0 ()


Modifying the content of a likelihood file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


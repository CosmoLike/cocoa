Interfacing the library with python
===================================

API - CMB
---------

The module clik contains the wrapper to the clik c library.
It contains only one object called ``clik``, which is initialized with a string containing the path to a likelihood file.

.. code-block:: python

    import clik
    
    clikid = clik.clik("clikidfile")
    

The ``has_cl``, ``lmax`` and parameter names array (see :ref:`querying`) can be queried by simpliy reading the ``has_cl``, ``lmax`` and ``extra_parameter_names`` attributes of the object

.. code-block:: python

    has_cl = clikid.has_cl
    print has_cl
    
A log likelihood is computed by calling the object with a list-like object (``tuple``, ``list`` of ``numpy.ndarray`` objects) containing the vector of parameters as described in :ref:`querying`.

.. code-block:: python

    loglkl = clikid(cl_and_pars)

The file ``click_example_py.py`` gives a simple example of the use of the python API. It is compiled and installed as :program:`clik_example_py`.

API - lensing
-------------

Similarly a lensing likelihood can be initialized by 

.. code-block:: python

    import clik
    
    clikid = clik.clik_lensing("clikidfile")

The ``lmax`` and parameter names array (see :ref:`querying`) can be queried by simpliy reading the ``lmax`` and ``extra_parameter_names`` attributes of the object.

The log likelihood is computed by calling 

.. code-block:: python

    loglkl = clikid(cl_and_pars)

``cl_and_pars`` must be a (lensing_lmax+1)+number_of_extra_parameters elements array. The first lmax+1 elements must be the clpp, the next the cltt. Optionnaly if providing only lmax+1+number_of_extra_parameters the likelihood will be computed using the fiducial cltt spectrum.

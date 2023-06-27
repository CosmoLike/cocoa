Interfacing the library with a f90 executable
===========================================

The following gives a description of the f90 API of the library, and how to correctly compile and link against it.

Compiling and linking
---------------------

The program `clik_f90-config` (installed in PREFIX/bin) spits out on the standard output the barbaric option and link line to give to your prefered c compiler when compiling and linking against the clik lib.

The file ``click_example_f90.f90`` gives a simple example of the use of the f90 API. It is compiled and installed as :program:`clik_example_f90`.


API - CMB likelihood
--------------------
All codes calling clik functions must

.. code-block:: fortran

    use clik
    
The library can initialize more than one likelihood. Likelihood are represented by a variable (in the following, named ``clikid``) of type ``type(clik_object)``.

Initialization
^^^^^^^^^^^^^^

The library must be initialized by calling 

.. c:function:: subroutine clik_init(clikid,hdffilepath)

  The subroutine sets the argument ``clikid``, which is of type ``type(clik_object)`` to a handle on an object containing the definition of the likelihood. It expects two arguments, ``hdffilepath`` a string containing the path to a likelihood file. In case of error, the library will only print out a message and force the calling program to exit.
  
Querying the likelihood object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: subroutine clik_get_has_cl(clikid,has_cl)

  This function fills the ``integer(kind=4), dimension(6):: has_cl`` array with flags describing which power spectra are needed by the likelihood compute function (see :ref:`querying`). The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. In case of error the program exit with an explaining message.


.. c:function:: subroutine clik_get_lmax(clikid,lmax)

  This function fills the array ``integer(kind=4), dimension(6):: lmax`` with the lmax value for each power spectra needed by the likelihood compute function (see :ref:`querying`). The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. In case of error the program exit with an explaining message in case of an error.

.. c:function:: subroutine clik_get_extra_parameter_names(clikid,names,numnames)

  This function sets ``integer::numnames`` to the number of nuisance parameters needed by the likelihood compute function (see :ref:`querying`) and fills with their names the array ``character(len=256), dimension(numnames)::names``. This array is allocated by the function and MUST be deallocated by the caller after use. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. In case of error the program exit with an explaining message.
    
Computing the log likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: real(kind=8) function clik_compute(clikid,cl_and_pars)
  
  This function returns the value of the log likelihood for the parameter vector ``cl_and_pars``. The content of this vector is described in :ref:`querying`. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`.  In case of error the program exit with an explaining message. This function can be called as many time as the user wants.
  
  
Cleanup
^^^^^^^

When a likelihood object is no more needed (i.e. when no more computation will be needed in the program), the memory it uses can be cleaned up calling

  .. c:function:: subroutine clik_cleanup(clikid)
  
  The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. 
 

API - lensing likelihood
------------------------
All codes calling clik functions must

.. code-block:: fortran

    use clik
    
The library can initialize more than one likelihood. Likelihood are represented by a variable (in the following, named ``clikid``) of type ``type(clik_object)``.

Testing whether a file contains a lensing likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can test whether a file contains a lensing likelihood by calling

.. c:function:: subroutine clik_try_lensing(hdffilepath, is_lensing);
	
	On return the logical argument ``is_lensing`` is set to true or false depending whether the ``hdffilepath`` argument points toward a lensing likelihood file. In case of error, the library will only print out a message and force the calling program to exit.


Initialization
^^^^^^^^^^^^^^

The library must be initialized by calling 

.. c:function:: subroutine clik_lensing_init(clikid,hdffilepath)

  The subroutine sets the argument ``clikid``, which is of type ``type(clik_object)`` to a handle on an object containing the definition of the likelihood. It expects two arguments, ``hdffilepath`` a string containing the path to a likelihood file. In case of error, the library will only print out a message and force the calling program to exit.
  
Querying the lensing likelihood object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: subroutine clik_get_lmax(clikid,lmax)

  On return the integer argument ``lmax`` take as value the lmax of both the clpp and cltt. The first argument of the function must be the return value from a previous call to :c:func:`clik_lensing_init`. In case of error the program exit with an explaining message in case of an error.

.. c:function:: subroutine clik_lensing_get_extra_parameter_names(clikid,names,numnames)

  This function sets ``integer::numnames`` to the number of nuisance parameters needed by the likelihood compute function (see :ref:`querying`) and fills with their names the array ``character(len=256), dimension(numnames)::names``. This array is allocated by the function and MUST be deallocated by the caller after use. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. In case of error the program exit with an explaining message.
    
Computing the log likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: real(kind=8) function clik_lensing_compute(clikid,cl_and_pars)
  
  This function returns the value of the log likelihood for the parameter vector ``cl_and_pars``. This vector must have 2*(lmax_lensing+1) + number_of_lensing_extra_parameters elements. They are first the lensing_lmax+1 values of clpp, then the lensing_lmax+1 values of the cltt, the the extra parameter values. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`.  In case of error the program exit with an explaining message. This function can be called as many time as the user wants.
  
  
Cleanup
^^^^^^^

When a likelihood object is no more needed (i.e. when no more computation will be needed in the program), the memory it uses can be cleaned up calling

  .. c:function:: subroutine clik_lensing_cleanup(clikid)
  
  The first argument of the function must be the return value from a previous call to :c:func:`clik_lensing_init`. 
 




    

 



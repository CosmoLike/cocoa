Interfacing the library with a c executable
===========================================

The following gives a description of the c API of the library, and how to correctly compile and link against it.

Compiling and linking
---------------------

The program `clik-config` (installed in PREFIX/bin) spits out on the standard output the barbaric option and link line to give to your prefered c compiler when compiling and linking against the clik lib.

The file ``click_example_c.c`` gives a simple example of the use of the c API. It is compiled and installed as :program:`clik_example_C`.

API - CMB likelihood
--------------------
All codes calling clik functions must

.. code-block:: c

    include "clik.h"
    

The library can initialize more than one likelihood. Likelihood are represented by a variable (in the following, named ``clikid``) of type ``clik_object*``.


Initialization
^^^^^^^^^^^^^^

The library must be initialized by calling 

.. c:function:: clik_object* clik_init(char* hdffilepath, error **err);

  The function returns a pointer on an object containing the definition of the likelihood. It expects two arguments, ``hdffilepath`` a string containing the path to a likelihood file, and ``err`` a c structure allowing error tracking. The error tracking system is provided by pmclib, please refer to its doc it for more info. If you don't which to use the error tracking system, set this argument to ``NULL``. In this case, the library will only print out a message and force the calling program to exit in case of an error.
  
Querying the likelihood object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void clik_get_has_cl(clik_object *clikid, int has_cl[6],error **err);

  This function fills the array ``has_cl`` with flags describing which power spectra are needed by the likelihood compute function (see :ref:`querying`). The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. The last argument allows error tracking. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error.


.. c:function:: void clik_get_lmax(clik_object *clikid, int lmax[6],error **err);

  This function fills the array ``lmax`` with the lmax value for each power spectra needed by the likelihood compute function (see :ref:`querying`). The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. The last argument allow to track errors. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error.

.. c:function:: int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error **err);

  This function returns the number of nuisance parameters needed by the likelihood compute function (see :ref:`querying`) and fills with their names the array ``*names``. This array is an array of parname, who are ``char[_pn_size]``. It is allocated by the function and MUST be deallocated by the caller after use. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. The last argument allow to track errors. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error.
    
Computing the log likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: double clik_compute(clik_object *clikid, double *cl_and_pars,error **err);
  
  This function returns the value of the log likelihood for the parameter vector ``cl_and_pars``. The content of this vector is desribed in :ref:`querying`. The first argument of the function must be the return value from a previous call to :c:func:`clik_init`. The last argument allow to track errors. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error. This function can be called as many time as the user wants.
  
  
Cleanup
^^^^^^^

When a likelihood object is no more needed (i.e. when no more computation will be needed in the program), the memory it uses can be cleaned up calling

  .. c:function:: void clik_cleanup(clik_object** pclikid);
  
  The first argument of the function must be the pointer on a variable containg the return value from a previous call to :c:func:`clik_init`. Upon return, the content of this variable will be changed to ``NULL``.
  
 
API - lensing likelihood
------------------------
All codes calling clik functions must

.. code-block:: c

    include "clik.h"
    

The library can initialize more than one likelihood. Likelihood are represented by a variable (in the following, named ``clikid``) of type ``clik_lensing_object*``.

Testing whether a file contains a lensing likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can test whether a file contains a lensing likelihood by calling

.. c:function:: int clik_try_lensing(char* hdffilepath, error **err);
	
	Return 1 or 0 depending if the file ``hdffilepath`` constains a lensing likelihood. If the file does not exist or cannot be read, an error us raised the usual way.


Initialization
^^^^^^^^^^^^^^

The lensing likelihood must be initialized by calling 

.. c:function:: clik_lensing_object* clik_lensing_init(char* hdffilepath, error **err);

  The function returns a pointer on an object containing the definition of the likelihood. It expects two arguments, ``hdffilepath`` a string containing the path to a lensing likelihood file, and ``err`` a c structure allowing error tracking. The error tracking system is provided by pmclib, please refer to its doc it for more info. If you don't which to use the error tracking system, set this argument to ``NULL``. In this case, the library will only print out a message and force the calling program to exit in case of an error.
  
Querying the lensing likelihood object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int clik_lensing_get_lmax(clik_lensing_object *clikid, error **err);

  This function returns the lmax value for both clpp and cltt. 

.. c:function:: int clik_get_lensing_extra_parameter_names(clik_lensing_object* clikid, parname **names, error **err);

  This function returns the number of nuisance parameters needed by the lensing likelihood compute function and fills with their names the array ``*names``. This array is an array of parname, who are ``char[_pn_size]``. It is allocated by the function and MUST be deallocated by the caller after use. The first argument of the function must be the return value from a previous call to :c:func:`clik_lensing_init`. The last argument allow to track errors. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error.
    
Computing the log likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: double clik_lensing_compute(clik_lensing_object *clikid, double *cl_and_pars,error **err);
  
  This function returns the value of the log likelihood for the parameter vector ``cl_and_pars``. This vector must have 2*(lmax_lensing+1) + number_of_lensing_extra_parameters elements. They are first the lensing_lmax+1 values of clpp, then the lensing_lmax+1 values of the cltt, the the extra parameter values.  The first argument of the function must be the return value from a previous call to :c:func:`clik_lensing_init`. The last argument allow to track errors. It can be left to ``NULL``, in which case no error tracking is performed and the program exit with an explaining message in case of an error. This function can be called as many time as the user wants.
  
  
Cleanup
^^^^^^^

When a lensing likelihood object is no more needed (i.e. when no more computation will be needed in the program), the memory it uses can be cleaned up calling

  .. c:function:: void clik_lensing_cleanup(clik_lensing_object** pclikid);
  
  The first argument of the function must be the pointer on a variable containg the return value from a previous call to :c:func:`clik_lensing_init`. Upon return, the content of this variable will be changed to ``NULL``.
  
 



    

 



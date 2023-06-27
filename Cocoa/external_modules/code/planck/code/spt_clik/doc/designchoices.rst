Design choices
==============

A likelihood is entirely defined by a likelihood file. 
------------------------------------------------------

The likelihood files (in fact directories, containing data as fits files and metadata as ascii files) contain all the info needed to define the likelihood. This encompass both the data needed to compute the likelihood,but also parameters describing the type of mathematical approximation used to compute this particular likelihood, and parameters describing the expected input of the likelihood. 

A likelihood file can also be the combinaison of various underlying likelihoods, the computation from the top likelihood will be the product of each of the underlying one. A tool is distributed to allow joining and disjoining likelihoods.
This means that a likelihood file containing the some of different likelihood approximation will have to contain the data for each of those approximation. This will translate into possibly huge file. And variations of those likelihoods will contains yet more copy of this huge data.
To solve partially this problem, it is possible for the likelihoods that need a lot of data to optionally refer to external files. This will be the case in particular for the WMAP and lowlike likelihoods.
In those cases, the data can either be included in the file (as described above) or installed in some directory, in which case the likelihood file simply refers to the path of this directory.The latter case improve the efficiency of the initialization of clik when using this kind of likelihood files.

In order to use the library, the first step is thus to initialize it with such a file. Functions or subroutine to perform this initialization are available in each language. Another function is provided to cleanup the memory at the end of the use of a given likelihood. Some likelihoods can be initialized within the same session, allowing to perform comparison between different likelihood approximation whithin the same run.


The library computes an approximation of the log likelihood
-----------------------------------------------------------

And nothing else. The library does not compute minus the log likelihood or a chi2 like value. Just the log likelihood.
A function is provided in each language to compute the log likelihood, given a set of parameters.

.. _querying:

The input of the compute function are multipoles of the power spectra and nuisance parameters.
----------------------------------------------------------------------------------------------

The ``compute`` function expect one single vector of double. This vector must contains power spectra one after the other, starting al l=0, in microK^2, and then the nuisance parameters. 
Since the likelihood is defined by the likelihood file the exact range of power spectra needed by this function can vary from likelihood approximation to likelihood approximations. The same holds for the nuisance parameters.
Functions are provided to query a likelihood file and obtain this info. More precisely, three such function are available.
 
    * clik_get_has_cl: retrieve an array of 6 flags, describing the power spectra needed (1 if needed 0 otherwise). The order is TT, EE, BB, TE, TB, EB. Thus if this function answers (1,1,0,1,0,0) it means that the input vector of the compute function must contain, in that order, the power spectra for TT, EE an TE.

    * clik_get_lmax: retrieve an array of 6 integer giving the lmax value for each power spectra. Order is same as above. -1 means that the spectra is not needed. Thus if this function answer (2000,1000,-1,2000,-1,-1). This means that the vector must contain, in that order, the first 2001 multipoles of the TT power spectra (from 0 to 2000 included) followed by the first 1001 multipoles of the EE power spectra (0 to 1000 included) and next by the first 2001 multipoles of the TE power spectra (0 to 2000 included). Thus the first 5003 elements of the parameter vectors are the values of different power spectra. Note that this is also the sum of the result array of clik_get_lmax plus 6. Isn't this fantastic ?

    * clik_get_parameter_names: returns the number of nuisance parameters and fills an array of string giving their names. For example, if this function returns 2, and ('sigma8', 'fwhm_error') it means that the last two elements of the parameter vector must be the value of sigma8 and fwhm_error, whatever those parameters mean.
    
To sumarize, the input vector of the compute function must be an array of N = Ncl + Nnuis doubles. Ncl being the sum + 6 of the return array of clik_get_lmax, and Nnuis is the return of clik_get_parameter_names. The power spectra must be the Ncl first elements of that array. They start at C0 and en up at some Clmax[i] (included), the ith elements of the return of clik_get_lmax. The ordering of the power spectra is always TT EE BB TE TB EB. The Cls are in microK^2. The nuisnce parameters are the Nnuis last element of the parameter vector. Their names are given by the return array of the function clik_get_parameter_names.

Pitfalls
--------

The function computes the log likelihood.

The Cls must be given in that order TT, EE, BB, TE, TB, EB. And, yes, the library expects C0 and C1 for each power spectra. 

The library expect power spectra and not l(l+1)Cl/2pi or other combination.

The library really wants power spectra in microK^2.


Lensing likelihood
------------------

The behaviour of the lensing likelihood is similar to the one of the CMB likelihoods except that it does not have nuisance parameters and expected the phiphi Cl as well as the TT Cl.

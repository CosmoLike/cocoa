# This file is part of EuclidEmulator2
# Copyright (c) 2018-2021 Mischa Knabenhans, Pedro Carrilho
#
# EuclidEmulator2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# EuclidEmulator2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

try:
    from classy import Class as _Class
except ImportError:
    print("\nClassy could not be found in your system.")
    print("Here are some suggestions:\n")
    print("\t -Download the Class from class-code.net and install it")
    print("\t  together with its wrapper classy (type 'make' instead of")
    print("\t  'make class'")
    print("\t -If you know that Class is installed on your system")
    print("\t  and yet classy could not be installed, try re-compiling")
    print("\t  Class with just ''make'' instead of ''make class''")
    print("NOTICE: Even without classy you can still use EuclidEmulator2")
    print("        to emulate boost factors. You won't be able to compute")
    print("        full power spectra, though.")


import cython
cimport cython
from libc.stdlib cimport free
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref

import sys as _sys
import numpy as np
from scipy.interpolate import CubicSpline as _CubicSpline
import warnings as _warnings

# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

# Import the C++ classes

cdef extern from "cosmo.h":
    cdef cppclass Cosmology:

        Cosmology(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s) except +



cdef extern from "emulator.h":
    cdef cppclass EuclidEmulator:

        double kvec[613];
        double Bvec[101][613];

        EuclidEmulator() except +

        void compute_nlc(Cosmology * csm, vector[double] redshift, int n_redshift);
        void write_nlc2file(const string &filename, vector[double] zvec, int n_redshift);


#Create new python classes for wrapping the c++ classes
cdef class PyCosmology:

     cdef Cosmology*cosm


     def __cinit__(self, double Omega_b , double Omega_m , double Sum_m_nu , double n_s , double h , double w_0 , double w_a , double A_s ):
        """Cython signature: void Cosmology(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s)"""
        assert isinstance(Omega_b, float), 'arg Omega_b wrong type'
        assert isinstance(Omega_m, float), 'arg Omega_m wrong type'
        assert isinstance(Sum_m_nu, float), 'arg Sum_m_nu wrong type'
        assert isinstance(n_s, float), 'arg n_s wrong type'
        assert isinstance(h, float), 'arg h wrong type'
        assert isinstance(w_0, float), 'arg w_0 wrong type'
        assert isinstance(w_a, float), 'arg w_a wrong type'
        assert isinstance(A_s, float), 'arg A_s wrong type'

        self.cosm = new Cosmology((<double>Omega_b), (<double>Omega_m), (<double>Sum_m_nu), (<double>n_s), (<double>h), (<double>w_0), (<double>w_a), (<double>A_s))

     def __dealloc__(self):
        del self.cosm


cdef class PyEuclidEmulator:

    cdef EuclidEmulator*ee2

    def __cinit__(self):
        self.ee2 = new EuclidEmulator()

    def compute_nlc(self,PyCosmology csm, redshift, n_redshift):
        self.ee2.compute_nlc((<Cosmology *> csm.cosm), redshift, n_redshift)

    def write_nlc2file(self,filename, zvec, n_redshift):
        self.ee2.write_nlc2file(<string>filename, zvec, n_redshift)

    def __dealloc__(self):
        del self.ee2

    # Attribute access
    @property
    def kvec(self):
        return self.ee2.kvec
    @kvec.setter
    def kvec(self, kvec):
        self.ee2.kvec = kvec
    @kvec.deleter
    def kvec(self):
        del self.kvec

    @property
    def Bvec(self):
        return self.ee2.Bvec
    @Bvec.setter
    def Bvec(self, Bvec):
        self.ee2.Bvec = Bvec
    @Bvec.deleter
    def Bvec(self):
        del self.Bvec



    #Extra functions to manipulate the results of the c++ computations

    ######################################################
    #################### Check input #####################
    ######################################################

    def check_param_range(self,par_dict): #, csm_index=0): #Only one cosmology for now
        """
        Checks if all parameters in the cosmology dictionary 'par_dict'
        passed to this function obey the limits set by the emulator.
        """

        om_b_range = [0.04, 0.06]
        om_m_range = [0.24, 0.40]
        m_nu_range = [0.00, 0.15]
        n_s_range  = [0.92, 1.00]
        h_range    = [0.61, 0.73]
        w_0_range  = [-1.3, -0.7]
        w_a_range  = [-0.7,  0.5]
        A_s_range  = [1.7e-9, 2.5e-9]

        om_b_not_in_range = om_b_range[0] > par_dict['Omega_b'] or\
                            om_b_range[1] < par_dict['Omega_b']

        om_m_not_in_range = om_m_range[0] > par_dict['Omega_m'] or\
                            om_m_range[1] < par_dict['Omega_m']

        m_nu_not_in_range = m_nu_range[0] > par_dict['m_ncdm'] or\
                            m_nu_range[1] < par_dict['m_ncdm']

        n_s_not_in_range = n_s_range[0] > par_dict['n_s'] or\
                           n_s_range[1] < par_dict['n_s']

        h_not_in_range = h_range[0] > par_dict['h'] or\
                         h_range[1] < par_dict['h']

        w_0_not_in_range = w_0_range[0] > par_dict['w0_fld'] or\
                           w_0_range[1] < par_dict['w0_fld']

        w_a_not_in_range = w_a_range[0] > par_dict['wa_fld'] or\
                           w_a_range[1] < par_dict['wa_fld']

        A_s_not_in_range = A_s_range[0] > par_dict['A_s'] or\
                           A_s_range[1] < par_dict['A_s']

        if om_b_not_in_range:
            raise ValueError("Parameter range violation: \nOmega_b is set to %f, but should be in the interval [0.04, 0.06]."
                             %(par_dict['Omega_b'] ))

        if om_m_not_in_range:
            raise ValueError("Parameter range violation: \nOmega_m is set to %f, but should be in the interval [0.24, 0.40]."
                             %(par_dict['Omega_m']))

        if m_nu_not_in_range:
            raise ValueError("Parameter range violation: \nm_ncdm is set to %f, but should be in the interval [0.00, 0.15]."
                             %(par_dict['m_ncdm']))

        if n_s_not_in_range:
            raise ValueError("Parameter range violation: \nn_s is set to %f, but should be in the interval [0.92, 1.00]."
                             %(par_dict['n_s']))

        if h_not_in_range:
            raise ValueError("Parameter range violation: \nh is set to %f, but should be in the interval [0.61, 0.73]."
                             %( par_dict['h']))

        if w_0_not_in_range:
            raise ValueError("Parameter range violation: \nw_0 is set to %f, but should be in the interval [-1.3, -0.7]."
                             %( par_dict['w0_fld']))

        if w_a_not_in_range:
            raise ValueError("Parameter range violation: \nw_a is set to %f, but should be in the interval [-0.7,  0.5]."
                             %( par_dict['wa_fld']))

        if A_s_not_in_range:
            raise ValueError("Parameter range violation: \nA_s is set to %e, but should be in the interval [1.7e-9, 2.5e-9]."
                             %( par_dict['A_s']))


    def convert_to_emu(self,class_pars_dict):
        """
        Signature:    convert_to_emu(class_pars_dict)

        Description:  Converts the set of parameters accepted by several
                      other codes into a set of parameters accepted by
                      EuclidEmulator2. Also checks the params exist.

        Input type:   python dictionary

        Output type:  python dictionary

        """
        if not isinstance(class_pars_dict, dict):
            raise TypeError("The cosmological parameters must be passed as a python dictionary.")

        if 'h' in class_pars_dict:
            h = class_pars_dict['h']
        elif 'hubble' in class_pars_dict:
            h = class_pars_dict['hubble']
        elif 'H0' in class_pars_dict:
            h = class_pars_dict['H0']/100.
        else:
            raise KeyError("Missing parameter h. Can't proceed.")

        if 'Omega_b' in class_pars_dict:
            Om_b = class_pars_dict['Omega_b']
        elif 'Omb' in class_pars_dict:
            Om_b = class_pars_dict['Omb']
        elif 'Omega_baryon' in class_pars_dict:
            Om_b = class_pars_dict['Omega_baryon']
        elif 'omega_b' in class_pars_dict:
            Om_b = class_pars_dict['omega_b']/h**2
        elif 'om_b' in class_pars_dict:
            Om_b = class_pars_dict['om_b']/h**2
        elif 'ombh2' in class_pars_dict:
            Om_b = class_pars_dict['ombh2']/h**2
        else:
            raise KeyError("Missing parameter Omega_b. Can't proceed.")

        # Currently only allowing this way of passing the neutrino mass
        # and only allowing one value to be passed.
        if 'm_ncdm' in class_pars_dict:
            m_ncdm = class_pars_dict['m_ncdm']
        elif 'm_nu' in class_pars_dict:
            m_ncdm = class_pars_dict['m_nu']
        elif 'mnu' in class_pars_dict:
            m_ncdm = class_pars_dict['mnu']
        elif 'neutrino_mass' in class_pars_dict:
            m_ncdm = class_pars_dict['neutrino_mass']
        else:
            print("Missing parameter m_nu. Will set to 0.")
            m_ncdm=0.0

        # Should give either Omega_m or Omega_cdm
        Om_cdm=0
        if 'Omega_m' in class_pars_dict:
            Om_m = class_pars_dict['Omega_m']
        elif 'Omm' in class_pars_dict:
            Om_m = class_pars_dict['Omm']
        elif 'Omega_matter' in class_pars_dict:
            Om_m = class_pars_dict['Omega_matter']
        elif 'omega_m' in class_pars_dict:
            Om_m = class_pars_dict['omega_m']/h**2
        elif 'om_m' in class_pars_dict:
            Om_m = class_pars_dict['om_m']/h**2
        elif 'ommh2' in class_pars_dict:
            Om_m = class_pars_dict['ommh2']/h**2
        elif 'Omega_cdm' in class_pars_dict:
            Om_cdm = class_pars_dict['Omega_cdm']
            Om_m = Om_b + Om_cdm
        elif 'Omc' in class_pars_dict:
            Om_cdm = class_pars_dict['Omc']
            Om_m = Om_b + Om_cdm
        elif 'omega_cdm' in class_pars_dict:
            Om_cdm = class_pars_dict['omega_cdm']/h**2
            Om_m = Om_b + Om_cdm
        elif 'omch2' in class_pars_dict:
            Om_cdm = class_pars_dict['omch2']/h**2
            Om_m = Om_b + Om_cdm
        else:
            raise KeyError("Missing parameter Omega_m or Omega_cdm. Can't proceed.")

        if 'n_s' in class_pars_dict:
            n_s = class_pars_dict['n_s']
        elif 'ns' in class_pars_dict:
            n_s = class_pars_dict['ns']
        else:
            raise KeyError("Missing parameter n_s. Can't proceed.")

        if 'A_s' in class_pars_dict:
            A_s = class_pars_dict['A_s']
        elif 'As' in class_pars_dict:
            A_s = class_pars_dict['As']
        elif 'ln10^{10}A_s' in class_pars_dict:
            A_s = np.exp(class_pars_dict['ln10^{10}A_s'])*1.0e-10
        else:
            raise KeyError("Missing parameter A_s or ln10^{10}A_s. Can't proceed.")

        # Using default values for DE params for LCDM case.
        if 'w0_fld' in class_pars_dict:
            w0_fld = class_pars_dict['w0_fld']
        elif 'w0' in class_pars_dict:
            w0_fld = class_pars_dict['w0']
        elif 'w_0' in class_pars_dict:
            w0_fld = class_pars_dict['w_0']
        elif 'w' in class_pars_dict:
            w0_fld = class_pars_dict['w']
        else:
            print("Missing parameter w0. Will set to -1.")
            w0_fld=-1.0

        if 'wa_fld' in class_pars_dict:
            wa_fld = class_pars_dict['wa_fld']
        elif 'wa' in class_pars_dict:
            wa_fld = class_pars_dict['wa']
        elif 'w_a' in class_pars_dict:
            wa_fld = class_pars_dict['w_a']
        else:
            print("Missing parameter wa. Will set to 0.")
            wa_fld=0.0

        # If Omega_cdm is given instead of Omega_m, add neutrino contribution using
        # approximation.
        if not(Om_cdm==0):
            Om_m=Om_m+m_ncdm/93.14/h**2

        emu_pars_dict = {'Omega_b': Om_b,
                         'Omega_m': Om_m,
                         'm_ncdm': m_ncdm,
                         'n_s': n_s,
                         'h': h,
                         'w0_fld': w0_fld,
                         'wa_fld': wa_fld,
                         'A_s': A_s,
                         'Omega_cdm': Om_cdm}

        return emu_pars_dict



    def get_boost(self,cosmo_par_in,redshifts,custom_kvec=None):

        if isinstance(redshifts, (int, float)):
            redshifts = np.asarray([redshifts])
        else:
            redshifts = np.asarray(redshifts)

        for z in redshifts:
            assert z <= 10.0 and z>=0.0, "EuclidEmulator2 allows only redshifts in the interval [0.0, 10.0]"

        #Check if all variables are passed and convert to emu dict
        cosmo_par = self.convert_to_emu(cosmo_par_in)
        #Check if all parameters are in range
        self.check_param_range(cosmo_par)

        cosmo=PyCosmology(cosmo_par['Omega_b'],
                          cosmo_par['Omega_m'],
                          cosmo_par['m_ncdm'],
                          cosmo_par['n_s'],
                          cosmo_par['h'],
                          cosmo_par['w0_fld'],
                          cosmo_par['wa_fld'],
                          cosmo_par['A_s'])

        len_redshifts = len(redshifts)
        self.compute_nlc(cosmo,redshifts,len_redshifts)

        k=np.asarray(self.kvec)
        logboost=np.reshape(self.Bvec[0:len_redshifts],(len_redshifts,len(k)))

        #Extrapolate for custom k-range
        kvals = k
        k_shape = kvals.shape

        do_extrapolate_above = False
        do_extrapolate_below = False
        if not(custom_kvec is None):

            if isinstance(custom_kvec, (int, float)):
                custom_kvec = np.asarray([custom_kvec])
            else:
                custom_kvec = np.asarray(custom_kvec)

            upper_mask = custom_kvec < max(kvals)
            lower_mask = custom_kvec > min(kvals)
            mask = [u and l for (u,l) in zip(lower_mask, upper_mask)]
            custom_k_within_range = custom_kvec[mask]
            custom_k_below = custom_kvec[[not(l) for l in lower_mask]]
            custom_k_above = custom_kvec[[not(u) for u in upper_mask]]

            if any(custom_kvec > max(kvals)):
                wrn_message = ("Warning:\nEuclidEmulator2 emulates the non-linear correction in \n"
                               "the interval [8.73e-3 h/Mpc, 9.41h/Mpc]. You are \n"
                               "requesting k modes beyond k_max = 9.41h/Mpc. \n"
                               "Higher k modes constantly extrapolated.")

                print(wrn_message)
                do_extrapolate_above = True

            if any(custom_kvec < min(kvals)):
                wrn_message = ("Warning:\nEuclidEmulator2 emulates the non-linear correction in \n"
                               "the interval [8.73e-3 h/Mpc, 9.41h/Mpc]. You are \n"
                               "requesting k modes below k_min = 8.73e-3 h/Mpc. \n"
                               "Lower k modes constantly extrapolated.")

                print(wrn_message)
                do_extrapolate_below = True

        len_kvals = len(kvals)


        bvals = {}
        for i in range(len_redshifts):
            tmp = logboost[i]
            if not(custom_kvec is None):
                bvals[i] = 10.0**_CubicSpline(np.log10(kvals),
                                              tmp.reshape(k_shape)
                                              )(np.log10(custom_k_within_range))

                #Extrapolate if necessary
                if do_extrapolate_below:
                    # below the k_min of EuclidEmulator2, we are in the linear regime where
                    # the boost factor is unity by construction
                    b_extrap = np.ones_like(custom_k_below)
                    bvals[i]= np.concatenate((b_extrap, bvals[i]))

                if do_extrapolate_above:
                    # We extrapolate by setting all b(k > k_max) to b(k_max)
                    b_extrap = bvals[i][-1] * np.ones_like(custom_k_above)
                    bvals[i] = np.concatenate((bvals[i], b_extrap))

            else:
                bvals[i] = 10.**tmp.reshape(k_shape)

        if not(custom_kvec is None):       # This could probably be done cleaner!
            kvals = custom_kvec

        return kvals, bvals


    def get_plin(self, emu_pars_dict, custom_kvec, redshifts):

        if _Class.__module__ not in _sys.modules:
            print("You have not imported neither classee nor classy.\n \
                   Computing linear power spectrum is hence not possible.")
            return None

        # Convert single redshift input argument to array
        if isinstance(redshifts, (int, float)):
            redshifts = np.asarray([redshifts])
        else:
            redshifts = np.asarray(redshifts)

        for z in redshifts:
            assert z <= 10.0 and z>=0.0, "EuclidEmulator2 allows only redshifts in the interval [0.0, 10.0]"

        # Convert single kvec input argument to array
        if isinstance(custom_kvec, (int, float)):
            custom_kvec = np.asarray([custom_kvec])
        else:
            custom_kvec = np.asarray(custom_kvec)

        # "Stringify" the input arrays to be understandable for classy.
        z_str=str(redshifts[0])

        if len(redshifts)>1:
          for i in range(1,len(redshifts)):
            z_str+=','+str(redshifts[i])

        # Convert the input dictionary into a Class-compatible dictionary
        cosmo_par = self.convert_to_emu(emu_pars_dict)

        # When no value for Omega_cdm is given, this is computed
        if cosmo_par['Omega_cdm']==0:
            cosmo_par['Omega_cdm']=cosmo_par['Omega_m']-cosmo_par['Omega_b']-cosmo_par['m_ncdm']/93.14/cosmo_par['h']**2

        # This parameter is eliminated as classy does not accept it.
        del cosmo_par['Omega_m']

        # Extend the input Class-compatible dictionary by the additional
        # information requested by classy.
        classy_pars = cosmo_par
        classy_pars['Omega_Lambda'] = 0.0
        classy_pars['output'] = 'mPk'
        classy_pars['P_k_max_1/Mpc'] = custom_kvec[-1]*cosmo_par['h']
        classy_pars['z_pk'] = z_str
        # Assuming a single massive neutrino with all the mass
        classy_pars['N_ur']=2.0308
        classy_pars['N_ncdm']=1


        # Create a "Class" instance called "cosmo" and run classy to compute
        # the cosmological quantities.
        cosmo = _Class()
        cosmo.set(classy_pars)
        cosmo.compute()

        # Convert k units: h/Mpc to 1/Mpc
        h = classy_pars['h']
        k_classy_arr = h*custom_kvec

        # Get shape of k vector
        k_shape = k_classy_arr.shape

        # Get power spectrum at tabulated z and k in units of Mpc^3
        linpower = {i:
                        np.array([cosmo.pk(k, z)*h*h*h
                                   for k in k_classy_arr]).reshape(k_shape)
                        for i, z in enumerate(redshifts)}

        return custom_kvec, linpower

    def get_pnonlin(self,emu_pars_dict, redshifts, custom_kvec=None):

        if _Class.__module__ not in _sys.modules:
            print("You have not imported neither classee nor classy.\n \
                   Emulating full power spectrum is hence not possible.")
            return None

        # Convert single redshift input argument to array
        if isinstance(redshifts, (int, float)):
            redshifts = np.asarray([redshifts])
        else:
            redshifts = np.asarray(redshifts)

        kvec, Bk = self.get_boost(emu_pars_dict, redshifts, custom_kvec)

        plin = self.get_plin(emu_pars_dict, kvec, redshifts)
        plin = plin[1]

        pnonlin = {}
        for i, z in enumerate(redshifts):
                pnonlin[i] = plin[i]*Bk[i]

        return kvec, pnonlin, plin, Bk

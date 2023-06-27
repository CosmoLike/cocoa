import os
import numpy as np
import ctypes as ct

from _plenslike import *

datadir = os.path.dirname(__file__) + "/data/"
pll     = ct.CDLL( os.path.dirname(__file__) + "/_plenslike.so")

# qest
class qest(ct.Structure):
    _fields_ = [ ("ntrm", ct.c_int),
                 ("lmax", ct.c_int),
                 ("s12L", ct.POINTER( ct.POINTER(ct.c_int) )),
                 ("w12L", ct.POINTER( ct.POINTER( ct.POINTER(ct.c_double) ))) ]

pll.fill_qe_resp.argtypes = [ ct.c_int, ct.POINTER(ct.c_double),
                              ct.POINTER(qest), ct.POINTER(qest),
                              ct.POINTER(ct.c_double), ct.c_int,
                              ct.POINTER(ct.c_double), ct.c_int ]

# mono
class plenslike_dat_mono(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl_fid",        ct.POINTER(ct.c_double)),
                 ("fl",            ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("al_inv",        ct.POINTER(ct.c_double)) ]

pll.load_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono), ct.c_char_p]
pll.free_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono) ]
pll.fill_qe_plm_resp_plm_mono.argtypes = [ ct.c_int, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                                           ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                                           ct.POINTER(ct.c_double), ct.POINTER(ct.c_double) ]
pll.calc_plenslike_mono.restype        = ct.c_double
pll.calc_plenslike_mono_renorm.restype = ct.c_double

class mono():
    def __init__(self, fname):
        print "plenslike:: loading mono likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_mono()
        pll.load_plenslike_dat_mono( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_mono(self.dat)

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmax )
        return pll.calc_plenslike_mono( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm(self, clpp, cltt, bl):
        assert( len(clpp) >= self.dat.lmax )
        assert( len(cltt) >= self.dat.lmax )
        assert( len(bl)   >= self.dat.lmax )

        return pll.calc_plenslike_mono_renorm( ct.byref(self.dat),
                                               clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               cltt.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               bl.ctypes.data_as(   ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt):
        bl = np.array( [self.dat.bl_fid[l] for l in xrange(0, self.dat.lmax+1)] )
        return self.calc_like_renorm(clpp, cltt, bl)

    def calc_clpp_renorm_cltt(self, clpp, cltt):
        resp = np.zeros( len(clpp) )
        
        pll.fill_qe_plm_resp_plm_mono( len(resp)-1, resp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                       self.dat.cltt_fid, self.dat.bl_fid, self.dat.fl,
                                       cltt.ctypes.data_as( ct.POINTER(ct.c_double) ), self.dat.bl_fid )

        return resp**2 / np.array( [ self.dat.al_inv[l] for l in xrange(0, len(clpp)) ] )**2 * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_mono_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins

# quad
class plenslike_dat_quad(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("lmaxt",         ct.c_int),
                 ("lmax1",         ct.c_int),
                 ("lmax2",         ct.c_int),
                 ("lmax3",         ct.c_int),
                 ("lmax4",         ct.c_int),
                 ("s4hat",         ct.c_double),
                 ("s4std",         ct.c_double),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("rl_inv",        ct.POINTER(ct.c_double)),
                 ("sl_fid",        ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl1n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl2n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl3n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl4n1_fid",     ct.POINTER(ct.c_double)),
                 ("fl1",           ct.POINTER(ct.c_double)),
                 ("fl2",           ct.POINTER(ct.c_double)),
                 ("fl3",           ct.POINTER(ct.c_double)),
                 ("fl4",           ct.POINTER(ct.c_double)),
                 ("qe12",          ct.POINTER(qest)),
                 ("qe34",          ct.POINTER(qest)) ]

pll.load_plenslike_dat_quad.argtypes   = [ ct.POINTER(plenslike_dat_quad), ct.c_char_p]
pll.free_plenslike_dat_quad.argtypes   = [ ct.POINTER(plenslike_dat_quad) ]

pll.fill_quad_resp_pp_cltt.argtypes    = [ ct.c_int, ct.POINTER(ct.c_double), ct.POINTER(plenslike_dat_quad), ct.POINTER(ct.c_double) ]

pll.calc_plenslike_quad.restype        = ct.c_double
pll.calc_plenslike_quad_renorm_cltt.restype = ct.c_double

class quad():
    def __init__(self, fname):
        print "plenslike:: loading quad likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_quad()
        pll.load_plenslike_dat_quad( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_quad(self.dat)

    def calc_qc_resp_pp_cltt(self, lmax, cltt):
        assert( len(cltt) >= self.dat.lmaxt )
        
        ret = np.zeros(lmax+1)
        pll.fill_quad_resp_pp_cltt( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                                    ct.byref(self.dat), cltt.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmax )
        return pll.calc_plenslike_quad( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt):
        assert( len(clpp) >= self.dat.lmax )
        assert( len(cltt) >= self.dat.lmaxt )

        return pll.calc_plenslike_quad_renorm_cltt( ct.byref(self.dat),
                                                    clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                                    cltt.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_clpp_renorm_cltt(self, clpp, cltt):
        resp = self.calc_qc_resp_pp_cltt( len(clpp)-1, cltt )
        return resp * np.array( [ self.dat.rl_inv[l] for l in xrange(0, len(clpp)) ] ) * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_quad_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins

# qecl
class plenslike_dat_qecl(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmaxcmb",       ct.c_int),
                 ("lmaxphi",       ct.c_int),
                 ("nqe",           ct.c_int),
                 ("nx",            ct.c_int),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("clee_fid",      ct.POINTER(ct.c_double)),
                 ("clte_fid",      ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("qlpp_fid",      ct.POINTER(ct.c_double)),
                 ("vlpp_inv",      ct.POINTER(ct.c_double)),
                 ("qlpp_inv",      ct.POINTER(ct.c_double)),
                 ("qes",           ct.POINTER(ct.POINTER(qest))),
                 ("qefs",          ct.POINTER(ct.c_int)),
                 ("qe12",          ct.POINTER(ct.c_int)),
                 ("qe34",          ct.POINTER(ct.c_int)) ]

pll.load_plenslike_dat_qecl.argtypes   = [ ct.POINTER(plenslike_dat_qecl), ct.c_char_p]
pll.free_plenslike_dat_qecl.argtypes   = [ ct.POINTER(plenslike_dat_qecl) ]

pll.calc_plenslike_qecl.restype        = ct.c_double
pll.calc_plenslike_qecl_renorm.restype = ct.c_double

class qecl():
    def __init__(self, fname):
        print "plenslike:: loading qecl likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_qecl()
        ret = pll.load_plenslike_dat_qecl( ct.byref(self.dat), fname )
        if ret != 0:
            self.dat = None
            assert(0)

    def __del__(self):
        if self.dat != None:
            pll.free_plenslike_dat_qecl(self.dat)

    def calc_qc_resp_pp(self, lmax, cltt, clee, clte):
        assert( len(cltt) >= self.dat.lmaxcmb )

        ret = np.zeros(lmax+1)
        pll.fill_qecl_resp_pp_cls( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   ct.byref(self.dat),
                                   cltt.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   clee.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   clte.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmaxpphi )
        return pll.calc_plenslike_qecl( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt, clee, clte):
        assert( len(clpp) >= self.dat.lmaxphi )
        assert( len(cltt) >= self.dat.lmaxcmb )
        assert( len(clee) >= self.dat.lmaxcmb )
        assert( len(clte) >= self.dat.lmaxcmb )

        return pll.calc_plenslike_qecl_renorm( ct.byref(self.dat),
                                               clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               cltt.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               clee.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               clte.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_clpp_renorm(self, clpp, cltt):
        resp = self.calc_qc_resp_pp_cls( len(clpp)-1, cltt )
        return resp * np.array( [ self.dat.rlpp_inv[l] for l in xrange(0, len(clpp)) ] ) * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_qecl_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins

# full
class plenslike_dat_full(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmaxcmb",       ct.c_int),
                 ("lmaxphi",       ct.c_int),
                 ("n1lqbins",      ct.c_int),
                 ("n1lpbins",      ct.c_int),
                 ("nqe",           ct.c_int),
                 ("nx",            ct.c_int),
                 ("s4fid",         ct.c_double),
                 ("s4std",         ct.c_double),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("clee_fid",      ct.POINTER(ct.c_double)),
                 ("clte_fid",      ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("qlpp_fid",      ct.POINTER(ct.c_double)),
                 ("qlss_fid",      ct.POINTER(ct.c_double)),
                 ("n1pp_fid",      ct.POINTER(ct.c_double)),
                 ("vlpp_inv",      ct.POINTER(ct.c_double)),
                 ("qlpp_inv",      ct.POINTER(ct.c_double)),
                 ("n1lqs",         ct.POINTER(ct.c_double)),
                 ("n1lps",         ct.POINTER(ct.c_double)),
                 ("mat_n1",        ct.POINTER(ct.c_double)),
                 ("qes",           ct.POINTER(ct.POINTER(qest))),
                 ("qefs",          ct.POINTER(ct.c_int)),
                 ("qe12",          ct.POINTER(ct.c_int)),
                 ("qe34",          ct.POINTER(ct.c_int)) ]

pll.load_plenslike_dat_full.argtypes   = [ ct.POINTER(plenslike_dat_full), ct.c_char_p ]
pll.free_plenslike_dat_full.argtypes   = [ ct.POINTER(plenslike_dat_full) ]

pll.calc_plenslike_full.restype        = ct.c_double
pll.calc_plenslike_full_renorm.restype = ct.c_double

class full():
    def __init__(self, fname):
        print "plenslike:: loading full likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_full()
        ret = pll.load_plenslike_dat_full( ct.byref(self.dat), fname )
        if ret != 0:
            self.dat = None
            assert(0)

    def __del__(self):
        if self.dat != None:
            pll.free_plenslike_dat_full(self.dat)

    def calc_n1(self, lmax, clpp):
        assert( len(clpp) >= self.dat.lmaxphi )

        ret = np.zeros(lmax+1)
        pll.fill_full_n1( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                          ct.byref(self.dat),
                          clpp.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_qc_resp_pp(self, lmax, cltt, clee, clte):
        assert( len(cltt) >= self.dat.lmaxcmb )

        ret = np.zeros(lmax+1)
        pll.fill_full_resp_pp_cls( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   ct.byref(self.dat),
                                   cltt.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   clee.ctypes.data_as(ct.POINTER(ct.c_double)),
                                   clte.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmaxpphi )
        return pll.calc_plenslike_full( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt, clee, clte):
        assert( len(clpp) >= self.dat.lmaxphi )
        assert( len(cltt) >= self.dat.lmaxcmb )
        assert( len(clee) >= self.dat.lmaxcmb )
        assert( len(clte) >= self.dat.lmaxcmb )

        return pll.calc_plenslike_full_renorm( ct.byref(self.dat),
                                               clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               cltt.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               clee.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               clte.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_clpp_renorm(self, clpp, cltt):
        resp = self.calc_qc_resp_pp_cls( len(clpp)-1, cltt )
        return resp * np.array( [ self.dat.rlpp_inv[l] for l in xrange(0, len(clpp)) ] ) * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_full_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins

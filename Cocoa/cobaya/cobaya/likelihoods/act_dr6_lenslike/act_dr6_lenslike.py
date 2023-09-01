import numpy as np
import warnings
from scipy.interpolate import interp1d
try:
    from cobaya.likelihoods.base_classes import InstallableLikelihood
except:
    InstallableLikelihood = object
import os
file_dir = os.path.abspath(os.path.dirname(__file__))
data_dir = f"{file_dir}/data/v1.1/"

variants =[x.strip() for x in  '''
act_baseline,
act_extended,
actplanck_baseline,
actplanck_extended,
act_polonly,
act_cibdeproj,
act_cinpaint
'''.strip().replace('\n','').split(',')]


# ================
# HELPER FUNCTIONS
# ================

def pp_to_kk(clpp,ell):
    return clpp * (ell*(ell+1.))**2. / 4.
    
def get_corrected_clkk(data_dict,clkk,cltt,clte,clee,clbb,suff=''):
    clkk_fid = data_dict['fiducial_cl_kk']
    cl_dict = {'tt':cltt,'te':clte,'ee':clee,'bb':clbb}
    N1_kk_corr = data_dict[f'dN1_kk{suff}'] @ (clkk-clkk_fid)
    dNorm = data_dict[f'dAL_dC{suff}']
    fid_norm = data_dict[f'fAL{suff}']
    N1_cmb_corr = 0.
    norm_corr = 0.
    for i,s in enumerate(['tt','ee','bb','te']):
        cldiff = (cl_dict[s]-data_dict[f'fiducial_cl_{s}'])
        N1_cmb_corr = N1_cmb_corr + (data_dict[f'dN1_{s}{suff}']@cldiff)
        c = - 2. * (dNorm[i] @ cldiff)
        if i==0:
            ls = np.arange(c.size)
        c[ls>=2] = c[ls>=2] / fid_norm[ls>=2]
        norm_corr = norm_corr + c
    nclkk = clkk + norm_corr*clkk_fid + N1_kk_corr + N1_cmb_corr
    return nclkk

def standardize(ls,cls,trim_lmax,lbuffer=2,extra_dims="y"):
    cstart = int(ls[0])
    diffs = np.diff(ls)
    if not(np.all(np.isclose(diffs,1.))): raise ValueError("Multipoles are not spaced by 1")
    if not(cstart<=2): raise ValueError("Multipoles start at value greater than 2")
    nlen = trim_lmax+lbuffer
    cend = nlen - cstart
    if extra_dims=="xyy":
        out = np.zeros((cls.shape[0],nlen,nlen))
        out[:,cstart:,cstart:] = cls[:,:cend,:cend]
    elif extra_dims=="yy":
        out = np.zeros((nlen,nlen))
        out[cstart:,cstart:] = cls[:cend,:cend]
    elif extra_dims=="xy":
        out = np.zeros((cls.shape[0],nlen))
        out[:,cstart:] = cls[:,:cend]
    elif extra_dims=="y":
        out = np.zeros(nlen)
        out[cstart:] = cls[:cend]
    else:
        raise ValueError
    return out

def get_limber_clkk_flat_universe(results,Pfunc,lmax,kmax,nz,zsrc=None):
    # Adapting code from Antony Lewis' CAMB notebook
    if zsrc is None:
        chistar = results.conformal_time(0)- results.tau_maxvis
    else:
        chistar = results.comoving_radial_distance(zsrc)
    chis = np.linspace(0,chistar,nz)
    zs=results.redshift_at_comoving_radial_distance(chis)
    dchis = (chis[2:]-chis[:-2])/2
    chis = chis[1:-1]
    zs = zs[1:-1]
    
    #Get lensing window function (flat universe)
    win = ((chistar-chis)/(chis**2*chistar))**2
    #Do integral over chi
    ls = np.arange(0,lmax+2, dtype=np.float64)
    cl_kappa=np.zeros(ls.shape)
    w = np.ones(chis.shape) #this is just used to set to zero k values out of range of interpolation
    for i, l in enumerate(ls[2:]):
        k=(l+0.5)/chis
        w[:]=1
        w[k<1e-4]=0
        w[k>=kmax]=0
        cl_kappa[i+2] = np.dot(dchis, w*Pfunc.P(zs, k, grid=False)*win/k**4)
    cl_kappa*= (ls*(ls+1))**2
    return cl_kappa

def get_camb_lens_obj(nz,kmax,zmax=None):
    import camb
    pars = camb.CAMBparams()
    # This cosmology is purely to go from chis->zs for limber integration;
    # the details do not matter
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    pars.InitPower.set_params(ns=0.965)
    results= camb.get_background(pars)
    nz = nz
    if zmax is None:
        chistar = results.conformal_time(0)- results.tau_maxvis
    else:
        chistar = results.comoving_radial_distance(zmax)
    chis = np.linspace(0,chistar,nz)
    zs=results.redshift_at_comoving_radial_distance(chis)
    cobj = {"CAMBdata": None,
            "Pk_interpolator": { "z": zs,
                                 "k_max": kmax,
                                 "nonlinear": True,
                                 "vars_pairs": ([["Weyl", "Weyl"]])}}
    return cobj


def parse_variant(variant):

    variant = variant.lower().strip()
    if variant not in variants: raise ValueError

    v = None
    if '_extended' in variant:
        baseline = False
    else:
        baseline = True
        if '_baseline' not in variant:
            v = variant.split('_')[-1]

    include_planck = True if 'actplanck' in variant else False
    return v,baseline,include_planck

# ==================            
# Generic likelihood
# ==================

"""
data_dict = load_data(data_directory) # pre-load data
# for each predicted spectra in chain
# cl_kk is CMB lensing convergence power spectrum (dimensionless, 
# no ell or 2pi factors)
# cl_tt, cl_ee, cl_te, cl_bb are lensed CMB power spectra
# (muK^2 units, no ell or 2pi factors)
lnlike = generic_lnlike(data_dict,cl_kk,cl_tt,cl_ee,cl_te,cl_bb)
This returns ln(Likelihood)
so for example,
chi_square = -2 lnlike
"""

def load_data(variant, ddir=data_dir,
              lens_only=False,
              apply_hartlap=True,like_corrections=True,mock=False,
              nsims_act=796,nsims_planck=400,trim_lmax=2998,scale_cov=None):
    """
    Given a data directory path, this function loads into a dictionary
    the data products necessary for evaluating the DR6 lensing likelihood.
    This includes:
    1. the ACT lensing bandpowers. Planck lensing bandpowers will be 
    appended if include_planck is True.
    2. the associated binning matrix to be applied to a theory curve
    3. the associated covariance matrix
    4. data products associated with applying likelihood corrections

    All these products will be standardized so that they apply
    to theory curves specified from L=0 to trim_lmax.

    A Hartlap correction will be applied to the covariance matrix
    corresponding to the lower of the number of simulations involved.
    
    """
    # TODO: review defaults


    v,baseline,include_planck = parse_variant(variant)
    

    # output data
    d = {}

    if lens_only and like_corrections: raise ValueError("Likelihood corrections should not be used in lens_only runs.")
    if not(lens_only) and not(like_corrections):
        warnings.warn("Neither using CMB-marginalized covariance matrix nor including likelihood corrections. Effective covariance may be underestimated.")

    d['include_planck'] = include_planck
    d['likelihood_corrections'] = like_corrections

    # Fiducial spectra
    if like_corrections:
        f_ls, f_tt, f_ee, f_bb, f_te = np.loadtxt(f"{ddir}/like_corrs/cosmo2017_10K_acc3_lensedCls.dat",unpack=True)
        f_tt = f_tt / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_ee = f_ee / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_bb = f_bb / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_te = f_te / (f_ls * (f_ls+1.)) * 2. * np.pi

        fd_ls, f_dd = np.loadtxt(f"{ddir}/like_corrs/cosmo2017_10K_acc3_lenspotentialCls.dat",unpack=True,usecols=[0,5])
        f_kk = f_dd * 2. * np.pi / 4.
        d['fiducial_cl_tt'] = standardize(f_ls,f_tt,trim_lmax)
        d['fiducial_cl_te'] = standardize(f_ls,f_te,trim_lmax)
        d['fiducial_cl_ee'] = standardize(f_ls,f_ee,trim_lmax)
        d['fiducial_cl_bb'] = standardize(f_ls,f_bb,trim_lmax)
        d['fiducial_cl_kk'] = standardize(fd_ls,f_kk,trim_lmax)

        
    # Return data bandpowers, covariance matrix and binning matrix
    if baseline:
        start = 2
        end = -6
    else:
        start = 2
        end = -3

    if v is None:
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act.txt')
    elif v=='cinpaint':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_cinpaint.txt')
    elif v=='polonly':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_polonly.txt')
    elif v=='cibdeproj':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_cibdeproj.txt')
    nbins_tot_act = y.size
    d['full_data_binned_clkk_act'] = y.copy()
    data_act = y[start:end].copy()
    d['data_binned_clkk'] = data_act
    nbins_act = data_act.size

    binmat = np.loadtxt(f'{ddir}/binning_matrix_act.txt')
    d['full_binmat_act'] = binmat.copy()
    pells = np.arange(binmat.shape[1])
    bcents = binmat@pells
    ls = np.arange(binmat.shape[1])
    d['binmat_act'] = standardize(ls,binmat[start:end,:],trim_lmax,extra_dims="xy")
    d['bcents_act'] = bcents.copy()

    if lens_only:
        if include_planck:
            if v not in [None,'cinpaint']: raise ValueError(f"Combination of {v} with Planck is not available")
            fcov = np.loadtxt(f'{ddir}/covmat_actplanck_cmbmarg.txt')
        else:
            if v=='cibdeproj':
                fcov = np.loadtxt(f"{ddir}/covmat_act_cibdeproj_cmbmarg.txt")
            elif v=='pol':
                fcov = np.loadtxt(f"{ddir}/covmat_act_polonly_cmbmarg.txt")
            else:
                fcov = np.loadtxt(f"{ddir}/covmat_act_cmbmarg.txt")
    else:
        if v not in [None,'cinpaint']: raise ValueError(f"Covmat for {v} without CMB marginalization is not available")
        if include_planck:
            fcov = np.loadtxt(f'{ddir}/covmat_actplanck.txt')
        else:
            fcov = np.loadtxt(f'{ddir}/covmat_act.txt')

    d['full_act_cov'] = fcov.copy()
    # Remove trailing bins from ACT part
    sel = np.s_[nbins_tot_act+end:nbins_tot_act]
    cov = np.delete(np.delete(fcov,sel,0),sel,1)
    # Remove leading bins from ACT part
    sel = np.s_[:start]
    cov = np.delete(np.delete(cov,sel,0),sel,1)
    
    # Test
    covmat = np.loadtxt(f'{ddir}/covmat_act.txt')
    covmat1 = covmat[start:end,start:end]
    cdiff = cov[:nbins_act,:nbins_act] - covmat1
    if not(np.all(np.isclose(cdiff,0))): raise ValueError

    if include_planck:
        data_planck = np.loadtxt(f'{ddir}/clkk_bandpowers_planck.txt')
        d['data_binned_clkk'] = np.append(d['data_binned_clkk'],data_planck)
        binmat = np.loadtxt(f'{ddir}/binning_matrix_planck.txt')
        pells = np.arange(binmat.shape[1])
        bcents = binmat@pells
        ls = np.arange(binmat.shape[1])
        d['binmat_planck'] = standardize(ls,binmat,trim_lmax,extra_dims="xy")
        d['bcents_planck'] = bcents.copy()

    if like_corrections:
        # Load matrices
        cmat = np.load(f"{ddir}/like_corrs/norm_correction_matrix_Lmin0_Lmax4000.npy")
        ls = np.arange(cmat.shape[1])
        d['dAL_dC'] = standardize(ls,cmat,trim_lmax,extra_dims="xyy")
        if include_planck:
            cmat = np.load(f"{ddir}/like_corrs/P18_norm_correction_matrix_Lmin0_Lmax3000.npy")
            ls = np.arange(cmat.shape[1])
            d['dAL_dC_planck'] = standardize(ls,cmat,trim_lmax,extra_dims="xyy")
            

        fAL_ls,fAL = np.loadtxt(f"{ddir}/like_corrs/n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt")
        d['fAL'] = standardize(fAL_ls,fAL,trim_lmax,extra_dims="y")
        if include_planck:
            fAL_ls,fAL = np.loadtxt(f"{ddir}/like_corrs/PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt")
            d['fAL_planck'] = standardize(fAL_ls,fAL,trim_lmax,extra_dims="y")

        for spec in ['kk','tt','ee','bb','te']:
            n1mat = np.loadtxt(f"{ddir}/like_corrs/N1der_{spec.upper()}_lmin600_lmax3000_full.txt")
            d[f'dN1_{spec}'] = standardize(fAL_ls,n1mat,trim_lmax,extra_dims="yy")
            if include_planck:
                n1mat = np.loadtxt(f"{ddir}/like_corrs/N1_planck_der_{spec.upper()}_lmin100_lmax2048.txt")
                d[f'dN1_{spec}_planck'] = standardize(fAL_ls,n1mat,trim_lmax,extra_dims="yy")

    nbins = d['data_binned_clkk'].size
    nsims = min(nsims_act,nsims_planck) if include_planck else nsims_act
    hartlap_correction = (nsims-nbins-2.)/(nsims-1.)
    if apply_hartlap:
        pass
        #warnings.warn("Hartlap correction to cinv: ", hartlap_correction)
    else:
        warnings.warn(f"Disabled Hartlap correction to cinv: {hartlap_correction}")
        hartlap_correction = 1.0
    if scale_cov is not None:
        warnings.warn(f"Covariance has been artificially scaled by: {scale_cov}")
        cov = cov * scale_cov
    d['cov'] = cov
    cinv = np.linalg.inv(cov) * hartlap_correction
    d['cinv'] = cinv

    if mock:
        mclpp = np.loadtxt(f"{self.ddir}/cls_default_dr6_accuracy.txt",usecols=[5])
        ls = np.arange(mclpp.size)
        mclkk = mclpp * 2. * np.pi / 4.
        self.clkk_data = self.binning_matrix @ mclkk[:self.kLmax]
    
    return d
    

def generic_lnlike(data_dict,ell_kk,cl_kk,ell_cmb,cl_tt,cl_ee,cl_te,cl_bb,trim_lmax = 2998):

    cl_kk = standardize(ell_kk,cl_kk,trim_lmax)
    cl_tt = standardize(ell_cmb,cl_tt,trim_lmax)
    cl_ee = standardize(ell_cmb,cl_ee,trim_lmax)
    cl_bb = standardize(ell_cmb,cl_bb,trim_lmax)
    cl_te = standardize(ell_cmb,cl_te,trim_lmax)
    
    d = data_dict
    cinv = d['cinv']
    clkk_act = get_corrected_clkk(data_dict,cl_kk,cl_tt,cl_te,cl_ee,cl_bb) if d['likelihood_corrections'] else cl_kk
    bclkk = d['binmat_act'] @ clkk_act
    if d['include_planck']:
        clkk_planck = get_corrected_clkk(data_dict,cl_kk,cl_tt,cl_te,cl_ee,cl_bb,'_planck') if d['likelihood_corrections'] else cl_kk
        bclkk = np.append(bclkk, d['binmat_planck'] @ clkk_planck)
    delta = d['data_binned_clkk'] - bclkk
    lnlike = -0.5 * np.dot(delta,np.dot(cinv,delta))
    return lnlike

    
# =================            
# Cobaya likelihood
# =================


class ACTDR6LensLike(InstallableLikelihood):

    lmax: int
    mock = False
    nsims_act = 792. # Number of sims used for covmat; used in Hartlap correction
    nsims_planck = 400. # Number of sims used for covmat; used in Hartlap correction
    no_like_corrections = False
    lens_only = False
    # Any ells above this will be discarded; likelihood must at least request ells up to this
    trim_lmax = 2998
    variant = "act_baseline"
    apply_hartlap = True
    # Limber integral parameters
    limber = False
    nz = 100
    kmax = 10
    scale_cov = None
    alens = False # Whether to divide the theory spectrum by Alens

    def initialize(self):
        if self.lens_only: self.no_like_corrections = True
        if self.lmax<self.trim_lmax: raise ValueError(f"An lmax of at least {self.trim_lmax} is required.")
        self.data = load_data(variant=self.variant,lens_only=self.lens_only,
                              like_corrections=not(self.no_like_corrections),apply_hartlap=self.apply_hartlap,
                              mock=self.mock,nsims_act=self.nsims_act,nsims_planck=self.nsims_planck,
                              trim_lmax=self.trim_lmax,scale_cov=self.scale_cov)
        
        if self.no_like_corrections:
            self.requested_cls = ["pp"]
        else:
            self.requested_cls = ["tt", "te", "ee", "bb", "pp"]

    def get_requirements(self):
        if self.no_like_corrections:
            ret = {'Cl': {'tt': self.lmax,'te': self.lmax,'ee': self.lmax,'pp':self.lmax}}
        else:
            ret = {'Cl': {'pp':self.lmax}}

        if self.limber:
            cobj = get_camb_lens_obj(self.nz,self.kmax,self.zmax)
            ret.update(cobj)
            
        return ret

    def logp(self, **params_values):
        cl = self.theory.get_Cl(ell_factor=False, units='FIRASmuK2')
        return self.loglike(cl, **params_values)

    def get_limber_clkk(self,**params_values):
        Pfunc = self.theory.get_Pk_interpolator(var_pair=("Weyl", "Weyl"), nonlinear=True, extrap_kmax=30.)
        results = self.provider.get_CAMBdata()
        return get_limber_clkk_flat_universe(results,Pfunc,self.trim_lmax,self.kmax,nz,zstar=None)

    def loglike(self, cl, **params_values):
        ell = cl['ell']
        Alens = 1
        if self.alens:
            Alens = self.theory.get_param('Alens')
        clpp = cl['pp'] / Alens
        if self.limber:
            cl_kk = self.get_limber_clkk( **params_values)
        else:
            cl_kk = pp_to_kk(clpp,ell)
        
        logp = generic_lnlike(self.data,ell,cl_kk,ell,cl['tt'],cl['ee'],cl['te'],cl['bb'],self.trim_lmax)
        self.log.debug(
            f"ACT-DR6-lensing-like lnLike value = {logp} (chisquare = {-2 * logp})")
        return logp

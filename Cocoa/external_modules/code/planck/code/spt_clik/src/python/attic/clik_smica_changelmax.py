#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import clik.hpy as hpy
import clik.smicahlp as hlp
import clik.parobject as php
import numpy as nm


def cutlminlmax(nlmin, nlmax, infile, outfile):
    """
    NAME
           cutlminlmax - cut l range

    SYNOPSIS
           cutlminlmax lmin_new lmax_new file_in file_out

    DESCRIPTION
           Cut the temperature component of a clik file to a new
           multipole range lmin <= l <= lmax, where lmin, lmax is
           adapted to the binning scheme of the file.
    """
    hpy.copyfile(infile,outfile)
    ff     = hpy.File(outfile,"r+")
    lmaxs  = ff["clik"].attrs["lmax"]
    lmin   = ff["clik/lkl_0"].attrs["lmin"]
    blmins = ff["clik/lkl_0/bin_lmin"][:]
    nbins = len(blmins)
    blmaxs = ff["clik/lkl_0/bin_lmax"][:]
    nlmin  = max(nlmin, lmin)
    nlmax  = min(nlmax, lmaxs[0])
    idmin  = blmins.searchsorted(nlmin-lmin)
    idmax  = blmaxs.searchsorted(nlmax-lmin)
    print lmin,lmaxs[0],nlmin,nlmax,idmin,idmax

    if idmin > idmax:
        print "No valid binning scheme found for input values" \
        + " lmin={lmin}, lmax={lmax}".format(lmin=nlmin, lmax=nlmax)
        ff.close()
        return

    print "Cut at lmin = %d"%(blmins[idmin]+lmin)
    print "Cut at lmax = %d"%(blmaxs[idmax]+lmin)
    nnlmin   = (blmins[idmin]+lmin)
    nnlmax   = (blmaxs[idmax]+lmin)
    lmaxs[0] = nnlmax
    ff["clik"].attrs["lmax"] = lmaxs
    ff["clik/lkl_0"].attrs["lmin"] = nnlmin
    ff["clik/lkl_0"].attrs["lmax"] = nnlmax
    zerobin = blmins[idmin]
    
    del ff["clik/lkl_0/bin_lmin"]
    ff["clik/lkl_0/bin_lmin"] = blmins[idmin:idmax+1]-zerobin
    del ff["clik/lkl_0/bin_lmax"]
    ff["clik/lkl_0/bin_lmax"] = blmaxs[idmin:idmax+1]-zerobin

    wl = ff["clik/lkl_0/bin_ws"][:]
    del ff["clik/lkl_0/bin_ws"]
    ff["clik/lkl_0/bin_ws"] = wl[nnlmin-lmin:nnlmax-lmin+1]
    ff["clik/lkl_0"].attrs["nbins"] = idmax-idmin+1

    rqh = ff["clik/lkl_0/Rq_hat"][:]
    nch = ff["clik/lkl_0"].attrs["m_channel_T"]
    rqh.shape=[-1,nch,nch]
    rqh = rqh[idmin:idmax+1,:,:]
    del ff["clik/lkl_0/Rq_hat"]
    ff["clik/lkl_0/Rq_hat"] = rqh.flatten()
    try :
        criterion = ff["clik/lkl_0/criterion_eig_norm"]
        del ff["clik/lkl_0/criterion_eig_norm"]
        ff["clik/lkl_0/criterion_eig_norm"] = criterion[idmin:idmax+1]
    except Exception,e:
        pass
    try :
        criterion = ff["clik/lkl_0/criterion_quad_mat"][:]
        del ff["clik/lkl_0/criterion_quad_mat"]
        sn = nm.sqrt(len(criterion)/nbins)
        criterion.shape=(nbins,sn,sn)
        ff["clik/lkl_0/criterion_quad_mat"] = (criterion[idmin:idmax+1]*1.).flat[:]
    except Exception,e:
        pass


    wq = ff["clik/lkl_0/wq"][:]
    del ff["clik/lkl_0/wq"]
    ff["clik/lkl_0/wq"] = wq[idmin:idmax+1]

    for i in range(1,ff["clik/lkl_0"].attrs["n_component"]):
        cmpt = "clik/lkl_0/component_%d"%i
        try:
            ff[cmpt].attrs["lmax"]=nnlmax
            ff[cmpt].attrs["lmin"]=nnlmin
        except Exception,e:
            pass
        try:
            rqh = ff[cmpt+"/Rq_0"][:]
            rqh.shape=[-1,nch,nch]
            rqh = rqh[idmin:idmax+1]
            del ff[cmpt+"/Rq_0"]
            ff[cmpt+"/Rq_0"] = rqh.flatten()
        except Exception,e:
            pass
        try:
            ng = nm.sum(ff[cmpt].attrs["ngcal"])
            blmn = nnlmin-lmin
            blmx = nnlmax-lmin
            if ff[cmpt].attrs["binned"]:
                blmn = idmin
                blmx = idmax
            tpl = ff[cmpt+"/gcaltpl"][:]
            del ff[cmpt+"/gcaltpl"]
            ff[cmpt+"/gcaltpl"] = tpl[blmn*ng:(blmx+1)*ng]
        except Exception,e:
            pass

    try:
        del ff["clik/check_param"]
        del ff["clik"].attrs["check_value"]
    except Exception,e:
        pass
    ff.close()


import sys
if __name__=="__main__":
    if len(sys.argv) != 5:
        print cutlminlmax.__doc__
    else:
        cutlminlmax(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],sys.argv[4])

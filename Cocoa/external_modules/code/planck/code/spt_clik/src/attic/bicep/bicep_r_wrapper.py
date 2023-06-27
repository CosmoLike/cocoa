#!/usr/bin/env python
# bicep_r_wrapper.py
# written by  D Barkats and K. Karkare. Adapted from code written by V Buza.
# This is a wrapper around functions to calculate the BICEP1/BICEP2 HL likelihood.
# The wrapper is a specific exemple that shows how to produce the HL likelihood versus r.
# This method is not intended to produce the results presented in the BICEP2 2014 I paper.
#
# This  wrapper will do the following
#      Loads up the bpwf and all the data products from the bicep2 data release page: http:bicepkeck.org/bicep2_2014_release
#      Loads in a BB theory spectrum from CAMB in this case BB r = 0.1
#      Loads the LCDM+lensing spectra from CAMB
#      Loops over r values at which to evaluate the likelihood
#          Scales fiducial BB by r, gets expvals
#          Computes the likelihood 
#      Saves the likelihood for that r
#
# Inputs:
#        - (hardwired) B1 or B2 data products
#        - (hardwired) B1 or B2 window functions 
#        - (hardwired) theory input spectrum
# Outputs:
#        - a text file output of the likelihood
#
# Calls:
#        - bicep_utils.py: all the utility functions.
#
# Current limitations:
#     - B2 data products not available for T, TE, TB so this only works for E, B, EB
#
#$Id: bicep_r_wrapper.py,v 1.1.2.6 2014/03/12 18:25:39 dbarkats Exp $ #

from optparse import OptionParser

# Execution starts here
if __name__ == '__main__':
    usage =   '''
    This is the top level function to produce a HL likelihood on r using BICEP data products.
    Refer to  http://lanl.arxiv.org/abs/1310.1422 for details on BICEP1.
    Refer to  http://bicepkeck.org/bicep2_2014_release for details on BICEP2
    For bug reports and feature requests, contact dbarkats [at] alma.cl or kkarkare [at] cfa.harvard.edu
    '''

    # Initialize parser
    parser = OptionParser()

    parser.add_option("-r", "--rlist", dest="rlist",
                      help="Specify a single value or range of tensor-to-scalar values r to calculate the likelihood for. Default is None and uses the hardwired rlist for the given experiment. Format is a single value of r or a comma-separated range: rmin,rmax,rstep (ie. -1,5,0.001)")

    parser.add_option("-f", "--fields", dest="fields",default = 'B',
                      help="Specify one of the following fields: T,E,B,TE,TB,EB,TEB,  B is default. Only E, B, and EB available now")

    parser.add_option("-e", "--exp", dest="experiment",default = 'bicep2',
                      help="Specify one of the following experiments to work with bicep1 or bicep2. Default is bicep2")
    
    (options, args) = parser.parse_args()

# Imports
import numpy as np
import sys
from numpy import linalg as LA

# BICEP-specific module
import bicep_util as bu  

# define hard-wired defaults
field = options.fields   # field to calculate likelihood on
if field != 'B' and field !='E' and field !='EB':
    print "only field = E, B, or EB are available for now"
    sys.exit()
exp = options.experiment
if  options.rlist == None:
    if exp == 'bicep1':
        rlist = np.arange(-1.2,5.01,0.01) # default r values to evaluate for B1
    elif exp == 'bicep2':
        rlist = np.arange(0,0.751,0.001) # default r values to evaluate for B2
else:
    orlist = options.rlist.split(',')
    if np.size(orlist) == 1:
        rlist = np.array([np.float(options.rlist)])
    elif np.size(orlist) == 3:
        rlist = np.arange(np.float(orlist[0]),np.float(orlist[1])+np.float(orlist[2]),np.float(orlist[2]))
    else:
        print 'Format is a single value of r or a comma separate range: rmin,rmax,rstep (ie. -1,5,0.001) '
        sys.exit()
        
# load the data products and bandpower window function from disk:
C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l =  bu.init(exp, field)

# Load the input r = 0.1 BB spectrum
if exp == 'bicep1':
    file_in = "camb_wmap5yr_withB_uK.txt"
elif exp == 'bicep2':
    file_in = "B2_3yr_camb_planck_withB_uK_20140314.txt"
(ell,BB_Cs_l) = bu.load_cmbfast(file_in)

# Load the LCDM spectra (contains lensing)
if exp == 'bicep1':
    lcdm = "camb_wmap5yr_lensed_uK.txt"
elif exp == 'bicep2':
    lcdm = "B2_3yr_camb_planck_lensed_uK_20140314.txt"
(ell1,inpmodLCDM_Cs_l) = bu.load_cmbfast(lcdm)
expvLCDM = bu.calc_expvals(ell1,inpmodLCDM_Cs_l,bpwf_l,bpwf_Cs_l)

# Initiating the logarithmic likelihood arrays
logLike = np.zeros(len(rlist))

# loop over r values
print "### Evaluating the HL likelihood..."
for i in range(0,len(rlist)):

    expv = expvLCDM.copy()
    
    r  =  rlist[i]
    # Scale the BB spectrum from fiducial model (r=0.1) and get its expectation value
    this_mod = BB_Cs_l.copy()
    this_mod[:,3] = np.multiply(this_mod[:,3],r/0.1)
    expvBB = bu.calc_expvals(ell,this_mod,bpwf_l,bpwf_Cs_l)
    
    # Add the scaled BB to the fixed LCDM + lensing BB 
    expv[:,3] = expv[:,3]+ expvBB[:,3]
    
    # Choose the relevant input spectra expecation value according to field
    if field == "T":
        C_l[:,0,0] = expv[:,0]
    elif field == 'E':
        C_l[:,0,0] = expv [:,2]
    elif field == 'B':
        C_l[:,0,0] = expv[:,3] 
    elif field == "EB":
        C_l[:,0,0] = expv[:,2]
        C_l[:,0,1] = expv[:,5]
        C_l[:,1,0] = expv[:,5]
        C_l[:,1,1] = expv[:,3]
    elif field == "TB":
        C_l[:,0,0] = expv[:,0]
        C_l[:,0,1] = expv[:,4]
        C_l[:,1,0] = expv[:,5]
        C_l[:,1,1] = expv[:,3]
    elif field == "TE":
        C_l[:,0,0] = expv[:,0]
        C_l[:,0,1] = expv[:,1]
        C_l[:,1,0] = expv[:,1]
        C_l[:,1,1] = expv[:,2]
    
    # add Noise bias N_l to expectation values.
    C_l = C_l + N_l

    logL = bu.evaluateLikelihood(C_l,C_l_hat,C_fl,M_inv)
    logLike[i]=logL

    if i == len(rlist)-1:
        sys.stdout.write("Likelihood: %4.3f %8.5f  \n" % (r, logLike[i]) )
    else:
        sys.stdout.write("Likelihood: %4.3f %8.5f  \r" % (r, logLike[i]) )
        sys.stdout.flush()

# Saves the logarithmic likelihood to file
bu.saveLikelihoodToText(rlist,logLike,field, exp = exp)


from getdist import MCSamples
from getdist import loadMCSamples
import getdist
import sys
import numpy as np
import getdist.plots as gplot
import os
import matplotlib
import subprocess

input_mcmc_chain  = "EXAMPLE_MCMC1"
number_chains_output = 8
max_number_chains_per_file = 1500

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
n_per_file = 2*max_number_chains_per_file
THIN_FACTOR = 0

while n_per_file > max_number_chains_per_file:
  THIN_FACTOR = THIN_FACTOR + 1

  readsamps=loadMCSamples("./" + input_mcmc_chain, settings={'ignore_rows':0.3})
  readsamps.thin(THIN_FACTOR)
  readsamps.saveAsText(input_mcmc_chain + '_THIN')

  if number_chains_output > 1:
    original_chain = np.loadtxt(input_mcmc_chain + '_THIN.txt')

    n_per_file = int(len(original_chain)/number_chains_output)
   
    for i in range(number_chains_output-1):
      np.savetxt(input_mcmc_chain + "_THIN." + str(i+1) + ".txt", original_chain[i*n_per_file:(i+1)*n_per_file])

    np.savetxt(input_mcmc_chain + "_THIN." + str(number_chains_output) + ".txt", original_chain[(number_chains_output-1)*n_per_file:])

    subprocess.run(["rm", input_mcmc_chain + '_THIN.txt'])

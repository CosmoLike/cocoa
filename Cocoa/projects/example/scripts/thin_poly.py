from getdist import MCSamples
from getdist import loadMCSamples
import getdist
import sys
import numpy as np
import getdist.plots as gplot
import os
import matplotlib
from sys import exit
import re
import random
import subprocess

input_poly_chain  = "EXAMPLE_POLY2"
number_chains_output = 8

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

input_file_name   = input_poly_chain + ".1.txt"
output_poly_chain = input_poly_chain + "_THIN"
output_file_name  = output_poly_chain + ".1.txt"

subprocess.run(["cp", input_poly_chain + ".input.yaml", output_poly_chain + ".input.yaml"])
subprocess.run(["cp", input_poly_chain + ".updated.yaml", output_poly_chain + ".updated.yaml"])

original_chain = np.loadtxt(input_file_name)
nsample = original_chain.shape[0]

# In polychord, max_weight should be 1
max_weight = max(original_chain[:,0])

for i in range(nsample):    
  if(original_chain[i,0]/max_weight < random.random()) :
    original_chain[i,0] = 0.0
  else :
    original_chain[i,0] = 1.0

original_chain = original_chain[original_chain[:,0]>0,:]

if number_chains_output > 1:
  n_per_file = int(len(original_chain)/number_chains_output)
 
  for i in range(number_chains_output-1):
    np.savetxt(input_poly_chain + "_THIN." + str(i+1) + ".txt", original_chain[i*n_per_file:(i+1)*n_per_file])

  np.savetxt(input_poly_chain + "_THIN." + str(number_chains_output) + ".txt", original_chain[(number_chains_output-1)*n_per_file:])
else:
  np.savetxt(output_file_name, original_chain)

# --------------------------------------------------------------
# TECHNICAL DETAIL - GETDIST NEED THE FIRST COMMENT LINE 
# ---------------------------------------------------------------

try :
  infile = open(input_file_name, "r")
except IOError :
  print("File not open")

try :
  lines = infile.readlines()
except :
  print("Failed Reading Lines of the File")
  exit(1)
finally : 
  infile.close()

# Need to write lines[0] at the top of each file
# Following https://stackoverflow.com/a/5917395
if number_chains_output > 1:
  for i in range(number_chains_output):
    with open(input_poly_chain + "_THIN." + str(i+1) + ".txt", 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(lines[0].rstrip('\r\n') + '\n' + content)
else:
  with open(output_file_name, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(lines[0].rstrip('\r\n') + '\n' + content)


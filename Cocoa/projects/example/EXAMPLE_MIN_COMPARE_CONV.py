#!/usr/bin/env python3

"""Plot how the minimum of EXAMPLE_EMUL_MINIMIZE1.py converges with its budget.

The convergence test runs EXAMPLE_EMUL_MINIMIZE1.py 25 times with growing
budgets; run i writes chains/EXAMPLE_EMUL_MIN_TEST_CONV<i>.txt, i = 0, ..., 24
(scripts/EXAMPLE_EMUL_MIN1_TEST_CONV.sbatch launches them as a job array). This
script reads the minimized chi2 of each run, the last column of its one-row
file, and plots |chi2_min(i) - chi2_min(24)|, the distance to the run with the
largest budget, against n_STW, the emcee steps per walker at each temperature,
with a log y axis from 1e-5 to 10. Where the curve flattens, more steps no
longer move the minimum.

The x values are computed, not read: (5000 + 2500 i)/(5*21) assumes that run i
spent 5000 + 2500 i likelihood evaluations over 5 temperatures and 21 walkers
(three per sampled parameter). The committed sbatch script passes --nstw
100 + 25 i instead, which does not follow this formula; the header line
"nswt (evals/Temp/walker)=..." of each run file records the value actually used.

Inputs: the 25 run files above, found through the ROOTDIR environment variable
(the Cocoa folder, exported by start_cocoa.sh).
Output: chains/EXAMPLE_PLOT_MIN_COMPARE_CONV.pdf in the project folder.

Run from the Cocoa folder after "source start_cocoa.sh":

    python ./projects/example/scripts/EXAMPLE_MIN_COMPARE_CONV.py
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

# Figure style: STIX fonts for text and mathematics, a light grid, and PDF output
# cropped to the drawn area.
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'

# Prefix of the convergence-test files: run i wrote <rt><i>.txt.
rt = os.environ['ROOTDIR']+"/projects/example/chains/EXAMPLE_EMUL_MIN_TEST_CONV"
plt.figure(figsize=(8, 5))
# Style lists; the one curve uses colors[0], markers[5] and linestyles[5]. A tuple
# such as (0,(3,1,1,1)) is a dash pattern: offset 0, then 3 points drawn, 1
# blank, 1 drawn, 1 blank.
colors = ['royalblue','lightcoral','black', 'purple']
markers = ['o', 's', '^', 'v', 'D', '*', 'x', 'P', '<', '>']
linestyles = ['solid',
              '-', 
              '--', 
              '-.', 
              ':', 
              (0,(3,1,1,1)), 
              (0,(5,2)), 
              (0,(1,1)), 
              (0,(3,5,1,5)), 
              (0,(1,10)), 
              (0,(5,1))]
# sz = number of convergence runs.
sz=25
# data [sz, 2]: the list comprehension builds one row per run i, holding the
# n_STW of the docstring formula and the minimized chi2 of the run (the last
# entry of the one-row file, which np.loadtxt returns as a 1D array).
data = np.array([[(5000+2500*i)/(5*21),np.loadtxt(f"{rt}{i}.txt")[-1]] for i in range(sz)])

# Distance of each minimum to the minimum of the largest-budget run, data[-1,1].
plt.plot(data[:,0], 
         abs(data[:,1]-data[-1,1]), 
         marker=markers[5],
         linestyle=linestyles[5],
         color=colors[0], 
         label="$\\Lambda$CDM, CMB+ACTL+DESI-Y3+DES-SN")

# Faint minor grid lines make the log scale easier to read.
ax = plt.gca()
ax.grid(True)
ax.grid(True, 
        which='minor', 
        color='grey', 
        linestyle='--', 
        linewidth=0.25, 
        alpha=0.1)
ax.minorticks_on()
ax.tick_params(axis='both', which='major',labelsize=15)
ax.tick_params(axis='both', which='minor',labelsize=15)
plt.yscale('log')
plt.ylim(1e-5, 10)
# Axis labels: n_STW and Delta chi2_min.
plt.xlabel("$n_{\\rm STW}$")
plt.ylabel("$\\Delta \\chi_{\\rm min}^2$")
ax.legend(fontsize=13, frameon=False)
plt.savefig(os.environ['ROOTDIR']+
            "/projects/example/chains/EXAMPLE_PLOT_MIN_COMPARE_CONV.pdf", 
            bbox_inches='tight')
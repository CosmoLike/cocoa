#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

# GENERAL PLOT OPTIONS
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

rt = os.environ['ROOTDIR']+"/projects/example/chains/EXAMPLE_EMUL_MIN_TEST_CONV"
plt.figure(figsize=(8, 5))
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
sz=25
data = np.array([[(5000+2500*i)/(5*21),np.loadtxt(f"{rt}{i}.txt")[-1]] for i in range(sz)])

plt.plot(data[:,0], 
         abs(data[:,1]-data[-1,1]), 
         marker=markers[5],
         linestyle=linestyles[5],
         color=colors[0], 
         label="$\\Lambda$CDM, CMB+ACTL+DESI-Y3+DES-SN")

# Get current axes
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
# Styling
plt.xlabel("$n_{\\rm STW}$")
plt.ylabel("$\\Delta \\chi_{\\rm min}^2$")
ax.legend(fontsize=13, frameon=False)
plt.savefig(os.environ['ROOTDIR']+
            "/projects/example/chains/EXAMPLE_PLOT_MIN_COMPARE_CONV.pdf", 
            bbox_inches='tight')
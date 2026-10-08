"""Plot the Omega_m, w, wa triangle of the w0wa chain EXAMPLE_EMUL_MCMC3.

EXAMPLE_EMUL_MCMC3.yaml runs Cobaya's Metropolis-Hastings MCMC on a w0wa model,
in which the dark-energy equation of state is w(a) = w + wa (1 - a): w is its
value today (w_0) and wa its slope in the scale factor (w_a). The YAML samples w
and w + wa and derives wa. The theory is emulated, and the data are Planck 2018
CMB, DES-Y5 supernovae, DESI DR2 BAO and ACT DR6 lensing.

GetDist, the chain-analysis package, draws the triangle plot: the marginalized
posterior of each parameter on the diagonal and the 68% and 95% contours of each
pair of parameters below it.

Input: the chain root ../chains/EXAMPLE_EMUL_MCMC3, a path relative to the
working folder, so the script runs from projects/example/scripts/.
Output: EXAMPLE_PLOT_MCMC.pdf in that folder; g.export() without a file name
names the figure after the script.

Run from the Cocoa folder after "source start_cocoa.sh":

    cd ./projects/example/scripts && python EXAMPLE_PLOT_MCMC.py
"""

import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

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

# Parameters of the triangle plot: Omega_m and the dark-energy parameters w (w_0)
# and wa (w_a).
parameter = [u'omegam', u'w', u'wa']
# Relative to the working folder: run the script from projects/example/scripts/.
chaindir=r'../chains/'

# GetDist settings. smooth_scale 0.35: Gaussian smoothing kernels 0.35 standard
# deviations wide, in one and two dimensions. ignore_rows 0.4: the first 40% of
# each chain is dropped as burn-in. range_confidence: the tail probability
# GetDist uses to choose each axis range.
analysissettings={'smooth_scale_1D':0.35,
                  'smooth_scale_2D':0.35,
                  'ignore_rows': u'0.4',
                  'range_confidence' : u'0.005'}


# GetDist plotter; the g.settings lines set fonts, line widths and the legend
# frame.
g=gplot.getSubplotPlotter(chain_dir=chaindir,analysis_settings=analysissettings,width_inch=4.5)
g.settings.axis_tick_x_rotation=65
g.settings.lw_contour = 1.2
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 13.0
g.settings.legend_fontsize = 13.5
g.settings.alpha_filled_add = 0.85
g.settings.lab_fontsize=15.5
g.legend_labels=False

# param_3d = None: no third parameter colors the points of the 2D panels.
param_3d = None
# One chain is plotted, so only the first entry of each style list is used.
g.triangle_plot(['../chains/EXAMPLE_EMUL_MCMC3'],
parameter,
plot_3d_with_param=param_3d,line_args=[
{'lw': 1.2,'ls': 'solid', 'color':'lightcoral'},
{'lw': 1.2,'ls': '--', 'color':'black'},
{'lw': 1.6,'ls': '-.', 'color': 'maroon'},
{'lw': 1.6,'ls': 'solid', 'color': 'indigo'},
],
contour_colors=['lightcoral','black','maroon','indigo'],
contour_ls=['solid','--','-.','solid'], 
contour_lws=[1.0,1.5,1.5,1.0],
filled=[True,False,False,True],
shaded=False,
legend_labels=[
'w0wa (Planck + DES-Y5 + BAO)',
],
legend_loc=(0.48, 0.80))

# Without a file name, export writes EXAMPLE_PLOT_MCMC.pdf (the script name) in
# the working folder.
g.export()
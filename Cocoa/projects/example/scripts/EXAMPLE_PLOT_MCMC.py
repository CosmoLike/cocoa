import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

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

parameter = [u'omegam', u'w', u'wa']
chaindir=r'../chains/'

analysissettings={'smooth_scale_1D':0.35,
                  'smooth_scale_2D':0.35,
                  'ignore_rows': u'0.4',
                  'range_confidence' : u'0.005'}


#GET DIST PLOT SETUP
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

param_3d = None
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

g.export()
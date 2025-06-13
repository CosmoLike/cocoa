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

parameter = [u'logA', u'ns', u'H0', u'omegabh2',  u'omegach2', u'omegaaxh2', u'logmx', u'tau', u'thetastar', u'chi2__BAO', u'chi2__CMB', u'chi2__SN']
chaindir=os.getcwd()

analysissettings={'smooth_scale_1D':0.35,'smooth_scale_2D':0.35,
'ignore_rows': u'0.4','range_confidence' : u'0.005'}

analysissettings2={'smooth_scale_1D':0.35,'smooth_scale_2D':0.35,
'ignore_rows': u'0.0','range_confidence' : u'0.005'}


root_chains = (
  '$ROOTDIR/projects/example/EXAMPLE_MCMC30'
)

# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + root_chains[0],settings=analysissettings)
p = samples.getParams()
samples.saveAsText('.VM_plot9_TMP1')
# --------------------------------------------------------------------------------

#GET DIST PLOT SETUP
g = gplot.getSubplotPlotter(
  chain_dir=chaindir,
  analysis_settings=analysissettings2,
  width_inch=5.0
)
#g.settings.axis_tick_x_rotation=65
g.settings.lw_contour = 1.2
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 13.5
g.settings.legend_fontsize = 13.0
g.settings.alpha_filled_add = 0.7
g.settings.lab_fontsize=15
g.legend_labels=False

param_3d = None
g.triangle_plot(
  [
   chaindir + '/.VM_plot9_TMP1',
  ],
parameter,
plot_3d_with_param=param_3d,line_args=[
  {'lw': 1.0,'ls': 'solid', 'color':'royalblue'},
  {'lw': 1.2,'ls': 'dashed', 'color':'lightcoral'},
  {'lw': 1.4,'ls': '-.', 'color':'grey'},
  {'lw': 1.6,'ls': 'dotted', 'color':'black'},
  {'lw': 1.0,'ls': 'dashdot', 'color':'purple'}
],
filled=[True,True,False,False,False],
contour_colors=['royalblue','lightcoral','grey','black', 'purple'],
contour_ls=['solid', 'dashed', '-.', 'dotted','dashdot'],
contour_lws=[1.0, 1.2, 1.4, 1.6, 1.0],
legend_labels=[
  'EXAMPLE MCMC 30', 
],
legend_loc=(0.425,0.725),
imax_shaded=0)


g.export()


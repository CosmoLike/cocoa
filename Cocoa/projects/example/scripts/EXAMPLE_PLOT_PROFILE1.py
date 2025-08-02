#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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

fig, ((ax1, ax2, ax3, ax4, ax5, ax11, ax13), (ax6, ax7, ax8, ax9, ax10, ax12, ax14)) = plt.subplots(2,7,sharex=True)
fig.set_figheight(12.1)
fig.set_figwidth(12.1*2*7/4.)
fig.subplots_adjust(hspace=0.05, wspace=0.275)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

data = np.loadtxt('../../scripts/profiles_saved/jul15/example_profile1_944_logmx.txt')

omegam=data[:,12]
H0 = data[:,11]
h2 = (H0/100.0)*(H0/100.0)
omegamh2=omegam*h2
omegax=data[:,8]/h2
omegal=(h2-omegamh2-data[:,8])/h2

h2 = (data[:,-2]/100.0)*(data[:,-2]/100.0)

ax1.plot(data[:,0], H0, 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0,
         label ="$n_{\\rm w}=5$, $n_{\\rm eval}^{\\rm T,w}=22,500$")
ax2.plot(data[:,0], omegam, 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax3.plot(data[:,0], omegax, 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax4.plot(data[:,0], omegal, 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax5.plot(data[:,0], data[:,3], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax11.plot(data[:,0], data[:,2], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax13.plot(data[:,0],data[:,7], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)


ax6.plot(data[:,0],data[:,14], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax7.plot(data[:,0],data[:,15], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax8.plot(data[:,0],data[:,16], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax9.plot(data[:,0], data[:,17], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax10.plot(data[:,0], data[:,18], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax12.plot(data[:,0],data[:,19], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)
ax14.plot(data[:,0],data[:,14]+data[:,15]+data[:,16]+data[:,17]+data[:,18]+data[:,19], 
         marker='D', c = 'black', linestyle=None, markersize=4,
         alpha=1.0,lw=1.0)

# --------------------------------------------------------------------------------
#data2 = np.loadtxt('../../scripts/profiles_saved/jul10/example_lcdm_min1_960.txt')
data2 = np.loadtxt('../../scripts/profiles_saved/jul15/example_lcdm_camb_min1_805.txt')

omegamh2=data2[3]+data2[4]+(0.06*(3.046/3)**0.75)/94.0708
h2 = (data2[7]/100.0)*(data2[7]/100.0)

ax1.axhline(y=data2[7],c="red",linewidth=2.5)
ax2.axhline(y=omegamh2/h2,c="red",linewidth=2.5)
ax3.axhline(y=0.0,c="red",linewidth=1.0)
ax4.axhline(y=1.0-omegamh2/h2,c="red",linewidth=2.5)
ax5.axhline(y=data2[1],c="red",linewidth=2.5)
ax11.axhline(y=data2[0],c="red",linewidth=2.5)
ax13.axhline(y=data2[5],c="red",linewidth=2.5)

ax6.axhline(y=data2[-7],c="red",linewidth=2.5)
ax7.axhline(y=data2[-6],c="red",linewidth=2.5)
ax8.axhline(y=data2[-5],c="red",linewidth=2.5)
ax9.axhline(y=data2[-4],c="red",linewidth=2.5)
ax10.axhline(y=data2[-3],c="red",linewidth=2.5)
ax12.axhline(y=data2[-2],c="red",linewidth=1.0)
ax14.axhline(y=data2[-7]+data2[-6]+data2[-5]+data2[-4]+data2[-3]+data2[-2],c="red",linewidth=2.5)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

for ax in (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12):
    ax.grid(True)                  # optional grid
    ax.grid(
            True,
            which='minor',
            color='black',
            linestyle='--',
            linewidth=0.25,
            alpha=0.1
        )
    ax.minorticks_on()
    ax.tick_params(
        axis='both',       # x and y
        which='major',     # major ticks
        labelsize=15       # <-- change this value
    )
    # minor ticks, if you’ve turned them on
    ax.tick_params(
        axis='both',
        which='minor',
        labelsize=15
    )

ax6.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax7.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax8.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax9.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax10.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax12.set_xlabel('$\\log(m_x)$',fontsize = 19)
ax14.set_xlabel('$\\log(m_x)$',fontsize = 19)

ax1.legend(fontsize=11.5, 
          loc = (0.025,0.1),
          frameon=False,
          borderpad=0.05,
          handlelength=0.5,
          scatteryoffsets=[0],
          markerscale=0.5)

ax1.set_ylabel('$H_{0}$',fontsize = 19)
ax2.set_ylabel('$\\Omega_{\\rm m}$',fontsize = 19)
ax3.set_ylabel('$\\Omega_{\\rm X}$',fontsize = 19)
ax4.set_ylabel('$\\Omega_{\\rm \\Lambda}$',fontsize = 19)
ax5.set_ylabel('$n_{\\rm s}$',fontsize = 19)
ax11.set_ylabel('$\\log(A_{\\rm s})$',fontsize = 19)
ax13.set_ylabel('$\\tau$',fontsize = 19)

ax6.set_ylabel('$\\chi^2_{\\rm CMB,TTTEEE}(l>30)$',fontsize = 19)
ax7.set_ylabel('$\\chi^2_{\\rm CMB,TT}(l<30)$',fontsize = 19)
ax8.set_ylabel('$\\chi^2_{\\rm CMB,EE}(l<30)$',fontsize = 19)

ax9.set_ylabel('$\\chi^2_{\\rm SN}$',fontsize = 19)
ax10.set_ylabel('$\\chi^2_{\\rm BAO}$',fontsize = 19)
ax12.set_ylabel('$\\chi^2_{\\phi \\phi}$',fontsize = 19)
ax14.set_ylabel('$\\chi^2_{like}$',fontsize = 19)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

plt.subplots_adjust(bottom=0.25, left = 0.2)
plt.savefig('plot3.pdf')
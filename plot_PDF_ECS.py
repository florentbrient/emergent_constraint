# -*- coding: utf-8 -*-
"""
Ponderate emergent constraint of ECS
Florent Brient
Created on Feb 07 2019
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

def colorYltoRed(nb):
   # Make color bar
   color1=np.linspace(1, 1, nb)
   color2=np.linspace(0.8, 0.0, nb)
   color3=np.linspace(0.0, 0.0, nb)
   colors=np.vstack((color1,color2,color3))
   black =np.array([0,0,0])# for obs
   colors=np.vstack((colors.conj().transpose(), black)).conj().transpose()
   return colors

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 0))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def openfilestat(file):
  f = open(file, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  names   = [ij[0] for ij in tab[1:]]
  mean    = [float(ij[1]) for ij in tab[1:]]
  std     = [float(ij[2]) for ij in tab[1:]]
  f.close()
  return tab,names,mean,std

file="data_ECS.txt"
tab,names,mean,std=openfilestat(file)

NCMIP = 2
NBECS = len(mean)-NCMIP


# Colors for PDF
colors=colorYltoRed(NBECS) # colors has the length of EC number
# Information for figure
typ = ['-','--','-.']
lw  = 5
fts = 25
xsize = 12.0
ysize = 10.0

x = np.linspace(1., 6., 100) # ECS range
#x = np.linspace(-3., 3, 100) # ECS range
allpdf = np.zeros(len(x))

# Name of figure
namefig="PDF_emergent_constraints"

CMIP5only = False
# all emergent constraints
listec = range(NCMIP,len(mean))
if CMIP5only:
  # CMIP5 emergent constraint
  listec = range(5,len(mean))
  namefig += '_CMIP5'


fig = plt.figure()
ax  = fig.add_subplot(111)

ax.plot(x, 100.*norm.pdf(x,loc=mean[0],scale=std[0]),'k--', lw=lw-1, label=names[0])
#pause
ax.plot(x, 100.*norm.pdf(x,loc=mean[1],scale=std[1]),'k-' , lw=lw-1, label=names[1])
for im in listec:
  pdfh  = 100.*norm.pdf(x,loc=mean[im],scale=std[im])/NBECS
  ax.plot(x, pdfh, lw=lw-2, label=names[im],color=colors[:,im-NCMIP])
  print norm.cdf(x,loc=mean[im],scale=std[im])[::4]
  allpdf += pdfh

ax.plot(x, allpdf,'b-', lw=lw, label='All EC')
ax.legend()

adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Equilibrium Climate Sensitivity (K)'
labely = 'Probability (%)'
ax.set_xlabel(labelx,fontsize=fts)
ax.set_ylabel(labely,fontsize=fts)
plt.xticks(size=fts)
plt.yticks(size=fts)
plt.tight_layout()
fig.set_size_inches(xsize, ysize)
fig.savefig(namefig + '.png')
fig.savefig(namefig + '.pdf')
plt.close()
#plt.show()






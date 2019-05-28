# -*- coding: utf-8 -*-
"""
Ponderate emergent constraint of ECS
Florent Brient
Created on Feb 07 2019
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import scipy.stats as stats
import tools as tl



def randomhist(nbsamp,mu,std,bins):
  samples = np.random.normal(mu, std, nbsamp)
  #histogram, bins2 = np.histogram(samples, bins=bins, density=True)
  #hist2=[float(ij) for ij in histogram]
  #hist3=hist2/np.sum(hist2)
  #bin_centers = 0.5*(bins2[1:] + bins2[:-1])
  return tl.makehist(samples,bins)

pathtxt="../text/"  
file   =pathtxt+"data_ECS.txt"
tab,names,mean,std=tl.opendataECS(file)
mean=np.array(mean)
std =np.array(std)

NCMIP = 2
NBECS = len(mean)-NCMIP


# Colors for PDF
colors=tl.colorYltoRed(NBECS) # colors has the length of EC number
# Information for figure
typ = ['-','--','-.']
lw  = 5
fts = 25
xsize = 12.0
ysize = 10.0

x = np.linspace(1., 6., 25) # ECS range
#x = np.linspace(-3., 3, 100) # ECS range
allpdf = np.zeros(len(x)-1)
nbsamp = 50000 # number of points to generate histogram

# Name of figure
namefig="PDF_emergent_constraints"
# path figure
pathfig="../figures/"

CMIP5only = False
# all emergent constraints
listec = range(NCMIP,len(mean))
if CMIP5only:
  # CMIP5 emergent constraint
  listec = range(5,len(mean))
  namefig += '_CMIP5'


fig = plt.figure()
ax  = fig.add_subplot(111)

bins = x
bin_centers,hist3=randomhist(nbsamp,mean[0],std[0],bins)
ax.plot(bin_centers,100.*hist3,'k--', lw=lw-1, label=names[0])

bin_centers,hist3=randomhist(nbsamp,mean[1],std[1],bins)
ax.plot(bin_centers,100.*hist3,'k-' , lw=lw-1, label=names[1])


#ax.plot(x, 100.*norm.pdf(x,loc=mean[0],scale=std[0]),'k--', lw=lw-1, label=names[0])
#plt.show()
#pause
#ax.plot(x, 100.*norm.pdf(x,loc=mean[1],scale=std[1]),'k-' , lw=lw-1, label=names[1])
for im in listec:
  #pdfh  = 100.*norm.pdf(x,loc=mean[im],scale=std[im])/NBECS
  #ax.plot(x, pdfh, lw=lw-2, label=names[im],color=colors[:,im-NCMIP])

  bin_centers,hist3 = randomhist(nbsamp,mean[im],std[im],bins)
  pdfh              = 100.*hist3/NBECS
  ax.plot(bin_centers,pdfh,'k-', lw=lw-2, label=names[im],color=colors[:,im-NCMIP])

  print pdfh.shape
  allpdf += pdfh


#ax.plot(x, allpdf,'b-', lw=lw, label='All EC')
ax.plot(bin_centers, allpdf,'b-', lw=lw, label='All EC')

# Sum of variances
meanec = np.mean(mean[listec])
stdec  = np.sqrt(np.mean(pow(std[listec],2.)) )
print meanec,stdec
bin_centers,hist3=randomhist(nbsamp,meanec,stdec,bins)
ax.plot(bin_centers, 100.*hist3,'g-', lw=lw, label='Sum of variances')
ax.legend()

# Kernel
#kernel_w  = stats.gaussian_kde(mean[listec],weights=std[listec])
#ECSpost   = kernel_w(bin_centers)
#print ECSpost
#ax.plot(bin_centers, 100.*ECSpost,'y-', lw=lw, label='Kernel')


tl.adjust_spines(ax, ['left', 'bottom'])
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
fig.savefig(pathfig+namefig + '.png')
fig.savefig(pathfig+namefig + '.pdf')
plt.close()
#plt.show()






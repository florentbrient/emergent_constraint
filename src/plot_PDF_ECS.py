# -*- coding: utf-8 -*-
"""
Ponderate emergent constraint of ECS
Florent Brient
Created on Feb 07, 2019
"""
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.stats import norm
import scipy.stats as stats
import tools as tl
mpl.rc('font',family='Helvetica')


# Open ECS values from CMIP models
pathtxt = "../text/"  
file    = pathtxt+"data_ECS.txt"
tab,names,mean,std=tl.opendataECS(file)
mean    = np.array(mean)
std     = np.array(std)

# Number of CMIPs (prior)
NCMIP = 2
# Number of emergent constraints
NBECS = len(mean)-NCMIP

# Information about the figure
#colors=tl.colorYltoRed(NBECS) # colors has the length of EC number
colorstart = [0.94,0.5,0.5]  # pink
colorend   = [0.5,0.94,0.94] # blue
colors     = tl.coloryourself(colorstart,colorend,NBECS)
typ        = ['-','--','-.']
lw         = 5
fts        = 25
xsize      = 12.0
ysize      = 10.0
# Name of figure
namefig="PDF_emergent_constraints"
# path figure
pathfig="../figures/"


x = np.linspace(1., 6., 100) # ECS range
allpdf  = np.zeros(len(x))
allpdf2 = np.ones(len(x))

# All emergent constraints
listec = range(NCMIP,len(mean))

CMIP5only = False
if CMIP5only:
  # CMIP5 emergent constraint
  listec = range(5,len(mean)) # check file
  namefig += '_CMIP5'

# Open figure
fig = plt.figure()
ax  = fig.add_subplot(111)

bins  = x
# CMIP distributions
colorCMIP = ['gray','k']
for ic in range(NCMIP):
  hist3 = norm.pdf(x,loc=mean[ic],scale=std[ic])
  ax.plot(bins,100.*hist3,'-',color=colorCMIP[ic], lw=lw-1, label=names[ic])

clevels = .9
ci_l,ci_u = tl.confidence_intervals(hist3,bins,clevels)
print  'Confidence interval CMIP5 : ',ci_l,ci_u

# Each emergent constraints
for im in listec:
  pdfh  = norm.pdf(x,loc=mean[im],scale=std[im])/NBECS
  # Distribution
  ax.plot(x, 100.*pdfh, lw=lw-2, label=names[im],color=colors[:,im-NCMIP], alpha=0.8)
  # Mode
  ax.plot(mean[im],0, marker='o',markersize=10,mew=0,mfc=colors[:,im-NCMIP],zorder=10,clip_on=False)
  allpdf  += pdfh
  allpdf2 *= NBECS*pdfh

###########
# Hereafter several methods have been tried to joint PDFs of emergent constraints
###########

values = mean[listec]

# Quantify the product of PDF (not used)
meanmean,meanstds,scalingfactor = tl.productPDF(values,std[listec])
prodec            = norm.pdf(x,loc=meanmean,scale=meanstds)

# Sum of variances (not used)
#print "listec : ",values
meanec = np.mean(values)
stdec  = np.sqrt(np.sum(pow(std[listec],2.)*pow(1./len(listec),2.)))
#print meanec,stdec,len(listec),std
#hist3  = norm.pdf(x,loc=meanec,scale=stdec)
#ax.plot(bins, 100.*hist3,'b-', lw=lw, label='Sum ECs')

# Sum of uncorrelated variables (not used)
stdec  = np.sqrt(np.sum(pow(std[listec],2.))/pow(len(listec),2.0))
#print meanec,stdec
#hist3  = norm.pdf(x,loc=meanec,scale=stdec)
#ax.plot(bins, 100.*hist3,'b-', lw=lw, label='Sum ECs')


# Weighted and unweighted distribution (used)
bdw    = 1.0 # kde.factor

# Unweighted kernel distribution
kernel = stats.gaussian_kde(values,bw_method=bdw);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g-', lw=lw, label='Unweighted ECs')
clevels = .9; ci_l,ci_u = tl.confidence_intervals(Z,bins,clevels)
print  'Confidence ECs unweighted: ',ci_l,ci_u,np.mean(values),bins[Z==np.max(Z)][0]
print kernel.covariance

# Weighted kernel distribution
# Test different weights
#weights1       = np.exp(std[listec]*-1.)
#weights1       = weights1/np.nansum(weights1)
#weights2       = np.exp(std[listec]*-1. - np.nanmax(std[listec]*-1.))
#weights2       = weights2/np.nansum(weights2);
weights3       = 1.0/(pow(std[listec],2.))
weights3       = weights3/np.nansum(weights3)
#weights4       = 1.0/(pow(std[listec],2.))/np.sum(1.0/(pow(std[listec],2.)))
#weights4       = weights4/np.nansum(weights4)

# Using 1/(sigma^2) as an relative optimal weight.
weights        = weights3

kernel = stats.gaussian_kde(values,bw_method=bdw,weights=weights);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g--', lw=lw, label='Weighted ECs')
clevels = .9; ci_l,ci_u = tl.confidence_intervals(Z,bins,clevels)
print  'Confidence ECs weighted : ',ci_l,ci_u,np.mean(values),bins[Z==np.max(Z)][0]


# Save figure
plt.legend(frameon=False,bbox_to_anchor=(0.63, 0.97, 0.42, .102),
           fontsize=fts*0.7)
tl.adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Equilibrium climate sensitivity ($^\circ$C)'
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






# -*- coding: utf-8 -*-
"""
Ponderate emergent constraint of ECS
Florent Brient
Created on Feb 07 2019
"""
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.stats import norm
import scipy.stats as stats
import tools as tl

mpl.rc('font',family='Helvetica')


#def randomhist(nbsamp,mu,std,bins):
#  samples = np.random.normal(mu, std, nbsamp)
#  return tl.makehist(samples,bins)

def productPDF(means,stds):
  Nb   = len(means)
  print Nb,(Nb-1)/2.
  meanstds = np.sqrt(1./(np.sum(1./(stds*stds))))
  meanmean = np.sum(means/(stds*stds))*np.power(meanstds,2.0)
  #scalingfactor = (1./np.power(2.*np.pi,(Nb-1)/2.)) \
  scalingfactor = (1./np.power(2.*np.pi,(Nb-1)/2.)) \
                 *(np.sqrt(np.power(meanstds,2.0)/np.prod(stds*stds))) \
                 *np.exp(-0.5*(np.sum(means*means/(stds*stds)) - meanmean*meanmean/np.power(meanstds,2.0)))
  print np.sum(stds)
  print scalingfactor
  return meanmean,meanstds,scalingfactor


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

x = np.linspace(1., 6., 50) # ECS range
allpdf  = np.zeros(len(x))#-1)
allpdf2 = np.ones(len(x))

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
hist3  = norm.pdf(x,loc=mean[0],scale=std[0])
ax.plot(bins,100.*hist3,'k--', lw=lw-1, label=names[0])

hist3  = norm.pdf(x,loc=mean[1],scale=std[1])
ax.plot(bins,100.*hist3,'k-' , lw=lw-1, label=names[1])

bx  = fig.add_subplot(111)
for im in listec:
  pdfh  = norm.pdf(x,loc=mean[im],scale=std[im])/NBECS
  bx.plot(x, 100.*pdfh, lw=lw-2, label=names[im],color=colors[:,im-NCMIP])
  print pdfh.shape
  allpdf  += pdfh
  allpdf2 *= NBECS*pdfh

meanmean,meanstds,scalingfactor = productPDF(mean[listec],std[listec])
prodec            = norm.pdf(x,loc=meanmean,scale=meanstds)#*meanstds #*scalingfactor
#ax.plot(bins, 100.*prodec,'y-', lw=lw, label='Prod EC')
#bx.legend(fontsize=fts/3)

# Wrong PDF
#ax.plot(bin_centers, allpdf,'b-', lw=lw, label='All EC')
print allpdf2
#ax.plot(bins, 10*1000*100*allpdf2,'g-', lw=lw, label='All EC')

# Sum of variances
meanec = np.mean(mean[listec])
stdec  = np.sqrt(np.mean(pow(std[listec],2.)) )
# Sum of uncorrelated variables
stdec  = np.sqrt(np.sum(pow(std[listec],2.))/pow(len(listec),2.0))
print meanec,stdec
hist3  = norm.pdf(x,loc=meanec,scale=stdec)
ax.plot(bins, 100.*hist3,'b-', lw=lw, label='Sum ECs')

values = mean[listec]
kernel = stats.gaussian_kde(values,bw_method=1.0);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g-', lw=lw, label='Kernel ECs')
weights        = np.exp(std[listec]*-1.)
w              = np.exp(std[listec]*-1. - np.nanmax(std[listec]*-1.))
w_model        = w/np.nansum(w);
print std[listec],weights,w_model
kernel = stats.gaussian_kde(values,bw_method=1.0,weights=w_model);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g--', lw=lw, label='Kernel ECs wght')


#ax.legend(fontsize=fts/2)
#plt.legend(bbox_to_anchor=(0.6, 0.70, 0.42, .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0., frameon=False,
#           fontsize=fts*0.7)
plt.legend(frameon=False,
           fontsize=fts*0.7)

tl.adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Equilibrium climate sensitivity (C)'
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






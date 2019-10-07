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
#colors=tl.colorYltoRed(NBECS) # colors has the length of EC number
colorstart = [0.94,0.5,0.5]
colorend   = [0.5,0.94,0.94]
colors     = tl.coloryourself(colorstart,colorend,NBECS)
# Information for figure
typ = ['-','--','-.']
lw  = 5
fts = 25
xsize = 12.0
ysize = 10.0

x = np.linspace(1., 6., 100) # ECS range
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
ax.plot(bins,100.*hist3,'-',color='gray', lw=lw-1, label=names[0])

hist3  = norm.pdf(x,loc=mean[1],scale=std[1])
ax.plot(bins,100.*hist3,'k-' , lw=lw-1, label=names[1])

clevels = .9
ci_l,ci_u = tl.confidence_intervals(hist3,bins,clevels)
print  'Confidence CMIP5 : ',ci_l,ci_u

for im in listec:
  pdfh  = norm.pdf(x,loc=mean[im],scale=std[im])/NBECS
  ax.plot(x, 100.*pdfh, lw=lw-2, label=names[im],color=colors[:,im-NCMIP], alpha=0.8)
  ax.plot(mean[im],0, marker='o',markersize=10,mew=0,mfc=colors[:,im-NCMIP],zorder=10,clip_on=False)
  print pdfh.shape
  allpdf  += pdfh
  allpdf2 *= NBECS*pdfh

meanmean,meanstds,scalingfactor = productPDF(mean[listec],std[listec])
prodec            = norm.pdf(x,loc=meanmean,scale=meanstds)#*meanstds #*scalingfactor

# Wrong PDF
#ax.plot(bin_centers, allpdf,'b-', lw=lw, label='All EC')
#print allpdf2
#ax.plot(bins, 10*1000*100*allpdf2,'g-', lw=lw, label='All EC')

# Sum of variances
print "listec : ",mean[listec]
meanec = np.mean(mean[listec])
stdec  = np.sqrt(np.sum(pow(std[listec],2.)*pow(1./len(listec),2.)))
print meanec,stdec,len(listec),std
#hist3  = norm.pdf(x,loc=meanec,scale=stdec)
#ax.plot(bins, 100.*hist3,'b-', lw=lw, label='Sum ECs')
#ax.legend()

# Sum of uncorrelated variables
#stdec  = np.sqrt(np.sum(pow(std[listec],2.))/pow(len(listec),2.0))
#print meanec,stdec
#hist3  = norm.pdf(x,loc=meanec,scale=stdec)
#ax.plot(bins, 100.*hist3,'b-', lw=lw, label='Sum ECs')

bdw    = 1.0
values = mean[listec]

# Unweighted kernel distribution
kernel = stats.gaussian_kde(values,bw_method=bdw);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g-', lw=lw, label='Unweighted ECs')
clevels = .9; ci_l,ci_u = tl.confidence_intervals(Z,bins,clevels)
print  'Confidence ECs : ',ci_l,ci_u,np.mean(values),bins[Z==np.max(Z)]
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

#print std[listec],weights1,weights2,weights3,weights4

# Using 1/(sigma^2) as an relative optimal weight.
weights        = weights3

#average        = np.average(values,weights=weights)
#print np.average(values),average
#print np.average((values-np.average(values))**2),np.average((values-average)**2, weights=weights)
#print np.sum((values-np.mean(values))**2)/pow(len(listec),1.0)

kernel = stats.gaussian_kde(values,bw_method=bdw,weights=weights);Z=kernel(bins)#;plt.plot(bins,Z);plt.show()
ax.plot(bins, 100.*Z,'g--', lw=lw, label='Weighted ECs')
clevels = .9; ci_l,ci_u = tl.confidence_intervals(Z,bins,clevels)
print  'Confidence ECs weighted : ',ci_l,ci_u,np.mean(values),bins[Z==np.max(Z)]


#ax.legend(fontsize=fts/2)
#plt.legend(bbox_to_anchor=(0.6, 0.70, 0.42, .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0., frameon=False,
#           fontsize=fts*0.7)
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
#plt.show()






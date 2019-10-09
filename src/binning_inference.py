# -*- coding: utf-8 -*-
"""
Binning inference from random relationship
Florent Brient
Created on April 17, 2019
"""
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import tools as tl
mpl.rc('font',family='Helvetica')

# Path of statistics files
pathtxt  = "../text/"
fileopen = pathtxt+"statistics_r2.0.10000.txt" # paper

# stats : slope,corr coef
# prior/post : mean,low66,high66,low90,high90
# Post1 : weighting
# Post2 : slope
stats,prior,post1,post2 = tl.openfilestat(fileopen)
print 'Number of emergent constraints : ', stats.shape[1]

# Binning range
bins_corr = np.linspace(0.5, 1., 50)
dbin_corr = bins_corr[1]-bins_corr[0]
bins_corr2, histcorr = tl.makehist(stats[1,:],bins_corr)

# Pdf of bins
pdf       = 100.*dbin_corr*histcorr
plt.plot(bins_corr2,pdf,'k-' , lw=5) #;plt.show()

# Plot prior and posterior distributions
bins_modes = np.linspace(2., 6., 30)
bins2, meanprior = tl.makehist(prior[0,:],bins_modes)
bins2, meanpost1 = tl.makehist(post1[0,:],bins_modes)
bins2, meanpost2 = tl.makehist(post2[0,:],bins_modes)
plt.plot(bins2,100.*meanprior,'k-',lw=5)
plt.plot(bins2,100.*meanpost1,'r-',lw=5)
plt.plot(bins2,100.*meanpost2,'b-',lw=5)
#plt.show()

bin_means = np.zeros((3,len(bins_corr)))*np.nan
bin_std   = np.zeros((3,len(bins_corr)))*np.nan
bin_low66 = np.zeros((3,len(bins_corr)))*np.nan
bin_high66= np.zeros((3,len(bins_corr)))*np.nan

# Values within bins
data      = stats[1,:]
digitized = np.digitize(data, bins_corr)

modes     = np.array([prior[0,:],post1[0,:],post2[0,:]])
low66     = np.array([prior[1,:]-prior[0,:],post1[1,:]-post1[0,:],post2[1,:]-post2[0,:]])
high66    = np.array([prior[2,:]-prior[0,:],post1[2,:]-post1[0,:],post2[2,:]-post2[0,:]])

# Fill the bins
for j in range(modes.shape[0]):
  for i in range(1, len(bins_corr)):
    if pdf[i-1]>0.5: # more than 0.5% in the bin
      #print 'histcorr[i-1]',histcorr[i-1]
      bin_means[j,i] = modes[j,digitized == i].mean()
      bin_std[j,i]   = modes[j,digitized == i].std()
      bin_low66[j,i] = low66[j,digitized == i].mean()
      bin_high66[j,i]= high66[j,digitized == i].mean()

# Plot figure
fig, ax = plt.subplots(figsize=(8, 6))
colors  = ['k','r','b']
labels  = ['Prior','Post (weight)','Post (slope)']
for ij in range(len(colors)):
  tl.plot_modes(bins_corr, bin_means[ij,:], bin_std[ij,:], color=colors[ij], ax=ax,alpha=0.5)
  ax.plot(bins_corr, bin_means[ij,:]+bin_high66[ij,:],color=colors[ij],linestyle='--',lw=1,label=labels[ij])
  ax.plot(bins_corr, bin_means[ij,:]+bin_low66[ij,:] ,color=colors[ij],linestyle='--',lw=1)

yaxis  = [2.0,6.0]
diffy  = (max(yaxis)-min(yaxis))/4.0
ypos   =  min(yaxis)#-diffy
#print ypos
ax.plot(bins_corr2,0.2*histcorr/max(histcorr)+ypos,'k',lw=2)

# Name of figure
namefig='modes_prior_post'
# path of the figures
pathfig="../figures/"
# Size of labels
fts = 25

# Informations figure
plt.legend(loc='upper left',frameon=False,fontsize=fts*0.7)
tl.adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Correlation coefficient'
labely = 'Modes'
xsize  = 12.0
ysize  = 10.0
xaxis  = [0.65,0.95] #boundary x-axis
ax.set_xlim(xaxis)
ax.set_ylim(yaxis)
ax.set_xlabel(labelx,fontsize=fts)
ax.set_ylabel(labely,fontsize=fts)
plt.xticks(size=fts)
plt.yticks(size=fts)
fig.set_size_inches(xsize, ysize)
fig.savefig(pathfig+namefig + '.png')
fig.savefig(pathfig+namefig + '.pdf')
plt.close()


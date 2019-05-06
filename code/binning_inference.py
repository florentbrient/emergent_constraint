# -*- coding: utf-8 -*-
"""
Binning inference from random relationship
Florent Brient
Created on Avril 17, 2019
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import os,sys
sys.path.append('/home/brientf/Documents/Articles/Emergent_Constraint/scipy-1.2.1/')
import scipy as sp
import scipy.stats as stats

def format1(value):
    return "%2.1f" % value
def format2(value):
    return "%4.2f" % value

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


def plot_modes(x2, y2, ci, color="#b9cfe7", ax=None, alpha=.5):
    if ax is None:
        ax = plt.gca()
    ax.plot(x2, y2, color=color, lw=5)
    ax.fill_between(x2, y2+ci, y2-ci, color=color, edgecolor="",alpha=alpha)
    return ax

def makehist(data,bins):
  histogram, bins2 = np.histogram(data, bins=bins, density=True)
  bin_centers = 0.5*(bins2[1:] + bins2[:-1])
  return bin_centers,histogram

def openfilestat(fileout):
  f       = open(fileout, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  n       = 3
  print tab[0],tab[1],tab[3]
  NB      = int(tab[n][-1])
  stats   = np.zeros((2,NB))
  prior   = np.zeros((5,NB))
  post1   = np.zeros((5,NB))
  post2   = np.zeros((5,NB))
  for ij in range(NB):
    print ij
    step      = n+(ij*5)+1
    print step,tab[step]
    line      = tab[step][0].split(',')[1:]
    print line
    stats[:,ij] = [float(line[ii]) for ii in range(2)]
    line      = tab[step+1][0].split(',')[1:]
    prior[:,ij] = [float(line[ii]) for ii in range(5)]
    line      = tab[step+2][0].split(',')[1:]
    post1[:,ij] = [float(line[ii]) for ii in range(5)]
    line      = tab[step+3][0].split(',')[1:]
    post2[:,ij] = [float(line[ii]) for ii in range(5)]
    print stats[:,ij]
  f.close()
  return stats,prior,post1,post2

pathtxt  = "../text/"
fileopen = pathtxt+"statistics_r3.0.10000.0.33.txt"
# stats : slope,corr coef
# prior/post : mean,low66,high66,low90,high90
stats,prior,post1,post2 = openfilestat(fileopen)
print type(stats)

bins_corr = np.linspace(0.5, 1., 50)
dbin_corr = bins_corr[1]-bins_corr[0]
bins_corr2, histcorr = makehist(stats[1,:],bins_corr)
print histcorr.shape
pdf       = 100.*dbin_corr*histcorr
plt.plot(bins_corr2,pdf,'k-' , lw=5);plt.show()

bins_modes = np.linspace(2., 6., 30)
bins2, meanprior = makehist(prior[0,:],bins_modes)
bins2, meanpost1 = makehist(post1[0,:],bins_modes)
bins2, meanpost2 = makehist(post2[0,:],bins_modes)
plt.plot(bins2,100.*meanprior,'k-',lw=5)
plt.plot(bins2,100.*meanpost1,'r-',lw=5)
plt.plot(bins2,100.*meanpost2,'b-',lw=5);plt.show()

data      = stats[1,:]
digitized = np.digitize(data, bins_corr)
#print digitized
bin_means = np.zeros((3,len(bins_corr)))*np.nan
bin_std   = np.zeros((3,len(bins_corr)))*np.nan
bin_low66 = np.zeros((3,len(bins_corr)))*np.nan
bin_high66= np.zeros((3,len(bins_corr)))*np.nan

modes     = np.array([prior[0,:],post1[0,:],post2[0,:]])
low66     = np.array([prior[1,:]-prior[0,:],post1[1,:]-post1[0,:],post2[1,:]-post2[0,:]])
high66    = np.array([prior[2,:]-prior[0,:],post1[2,:]-post1[0,:],post2[2,:]-post2[0,:]])
print modes.shape
for j in range(modes.shape[0]):
  for i in range(1, len(bins_corr)):
    if pdf[i-1]>0.5:
      print 'histcorr[i-1]',histcorr[i-1]
      bin_means[j,i] = modes[j,digitized == i].mean()
      bin_std[j,i]   = modes[j,digitized == i].std()
      bin_low66[j,i] = low66[j,digitized == i].mean()
      bin_high66[j,i]= high66[j,digitized == i].mean()

fig, ax = plt.subplots(figsize=(8, 6))
colors  = ['k','r','b']
for ij in range(len(colors)):
  plot_modes(bins_corr, bin_means[ij,:], bin_std[ij,:], color=colors[ij], ax=ax,alpha=0.5)
  print bin_high66[ij,:]
  ax.plot(bins_corr, bin_means[ij,:]+bin_high66[ij,:],color=colors[ij],linestyle='--',lw=1)
  ax.plot(bins_corr, bin_means[ij,:]+bin_low66[ij,:] ,color=colors[ij],linestyle='--',lw=1)
  #ax.fill_between(bins_corr, bin_means[ij,:]+bin_high66[ij,:], bin_means[ij,:]+bin_low66[ij,:], color=colors[ij], edgecolor="",alpha=0.2)

yaxis  = [2.0,6.0]
diffy  = (max(yaxis)-min(yaxis))/4.0
ypos   =  min(yaxis)#-diffy
print ypos
ax.plot(bins_corr2,0.2*histcorr/max(histcorr)+ypos,'k',lw=2)

# Name of figure
namefig='modes_prior_post'
# path figure
pathfig="../figures/"

adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Correlation Coefficient'
labely = 'Modes'
fts    = 25
xsize  = 12.0
ysize  = 10.0
xaxis  = [0.65,0.95] #[min(bins_corr),max(bins_corr)])
ax.set_xlim(xaxis)
ax.set_ylim(yaxis)
ax.set_xlabel(labelx,fontsize=fts)
ax.set_ylabel(labely,fontsize=fts)
#plt.xticks([0,1,2,3],size=fts)
plt.xticks(size=fts)
plt.yticks(size=fts)
#plt.tight_layout()
fig.set_size_inches(xsize, ysize)
fig.savefig(pathfig+namefig + '.png')
fig.savefig(pathfig+namefig + '.pdf')
plt.close()

#plt.plot(bin_means[0,:],'k');plt.plot(bin_means[0,:]-bin_std[0,:],'k');plt.plot(bin_means[0,:]+bin_std[0,:],'k')
#plt.plot(bin_means[1,:],'r');plt.plot(bin_means[1,:]-bin_std[1,:],'r');plt.plot(bin_means[1,:]+bin_std[1,:],'r')
#plt.plot(bin_means[2,:],'b');plt.plot(bin_means[2,:]-bin_std[2,:],'b');plt.plot(bin_means[2,:]+bin_std[2,:],'b')
#plt.show()



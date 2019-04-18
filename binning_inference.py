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


fileopen = 'statistics_save.txt'
# stats : slope,corr coef
# prior/post : mean,low66,high66,low90,high90
stats,prior,post1,post2 = openfilestat(fileopen)

bins = np.linspace(0.6, 1., 20)
bins2, histogram = makehist(stats[1,:],bins)
#plt.plot(bins2[:-1],100.*histogram,'k-' , lw=5);plt.show()

bins = np.linspace(2., 6., 30)
bins2, meanprior = makehist(prior[0,:],bins)
bins2, meanpost1 = makehist(post1[0,:],bins)
bins2, meanpost2 = makehist(post2[0,:],bins)
plt.plot(bins2,100.*meanprior,'k-',lw=5)
plt.plot(bins2,100.*meanpost1,'r-',lw=5)
plt.plot(bins2,100.*meanpost2,'b-',lw=5);plt.show()

bins      = np.linspace(0.6, 1., 20)
data      = stats[1,:]
digitized = np.digitize(data, bins)
print digitized
data2     = post2[0,:]
bin_means = [data2[digitized == i].mean() for i in range(1, len(bins))]
plt.plot(bin_means);plt.show()



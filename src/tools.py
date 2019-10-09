# -*- coding: utf-8 -*-
"""
Tools for emergent constraints and inference
Florent Brient
Created on May 28, 2019
"""

import numpy as np
from matplotlib import pyplot as plt
import os,sys
sys.path.append('/home/brientf/Documents/Articles/Emergent_Constraint/scipy-1.2.1/')
import scipy as sp
import scipy.stats as stats

# Make histogram
def makehist(data,bins):
  histogram, bins2 = np.histogram(data, bins=bins, density=True)
  bin_centers = 0.5*(bins2[1:] + bins2[:-1])
  return bin_centers,histogram

def openBS16(fileout):
  f       = open(fileout, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  tab     = np.array(tab); NM=len(tab)
  yall      = [float(tab[xx][4]) for xx in range(NM)]
  xall      = [float(tab[xx][8]) for xx in range(NM)]
  sigma_mod = [float(tab[xx][10]) for xx in range(NM)]
  print yall,xall,sigma_mod
  return yall,xall,sigma_mod

# Open ASCII file
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
    #print ij
    step      = n+(ij*5)+1
    #print step,tab[step]
    line      = tab[step][0].split(',')[1:]
    #print line
    stats[:,ij] = [float(line[ii]) for ii in range(2)]
    line      = tab[step+1][0].split(',')[1:]
    prior[:,ij] = [float(line[ii]) for ii in range(5)]
    line      = tab[step+2][0].split(',')[1:]
    post1[:,ij] = [float(line[ii]) for ii in range(5)]
    line      = tab[step+3][0].split(',')[1:]
    post2[:,ij] = [float(line[ii]) for ii in range(5)]
    print 'slope, corrcoef : ',stats[:,ij]
  f.close()
  return stats,prior,post1,post2


def opendataECS(file):
  f = open(file, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  names   = [ij[0] for ij in tab[1:]]
  mean    = [float(ij[1]) for ij in tab[1:]]
  std     = [float(ij[2]) for ij in tab[1:]]
  f.close()
  return tab,names,mean,std

# Create PDF distribution
def productPDF(means,stds):
  Nb   = len(means)
  #print Nb,(Nb-1)/2.
  meanstds = np.sqrt(1./(np.sum(1./(stds*stds))))
  meanmean = np.sum(means/(stds*stds))*np.power(meanstds,2.0)
  scalingfactor = (1./np.power(2.*np.pi,(Nb-1)/2.)) \
                 *(np.sqrt(np.power(meanstds,2.0)/np.prod(stds*stds))) \
                 *np.exp(-0.5*(np.sum(means*means/(stds*stds)) - meanmean*meanmean/np.power(meanstds,2.0)))
  #print np.sum(stds)
  #print scalingfactor
  return meanmean,meanstds,scalingfactor


########################
# Confidence intervals #
######################## 
def confidence_intervals(dpdf,x,clevels):
#CONFIDENCE_INTERVALS Confidence intervals from given pdf
#   calculates the lower bounds ci_l and upper bounds ci_u for the confidence intervals with
#   confidence levels clevels (clevels can be a vector). Inputs are a
#   discrete pdf (dpdf) given at points x.
#   return ci_l, ci_u
 
  dcdf           = sp.integrate.cumtrapz(dpdf,x,initial=0);
 
  # posterior confidence intervals
  lc             = (1-clevels)/2;     # lower bound probabilities
  uc             = lc + clevels;       # upper bound probabilities
  ci_l           = np.interp(lc, dcdf, x) #, 'pchip');
  ci_u           = np.interp(uc, dcdf, x) #, 'pchip');
  return ci_l,ci_u 


def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b)

###################### 
#    For figures     # 
###################### 

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

def colorYltoRed(nb):
   # Make color bar
   color1=np.linspace(1, 1, nb)
   color2=np.linspace(0.8, 0.0, nb)
   color3=np.linspace(0.0, 0.0, nb)
   colors=np.vstack((color1,color2,color3))
   black =np.array([0,0,0])# for obs
   colors=np.vstack((colors.conj().transpose(), black)).conj().transpose()
   return colors

def coloryourself(start,end,nb):
   color1=np.linspace(start[0],end[0],nb)
   color2=np.linspace(start[1],end[1],nb)
   color3=np.linspace(start[2],end[2],nb)
   colors=np.vstack((color1,color2,color3))
   return colors

def plot_modes(x2, y2, ci, color="#b9cfe7", ax=None, alpha=.5):
    if ax is None:
        ax = plt.gca()
    ax.plot(x2, y2, color=color, lw=5)
    ax.fill_between(x2, y2+ci, y2-ci, color=color, edgecolor="",alpha=alpha)
    return ax


def plot_ci_manual(t, s_err, n, x, x2, y2, ax=None):
    """Return an axes of confidence bands using a simple approach.

    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}

    References
    ----------
    .. [1]: M. Duarte.  "Curve fitting," JUpyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb

    """
    if ax is None:
        ax = plt.gca()
    ci = t*s_err*np.sqrt(1/n + (x2-np.mean(x))**2/np.sum((x-np.mean(x))**2))
    ax.fill_between(x2, y2+ci, y2-ci, color="#b9cfe7", edgecolor="",alpha=.5)
    return ax


def plot_ci_bootstrap(xs, ys, resid, nboot=2000, ax=None):
    """Return an axes of confidence bands using a bootstrap approach.

    Notes
    -----
    The bootstrap approach iteratively resampling residuals.
    It plots `nboot` number of straight lines and outlines the shape of a band.
    The density of overlapping lines indicates improved confidence.

    Returns
    -------
    ax : axes
        - Cluster of lines
        - Upper and Lower bounds (high and low) (optional)  Note: sensitive to outliers

    References
    ----------
    .. [1] J. Stults. "Visualizing Confidence Intervals", Various Consequences.
       http://www.variousconsequences.com/2010/02/visualizing-confidence-intervals.html

    """ 
    if ax is None:
        ax = plt.gca()
    bootindex = sp.random.randint
    for _ in range(nboot):
        resamp_resid = resid[bootindex(0, len(resid)-1, len(resid))]
        # Make coeffs of for polys
        pc = sp.polyfit(xs, ys + resamp_resid, 1)                   
        # Plot bootstrap cluster
        ax.plot(xs, sp.polyval(pc, xs), "b-", linewidth=2, alpha=3.0/float(nboot))
    return ax







# -*- coding: utf-8 -*-
"""
Create artificial emergent constraint and inference
Florent Brient
Created on Feb 12 2019
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import os,sys
sys.path.append('/home/brientf/Documents/Articles/Emergent_Constraint/scipy-1.2.1/')
import scipy as sp
import scipy.stats as stats


def confidence_intervals(dpdf,x,clevels):
#CONFIDENCE_INTERVALS Confidence intervals from given pdf
#   calculates the lower bounds ci_l and upper bounds ci_u for the confidence intervals with
#   confidence levels clevels (clevels can be a vector). Inputs are a
#   discrete pdf (dpdf) given at points x.
#   return ci_l, ci_u = 
 
  dcdf           = sp.integrate.cumtrapz(dpdf,x,initial=0);
 
  # posterior confidence intervals
  lc             = (1-clevels)/2;     # lower bound probabilities
  uc             = lc + clevels;       # upper bound probabilities
  ci_l           = np.interp(lc, dcdf, x) #, 'pchip');
  ci_u           = np.interp(uc, dcdf, x) #, 'pchip');
  return ci_l,ci_u 


# Modeling with Numpy
def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b)


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
    ax.fill_between(x2, y2+ci, y2-ci, color="#b9cfe7", edgecolor="")
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



# Number of models
NB=29
# Predictor min,max
MM=[0,3]
# Predictand min,max
ECS=[1.5,4.5]

# randomness of slope
data=np.random.random(NB)
#data=MM[0]+data*(MM[1]-MM[0])
data=data-np.mean(data)

# strength randomness
rdm=2.0

# slope y=ax+b
aa = (ECS[1]-ECS[0])/(MM[1]-MM[0])
bb = ECS[0]

# Perfect slope
# Equidistant between models
xe = np.linspace(min(MM), max(MM), NB)
# Random distance between models
xx = np.random.random(NB) # between 0 and 1
xx = min(MM)+xx*(max(MM)-min(MM))

y0 = aa*xx+bb
print xx,y0

# Unperfect slope
y1 = y0+data*rdm

# Correlation coefficient
corr=np.corrcoef(xx,y1)[0,1]

# Confidence interval slope
p, cov  = np.polyfit(xx, y1, 1, cov=True) 
y_model = equation(p, xx)
   # Statistics
n       = y1.size                                           # number of observations
m       = p.size                                            # number of parameters
dof     = n - m                                             # degrees of freedom
t       = stats.t.ppf(0.975, n - m)                         # used for CI and PI bands
t       = stats.t.ppf(0.95, n - m)                          # used for CI and PI bands
   # Estimates of Error in Data/Model
resid    = y1 - y_model                           
chi2     = np.sum((resid/y_model)**2)                       # chi-squared; estimates error in data
chi2_red = chi2/(dof)                                       # reduced chi-squared; measures goodness of fit
s_err    = np.sqrt(np.sum(resid**2)/(dof))                  # standard deviation of the error

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(xx,y1)
ax.plot(xx,y_model,"-", color="0.1", linewidth=1.5, alpha=0.5, label="Fit") 
x2 = np.linspace(np.min(xx), np.max(xx), 100)
y2 = equation(p, x2)
plot_ci_manual(t, s_err, n, xx, x2, y2, ax=ax)
#plot_ci_bootstrap(xx, y1, resid, ax=ax, nboot=1000)

# Prediction Interval
pi = t*s_err*np.sqrt(1+1/n+(x2-np.mean(xx))**2/np.sum((xx-np.mean(xx))**2))   
ax.fill_between(x2, y2+pi, y2-pi, color="None", linestyle="--")
ax.plot(x2, y2-pi, "--", color="0.5", label="90% Prediction Limits")
ax.plot(x2, y2+pi, "--", color="0.5")


# Observations
obsmean = 0.66*max(MM)
obssigma= 0.3
#c       = np.random.normal(obsmean, obssigma, len(xx))
#hist,bin_edges = np.histogram(c, 20)
xee     = np.linspace(min(MM), max(MM), 10*NB)
obspdf  = norm.pdf(xee,loc=obsmean,scale=obssigma)
obspdf  = obspdf/np.max(obspdf)

# Inference with confidence interval of the curve
nbboot  = 10000
sigma   = s_err                         # Standard deviation of the error
yinfer  = np.zeros(nbboot)
bootindex = sp.random.randint
for ij in range(nbboot):
  idx = bootindex(0, NB-1, NB)
  #resamp_resid = resid[bootindex(0, len(resid)-1, len(resid))]
  # Make coeffs of for polys
  pc = sp.polyfit(xx[idx], y1[idx], 1) # error in xx?
  yinfer[ij]  = pc[0]*obsmean + pc[1] + sigma*np.random.random(1)[0]  # prediction inference

# Confidence interval of yinfer
yimean  = np.mean(yinfer)
yistd   = np.std(yinfer)

yi66    = [yimean-yistd,yimean+yistd]
yi90    = [yimean-2.0*yistd,yimean+2.0*yistd]
xi      = min(xx)-0.1
plt.scatter(xi,yimean,color='g')
plt.plot([xi,xi],yi66,lw=2,color='g')
plt.plot([xi,xi],yi90,lw=1,color='g')

# Kullback–Leibler divergence
log_llh   = np.zeros(NB)
sigma_mod = np.ones(NB)*obssigma # same sigma for each model (sigma=sigmaobs)
for ij in range(NB):
  model_pdf   = norm.pdf(xee,loc=xx[ij],scale=sigma_mod[ij])
  model_pdf   = model_pdf/np.max(model_pdf)
  log_llh[ij] = np.trapz(xee, obspdf * np.log(obspdf / model_pdf))
  print log_llh[ij],xx[ij],y1[ij],sigma_mod[ij]
  #plt.plot(xee, obspdf, 'g', xee, model_pdf, 'r');plt.show()

w              = np.exp(log_llh - np.nanmax(log_llh));
w_model        = w/np.nansum(w);


yee       = np.linspace(min(ECS)-1.5, max(ECS)+1.5, 10*NB)
idx=np.argsort(y1)
kernel    = stats.gaussian_kde(y1[idx])
ECSprior  = kernel(yee)
kernel_w  = stats.gaussian_kde(y1[idx],weights=w_model[idx])
ECSpost   = kernel_w(yee)

priormax           = yee[ECSprior==max(ECSprior)]
priorl90,prioru90  = confidence_intervals(ECSprior,yee,.9)
priorl66,prioru66  = confidence_intervals(ECSprior,yee,.66)
postmax            = yee[ECSpost==max(ECSpost)]
postl90,postu90    = confidence_intervals(ECSpost,yee,.9)
postl66,postu66    = confidence_intervals(ECSpost,yee,.66)

# plot Prior
xi      = min(xx)-0.5
plt.scatter(xi,priormax,color='k')
plt.plot([xi,xi],[priorl66,prioru66],lw=2,color='k')
plt.plot([xi,xi],[priorl90,prioru90],lw=1,color='k')
# plot Post
xi      = min(xx)-0.3
plt.scatter(xi,postmax,color='r')
plt.plot([xi,xi],[postl66,postu66],lw=2,color='r')
plt.plot([xi,xi],[postl90,postu90],lw=1,color='r')


plt.plot(xee,obspdf,'g')
plt.plot([obsmean], [0], marker='o', markersize=3, color="green")

namefig='filename'
adjust_spines(ax, ['left', 'bottom'])
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
labelx = 'Predictor (-)'
labely = 'Equilibrium Climate Sensitivity (K)'
fts    = 25
xsize  = 12.0
ysize  = 10.0
ax.set_ylim([0,max(ECS)+1.5])
ax.set_xlabel(labelx,fontsize=fts)
ax.set_ylabel(labely,fontsize=fts)
plt.xticks([0,1,2,3],size=fts)
plt.yticks(size=fts)
#plt.tight_layout()
fig.set_size_inches(xsize, ysize)
fig.savefig(namefig + '.png')
fig.savefig(namefig + '.pdf')
plt.close()















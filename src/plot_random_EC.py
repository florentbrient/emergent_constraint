# -*- coding: utf-8 -*-
"""
Create artificial emergent constraint and calculate inference
Florent Brient
Created on Feb 12 2019
"""

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.stats import norm
import os,sys
sys.path.append('/home/brientf/Documents/Articles/Emergent_Constraint/scipy-1.2.1/')
import scipy as sp
import scipy.stats as stats
<<<<<<< HEAD
import tools as tl
=======
mpl.rc('font',family='Helvetica')
>>>>>>> at_home

def format1(value):
    return "%2.1f" % value
def format2(value):
    return "%4.2f" % value
def format3(value):
    return "%7s" % value


# if makerandom, create artificial relationship
# otherwize, upload data
makerandom=0
# makefigure
makefigure=1

if makerandom:
  # Number of models
  NB=29
  # Predictor min,max
  MM=[0,3]
  # Predictand min,max
  ECS=[2.0,5.0]
  # strength randomness of the relationship
  rdm=2.0 #1.5

  # slope y=ax+b
  aa = (ECS[1]-ECS[0])/(MM[1]-MM[0])
  bb = ECS[0]

  # Perfect slope
  # Equidistant between models
  xe = np.linspace(min(MM), max(MM), NB)
  # Number of random set (default = 1)
  NR = 1 #10000
  # randomness of slope
  data = np.random.random(NB)
  data = data-np.mean(data)
  # Random distance between models
  xall = np.zeros((NR,NB))
  for ii in range(NR):
    xall[ii,:] = np.random.random(NB) # between 0 and 1
    xall[ii,:] = min(MM)+xall[ii,:]*(max(MM)-min(MM))

  y0   = aa*xall+bb
  #print aa,bb#,xx,y0

  # Unperfect slope
  yall = y0+data*rdm 

  # Observations
  obsmean = 0.66*max(MM) #default
  obssigma= 0.3
  #xplot,obspdf = makehist(obsmean,obssigma,xe)

  #xplot  = np.linspace(min(MM), max(MM), 10*NB)
  diff         = (max(MM)-min(MM))/2.0
  xplot        = np.linspace(min(MM)-diff, max(MM)+diff, 10*NB)

  #plt.plot(xplot,100.*obspdf);plt.show()
  # Make observation distribution
  obspdf  = norm.pdf(xplot,loc=obsmean,scale=obssigma)
  obspdf  = obspdf*obssigma
  #obspdf  = obspdf/np.mean(obspdf)


  # Models
  sigma_mod = np.ones(NB)*obssigma # same sigma for each model (sigma=sigmaobs)
  model_pdf_all = np.zeros((NR,NB,len(obspdf)))
  for ii in range(NR):
    for ij in range(NB):
      pdf                    = norm.pdf(xplot,loc=xall[ii,ij],scale=sigma_mod[ij])
      model_pdf_all[ii,ij,:] = pdf*sigma_mod[ij] #/np.mean(pdf)
  
else:
  print 'import your data --  not ready yet'
  diropen  = '../text/'
  fileopen = diropen+'data_cloud_Brient_Schneider.txt'
  yall,xall,sigma_mod     = tl.openBS16(fileopen)
  NR = 1; NB = len(yall)
  yall = np.reshape(yall,(NR,NB))
  xall = np.reshape(xall,(NR,NB))
  #sigma_mod = np.reshape(sigma_mod,(NR,NB))
  #print sigma_mod.shape
  minMM,maxMM   = np.min(xall),np.max(xall)
  diff          = (maxMM-minMM)/2.0
  print 'aa ',minMM,maxMM,diff
  xplot         = np.linspace(minMM-diff, maxMM+diff, 10*NB)
  print len(xplot)
  NX            = len(xplot)
  model_pdf_all = np.zeros((NR,NB,NX))
  for ii in range(NR):
    for ij in range(NB):
      print xall[ii,ij],sigma_mod[ij]
      pdf                    = norm.pdf(xplot,loc=xall[ii,ij],scale=sigma_mod[ij])
      model_pdf_all[ii,ij,:] = pdf*sigma_mod[ij] #/np.mean(pdf)

  # Make observation distribution
  obsmean = -0.96
  obssigma= 0.22
  obspdf  = norm.pdf(xplot,loc=obsmean,scale=obssigma)
  obspdf  = obspdf*obssigma
  #exit(1)

# Calculate correlation coefficient, statistical inference

pathtxt  = "../text/"
filesave = pathtxt+"statistics.txt"
f        = open(filesave, 'wb')
f.write('Statistics for random relationship\n')
f.write('Number of models: '+str(NB)+'\n')
if makerandom:
  f.write('Randomness slope: '+format1(rdm)+'\n')
  f.write('Number of set   : '+str(NR)+'\n')
textformat0 = "stats,{slope},{r2}\n"
textformat1 = "{typ},{mode},{low66},{high66},{low90},{high90}\n"

for ii in range(NR):
  print ii
  xx = xall[ii,:]
  y1 = yall[ii,:]
  model_pdf = model_pdf_all[ii,:,:]

  # Correlation coefficient
  corr=np.corrcoef(xx,y1)[0,1]

  #### Confidence interval slope
  p, cov  = np.polyfit(xx, y1, 1, cov=True) 
  #print 'p ',p
  y_model = tl.equation(p, xx)
    # Statistics
  n       = y1.size                                           # number of observations
  m       = p.size                                            # number of parameters
  dof     = n - m                                             # degrees of freedom
  #t       = stats.t.ppf(0.975, n - m)                         # used for CI and PI bands
  t       = stats.t.ppf(0.95, n - m)                          # used for CI and PI bands
     # Estimates of Error in Data/Model
  resid    = y1 - y_model                           
  chi2     = np.sum((resid/y_model)**2)                       # chi-squared; estimates error in data
  chi2_red = chi2/(dof)                                       # reduced chi-squared; measures goodness of fit
  s_err    = np.sqrt(np.sum(resid**2)/(dof))                  # standard deviation of the error

  # Inference with confidence interval of the curve
  nbboot  = 10000                         # number of bootstrap
  sigma   = s_err                         # Standard deviation of the error
  yinfer  = np.zeros(nbboot)
  bootindex = sp.random.randint
  for ij in range(nbboot):
    idx = bootindex(0, NB-1, NB)
    #resamp_resid = resid[bootindex(0, len(resid)-1, len(resid))]
    # Make coeffs of for polys
    pc = sp.polyfit(xx[idx], y1[idx], 1) # error in xx?
    yinfer[ij]  = pc[0]*obsmean + pc[1] + sigma*np.random.randn()  # prediction inference

  # Confidence interval of yinfer
  yimean  = np.mean(yinfer)
  yistd   = np.std(yinfer)
  yi66    = [yimean-yistd,yimean+yistd]
  yi90    = [yimean-2.0*yistd,yimean+2.0*yistd]

  # Kullback–Leibler divergence
  log_llh   = np.zeros(NB)
  for ij in range(NB):
    #plt.plot(xplot, obspdf, 'g', xplot, model_pdf[ij,:], 'r');plt.show()
    log_llh[ij] = np.trapz(xplot, obspdf * np.log(obspdf / model_pdf[ij,:]))
    #print log_llh[ij],xx[ij],y1[ij],sigma_mod[ij]

  w              = np.exp(log_llh - np.nanmax(log_llh));
  w_model        = w/np.nansum(w);

  yee       = np.linspace(min(y1), max(y1), NB) #np.linspace(min(ECS)-1.5, max(ECS)+1.5, 10*NB)
  idx       = np.argsort(y1)
  kernel    = stats.gaussian_kde(y1[idx])
  ECSprior  = kernel(yee)
  kernel_w  = stats.gaussian_kde(y1[idx],weights=w_model[idx])
  ECSpost   = kernel_w(yee)

  priormax           = yee[ECSprior==max(ECSprior)]
  priorl90,prioru90  = tl.confidence_intervals(ECSprior,yee,.9)
  priorl66,prioru66  = tl.confidence_intervals(ECSprior,yee,.66)
  postmax            = yee[ECSpost==max(ECSpost)]
  postl90,postu90    = tl.confidence_intervals(ECSpost,yee,.9)
  postl66,postu66    = tl.confidence_intervals(ECSpost,yee,.66)

  #textformat = "{typ}:  {mode},{low66},{high66},{low90},{high90}"
  f.write(textformat0.format(slope=format2(p[0]),r2=format2(corr)))
  f.write(textformat1.format(typ='Prior',mode=format2(priormax),low66=format2(priorl66)
   ,high66=format2(prioru66),low90=format2(priorl90),high90=format2(prioru90)))
  f.write(textformat1.format(typ='Post1',mode=format2(postmax),low66=format2(postl66)
   ,high66=format2(postu66),low90=format2(postl90),high90=format2(postu90)))
  f.write(textformat1.format(typ='Post2',mode=format2(yimean),low66=format2(yi66[0])
   ,high66=format2(yi66[1]),low90=format2(yi90[0]),high90=format2(yi90[1])))
  f.write("***\n")
f.close()

if makefigure:

  fig, ax = plt.subplots(figsize=(8, 6))
  ax.scatter(xx,y1,color='k')
  ax.plot(xx,y_model,"-", color="0.1", linewidth=1.5, alpha=0.5, label="Fit") 
  x2 = np.linspace(np.min(xx), np.max(xx), 100)
  y2 = tl.equation(p, x2)
  tl.plot_ci_manual(t, s_err, n, xx, x2, y2, ax=ax)
  #tl.plot_ci_bootstrap(xx, y1, resid, ax=ax, nboot=1000)

  # Prediction Interval
  pi = t*s_err*np.sqrt(1+1/n+(x2-np.mean(xx))**2/np.sum((xx-np.mean(xx))**2))
  #ax.fill_between(x2, y2+pi, y2-pi, color="None", linestyle="--")
  ax.plot(x2, y2-pi, "--", color="0.5", label="90% Prediction Limits")
  ax.plot(x2, y2+pi, "--", color="0.5")

  # plot Prior
  xi      = min(xx)-0.5
  plt.plot([xi,xi],[priorl66,prioru66],lw=3,color='k')
  plt.plot([xi,xi],[priorl90,prioru90],lw=1,color='k')
  plt.plot([xi],priormax, marker='o',markersize=6,color='k',markeredgecolor='None')

  # plot Post
  xi      = min(xx)-0.3
  plt.plot([xi,xi],[postl66,postu66],lw=3,color='r')
  plt.plot([xi,xi],[postl90,postu90],lw=1,color='r')
  plt.plot([xi],postmax, marker='o',markersize=6,color='r',markeredgecolor='None')

  xi      = min(xx)-0.1
  #plt.scatter(xi,yimean,color='g')
  plt.plot([xi,xi],yi66,lw=3,color='b')
  plt.plot([xi,xi],yi90,lw=2,color='b')
  plt.plot([xi],yimean,marker='o',markersize=6,color='b',markeredgecolor='None')

  diffx  = (max(xx)-min(xx))/4.0
  diffy  = (max(y1)-min(y1))/4.0
  ypos   =  min(y1)-diffy
  plt.plot(xplot,obspdf/max(obspdf)+ypos,'g',lw=2)
  plt.plot([obsmean], ypos, marker='o', markersize=8, color="green")

  if makerandom:
    title = 'Relationship preditor/predictand (randomness={rdm:1.1f}$\sigma$,r={corr:02.2f})'.format(rdm=rdm,corr=corr)
  else: 
    title = 'Not random'
  plt.title(title)

  # Name of figure
  namefig='filename'
  # path figure
  pathfig="../figures/"

  tl.adjust_spines(ax, ['left', 'bottom'])
  ax.get_yaxis().set_tick_params(direction='out')
  ax.get_xaxis().set_tick_params(direction='out')
  labelx = 'Predictor A (-)'
  labely = 'Predictand B (-)'
  fts    = 25
  xsize  = 12.0
  ysize  = 10.0
  ax.set_xlim([min(xx)-diffx,max(xx)+diffx])
  ax.set_ylim([ypos,max(y1)+diffy])
  ax.set_xlabel(labelx,fontsize=fts)
  ax.set_ylabel(labely,fontsize=fts)
  plt.xticks([0,1,2,3],size=fts)
  plt.yticks(size=fts)
  #plt.tight_layout()
  fig.set_size_inches(xsize, ysize)
  fig.savefig(pathfig+namefig + '.png')
  fig.savefig(pathfig+namefig + '.pdf')
  plt.close()















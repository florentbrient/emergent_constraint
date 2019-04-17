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
    step      = n+(ij*5)+1
    print step,tab[step]
    line      = tab[step][0].split(',')
    print line
    stats[:,ij] = [float(line[ij+1]) for ij in range(2)]
    line      = tab[step+1][0].split(',')
    prior[:,ij] = [float(line[ij+1]) for ij in range(5)]
  f.close()
  return tab,nametab,Nmin,Nval,ijmin


fileopen = 'statistics.txt'
openfilestat(fileopen)


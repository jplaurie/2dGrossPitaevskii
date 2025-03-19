
import numpy as np
from math import *

from matplotlib import pyplot as plt
#from matplotlib.pyplot import *
#from matplotlib.ticker import MultipleLocator
from matplotlib import animation,cm,rc,colors
#from matplotlib  import ticker as mticker

#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource, LogNorm
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import scipy as scp
from scipy import io
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit

import seaborn as sns


##from time import time as ctime
##import os,glob,subprocess,shutil
##import warnings; warnings.simplefilter('ignore')
##import timeit,itertools

##import sys
#print(sys.executable)
#%%
def newfig(a=1,b=1,figheight=12,aspectratio=0.75,fontsize=34,**kwargs):
    rc('legend', frameon=False,fancybox=False,shadow=False,fontsize=28,loc='best')
    rc('lines', linewidth=1)
    font = {'family':'serif','size':fontsize}
    rc('font',**font)
    rc('text', usetex=True)
    rc('xtick',labelsize=fontsize)
    rc('ytick',labelsize=fontsize)
    rc('savefig',format='pdf')
    return plt.subplots(a,b,figsize=(b*figheight,a*figheight*aspectratio),clear=True,tight_layout=True,**kwargs)


import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np
import os
import sys

N = int(sys.argv[1])
results_b = np.loadtxt('g_vector%i.dat'%N)
results_c = np.loadtxt('s_vector%i.dat'%N)
x = np.linspace(0,1,(10**N)+1)

def u(parameter):
    u = 1 - (1 - np.exp(-10))*parameter - np.exp(-10*parameter)
    return u

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(x,u(x),'g',x,results_b,'c',x,results_c,'r--')
plt.title('Analytic vs numerical solution with N = %i'%10**N,fontsize=18)
plt.legend(['Closed form solution','General algorithm solution','Specialized algorithm solution'])
plt.xlabel(r'\textit{$x_i$}', fontsize=16)
plt.ylabel(r'\textit{$u_i$ and $v_i$ values}',fontsize=16)
plt.show()

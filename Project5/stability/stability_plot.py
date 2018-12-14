import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display

N = 40

results_a = np.loadtxt('stability_euler1.dat')
results_aa = np.loadtxt('stability_leap1.dat')

results_b = np.loadtxt('stability_euler2.dat')
results_bb = np.loadtxt('stability_leap2.dat')

results_c = np.loadtxt('stability_euler3.dat')
results_cc = np.loadtxt('stability_leap3.dat')

results_d = np.loadtxt('stability_euler4.dat')
results_dd = np.loadtxt('stability_leap4.dat')


x = np.linspace(0,1,(N)+1)
def u(parameter):
    u = np.sin(np.pi*4*parameter)
    return u


plt.figure(figsize=(10,10))
plt.suptitle('Comparison between forward-Euler and leapfrog time schemes',fontsize=15)
plt.subplot(221)
plt.plot(x,results_d,'r--',zorder=1)
plt.plot(x,results_dd,zorder=1,alpha = 0.8)
plt.legend(['Forward-Euler scheme','Leapfrog scheme'])
plt.title(r'dt = 0.01')
plt.ylim(-0.40,1.80)
plt.xlabel(r'$x$ (west-east)',fontsize=12)
plt.ylabel(r'Streamfunction $\psi (x,t)$',fontsize=12)
plt.grid(True)
plt.subplot(222)
plt.plot(x,results_a,'r--',zorder=1)
plt.plot(x,results_aa,zorder=1,alpha = 0.8)
plt.legend(['Forward-Euler scheme','Leapfrog scheme'])
plt.title(r'dt = 1.0')
plt.xlabel(r'$x$ (west-east)',fontsize=12)
plt.ylim(-0.40,1.80)
plt.grid(True)
plt.savefig('comparison.pdf',bbox_inches='tight')

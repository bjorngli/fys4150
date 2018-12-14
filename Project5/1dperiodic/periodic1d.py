import matplotlib.pyplot as plt
import numpy as np
import os
import sys

N = 40

results_b = np.loadtxt('periodic1d.dat')

x = np.linspace(0,1,(N)+1)
def u(parameter):
    u = np.sin(np.pi*4*parameter)
    return u

#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('text', usetex=True)

#plt.plot(x,u(x))
#plt.plot(x,u(x),'b',zorder=2,alpha=0.4)
#plt.plot(x,results_b[3,:],'r--',zorder=1)#,x,results_c)
plt.plot(x,results_b,'r--',zorder=1)#,x,results_c)
plt.title('Analytic and numerical solution of the Poisson eq.')
plt.legend(['Closed form solution u(x)','General algorithm solution v(x)'])#,'Simplified algorithm'])
plt.xlabel(r'$x$',fontsize=12)
plt.ylabel(r'$v(x)$' '  '  r'and ' '  '  r'$u(x)$',fontsize=12)
plt.grid(True)
plt.savefig('Num_plot%i.png'%N)
plt.show()

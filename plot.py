import matplotlib.pyplot as plt
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

#plt.plot(x,u(x))
plt.plot(x,u(x),x,results_b,x,results_c)
plt.title('Values of vector v with N = %i'%N)
plt.legend(['Closed form solution','General algorithm solution','Simplified algorithm'])
#plt.xlabel('')
#plt.ylabel('')
plt.show()

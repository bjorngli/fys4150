import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import scipy.optimize as scipy

matrix_dim = np.loadtxt('dimensions.dat')
iterations = np.loadtxt('iterations.dat')

def func(x, a, b):
    return a*x**2 + b

popt,pcov = scipy.curve_fit(func,  matrix_dim,  iterations,  p0=(0, 0))
a,b = popt

plt.plot(matrix_dim,iterations,'b',zorder=2,alpha=0.7)
plt.plot(matrix_dim,func(matrix_dim,*popt),'r--',alpha=1)

plt.title('Iterations needed as a function of matrix dimension')
#plt.legend(['Recorded number of iterations','Fitted curve {:4.0f}*exp({:0.4f}x) {:4.0f}'.format(a,b)])
plt.legend(['Recorded number of iterations','Fitted curve {:4.0f}*x**2 {:4.0f}'.format(a,b)])
plt.xlabel(r'$Matrix dimension$',fontsize=12)
plt.ylabel(r'$Iterations$',fontsize=12)
plt.grid(True)
plt.savefig('dim_vs_iter.png')
plt.show()

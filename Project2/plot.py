import matplotlib.pyplot as plt
import numpy as np
import os
import sys

matrix_dim = np.loadtxt('dimensions.dat')
iterations = np.loadtxt('iterations.dat')


plt.plot(matrix_dim,iterations,'b',zorder=2,alpha=0.4)

plt.title('Iterations needed as a function of matrix dimension')
plt.legend(['Iterations needed as a function of matrix dimension'])
plt.xlabel(r'$Matrix dimension$',fontsize=12)
plt.ylabel(r'$Iterations$',fontsize=12)
plt.grid(True)
plt.savefig('dim_vs_iter.png')
plt.show()

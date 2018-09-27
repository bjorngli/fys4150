import matplotlib.pyplot as plt
import numpy as np
import os
import sys

eigen_values = np.loadtxt('one_electron.dat')
eigen_values = np.sort(eigen_values)

print(eigen_values)
print(eigen_values[0],eigen_values[1],eigen_values[2],eigen_values[3])

# plt.title('Iterations needed as a function of matrix dimension')
# plt.legend(['Iterations needed as a function of matrix dimension'])
# plt.xlabel(r'$Matrix dimension$',fontsize=12)
# plt.ylabel(r'$Iterations$',fontsize=12)
# plt.grid(True)
# plt.savefig('dim_vs_iter.png')
# plt.show()

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Load in data for n = 200 and rho_max = 10
eigen_values1 = np.loadtxt('two_electron1.dat')
eigen_vectors1 = np.loadtxt('two_outfile1.txt')

eigen_values2 = np.loadtxt('two_electron2.dat')
eigen_vectors2 = np.loadtxt('two_outfile2.txt')

eigen_values3 = np.loadtxt('two_electron3.dat')
eigen_vectors3 = np.loadtxt('two_outfile3.txt')

eigen_values4 = np.loadtxt('two_electron4.dat')
eigen_vectors4 = np.loadtxt('two_outfile4.txt')

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute1 = eigen_values1.argsort()
eigen_values1 = eigen_values1[permute1]
eigen_vectors1 = eigen_vectors1[:,permute1]

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute2 = eigen_values2.argsort()
eigen_values2 = eigen_values2[permute2]
eigen_vectors2 = eigen_vectors2[:,permute2]

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute3 = eigen_values3.argsort()
eigen_values3 = eigen_values3[permute3]
eigen_vectors3 = eigen_vectors3[:,permute3]

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute4 = eigen_values4.argsort()
eigen_values4 = eigen_values4[permute4]
eigen_vectors4 = eigen_vectors4[:,permute4]


# Print statements for testing
#print(eigen_vectors[0])
print(eigen_values1[0])
print(eigen_values2[0])
print(eigen_values3[0])
print(eigen_values4[0])
#print(eigen_values[0],eigen_values[1],eigen_values[2],eigen_values[3])


# Plot of the three lowest lying eigenstates
rho_max = 5
n = 200
r = np.linspace(0,rho_max,n)

FirstEigvector1 = eigen_vectors1[:,0]
FirstEigvector2 = eigen_vectors2[:,0]
FirstEigvector3 = eigen_vectors3[:,0]
FirstEigvector4 = eigen_vectors4[:,0]

plt.figure()


#area = np.trapz(FirstEigvector**2, dx=(10/200))
#print(area)
#sum1 = np.sum(FirstEigvector)
#print(sum1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(r, FirstEigvector1**2,label=r'$\omega_r=0.01$')
plt.plot(r, FirstEigvector2**2,label=r'$\omega_r=0.5$')
plt.plot(r, FirstEigvector3**2,label=r'$\omega_r=1$')
plt.plot(r, FirstEigvector4**2,label=r'$\omega_r=5$')
plt.legend()
plt.xlabel(r'$\rho$',fontsize=14)
plt.ylabel(r'Radial probability $|\psi(\rho)|^2$',fontsize=14)
plt.suptitle(r'Radial probability distributions at ground state',fontsize=16)
plt.title(r' with different oscillator frequencies $\omega_r$',fontsize=16)
plt.savefig('eigenvector2.pdf')
plt.show()

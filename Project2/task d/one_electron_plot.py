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
eigen_values = np.loadtxt('one_electron.dat')
eigen_vectors = np.loadtxt('one_outfile.txt')

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute = eigen_values.argsort()
eigen_values = eigen_values[permute]
eigen_vectors = eigen_vectors[:,permute]


# Print statements for testing
#print(eigen_vectors[0])
print(eigen_values[0],eigen_values[1],eigen_values[2],eigen_values[3])


# Plot of the three lowest lying eigenstates
rho_max = 5
n = 200
r = np.linspace(0,rho_max,n)

FirstEigvector = eigen_vectors[:,0]
SecondEigvector = eigen_vectors[:,1]
ThirdEigvector = eigen_vectors[:,2]
plt.figure(1,figsize=[8,8])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.subplot(311)
plt.title(r'Radial probability distribution for the eigenstates with l=0',fontsize=16)
plt.plot(r, FirstEigvector**2 ,'b-',label=r'$\lambda_1=$%7.4f'%eigen_values[0])
plt.xlim((0,5))
plt.xlabel(r'$\rho$',fontsize=14)
plt.ylabel(r'Radial probability $|\psi(\rho)|^2$',fontsize=14)
#plt.grid(True)
plt.legend()
plt.subplot(312)
plt.plot(r, SecondEigvector**2 ,'g-',label=r'$\lambda_2=$%7.4f'%eigen_values[1])
plt.xlim((0,5))
plt.xlabel(r'$\rho$',fontsize=14)
plt.ylabel(r'Radial probability $|\psi(\rho)|^2$',fontsize=14)
#plt.grid(True)
plt.legend()
plt.subplot(313)
plt.plot(r, ThirdEigvector**2 ,'r-',label=r'$\lambda_3=$%7.4f'%eigen_values[2])
plt.xlim((0,5))
plt.xlabel(r'$\rho$',fontsize=14)
plt.ylabel(r'Radial probability $|\psi(\rho)|^2$',fontsize=14)
#plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('eigenvector1.pdf')
plt.show()

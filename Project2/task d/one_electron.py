import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd

# Load in data for n = 200 and rho_max = 10
eigen_values = np.loadtxt('one_electron200.dat')
eigen_vectors = np.loadtxt('outfile200.txt')

# Sort eigen values. And change placement of eigen vectors to correspond to the sorted eigen values.
permute = eigen_values.argsort()
eigen_values = eigen_values[permute]
eigen_vectors = eigen_vectors[:,permute]


# Print statements for testing
#print(eigen_vectors[0])
#print(eigen_values[0],eigen_values[1],eigen_values[2],eigen_values[3])


# Plot of the three lowest lying eigenstates
rho_max = 10
n = 200
r = np.linspace(0,rho_max,n)

FirstEigvector = eigen_vectors[:,0]
SecondEigvector = eigen_vectors[:,1]
ThirdEigvector = eigen_vectors[:,2]
plt.figure(1,figsize=[8,10])

plt.subplot(311)
plt.plot(r, FirstEigvector**2 ,'b-')
plt.axis([0,4.6,0.0, 0.05])
plt.xlabel(r'$r$')
plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
plt.title(r'Radial probability distributions for first lowest-lying states')
plt.grid(True)
plt.subplot(312)
plt.plot(r, SecondEigvector**2 ,'g-')
plt.axis([0,4.6,0.0, 0.05])
plt.xlabel(r'$r$')
plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
plt.title(r'Radial probability distributions for second lowest-lying states')
plt.grid(True)
plt.subplot(313)
plt.plot(r, ThirdEigvector**2 ,'r-')
plt.axis([0,4.6,0.0, 0.05])
plt.xlabel(r'$r$')
plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
plt.title(r'Radial probability distributions for third lowest-lying states')
plt.grid(True)
plt.tight_layout()
plt.savefig('eigenvector1.pdf')
plt.show()

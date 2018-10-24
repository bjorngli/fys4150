import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
energy = pd.read_csv('energy_error.dat',sep='\s+',header=None,names=['error'])
momentum = pd.read_csv('momentum_error.dat',sep='\s+',header=None,names=['error'])
dist = pd.read_csv('dist_error.dat',sep='\s+',header=None,names=['error'])
dist_euler = pd.read_csv('dist_error_euler.dat',sep='\s+',header=None,names=['error'])
error_line = np.repeat(1e-5, len(energy))

#energy_mom_euler = pd.read_csv('euler_energy.dat',sep='\s+',header=None,names=['energy','momentum'])

# #Energy
# plt.plot(np.linspace(10/20,10/(len(energy)+20),len(energy)),energy['error'],label='Energy error')
# plt.plot(np.linspace(10/20,10/(len(energy)+20),len(energy)),error_line,label='Accepted error 1e-5')
# plt.ylabel('Maximum difference from starting energy')
# plt.xlabel('Time step deltaT')
# plt.legend()
# plt.grid(True)
# plt.show()
plt.figure(figsize=(7,6))
#Energy
plt.subplot(121)

plt.plot(np.linspace(20,(len(energy)+20),len(energy)),energy['error'],label='Energy error')
plt.plot(np.linspace(20,(len(energy)+20),len(energy)),error_line,label='Accepted error 1e-5')
plt.ylabel(r'Maximum difference from starting energy',size=12)
plt.xlabel(r'Number of integration points',size=12)
#plt.legend()
plt.grid(True)
plt.savefig('energy.pdf',bbox_inches='tight')
plt.show()

#
#Momentum
plt.subplot(122)
plt.plot(np.linspace(20,(len(momentum)+20),len(momentum)),np.log10(momentum['error']),label='Momentum error')
plt.plot(np.linspace(20,(len(momentum)+20),len(momentum)),np.log10(error_line),label='Accepted error 1e-5')
plt.ylabel(r'$log_{10}$ of max difference from starting angular momentum',size=12)
plt.xlabel(r'Number of integration points',size=12)
#plt.legend()
plt.grid(True)
plt.savefig('momentum.pdf',bbox_inches='tight')
plt.show()
#


plt.figure(figsize=(8,4))
# Distance euler
plt.subplot(121)
plt.plot(np.linspace(20,(len(dist_euler)+20),len(dist_euler)),dist_euler['error'],label='Distance error')
plt.ylabel(r'Displacement in distance',size=12)
plt.xlabel(r'Number of integration points',size=12)
plt.title('Verlet method: Distance between start and end position for 1 year. Ideally zero.')
#plt.legend()
plt.grid(True)

# Distance verlet
plt.subplot(122)
plt.plot(np.linspace(20,(len(dist)+20),len(dist)),dist['error'],label='Distance error')
plt.ylabel('Displacement in distance')
plt.xlabel(r'Number of integration points',size=12)
plt.title('Verlet method: Distance between start and end position for 1 year. Ideally zero.')
#plt.legend()
plt.grid(True)
plt.savefig('dist.pdf', bbox_inches='tight')
plt.show()

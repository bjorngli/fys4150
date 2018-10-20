import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

energy = pd.read_csv('energy_error.dat',sep='\s+',header=None,names=['error'])
momentum = pd.read_csv('momentum_error.dat',sep='\s+',header=None,names=['error'])
dist = pd.read_csv('dist_error.dat',sep='\s+',header=None,names=['error'])
dist_euler = pd.read_csv('dist_error_euler.dat',sep='\s+',header=None,names=['error'])
error_line = np.repeat(1e-5, len(energy))

fig=plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

# Energy
#plt.plot(range(len(energy)),energy['error'],label='Energy error')
#plt.plot(range(len(energy)),error_line,label='Error line')
# plt.ylabel('Maximum difference from starting energy')
# plt.xlabel('Number of integration points')
# plt.legend()
# plt.grid(True)
# plt.show()

# Momentum
# plt.plot(range(len(momentum)),np.log(momentum['error']),label='Momentum error')
# plt.plot(range(len(energy)),np.log(error_line),label='Error line')
# plt.ylabel('Logarithm of maximum difference from starting angular momentum')
# plt.xlabel('Number of integration points')
# plt.legend()
# plt.grid(True)
# plt.show()

# Distance verlet
plt.plot(range(10,401),dist['error'],label='Distance error')
plt.ylabel('Displacement in distance')
plt.xlabel('Number of integration points')
plt.title('Distance between start and end position for 1 year. Ideally zero.')
plt.legend()
plt.grid(True)
plt.show()

# Distance euler
plt.plot(range(10,401),dist_euler['error'],label='Distance error')
plt.ylabel('Displacement in distance')
plt.xlabel('Number of integration points')
plt.title('Distance between start and end position for 1 year. Ideally zero.')
plt.legend()
plt.grid(True)
plt.show()

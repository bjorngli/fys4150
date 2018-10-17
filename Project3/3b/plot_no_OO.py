import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#df = pd.read_csv('cpp.dat')

x_euler = np.loadtxt('xpos_euler_noOO.dat')
y_euler = np.loadtxt('ypos_euler_noOO.dat')
x_verlet = np.loadtxt('xpos_verlet_noOO.dat')
y_verlet = np.loadtxt('ypos_verlet_noOO.dat')

plt.plot(x_euler,y_euler)
#plt.plot(x_verlet,y_verlet)
plt.show()

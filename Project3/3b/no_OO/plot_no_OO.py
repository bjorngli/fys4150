import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

x_euler = np.loadtxt('xpos_euler_noOO.dat')
y_euler = np.loadtxt('ypos_euler_noOO.dat')
x_verlet = np.loadtxt('xpos_verlet_noOO.dat')
y_verlet = np.loadtxt('ypos_verlet_noOO.dat')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(8,6))
plt.subplot(121)
a=plt.plot(x_euler,y_euler,'b')
plt.title(r'\textbf{Euler}',fontsize=16)
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
plt.setp(a,linewidth=2)
plt.axis('equal')
plt.plot([0],[0],marker='o',markersize=10,color='y')
plt.subplot(122)
plt.title(r'\textbf{Verlet}',fontsize=16)
plt.xlabel(r'x position (AU)',fontsize=14)
#plt.ylabel(r'y position (AU)',fontsize=14)
b=plt.plot(x_verlet,y_verlet,'b')
plt.setp(b,linewidth=2)
plt.axis('equal')
plt.plot([0],[0],marker='o',markersize=10,color='y')
plt.savefig('eulervsverlet.pdf',bbox_inches='tight')
plt.show()

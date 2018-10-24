import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

planets = ['Earth','Jupiter','Mars','Venus','Saturn','Mercury','Uranus','Neptune','Pluto']

earth = pd.read_csv('Earth.dat',sep='\\s+',header=None,names=['x','y','z'])
jupiter = pd.read_csv('Jupiter.dat',sep='\\s+',header=None,names=['x','y','z'])
mars = pd.read_csv('Mars.dat',sep='\\s+',header=None,names=['x','y','z'])
venus = pd.read_csv('Venus.dat',sep='\\s+',header=None,names=['x','y','z'])
saturn = pd.read_csv('Saturn.dat',sep='\\s+',header=None,names=['x','y','z'])
mercury = pd.read_csv('Mercury.dat',sep='\\s+',header=None,names=['x','y','z'])
uranus = pd.read_csv('Uranus.dat',sep='\\s+',header=None,names=['x','y','z'])
neptune = pd.read_csv('Neptune.dat',sep='\\s+',header=None,names=['x','y','z'])
pluto = pd.read_csv('Pluto.dat',sep='\\s+',header=None,names=['x','y','z'])

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
ax.set_zlabel(r'z position (AU)',fontsize=14)
plt.plot(earth['x'],earth['y'],earth['z'])
plt.plot(jupiter['x'],jupiter['y'],jupiter['z'])
plt.plot(mars['x'],mars['y'],mars['z'])
plt.plot(venus['x'],venus['y'],venus['z'])
plt.plot(saturn['x'],saturn['y'],saturn['z'])
plt.plot(mercury['x'],mercury['y'],mercury['z'])
plt.plot(uranus['x'],uranus['y'],uranus['z'])
plt.plot(neptune['x'],neptune['y'],neptune['z'])
plt.plot(pluto['x'],pluto['y'],pluto['z'])
plt.grid(True)
plt.savefig('solar.pdf',bbox_inches='tight')
plt.show()

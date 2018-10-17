import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

earth = pd.read_csv('Earth.dat',sep='\s+',header=None,names=['x','y','z'])
jupiter = pd.read_csv('Jupiter.dat',sep='\s+',header=None,names=['x','y','z'])
mars = pd.read_csv('Mars.dat',sep='\s+',header=None,names=['x','y','z'])
venus = pd.read_csv('Venus.dat',sep='\s+',header=None,names=['x','y','z'])
saturn = pd.read_csv('Saturn.dat',sep='\s+',header=None,names=['x','y','z'])
mercury = pd.read_csv('Mercury.dat',sep='\s+',header=None,names=['x','y','z'])
uranus = pd.read_csv('Uranus.dat',sep='\s+',header=None,names=['x','y','z'])
neptune = pd.read_csv('Neptune.dat',sep='\s+',header=None,names=['x','y','z'])
pluto = pd.read_csv('Pluto.dat',sep='\s+',header=None,names=['x','y','z'])

fig=plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
plt.plot(earth['x'],earth['y'],earth['z'],label='Earth')
plt.plot(jupiter['x'],jupiter['y'],jupiter['z'], label='Jupiter')
plt.plot(mars['x'],mars['y'],mars['z'],label='Mars')
plt.plot(venus['x'],venus['y'],venus['z'], label='Venus')
plt.plot(saturn['x'],saturn['y'],saturn['z'], label='Saturn')
plt.plot(mercury['x'],mercury['y'],mercury['z'], label='Mercury')
plt.plot(uranus['x'],uranus['y'],uranus['z'], label='Uranus')
plt.plot(neptune['x'],neptune['y'],neptune['z'], label='Neptune')
#plt.plot(pluto['x'],pluto['y'],pluto['z'], label='Pluto')
ax.set_xlabel('x [AU]',size=14)
ax.set_ylabel('y [AU]',size=14)
ax.set_zlabel('z [AU]',size=14)
ax.set_zlim(-2,2) 
ax.legend()
ax.grid(True)
plt.show()

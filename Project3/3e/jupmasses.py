import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

sun=pd.read_csv('3c_Sun.dat',sep='\s+',header=None,names=['x','y','z'])
sun10=pd.read_csv('3c_Sun10.dat',sep='\s+',header=None,names=['x','y','z'])
sun1000=pd.read_csv('3c_Sun1000.dat',sep='\s+',header=None,names=['x','y','z'])
earth = pd.read_csv('3c_Earth.dat',sep='\s+',header=None,names=['x','y','z'])
earth10 = pd.read_csv('3c_Earth10.dat',sep='\s+',header=None,names=['x','y','z'])
earth1000 = pd.read_csv('3c_Earth1000.dat',sep='\s+',header=None,names=['x','y','z'])
jupiter = pd.read_csv('3c_Jupiter.dat',sep='\s+',header=None,names=['x','y','z'])
jupiter10 = pd.read_csv('3c_Jupiter10.dat',sep='\s+',header=None,names=['x','y','z'])
jupiter1000 = pd.read_csv('3c_Jupiter1000.dat',sep='\s+',header=None,names=['x','y','z'])

#fig=plt.figure(figsize=(8,8))
'''
ax = fig.add_subplot(111, projection='3d')
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
ax.set_zlabel(r'z position (AU)',fontsize=14)
plt.plot(sun1000['x'],sun1000['y'],sun1000['z'],color='g',label='Sun')
plt.plot(earth1000['x'],earth1000['y'],earth1000['z'],color='b',label='Earth')
plt.plot(jupiter1000['x'],jupiter1000['y'],jupter1000['z'],color='r',label='Jupiter')
plt.legend()
plt.grid(True)
plt.savefig('jup10002d.pdf',bbox_inches='tight')
plt.show()
'''
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
#ax.set_zlabel(r'z position (AU)',fontsize=14)
plt.plot(sun['x'],sun['y'],color='g',label='Sun')
plt.plot(earth['x'],earth['y'],color='b',label='Earth')
plt.plot(jupiter['x'],jupiter['y'],color='r',label='Jupiter')
plt.legend()
plt.grid(True)
plt.savefig('jup.pdf',bbox_inches='tight')
plt.show()

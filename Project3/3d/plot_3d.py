import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

earth1 = pd.read_csv('earth_2.00.dat',sep='\s+',header=None,names=['x','y','z'])
earth2 = pd.read_csv('earth_2.25.dat',sep='\s+',header=None,names=['x','y','z'])
earth3 = pd.read_csv('earth_2.50.dat',sep='\s+',header=None,names=['x','y','z'])
earth4 = pd.read_csv('earth_2.75.dat',sep='\s+',header=None,names=['x','y','z'])
earth5 = pd.read_csv('earth_3.00.dat',sep='\s+',header=None,names=['x','y','z'])




fig=plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
plt.plot(earth5['x'],earth5['y'],label=r'$\beta$ = 3.00')
plt.plot(earth4['x'],earth4['y'],label=r'$\beta$ = 2.75')
plt.plot(earth3['x'],earth3['y'],label=r'$\beta$ = 2.50')
plt.plot(earth2['x'],earth2['y'],label=r'$\beta$ = 2.25')
plt.plot(earth1['x'],earth1['y'],label=r'$\beta$ = 2.00')
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig('beta.pdf',bbox_inches='tight')
plt.show()

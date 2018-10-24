import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Sun_escape = pd.read_csv('Sun_escape.dat',sep='\s+',header=None,names=['x','y','z'])
Planet_escape_88 = pd.read_csv('Planet_escape_88.dat',sep='\s+',header=None,names=['x','y','z'])
Planet_escape_87 = pd.read_csv('Planet_escape_87.dat',sep='\s+',header=None,names=['x','y','z'])

pax_88=Planet_escape_88['x'][Planet_escape_88.index[-1]]
pay_88=Planet_escape_88['y'][Planet_escape_88.index[-1]]

pax_87=Planet_escape_87['x'][Planet_escape_87.index[-1]]
pay_87=Planet_escape_87['y'][Planet_escape_87.index[-1]]


sux=Sun_escape['x'][Sun_escape.index[-1]]
suy=Sun_escape['y'][Sun_escape.index[-1]]

plt.figure(figsize=(8,4))
plt.subplot(121)
plt.plot(Sun_escape['x'],Sun_escape['y'],'yo',label=r'Sun')#,Sun_escape['z']
plt.plot(Planet_escape_88['x'],Planet_escape_88['y'],label=r'$v=8.88 AU/year$')#,Planet_escape['z']
plt.plot([pax_88],[pay_88],'bo', label=r"Earth")
plt.xlabel('x position (AU)',size=12)
plt.ylabel('y position (AU)',size=12)
plt.legend()
plt.grid(True)

plt.subplot(122)
plt.plot(Sun_escape['x'],Sun_escape['y'],'yo',label=r'Sun')#,Sun_escape['z']
plt.plot(Planet_escape_87['x'],Planet_escape_87['y'],label=r'$v=8.87 AU/year$')#,Planet_escape['z']
plt.plot([pax_87],[pay_87],'bo', label=r"Earth")
plt.xlabel('x position (AU)',size=12)
plt.legend()
plt.grid(True)
#plt.savefig('esc.pdf',bbox_inches='tight')
plt.show()

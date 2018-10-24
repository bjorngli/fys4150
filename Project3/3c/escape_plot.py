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

earth = pd.read_csv('3b_Earth_verlet.dat',sep='\s+',header=None,names=['x','y','z'])
orbit = pd.read_csv('3b_Earth_verlet1.dat',sep='\s+',header=None,names=['x','y','z'])
sun = pd.read_csv('3b_Sun.dat',sep='\s+',header=None,names=['x','y','z'])

fig=plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)#, projection='3d')
plt.plot(earth['x'],earth['z'],'r')#,earth['z'],'r')
plt.plot(orbit['x'],orbit['y'],'b--')#,orbit['z'],'b--')
plt.plot(sun['x'],sun['y'],'yo')#,sun['z'],'y')
#p = plt.Circle((0., 0.), 0.1,color='y')
#ax.add_patch(p)
#art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")
plt.title(r'\textbf{Earth escaping its orbit}',fontsize=16)
plt.xlabel(r'x position (AU)',fontsize=14)
plt.ylabel(r'y position (AU)',fontsize=14)
plt.legend([r'$v=\sqrt{8}\pi  AU/yr$',r'$v=2\pi  AU/yr$'],fontsize=14)
#ax.set_zlabel(r'z position (AU)')
plt.axis('equal')
#plt.savefig('earth_escape1.pdf',bbox_inches='tight')
plt.show()

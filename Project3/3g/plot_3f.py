import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

earth = pd.read_csv('3f_pos.dat',sep='\s+',header=None,names=['x','y'])

print(max(np.arctan(earth['y']/earth['x'])))

plt.plot(np.linspace(0,100,len(earth['y'])),np.arctan(earth['y']/earth['x']))
plt.xlabel(r'Earth years',fontsize=14)
plt.ylabel(r'Angle',fontsize=14)
#plt.title(r"\textbf{Mercury's Perihelion Precession}")
plt.grid(True)
plt.savefig('perihelion.pdf',bbox_inches='tight')
plt.show()

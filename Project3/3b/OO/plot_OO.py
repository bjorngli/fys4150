import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

earth = pd.read_csv('3b_Earth_verlet.dat',sep='\s+',header=None,names=['x','y','z'])


fig=plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
plt.plot(earth['x'],earth['y'],earth['z'])
plt.grid(True)
plt.show()

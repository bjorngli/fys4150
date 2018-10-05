import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('cpp.dat')

x = np.loadtxt('xposition_a.dat')
y = np.loadtxt('yposition_a.dat')

plt.plot(x,y)
plt.show()

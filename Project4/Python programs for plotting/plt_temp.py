import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df=pd.read_csv('counter.dat',header=None,sep='\s+',names=(['Accepted']),index_col = 0)

E1=pd.to_numeric(df['Accepted'])
x = np.linspace(1,5.1,len(E1))
E1_new = E1 / E1.index

plt.title(r'\textbf{Fraction of accepted cycles as a function of T}',fontsize=14)
plt.plot(x,E1_new)
plt.grid(True)
plt.axvline(x=2.269, color='k', linestyle='--')
plt.ylabel(r'Fraction of accepted MC cycles',fontsize=14)
plt.xlabel(r'Temperature',fontsize=14)
plt.savefig('temp.pdf',bbox_inches='tight')
plt.show()

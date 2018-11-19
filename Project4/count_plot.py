import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df=pd.read_csv('count1.dat',header=None,sep='\s+',names=(['mcs','frac']))
df1=pd.read_csv('count2.dat',header=None,sep='\s+',names=(['mcs','frac']))

x=df['mcs']


plt.plot(x,df['frac'],x,df1['frac'])
plt.title(r'\textbf{Number of accepted states for T=1.0 and T=2.4}',fontsize=15)
plt.xlabel(r'Number of states',fontsize=14)
plt.ylabel(r'Number of accepted states',fontsize=14)
plt.legend([r'T=1.0',r'T=2.4'],fontsize=13)
plt.savefig('count.pdf',bbox_inches='tight')
plt.show()

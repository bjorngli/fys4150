import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df=pd.read_csv('mc_calc_random.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))
df1=pd.read_csv('mc_calc_uniform.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))

#file1=df[:(len(df)/2)]
#file2=df[(len(df)/2):]
E1=pd.to_numeric(df['E'][1.0])
E2=pd.to_numeric(df['E'][2.4])
absM1=pd.to_numeric(df['absM'][1.0])
absM2=pd.to_numeric(df['absM'][2.4])
M1=pd.to_numeric(df['M'][1.0])
M2=pd.to_numeric(df['M'][2.4])

E1u=pd.to_numeric(df1['E'][1.0])
E2u=pd.to_numeric(df1['E'][2.4])
absM1u=pd.to_numeric(df1['absM'][1.0])
absM2u=pd.to_numeric(df1['absM'][2.4])
M1u=pd.to_numeric(df1['M'][1.0])
M2u=pd.to_numeric(df1['M'][2.4])
x=np.linspace(1,10000,len(df)/2)

fig=plt.figure(figsize=(8,6))

plt.subplot(221)
plt.title(r'\textbf{Random lattice}',fontsize=14)
plt.plot(x,E1,x,E2)
plt.grid(True)
plt.ylim(-2.1,-0.7)
plt.legend([r'T=1.0',r'T=2.4'])
plt.ylabel(r'$\langle E \rangle$',fontsize=14)

plt.subplot(222)
plt.title(r'\textbf{All spins up lattice}',fontsize=14)
plt.plot(x,E1u,x,E2u)
plt.grid(True)
plt.ylim(-2.1,-0.7)
plt.legend([r'T=1.0',r'T=2.4'])

plt.subplot(223)
plt.plot(x,absM1,x,absM2)
plt.grid(True)
plt.ylim(0.0,1.1)
plt.legend([r'T=1.0',r'T=2.4'])
plt.xlabel(r'Number of Monte Carlo cycles',fontsize=13)
plt.ylabel(r'$\langle |{\mathcal{M}}| \rangle$',fontsize=14)

plt.subplot(224)
plt.plot(x,absM1u,x,absM2u)
plt.grid(True)
plt.ylim(0.0,1.1)
plt.legend([r'T=1.0',r'T=2.4'])
plt.xlabel(r'Number of Monte Carlo cycles',fontsize=13)
plt.savefig('eandm.pdf',bbox_inches='tight')
plt.show()

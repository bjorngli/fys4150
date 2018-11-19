import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df20=pd.read_csv('mc_calc20.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))
df40=pd.read_csv('mc_calc40.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))
df60=pd.read_csv('mc_calc60.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))
df80=pd.read_csv('mc_calc80.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))
df100=pd.read_csv('mc_calc100.dat',header=None,sep='\s+',names=(['E','Cv','M','chi','absM']))

x=np.linspace(2.0,2.4,len(df20))
E20=pd.to_numeric(df20['E'])
Cv20=pd.to_numeric(df20['Cv'])
chi20=pd.to_numeric(df20['chi'])
absM20=pd.to_numeric(df20['absM'])

E40=pd.to_numeric(df40['E'])
Cv40=pd.to_numeric(df40['Cv'])
chi40=pd.to_numeric(df40['chi'])
absM40=pd.to_numeric(df40['absM'])

E60=pd.to_numeric(df60['E'])
Cv60=pd.to_numeric(df60['Cv'])
chi60=pd.to_numeric(df60['chi'])
absM60=pd.to_numeric(df60['absM'])

E80=pd.to_numeric(df80['E'])
Cv80=pd.to_numeric(df80['Cv'])
chi80=pd.to_numeric(df80['chi'])
absM80=pd.to_numeric(df80['absM'])

E100=pd.to_numeric(df100['E'])
Cv100=pd.to_numeric(df100['Cv'])
chi100=pd.to_numeric(df100['chi'])
absM100=pd.to_numeric(df100['absM'])
print('t_cv:',Cv20.idxmax(),Cv40.idxmax(),Cv60.idxmax(),Cv80.idxmax(),Cv100.idxmax(),'t_chi:',chi20.idxmax(),chi40.idxmax(),chi60.idxmax(),chi80.idxmax(),chi100.idxmax())

fig=plt.figure(figsize=(8,8))

plt.subplot(411)
plt.title(r'\textbf{Expectation values for lattices of size L $\times$ L}',fontsize=15)
plt.plot(x,E20,x,E40,x,E60,x,E80,x,E100)
plt.axvline(x=2.269, color='k', linestyle='--')
plt.grid(True)
plt.xlim(2.0,2.4)
plt.legend([r'L=20',r'L=40',r'L=60',r'L=80',r'L=100'])
plt.ylabel(r'$\langle E \rangle$',fontsize=14)

plt.subplot(412)
plt.plot(x,absM20,x,absM40,x,absM60,x,absM80,x,absM100)
plt.axvline(x=2.269, color='k', linestyle='--')
plt.grid(True)
plt.xlim(2.0,2.4)
#plt.legend(['L=20','L=40','L=60','L=80','L=100'])
plt.ylabel(r'$\langle |{\mathcal{M}}| \rangle$',fontsize=14)

plt.subplot(413)
plt.plot(x,Cv20,x,Cv40,x,Cv60,x,Cv80,x,Cv100)
plt.axvline(x=2.269, color='k', linestyle='--')
plt.grid(True)
plt.xlim(2.0,2.4)
#plt.legend(['L=20','L=40','L=60','L=80','L=100'])
plt.ylabel(r'$\langle C_V \rangle$',fontsize=14)

plt.subplot(414)
plt.plot(x,chi20,x,chi40,x,chi60,x,chi80,x,chi100)
plt.axvline(x=2.269, color='k', linestyle='--')
plt.grid(True)
plt.xlim(2.0,2.4)
#plt.legend(['L=20','L=40','L=60','L=80','L=100'])
plt.xlabel(r'Temperature',fontsize=14)
plt.ylabel(r'$\langle \chi \rangle$',fontsize=14)
plt.savefig('last.pdf',bbox_inches='tight')
plt.show()

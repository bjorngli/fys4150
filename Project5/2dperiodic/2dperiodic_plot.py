import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n = 41

df = pd.read_csv('2dperiodic_0.dat',sep='\s+',header=None)
df1 = pd.read_csv('2dperiodic_150.dat',sep='\s+',header=None)

df2 = df.iloc[:,0:n].transpose()
df3 = df.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)

for i in range(1,41):
    df_temp = df.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    df2 = pd.concat([df2,df_temp],ignore_index=True,axis=1,join_axes=[df2.index])



dff2 = df1.iloc[:,0:n].transpose()
dff3 = df1.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)

for i in range(1,41):
    dff_temp = df1.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    dff2 = pd.concat([dff2,dff_temp],ignore_index=True,axis=1,join_axes=[dff2.index])


X=df2.columns.values
Y=df2.index.values
Z=df2.values
x,y=np.meshgrid(X, Y)

X2=dff2.columns.values
Y2=dff2.index.values
Z2=dff2.values
x2,y2=np.meshgrid(X2, Y2)


plt.figure(figsize=(10,5))
plt.suptitle(r'Contour field for $\psi (x,y,t)$ at t = 0 and t = $150$',fontsize=15)
plt.subplot(121)
plt.style.use("ggplot")
CS = plt.contourf(x, y, Z.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (west-east)', fontsize = 13)
plt.ylabel('y (south-north)', fontsize = 13)

plt.subplot(122)
plt.style.use("ggplot")
CS = plt.contourf(x2, y2, Z2.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (west-east)', fontsize = 13)
plt.ylabel('y (south-north)', fontsize = 13)
plt.savefig('2dperiodic.pdf',bbox_inches='tight')

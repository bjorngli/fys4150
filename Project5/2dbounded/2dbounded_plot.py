import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n = 41

df1 = pd.read_csv('2dbounded_0.dat',sep='\s+',header=None)
df2 = pd.read_csv('2dbounded_50.dat',sep='\s+',header=None)
df3 = pd.read_csv('2dbounded_100.dat',sep='\s+',header=None)
df4 = pd.read_csv('2dbounded_150.dat',sep='\s+',header=None)

df_2 = df1.iloc[:,0:n].transpose()
df_3 = df1.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)
dftest = pd.concat([df_2,df_3],ignore_index=True,axis=1,join_axes=[df_2.index])


for i in range(1,41):
    df_temp = df1.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    df_2 = pd.concat([df_2,df_temp],ignore_index=True,axis=1,join_axes=[df_2.index])

dff_2 = df2.iloc[:,0:n].transpose()
dff_3 = df2.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)
dfftest = pd.concat([dff_2,dff_3],ignore_index=True,axis=1,join_axes=[dff_2.index])


for i in range(1,41):
    dff_temp = df2.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    dff_2 = pd.concat([dff_2,dff_temp],ignore_index=True,axis=1,join_axes=[dff_2.index])

dfff_2 = df3.iloc[:,0:n].transpose()
dfff_3 = df3.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)
dffftest = pd.concat([dfff_2,dfff_3],ignore_index=True,axis=1,join_axes=[dfff_2.index])


for i in range(1,41):
    dfff_temp = df3.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    dfff_2 = pd.concat([dfff_2,dfff_temp],ignore_index=True,axis=1,join_axes=[dfff_2.index])

dffff_2 = df4.iloc[:,0:n].transpose()
dffff_3 = df4.iloc[:,n*1:n*1+n].transpose().reset_index().drop(['index'],axis=1)
dfffftest = pd.concat([dffff_2,dffff_3],ignore_index=True,axis=1,join_axes=[dffff_2.index])


for i in range(1,41):
    dffff_temp = df4.iloc[:,n*i:n*i+n].transpose().reset_index().drop(['index'],axis=1)
    dffff_2 = pd.concat([dffff_2,dffff_temp],ignore_index=True,axis=1,join_axes=[dffff_2.index])


X = np.linspace(0,1,41)
Y = np.linspace(0,1,41)



plt.figure(figsize=(10,10))
plt.suptitle(r'Contour plot of the two-dimensional sinus wave',fontsize=15)
plt.subplot(221)
plt.style.use("ggplot")
Z=df_2.values
x,y=np.meshgrid(X, Y)
CS = plt.contourf(x, y, Z.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.clim(-2,2)
plt.title(r'$\psi(x,y,0) = sin(\pi y)sin(4\pi x)$')
plt.ylabel('y (south-north)')
plt.xlabel('x (west-east)')

plt.subplot(222)
plt.style.use("ggplot")
Z=dff_2.values
x,y=np.meshgrid(X, Y)
CS = plt.contourf(x, y, Z.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.clim(-2,2)
plt.title(r'$\psi(x,y,50) = sin(\pi y)sin(4\pi x)$')
plt.ylabel('y (south-north)')
plt.xlabel('x (west-east)')

plt.subplot(223)
plt.style.use("ggplot")
Z=dfff_2.values
x,y=np.meshgrid(X, Y)
CS = plt.contourf(x, y, Z.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.clim(-2,2)
plt.title(r'$\psi(x,y,100) = sin(\pi y)sin(4\pi x)$')
plt.ylabel('y (south-north)')
plt.xlabel('x (west-east)')

plt.subplot(224)
plt.style.use("ggplot")
Z=dffff_2.values
x,y=np.meshgrid(X, Y)
CS = plt.contourf(x, y, Z.transpose(), 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.clim(-2,2)
plt.title(r'$\psi(x,y,150) = sin(\pi y)sin(4\pi x)$')
plt.ylabel('y (south-north)')
plt.xlabel('x (west-east)')

plt.savefig('2dbounded_all.pdf',bbox_inches='tight')

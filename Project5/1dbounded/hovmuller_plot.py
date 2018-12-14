import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('HovmullerBoundedSine.dat',header=None,sep='\s+')
df2 = pd.read_csv('HovmullerBoundedGaus1.dat',header=None,sep='\s+')
df3 = pd.read_csv('HovmullerBoundedGaus2.dat',header=None,sep='\s+')

# Hovmuller Bounded sine
Z=df.values

x = np.linspace(0,1,41)
t = np.linspace(0,150,7501)

plt.figure(4)
plt.style.use("ggplot")
fig = plt.figure(figsize = (9,10))
CS = plt.contourf(x, t, Z, 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (west-east)', fontsize = 13)
plt.ylabel('Time, t', fontsize = 13)
plt.title(r'Hovmüller diagram of $\psi(x, t)$ from I.C. wave $\psi(x, 0) = sin(4 \pi x)$')
plt.savefig('HovmullerBoundedSine.pdf',bbox_inches='tight')

# Hovmuller Bounded Gaus sigma 0.1
Z=df2.values

x = np.linspace(0,1,41)
t = np.linspace(0,150,7501)

plt.figure(4)
plt.style.use("ggplot")
fig = plt.figure(figsize = (9,10))
CS = plt.contourf(x, t, Z, 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (east-west)', fontsize = 13)
plt.ylabel('Time, t', fontsize = 13)
plt.title(r'Hovmüller diagram of $\psi(x, t)$ from I.C. wave $\psi(x, 0) = e^{-(\frac{x-0.5}{\sigma})^2}$ for $\sigma$ = 0.1')
plt.savefig('HovmullerBoundedGaus1.pdf')


# Hovmuller Bounded Gaus sigma 0.2
Z=df3.values

x = np.linspace(0,1,41)
t = np.linspace(0,150,7501)

plt.figure(4)
plt.style.use("ggplot")
fig = plt.figure(figsize = (9,10))
CS = plt.contourf(x, t, Z, 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (east-west)', fontsize = 13)
plt.ylabel('Time, t', fontsize = 13)
plt.title(r'Hovmüller diagram of $\psi(x, t)$ from I.C. wave $\psi(x, 0) = e^{-(\frac{x-0.5}{\sigma})^2}$ for $\sigma$ = 0.2')
plt.savefig('HovmullerBoundedGaus2.pdf')


# Hovmuller Bounded Gaus sigma 0.1 and sigma 0.2
Z=df2.values
Z1 = df3.values

x = np.linspace(0,1,41)
t = np.linspace(0,150,7501)

plt.figure(figsize=(10,5))
plt.suptitle(r'Hovmüller diagram of $\psi(x, t)$ from I.C. wave $\psi(x, 0) = e^{-(\frac{x-0.5}{\sigma})^2}$ for $\sigma = 0.1$ and $0.2$')
plt.subplot(121)
plt.style.use("ggplot")
CS = plt.contourf(x, t, Z, 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (west-east)', fontsize = 13)
plt.ylabel('Time, t', fontsize = 13)

plt.subplot(122)
plt.style.use("ggplot")
CS = plt.contourf(x, t, Z1, 20, cmap = plt.cm.BrBG)
plt.colorbar(CS, orientation = "horizontal")
plt.xlabel('x (west-east)', fontsize = 13)
plt.ylabel('Time, t', fontsize = 13)
plt.savefig('HovmullerBoundedGausBoth.pdf',bbox_inches='tight')

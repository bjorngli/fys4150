import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np

h_list = np.zeros(7)
errors = np.zeros(7)
def u(x):
    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return u

for i in range(1,8):
    N=10**i
    h = 1./(N+1)
    h_list[i-1]=h
    x = np.linspace(h,1-h,N-1)
    results = np.loadtxt('s_vector%i.dat'%(i))
    errors[i-1] = max(abs((results[1:N]-u(x))/u(x)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.loglog(h_list,errors,'o-')
plt.title('Relative error for different values of h',fontsize=18)
plt.legend()
plt.xlabel(r'\textit{$log_{10}(h)$}',fontsize=16)
plt.ylabel(r'\textit{$max(log_{10}(|v_i-u_i/u_i|)$',fontsize=16)
plt.show()

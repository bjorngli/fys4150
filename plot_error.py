import matplotlib.pyplot as plt
import numpy as np
import os
import sys

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
print(errors)
plt.loglog(h_list,errors,'o-')
plt.title('Relative error for different values of h')
plt.legend()
plt.xlabel('log10(h)')
plt.ylabel('log10(v-u/u)')
plt.show()

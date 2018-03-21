import time
from sympy import *
c=cos
s=sin

import numpy as np
import matplotlib.pyplot as plt

def W_i(t_i,t_i_min,t_i_max,tau_i,W_o):
    if t_i<=t_i_min:
        return W_o
    elif t_i<=t_i_min+tau_i:
        return W_o/2.0*(1+cos(pi*(t_i-t_i_min)/tau_i))
    elif t_i<t_i_max-tau_i:
        return 0
    elif t_i<t_i_max:
        return W_o/2.0*(1+cos(pi*(t_i_max-t_i)/tau_i))
    else:
        return W_o
        
#theta 1---------------------------------------------------------------------
plt.ylabel('Weight')
plt.xlabel('Joint Angle')
plt.ylim(-0.01,3.1)
plt.xlim(-3.14,3.14)

t_list=[]
W_list=[]
t=-3.14

t_min=-2
t_max=2
tau=(t_max-t_min)/4.0

delta=6.28/1000
for i in range(1000):
    t=t+delta
    t_list.append(t)
    W_list.append(W_i(t,t_min,t_max,tau,3))
    print(t,W_i(t,t_min,t_max,tau,3))

plt.plot(t_list,W_list)
plt.show()

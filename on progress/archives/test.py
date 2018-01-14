import numpy as np
import matplotlib.pyplot as plt

t1_min=5
t1_max=10
t1_min_list=[]
t1_max_list=[]
for i in range(300):
    t1_min_list.append(t1_min)
for i in range(300):
    t1_max_list.append(t1_max)
#Plotting
x=np.arange(0, 300, 1)
plt.figure(1)

plt.subplot(111)
plt.ylim(0,15)
plt.ylabel('joint angle #1')
plt.plot(x,t1_min_list,'r--')
plt.plot(x,t1_max_list,'r--')

plt.show()

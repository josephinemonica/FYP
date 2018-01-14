x=np.arange(0, 300, 1)
plt.figure(1)

#theta 1---------------------------------------------------------------------
plt.subplot(411)
plt.ylabel('joint angle #1')
plt.ylim(float(t1_min)-0.1,float(t1_max)+0.1)
plt.plot(t1_list)

#theta 1 limit
t1_min_list=[]
t1_max_list=[]
t1_relax_low_list=[]
t1_relax_up_list=[]
for i in range(300):
    t1_min_list.append(t1_min)
for i in range(300):
    t1_max_list.append(t1_max)
for i in range(300):
    t1_relax_low_list.append(t1_min+tau1)
for i in range(300):
    t1_relax_up_list.append(t1_max-tau1)
plt.plot(x,t1_min_list,'r--')
plt.plot(x,t1_max_list,'r--')
plt.plot(x,t1_relax_low_list,'y--')
plt.plot(x,t1_relax_up_list,'y--')

#theta 2---------------------------------------------------------------------
plt.subplot(412)
plt.ylabel('joint angle #2')
plt.ylim(float(t2_min)-0.1,float(t2_max)+0.1)
plt.plot(t2_list)

#theta 2 limit
t2_min_list=[]
t2_max_list=[]
t2_relax_low_list=[]
t2_relax_up_list=[]
for i in range(300):
    t2_min_list.append(t2_min)
for i in range(300):
    t2_max_list.append(t2_max)
for i in range(300):
    t2_relax_low_list.append(t2_min+tau2)
for i in range(300):
    t2_relax_up_list.append(t2_max-tau2)
plt.plot(x,t2_min_list,'r--')
plt.plot(x,t2_max_list,'r--')
plt.plot(x,t2_relax_low_list,'y--')
plt.plot(x,t2_relax_up_list,'y--')

#theta 3---------------------------------------------------------------------
plt.subplot(413)
plt.ylabel('joint angle #3')
plt.ylim(float(t3_min)-0.1,float(t3_max)+0.1)
plt.plot(t3_list)

#theta 3 limit
t3_min_list=[]
t3_max_list=[]
t3_relax_low_list=[]
t3_relax_up_list=[]
for i in range(300):
    t3_min_list.append(t3_min)
for i in range(300):
    t3_max_list.append(t3_max)
for i in range(300):
    t3_relax_low_list.append(t3_min+tau3)
for i in range(300):
    t3_relax_up_list.append(t3_max-tau3)
plt.plot(x,t3_min_list,'r--')
plt.plot(x,t3_max_list,'r--')
plt.plot(x,t3_relax_low_list,'y--')
plt.plot(x,t3_relax_up_list,'y--')

#theta 4---------------------------------------------------------------------
plt.subplot(414)
plt.ylabel('joint angle #4')
plt.ylim(-3.15,3.15)
plt.plot(t4_list)


#Compute Jacobian matrix of the UAV+arm system
import numpy as np

import math
c=math.cos
s=math.sin
####################################################################################################################################
#System parameters
xb=0.11316
l1=0.0695
l2=0.17
l3=0.07025
l4=0.025
alpha1=np.pi/2
alpha2=0
alpha3=np.pi/2
alpha4=0
t1=1
t2=2
t3=3
t4=4
a1=0
a2=l2
a3=0
a4=0
b1=l1
b2=0
b3=0
b4=l3+l4
DH_parameters=[[a1,b1,alpha1,t1],[a2,b2,alpha2,t2],[a3,b3,alpha3,t3],[a4,b4,alpha4,t4]]
n=4 #number of DOF for arms

#BASE
theta=0
phi=0
psi=0
x=0
y=0
z=0
p_b=np.matrix([[x],[y],[z]])

a_vec_list=[]
###################################################################################################################

#Computing rotation matrix Q
def Q(alpha,theta):
    mu=s(alpha)
    lam=c(alpha)
    return np.matrix([[c(theta),-lam*s(theta),mu*s(theta)],
                    [s(theta),lam*c(theta),-mu*c(theta)],
                    [0,mu,lam]])

#Computing a_vector, vector from O_i to O_i+1
def a_vector(a,b,theta):
    return np.matrix([[a*c(theta)],[a*s(theta)],[b]])
#########################################################################################################    
#compute P and Q
Q_list=[]                   #list of rotation matrix Q, will be reused in the second part of jacobian computation
P_list=[]                   #list of P, will be reused in the second part of jacobian computation
    
#initialize P_0 with identity
P_dummy=np.identity(3)
for i in range(n):              #compute all Q and P
    a=DH_parameters[i][0]         
    b=DH_parameters[i][1]             
    alpha=DH_parameters[i][2]        
    theta=DH_parameters[i][3]
    Q_dummy=Q(alpha,theta)
    Q_list.append(Q_dummy)
    P_dummy=np.dot(P_dummy,Q_dummy)
    P_list.append(P_dummy)
###########################################################################################################
#compute Jacobian Matrix        
#initialize matrix: 6xn
J=np.zeros(shape=(6,n))
####evaluation of upper jacobian
#fill first column of upper half jacobian
J[0][0]=0
J[1][0]=0
J[2][0]=1
     
for i in range(1,n):                    # i= 1 to n-1
    #ith column of J = last column(z component) of P_i
    P_dummy=P_list[i-1]
    J[0][i]=P_dummy.item((0,2))
    J[1][i]=P_dummy.item((1,2))
    J[2][i]=P_dummy.item((2,2))
    
####evaluation of lower jacobian
    
#initialize
a=DH_parameters[n-1][0]      
b=DH_parameters[n-1][1]     
alpha=DH_parameters[n-1][2]
theta=DH_parameters[n-1][3]
a_vec=a_vector(a,b,theta)   #a_n #3x1
a_vec_list.insert(0,a_vec)
r_vec=a_vec                 #r_n #3x1
lower=np.dot(P_list[n-2], (np.cross([0,0,1],r_vec.T).T))    #P_n-1 ([0,0,1]xr_vec_n)
    
#fill last column of lower half jacobian
J[3][n-1]=lower[0]
J[4][n-1]=lower[1]
J[5][n-1]=lower[2]
    
for i in range(n-1,0,-1):                   #i= n-1 to 1
    a=DH_parameters[i-1][0]                 #a_n-1 to a_1
    b=DH_parameters[i-1][1]                 #b_n-1 to b_1
    alpha=DH_parameters[i-1][2]             #alpha_n-1 to alpha_1
    theta=DH_parameters[i-1][3]             #theta_n-1 to theta_1
    a_vec=a_vector(a,b,theta)               #
    a_vec_list.insert(0,a_vec)
    Q_dummy=Q_list[i-1]                     #Q_list[i-1]=Q_i..  Q_n-1 to Q_1
       
    P_dummy=P_list[i-2]                     #P_list[i-2]=P_i-1..P_n-2 to P_0, but P_0 =identity
    if(i-2==-1):
        P_dummy=np.identity(3)
    r_vec=a_vec+np.dot(Q_dummy,r_vec) #r_n-1 to r_1
        
    lower=np.dot(P_dummy, (np.cross([0,0,1],r_vec.T)).T)
    J[3][i-1]=lower[0]
    J[4][i-1]=lower[1]
    J[5][i-1]=lower[2]
#J is found
######################################################################################################################
#Rotation matrix from world to base<=> coor transform from base to world
Rb=np.matrix([[c(theta)*c(psi),   s(phi)*s(theta)*c(psi)-c(theta)*s(psi),  s(psi)*s(phi)+c(phi)*s(theta)*c(psi)],
            [c(theta)*s(psi),   c(psi)*c(phi)+s(phi)*s(theta)*s(psi),   c(phi)*s(theta)*s(psi)-s(phi)*c(psi)],
            [-s(theta),         c(theta)*s(phi),                        c(theta)*c(phi)                     ]])
######################################################################################################################
#orientation of EE to world
#Rb Q1 Q2 Q3 ... Qn= Rb Pn
Pn=P_list[len(P_list)-1]
Q_ee=np.dot(Rb,Pn)

#position of EE to world
p_ee_b=np.matrix([[0],[0],[0]])   #position EE w.r.t base
for i in range(n):
    P_dummy=P_list[i-1]                     #P_list[i-2]=P_i-1..P_n-2 to P_0, but P_0 =identity

    a_vec_dummy=a_vec_list[i]

    if(i==0):
        P_dummy=np.identity(3)
    p_ee_b=p_ee_b+np.dot(P_dummy,a_vec_dummy)
    
p_ee=p_b+np.dot(Rb,(np.matrix([xb,0,0])).T)+np.dot(Rb,(np.matrix(p_ee_b)))
#p_ee computed
###################################################################################################################
print(p_ee)
print(Q_ee)

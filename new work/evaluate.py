from sympy import *
c=cos
s=sin

import numpy as np
import matplotlib.pyplot as plt
x,y,z,theta,phi,psi,t1,t2,t3,t4= symbols('x y z theta phi psi t1 t2 t3 t4') #UAV and arm postures
phi=0;theta=0;#underactuated
l1=0.0695; l2=0.17; l3=0.07025; l4=0.025        
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0    #DH params
xshift=0.11316; yshift=0; zshift=0  #shift from base to first link origin
shift=Matrix([[xshift],[yshift],[zshift]])
a1=0; a2=l2; a3=0; a4=0             #DH params
b1=l1; b2=0; b3=0; b4=l3+l4         #DH params

#================================================================================================================
#Compute rotation matrix for euler angles psi,theta, and phi
def R(psi,theta,phi):
    return Matrix([
        [cos(psi)*cos(theta), sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi), sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)], 
        [sin(psi)*cos(theta), sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi), -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)], 
        [-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)]])

#Computing rotation matrix Q of a link frame
#params: alpha, theta -> DH params
def Q(alpha,theta): 
    mu=s(alpha)
    lam=c(alpha)
    return Matrix([[c(theta),-lam*s(theta),mu*s(theta)],
                    [s(theta),lam*c(theta),-mu*c(theta)],
                    [0,mu,lam]])

#Computing a_vector, vector from O_i to O_i+1
def a_vector(a,b,theta):
    return Matrix([[a*c(theta)],[a*s(theta)],[b]])
    
#Computing cross product matrix of 3x1 vector
def CPM(v):
    return Matrix([[0,      -v[2],  v[1]],
                   [v[2],   0,     -v[0]],
                   [-v[1],  v[0],   0]

                   ])
                   
#Computing vector of a matrix of 3x3 matrix; opposite of CPM
def vect(A):
    return 0.5*Matrix([ [A[7]-A[5]],
                    [A[2]-A[6]],
                    [A[3]-A[1]]])
                    
#compute ZYX euler angles given a 3x3rotation matrix   
#https://www.geometrictools.com/Documentation/EulerAngles.pdf          
def euler_angles(R):
    if(R[6]<1):
        if(R[6]>-1):
            theta=asin(-R[6])
            psi=atan2(R[3],R[0])
            phi=atan2(R[7],R[8])
        else:
            theta=pi/2
            psi=-atan2(-R[5],R[4])
            phi=0
    else:
        theta=-pi/2
        psi=atan2(-R[5],R[4])#atan2(R[1][0],R[1][1])
        phi=0

    return Matrix([[psi],[theta],[phi]])
    
#T matrix to transfrom rates of zyx euler angles to angular velocity. omega=T d/dt(psi theta phi)
def T(psi,theta,phi):
    return Matrix([
                    [0, -s(psi),    c(psi)*c(theta) ],
                    [0, c(psi),     s(psi)*c(theta) ],
                    [1, 0,          -s(theta)       ]])  

#=======================================================================================================

#Q1, Q2, Q3, Q4
Q1=Q(alpha1,t1)
Q2=Q(alpha2,t2)
Q3=Q(alpha3,t3)
Q4=Q(alpha4,t4)
#P1, P2, P3, P4
P1=Q1
P2=P1*Q2
P3=P2*Q3
P4=P3*Q4
#vectors: a1_vec, a2_vec, a3_vec, a4_vec
a1_vec=a_vector(a1,b1,t1)
a2_vec=a_vector(a2,b2,t2)
a3_vec=a_vector(a3,b3,t3)
a4_vec=a_vector(a4,b4,t4)
#Rb: Rotation matrix from world to base. Euler angles ZYX : psi theta phi
#Rb=Q_psi*Q_theta*Q_phi
Rb=Matrix([
[cos(psi)*cos(theta), sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi), sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)], 
[sin(psi)*cos(theta), sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi), -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)], 
[-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)]])

Rs=Matrix([
    [1,0,0],
    [0,-1,0],
    [0,0,-1]])

#p_eb= position of ee w.r.t base in base's frame
p_eb=a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec+shift



#*********************************************************
p_b=Matrix([[x],[y],[z]])    #base position
p_1=p_b+Rb*Rs*(shift)
p_2=p_b+Rb*Rs*(shift+a1_vec)
p_3=p_b+Rb*Rs*(shift+a1_vec+P1*a2_vec)

#cylinder axis unit vectors in world frame
e_1=Rb*Rs*Matrix([[0],[0],[1]])
e_2=Rb*Rs*P1*Matrix([[0],[0],[1]])
e_3=Rb*Rs*P2*Matrix([[0],[0],[1]])


#obstacle info
x_obs,y_obs,z_obs,R_obs=symbols('x_obs y_obs z_obs R_obs')
p_obs=Matrix([[x_obs],[y_obs],[z_obs]])
        
p_obs_1=p_obs-p_1
p_obs_2=p_obs-p_2
p_obs_3=p_obs-p_3

#h_1=p_1+((p_obs_1.transpose()*e_1)[0])*e_1
#h_1=p_1
h_1=p_1+l1*e_1
J_1=diff(h_1,x)
J_1=J_1.row_join(diff(h_1,y))
J_1=J_1.row_join(diff(h_1,z))
J_1=J_1.row_join(diff(h_1,psi))
J_1=J_1.row_join(diff(h_1,t1))
J_1=J_1.row_join(diff(h_1,t2))
J_1=J_1.row_join(diff(h_1,t3))
J_1=J_1.row_join(diff(h_1,t4))
u_1=(h_1-p_obs)/sqrt(((h_1-p_obs).transpose()*(h_1-p_obs))[0])
J_1=-u_1.transpose()*J_1
#print(J_1)


#h_2=p_2+((p_obs_2.transpose()*e_2)[0])*e_2
#h_2=p_2
h_2=p_2+l2*e_2
J_2=diff(h_2,x)
J_2=J_2.row_join(diff(h_2,y))
J_2=J_2.row_join(diff(h_2,z))
J_2=J_2.row_join(diff(h_2,psi))
J_2=J_2.row_join(diff(h_2,t1))
J_2=J_2.row_join(diff(h_2,t2))
J_2=J_2.row_join(diff(h_2,t3))
J_2=J_2.row_join(diff(h_2,t4))
u_2=(h_2-p_obs)/sqrt(((h_2-p_obs).transpose()*(h_2-p_obs))[0])
J_2=-u_2.transpose()*J_2
#print(J_2)

#h_3=p_3+((p_obs_3.transpose()*e_3)[0])*e_3
#h_3=p_3
h_3=p_3+l3*e_3
J_3=diff(h_3,x)
J_3=J_3.row_join(diff(h_3,y))
J_3=J_3.row_join(diff(h_3,z))
J_3=J_3.row_join(diff(h_3,psi))
J_3=J_3.row_join(diff(h_3,t1))
J_3=J_3.row_join(diff(h_3,t2))
J_3=J_3.row_join(diff(h_3,t3))
J_3=J_3.row_join(diff(h_3,t4))
u_3=(h_3-p_obs)/sqrt(((h_3-p_obs).transpose()*(h_3-p_obs))[0])
J_3=-u_3.transpose()*J_3
print(J_3)

J_oa=J_1.col_join(J_2)
J_oa=J_oa.col_join(J_3)

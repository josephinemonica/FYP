import time
from sympy import *

c=cos
s=sin

import numpy as np
import matplotlib.pyplot as plt
#================================================================================================================
Rs=Matrix([
    [1,0,0],
    [0,-1,0],
    [0,0,-1]])
    
#Variables
x,y,z,theta,phi,psi,t1,t2,t3,t4= symbols('x y z theta phi psi t1 t2 t3 t4') #UAV and arm postures
l1=0.0695; l2=0.17; l3=0.07025; l4=0.025        
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0    #DH params
xshift=0.11316; yshift=0; zshift=0  #shift from base to first link origin
shift=Matrix([[xshift],[yshift],[zshift]])
a1=0; a2=l2; a3=0; a4=0             #DH params
b1=l1; b2=0; b3=0; b4=l3+l4         #DH params
#=========for plotting===========================================================================================
t1_list=[]
t2_list=[]
t3_list=[]
t4_list=[]
psi_list=[]
theta_list=[]
phi_list=[]
psi_dot_list_d=[] #desired
theta_dot_list_d=[] #desired
phi_dot_list_d=[] #desired
psi_dot_list_actual=[] #actual from the algo
theta_dot_list_actual=[]
phi_dot_list_actual=[]
x_b_list=[]
y_b_list=[]
z_b_list=[]
psi_b_list=[]
x_list=[]
y_list=[]
z_list=[]
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
#==========JOINT LIMITS
t1_min=-pi/6
t1_max=pi/6
t2_min=-pi+pi
t2_max=0+pi
t3_min=0-pi/2
t3_max=pi-pi/2
tau1=(t1_max-t1_min)/4.0
tau2=(t2_max-t2_min)/4.0
tau3=(t3_max-t3_min)/4.0

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
W1,W2,W3=symbols('W1 W2 W3')
W_JL=Matrix([[0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [0,0,0,0,W1,0,0,0],
             [0,0,0,0,0,W2,0,0],
             [0,0,0,0,0,0,W3,0],
             [0,0,0,0,0,0,0,0]])
# W_JL=W_JL.evalf(subs={W1:W_i(t1,t1_min,t1_max,tau1,5),W2:W_i(t2,t2_min,t2_max,tau2,5),W3:W_i(t3,t3_min,t3_max,tau3,5)})
#=======================================================================================================
#Compute the total Jacobian of UAV+arm system
#NOTE: t4 does not affect the Jacobian for our robot actually
def Jacobian(t1,t2,t3,t4,psi):

    phi=0;theta=0   #underactuated parameters
    #params
    l1=0.0695; l2=0.17; l3=0.07025; l4=0.025
    alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0
    a1=0; a2=l2; a3=0; a4=0
    b1=l1; b2=0; b3=0; b4=l3+l4
    shift=Matrix([[0.11316],[0],[0]])
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
    
    #p_eb= position of ee w.r.t base in base's frame
    p_eb=a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec+shift

    #Jacobian of EE w.r.t base in base's frame
    ###CHANGED
    J_eb=Matrix([
    [-l2*c(t2)*s(t1)-(l3+l4)*s(t2+t3)*s(t1),    -l2*c(t1)*s(t2)+(l3+l4)*c(t1)*c(t2+t3)      ,(l3+l4)*c(t1)*c(t2+t3) ,0],
    [l2*c(t1)*c(t2)+(l3+l4)*s(t2+t3)*c(t1),     -l2*s(t1)*s(t2)+(l3+l4)*s(t1)*c(t2+t3)      ,(l3+l4)*s(t1)*c(t2+t3) ,0],
    [0,                                         l2*c(t2)+(l3+l4)*s(t2+t3)                   ,(l3+l4)*s(t2+t3)       ,0],
    [0,s(t1)    ,s(t1)    ,c(t1)*s(t2+t3)],
    [0,-c(t1)   ,-c(t1)   ,s(t1)*s(t2+t3)],
    [1,0        ,0        ,-c(t2+t3)     ]])
    
    #construct total jacobian such that: J [xb. yb. zb. psi. t1. t2. t3. t4.]^T =[d/dt(p_e) omega]^T
    F_left=eye(3)
    F_left=F_left.col_join(zeros(3,3))
    F_right=-CPM(Rb*Rs*p_eb)*Matrix([[0],[0],[1]])
    F_right=F_right.col_join(Matrix([[0],[0],[1]]))
    F=F_left.row_join(F_right)

    G_left=Rb*Rs
    G_left=G_left.col_join(zeros(3,3))
    G_right=zeros(3,3)
    G_right=G_right.col_join(Rb*Rs)
    G=G_left.row_join(G_right)
    G=G*J_eb
    J=F.row_join(G)

    return J.evalf()
#=======================================================================================================
#FK 
#@param joint posture [xb yb zb psi t1 t2 t3 t4]^T
#@return pose of EE for given q. [x y z psi theta phi]^T
def FK(q):
    x=q[0];y=q[1];z=q[2];psi=q[3];t1=q[4];t2=q[5];t3=q[6];t4=q[7]
    
    l1=0.0695; l2=0.17; l3=0.07025; l4=0.025
    alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0
    a1=0; a2=l2; a3=0; a4=0
    b1=l1; b2=0; b3=0; b4=l3+l4
    shift=Matrix([[0.11316],[0],[0]])
    theta=0;phi=0   #underactuated variables
    p_b=Matrix([[x],[y],[z]])    #base position
    
    #Rb
    Rb=Matrix([
    [cos(psi)*cos(theta), sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi), sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)], 
    [sin(psi)*cos(theta), sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi), -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)], 
    [-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)]])

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
    
    #position
    p=p_b+Rb*Rs*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec)
    #orientation
    Q_e=Rb*Rs*P4
    
    #compute orientation angles
    ori=euler_angles(Q_e)
    
    #complete position and orientation result
    p_result=p.col_join(ori)
    return p_result.evalf()
#=======================================================================================================    
#IK solver
#@param x_d: desired goal position [x y z]^T 
#@param z_d: desired goal orientation [psi theta phi]^T
#@param q_1: initial joint postures [xb yb zb psi t1 t2 t3 t4]^T
#@param T: required time to move from x_1 to x_d
#return: joint postures q that will generate x_d
def IK(x_d,z_d,q_1,Time):
    N=300
    q_k=q_1     #initialize q_k
    x_k=FK(q_1)[0:3,0:1] #initialize x_k
    z_k=FK(q_1)[3:6,0:1] #initialize z_k
    
    dt=Time*1.0/N  #integration step time #multiply 1.0 to avoid integer division error
    
    #Weighting matrices
    We=3.0*eye(3)   #main task
    Wv=0.1*eye(8)   #singularity avoidance
    Wc=3.0*eye(3)   #additional task
    C=2
    for k in range(1,N+1):  #k=1,2,3, ... ,N

        #calculate planned velocity
        x_k_dot=C*(x_d-x_k)/((N+1-k)*dt)
        z_k_dot=C*(z_d-z_k)/((N+1-k)*dt)
        
        #find joint rates that generate x_k_dot
        #use configuration control method
        Jac=Jacobian(q_k[4],q_k[5],q_k[6],q_k[7],q_k[3])  #compute the jacobian at iteration k. Jacobian(t1,t2,t3,t4,psi)
        Je=Jac[0:3,0:8] #for position-> upper half of total Jacobian-> Je qdot=xdot. x=pose
        Jc=(T(z_k[0],z_k[1],z_k[2]))**-1 *Jac[3:6,0:8]  #for orientation-> lower half of total Jacobian. Jc qdot=zdot. z=euler angles
        
        #joint limit
        t1=q_k[4]
        t2=q_k[5]
        t3=q_k[6]
        W_JL_calculated=W_JL.evalf(subs={W1:W_i(t1,t1_min,t1_max,tau1,100),W2:W_i(t2,t2_min,t2_max,tau2,100),W3:W_i(t3,t3_min,t3_max,tau3,100)})

        q_k_dot=(Je.transpose()*We*Je + Jc.transpose()*Wc*Jc+ Wv +W_JL_calculated)**-1 * (Je.transpose()*We*x_k_dot+ Jc.transpose()*Wc*z_k_dot)

        #Integrate to find the next q_k
        q_k=q_k+q_k_dot*dt
        
        #Compute the new pose
        FK_result=FK(q_k)
        x_k=FK_result[0:3,0:1]
        z_k=FK_result[3:6,0:1]
        
        #put data on the list for plotting
        t1_list.append(q_k[4])
        t2_list.append(q_k[5])
        t3_list.append(q_k[6])
        t4_list.append(q_k[7])
        
        x_b_list.append(q_k[0])
        y_b_list.append(q_k[1])
        z_b_list.append(q_k[2])
        psi_b_list.append(q_k[3])
        
        psi_list.append(z_k[0])
        theta_list.append(z_k[1])
        phi_list.append(z_k[2])
        
        x_list.append(x_k[0])
        y_list.append(x_k[1])
        z_list.append(x_k[2])
        
        psi_dot_list_d.append(z_k_dot[0])
        theta_dot_list_d.append(z_k_dot[1])
        phi_dot_list_d.append(z_k_dot[2])
        
        z_dot_actual=Jc*q_k_dot
        psi_dot_list_actual.append(z_dot_actual[0])
        theta_dot_list_actual.append(z_dot_actual[1])
        phi_dot_list_actual.append(z_dot_actual[2])

    return q_k  #return the final joint postures i.e. joint postures that will generate x_d, z_d
    
#=======================================================================================================        
#==========JOINT LIMITS
#t1_min=-pi/6
#t1_max=pi/6
#t2_min=-pi+pi
#t2_max=0+pi
#t3_min=0-pi/2
#t3_max=pi-pi/2

begin=time.time()
q_1=Matrix([[0],[0],[0],[0],[0],[1.5],[0],[0]])   #initial posture
pose_goal=Matrix([[1],[1],[0.4],[1],[2],[3]])   #pose goal

pose_goal[3:6,0:1]=euler_angles(R(pose_goal[3],pose_goal[4],pose_goal[5])).evalf() #TODO

q_computed=IK(pose_goal[0:3,0:1],pose_goal[3:6,0:1],q_1,100)    #compute posture with IK
print('q computed by IK algo is: ',q_computed)

pose_computed=FK(q_computed)
print('the q computed would result to pose: ',pose_computed)

pose_error=pose_goal-pose_computed  #error ( goal- IK solution)
print('the error in pose is: ', pose_error)

print('-------------aaaaaaaaaaaaaa')
print(x_list)
print('-------------aaaaaaaaaaaaaa')
print(y_list)
print('-------------aaaaaaaaaaaaaa')
print(z_list)
print('-------------aaaaaaaaaaaaaa')

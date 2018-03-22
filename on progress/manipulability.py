from sympy import *
c=cos
s=sin
import time
import matplotlib.pyplot as plt

import numpy as np
from scipy.linalg import svdvals    #for SVD (Singular value Decomposition)
from mpl_toolkits.mplot3d import Axes3D # 3d plot
#---------------------------------------------------------------------------------------------------------------------------#
#rotation matrix between frame 1 of the arm with frame xyz of the world
Rs=Matrix([
    [1,0,0],
    [0,-1,0],
    [0,0,-1]])
    
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
#---------------------------------------------------------------------------------------------------------------------------#    
#MANIPULABILITY
#manipulability 
#@param joint posture [xb yb zb psi t1 t2 t3 t4]^T
#@return the manipulability at the given posture
def m(q):
    
    #Jacobian matrix
    Jac=Jacobian(q[4],q[5],q[6],q[7],q[3])  #compute the jacobian at iteration k. Jacobian(t1,t2,t3,t4,psi)
    print(Jac*Jac.transpose())
    #manipulability
    manip=(Jac*Jac.transpose()).det()
    
    #TODO uncomment me
    if abs(manip)<5* 10**-15:   #for very small number there can be error that results the umber to be negative because limit of number accuracy when calculating
        return 0
        
    return sqrt(manip)
    
#---------------------------------------------------------------------------------------------------------------------------#
#MANIPULABILITY using SVD
#manipulability 
#@param joint posture [xb yb zb psi t1 t2 t3 t4]^T
#@return the manipulability at the given posture
def m_svd(q):
    
    #Jacobian matrix
    Jac=Jacobian(q[4],q[5],q[6],q[7],q[3])  #compute the jacobian at iteration k. Jacobian(t1,t2,t3,t4,psi)

    Jac_formatted = np.array([
                            [float(Jac[0,0]),float(Jac[0,1]),float(Jac[0,2]),float(Jac[0,3]),float(Jac[0,4]),float(Jac[0,5]),float(Jac[0,6]),float(Jac[0,7])],
                            [float(Jac[1,0]),float(Jac[1,1]),float(Jac[1,2]),float(Jac[1,3]),float(Jac[1,4]),float(Jac[1,5]),float(Jac[1,6]),float(Jac[1,7])],
                            [float(Jac[2,0]),float(Jac[2,1]),float(Jac[2,2]),float(Jac[2,3]),float(Jac[2,4]),float(Jac[2,5]),float(Jac[2,6]),float(Jac[2,7])],
                            [float(Jac[3,0]),float(Jac[3,1]),float(Jac[3,2]),float(Jac[3,3]),float(Jac[3,4]),float(Jac[3,5]),float(Jac[3,6]),float(Jac[3,7])],
                            [float(Jac[4,0]),float(Jac[4,1]),float(Jac[4,2]),float(Jac[4,3]),float(Jac[4,4]),float(Jac[4,5]),float(Jac[4,6]),float(Jac[4,7])],
                            [float(Jac[5,0]),float(Jac[5,1]),float(Jac[5,2]),float(Jac[5,3]),float(Jac[5,4]),float(Jac[5,5]),float(Jac[5,6]),float(Jac[5,7])]
                            ])
    vals=svdvals(Jac_formatted) # compute Singular values
    
    #compute the product
    manip=1
    for val in vals:
        manip=manip*val
    
    return manip
#---------------------------------------------------------------------------------------------------------------------------#
#manipulability of ARM only
#DETERMINANT method
def mm(t1,t2,t3):
    J_eb=Matrix([
    [-l2*c(t2)*s(t1)-(l3+l4)*s(t2+t3)*s(t1),    -l2*c(t1)*s(t2)+(l3+l4)*c(t1)*c(t2+t3)      ,(l3+l4)*c(t1)*c(t2+t3) ,0],
    [l2*c(t1)*c(t2)+(l3+l4)*s(t2+t3)*c(t1),     -l2*s(t1)*s(t2)+(l3+l4)*s(t1)*c(t2+t3)      ,(l3+l4)*s(t1)*c(t2+t3) ,0],
    [0,                                         l2*c(t2)+(l3+l4)*s(t2+t3)                   ,(l3+l4)*s(t2+t3)       ,0],
    [0,s(t1)    ,s(t1)    ,c(t1)*s(t2+t3)],
    [0,-c(t1)   ,-c(t1)   ,s(t1)*s(t2+t3)],
    [1,0        ,0        ,-c(t2+t3)     ]])
    return sqrt((J_eb*J_eb.transpose()).det()) 
#---------------------------------------------------------------------------------------------------------------------------#  
def vary_t1():
    t1_list=np.arange(0,2*pi+0.1,0.1)
    m_t1=[]
    for t1 in t1_list:
        q_test=Matrix([[1],[100],[1],[2],[t1],[1+pi],[0-pi/2],[0]]) 
        man=m(q_test)
        m_t1.append(man)
 
    plt.figure(1)   
    plt.ylim(float(min(m_t1))-0.1,float(max(m_t1))+0.1)
    plt.plot(t1_list,m_t1)
    plt.show()

#=========================
def vary_t2():
    t2_list=np.arange(0,2*pi+0.1,0.1)
    m_t2=[]
    for t2 in t2_list:
        q_test=Matrix([[1],[100],[1],[2],[0],[t2],[0-pi/2],[0]]) 
        man=m(q_test)
        m_t2.append(man)
    plt.figure(1)   
    plt.ylim(float(min(m_t2))-0.1,float(max(m_t2))+0.1)
    plt.plot(t2_list,m_t2)
    plt.show()  

#=========================
def vary_t3():
    t3_list=np.arange(0,2*pi+0.1,0.1)
    m_t3=[]
    for t3 in t3_list:
        
        #man const=0
  #      q_test=Matrix([ [ 0.970961285262131],[  1.02132778481411],[ 0.694669112087531],[ 0.679252605172343],[-0.323322198466938[  2.27175208891396],[ t3],[  1.00257187369481]])
        q_test=Matrix([[1],[100],[1],[2],[0],[1+pi],[t3],[0]]) 
        man=m(q_test)
        m_t3.append(man)
    plt.figure(1)   
    plt.ylim(float(min(m_t3))-0.1,float(max(m_t3))+0.1)
    plt.plot(t3_list,m_t3)
    plt.show()
#=========================    
def vary_t2t3():
    t2_list = np.array([
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),
    np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100),np.linspace(0, 2 * np.pi, 100)
    ])
    
    t3_list=np.transpose(t2_list)

    j=0   
    man_list=np.zeros((100,100)) #TODO CHANGE size
    #t2_list, t3_list = np.meshgrid(t2_list, t3_list)
    for i in range(len(t2_list)):
        for k in range(len(t2_list[0])):
            print(i,k)
            #print(j)
            #j=j+1
            q_test=Matrix([[1],[100],[1],[2],[0],[t2_list[i,k]],[t3_list[i,k]],[0]])
            man=m(q_test)
            print(man)
            man_list[i,k]=man
        #man#_list[i]=m_child
    #print(man_list)
    print(man_list)
    fig = plt.figure()
    ax = fig.gca(projection='3d')  
    ax.plot_surface(t2_list, t3_list, man_list, color='b')
    ax.set_xlabel('t2')
    ax.set_ylabel('t3')
    ax.set_zlabel('manipulability')
    plt.show()
#********************************************************************************************************************************#
#vary_t2()
vary_t3()
#vary_t2t3()
#vary_t2t3()


#begin=time.time()
#for i in range(10):
#    q_test=Matrix([[1],[100],[1],[2],[2.2],[0+pi],[0-pi/2],[0]]) 
#    print(m(q_test))
#    #print(m_svd(q_test))
#end=time.time()
#print(end-begin)

from sympy import *
c=cos
s=sin

######################################################################################################################################
#Computing rotation matrix Q
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
                   
#Computing vector of a matrix of 3x3 matrix
def vect(A):
    return 0.5*Matrix([ [A[7]-A[5]],
                    [A[2]-A[6]],
                    [A[3]-A[1]]])
#compute ZYZ euler angles given a 3x3rotation matrix             
def euler_angles(R):
    if(R[8]<1):    #if(R[2][2]<1):
        if(R[8]>-1):
      
            theta=acos(R[8])
            phi=atan2(R[5],R[2])#atan2(R[1][2],R[0][2]) #CHANGE
            psi=atan2(R[7],-R[6])#atan2(R[2][1],-R[2][0])
        elif(R[8]==-1):

            theta=pi
            phi=-atan2(R[3],R[4])#-atan2(R[1][0],R[1][1])
            psi=0
    else:

        theta=0
        phi=atan2(R[3],R[4])#atan2(R[1][0],R[1][1])
        psi=0
    return Matrix([[phi],[theta],[psi]])
    #https://www.geometrictools.com/Documentation/EulerAngles.pdf
  
##################################################################################################################################
#variables 
x,y,z,theta,phi,psi,t1,t2,t3,t4= symbols('x y z theta phi psi t1 t2 t3 t4')
l1,l2,l3,l4,alpha1,alpha2,alpha3,alpha4=symbols('l1 l2 l3 l4 alpha1 alpha2 alpha3 alpha4')  #known from system
xshift, yshift, zshift=symbols('xshift yshift zshift')                                      #known from system
a1,a2,a3,a4,b1,b2,b3,b4=symbols('a1 a2 a3 a4 b1 b2 b3 b4')                                  
xd,yd,zd,thetad,psid,phid=symbols('xd yd zd thetad psid phid')                              #goal pose given
rho1,rho2=symbols('rho1 rho2')  #redundancy parameters
shift=symbols('shift')

l1=0.0695; l2=0.17; l3=0.07025; l4=0.025
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0
xshift=0.11316; yshift=0; zshift=0
a1=0; a2=l2; a3=0; a4=0
b1=l1; b2=0; b3=0; b4=l3+l4
shift=Matrix([[xshift],[yshift],[zshift]])

#underactuated parameters of UAV
theta=0
psi=0
############## specify GOAL POSE ###################################
xd=10
yd=2
zd=3
thetad=1
psid=1
phid=1
####################################################################

################# redundancy parameter##############################

####################################################################
p_b=Matrix([[x],[y],[z]])    #base position

#Rb
Rb= Matrix([[c(phi)*c(theta)*c(psi)-s(phi)*s(psi),  -c(phi)*c(theta)*s(psi)-s(phi)*s(psi),  c(phi)*s(theta)],
            [s(phi)*c(theta)*c(psi)+c(phi)*s(psi),  -s(phi)*c(theta)*s(psi)+c(phi)*s(psi),  s(phi)*s(theta)],
            [-s(phi)*c(psi),                        s(phi)*s(psi),                          c(theta)        ]])
            
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

#JACOBIAN###############################
#p_eb= position of ee w.r.t base in base's frame
p_eb=a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec+shift

#T : matrix relating (rate of  euler angles) to  angular velocity
T=Matrix([
        [0,-s(phi),c(phi)*c(theta)],
        [0,c(phi), s(phi)*s(theta)],
        [1,0,c(theta)]])
        
#Jacobian of EE w.r.t base in base's frame
J_eb=Matrix([
    [0,s(t1)    ,s(t1)    ,c(t1)*s(t2+t3)],
    [0,-c(t1)   ,-c(t1)   ,s(t1)*s(t2+t3)],
    [1,0        ,0        ,-c(t2+t3)     ],
    [-l2*c(t2)*s(t1)-(l3+l4)*s(t2+t3)*s(t1),    -l2*c(t1)*s(t2)+(l3+l4)*c(t1)*c(t2+t3)      ,(l3+l4)*c(t1)*c(t2+t3) ,0],
    [l2*c(t1)*c(t2)+(l3+l4)*s(t2+t3)*c(t1),     -l2*s(t1)*s(t2)+(l3+l4)*s(t1)*c(t2+t3)      ,(l3+l4)*s(t1)*c(t2+t3) ,0],
    [0,                                         l2*c(t2)+(l3+l4)*s(t2+t3)                   ,(l3+l4)*s(t2+t3)       ,0]])
    
F_left=eye(3)
F_left=F_left.col_join(zeros(3,3))
F_right=-CPM(Rb*p_eb)*T
F_right=F_right.col_join(T)

#TODO slice the last 2 columns of F
F_right=Matrix([[F_right[0]],[F_right[3]],[F_right[6]],[F_right[9]],[F_right[12]],[F_right[15]]])

F=F_left.row_join(F_right)



G_left=Rb
G_left=G_left.col_join(zeros(3,3))
G_right=zeros(3,3)
G_right=G_right.col_join(Rb)
G=G_left.row_join(G_right)
G=G*J_eb
Jacobian=F.row_join(G)

def J(t1,t2,t3,t4,phi):
    return Matrix([[1, 0, 0, -(0.09525*sin(t1)*sin(t2)*cos(t3) + 0.09525*sin(t1)*sin(t3)*cos(t2) + 0.17*sin(t1)*cos(t2))*cos(phi) - (0.09525*sin(t2)*cos(t1)*cos(t3) + 0.09525*sin(t3)*cos(t1)*cos(t2) + 0.17*cos(t1)*cos(t2) + 0.11316)*sin(phi), cos(phi), sin(t1)*cos(phi), sin(t1)*cos(phi), sin(t2 + t3)*cos(phi)*cos(t1) - cos(phi)*cos(t2 + t3)], [0, 1, 0, (0.09525*sin(t2)*sin(t3) + 0.17*sin(t2) - 0.09525*cos(t2)*cos(t3) + 0.0695)*cos(phi) + (0.09525*sin(t2)*cos(t1)*cos(t3) + 0.09525*sin(t3)*cos(t1)*cos(t2) + 0.17*cos(t1)*cos(t2) + 0.11316)*cos(phi), 0, sin(phi)*sin(t1) - cos(phi)*cos(t1), sin(phi)*sin(t1) - cos(phi)*cos(t1), sin(phi)*sin(t2 + t3)*cos(t1) + sin(t1)*sin(t2 + t3)*cos(phi)], [0, 0, 1, 0, 1, -sin(phi)*sin(t1), -sin(phi)*sin(t1), -sin(phi)*sin(t2 + t3)*cos(t1) - cos(t2 + t3)], [0, 0, 0, 0, (-0.09525*sin(t1)*sin(t2 + t3) - 0.17*sin(t1)*cos(t2))*cos(phi), (-0.17*sin(t2)*cos(t1) + 0.09525*cos(t1)*cos(t2 + t3))*cos(phi) + (0.09525*sin(t2 + t3) + 0.17*cos(t2))*cos(phi), 0.09525*sin(t2 + t3)*cos(phi) + 0.09525*cos(phi)*cos(t1)*cos(t2 + t3), 0], [0, 0, 0, 0, (-0.09525*sin(t1)*sin(t2 + t3) - 0.17*sin(t1)*cos(t2))*sin(phi) + (0.09525*sin(t2 + t3)*cos(t1) + 0.17*cos(t1)*cos(t2))*cos(phi), (-0.17*sin(t1)*sin(t2) + 0.09525*sin(t1)*cos(t2 + t3))*cos(phi) + (-0.17*sin(t2)*cos(t1) + 0.09525*cos(t1)*cos(t2 + t3))*sin(phi), 0.09525*sin(phi)*cos(t1)*cos(t2 + t3) + 0.09525*sin(t1)*cos(phi)*cos(t2 + t3), 0], [0, 0, 0, 1, -(-0.09525*sin(t1)*sin(t2 + t3) - 0.17*sin(t1)*cos(t2))*sin(phi), -(-0.17*sin(t2)*cos(t1) + 0.09525*cos(t1)*cos(t2 + t3))*sin(phi) + 0.09525*sin(t2 + t3) + 0.17*cos(t2), -0.09525*sin(phi)*cos(t1)*cos(t2 + t3) + 0.09525*sin(t2 + t3), 0]])
####################################################################################################################
#q=x y z phi t1 t2 t3 t4
def FK(q):
    x=q[0];y=q[1];z=q[2];phi=q[3];t1=q[4];t2=q[5];t3=q[6];t4=q[7]
    xshift=0.11316; yshift=0; zshift=0
    shift=Matrix([[xshift],[yshift],[zshift]])
    theta=0;psi=0   #underactuated variables
    p_b=Matrix([[x],[y],[z]])    #base position
    #Rb
    Rb= Matrix([[c(phi)*c(theta)*c(psi)-s(phi)*s(psi),  -c(phi)*c(theta)*s(psi)-s(phi)*c(psi),  c(phi)*s(theta)],
            [s(phi)*c(theta)*c(psi)+c(phi)*s(psi),  -s(phi)*c(theta)*s(psi)+c(phi)*s(psi),  s(phi)*s(theta)],
            [-s(phi)*c(psi),                        s(phi)*s(psi),                          c(theta)        ]])
            
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
    p=p_b+Rb*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec)
    #orientation
    Q_e=Rb*P4
    q_e=vect(Q_e)
    #compute orientation angles
    ori=euler_angles(Q_e) #TODO
    
    #complete position and orientation result
    p_result=p.col_join(ori)
    return p_result
###################################################################################3
    
#Weight matrices
We= 1.*eye(6)
Ws= 0.01*eye(8)

T=2     #time required to reach goal
n=100    #number of iterations
xi=0;yi=0;zi=0;phii=0;t1i=0;t2i=0;t3i=0;t4i=1
q_current=Matrix([[xi],[yi],[zi],[phii],[t1i],[t2i],[t3i],[t4i]])    #initial q
p_goal=Matrix([[xd], [yd], [zd], [thetad], [psid], [phid]])
p_current=FK(q_current) #6x1
for k in range(n):

    #p_dot
    p_dot=(p_goal-p_current)/(T-k*T/n)  #6x1

    #q_dot
    t1=q_current[4];t2=q_current[5];t3=q_current[6];t4=q_current[7];phi=q_current[4];
    Jac=J(t1,t2,t3,t4,phi).evalf()

    det=(Jac.transpose()*We*Jac+Ws).det()
    q_dot=det*(Jac.transpose()*We*p_dot).evalf()    #8x1

    #update q_current
    q_current=q_current+q_dot*(T*1.0/n)

    #p_current
    p_current=FK(q_current) #6x1
    

    
print(q_current)
print(FK(q_current).evalf())
#JJt=J*J.transpose()
#print(JJt)   

############################################################################################
#man=JJt.det()
#print(man)

#EQUATIONS
###position###
#position_eqn=p_goal-p_b-Rb*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec)

###orientation###
#desired orientation
#Q_goal=Matrix([[c(phid)*c(thetad)*c(psid)-s(phid)*s(psid),  -c(phid)*c(thetad)*s(psid)-s(phid)*s(psid),  c(phid)*c(thetad)],
#           [s(phid)*c(thetad)*c(psid)+c(phid)*s(psid),  -s(phid)*c(thetad)*s(psid)+c(phid)*c(psid),  s(phid)*s(thetad)],
#            [-s(phid)*c(psid),                        s(phid)*s(psid),                          c(thetad)        ]])
#q_goal=vect(Q_goal)
#ee true orientation
#Q_e=Rb*P4
#q_e=vect(Q_e)

#orientation_eqn=q_goal-q_e

#variables to solve
#var=[x,y,z,t1,t3,t4]    #phi,t2-> redundancy
#eqns=[position_eqn[0],position_eqn[1],position_eqn[2],orientation_eqn[0],orientation_eqn[1],orientation_eqn[2]]
#var=[t1,t3]    #phi,t2-> redundancy
#eqns=[orientation_eqn[0],orientation_eqn[1]]
#print(solve(eqns,var))





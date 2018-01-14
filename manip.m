function [ ] =manip()

syms x y z theta phi psi t1 t2 t3 t4
l1=0.0695; l2=0.17; l3=0.07025; l4=0.025
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0
a1=0; a2=l2; a3=0; a4=0
b1=l1; b2=0; b3=0; b4=l3+l4
theta=0;psi=0   %underactuated params of UAV
shift=[0.11316;0;0]

%redundancy
t1=0

%Q1, Q2, Q3, Q4
Q1=Q(alpha1,t1)
Q2=Q(alpha2,t2)
Q3=Q(alpha3,t3)
Q4=Q(alpha4,t4)
%P1, P2, P3, P4
P1=Q1
P2=P1*Q2
P3=P2*Q3
P4=P3*Q4

T=[0,-sin(phi),cos(phi)*cos(theta);
0,cos(phi), sin(phi)*sin(theta);
1,0,cos(theta)]

%vectors: a1_vec, a2_vec, a3_vec, a4_vec
a1_vec=a_vector(a1,b1,t1)
a2_vec=a_vector(a2,b2,t2)
a3_vec=a_vector(a3,b3,t3)
a4_vec=a_vector(a4,b4,t4)

%p_eb= position of ee w.r.t base in base's frame
p_eb=a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec+shift

%Rb rotation matrix of base frame
Rb= [cos(phi)*cos(theta)*cos(psi)-sin(phi)*sin(psi),    -cos(phi)*cos(theta)*sin(psi)-sin(phi)*sin(psi),    cos(phi)*cos(theta);
    sin(phi)*cos(theta)*cos(psi)+cos(phi)*sin(psi),     -sin(phi)*cos(theta)*sin(psi)+cos(phi)*cos(psi),    sin(phi)*sin(theta);
    -sin(phi)*cos(psi),                                 sin(phi)*sin(psi),                                  cos(theta)]
    
J_eb=[0,sin(t1)    ,sin(t1)    ,cos(t1)*sin(t2+t3);
0,-cos(t1)   ,-cos(t1)   ,sin(t1)*sin(t2+t3);
1,0        ,0        ,-cos(t2+t3);
-l2*cos(t2)*sin(t1)-(l3+l4)*sin(t2+t3)*sin(t1),    -l2*cos(t1)*sin(t2)+(l3+l4)*cos(t1)*cos(t2+t3)      ,(l3+l4)*cos(t1)*cos(t2+t3) ,0;
l2*cos(t1)*cos(t2)+(l3+l4)*sin(t2+t3)*cos(t1),     -l2*sin(t1)*sin(t2)+(l3+l4)*sin(t1)*cos(t2+t3)      ,(l3+l4)*sin(t1)*cos(t2+t3) ,0;
0,                                         l2*cos(t2)+(l3+l4)*sin(t2+t3)                   ,(l3+l4)*sin(t2+t3)       ,0]

F_left=eye(3)
F_left=vertcat(F_left,zeros(3,3))
F_right=-CPM(Rb*p_eb)*T
F_right=vertcat(F_right,T)
F=horzcat(F_left,F_right)
G_left=Rb
G_left=vertcat(G_left,zeros(3,3))
G_right=zeros(3,3)
G_right=vertcat(G_right,Rb)
G=horzcat(G_left,G_right)
G=G*J_eb
J=horzcat(F,G)
JJt=J*transpose(J) 
solve(diff(det(JJt) ,t2),t2)
end

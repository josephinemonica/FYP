function [ q ] =IK_solver(x,y,z,psi,t1,t2,t3,t4)
%constant shift from UAV to base
shift=[0.11316;0;0];
%xshift=0.11316; yshift=0; zshift=0

%DH constant parameters
%syms l1 l2 l3 l4 alpha1 alpha2 alpha3 alpha4

%DH params
%syms a1 a2 a3 a4 b1 b2 b3 b4

l1=0.0695; l2=0.17; l3=0.07025; l4=0.025;
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0;
a1=0; a2=l2; a3=0; a4=0;
b1=l1; b2=0; b3=0; b4=l3+l4;
theta=0;phi=0;   %underactuated params of UAV

p_b=[x;y;z];         %base position


%Rb rotation matrix of base frame
Rb=[cos(psi)*cos(theta), sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi), sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi),
sin(psi)*cos(theta), sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi), -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi), 
-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)] ; 
    
%Q1, Q2, Q3, Q4
Q1=Q(alpha1,t1);
Q2=Q(alpha2,t2);
Q3=Q(alpha3,t3);
Q4=Q(alpha4,t4);

%P1, P2, P3, P4
P1=Q1;
P2=P1*Q2;
P3=P2*Q3;
P4=P3*Q4;

%vectors: a1_vec, a2_vec, a3_vec, a4_vec
a1_vec=a_vector(a1,b1,t1);
a2_vec=a_vector(a2,b2,t2);
a3_vec=a_vector(a3,b3,t3);
a4_vec=a_vector(a4,b4,t4);

%position equations
p_result=p_b+Rb*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec);
%orientation equations
Q_e=Rb*P4;   %ee true orientation
disp(vpa(Q_e))
ori=euler_angles(Q_e);
disp(vpa(p_result))
disp(vpa(ori))

end


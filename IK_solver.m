function [ q ] =IK_solver(xd,yd,zd,psid,thetad,phid)
t2=0.1;psi=0.1;
%actuating parameters of UAV and arm
syms x y z theta phi t1 t3 t4;

%constant shift from UAV to base
shift=[0.11316;0;0];

%DH params
%syms a1 a2 a3 a4 b1 b2 b3 b4
%syms l1 l2 l3 l4 alpha1 alpha2 alpha3 alpha4
l1=0.0695; l2=0.17; l3=0.07025; l4=0.025;
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0;
a1=0; a2=l2; a3=0; a4=0;
b1=l1; b2=0; b3=0; b4=l3+l4;
theta=0;phi=0;   %underactuated params of UAV

%Solve IK of a system
%input: Goal pose
%---------------------------------------------------------------------------------------------
p_b=[x;y;z];         %base position
p_goal=[xd;yd;zd];   %goal position, GIVEN
                    %goal orientation, GIVEN
Q_goal=[cos(psid)*cos(thetad), sin(phid)*sin(thetad)*cos(psid) - sin(psid)*cos(phid), sin(phid)*sin(psid) + sin(thetad)*cos(phid)*cos(psid),
sin(psid)*cos(thetad), sin(phid)*sin(psid)*sin(thetad) + cos(phid)*cos(psid), -sin(phid)*cos(psid) + sin(psid)*sin(thetad)*cos(phid), 
-sin(thetad), sin(phid)*cos(thetad), cos(phid)*cos(thetad)] ; 
q_goal=vect(Q_goal);

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
position_eqn=p_goal-p_b-Rb*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec);
%orientation equations
Q_e=Rb*P4;   %ee true orientation
q_e=vect(Q_e);
orientation_eqn=q_goal-q_e;

%-------------------------------------------------------------------------------
p_eb=a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec+shift;
J_eb=[0,sin(t1)    ,sin(t1)    ,cos(t1)*sin(t2+t3);
0,-cos(t1)   ,-cos(t1)   ,sin(t1)*sin(t2+t3);
1,0        ,0        ,-cos(t2+t3);
-l2*cos(t2)*sin(t1)-(l3+l4)*sin(t2+t3)*sin(t1),    -l2*cos(t1)*sin(t2)+(l3+l4)*cos(t1)*cos(t2+t3)      ,(l3+l4)*cos(t1)*cos(t2+t3) ,0;
l2*cos(t1)*cos(t2)+(l3+l4)*sin(t2+t3)*cos(t1),     -l2*sin(t1)*sin(t2)+(l3+l4)*sin(t1)*cos(t2+t3)      ,(l3+l4)*sin(t1)*cos(t2+t3) ,0;
0,                                         l2*cos(t2)+(l3+l4)*sin(t2+t3)                   ,(l3+l4)*sin(t2+t3)       ,0];

F_left=eye(3);
F_left=vertcat(F_left,zeros(3,3));
F_right=-CPM(Rb*p_eb)*[0;0;1];
F_right=vertcat(F_right,[0;0;1]);
F=horzcat(F_left,F_right);
G_left=Rb;
G_left=vertcat(G_left,zeros(3,3));
G_right=zeros(3,3);
G_right=vertcat(G_right,Rb);
G=horzcat(G_left,G_right);
G=G*J_eb;
J=horzcat(F,G);

%---------------------------------------------------------------------------
%SOLVE


var=[x y z t1 t3 t4];
eqns=[position_eqn(1),position_eqn(2),position_eqn(3),orientation_eqn(1),orientation_eqn(2),orientation_eqn(3)];
S=vpasolve(eqns,var);



disp(S.x)
disp(S.y)
disp(S.z)
disp(psi)
disp(vpa(mod(S.t1,2*pi)))
disp(t2)
disp(vpa(mod(S.t3,2*pi)))
disp(vpa(mod(S.t4,2*pi)))

Q_e=Rb*P4;   %ee true orientation

x=S.x;
y=S.y;
z=S.z;
t1=mod(S.t1,2*pi);
t3=mod(S.t3,2*pi);
t4=mod(S.t4,2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
position_eqn=p_goal-p_b-Rb*(shift+a1_vec+P1*a2_vec+P2*a3_vec+P3*a4_vec);
%orientation equations
Q_e=Rb*P4;   %ee true orientation
disp(vpa((Q_goal)))
disp(vpa((Q_e)))
end


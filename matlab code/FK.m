function [ y ] = FK( pos )

%Assign values of the pose accordingly to each variable
x=pos(1,1);
y=pos(2,1);
z=pos(3,1);
phi=pos(4,1);
theta=pos(5,1);
psi=pos(6,1);
theta1=pos(7,1);
theta2=pos(8,1);
theta3=pos(9,1);
theta4=pos(10,1);
l1=         0.04;            %%
l2=0.17;
l3=0.07025;
l4=0.025;
xb=   0.1;   %%
zb=   0.02;    %%

%Values that have a '%' sign on the right are not the exact values


%Transformation matrix from one frame to another
T4_e = [1 0 0 0; 
    0 1 0 0;
    0 0 1 l4; 
    0 0 0 1];
T4_3 = [cos(theta4) -sin(theta4) 0 0;
    0 0 -1 -l3;
    sin(theta4) cos(theta4) 0 0;
    0 0 0 1];
T2_3 = [-sin(theta3) -cos(theta3) 0 l2;
    cos(theta3) -sin(theta3) 0 0;
    0 0 1 0;
    0 0 0 1];
T1_2 = [cos(theta2) -sin(theta2) 0 0;
    0 0 -1 0;
    sin(theta2) cos(theta2) 0 0;
    0 0 0 1];
T0_1 = [cos(theta1) -sin(theta1) 0 0;
    sin(theta1) cos(theta1) 0 0;
    0 0 1 l1;
    0 0 0 1];
Tb_0 = [1 0 0 xb;
    0 1 0 0;
    0 0 1 zb;
    0 0 0 1];
%Transformation between the end effector and the base
Tb_e = Tb_0 * T0_1 * T1_2 * T2_3 * T4_3;

Rb = [ cos(psi)*cos(theta) sin(phi)*sin(theta)*cos(psi)-sin(psi)*cos(phi) sin(psi)*sin(phi)+cos(psi)*sin(theta)*cos(phi);
    sin(psi)*cos(theta) cos(psi)*cos(phi)+sin(psi)*sin(theta)*sin(phi) sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
    -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];

Tb(1:3,1:3)=Rb(1:3,1:3);
Tb(1:4,4) = [x;y;z;1];
Tb(4,1:3) = [0 0 0];

%Transformation between the end effector and the world frane
T = Tb*Tb_e;

R = T(1:3,1:3);

%Position of the end effector
P = T(1:3,4);

y(1:3,1)= P;

%Orientation of the end effector in X,Y,Z
eulZYX = rotm2eul(R);
y(4,1)= eulZYX(1,3) ;
y(5,1)= eulZYX(1,2) ;
y(6,1)= eulZYX(1,1) ;


end


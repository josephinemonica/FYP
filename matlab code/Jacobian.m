function [ J ] = Jacobian( a )
 
%Values that have a '%' sign on the right are not the exact values
xb= 0.1;  %%
zb=  0.02; %%
l1= 0.04;  %%
l2=0.17;
l3=0.07025;
l4=0.025;
phi=a(1,1);
theta=a(2,1);
psi=a(3,1);
theta1=a(4,1);
theta2=a(5,1);
theta3=a(6,1);
theta4=a(7,1);

%Jacobian of the end effector wrt the base
Jb_eb = [-l3*sin(theta1)*cos(theta2+theta3)-l4*sin(theta1)*cos(theta2+theta3)-l2*sin(theta1)*cos(theta2) -l3*cos(theta1)*sin(theta2+theta3)-l4*cos(theta1)*sin(theta2+theta3)-l2*sin(theta2)*cos(theta1) -l3*cos(theta1)*sin(theta2+theta3)-l4*cos(theta1)*sin(theta2+theta3) 0;
    l3*cos(theta1)*cos(theta2+theta3)+l4*cos(theta1)*cos(theta2+theta3)+l2*cos(theta1)*cos(theta2) -l3*sin(theta1)*sin(theta2+theta3)-l4*sin(theta1)*sin(theta2+theta3)-l2*sin(theta1)*sin(theta2) -l3*sin(theta1)*sin(theta2+theta3)-l4*sin(theta1)*sin(theta2+theta3) 0;
    0 l3*cos(theta2+theta3)+l4*cos(theta2+theta3)+l2*cos(theta2) l3*cos(theta2+theta3)+l4*cos(theta2+theta3)+l2*cos(theta2) 0;
    0 sin(theta1) sin(theta1) cos(theta1)*cos(theta2+theta3);
    0 -cos(theta1) -cos(theta1) sin(theta1)*cos(theta2+theta3);
    1 0 0 sin(theta2+theta3)];

%Rotation matrix between the UAV and the world frame
Rb = [ cos(psi)*cos(theta) sin(phi)*sin(theta)*cos(psi)-sin(psi)*cos(phi) sin(psi)*sin(phi)+cos(psi)*sin(theta)*cos(phi);
    sin(psi)*cos(theta) cos(psi)*cos(phi)+sin(psi)*sin(theta)*sin(phi) sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
    -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];

%Jacobian of the end effector wrt the world frame
Jeb_1(1:3,1:3) = Rb(1:3,1:3);
Jeb_1(1:3,4:6) = zeros(3);
Jeb_1(4:6,4:6) = Rb(1:3,1:3);
Jeb_1(4:6,1:3) = zeros(3);
Jeb = Jeb_1 * Jb_eb;

%Position of the end effector wrt the base
pb_eb = [xb - l3*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + l2*cos(theta1)*cos(theta2);
      l2*cos(theta2)*sin(theta1) - l3*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1));
      l1 + zb + l3*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)) + l2*sin(theta2)];
Rbpb = Rb*pb_eb;

%Skew matrix of Rbpb
eve = [0 -Rbpb(3) Rbpb(2);
    Rbpb(3) 0 -Rbpb(1);
    -Rbpb(2) Rbpb(1) 0];

%Jacobian of the UAV
Jb(1:3,4:6) = eve(1:3,1:3);
Jb(1:3,1:3) = eye(3);
Jb(4:6,4:6) = eye(3);
Jb(4:6,1:3) = zeros(3);

%Tb - Transformation between angular velocity and time derivative of Euler
%angles

Tb = [0 -sin(psi) cos(psi)*cos(theta);
    0 cos(psi) sin(psi)*cos(theta);
    1 0 -sin(theta)];
Qb(4:6,4:6) = Tb(1:3,1:3);
Qb(1:3,1:3) = eye(3);
Qb(1:3,4:6) = zeros(3);
Qb(4:6,1:3) = zeros(3);

%Jacobian of the UAV (angle in terms of time derivative of Euler angles)
J_before=Jb*Qb;
J(1:6,1:3) = J_before(1:6,1:3) ;
J(1:6,4) = J_before(1:6,6) ;

%The complete Jacobian
J(1:6,5:8)=Jeb(1:6,1:4);


end


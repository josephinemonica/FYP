%Parametes of the system

%actuating parameters of UAV and arm
syms x y z theta phi psi t1 t2 t3 t4

%constant shift from UAV to base
syms shift
syms xshift yshift zshift
shift=[xshift;yshift;zshift]

%DH constant parameters
syms l1 l2 l3 l4 alpha1 alpha2 alpha3 alpha4

%DH params
syms a1 a2 a3 a4 b1 b2 b3 b4

l1=0.0695; l2=0.17; l3=0.07025; l4=0.025
alpha1=pi/2; alpha2=0; alpha3=pi/2; alpha4=0
xshift=0.11316; yshift=0; zshift=0
a1=0; a2=l2; a3=0; a4=0
b1=l1; b2=0; b3=0; b4=l3+l4
theta=0;psi=0   %underactuated params of UAV

%redundancy params
phi=0;t2=0


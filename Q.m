function [ Q ] = Q( alpha,theta )
%compute rotation matrix
%alpha, theta : DH parameters
mu=sin(alpha);
lam=cos(alpha);
Q=[ cos(theta),-lam*sin(theta),mu*sin(theta);
    sin(theta),lam*cos(theta),-mu*cos(theta);
    0,mu,lam];
end


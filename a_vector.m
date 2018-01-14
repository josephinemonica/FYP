function [ a_vector ] = a_vector( a,b,theta)
%Compute a_vector, vector from O_i to O_i+1
a_vector=[  a*cos(theta);
            a*sin(theta);
            b];
end


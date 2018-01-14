function [ vect ] = vect(A)
%Compute vector of a matrix of 3x3 matrix
vect=0.5*[A(3,2)-A(2,3);
        A(1,3)-A(3,1);
        A(2,1)-A(1,2)];
end


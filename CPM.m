function [ CPM ] = CPM(v)
%Compute cross product matrix of 3x1 vector
CPM=[0, -v(3),  v(2);
    v(3),   0,  -v(1);
    -v(2),  v(1),   0];
end


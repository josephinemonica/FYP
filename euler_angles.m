function [ angles ] =euler_angles(R)
%compute ZYX euler angles given a 3x3rotation matrix   
if (R(3,1)<1)
    if(R(3,1)>-1)
        theta=asin(-R(3,1));
        psi=atan2(R(2,1),R(1,1));
        phi=atan2(R(3,2),R(3,3));
    else
        theta=pi/2;
        psi=-atan2(-R(2,3),R(2,2));
        phi=0;
    end
else
    theta=-pi/2;
    psi=atan2(-R(2,3),R(2,2));
    phi=0;
end
angles= [psi;theta;phi];

end


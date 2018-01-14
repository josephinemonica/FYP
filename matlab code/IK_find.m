%All lengths are in meters and angles are in radians
x = [];
desired = [];

%Get the desired end effector pose
prompt = 'Enter the desired pose of the end effector:\n';
for i=1:6
    desired(i,1) = input(prompt); 
end

%Get the time duration in which the pose is expected to be achieved
prompt = 'Enter the time duration in which the pose must be reached:\n';
t=input(prompt);

%The initial pose of the UAV and the base are recorded
prompt = 'Enter the initial base position coordinates and joint values one by one:\n';
for i=1:10
    x(i,1) = input(prompt); 
end

for k=1:10

%Compute the position of the end effector (Forward Kinematics)
y = FK(x);

%Find the velocity of the end effector
%fprintf('aaa');
vel = (10/(t*(11-k))).*(desired-y);
%disp(vel)
%fprintf('bbb')
%Weight matrices
We= 1.*eye(6);
Ws= 0.1.*eye(8);


Je=Jacobian(x(4:10,1));

Jacobian1=Je';
fprintf('aaa');
disp(Jacobian1)
fprintf('bbb')

fordet = det( (Jacobian1*We*Je) +  Ws);


oldpos = x(1:3,1);
oldpos(4:8,1) = x(6:10,1);

newvel = (Jacobian1*We*vel);
%Velocity of the joints
newvel = fordet.*newvel;

%Numerical integration to find the new position
pos = (t/10).*newvel + oldpos;

x(1:4,1)=pos(1:4,1);
x(7:10,1)=pos(5:8,1);
end
disp('The required base coordinates and joint angles are:');
disp(pos);

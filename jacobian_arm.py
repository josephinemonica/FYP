#Compute Jacobian matrix of the UAV+arm system
import numpy as np

import math
c=math.cos
s=math.sin

#System parameters
xb=0.11316
l1=0.0695
l2=0.17
l3=0.07025
l4=0.025
t1=1
t2=2
t3=3
t4=4
#angles of base
theta=0
phi=0
psi=0
###################
def jacobian(t1,t2,t3,t4,l1,l2,l3,l4):

    J=np.matrix([
    [0,s(t1)    ,s(t1)    ,c(t1)*s(t2+t3)],
    [0,-c(t1)   ,-c(t1)   ,s(t1)*s(t2+t3)],
    [1,0        ,0        ,-c(t2+t3)     ],
    [-l2*c(t2)*s(t1)-(l3+l4)*s(t2+t3)*s(t1),    -l2*c(t1)*s(t2)+(l3+l4)*c(t1)*c(t2+t3)      ,(l3+l4)*c(t1)*c(t2+t3) ,0],
    [l2*c(t1)*c(t2)+(l3+l4)*s(t2+t3)*c(t1),     -l2*s(t1)*s(t2)+(l3+l4)*s(t1)*c(t2+t3)      ,(l3+l4)*s(t1)*c(t2+t3) ,0],
    [0,                                         l2*c(t2)+(l3+l4)*s(t2+t3)                   ,(l3+l4)*s(t2+t3)       ,0]
                ])
    return J
                
print jacobian(t1,t2,t3,t4,l1,l2,l3,l4)

#Rotation matrix from world to base<=> coor transform from base to world
Rb=np.matrix([
[c(theta)*c(psi),   s(phi)*s(theta)*c(psi)-c(theta)*s(psi),  s(psi)*s(phi)+c(phi)*s(theta)*c(psi)],
[c(theta)*s(psi),   c(psi)*c(phi)+s(phi)*s(theta)*s(psi),   c(phi)*s(theta)*s(psi)-s(phi)*c(psi)],
[-s(theta),         c(theta)*s(phi),                        c(theta)*c(phi)                     ]
])

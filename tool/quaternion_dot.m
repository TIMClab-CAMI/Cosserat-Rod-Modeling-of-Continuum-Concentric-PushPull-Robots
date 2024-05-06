function dQ = quaternion_dot(Q,Omega)

% 
R = quat2rot(Q);
Omega = R*Omega;

norme_Q=sqrt(Q(1)^2+Q(2)^2+Q(3)^2+Q(4)^2);
Q(1)=Q(1)/norme_Q;
Q(2)=Q(2)/norme_Q;
Q(3)=Q(3)/norme_Q;
Q(4)=Q(4)/norme_Q;
    

A_Omega = [       0, -Omega(1), -Omega(2), -Omega(3)
           Omega(1),         0, -Omega(3),  Omega(2)
           Omega(2),  Omega(3),         0, -Omega(1)
           Omega(3), -Omega(2),  Omega(1),        0];


dQ  = 1/2*A_Omega*Q;

